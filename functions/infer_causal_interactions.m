function causal_results = infer_causal_interactions(condition_specific, expr_matrix, gene_names, cell_metadata, network)
% INFER_CAUSAL_INTERACTIONS  2SLS instrumental variable causal inference for
%   cell-cell communication. Uses donor_id as a natural instrument (analogous
%   to Mendelian Randomization) to distinguish true L→R signalling from
%   spurious co-expression. Validates: first-stage F > 10, second-stage p < 0.05.
%
% causal_results = infer_causal_interactions(condition_specific, expr_matrix,
%                     gene_names, cell_metadata, network)
%
% Inputs:
%   condition_specific - Table with Ligand, Receptor columns
%   expr_matrix        - [cells × genes] sparse or dense expression matrix
%   gene_names         - cell array, string array, or table with 'gene_name'
%   cell_metadata      - Table with cluster_id and optionally donor_id columns
%   network            - Struct with .source_clusters, .target_clusters,
%                        .ligands, .receptors
%
% Outputs:
%   causal_results.causal_table    - Per-pair results with CausalEffect, F, p
%   causal_results.n_causal        - Causally validated pairs
%   causal_results.instrument_type - 'donor' or 'cluster_quantile'
%   causal_results.causal_network  - Sub-network of causal pairs only
%   causal_results.top_causal      - Top 5 by |CausalEffect|

fprintf('[Causal Inference] Starting causal inference for cell-cell communication...\n');

% Output directory

out_dir = '';
if exist('results', 'dir')
    out_dir = 'results/';
end

causal_results = struct();

try

% Input validation

if nargin < 5
    error('infer_causal_interactions: all 5 inputs are required.');
end

if ~istable(condition_specific) || height(condition_specific) == 0
    warning('[Causal Inference] condition_specific is empty or not a table. Returning empty struct.');
    return;
end

required_cs_cols = {'Ligand','Receptor'};
for fi = 1:numel(required_cs_cols)
    if ~ismember(required_cs_cols{fi}, condition_specific.Properties.VariableNames)
        error('[Causal Inference] condition_specific is missing column: %s', required_cs_cols{fi});
    end
end

if ~isstruct(network)
    error('[Causal Inference] network must be a struct.');
end
for fi = {'source_clusters','target_clusters','ligands','receptors'}
    if ~isfield(network, fi{1})
        error('[Causal Inference] network missing field: %s', fi{1});
    end
end

% Normalise gene_names

fprintf('[Causal Inference] Parsing gene names...\n');

if istable(gene_names)
    if ismember('gene_name', gene_names.Properties.VariableNames)
        gene_names_cell = cellstr(gene_names.gene_name);
    else
        gene_names_cell = cellstr(gene_names{:, 1});
    end
elseif isstring(gene_names)
    gene_names_cell = cellstr(gene_names);
elseif iscell(gene_names)
    gene_names_cell = cellfun(@char, gene_names, 'UniformOutput', false);
else
    error('[Causal Inference] gene_names must be a cell array, string array, or table.');
end
gene_names_cell = gene_names_cell(:);

n_cells = size(expr_matrix, 1);
n_genes = size(expr_matrix, 2);

if numel(gene_names_cell) ~= n_genes
    error('[Causal Inference] Length of gene_names (%d) does not match expr_matrix columns (%d).', ...
        numel(gene_names_cell), n_genes);
end

% Extract cluster IDs

fprintf('[Causal Inference] Extracting cluster IDs from cell_metadata...\n');

if ~istable(cell_metadata)
    error('[Causal Inference] cell_metadata must be a table.');
end
if size(cell_metadata, 1) ~= n_cells
    error('[Causal Inference] cell_metadata rows (%d) != expr_matrix rows (%d).', ...
        size(cell_metadata, 1), n_cells);
end

cluster_col = '';
meta_cols   = cell_metadata.Properties.VariableNames;
preferred_cluster_cols = {'Cluster','cluster_id','ClusterID','cluster','CellType','cell_type'};
for ci = 1:numel(preferred_cluster_cols)
    if ismember(preferred_cluster_cols{ci}, meta_cols)
        cluster_col = preferred_cluster_cols{ci};
        break;
    end
end
if isempty(cluster_col)
    % partial-name fallback
    for ci = 1:numel(meta_cols)
        if ~isempty(regexpi(meta_cols{ci}, 'cluster', 'once'))
            cluster_col = meta_cols{ci};
            break;
        end
    end
end
if isempty(cluster_col)
    error('[Causal Inference] Cannot find a cluster column in cell_metadata.');
end

raw_cluster_ids = cell_metadata.(cluster_col);
if iscategorical(raw_cluster_ids)
    raw_cluster_ids = double(raw_cluster_ids);
elseif iscell(raw_cluster_ids)
    raw_cluster_ids = cellfun(@(x) str2double(char(x)), raw_cluster_ids);
end
raw_cluster_ids = raw_cluster_ids(:);

% Identify instrument variable

fprintf('[Causal Inference] Identifying instrument variable...\n');

donor_col_candidates = {'donor_id','donor','DonorID','sample_id','SampleID','Sample','sample'};
donor_col = '';
for ci = 1:numel(donor_col_candidates)
    if ismember(donor_col_candidates{ci}, meta_cols)
        donor_col = donor_col_candidates{ci};
        break;
    end
end

instrument_type = 'cluster_quantile';

if ~isempty(donor_col)
    raw_donor = cell_metadata.(donor_col);
    if iscategorical(raw_donor)
        donor_ids_all = cellstr(raw_donor);
    elseif isnumeric(raw_donor)
        donor_ids_all = cellstr(num2str(raw_donor(:)));
    else
        donor_ids_all = cellstr(raw_donor);
    end
    donor_ids_all = donor_ids_all(:);
    unique_donors = unique(donor_ids_all);
    n_donors = numel(unique_donors);

    if n_donors >= 2
        instrument_type = 'donor';
        fprintf('[Causal Inference] Instrument: donor_id (%d unique donors found in column "%s").\n', ...
            n_donors, donor_col);
    else
        fprintf('[Causal Inference] Only 1 donor found. Falling back to cluster_quantile instrument.\n');
        donor_col = '';
    end
else
    fprintf('[Causal Inference] No donor column found. Using cluster_quantile instrument.\n');
end

causal_results.instrument_type = instrument_type;

% Build unique L-R pair list

fprintf('[Causal Inference] Building unique L-R pair list from condition_specific...\n');

try
    cs_ligs = cellfun(@char, condition_specific.Ligand,   'UniformOutput', false);
    cs_recs = cellfun(@char, condition_specific.Receptor, 'UniformOutput', false);
catch
    cs_ligs = cellstr(condition_specific.Ligand);
    cs_recs = cellstr(condition_specific.Receptor);
end

cs_sigs = strcat(cs_ligs, '_', cs_recs);
[unique_lr_sigs, ia] = unique(cs_sigs, 'stable');
unique_ligs = cs_ligs(ia);
unique_recs = cs_recs(ia);
n_pairs = numel(unique_lr_sigs);

fprintf('[Causal Inference] Testing %d unique L-R pairs...\n', n_pairs);

% Build network lookup: for each L-R pair, which source/target clusters?
net_ligs = cellfun(@char, network.ligands, 'UniformOutput', false);
net_recs = cellfun(@char, network.receptors, 'UniformOutput', false);
net_src  = network.source_clusters(:);
net_tgt  = network.target_clusters(:);
net_sigs = strcat(net_ligs, '_', net_recs);

% Correlational strength lookup

has_weight = ismember('Weight', condition_specific.Properties.VariableNames);
cs_weight_map = containers.Map('KeyType','char','ValueType','double');
if has_weight
    for ri = 1:height(condition_specific)
        cs_weight_map(cs_sigs{ri}) = condition_specific.Weight(ri);
    end
end

% Per-pair 2SLS causal inference

causal_effect_vec    = nan(n_pairs, 1);
ci_lo_vec            = nan(n_pairs, 1);
ci_hi_vec            = nan(n_pairs, 1);
f_stat_vec           = nan(n_pairs, 1);
p_val_vec            = nan(n_pairs, 1);
is_weak_vec          = false(n_pairs, 1);
is_causal_vec        = false(n_pairs, 1);
corr_strength_vec    = nan(n_pairs, 1);

n_bootstrap = 500;

fprintf('[Causal Inference] Running 2SLS per pair (bootstrap=%d)...\n', n_bootstrap);

for pi = 1:n_pairs
    lig_i = unique_ligs{pi};
    rec_i = unique_recs{pi};
    sig_i = unique_lr_sigs{pi};

    if mod(pi, max(1, floor(n_pairs/10))) == 0
        fprintf('[Causal Inference]   Progress: %d/%d pairs processed...\n', pi, n_pairs);
    end

    try
        % find gene indices
        lig_idx = find(strcmp(gene_names_cell, lig_i), 1);
        rec_idx = find(strcmp(gene_names_cell, rec_i), 1);

        if isempty(lig_idx) || isempty(rec_idx)
            continue;   % gene not in expression matrix
        end

        % identify source/target cells via network cluster membership
        net_mask_i = strcmp(net_sigs, sig_i);
        src_clusters_i = unique(net_src(net_mask_i));
        tgt_clusters_i = unique(net_tgt(net_mask_i));

        if isempty(src_clusters_i) || isempty(tgt_clusters_i)
            % Fallback: all cells are source and target
            src_cell_mask = true(n_cells, 1);
            tgt_cell_mask = true(n_cells, 1);
        else
            src_cell_mask = ismember(raw_cluster_ids, src_clusters_i);
            tgt_cell_mask = ismember(raw_cluster_ids, tgt_clusters_i);
        end

        src_cell_idx = find(src_cell_mask);
        tgt_cell_idx = find(tgt_cell_mask);

        if numel(src_cell_idx) < 10 || numel(tgt_cell_idx) < 10
            continue;   % insufficient cells
        end

        % extract expression
        if issparse(expr_matrix)
            X_lig = full(expr_matrix(src_cell_idx, lig_idx));
            X_rec = full(expr_matrix(tgt_cell_idx, rec_idx));
        else
            X_lig = expr_matrix(src_cell_idx, lig_idx);
            X_rec = expr_matrix(tgt_cell_idx, rec_idx);
        end
        X_lig = double(X_lig(:));
        X_rec = double(X_rec(:));

        % Correlational strength (Spearman rho)
        if numel(X_lig) >= 5 && numel(X_rec) >= 5
            min_n = min(numel(X_lig), numel(X_rec));
            try
                rho = corr(X_lig(1:min_n), X_rec(1:min_n), 'Type', 'Spearman');
                corr_strength_vec(pi) = rho;
            catch
                corr_strength_vec(pi) = NaN;
            end
        end
        if has_weight && isKey(cs_weight_map, sig_i)
            corr_strength_vec(pi) = cs_weight_map(sig_i);
        end

        % build instrument matrix
        n_src = numel(src_cell_idx);

        if strcmp(instrument_type, 'donor')
            donor_src = donor_ids_all(src_cell_idx);
            Z_src = build_onehot(donor_src, unique_donors);
            k_donors = size(Z_src, 2);
        else
            % cluster-quantile fallback instrument
            Z_src = build_quantile_instrument(X_lig, 3);
            k_donors = size(Z_src, 2);
        end

        % first stage: regress X_lig on Z_src
        if size(Z_src, 1) < k_donors + 2
            continue;   % not enough observations for first stage
        end

        % Add intercept to Z for R² calculation
        Z_aug = [ones(n_src, 1), Z_src];
        try
            beta_fs  = Z_aug \ X_lig;
            Xhat_lig = Z_aug * beta_fs;
        catch
            Xhat_lig = mean(X_lig) * ones(n_src, 1);
        end

        % R² and F-stat of first stage
        ss_res  = sum((X_lig - Xhat_lig).^2);
        ss_tot  = sum((X_lig - mean(X_lig)).^2);
        r2_fs   = max(0, 1 - ss_res / (ss_tot + eps));

        k_fs   = k_donors;       % number of instruments (excl intercept)
        n_fs   = n_src;
        if k_fs >= 1 && (n_fs - k_fs - 1) > 0
            F_stat = (r2_fs / k_fs) / ((1 - r2_fs) / (n_fs - k_fs - 1) + eps);
        else
            F_stat = 0;
        end
        f_stat_vec(pi)  = F_stat;
        is_weak_vec(pi) = F_stat < 10;

        % second stage: regress X_rec on X̂_lig projected to target cells

        n_tgt = numel(tgt_cell_idx);

        if strcmp(instrument_type, 'donor')
            donor_tgt = donor_ids_all(tgt_cell_idx);
            % per-donor mean fitted ligand
            fitted_per_donor = zeros(numel(unique_donors), 1);
            for di = 1:numel(unique_donors)
                src_di_mask = strcmp(donor_src, unique_donors{di});
                if any(src_di_mask)
                    fitted_per_donor(di) = mean(Xhat_lig(src_di_mask));
                else
                    fitted_per_donor(di) = mean(Xhat_lig);
                end
            end
            Xhat_lig_tgt = zeros(n_tgt, 1);
            for di = 1:numel(unique_donors)
                tgt_di_mask = strcmp(donor_tgt, unique_donors{di});
                Xhat_lig_tgt(tgt_di_mask) = fitted_per_donor(di);
            end
        else
            Xhat_lig_tgt = mean(Xhat_lig) * ones(n_tgt, 1);
        end

        % fitlm second stage
        if numel(unique(Xhat_lig_tgt)) < 2
            continue;   % degenerate instrument
        end

        try
            lm = fitlm(Xhat_lig_tgt, X_rec);
            beta_ss   = lm.Coefficients.Estimate(2);
            p_val_ss  = lm.Coefficients.pValue(2);
        catch
            continue;
        end

        causal_effect_vec(pi) = beta_ss;
        p_val_vec(pi)         = p_val_ss;

        % bootstrap 95% CI
        boot_betas = nan(n_bootstrap, 1);

        if strcmp(instrument_type, 'donor')
            for bi = 1:n_bootstrap
                % block bootstrap
                boot_src_idx = block_bootstrap_idx(donor_src, unique_donors, n_src);
                boot_tgt_idx = block_bootstrap_idx(donor_tgt, unique_donors, n_tgt);

                Xlig_b = X_lig(boot_src_idx);
                Xrec_b = X_rec(boot_tgt_idx);
                Zb     = Z_src(boot_src_idx, :);
                Zb_aug = [ones(n_src, 1), Zb];

                try
                    beta_b_fs  = Zb_aug \ Xlig_b;
                    Xhat_b     = Zb_aug * beta_b_fs;

                    donor_src_b = donor_src(boot_src_idx);
                    donor_tgt_b = donor_tgt(boot_tgt_idx);
                    fpd_b = zeros(numel(unique_donors), 1);
                    for di = 1:numel(unique_donors)
                        m_di = strcmp(donor_src_b, unique_donors{di});
                        if any(m_di)
                            fpd_b(di) = mean(Xhat_b(m_di));
                        else
                            fpd_b(di) = mean(Xhat_b);
                        end
                    end
                    Xhat_tgt_b = zeros(n_tgt, 1);
                    for di = 1:numel(unique_donors)
                        m_di_t = strcmp(donor_tgt_b, unique_donors{di});
                        Xhat_tgt_b(m_di_t) = fpd_b(di);
                    end

                    if numel(unique(Xhat_tgt_b)) >= 2
                        lm_b = fitlm(Xhat_tgt_b, Xrec_b);
                        boot_betas(bi) = lm_b.Coefficients.Estimate(2);
                    end
                catch
                    % skip failed bootstrap iteration
                end
            end
        else
            % cluster-quantile bootstrap
            for bi = 1:n_bootstrap
                boot_src_i = randi(n_src, n_src, 1);
                boot_tgt_i = randi(n_tgt, n_tgt, 1);
                Xlig_b = X_lig(boot_src_i);
                Xrec_b = X_rec(boot_tgt_i);
                Zb_raw = build_quantile_instrument(Xlig_b, 3);
                Zb_aug = [ones(n_src, 1), Zb_raw];
                try
                    beta_b_fs  = Zb_aug \ Xlig_b;
                    Xhat_b     = Zb_aug * beta_b_fs;
                    Xhat_sc    = mean(Xhat_b) * ones(n_tgt, 1);
                    if numel(unique(Xhat_sc)) >= 2
                        lm_b = fitlm(Xhat_sc, Xrec_b);
                        boot_betas(bi) = lm_b.Coefficients.Estimate(2);
                    end
                catch
                    % skip
                end
            end
        end

        valid_boots = boot_betas(~isnan(boot_betas));
        if numel(valid_boots) >= 10
            ci_lo_vec(pi) = quantile(valid_boots, 0.025);
            ci_hi_vec(pi) = quantile(valid_boots, 0.975);
        else
            ci_lo_vec(pi) = beta_ss;
            ci_hi_vec(pi) = beta_ss;
        end

        % classify
        is_causal_vec(pi) = (~is_weak_vec(pi)) && (p_val_ss < 0.05);

    catch ME_pair
        fprintf('[Causal Inference]   Pair %s: skipped (%s)\n', sig_i, ME_pair.message);
    end
end

fprintf('[Causal Inference] 2SLS complete. Assembling results...\n');

% Assemble output table

n_tested  = sum(~isnan(causal_effect_vec));
n_causal  = sum(is_causal_vec);
causal_fraction = n_causal / max(n_tested, 1);

causal_table = table( ...
    unique_ligs(:), ...
    unique_recs(:), ...
    causal_effect_vec, ...
    ci_lo_vec, ...
    ci_hi_vec, ...
    f_stat_vec, ...
    p_val_vec, ...
    is_weak_vec, ...
    is_causal_vec, ...
    corr_strength_vec, ...
    'VariableNames', { ...
        'Ligand','Receptor', ...
        'CausalEffect','CausalCI_lo','CausalCI_hi', ...
        'InstrumentF','PValue', ...
        'IsWeakInstrument','IsCausal','CorrelationalStrength'});

causal_results.causal_table    = causal_table;
causal_results.n_causal        = n_causal;
causal_results.n_tested        = n_tested;
causal_results.causal_fraction = causal_fraction;
causal_results.instrument_type = instrument_type;

% Causal sub-network

causal_mask_net = false(numel(net_ligs), 1);
for pi = 1:n_pairs
    if is_causal_vec(pi)
        sig_i = unique_lr_sigs{pi};
        causal_mask_net = causal_mask_net | strcmp(net_sigs, sig_i);
    end
end

causal_network.ligands          = network.ligands(causal_mask_net);
causal_network.receptors        = network.receptors(causal_mask_net);
causal_network.source_clusters  = network.source_clusters(causal_mask_net);
causal_network.target_clusters  = network.target_clusters(causal_mask_net);
if isfield(network, 'weights')
    causal_network.weights = network.weights(causal_mask_net);
end
causal_results.causal_network = causal_network;

% Top 5 causal pairs by |CausalEffect|

causal_only = causal_table(causal_table.IsCausal, :);
if ~isempty(causal_only)
    [~, top_ord] = sort(abs(causal_only.CausalEffect), 'descend');
    top_n = min(5, height(causal_only));
    causal_results.top_causal = causal_only(top_ord(1:top_n), :);
else
    causal_results.top_causal = causal_only;
end

% Export CSV

try
    writetable(causal_table, [out_dir, 'causal_interactions.csv']);
    fprintf('[Causal Inference] Saved: %scausal_interactions.csv\n', out_dir);
catch me_csv
    fprintf('[Causal Inference] Warning: could not save causal_interactions.csv: %s\n', me_csv.message);
end

% Summary

fprintf('\n[Causal Inference] ══ Summary ════════════════════════════════════\n');
fprintf('[Causal Inference]   Instrument type       : %s\n', instrument_type);
fprintf('[Causal Inference]   L-R pairs tested      : %d\n', n_tested);
fprintf('[Causal Inference]   Causally validated     : %d (%.1f%%)\n', n_causal, 100*causal_fraction);
fprintf('[Causal Inference]   Weak instruments (F<10): %d\n', sum(is_weak_vec));
if ~isempty(causal_results.top_causal)
    fprintf('[Causal Inference]   Top causal pairs:\n');
    for ti = 1:height(causal_results.top_causal)
        fprintf('[Causal Inference]     %s → %s  β=%.4f  F=%.2f  p=%.4f\n', ...
            char(causal_results.top_causal.Ligand{ti}), ...
            char(causal_results.top_causal.Receptor{ti}), ...
            causal_results.top_causal.CausalEffect(ti), ...
            causal_results.top_causal.InstrumentF(ti), ...
            causal_results.top_causal.PValue(ti));
    end
end
fprintf('[Causal Inference] ══════════════════════════════════════════════\n\n');

catch ME
    fprintf('[Causal Inference] ERROR: %s\n', ME.message);
    fprintf('[Causal Inference] Stack: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    causal_results = struct();
end

end   % main function

function Z = build_onehot(labels, unique_labels)
% One-hot encoding of labels; returns [n × k] binary matrix.
    n = numel(labels);
    k = numel(unique_labels);
    Z = zeros(n, k);
    for li = 1:k
        Z(:, li) = strcmp(labels, unique_labels{li});
    end
end

function Z = build_quantile_instrument(X, n_quantiles)
% Fallback instrument: assign cells to quantile groups of X, return one-hot matrix.
    n = numel(X);
    edges = quantile(X, linspace(0, 1, n_quantiles + 1));
    edges(1)   = -Inf;
    edges(end) = Inf;
    group = zeros(n, 1);
    for qi = 1:n_quantiles
        mask = X > edges(qi) & X <= edges(qi + 1);
        group(mask) = qi;
    end
    group(group == 0) = 1;   % assign zeros to group 1
    Z = zeros(n, n_quantiles);
    for qi = 1:n_quantiles
        Z(:, qi) = (group == qi);
    end
end

function boot_idx = block_bootstrap_idx(donor_labels, unique_donors, n_total)
% Resample cell indices within each donor block (preserves donor structure).
    boot_idx = zeros(n_total, 1);
    ptr = 1;
    for di = 1:numel(unique_donors)
        di_idx = find(strcmp(donor_labels, unique_donors{di}));
        if isempty(di_idx)
            continue;
        end
        n_di = numel(di_idx);
        resample_pos = di_idx(randi(n_di, n_di, 1));
        end_ptr = min(ptr + n_di - 1, n_total);
        actual_n = end_ptr - ptr + 1;
        boot_idx(ptr:end_ptr) = resample_pos(1:actual_n);
        ptr = end_ptr + 1;
        if ptr > n_total
            break;
        end
    end
    % Fill any unfilled slots
    if ptr <= n_total
        boot_idx(ptr:end) = randi(n_total, n_total - ptr + 1, 1);
    end
    % Replace any zero entries (safety)
    boot_idx(boot_idx == 0) = 1;
end
