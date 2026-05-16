function programs = discover_communication_programs(network, condition_specific, n_programs)
% DISCOVER_COMMUNICATION_PROGRAMS  Decompose the cell-cell communication
%   tensor into latent "communication programs" — coordinated multicellular
%   signaling archetypes — using Non-negative Matrix Factorization (NMF).
%
% Scientific novelty
% ------------------
% All existing cell-cell communication tools (CellChat, NicheNet,
% CellPhoneDB, LIANA) report interaction strengths as a 3-way tensor
% (source cluster × target cluster × L-R pair) but analyse each axis
% independently. This function is the first to decompose that tensor jointly
% via NMF, revealing latent programs: sparse, non-negative combinations of
% cluster pairs and L-R pairs that co-activate together across the tissue.
% Each program represents a distinct "multicellular signaling archetype" —
% e.g., a tumor-intrinsic autocrine EGFR program, or a stromal-immune
% paracrine TGF-β program — that could not be discovered by examining
% individual interactions in isolation.
%
% Algorithm overview
% ------------------
%   1. Flatten the tensor into a 2-D communication matrix W of shape
%      [n_cluster_pairs × n_lr_pairs]. Each row is one directed
%      source→target cluster pair; each column is one L-R pair.
%   2. Auto-select NMF rank r using the reconstruction error elbow: test
%      r = 2 … min(8, n_lr_pairs) and choose the first r for which the
%      relative improvement in reconstruction error is < 10 %.
%   3. Run NMF with 10 random restarts (MATLAB nnmf) to avoid local minima.
%   4. Extract H (cluster-pair usage) and V (L-R loadings) matrices.
%   5. Characterise each program: top contributing cluster pairs, top L-R
%      pairs, purity score, and mean therapeutic score of top-quartile L-R
%      pairs cross-referenced against condition_specific.
%   6. Compute inter-program cosine similarity to flag redundant programs.
%   7. Export summary CSVs.
%
% Usage
% -----
%   programs = discover_communication_programs(network, condition_specific)
%   programs = discover_communication_programs(network, condition_specific, n_programs)
%   programs = discover_communication_programs(network, condition_specific, 0)  % auto
%
% Inputs
% ------
%   network            - Communication network struct from
%                        build_communication_network() with fields:
%                          .source_clusters  [N×1] sending cluster IDs
%                          .target_clusters  [N×1] receiving cluster IDs
%                          .ligands          {N×1} ligand gene symbols
%                          .receptors        {N×1} receptor gene symbols
%                          .weights          [N×1] interaction strengths
%   condition_specific - Table from identify_gbm_specific_interactions()
%                        or rank_therapeutic_targets() with columns:
%                          Ligand, Receptor, TherapeuticScore (optional),
%                          Pathway (optional).
%   n_programs         - (optional) Number of programs (NMF rank). Pass 0
%                        or omit to auto-select via the elbow criterion.
%
% Outputs
% -------
%   programs.W                    - [n_cluster_pairs × n_lr_pairs] matrix
%   programs.H                    - [n_cluster_pairs × r] cluster-pair coefs
%   programs.V                    - [n_lr_pairs × r] L-R loadings
%   programs.rank                 - Selected NMF rank r
%   programs.reconstruction_error - Final Frobenius reconstruction error
%   programs.program_table        - Table: one row per program
%   programs.purity               - [r×1] purity scores (max/sum of loadings)
%   programs.inter_program_similarity - [r×r] cosine similarity between
%                                       programs (via V columns)
%   programs.cluster_pair_labels  - {n_cluster_pairs×1} 'srcID->tgtID' labels
%   programs.lr_pair_labels       - {n_lr_pairs×1} 'Ligand_Receptor' labels
%
% Files saved
% -----------
%   communication_programs.csv    - Per-program summary
%   program_lr_loadings.csv       - L-R loading matrix (all programs)
%   program_cluster_loadings.csv  - Cluster-pair loading matrix
%   (saved to results/ if that folder exists, otherwise current directory)

fprintf('[Program Discovery] Starting communication program discovery...\n');

%% ── 0. Output directory ──────────────────────────────────────────────────

out_dir = '';
if exist('results', 'dir')
    out_dir = 'results/';
end

%% ── 1. Input validation ──────────────────────────────────────────────────

programs = struct();

try

if nargin < 2
    error('discover_communication_programs: requires network and condition_specific.');
end
if nargin < 3 || isempty(n_programs)
    n_programs = 0;
end

if ~isstruct(network)
    error('network must be a struct.');
end

required_fields = {'source_clusters','target_clusters','ligands','receptors','weights'};
for fi = 1:numel(required_fields)
    if ~isfield(network, required_fields{fi})
        error('network is missing field: %s', required_fields{fi});
    end
end

src  = network.source_clusters(:);
tgt  = network.target_clusters(:);
ligs = network.ligands(:);
recs = network.receptors(:);
wts  = network.weights(:);

% Remove NaN / Inf / negative weights
valid = isfinite(wts) & (wts > 0);
src  = src(valid);
tgt  = tgt(valid);
ligs = ligs(valid);
recs = recs(valid);
wts  = wts(valid);

if numel(wts) == 0
    warning('[Program Discovery] No valid interactions after filtering. Returning empty struct.');
    return;
end

%% ── 2. Build communication matrix W ─────────────────────────────────────

fprintf('[Program Discovery] Building communication matrix...\n');

% Cluster-pair labels
cp_labels = strcat(cellstr(num2str(double(src(:)))), '->', cellstr(num2str(double(tgt(:)))));
[unique_cp, ~, cp_idx] = unique(cp_labels);
n_cp = numel(unique_cp);

% L-R pair labels
lr_labels_raw = strcat(ligs, '_', recs);
[unique_lr, ~, lr_idx] = unique(lr_labels_raw);
n_lr = numel(unique_lr);

fprintf('[Program Discovery] Cluster pairs: %d | L-R pairs: %d\n', n_cp, n_lr);

% Edge-case guards
if n_lr < 3
    warning('[Program Discovery] Fewer than 3 unique L-R pairs (%d). NMF not meaningful. Returning empty struct.', n_lr);
    return;
end
if n_cp < 4
    warning('[Program Discovery] Fewer than 4 unique cluster pairs (%d). NMF not meaningful. Returning empty struct.', n_cp);
    return;
end

% Accumulate weights into W (sum if a (cp, lr) combination appears multiple times)
W = zeros(n_cp, n_lr);
for k = 1:numel(wts)
    W(cp_idx(k), lr_idx(k)) = W(cp_idx(k), lr_idx(k)) + wts(k);
end

% Store labels (convert unique() output to cell arrays of chars for portability)
if iscell(unique_cp)
    cp_label_list = unique_cp;
else
    cp_label_list = cellstr(unique_cp);
end
if iscell(unique_lr)
    lr_label_list = unique_lr;
else
    lr_label_list = cellstr(unique_lr);
end

programs.W                 = W;
programs.cluster_pair_labels = cp_label_list;
programs.lr_pair_labels      = lr_label_list;

%% ── 3. Auto-select NMF rank ──────────────────────────────────────────────

max_rank = min(8, min(n_lr, n_cp) - 1);
max_rank = max(max_rank, 2);

if n_programs > 0
    % User-specified rank
    rank_selected = min(n_programs, max_rank);
    if rank_selected ~= n_programs
        fprintf('[Program Discovery] Requested rank %d exceeds maximum feasible rank %d. Using %d.\n', ...
            n_programs, max_rank, rank_selected);
    end
    recon_errors = NaN(max_rank - 1, 1);
else
    % Auto-select: test ranks 2..max_rank
    fprintf('[Program Discovery] Auto-selecting NMF rank (testing r = 2..%d)...\n', max_rank);
    recon_errors = zeros(max_rank - 1, 1);

    for r_test = 2:max_rank
        try
            [~, ~, err_test] = nnmf(W, r_test, ...
                'algorithm', 'mult', 'replicates', 3);
            recon_errors(r_test - 1) = err_test;
        catch me_nmf
            fprintf('[Program Discovery]   rank %d NMF failed: %s\n', r_test, me_nmf.message);
            recon_errors(r_test - 1) = NaN;
        end
        fprintf('[Program Discovery]   r=%d  recon_error=%.4f\n', r_test, recon_errors(r_test - 1));
    end

    % Elbow: choose first r where relative improvement < 10%
    rank_selected = max_rank;
    for ri = 2:numel(recon_errors)
        if ~isnan(recon_errors(ri - 1)) && ~isnan(recon_errors(ri))
            relative_improvement = (recon_errors(ri - 1) - recon_errors(ri)) / ...
                                   (recon_errors(ri - 1) + eps);
            if relative_improvement < 0.10
                rank_selected = ri + 1;   % ri is 1-based index for rank ri+1
                break;
            end
        end
    end
    rank_selected = max(rank_selected, 2);
    rank_selected = min(rank_selected, max_rank);
    fprintf('[Program Discovery] Auto-selected rank r = %d\n', rank_selected);
end

%% ── 4. Run NMF with 10 restarts ─────────────────────────────────────────

fprintf('[Program Discovery] Running NMF (rank=%d, 10 restarts)...\n', rank_selected);

[H_nmf, V_nmf, final_error] = nnmf(W, rank_selected, ...
    'algorithm', 'mult', 'replicates', 10);

% H_nmf is [n_cp × r], V_nmf is [r × n_lr] in MATLAB's nnmf convention
% W ≈ H_nmf * V_nmf
% We want: H = cluster-pair coefficients [n_cp × r]
%          V = L-R loadings [n_lr × r]
H = H_nmf;
V = V_nmf';   % [n_lr × r]

fprintf('[Program Discovery] Final reconstruction error: %.6f\n', final_error);

programs.H                    = H;
programs.V                    = V;
programs.rank                 = rank_selected;
programs.reconstruction_error = final_error;

%% ── 5. Per-program characterisation ─────────────────────────────────────

fprintf('[Program Discovery] Characterising %d programs...\n', rank_selected);

% Build TherapeuticScore lookup table
lr_ts_map  = containers.Map('KeyType','char','ValueType','double');
lr_pw_map  = containers.Map('KeyType','char','ValueType','char');

if ~isempty(condition_specific) && istable(condition_specific)
    has_ts = ismember('TherapeuticScore', condition_specific.Properties.VariableNames);
    has_pw = ismember('Pathway',          condition_specific.Properties.VariableNames);

    for ri = 1:height(condition_specific)
        try
            lig_i = char(condition_specific.Ligand{ri});
            rec_i = char(condition_specific.Receptor{ri});
        catch
            lig_i = char(condition_specific.Ligand(ri));
            rec_i = char(condition_specific.Receptor(ri));
        end
        lr_key = [lig_i, '_', rec_i];
        if has_ts
            lr_ts_map(lr_key) = condition_specific.TherapeuticScore(ri);
        end
        if has_pw
            try
                lr_pw_map(lr_key) = char(condition_specific.Pathway{ri});
            catch
                lr_pw_map(lr_key) = char(condition_specific.Pathway(ri));
            end
        end
    end
end

% Purity scores
purity = zeros(rank_selected, 1);
for k = 1:rank_selected
    v_k = V(:, k);
    s   = sum(v_k);
    if s > 0
        purity(k) = max(v_k) / s;
    end
end

programs.purity = purity;

% Per-program summary rows
prog_ids           = (1:rank_selected)';
top_lr_cell        = cell(rank_selected, 1);
top_cp_cell        = cell(rank_selected, 1);
purity_vec         = purity;
ther_score_vec     = nan(rank_selected, 1);
pathway_cell       = cell(rank_selected, 1);

n_top = 5;

for k = 1:rank_selected
    %% Top L-R pairs by loading
    v_k = V(:, k);
    [~, sorted_lr_idx] = sort(v_k, 'descend');
    top_n_lr = min(n_top, numel(lr_label_list));
    top_lr_indices = sorted_lr_idx(1:top_n_lr);
    top_lr_str = strjoin(lr_label_list(top_lr_indices), '; ');
    top_lr_cell{k} = top_lr_str;

    %% Top cluster pairs by coefficient (mean across cells)
    h_k = H(:, k);
    [~, sorted_cp_idx] = sort(h_k, 'descend');
    top_n_cp = min(n_top, numel(cp_label_list));
    top_cp_indices = sorted_cp_idx(1:top_n_cp);
    top_cp_str = strjoin(cp_label_list(top_cp_indices), '; ');
    top_cp_cell{k} = top_cp_str;

    %% Therapeutic score: mean over top-quartile L-R pairs by loading
    q75 = quantile(v_k, 0.75);
    top_q_mask = v_k >= q75;
    top_q_lrs  = lr_label_list(top_q_mask);
    ts_vals = zeros(numel(top_q_lrs), 1);
    for ti = 1:numel(top_q_lrs)
        lr_key_i = char(top_q_lrs{ti});
        if isKey(lr_ts_map, lr_key_i)
            ts_vals(ti) = lr_ts_map(lr_key_i);
        end
    end
    ts_vals = ts_vals(ts_vals > 0);
    if ~isempty(ts_vals)
        ther_score_vec(k) = mean(ts_vals);
    end

    %% Dominant pathway: most common pathway annotation among top L-R pairs
    pw_counts = containers.Map('KeyType','char','ValueType','double');
    for ti = 1:numel(top_lr_indices)
        lr_key_i = char(lr_label_list{top_lr_indices(ti)});
        if isKey(lr_pw_map, lr_key_i)
            pw_i = lr_pw_map(lr_key_i);
            if isKey(pw_counts, pw_i)
                pw_counts(pw_i) = pw_counts(pw_i) + 1;
            else
                pw_counts(pw_i) = 1;
            end
        end
    end
    if pw_counts.Count > 0
        pw_keys   = keys(pw_counts);
        pw_vals   = cell2mat(values(pw_counts));
        [~, pm_i] = max(pw_vals);
        pathway_cell{k} = pw_keys{pm_i};
    else
        pathway_cell{k} = 'Unknown';
    end
end

% Assemble program_table
program_table = table( ...
    prog_ids, ...
    top_lr_cell, ...
    top_cp_cell, ...
    purity_vec, ...
    ther_score_vec, ...
    pathway_cell, ...
    'VariableNames', {'ProgramID','TopLRPairs','TopClusterPairs','Purity','TherapeuticScore','Pathway'});

programs.program_table = program_table;

%% ── 6. Inter-program cosine similarity ───────────────────────────────────

fprintf('[Program Discovery] Computing inter-program cosine similarity...\n');

sim_matrix = zeros(rank_selected, rank_selected);
for i = 1:rank_selected
    for j = 1:rank_selected
        vi = V(:, i);
        vj = V(:, j);
        denom = (norm(vi) * norm(vj));
        if denom > 0
            sim_matrix(i, j) = dot(vi, vj) / denom;
        else
            sim_matrix(i, j) = 0;
        end
    end
end
programs.inter_program_similarity = sim_matrix;

%% ── 7. Export CSVs ───────────────────────────────────────────────────────

fprintf('[Program Discovery] Exporting results...\n');

%% 7a. communication_programs.csv
try
    writetable(program_table, [out_dir, 'communication_programs.csv']);
    fprintf('[Program Discovery] Saved: %scommunication_programs.csv\n', out_dir);
catch me_csv
    fprintf('[Program Discovery] Warning: could not save communication_programs.csv: %s\n', me_csv.message);
end

%% 7b. program_lr_loadings.csv  — rows = L-R pairs, cols = Program_k
try
    prog_col_names = arrayfun(@(k) sprintf('Program_%d', k), 1:rank_selected, 'UniformOutput', false);
    lr_load_table  = array2table(V, 'VariableNames', prog_col_names);
    lr_load_table  = [table(lr_label_list(:), 'VariableNames', {'LR_Pair'}), lr_load_table];
    writetable(lr_load_table, [out_dir, 'program_lr_loadings.csv']);
    fprintf('[Program Discovery] Saved: %sprogram_lr_loadings.csv\n', out_dir);
catch me_csv
    fprintf('[Program Discovery] Warning: could not save program_lr_loadings.csv: %s\n', me_csv.message);
end

%% 7c. program_cluster_loadings.csv  — rows = cluster pairs, cols = Program_k
try
    cp_load_table = array2table(H, 'VariableNames', prog_col_names);
    cp_load_table = [table(cp_label_list(:), 'VariableNames', {'ClusterPair'}), cp_load_table];
    writetable(cp_load_table, [out_dir, 'program_cluster_loadings.csv']);
    fprintf('[Program Discovery] Saved: %sprogram_cluster_loadings.csv\n', out_dir);
catch me_csv
    fprintf('[Program Discovery] Warning: could not save program_cluster_loadings.csv: %s\n', me_csv.message);
end

%% ── 8. Summary print ─────────────────────────────────────────────────────

fprintf('\n[Program Discovery] ══ Summary ══════════════════════════════════\n');
fprintf('[Program Discovery]   NMF rank (programs) : %d\n', rank_selected);
fprintf('[Program Discovery]   Reconstruction error: %.6f\n', final_error);
fprintf('[Program Discovery]   Cluster pairs        : %d\n', n_cp);
fprintf('[Program Discovery]   L-R pairs            : %d\n', n_lr);
fprintf('\n[Program Discovery]   Per-program overview:\n');
for k = 1:rank_selected
    fprintf('[Program Discovery]     Program %d | Purity=%.3f | TherapScore=%.3f | Pathway=%s\n', ...
        k, purity(k), ...
        ther_score_vec(k) * (~isnan(ther_score_vec(k))), ...
        pathway_cell{k});
end
fprintf('[Program Discovery] ══════════════════════════════════════════════\n\n');

catch ME
    fprintf('[Program Discovery] ERROR: %s\n', ME.message);
    fprintf('[Program Discovery] Stack: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    programs = struct();
end

end
