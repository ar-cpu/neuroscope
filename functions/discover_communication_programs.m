function programs = discover_communication_programs(network, condition_specific, n_programs)
% DISCOVER_COMMUNICATION_PROGRAMS  Decompose the cell-cell communication
%   tensor into latent programs via Non-negative Matrix Factorization (NMF).
%
% Usage:
%   programs = discover_communication_programs(network, condition_specific)
%   programs = discover_communication_programs(network, condition_specific, n_programs)
%   programs = discover_communication_programs(network, condition_specific, 0)  % auto-rank
%
% Inputs:
%   network            - Network struct with fields: source_clusters,
%                        target_clusters, ligands, receptors, weights.
%   condition_specific - Table with columns Ligand, Receptor, and optionally
%                        TherapeuticScore, Pathway.
%   n_programs         - NMF rank. Pass 0 or omit to auto-select via elbow.
%
% Outputs:
%   programs.W, H, V              - Communication matrix and NMF factors
%   programs.rank                 - Selected NMF rank
%   programs.reconstruction_error - Final Frobenius error
%   programs.program_table        - Per-program summary (purity, pathway, score)
%   programs.purity               - [r×1] purity scores
%   programs.inter_program_similarity - [r×r] cosine similarity between programs

fprintf('[Program Discovery] Starting communication program discovery...\n');

out_dir = '';
if exist('results', 'dir')
    out_dir = 'results/';
end

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

fprintf('[Program Discovery] Building communication matrix...\n');

cp_labels = strcat(cellstr(num2str(double(src(:)))), '->', cellstr(num2str(double(tgt(:)))));
[unique_cp, ~, cp_idx] = unique(cp_labels);
n_cp = numel(unique_cp);

lr_labels_raw = strcat(ligs, '_', recs);
[unique_lr, ~, lr_idx] = unique(lr_labels_raw);
n_lr = numel(unique_lr);

fprintf('[Program Discovery] Cluster pairs: %d | L-R pairs: %d\n', n_cp, n_lr);

if n_lr < 3
    warning('[Program Discovery] Fewer than 3 unique L-R pairs (%d). NMF not meaningful.', n_lr);
    return;
end
if n_cp < 4
    warning('[Program Discovery] Fewer than 4 unique cluster pairs (%d). NMF not meaningful.', n_cp);
    return;
end

W = zeros(n_cp, n_lr);
for k = 1:numel(wts)
    W(cp_idx(k), lr_idx(k)) = W(cp_idx(k), lr_idx(k)) + wts(k);
end

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

programs.W                   = W;
programs.cluster_pair_labels = cp_label_list;
programs.lr_pair_labels      = lr_label_list;

max_rank = max(min(8, min(n_lr, n_cp) - 1), 2);

if n_programs > 0
    rank_selected = min(n_programs, max_rank);
    if rank_selected ~= n_programs
        fprintf('[Program Discovery] Requested rank %d exceeds max feasible %d. Using %d.\n', ...
            n_programs, max_rank, rank_selected);
    end
    recon_errors = NaN(max_rank - 1, 1);
else
    fprintf('[Program Discovery] Auto-selecting NMF rank (testing r = 2..%d)...\n', max_rank);
    recon_errors = zeros(max_rank - 1, 1);

    for r_test = 2:max_rank
        try
            [~, ~, err_test] = nnmf(W, r_test, 'algorithm', 'mult', 'replicates', 3);
            recon_errors(r_test - 1) = err_test;
        catch me_nmf
            fprintf('[Program Discovery]   rank %d NMF failed: %s\n', r_test, me_nmf.message);
            recon_errors(r_test - 1) = NaN;
        end
        fprintf('[Program Discovery]   r=%d  recon_error=%.4f\n', r_test, recon_errors(r_test - 1));
    end

    % First r where relative reconstruction improvement falls below 10%
    rank_selected = max_rank;
    for ri = 2:numel(recon_errors)
        if ~isnan(recon_errors(ri - 1)) && ~isnan(recon_errors(ri))
            relative_improvement = (recon_errors(ri - 1) - recon_errors(ri)) / ...
                                   (recon_errors(ri - 1) + eps);
            if relative_improvement < 0.10
                rank_selected = ri + 1;
                break;
            end
        end
    end
    rank_selected = max(min(rank_selected, max_rank), 2);
    fprintf('[Program Discovery] Auto-selected rank r = %d\n', rank_selected);
end

fprintf('[Program Discovery] Running NMF (rank=%d, 10 restarts)...\n', rank_selected);

[H_nmf, V_nmf, final_error] = nnmf(W, rank_selected, 'algorithm', 'mult', 'replicates', 10);

% MATLAB nnmf: W ≈ H_nmf * V_nmf; transpose V to get [n_lr × r] loadings
H = H_nmf;
V = V_nmf';

fprintf('[Program Discovery] Final reconstruction error: %.6f\n', final_error);

programs.H                    = H;
programs.V                    = V;
programs.rank                 = rank_selected;
programs.reconstruction_error = final_error;

fprintf('[Program Discovery] Characterising %d programs...\n', rank_selected);

lr_ts_map = containers.Map('KeyType','char','ValueType','double');
lr_pw_map = containers.Map('KeyType','char','ValueType','char');

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

purity = zeros(rank_selected, 1);
for k = 1:rank_selected
    v_k = V(:, k);
    s   = sum(v_k);
    if s > 0
        purity(k) = max(v_k) / s;
    end
end

programs.purity = purity;

top_lr_cell    = cell(rank_selected, 1);
top_cp_cell    = cell(rank_selected, 1);
ther_score_vec = nan(rank_selected, 1);
pathway_cell   = cell(rank_selected, 1);
n_top = 5;

for k = 1:rank_selected
    v_k = V(:, k);
    [~, sorted_lr_idx] = sort(v_k, 'descend');
    top_lr_indices = sorted_lr_idx(1:min(n_top, numel(lr_label_list)));
    top_lr_cell{k} = strjoin(lr_label_list(top_lr_indices), '; ');

    h_k = H(:, k);
    [~, sorted_cp_idx] = sort(h_k, 'descend');
    top_cp_indices = sorted_cp_idx(1:min(n_top, numel(cp_label_list)));
    top_cp_cell{k} = strjoin(cp_label_list(top_cp_indices), '; ');

    % Therapeutic score: mean over top-quartile L-R pairs by loading
    q75 = quantile(v_k, 0.75);
    top_q_lrs = lr_label_list(v_k >= q75);
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

    % Dominant pathway: most frequent annotation among top L-R pairs
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
        pw_keys = keys(pw_counts);
        pw_vals = cell2mat(values(pw_counts));
        [~, pm_i] = max(pw_vals);
        pathway_cell{k} = pw_keys{pm_i};
    else
        pathway_cell{k} = 'Unknown';
    end
end

program_table = table( ...
    (1:rank_selected)', ...
    top_lr_cell, ...
    top_cp_cell, ...
    purity, ...
    ther_score_vec, ...
    pathway_cell, ...
    'VariableNames', {'ProgramID','TopLRPairs','TopClusterPairs','Purity','TherapeuticScore','Pathway'});

programs.program_table = program_table;

fprintf('[Program Discovery] Computing inter-program cosine similarity...\n');

sim_matrix = zeros(rank_selected, rank_selected);
for i = 1:rank_selected
    for j = 1:rank_selected
        vi = V(:, i);
        vj = V(:, j);
        denom = norm(vi) * norm(vj);
        if denom > 0
            sim_matrix(i, j) = dot(vi, vj) / denom;
        end
    end
end
programs.inter_program_similarity = sim_matrix;

fprintf('[Program Discovery] Exporting results...\n');

try
    writetable(program_table, [out_dir, 'communication_programs.csv']);
    fprintf('[Program Discovery] Saved: %scommunication_programs.csv\n', out_dir);
catch me_csv
    fprintf('[Program Discovery] Warning: could not save communication_programs.csv: %s\n', me_csv.message);
end

try
    prog_col_names = arrayfun(@(k) sprintf('Program_%d', k), 1:rank_selected, 'UniformOutput', false);
    lr_load_table  = array2table(V, 'VariableNames', prog_col_names);
    lr_load_table  = [table(lr_label_list(:), 'VariableNames', {'LR_Pair'}), lr_load_table];
    writetable(lr_load_table, [out_dir, 'program_lr_loadings.csv']);
    fprintf('[Program Discovery] Saved: %sprogram_lr_loadings.csv\n', out_dir);
catch me_csv
    fprintf('[Program Discovery] Warning: could not save program_lr_loadings.csv: %s\n', me_csv.message);
end

try
    cp_load_table = array2table(H, 'VariableNames', prog_col_names);
    cp_load_table = [table(cp_label_list(:), 'VariableNames', {'ClusterPair'}), cp_load_table];
    writetable(cp_load_table, [out_dir, 'program_cluster_loadings.csv']);
    fprintf('[Program Discovery] Saved: %sprogram_cluster_loadings.csv\n', out_dir);
catch me_csv
    fprintf('[Program Discovery] Warning: could not save program_cluster_loadings.csv: %s\n', me_csv.message);
end

fprintf('\n[Program Discovery] NMF rank: %d | Recon error: %.6f | CPs: %d | LR pairs: %d\n', ...
    rank_selected, final_error, n_cp, n_lr);
fprintf('[Program Discovery] Per-program overview:\n');
for k = 1:rank_selected
    fprintf('[Program Discovery]   Program %d | Purity=%.3f | TherapScore=%.3f | Pathway=%s\n', ...
        k, purity(k), ...
        ther_score_vec(k) * (~isnan(ther_score_vec(k))), ...
        pathway_cell{k});
end

catch ME
    fprintf('[Program Discovery] ERROR: %s\n', ME.message);
    fprintf('[Program Discovery] Stack: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    programs = struct();
end

end
