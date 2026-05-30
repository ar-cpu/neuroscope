function network = build_communication_network(expr_matrix, genes_table, clusters_table, cellphonedb, dataset_name, varargin)
% BUILD_COMMUNICATION_NETWORK  Construct a weighted directed ligand-receptor
%   communication network from single-cell expression data.
%
% Usage:
%   network = build_communication_network(expr_matrix, genes_table, clusters_table, cellphonedb, dataset_name)
%
% Inputs:
%   expr_matrix    - [cells × genes] sparse or dense expression matrix (raw/normalized counts)
%   genes_table    - Table with 'gene_name' column, or cell array of gene symbols
%   clusters_table - Table with 'Cluster' column (numeric IDs, 1-based)
%   cellphonedb    - CellPhoneDB struct or {N×2} cell array of [ligand, receptor] pairs
%   dataset_name   - String label used in progress messages
%
% Outputs:
%   network.ligands, .receptors, .weights, .source_clusters, .target_clusters,
%   .mi_scores, .significance, .raw_strength, .channel_capacity, .ci_lower,
%   .ci_upper, .mi_mode  — see field names for descriptions

fprintf('Network construction: %s\n', dataset_name);

gene_names = genes_table.gene_name;
cluster_ids = clusters_table.Cluster;
if iscategorical(cluster_ids)
    cluster_ids = double(cluster_ids);
end

unique_clusters = unique(cluster_ids);
n_clusters = length(unique_clusters);
n_cells = size(expr_matrix, 1);
n_genes = length(gene_names);

% Look for donor/sample column to enable cross-sample KSG MI.
sample_id_col = '';
for cand = {'donor_id', 'sample_id', 'donor', 'sample', 'SubjectID', 'individual'}
    if ismember(cand{1}, clusters_table.Properties.VariableNames)
        sample_id_col = cand{1};
        break;
    end
end

if ~isempty(sample_id_col)
    raw_sample_ids = clusters_table.(sample_id_col);
    if isnumeric(raw_sample_ids) || islogical(raw_sample_ids)
        sample_ids = arrayfun(@num2str, raw_sample_ids, 'UniformOutput', false);
    else
        sample_ids = cellstr(raw_sample_ids);
    end
    unique_samples = unique(sample_ids);
    n_samples = length(unique_samples);
    fprintf('  Cells: %d | Genes: %d | Clusters: %d | Samples: %d\n', ...
            n_cells, n_genes, n_clusters, n_samples);
    fprintf('  MI mode: cross-sample (paired per-sample means)\n\n');
else
    sample_ids   = [];
    n_samples    = 0;
    unique_samples = {};
    fprintf('  Cells: %d | Genes: %d | Clusters: %d\n', n_cells, n_genes, n_clusters);
    fprintf('  WARNING: no donor_id/sample_id column found. Falling back to Spearman MI.\n\n');
end

MAX_CELLS_PER_CLUSTER = 500;

fprintf('Computing cluster mean expression...\n');

cluster_means = zeros(n_clusters, n_genes);
cluster_std = zeros(n_clusters, n_genes);

for k = 1:n_clusters
    cells_in_cluster = (cluster_ids == unique_clusters(k));
    n_in_cluster = sum(cells_in_cluster);

    if n_in_cluster == 0
        continue;
    end

    if n_in_cluster > MAX_CELLS_PER_CLUSTER
        cell_indices = find(cells_in_cluster);
        sampled = cell_indices(randsample(n_in_cluster, MAX_CELLS_PER_CLUSTER));
        cluster_expr = full(expr_matrix(sampled, :));
    else
        cluster_expr = full(expr_matrix(cells_in_cluster, :));
    end

    if size(cluster_expr, 1) > 1
        cluster_means(k, :) = mean(cluster_expr, 1);
        cluster_std(k, :) = std(cluster_expr, 0, 1);
    else
        cluster_means(k, :) = cluster_expr;
        cluster_std(k, :) = 0.1;
    end

    if mod(k, 50) == 0
        fprintf('  Cluster %d/%d\n', k, n_clusters);
    end
end

% GP smoothing: squared-exponential kernel over expression space.
fprintf('Applying GP smoothing (sigma=0.5)...\n');

sigma_gp = 0.5;
cluster_means_smooth = zeros(size(cluster_means));

batch_size = 100;
for start_idx = 1:batch_size:n_clusters
    end_idx = min(start_idx + batch_size - 1, n_clusters);
    for k = start_idx:end_idx
        distances = sqrt(sum((cluster_means - cluster_means(k, :)).^2, 2));
        weights = exp(-distances.^2 / (2 * sigma_gp^2));
        weights = weights / sum(weights);
        cluster_means_smooth(k, :) = sum(cluster_means .* weights, 1);
    end
end

cluster_means = real(cluster_means_smooth);
fprintf('  Smoothing complete\n');

% Shannon entropy specificity: 1 = cluster-exclusive, 0 = ubiquitous.
fprintf('Computing gene expression specificity...\n');

gene_specificity = zeros(1, n_genes);
log2_n_clusters = log2(max(n_clusters, 2));

for g = 1:n_genes
    probs = cluster_means(:, g);
    total = sum(probs);
    if total < 1e-10
        gene_specificity(g) = 0;
        continue;
    end
    probs = probs / total;
    probs = probs(probs > 0);
    H = -sum(probs .* log2(probs));
    gene_specificity(g) = 1 - H / log2_n_clusters;
end

fprintf('  Mean gene specificity: %.3f\n', mean(gene_specificity));

fprintf('Loading L-R database...\n');

lr_pairs = get_cellphonedb_interactions(cellphonedb, gene_names);

if isempty(lr_pairs)
    warning('No pairs from CellPhoneDB - using curated fallback list');
    lr_pairs = get_known_ligand_receptor_pairs();
end

n_lr_pairs = size(lr_pairs, 1);
fprintf('  L-R pairs to test: %d\n', n_lr_pairs);

fprintf('Computing cluster proximity...\n');

cluster_similarity = corr(cluster_means');
cluster_similarity = max(0, cluster_similarity);

% One-step normalized graph Laplacian diffusion.
D = diag(sum(cluster_similarity, 2));
diffusion_matrix = eye(n_clusters) + 0.3 * (D \ (cluster_similarity - D));

fprintf('Testing interactions (%d pairs x %d cluster pairs)...\n', ...
    n_lr_pairs, n_clusters^2);

network.ligands           = {};
network.receptors         = {};
network.weights           = [];
network.source_clusters   = [];
network.target_clusters   = [];
network.mi_scores         = [];
network.raw_strength      = [];
network.channel_capacity  = [];   % normalized MI I(X;Y)/H(X) in [0,1]; primary score
network.ci_lower          = [];   % bootstrap 2.5th percentile of channel_capacity
network.ci_upper          = [];   % bootstrap 97.5th percentile of channel_capacity
network.mi_mode           = {};   % 'ksg_cross_sample' or 'spearman_copula'

expression_threshold = 0.1;
if ~isempty(varargin) && isnumeric(varargin{1})
    expression_threshold = varargin{1};
end

interaction_count = 0;
tested_count = 0;
progress_interval = max(1, floor(n_lr_pairs / 20));

gene_idx_map = containers.Map(gene_names, 1:length(gene_names));

% Primary score is channel capacity (KSG MI / H(ligand)), not geometric mean.
% This tests whether sender ligand variation predicts receiver receptor variation,
% rather than simply whether both genes are co-expressed.
capacity_threshold = 0.01;   % minimum normalized MI to call an interaction
MIN_SAMPLES_KSG    = 10;     % below this, Spearman fallback is used

for lr_idx = 1:n_lr_pairs
    if mod(lr_idx, progress_interval) == 0
        fprintf('  %d/%d pairs, %d interactions found\n', ...
            lr_idx, n_lr_pairs, interaction_count);
    end

    ligand = lr_pairs{lr_idx, 1};
    receptor = lr_pairs{lr_idx, 2};

    if ~isKey(gene_idx_map, ligand) || ~isKey(gene_idx_map, receptor)
        continue;
    end

    lig_idx = gene_idx_map(ligand);
    rec_idx = gene_idx_map(receptor);

    if max(cluster_means(:, lig_idx)) < expression_threshold || ...
       max(cluster_means(:, rec_idx)) < expression_threshold
        continue;
    end

    lig_spec = gene_specificity(lig_idx);
    rec_spec = gene_specificity(rec_idx);
    specificity_weight = sqrt(max(lig_spec, 0) * max(rec_spec, 0));

    for src = 1:n_clusters
        lig_expr = cluster_means(src, lig_idx);

        if lig_expr < expression_threshold
            continue;
        end

        for tgt = 1:n_clusters
            rec_expr = cluster_means(tgt, rec_idx);

            if rec_expr < expression_threshold
                continue;
            end

            tested_count = tested_count + 1;

            % Cross-sample KSG MI: X_i = mean ligand in sender cells from donor i,
            % Y_i = mean receptor in receiver cells from donor i. This is a valid
            % pairing; single-cell pairings would give MI ≈ 0 (unpaired populations).
            if n_samples >= MIN_SAMPLES_KSG
                X_samp = zeros(n_samples, 1);
                Y_samp = zeros(n_samples, 1);
                for s = 1:n_samples
                    sid = unique_samples{s};
                    src_in_s = (cluster_ids == unique_clusters(src)) & strcmp(sample_ids, sid);
                    tgt_in_s = (cluster_ids == unique_clusters(tgt)) & strcmp(sample_ids, sid);
                    if any(src_in_s)
                        X_samp(s) = mean(full(expr_matrix(src_in_s, lig_idx)));
                    end
                    if any(tgt_in_s)
                        Y_samp(s) = mean(full(expr_matrix(tgt_in_s, rec_idx)));
                    end
                end
                valid = (X_samp > 0) | (Y_samp > 0);
                X_samp = X_samp(valid);
                Y_samp = Y_samp(valid);

                if numel(X_samp) < 5
                    continue;  % too few valid donors for KSG
                end

                [capacity, mi_score, ~, ~, ci_lower, ci_upper] = estimate_channel_capacity( ...
                    X_samp, Y_samp, 'n_perm', 0, 'n_bootstrap', 500);

            else
                % Spearman copula fallback: MI ≈ -0.5·log(1-ρ²), computed on
                % cluster-mean vectors (valid pairing: both properties of the same cluster).
                X_vec = cluster_means(:, lig_idx);
                Y_vec = cluster_means(:, rec_idx);
                rho = corr(X_vec, Y_vec, 'type', 'Spearman');
                rho = max(-0.9999, min(0.9999, rho));
                mi_score = max(0, -0.5 * log(1 - rho^2) / log(2));
                H_x = max(1e-10, 0.5 * log(2*pi*exp(1)*var(X_vec) + 1e-12) / log(2));
                capacity = max(0, min(1, mi_score / H_x));
                ci_lower = NaN;
                ci_upper = NaN;
            end

            if capacity < capacity_threshold
                continue;
            end

            geometric_mean  = sqrt(lig_expr * rec_expr);
            expression_prior = geometric_mean / (geometric_mean + 0.5);   % sigmoid in [0,1]
            spatial_weight  = diffusion_matrix(src, tgt);

            strength = 0.55 * capacity ...
                     + 0.20 * expression_prior ...
                     + 0.15 * specificity_weight ...
                     + 0.10 * spatial_weight;

            interaction_count = interaction_count + 1;
            network.ligands{end+1}          = ligand;
            network.receptors{end+1}        = receptor;
            network.weights(end+1)          = strength;
            network.raw_strength(end+1)     = geometric_mean;
            network.source_clusters(end+1)  = unique_clusters(src);
            network.target_clusters(end+1)  = unique_clusters(tgt);
            network.mi_scores(end+1)        = mi_score;
            network.channel_capacity(end+1) = capacity;
            network.ci_lower(end+1)         = ci_lower;
            network.ci_upper(end+1)         = ci_upper;
            if n_samples >= MIN_SAMPLES_KSG
                network.mi_mode{end+1} = 'ksg_cross_sample';
            else
                network.mi_mode{end+1} = 'spearman_copula';
            end
        end
    end
end

network.ligands          = network.ligands(:);
network.receptors        = network.receptors(:);
network.weights          = network.weights(:);
network.raw_strength     = network.raw_strength(:);
network.source_clusters  = network.source_clusters(:);
network.target_clusters  = network.target_clusters(:);
network.mi_scores        = network.mi_scores(:);
network.channel_capacity = network.channel_capacity(:);
network.ci_lower         = network.ci_lower(:);
network.ci_upper         = network.ci_upper(:);
network.mi_mode          = network.mi_mode(:);
network.n_samples_used   = n_samples;

fprintf('  Tested: %d | Detected: %d (%.1f%%)\n', ...
    tested_count, interaction_count, 100*interaction_count/max(1,tested_count));
if n_samples >= MIN_SAMPLES_KSG
    fprintf('  MI method: KSG cross-sample (n_samples=%d)\n', n_samples);
else
    fprintf('  MI method: Spearman copula fallback — add donor_id for KSG MI\n');
end
if ~isempty(network.channel_capacity)
    fprintf('  Mean norm. MI: %.4f\n', mean(network.channel_capacity));
end

fprintf('Permutation testing (%d permutations)...\n', 50);

n_permutations = 50;
n_interactions = length(network.weights);

p_value_counters = zeros(n_interactions, 1);

gene_idx_cache = containers.Map();
for i = 1:n_interactions
    key = sprintf('%s_%s', network.ligands{i}, network.receptors{i});
    if ~isKey(gene_idx_cache, key)
        gene_idx_cache(key) = [gene_idx_map(network.ligands{i}), gene_idx_map(network.receptors{i})];
    end
end

max_perm_interactions = min(100000, n_interactions);
if n_interactions > max_perm_interactions
    perm_sample_idx = randsample(n_interactions, max_perm_interactions);
else
    perm_sample_idx = 1:n_interactions;
end

unique_genes_needed = unique([network.ligands; network.receptors]);
genes_needed_idx = zeros(length(unique_genes_needed), 1);
for g = 1:length(unique_genes_needed)
    genes_needed_idx(g) = gene_idx_map(unique_genes_needed{g});
end

if n_interactions == 0
    fprintf('  No interactions to permutation-test.\n');
    network.significance = zeros(0, 1);
    return;
end

cluster_id_to_row = containers.Map(num2cell(double(unique_clusters(:)')), ...
                                   num2cell(1:n_clusters));

for perm = 1:n_permutations
    perm_clusters = cluster_ids(randperm(length(cluster_ids)));

    perm_cluster_means = zeros(n_clusters, length(genes_needed_idx));
    for k = 1:n_clusters
        perm_cells = (perm_clusters == unique_clusters(k));
        if sum(perm_cells) > 0
            if sum(perm_cells) > MAX_CELLS_PER_CLUSTER
                cell_idx = find(perm_cells);
                cell_idx = cell_idx(randsample(length(cell_idx), MAX_CELLS_PER_CLUSTER));
                perm_cluster_means(k, :) = mean(full(expr_matrix(cell_idx, genes_needed_idx)), 1);
            else
                perm_cluster_means(k, :) = mean(full(expr_matrix(perm_cells, genes_needed_idx)), 1);
            end
        end
    end

    gene_idx_lookup = containers.Map(genes_needed_idx, 1:length(genes_needed_idx));

    for s_idx = 1:length(perm_sample_idx)
        i = perm_sample_idx(s_idx);

        key = sprintf('%s_%s', network.ligands{i}, network.receptors{i});
        indices = gene_idx_cache(key);

        lig_perm_idx = gene_idx_lookup(indices(1));
        rec_perm_idx = gene_idx_lookup(indices(2));

        src_row = cluster_id_to_row(double(network.source_clusters(i)));
        tgt_row = cluster_id_to_row(double(network.target_clusters(i)));

        perm_strength = sqrt(perm_cluster_means(src_row, lig_perm_idx) * ...
                             perm_cluster_means(tgt_row, rec_perm_idx));

        if perm_strength >= network.raw_strength(i)
            p_value_counters(i) = p_value_counters(i) + 1;
        end
    end

    if mod(perm, 10) == 0
        fprintf('  Permutation %d/%d\n', perm, n_permutations);
    end
end

network.significance = p_value_counters / n_permutations;

% Non-sampled interactions are conservatively assigned p=1 (not tested).
if length(perm_sample_idx) < n_interactions
    non_sampled = setdiff(1:n_interactions, perm_sample_idx);
    network.significance(non_sampled) = 1.0;
    fprintf('  Note: %d interactions not permutation-tested (set p=1.0)\n', ...
        length(non_sampled));
end

n_significant = sum(network.significance < 0.05);
fprintf('  Significant (p<0.05): %d/%d (%.1f%%)\n', ...
    n_significant, n_interactions, 100*n_significant/n_interactions);

fprintf('Computing topology...\n');

all_genes = unique([network.ligands; network.receptors]);
n_nodes = length(all_genes);

network.topology_metrics = struct();
network.topology_metrics.n_nodes = n_nodes;
network.topology_metrics.n_edges = n_interactions;
network.topology_metrics.density = n_interactions / max(1, n_nodes * (n_nodes - 1));

if n_interactions <= 100000
    adj_matrix = zeros(n_nodes, n_nodes);
    gene_to_idx = containers.Map(all_genes, 1:n_nodes);

    for i = 1:n_interactions
        lig_node = gene_to_idx(network.ligands{i});
        rec_node = gene_to_idx(network.receptors{i});
        adj_matrix(lig_node, rec_node) = network.weights(i);
    end

    network.topology_metrics.avg_degree = mean(sum(adj_matrix > 0, 2));
    G = digraph(adj_matrix > 0);
    network.topology_metrics.n_components = length(unique(conncomp(G)));

    out_degree = sum(adj_matrix > 0, 2);
    in_degree = sum(adj_matrix > 0, 1)';
    total_degree = out_degree + in_degree;
    hub_threshold = prctile(total_degree, 90);
    network.topology_metrics.top_hubs = all_genes(total_degree > hub_threshold);
else
    network.topology_metrics.avg_degree = 2 * n_interactions / n_nodes;
    network.topology_metrics.n_components = 1;
    network.topology_metrics.top_hubs = all_genes(1:min(10, n_nodes));
end

fprintf('\nNetwork complete: %s\n', dataset_name);
fprintf('  Interactions: %d | Significant: %d | Ligands: %d | Receptors: %d\n', ...
    n_interactions, n_significant, ...
    length(unique(network.ligands)), length(unique(network.receptors)));

end

function mi = compute_gaussian_mi(X, Y)
% Spearman copula MI: MI(X,Y) ≈ -0.5*log(1-rho^2). More robust than Pearson
% version for zero-inflated scRNA-seq counts.

n_x = length(X);
n_y = length(Y);

if n_x < 3 || n_y < 3
    mi = 0;
    return;
end

max_sample = min(100, min(n_x, n_y));

if n_x > max_sample
    X = X(randsample(n_x, max_sample));
end

if n_y > max_sample
    Y = Y(randsample(n_y, max_sample));
end

if length(X) == length(Y)
    rho = corr(X, Y, 'Type', 'Spearman');
    mi = max(0, -0.5 * log(1 - rho^2 + 1e-10));
    mi = min(1, mi);
else
    mi = 0;
end

end
