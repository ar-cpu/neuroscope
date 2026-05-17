function comm_profiles = compute_cell_communication_profiles(expr_matrix, gene_names, network, cell_metadata, varargin)
% COMPUTE_CELL_COMMUNICATION_PROFILES  Per-cell sending/receiving signal
%   profiles and communication-based cell state discovery.
%
% Usage:
%   comm_profiles = compute_cell_communication_profiles( ...
%       expr_matrix, gene_names, network, cell_metadata)
%   comm_profiles = compute_cell_communication_profiles( ...
%       expr_matrix, gene_names, network, cell_metadata, ...
%       'n_pca_dims', 30, 'n_comm_clusters_factor', 3, 'batch_size', 1000)
%
% Inputs:
%   expr_matrix   - [n_cells x n_genes] log-normalized expression matrix (dense double preferred)
%   gene_names    - gene names matching columns: table with 'gene_name', cell array, or struct
%   network       - struct from build_communication_network with fields .ligands, .receptors
%   cell_metadata - table with columns: CellID, Cluster
%
% Optional parameters:
%   'n_pca_dims'             - PCA dimensions (default 20)
%   'n_comm_clusters_factor' - k = n_transcriptomic_clusters * factor (default 2)
%   'batch_size'             - cells per batch (default 1000)
%   'min_expr_threshold'     - min log1p expression to count as expressed (default 0.1)
%   'boundary_quantile'      - top fraction flagged as boundary cells (default 0.15)
%   'verbose'                - print progress (default true)
%
% Output struct fields:
%   .cell_profiles    [n_cells x (n_lig+n_rec)]  per-cell send/recv deviations
%   .ligand_names / .receptor_names  genes found in expr_matrix
%   .pca_coords       [n_cells x n_pca_dims]
%   .pca_variance_explained  [n_pca_dims x 1]
%   .comm_clusters    [n_cells x 1]  communication-state assignments
%   .boundary_cells   [n_boundary x 1]  cells near state boundaries
%   .dominant_interactions  table  top L-R features per comm state
%   .cluster_means_profiles [n_comm_clusters x (n_lig+n_rec)]

p = inputParser();
p.addRequired('expr_matrix');
p.addRequired('gene_names');
p.addRequired('network');
p.addRequired('cell_metadata');
p.addParameter('n_pca_dims',             20,    @(x) isnumeric(x) && x > 0);
p.addParameter('n_comm_clusters_factor',  2,    @(x) isnumeric(x) && x > 0);
p.addParameter('batch_size',           1000,    @(x) isnumeric(x) && x > 0);
p.addParameter('min_expr_threshold',    0.1,    @(x) isnumeric(x) && x >= 0);
p.addParameter('boundary_quantile',    0.15,    @(x) isnumeric(x) && x > 0 && x < 1);
p.addParameter('verbose',             true,    @islogical);
p.parse(expr_matrix, gene_names, network, cell_metadata, varargin{:});

n_pca_dims     = p.Results.n_pca_dims;
cluster_factor = p.Results.n_comm_clusters_factor;
batch_size     = p.Results.batch_size;
min_expr       = p.Results.min_expr_threshold;
boundary_q     = p.Results.boundary_quantile;
verbose        = p.Results.verbose;

% Resolve gene name format
if istable(gene_names)
    if ismember('gene_name', gene_names.Properties.VariableNames)
        gn = gene_names.gene_name;
    else
        gn = gene_names{:,1};
    end
elseif iscell(gene_names)
    gn = gene_names;
elseif isstruct(gene_names) && isfield(gene_names, 'gene_name')
    gn = gene_names.gene_name;
else
    error('compute_cell_communication_profiles: unrecognized gene_names format.');
end
if ~iscell(gn)
    gn = cellstr(num2str(gn));
end
gn = gn(:);

[n_cells, n_genes] = size(expr_matrix);

if verbose
    fprintf('compute_cell_communication_profiles:\n');
    fprintf('  Input: %d cells x %d genes\n', n_cells, n_genes);
end

% Extract cluster IDs
cluster_ids = cell_metadata.Cluster;
if iscategorical(cluster_ids)
    cluster_ids = double(cluster_ids);
elseif iscell(cluster_ids)
    [~, ~, cluster_ids] = unique(cluster_ids);
end
cluster_ids = double(cluster_ids(:));

unique_clusters = unique(cluster_ids);
n_transcr_clust = length(unique_clusters);

if verbose
    fprintf('  Transcriptomic clusters: %d\n', n_transcr_clust);
end

% Map ligands and receptors to gene indices
all_ligands   = unique(network.ligands);
all_receptors = unique(network.receptors);

[~, lig_in_net, lig_in_mat] = intersect(all_ligands,   gn, 'stable');
[~, rec_in_net, rec_in_mat] = intersect(all_receptors, gn, 'stable');

ligand_names_found   = all_ligands(lig_in_net);
receptor_names_found = all_receptors(rec_in_net);

n_lig = length(ligand_names_found);
n_rec = length(receptor_names_found);
n_lr  = n_lig + n_rec;

if verbose
    fprintf('  Ligands: %d in network, %d in matrix\n', length(all_ligands), n_lig);
    fprintf('  Receptors: %d in network, %d in matrix\n', length(all_receptors), n_rec);
    fprintf('  Profile dimension: %d\n', n_lr);
end

if n_lr == 0
    error('compute_cell_communication_profiles: no ligands or receptors found. Check gene name format.');
end

% Step 1: Cluster mean expression for L/R genes
if verbose, fprintf('\n--- Step 1: Cluster mean expression ---\n'); end

lig_cols = lig_in_mat;
rec_cols = rec_in_mat;

cluster_mean_lig = zeros(n_transcr_clust, n_lig);
cluster_mean_rec = zeros(n_transcr_clust, n_rec);

for k = 1:n_transcr_clust
    clust_k   = unique_clusters(k);
    cell_mask = (cluster_ids == clust_k);
    if sum(cell_mask) == 0, continue; end

    if issparse(expr_matrix)
        lig_block = full(expr_matrix(cell_mask, lig_cols));
        rec_block = full(expr_matrix(cell_mask, rec_cols));
    else
        lig_block = expr_matrix(cell_mask, lig_cols);
        rec_block = expr_matrix(cell_mask, rec_cols);
    end

    cluster_mean_lig(k, :) = mean(lig_block, 1);
    cluster_mean_rec(k, :) = mean(rec_block, 1);
end

[~, cluster_id_idx] = ismember(cluster_ids, unique_clusters);
cell_cluster_mean_lig = cluster_mean_lig(cluster_id_idx, :);
cell_cluster_mean_rec = cluster_mean_rec(cluster_id_idx, :);

% Step 2: Per-cell communication profiles (batched)
% send_score = expr / cluster_mean; zero below min_expr; capped at 10x
if verbose
    fprintf('\n--- Step 2: Per-cell communication profiles ---\n');
    fprintf('  Processing %d cells in batches of %d...\n', n_cells, batch_size);
end

t_profiles = tic;
epsilon   = 1e-6;
max_ratio = 10.0;

cell_profiles = zeros(n_cells, n_lr, 'single');
n_batches = ceil(n_cells / batch_size);

for b = 1:n_batches
    idx_start = (b-1)*batch_size + 1;
    idx_end   = min(b*batch_size, n_cells);
    batch_idx = idx_start:idx_end;

    if issparse(expr_matrix)
        lig_expr_batch = full(expr_matrix(batch_idx, lig_cols));
        rec_expr_batch = full(expr_matrix(batch_idx, rec_cols));
    else
        lig_expr_batch = expr_matrix(batch_idx, lig_cols);
        rec_expr_batch = expr_matrix(batch_idx, rec_cols);
    end
    lig_expr_batch = double(lig_expr_batch);
    rec_expr_batch = double(rec_expr_batch);

    cm_lig_batch = cell_cluster_mean_lig(batch_idx, :);
    cm_rec_batch = cell_cluster_mean_rec(batch_idx, :);

    send_scores = min(lig_expr_batch ./ max(cm_lig_batch, epsilon), max_ratio);
    recv_scores = min(rec_expr_batch ./ max(cm_rec_batch, epsilon), max_ratio);

    % Cells not expressing a gene get score 0, not ~1
    send_scores(lig_expr_batch < min_expr) = 0;
    recv_scores(rec_expr_batch < min_expr) = 0;

    cell_profiles(batch_idx, :) = single([send_scores, recv_scores]);

    if verbose && mod(b, 50) == 0
        fprintf('  Batch %d/%d (cells %d-%d)\n', b, n_batches, idx_start, idx_end);
    end
end

if verbose
    fprintf('  Profiles computed in %.1f s\n', toc(t_profiles));
    mean_send = mean(cell_profiles(:, 1:n_lig), 'all');
    mean_recv = mean(cell_profiles(:, n_lig+1:end), 'all');
    fprintf('  Mean sending score: %.3f | Mean receiving score: %.3f\n', mean_send, mean_recv);
    above_avg_send = mean(max(cell_profiles(:, 1:n_lig), [], 2) > 1.5);
    fprintf('  Cells with strong sending (max>1.5): %.1f%%\n', above_avg_send*100);
end

% Step 3: PCA on communication profiles
if verbose, fprintf('\n--- Step 3: PCA on communication profiles ---\n'); end

t_pca = tic;

profile_double   = double(cell_profiles);
profile_mean     = mean(profile_double, 1);
profile_centered = bsxfun(@minus, profile_double, profile_mean);

n_pca_actual = min(n_pca_dims, min(n_cells, n_lr) - 1);

if n_cells > n_lr
    cov_matrix = (profile_centered' * profile_centered) / (n_cells - 1);
    [V, D] = eigs(cov_matrix, n_pca_actual, 'largestabs');
    eigenvalues = diag(D);
    pca_coords  = profile_centered * V;
else
    [~, S, V]   = svds(profile_centered, n_pca_actual);
    eigenvalues = diag(S).^2 / (n_cells - 1);
    pca_coords  = profile_centered * V;
end

total_var   = sum(var(profile_double, 0, 1));
var_explained = eigenvalues / total_var;

if verbose
    fprintf('  PCA dimensions: %d\n', n_pca_actual);
    cumvar = cumsum(var_explained);
    for pc = [1, 2, 3, 5, 10, min(20, n_pca_actual)]
        if pc <= n_pca_actual
            fprintf('    PC1-%2d: %.1f%%\n', pc, cumvar(pc)*100);
        end
    end
end
fprintf('  PCA complete (%.1f s)\n', toc(t_pca));

% Step 4: K-means on communication manifold
% k = n_transcriptomic_clusters * cluster_factor; bounded by cell count
if verbose, fprintf('\n--- Step 4: Communication state discovery ---\n'); end

t_kmeans = tic;

n_comm_clusters = n_transcr_clust * cluster_factor;
n_comm_clusters = min(n_comm_clusters, floor(n_cells / 10));
n_comm_clusters = max(n_comm_clusters, 2);

if verbose
    fprintf('  k-means k=%d (%d clusters x %d factor)...\n', ...
        n_comm_clusters, n_transcr_clust, cluster_factor);
end

try
    [comm_clusters, centroids, ~] = kmeans(pca_coords, n_comm_clusters, ...
        'Replicates', 5, 'MaxIter', 300, 'Display', 'off', 'Distance', 'sqeuclidean');
catch ME
    fprintf('  K-means warning: %s. Using single replicate.\n', ME.message);
    [comm_clusters, centroids, ~] = kmeans(pca_coords, n_comm_clusters, ...
        'Replicates', 1, 'MaxIter', 200, 'Display', 'off');
end

cluster_sizes_comm = histcounts(comm_clusters, 1:n_comm_clusters+1);

if verbose
    fprintf('  K-means complete (%.1f s)\n', toc(t_kmeans));
    fprintf('  Cluster sizes min/median/max: %d / %.0f / %d\n', ...
        min(cluster_sizes_comm), median(cluster_sizes_comm), max(cluster_sizes_comm));
end

% Step 5: Boundary cells (near state transition)
% Boundary score = dist_to_assigned / (dist_to_assigned + dist_to_nearest_other)
if verbose, fprintf('\n--- Step 5: Boundary cell detection ---\n'); end

t_bound = tic;

dist_to_assigned      = zeros(n_cells, 1);
dist_to_nearest_other = zeros(n_cells, 1);

for b = 1:n_batches
    idx_start  = (b-1)*batch_size + 1;
    idx_end    = min(b*batch_size, n_cells);
    batch_idx  = idx_start:idx_end;

    pca_batch    = pca_coords(batch_idx, :);
    labels_batch = comm_clusters(batch_idx);
    dist_to_all  = pdist2(pca_batch, centroids);

    for j = 1:length(batch_idx)
        dist_to_assigned(batch_idx(j)) = dist_to_all(j, labels_batch(j));
        other_dists = dist_to_all(j, :);
        other_dists(labels_batch(j)) = Inf;
        dist_to_nearest_other(batch_idx(j)) = min(other_dists);
    end
end

boundary_score = dist_to_assigned ./ (dist_to_assigned + dist_to_nearest_other);
boundary_score(~isfinite(boundary_score)) = 0;

boundary_threshold = quantile(boundary_score, 1 - boundary_q);
boundary_cells = find(boundary_score >= boundary_threshold);

if verbose
    fprintf('  Boundary cells: %d (%.1f%%), threshold=%.3f (%.1f s)\n', ...
        length(boundary_cells), 100*length(boundary_cells)/n_cells, ...
        boundary_threshold, toc(t_bound));
end

% Step 6: Dominant interactions per communication state
if verbose, fprintf('\n--- Step 6: Dominant interactions ---\n'); end

cluster_means_profiles = zeros(n_comm_clusters, n_lr);
for k = 1:n_comm_clusters
    mask_k = (comm_clusters == k);
    if sum(mask_k) > 0
        cluster_means_profiles(k, :) = mean(double(cell_profiles(mask_k, :)), 1);
    end
end

all_feature_names = [ligand_names_found; receptor_names_found];
feature_type = [repmat({'ligand'}, n_lig, 1); repmat({'receptor'}, n_rec, 1)];

dominant_rows = {};
TOP_N = 10;

for k = 1:n_comm_clusters
    profile_k = cluster_means_profiles(k, :);
    [~, top_idx] = sort(profile_k, 'descend');
    top_idx = top_idx(1:min(TOP_N, length(top_idx)));
    for j = 1:length(top_idx)
        feat_idx = top_idx(j);
        dominant_rows{end+1} = {k, all_feature_names{feat_idx}, ...
            feature_type{feat_idx}, profile_k(feat_idx), j};
    end
end

if ~isempty(dominant_rows)
    dom_matrix = vertcat(dominant_rows{:});
    dominant_interactions = table( ...
        cell2mat(dom_matrix(:,1)), dom_matrix(:,2), dom_matrix(:,3), ...
        cell2mat(dom_matrix(:,4)), cell2mat(dom_matrix(:,5)), ...
        'VariableNames', {'CommCluster','Gene','Type','MeanScore','Rank'});
else
    dominant_interactions = table();
end

if verbose
    fprintf('  Dominant interaction table: %d rows\n', height(dominant_interactions));
    for k = 1:min(5, n_comm_clusters)
        mask_k = (dominant_interactions.CommCluster == k) & ...
                 strcmp(dominant_interactions.Type, 'ligand') & ...
                 (dominant_interactions.Rank == 1);
        if any(mask_k)
            row = find(mask_k, 1);
            fprintf('    State %2d: top ligand = %-12s (score=%.3f)\n', ...
                k, dominant_interactions.Gene{row}, dominant_interactions.MeanScore(row));
        end
    end
end

% Assemble output
comm_profiles = struct();
comm_profiles.cell_profiles           = cell_profiles;
comm_profiles.ligand_names            = ligand_names_found;
comm_profiles.receptor_names          = receptor_names_found;
comm_profiles.pca_coords              = pca_coords;
comm_profiles.pca_variance_explained  = var_explained;
comm_profiles.comm_clusters           = comm_clusters;
comm_profiles.n_comm_clusters         = n_comm_clusters;
comm_profiles.boundary_cells          = boundary_cells;
comm_profiles.boundary_score          = boundary_score;
comm_profiles.dominant_interactions   = dominant_interactions;
comm_profiles.cluster_means_profiles  = cluster_means_profiles;
comm_profiles.profile_feature_names   = all_feature_names;
comm_profiles.n_lig                   = n_lig;
comm_profiles.n_rec                   = n_rec;

if exist('results', 'dir')
    try
        comm_cluster_table = table((1:n_cells)', comm_clusters, boundary_score, ...
            'VariableNames', {'cell_index','comm_cluster','boundary_score'});
        writetable(comm_cluster_table, fullfile('results', 'communication_clusters.csv'));
        if height(dominant_interactions) > 0
            writetable(dominant_interactions, fullfile('results', 'dominant_interactions.csv'));
        end
        if verbose, fprintf('  Progress CSVs saved to results/\n'); end
    catch
    end
end

if verbose
    fprintf('\n=== compute_cell_communication_profiles summary ===\n');
    fprintf('  Cells: %d | Features: %d (%d lig + %d rec)\n', n_cells, n_lr, n_lig, n_rec);
    fprintf('  PCA dims: %d | Comm states: %d | Boundary cells: %d (%.1f%%)\n', ...
        n_pca_actual, n_comm_clusters, length(boundary_cells), 100*length(boundary_cells)/n_cells);
    fprintf('compute_cell_communication_profiles done.\n\n');
end

end
