function comm_profiles = compute_cell_communication_profiles(expr_matrix, gene_names, network, cell_metadata, varargin)
% COMPUTE_CELL_COMMUNICATION_PROFILES  Per-cell sending/receiving signal
%   profiles and discovery of communication-based cell states.
%
% Scientific rationale
% --------------------
% Every existing cell-cell communication tool — CellChat, NicheNet, CellPhoneDB,
% LIANA, Tensor-cell2cell — aggregates single-cell data to CLUSTER MEANS before
% computing ligand-receptor co-expression scores. This destroys within-cluster
% heterogeneity. The critical biological reality is:
%
%   Cells within a "cluster" are NOT identical communicators.
%
% Consider a cluster of 500 GABAergic interneurons: some may be in an active
% signaling state (high ligand expression), while others are in a quiescent
% state (low ligand expression). The cluster mean averages these together into
% a single fictitious "average cell" that may not correspond to any actual cell.
%
% This function introduces COMMUNICATION PROFILES — per-cell vectors that
% capture each cell's deviation from its cluster's average signaling state:
%
%   send_score_i(L) = expr_i(L) / mean_cluster(L)   for each ligand L
%   recv_score_i(R) = expr_i(R) / mean_cluster(R)   for each receptor R
%
% A cell with send_score >> 1 is an above-average ligand sender.
% A cell with recv_score >> 1 is an above-average receptor receiver.
%
% By performing PCA and k-means clustering on these deviation profiles, we
% discover COMMUNICATION STATES that are invisible to cluster-level analysis.
% A "communication state" is a group of cells — potentially spanning multiple
% transcriptomic clusters — that share the same signaling repertoire.
%
% This represents a fundamentally new level of resolution for cell-cell
% communication analysis and is the core novel contribution of NEUROSCOPE.
%
% Computational approach
% ----------------------
% For 222K cells with ~3K HVGs, we process in batches of 1000 cells to
% avoid constructing the full pairwise similarity matrix (222K x 222K = 200GB).
% PCA on [n_cells x (n_lig+n_rec)] communication profile matrix is tractable
% since n_lig+n_rec << n_cells.
%
% Usage:
%   comm_profiles = compute_cell_communication_profiles( ...
%       expr_matrix, gene_names, network, cell_metadata)
%
%   comm_profiles = compute_cell_communication_profiles( ...
%       expr_matrix, gene_names, network, cell_metadata, ...
%       'n_pca_dims', 30, 'n_comm_clusters_factor', 3, 'batch_size', 1000)
%
% Inputs:
%   expr_matrix   - [n_cells x n_genes] log-normalized expression matrix.
%                   This should be the OUTPUT of preprocess_expression (HVG-
%                   filtered, log1p-normalized). Dense double preferred.
%   gene_names    - Gene names matching expr_matrix columns. Accepts:
%                     (a) table with 'gene_name' column
%                     (b) cell array of strings
%                     (c) struct with 'gene_name' field
%   network       - Communication network struct from build_communication_network.
%                   Must have fields: .ligands {N×1}, .receptors {N×1}.
%   cell_metadata - Table with columns: CellID, Cluster.
%
% Optional name-value parameters:
%   'n_pca_dims'               - PCA dimensions for communication manifold
%                                (default: 20).
%   'n_comm_clusters_factor'   - k for k-means = n_transcriptomic_clusters *
%                                this factor (default: 2). Allows discovery of
%                                communication states within clusters.
%   'batch_size'               - Cells processed per batch for memory efficiency
%                                (default: 1000).
%   'min_expr_threshold'       - Minimum expression (log1p units) to consider a
%                                cell as "expressing" a gene (default: 0.1).
%   'boundary_quantile'        - Quantile threshold for boundary cell detection
%                                (default: 0.15 = bottom 15% in cluster distance).
%   'verbose'                  - Print progress (default: true).
%
% Output struct fields:
%   .cell_profiles    [n_cells x (n_lig+n_rec)]  per-cell send/recv deviations
%   .ligand_names     {n_lig x 1}  ligand genes found in expr_matrix
%   .receptor_names   {n_rec x 1}  receptor genes found in expr_matrix
%   .pca_coords       [n_cells x n_pca_dims]  communication manifold PCA coords
%   .pca_variance_explained  [n_pca_dims x 1]  fraction of variance per PC
%   .comm_clusters    [n_cells x 1]  communication-based cluster assignments
%   .n_comm_clusters  scalar  number of communication states discovered
%   .boundary_cells   [n_boundary x 1]  cell indices at state transitions
%   .dominant_interactions  table  per comm-cluster top L-R interactions
%   .cluster_means_profiles [n_comm_clusters x (n_lig+n_rec)]  centroid profiles

%% --- Parse inputs --------------------------------------------------------

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

n_pca_dims          = p.Results.n_pca_dims;
cluster_factor      = p.Results.n_comm_clusters_factor;
batch_size          = p.Results.batch_size;
min_expr            = p.Results.min_expr_threshold;
boundary_q          = p.Results.boundary_quantile;
verbose             = p.Results.verbose;

%% --- Resolve gene name format -------------------------------------------

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

%% --- Extract cluster information ----------------------------------------

cluster_ids = cell_metadata.Cluster;
if iscategorical(cluster_ids)
    cluster_ids = double(cluster_ids);
elseif iscell(cluster_ids)
    [~, ~, cluster_ids] = unique(cluster_ids);
end
cluster_ids = double(cluster_ids(:));

unique_clusters  = unique(cluster_ids);
n_transcr_clust  = length(unique_clusters);

if verbose
    fprintf('  Transcriptomic clusters: %d\n', n_transcr_clust);
end

%% --- Map ligands and receptors to gene indices --------------------------

% Get unique ligands and receptors from network
all_ligands   = unique(network.ligands);
all_receptors = unique(network.receptors);

% Find which ones are present in our expression matrix
[~, lig_in_net, lig_in_mat] = intersect(all_ligands, gn, 'stable');
[~, rec_in_net, rec_in_mat] = intersect(all_receptors, gn, 'stable');

ligand_names_found   = all_ligands(lig_in_net);
receptor_names_found = all_receptors(rec_in_net);

n_lig = length(ligand_names_found);
n_rec = length(receptor_names_found);
n_lr  = n_lig + n_rec;

if verbose
    fprintf('  Ligands in network: %d | Found in expr_matrix: %d\n', ...
        length(all_ligands), n_lig);
    fprintf('  Receptors in network: %d | Found in expr_matrix: %d\n', ...
        length(all_receptors), n_rec);
    fprintf('  Communication profile dimension: %d (ligands) + %d (receptors) = %d\n', ...
        n_lig, n_rec, n_lr);
end

if n_lr == 0
    error('compute_cell_communication_profiles: no ligands or receptors found in expression matrix. Check gene name format.');
end

%% =========================================================================
%% Step 1: Compute per-cluster mean expression for L/R genes
%% =========================================================================
%
% We need the cluster-mean expression of each ligand and receptor so we can
% compute per-cell DEVIATIONS from the cluster average. This is what makes
% communication profiles novel — they capture within-cluster heterogeneity.

if verbose
    fprintf('\n--- Step 1: Computing cluster mean expression ---\n');
end

% Extract L/R columns from expression matrix
% For large matrices, do this in pieces
lig_cols = lig_in_mat;   % indices into expr_matrix columns
rec_cols = rec_in_mat;

% Compute cluster means for ligands
cluster_mean_lig = zeros(n_transcr_clust, n_lig);  % [k x n_lig]
cluster_mean_rec = zeros(n_transcr_clust, n_rec);  % [k x n_rec]

for k = 1:n_transcr_clust
    clust_k = unique_clusters(k);
    cell_mask = (cluster_ids == clust_k);
    n_in_cluster = sum(cell_mask);

    if n_in_cluster == 0
        continue;
    end

    % Extract ligand columns for this cluster
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

% Build per-cell cluster mean arrays (needed for deviation computation)
% cluster_id_idx maps cell -> row in cluster_mean arrays
[~, cluster_id_idx] = ismember(cluster_ids, unique_clusters);

% Per-cell cluster means: [n_cells x n_lig] and [n_cells x n_rec]
cell_cluster_mean_lig = cluster_mean_lig(cluster_id_idx, :);  % [n_cells x n_lig]
cell_cluster_mean_rec = cluster_mean_rec(cluster_id_idx, :);  % [n_cells x n_rec]

if verbose
    fprintf('  Cluster mean matrices computed: %d clusters x (%d lig + %d rec)\n', ...
        n_transcr_clust, n_lig, n_rec);
end

%% =========================================================================
%% Step 2: Compute per-cell communication profiles in batches
%% =========================================================================
%
% For each cell i:
%   send_score_i(L) = expr_i(L) / max(cluster_mean(L), epsilon)
%   recv_score_i(R) = expr_i(R) / max(cluster_mean(R), epsilon)
%
% Division by epsilon prevents extreme ratios when cluster mean is near 0.
% Scores are capped at 10x to prevent outlier cells from dominating PCA.
%
% The concatenated [send_scores; recv_scores] forms the communication
% profile: a fingerprint of how this cell's signaling state deviates from
% its transcriptomic cluster's average.

if verbose
    fprintf('\n--- Step 2: Computing per-cell communication profiles ---\n');
    fprintf('  Processing %d cells in batches of %d...\n', n_cells, batch_size);
end

t_profiles = tic;
epsilon = 1e-6;   % Prevents division by zero; also handles near-zero cluster means
max_ratio = 10.0; % Cap ratio at 10x to prevent outlier dominance

cell_profiles = zeros(n_cells, n_lr, 'single');  % [n_cells x n_lr]

n_batches = ceil(n_cells / batch_size);
for b = 1:n_batches
    idx_start = (b-1)*batch_size + 1;
    idx_end   = min(b*batch_size, n_cells);
    batch_idx = idx_start:idx_end;
    n_batch   = length(batch_idx);

    % Extract expression for this batch
    if issparse(expr_matrix)
        lig_expr_batch = full(expr_matrix(batch_idx, lig_cols));
        rec_expr_batch = full(expr_matrix(batch_idx, rec_cols));
    else
        lig_expr_batch = expr_matrix(batch_idx, lig_cols);
        rec_expr_batch = expr_matrix(batch_idx, rec_cols);
    end
    lig_expr_batch = double(lig_expr_batch);
    rec_expr_batch = double(rec_expr_batch);

    % Get corresponding cluster means for this batch
    cm_lig_batch = cell_cluster_mean_lig(batch_idx, :);  % [n_batch x n_lig]
    cm_rec_batch = cell_cluster_mean_rec(batch_idx, :);  % [n_batch x n_rec]

    % Compute deviation ratios
    send_scores = lig_expr_batch ./ max(cm_lig_batch, epsilon);
    recv_scores = rec_expr_batch ./ max(cm_rec_batch, epsilon);

    % Cap extreme ratios
    send_scores = min(send_scores, max_ratio);
    recv_scores = min(recv_scores, max_ratio);

    % Zero out scores for genes below expression threshold
    % (cells not expressing a gene at all should have score 0, not ~1)
    send_scores(lig_expr_batch < min_expr) = 0;
    recv_scores(rec_expr_batch < min_expr) = 0;

    % Concatenate into profile vector
    cell_profiles(batch_idx, :) = single([send_scores, recv_scores]);

    if verbose && mod(b, 50) == 0
        fprintf('  Batch %d/%d (cells %d-%d)\n', b, n_batches, idx_start, idx_end);
    end
end

if verbose
    fprintf('  Profiles computed in %.1f s\n', toc(t_profiles));
    % Profile statistics
    mean_send = mean(cell_profiles(:, 1:n_lig), 'all');
    mean_recv = mean(cell_profiles(:, n_lig+1:end), 'all');
    fprintf('  Mean sending score: %.3f (expected ~1.0 for cluster-average cells)\n', mean_send);
    fprintf('  Mean receiving score: %.3f\n', mean_recv);

    % Fraction of cells with above-average signaling
    above_avg_send = mean(max(cell_profiles(:, 1:n_lig), [], 2) > 1.5);
    fprintf('  Cells with strong sending (max send>1.5): %.1f%%\n', above_avg_send*100);
end

%% =========================================================================
%% Step 3: PCA on communication profiles (the "communication manifold")
%% =========================================================================
%
% PCA on the [n_cells x n_lr] communication profile matrix.
% n_lr is typically 200-600 dimensions (unique ligands + receptors).
% PCA is tractable directly on this matrix.
%
% The top PCs capture the dominant axes of variation in how cells deviate
% from their cluster's average signaling. PC1 might separate "sender-dominant"
% from "receiver-dominant" cells. PC2 might separate cells using different
% ligand families, etc.

if verbose
    fprintf('\n--- Step 3: PCA on communication profiles ---\n');
end

t_pca = tic;

% Center the profiles before PCA (zero-mean per feature)
profile_double = double(cell_profiles);
profile_mean   = mean(profile_double, 1);
profile_centered = bsxfun(@minus, profile_double, profile_mean);

% PCA via SVD (economy mode — only compute top n_pca_dims components)
% For n_cells >> n_lr, we can use the covariance matrix approach for speed
n_pca_actual = min(n_pca_dims, min(n_cells, n_lr) - 1);

if n_cells > n_lr
    % Standard PCA via covariance matrix (faster when n_cells >> n_lr)
    % C = X'X / (n-1) is [n_lr x n_lr], feasible
    cov_matrix = (profile_centered' * profile_centered) / (n_cells - 1);
    [V, D] = eigs(cov_matrix, n_pca_actual, 'largestabs');  % top eigenvectors
    eigenvalues = diag(D);
    pca_coords = profile_centered * V;  % [n_cells x n_pca_actual]
else
    % Gram matrix approach when n_cells < n_lr (shouldn't happen in practice)
    [~, S, V] = svds(profile_centered, n_pca_actual);
    eigenvalues = diag(S).^2 / (n_cells - 1);
    pca_coords = profile_centered * V;
end

% Variance explained per PC
total_var = sum(var(profile_double, 0, 1));  % total variance in profiles
var_explained = eigenvalues / total_var;

if verbose
    fprintf('  PCA dimensions: %d\n', n_pca_actual);
    fprintf('  Variance explained:\n');
    cumvar = cumsum(var_explained);
    for pc = [1, 2, 3, 5, 10, min(20, n_pca_actual)]
        if pc <= n_pca_actual
            fprintf('    PC1-%2d: %.1f%%\n', pc, cumvar(pc)*100);
        end
    end
end

fprintf('  PCA complete (%.1f s)\n', toc(t_pca));

%% =========================================================================
%% Step 4: K-means clustering on communication manifold
%% =========================================================================
%
% We cluster cells on their PCA communication coordinates to find
% "communication states" — groups of cells that share the same signaling
% identity, regardless of their transcriptomic cluster assignment.
%
% k = n_transcr_clusters * cluster_factor (default 2x)
% This lets us find at least 2 communication sub-states per transcriptomic
% cluster. The actual number could be fewer if clusters have uniform signaling.
%
% After k-means, we identify "boundary cells" — cells that are close to
% the centroid of a DIFFERENT communication state than their assigned one.
% These are likely cells mid-way through a communication state transition.

if verbose
    fprintf('\n--- Step 4: Communication state discovery (k-means) ---\n');
end

t_kmeans = tic;

n_comm_clusters = n_transcr_clust * cluster_factor;
n_comm_clusters = min(n_comm_clusters, floor(n_cells / 10));  % at least 10 cells per cluster
n_comm_clusters = max(n_comm_clusters, 2);                    % at least 2 states

if verbose
    fprintf('  k-means with k=%d communication states...\n', n_comm_clusters);
    fprintf('  (Rationale: %d transcriptomic clusters x %d factor)\n', ...
        n_transcr_clust, cluster_factor);
end

% Use k-means with multiple restarts for stability
try
    [comm_clusters, centroids, sumd] = kmeans(pca_coords, n_comm_clusters, ...
        'Replicates', 5, 'MaxIter', 300, 'Display', 'off', 'Distance', 'sqeuclidean');
catch ME
    fprintf('  K-means warning: %s. Using single replicate.\n', ME.message);
    [comm_clusters, centroids, sumd] = kmeans(pca_coords, n_comm_clusters, ...
        'Replicates', 1, 'MaxIter', 200, 'Display', 'off');
end

% Cluster sizes
cluster_sizes_comm = histcounts(comm_clusters, 1:n_comm_clusters+1);

if verbose
    fprintf('  K-means complete (%.1f s)\n', toc(t_kmeans));
    fprintf('  Communication cluster sizes (min/median/max): %d / %.0f / %d\n', ...
        min(cluster_sizes_comm), median(cluster_sizes_comm), max(cluster_sizes_comm));
end

%% =========================================================================
%% Step 5: Identify boundary cells (state transition cells)
%% =========================================================================
%
% Boundary cells are cells that:
%   - Are assigned to communication cluster C_i
%   - BUT their distance to another cluster C_j's centroid is close to their
%     distance to C_i's centroid (ratio > boundary_q threshold)
%
% These cells likely represent transitional states between two communication
% programs. They are scientifically interesting because they may reveal the
% mechanisms of communication state switching.

if verbose
    fprintf('\n--- Step 5: Boundary cell detection ---\n');
end

t_bound = tic;

% For each cell, compute distance to assigned centroid and to nearest OTHER centroid
dist_to_assigned  = zeros(n_cells, 1);
dist_to_nearest_other = zeros(n_cells, 1);

% Process in batches to avoid huge distance matrix
for b = 1:n_batches
    idx_start = (b-1)*batch_size + 1;
    idx_end   = min(b*batch_size, n_cells);
    batch_idx = idx_start:idx_end;

    pca_batch = pca_coords(batch_idx, :);  % [n_batch x n_pca]
    labels_batch = comm_clusters(batch_idx);

    % Distance to each centroid: [n_batch x k]
    dist_to_all = pdist2(pca_batch, centroids);

    % Distance to assigned centroid
    for j = 1:length(batch_idx)
        dist_to_assigned(batch_idx(j)) = dist_to_all(j, labels_batch(j));
        % Distance to nearest other centroid
        other_dists = dist_to_all(j, :);
        other_dists(labels_batch(j)) = Inf;
        dist_to_nearest_other(batch_idx(j)) = min(other_dists);
    end
end

% Boundary score: how close is the cell to another cluster relative to its own?
% Score near 1.0 = equidistant (boundary), score near 0 = well inside cluster
boundary_score = dist_to_assigned ./ (dist_to_assigned + dist_to_nearest_other);
boundary_score(~isfinite(boundary_score)) = 0;

% Boundary cells: those in top boundary_q fraction of boundary_score
boundary_threshold = quantile(boundary_score, 1 - boundary_q);
boundary_cells = find(boundary_score >= boundary_threshold);

if verbose
    fprintf('  Boundary cells identified: %d (%.1f%% of total)\n', ...
        length(boundary_cells), 100*length(boundary_cells)/n_cells);
    fprintf('  Boundary score threshold: %.3f\n', boundary_threshold);
    fprintf('  Boundary detection complete (%.1f s)\n', toc(t_bound));
end

%% =========================================================================
%% Step 6: Dominant interactions per communication state
%% =========================================================================
%
% For each communication-based cluster (state), identify the L-R pairs that
% most strongly characterize that state. We compute the mean communication
% profile for each state and identify the top-scoring ligand-receptor features.

if verbose
    fprintf('\n--- Step 6: Dominant interactions per communication state ---\n');
end

% Compute cluster-mean communication profiles
cluster_means_profiles = zeros(n_comm_clusters, n_lr);
for k = 1:n_comm_clusters
    mask_k = (comm_clusters == k);
    if sum(mask_k) > 0
        cluster_means_profiles(k, :) = mean(double(cell_profiles(mask_k, :)), 1);
    end
end

% Build a table of dominant interactions per cluster
all_feature_names = [ligand_names_found; receptor_names_found];
feature_type = [repmat({'ligand'}, n_lig, 1); repmat({'receptor'}, n_rec, 1)];

dominant_rows = {};
TOP_N = 10;  % top 10 features per communication state

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

% Convert to table
if ~isempty(dominant_rows)
    dom_matrix = vertcat(dominant_rows{:});
    dominant_interactions = table( ...
        cell2mat(dom_matrix(:,1)), ...
        dom_matrix(:,2), ...
        dom_matrix(:,3), ...
        cell2mat(dom_matrix(:,4)), ...
        cell2mat(dom_matrix(:,5)), ...
        'VariableNames', {'CommCluster', 'Gene', 'Type', 'MeanScore', 'Rank'});
else
    dominant_interactions = table();
end

if verbose
    fprintf('  Dominant interaction table: %d rows\n', height(dominant_interactions));
    % Show top interaction for each of first 5 comm clusters
    fprintf('  Top ligand/receptor per communication state (first 5):\n');
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

%% =========================================================================
%% Assemble output struct
%% =========================================================================

comm_profiles = struct();
comm_profiles.cell_profiles           = cell_profiles;         % [n_cells x n_lr] single
comm_profiles.ligand_names            = ligand_names_found;    % {n_lig x 1}
comm_profiles.receptor_names          = receptor_names_found;  % {n_rec x 1}
comm_profiles.pca_coords              = pca_coords;            % [n_cells x n_pca]
comm_profiles.pca_variance_explained  = var_explained;         % [n_pca x 1]
comm_profiles.comm_clusters           = comm_clusters;         % [n_cells x 1]
comm_profiles.n_comm_clusters         = n_comm_clusters;       % scalar
comm_profiles.boundary_cells          = boundary_cells;        % [n_boundary x 1]
comm_profiles.boundary_score          = boundary_score;        % [n_cells x 1]
comm_profiles.dominant_interactions   = dominant_interactions; % table
comm_profiles.cluster_means_profiles  = cluster_means_profiles; % [k x n_lr]
comm_profiles.profile_feature_names   = all_feature_names;    % {n_lr x 1}
comm_profiles.n_lig                   = n_lig;
comm_profiles.n_rec                   = n_rec;

% Save progress CSVs
if exist('results', 'dir')
    try
        % Communication cluster assignments
        comm_cluster_table = table((1:n_cells)', comm_clusters, boundary_score, ...
            'VariableNames', {'cell_index', 'comm_cluster', 'boundary_score'});
        writetable(comm_cluster_table, fullfile('results', 'communication_clusters.csv'));

        % Dominant interactions
        if height(dominant_interactions) > 0
            writetable(dominant_interactions, fullfile('results', 'dominant_interactions.csv'));
        end
        if verbose
            fprintf('  Progress CSVs saved to results/\n');
        end
    catch
        % Non-fatal
    end
end

if verbose
    fprintf('\n=== compute_cell_communication_profiles summary ===\n');
    fprintf('  Cells:                  %d\n', n_cells);
    fprintf('  Communication features: %d (%d lig + %d rec)\n', n_lr, n_lig, n_rec);
    fprintf('  PCA dimensions:         %d (%.1f%% variance)\n', ...
        n_pca_actual, cumsum(var_explained(end))*100);
    fprintf('  Communication states:   %d (vs. %d transcriptomic clusters)\n', ...
        n_comm_clusters, n_transcr_clust);
    fprintf('  Boundary cells:         %d (%.1f%%)\n', ...
        length(boundary_cells), 100*length(boundary_cells)/n_cells);
    fprintf('compute_cell_communication_profiles done.\n\n');
end

end
