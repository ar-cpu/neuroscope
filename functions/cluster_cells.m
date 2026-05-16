function clusters = cluster_cells(expr_matrix, metadata, dataset_name)
% CLUSTER_CELLS  Partition cells into transcriptionally coherent communities
%                using a graph-based Louvain pipeline that mirrors the approach
%                implemented in Seurat (Butler et al. 2018) and Scanpy
%                (Wolf et al. 2018).
%
% Pipeline overview
% -----------------
%   Step 1 — Intelligent subsampling: datasets exceeding 15,000 cells are
%             randomly subsampled for graph construction. Remaining cells are
%             assigned to their nearest centroid in PCA space.
%   Step 2 — Highly variable gene (HVG) selection: genes are ranked by the
%             variance-to-mean dispersion ratio. The top 30% (capped at 1,000)
%             are retained to focus the embedding on informative features.
%   Step 3 — PCA: HVG expression is centered, unit-scaled, and decomposed.
%             The number of principal components is chosen to explain ≥80% of
%             variance, bounded to [15, 30].
%   Step 4 — Shared nearest-neighbor (SNN) graph: a k-NN graph (k=30) is
%             built in PCA space (chunked Euclidean distance). Edge weights
%             are Jaccard similarities of overlapping neighborhoods, yielding a
%             topology that is more robust to local density variation than raw
%             distance.
%   Step 5 — Louvain community detection: the SNN graph is thresholded at its
%             50th-percentile edge weight to reduce noise, then subjected to
%             greedy modularity maximization. The number of communities is
%             determined automatically.
%   Step 6 — Silhouette validation: average silhouette score over the sampled
%             cells is reported as a data-driven quality metric.
%
% Usage:
%   clusters = cluster_cells(expr_matrix, metadata, dataset_name)
%
% Inputs:
%   expr_matrix  - [cells × genes] sparse or dense normalized expression matrix.
%   metadata     - Table of cell-level metadata (unused internally; retained for
%                  API compatibility with the NEUROSCOPE pipeline).
%   dataset_name - String label used in progress messages.
%
% Outputs:
%   clusters - Table with columns Cell_ID (1-based integer) and Cluster
%              (1-based integer community label). All cells receive an
%              assignment.
%
% References:
%   Butler et al. (2018) Integrating single-cell transcriptomic data across
%     different conditions, technologies, and species. Nature Biotechnology 36:411-420.
%   Wolf et al. (2018) SCANPY: large-scale single-cell gene expression data
%     analysis. Genome Biology 19:15.
%   Blondel et al. (2008) Fast unfolding of communities in large networks.
%     Journal of Statistical Mechanics: Theory and Experiment 2008:P10008.

fprintf('      Graph-based clustering (Louvain + SNN)\n\n');

n_cells = size(expr_matrix, 1);
n_genes = size(expr_matrix, 2);

fprintf('      Dataset: %s (%d cells x %d genes)\n', dataset_name, n_cells, n_genes);

%% Step 1: Subsampling for memory-efficient graph construction
max_sample = 15000;
if n_cells > max_sample
    fprintf('      Subsampling to %d cells for graph construction...\n', max_sample);
    sample_idx = randsample(n_cells, max_sample);
    expr_sample = expr_matrix(sample_idx, :);
else
    expr_sample = expr_matrix;
    sample_idx = 1:n_cells;
end

%% Step 2: Highly variable gene selection (dispersion-based)
fprintf('      Selecting highly variable genes (HVG)...\n');

% Convert to full for computation
if issparse(expr_sample)
    expr_full = full(expr_sample);
else
    expr_full = expr_sample;
end

% Calculate variance-to-mean ratio (dispersion)
gene_means = mean(expr_full, 1);
gene_vars = var(expr_full, 0, 1);

% Avoid division by zero
gene_means(gene_means == 0) = 1e-10;
dispersion = gene_vars ./ gene_means;

% Select top variable genes
[~, var_idx] = sort(dispersion, 'descend');
n_hvg = min(1000, floor(n_genes * 0.3)); % Top 30% or 1000 genes
hvg_idx = var_idx(1:n_hvg);

expr_hvg = expr_full(:, hvg_idx);
fprintf('      Selected %d highly variable genes\n', n_hvg);

%% Step 3: PCA dimensionality reduction
fprintf('      Running PCA for dimensionality reduction...\n');

% Center and scale
expr_centered = expr_hvg - mean(expr_hvg, 1);
expr_scaled = expr_centered ./ (std(expr_hvg, 0, 1) + 1e-10);
expr_scaled(isnan(expr_scaled)) = 0;

% PCA
[coeff, score, ~, ~, explained] = pca(expr_scaled, 'NumComponents', 50, 'Algorithm', 'svd');

% Select PCs explaining >80% variance
cumvar = cumsum(explained);
n_pcs = find(cumvar >= 80, 1, 'first');
n_pcs = max(15, min(n_pcs, 30)); % Between 15-30 PCs

score_pcs = score(:, 1:n_pcs);
fprintf('      Using %d PCs (%.1f%% variance explained)\n', n_pcs, cumvar(n_pcs));

%% Step 4: Shared nearest-neighbor (SNN) graph construction
fprintf('      Building SNN graph (k=30 nearest neighbors)...\n');

k_neighbors = 30;
n_sample = size(score_pcs, 1);

% Build k-NN graph using Euclidean distance
fprintf('      Computing pairwise distances...\n');

% Use pdist2 in chunks to save memory
chunk_size = 1000;
n_chunks = ceil(n_sample / chunk_size);
knn_indices = zeros(n_sample, k_neighbors);
knn_distances = zeros(n_sample, k_neighbors);

for chunk = 1:n_chunks
    start_idx = (chunk-1)*chunk_size + 1;
    end_idx = min(chunk*chunk_size, n_sample);
    chunk_cells = start_idx:end_idx;
    
    % Compute distances for this chunk
    dists = pdist2(score_pcs(chunk_cells, :), score_pcs, 'euclidean');
    
    % Find k nearest neighbors
    [sorted_dists, sorted_idx] = sort(dists, 2, 'ascend');
    knn_indices(chunk_cells, :) = sorted_idx(:, 2:k_neighbors+1); % Skip self
    knn_distances(chunk_cells, :) = sorted_dists(:, 2:k_neighbors+1);
    
    if mod(chunk, 10) == 0
        fprintf('      Progress: %d/%d chunks\n', chunk, n_chunks);
    end
end

% Convert to SNN graph (Jaccard similarity of neighborhoods)
fprintf('      Computing shared nearest neighbor weights...\n');

% Build adjacency matrix
snn_graph = zeros(n_sample, n_sample);

for i = 1:n_sample
    neighbors_i = knn_indices(i, :);
    
    for j_idx = 1:k_neighbors
        j = neighbors_i(j_idx);
        if j > i % Only compute upper triangle
            neighbors_j = knn_indices(j, :);
            
            % Jaccard similarity
            shared = length(intersect(neighbors_i, neighbors_j));
            union_size = 2*k_neighbors - shared;
            jaccard = shared / union_size;
            
            snn_graph(i, j) = jaccard;
            snn_graph(j, i) = jaccard;
        end
    end
    
    if mod(i, 500) == 0
        fprintf('      SNN progress: %d/%d cells\n', i, n_sample);
    end
end

fprintf('      SNN graph constructed (%.2f%% edges)\n', 100*nnz(snn_graph)/numel(snn_graph));

%% Step 5: Louvain community detection
fprintf('      Running Louvain algorithm for community detection...\n');

% Threshold graph to make sparse
threshold = prctile(snn_graph(snn_graph > 0), 50); % Keep top 50%
snn_sparse = sparse(snn_graph .* (snn_graph > threshold));

% Louvain clustering via modularity optimization
cluster_labels_sample = louvain_clustering(snn_sparse);
n_clusters = length(unique(cluster_labels_sample));

fprintf('      Detected %d communities via Louvain\n', n_clusters);

%% Step 6: Silhouette-score validation
fprintf('      Validating clusters with silhouette analysis...\n');

silhouette_vals = silhouette(score_pcs, cluster_labels_sample, 'euclidean');
avg_silhouette = mean(silhouette_vals);

fprintf('      Average silhouette score: %.3f ', avg_silhouette);
if avg_silhouette > 0.5
    fprintf('(Excellent)\n');
elseif avg_silhouette > 0.3
    fprintf('(Good)\n');
else
    fprintf('(Fair)\n');
end

%% Step 7: Assign all cells to nearest community centroid
fprintf('      Assigning all cells to nearest cluster centroids...\n');

% Compute cluster centroids in PCA space
cluster_centroids = zeros(n_clusters, n_pcs);
for k = 1:n_clusters
    cluster_centroids(k, :) = mean(score_pcs(cluster_labels_sample == k, :), 1);
end

% Assign remaining cells
cluster_labels_all = zeros(n_cells, 1);
cluster_labels_all(sample_idx) = cluster_labels_sample;

if length(sample_idx) < n_cells
    remaining = setdiff(1:n_cells, sample_idx);
    
    % Project remaining cells to PCA space
    expr_remaining = full(expr_matrix(remaining, hvg_idx));
    expr_remaining = (expr_remaining - mean(expr_hvg, 1)) ./ (std(expr_hvg, 0, 1) + 1e-10);
    expr_remaining(isnan(expr_remaining)) = 0;
    score_remaining = expr_remaining * coeff(:, 1:n_pcs);
    
    % Assign to nearest centroid
    dists_to_centroids = pdist2(score_remaining, cluster_centroids, 'euclidean');
    [~, nearest_cluster] = min(dists_to_centroids, [], 2);
    cluster_labels_all(remaining) = nearest_cluster;
    
    fprintf('      Assigned %d remaining cells\n', length(remaining));
end

%% OUTPUT
cell_ids = (1:n_cells)';
clusters = table(cell_ids, cluster_labels_all, 'VariableNames', {'Cell_ID', 'Cluster'});

% Summary
fprintf('\n      Clustering complete: %d clusters, silhouette=%.3f\n\n', n_clusters, avg_silhouette);

end

%% Local function: greedy Louvain modularity optimization

function communities = louvain_clustering(adjacency_matrix)
% Greedy modularity optimization following the first phase of the Louvain
% algorithm (Blondel et al. 2008). Each node is iteratively reassigned to
% the neighboring community that yields the greatest modularity gain until
% no further improvement is possible or max_iter is reached.

n = size(adjacency_matrix, 1);
communities = (1:n)'; % Initialize each node in own community

% Total edge weight
m = full(sum(adjacency_matrix(:))) / 2;
if m == 0
    return; % No edges
end

% Node degrees
degrees = full(sum(adjacency_matrix, 2));

% Iterative optimization
max_iter = 20;
improvement = true;
iter = 0;

while improvement && iter < max_iter
    iter = iter + 1;
    improvement = false;
    
    % Try moving each node to neighboring communities
    for i = randperm(n) % Random order
        current_comm = communities(i);
        
        % Find neighboring communities
        neighbors = find(adjacency_matrix(i, :) > 0);
        neighbor_comms = unique(communities(neighbors));
        neighbor_comms(neighbor_comms == current_comm) = [];
        
        if isempty(neighbor_comms)
            continue;
        end
        
        % Calculate modularity gain for each move
        best_gain = 0;
        best_comm = current_comm;
        
        for new_comm = neighbor_comms'
            % Modularity gain formula
            k_i_in = sum(adjacency_matrix(i, communities == new_comm));
            sum_tot = sum(degrees(communities == new_comm));
            
            delta_q = (k_i_in - degrees(i) * sum_tot / (2*m)) / m;
            
            if delta_q > best_gain
                best_gain = delta_q;
                best_comm = new_comm;
            end
        end
        
        % Move if beneficial
        if best_comm ~= current_comm
            communities(i) = best_comm;
            improvement = true;
        end
    end
end

% Renumber communities consecutively
[~, ~, communities] = unique(communities);

end