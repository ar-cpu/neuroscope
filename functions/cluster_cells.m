function clusters = cluster_cells(expr_matrix, metadata, dataset_name)
% CLUSTER_CELLS  Partition cells into communities using a Louvain + SNN
%   pipeline (Seurat/Scanpy-style): HVG selection → PCA → k-NN → SNN graph
%   → Louvain → silhouette validation → centroid assignment for held-out cells.
%
% Usage:
%   clusters = cluster_cells(expr_matrix, metadata, dataset_name)
%
% Inputs:
%   expr_matrix  - [cells × genes] sparse or dense normalized expression matrix.
%   metadata     - Table of cell metadata (unused; kept for API compatibility).
%   dataset_name - Label used in progress messages.
%
% Outputs:
%   clusters - Table with columns Cell_ID (1-based) and Cluster (1-based).

fprintf('      Graph-based clustering (Louvain + SNN)\n\n');

n_cells = size(expr_matrix, 1);
n_genes = size(expr_matrix, 2);

fprintf('      Dataset: %s (%d cells x %d genes)\n', dataset_name, n_cells, n_genes);

% Step 1: Subsample large datasets for memory-efficient graph construction
max_sample = 15000;
if n_cells > max_sample
    fprintf('      Subsampling to %d cells for graph construction...\n', max_sample);
    sample_idx = randsample(n_cells, max_sample);
    expr_sample = expr_matrix(sample_idx, :);
else
    expr_sample = expr_matrix;
    sample_idx = 1:n_cells;
end

% Step 2: HVG selection by variance-to-mean dispersion ratio
fprintf('      Selecting highly variable genes (HVG)...\n');

if issparse(expr_sample)
    expr_full = full(expr_sample);
else
    expr_full = expr_sample;
end

gene_means = mean(expr_full, 1);
gene_vars = var(expr_full, 0, 1);
gene_means(gene_means == 0) = 1e-10;
dispersion = gene_vars ./ gene_means;

[~, var_idx] = sort(dispersion, 'descend');
n_hvg = min(1000, floor(n_genes * 0.3));
hvg_idx = var_idx(1:n_hvg);

expr_hvg = expr_full(:, hvg_idx);
fprintf('      Selected %d highly variable genes\n', n_hvg);

% Step 3: PCA — retain PCs explaining >=80% variance, bounded to [15, 30]
fprintf('      Running PCA for dimensionality reduction...\n');

expr_centered = expr_hvg - mean(expr_hvg, 1);
expr_scaled = expr_centered ./ (std(expr_hvg, 0, 1) + 1e-10);
expr_scaled(isnan(expr_scaled)) = 0;

[coeff, score, ~, ~, explained] = pca(expr_scaled, 'NumComponents', 50, 'Algorithm', 'svd');

cumvar = cumsum(explained);
n_pcs = find(cumvar >= 80, 1, 'first');
n_pcs = max(15, min(n_pcs, 30));

score_pcs = score(:, 1:n_pcs);
fprintf('      Using %d PCs (%.1f%% variance explained)\n', n_pcs, cumvar(n_pcs));

% Step 4: SNN graph (k=30) with Jaccard edge weights
fprintf('      Building SNN graph (k=30 nearest neighbors)...\n');

k_neighbors = 30;
n_sample = size(score_pcs, 1);

fprintf('      Computing pairwise distances...\n');

chunk_size = 1000;
n_chunks = ceil(n_sample / chunk_size);
knn_indices = zeros(n_sample, k_neighbors);
knn_distances = zeros(n_sample, k_neighbors);

for chunk = 1:n_chunks
    start_idx = (chunk-1)*chunk_size + 1;
    end_idx = min(chunk*chunk_size, n_sample);
    chunk_cells = start_idx:end_idx;

    dists = pdist2(score_pcs(chunk_cells, :), score_pcs, 'euclidean');

    [sorted_dists, sorted_idx] = sort(dists, 2, 'ascend');
    knn_indices(chunk_cells, :) = sorted_idx(:, 2:k_neighbors+1);   % skip self
    knn_distances(chunk_cells, :) = sorted_dists(:, 2:k_neighbors+1);

    if mod(chunk, 10) == 0
        fprintf('      Progress: %d/%d chunks\n', chunk, n_chunks);
    end
end

fprintf('      Computing shared nearest neighbor weights...\n');

snn_graph = zeros(n_sample, n_sample);

for i = 1:n_sample
    neighbors_i = knn_indices(i, :);

    for j_idx = 1:k_neighbors
        j = neighbors_i(j_idx);
        if j > i
            neighbors_j = knn_indices(j, :);
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

% Step 5: Louvain community detection on thresholded SNN graph
fprintf('      Running Louvain algorithm for community detection...\n');

threshold = prctile(snn_graph(snn_graph > 0), 50);
snn_sparse = sparse(snn_graph .* (snn_graph > threshold));

cluster_labels_sample = louvain_clustering(snn_sparse);
n_clusters = length(unique(cluster_labels_sample));

fprintf('      Detected %d communities via Louvain\n', n_clusters);

% Step 6: Silhouette score as cluster quality metric
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

% Step 7: Assign held-out cells to nearest centroid in PCA space
fprintf('      Assigning all cells to nearest cluster centroids...\n');

cluster_centroids = zeros(n_clusters, n_pcs);
for k = 1:n_clusters
    cluster_centroids(k, :) = mean(score_pcs(cluster_labels_sample == k, :), 1);
end

cluster_labels_all = zeros(n_cells, 1);
cluster_labels_all(sample_idx) = cluster_labels_sample;

if length(sample_idx) < n_cells
    remaining = setdiff(1:n_cells, sample_idx);

    expr_remaining = full(expr_matrix(remaining, hvg_idx));
    expr_remaining = (expr_remaining - mean(expr_hvg, 1)) ./ (std(expr_hvg, 0, 1) + 1e-10);
    expr_remaining(isnan(expr_remaining)) = 0;
    score_remaining = expr_remaining * coeff(:, 1:n_pcs);

    dists_to_centroids = pdist2(score_remaining, cluster_centroids, 'euclidean');
    [~, nearest_cluster] = min(dists_to_centroids, [], 2);
    cluster_labels_all(remaining) = nearest_cluster;

    fprintf('      Assigned %d remaining cells\n', length(remaining));
end

cell_ids = (1:n_cells)';
clusters = table(cell_ids, cluster_labels_all, 'VariableNames', {'Cell_ID', 'Cluster'});

fprintf('\n      Clustering complete: %d clusters, silhouette=%.3f\n\n', n_clusters, avg_silhouette);

end

function communities = louvain_clustering(adjacency_matrix)
% Greedy modularity optimization — phase 1 of Blondel et al. (2008).

n = size(adjacency_matrix, 1);
communities = (1:n)';

m = full(sum(adjacency_matrix(:))) / 2;
if m == 0
    return;
end

degrees = full(sum(adjacency_matrix, 2));

max_iter = 20;
improvement = true;
iter = 0;

while improvement && iter < max_iter
    iter = iter + 1;
    improvement = false;

    for i = randperm(n)
        current_comm = communities(i);

        neighbors = find(adjacency_matrix(i, :) > 0);
        neighbor_comms = unique(communities(neighbors));
        neighbor_comms(neighbor_comms == current_comm) = [];

        if isempty(neighbor_comms)
            continue;
        end

        best_gain = 0;
        best_comm = current_comm;

        for new_comm = neighbor_comms'
            k_i_in = sum(adjacency_matrix(i, communities == new_comm));
            sum_tot = sum(degrees(communities == new_comm));
            delta_q = (k_i_in - degrees(i) * sum_tot / (2*m)) / m;

            if delta_q > best_gain
                best_gain = delta_q;
                best_comm = new_comm;
            end
        end

        if best_comm ~= current_comm
            communities(i) = best_comm;
            improvement = true;
        end
    end
end

[~, ~, communities] = unique(communities);

end
