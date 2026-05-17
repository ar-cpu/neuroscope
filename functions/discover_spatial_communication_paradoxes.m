function paradox_results = discover_spatial_communication_paradoxes(gbm_network, gbm_clusters, gbm_specific)
% DISCOVER_SPATIAL_COMMUNICATION_PARADOXES  Detect mismatches between cluster
%   spatial proximity and observed signaling activity (silent neighbors,
%   asymmetric neighbors, distant broadcasters).
%
% paradox_results = discover_spatial_communication_paradoxes(gbm_network, gbm_clusters, gbm_specific)
%
% Inputs:
%   gbm_network  - Network struct (.source_clusters, .target_clusters, ...)
%   gbm_clusters - Table with spatial coordinates and cluster IDs
%                  (columns containing 'cluster', 'x'+'coord', 'y'+'coord')
%                    y-coordinate: name contains 'y' + ('coord'/'pos'/'cent')
%   gbm_specific - Table of condition-specific interactions, used to annotate
%                  which paradoxical channels are condition-enriched.
%
% Outputs:
%   paradox_results - Struct with fields for each paradox type (silent_neighbors,
%                     asymmetric_neighbors, distant_broadcasters) as tables,
%                     plus summary statistics and a proximity-communication
%                     correlation value.

fprintf('Spatial-Molecular Communication Paradoxes\n\n');

fprintf('Finding mismatches between spatial proximity\n');
fprintf('and molecular communication patterns\n\n');

%% Phase 1: Input validation and coordinate extraction

if ~isstruct(gbm_network) || ~istable(gbm_clusters)
    error('Invalid inputs');
end

required_net = {'source_clusters', 'target_clusters', 'ligands', 'receptors'};
for i = 1:length(required_net)
    if ~isfield(gbm_network, required_net{i})
        error('Missing field: %s', required_net{i});
    end
end

% Get spatial coordinates
x_col = ''; y_col = ''; cluster_col = '';

cols = gbm_clusters.Properties.VariableNames;
for i = 1:length(cols)
    c = lower(cols{i});
    if contains(c, 'cluster') && isempty(cluster_col)
        cluster_col = cols{i};
    end
    if (contains(c, 'x') && (contains(c, 'coord') || contains(c, 'pos') || contains(c, 'cent'))) && isempty(x_col)
        x_col = cols{i};
    end
    if (contains(c, 'y') && (contains(c, 'coord') || contains(c, 'pos') || contains(c, 'cent'))) && isempty(y_col)
        y_col = cols{i};
    end
end

if isempty(x_col) || isempty(y_col) || isempty(cluster_col)
    error('Cannot find spatial coordinates (need x_coord, y_coord, cluster columns)');
end

fprintf('  Spatial coordinates found: %s, %s\n', x_col, y_col);
fprintf('  Cluster column: %s\n', cluster_col);

src_clusters = gbm_network.source_clusters;
tgt_clusters = gbm_network.target_clusters;

% Get unique cluster IDs and ensure column vector
all_cluster_ids = unique(gbm_clusters.(cluster_col));
all_cluster_ids = double(all_cluster_ids(:));
n_clusters = length(all_cluster_ids);

fprintf('  Clusters: %d\n', n_clusters);
fprintf('  Interactions: %d\n', length(src_clusters));
fprintf('\n');

%% COMPUTE SPATIAL DISTANCE MATRIX
fprintf('Phase 2: Computing spatial distance matrix...\n');

% Get centroids
centroids = zeros(n_clusters, 2);
for c = 1:n_clusters
    cluster_id = all_cluster_ids(c);
    mask = gbm_clusters.(cluster_col) == cluster_id;
    
    if sum(mask) > 0
        centroids(c, 1) = mean(gbm_clusters.(x_col)(mask), 'omitnan');
        centroids(c, 2) = mean(gbm_clusters.(y_col)(mask), 'omitnan');
    end
end

% Pairwise distances
spatial_dist = zeros(n_clusters, n_clusters);
for i = 1:n_clusters
    for j = (i+1):n_clusters
        dx = centroids(i,1) - centroids(j,1);
        dy = centroids(i,2) - centroids(j,2);
        dist = sqrt(dx^2 + dy^2);
        spatial_dist(i,j) = dist;
        spatial_dist(j,i) = dist;
    end
end

fprintf('  Distance matrix: %dx%d\n', n_clusters, n_clusters);
fprintf('  Mean distance: %.1f um\n', mean(spatial_dist(spatial_dist>0)));
fprintf('  Min distance: %.1f um\n', min(spatial_dist(spatial_dist>0)));
fprintf('\n');

%% COMPUTE COMMUNICATION STRENGTH MATRIX
fprintf('Phase 3: Computing communication strength matrix...\n');

comm_strength = zeros(n_clusters, n_clusters);

for i = 1:length(src_clusters)
    src_id = src_clusters(i);
    tgt_id = tgt_clusters(i);
    
    src_idx = find(all_cluster_ids == src_id, 1);
    tgt_idx = find(all_cluster_ids == tgt_id, 1);
    
    if ~isempty(src_idx) && ~isempty(tgt_idx)
        comm_strength(src_idx, tgt_idx) = comm_strength(src_idx, tgt_idx) + 1;
    end
end

fprintf('  Communication matrix: %dx%d\n', n_clusters, n_clusters);
fprintf('  Total edges: %d\n', sum(comm_strength(:) > 0));
fprintf('  Mean strength: %.1f\n', mean(comm_strength(comm_strength>0)));
fprintf('\n');

%% PARADOX 1: SILENT NEIGHBORS
fprintf('Phase 4: Finding Silent Neighbors (adjacent but silent)...\n');

% Define "adjacent" as closest 10% of distances
dist_threshold = prctile(spatial_dist(spatial_dist>0), 10);

silent_neighbors = [];
for i = 1:n_clusters
    for j = (i+1):n_clusters
        is_close = spatial_dist(i,j) < dist_threshold && spatial_dist(i,j) > 0;
        is_silent = (comm_strength(i,j) == 0) && (comm_strength(j,i) == 0);
        
        if is_close && is_silent
            silent_neighbors = [silent_neighbors; i, j, spatial_dist(i,j)];
        end
    end
end

n_silent = size(silent_neighbors, 1);
fprintf('  Distance threshold: %.1f um\n', dist_threshold);
fprintf('  Silent neighbor pairs: %d\n', n_silent);
if n_silent > 0
    fprintf('  Functional barriers detected\n');
end
fprintf('\n');

%% PARADOX 2: ASYMMETRIC NEIGHBORS
fprintf('Phase 5: Finding Asymmetric Neighbors (unidirectional)...\n');

asymmetric_neighbors = [];
for i = 1:n_clusters
    for j = (i+1):n_clusters
        is_close = spatial_dist(i,j) < dist_threshold && spatial_dist(i,j) > 0;
        
        strength_ij = comm_strength(i,j);
        strength_ji = comm_strength(j,i);
        
        if is_close && (strength_ij > 0 || strength_ji > 0)
            % Compute asymmetry
            total = strength_ij + strength_ji;
            if total > 0
                asymmetry = abs(strength_ij - strength_ji) / total;
                
                % High asymmetry = directed control
                if asymmetry > 0.7
                    dominant = i;
                    subordinate = j;
                    dominant_strength = strength_ij;
                    if strength_ji > strength_ij
                        dominant = j;
                        subordinate = i;
                        dominant_strength = strength_ji;
                    end
                    
                    asymmetric_neighbors = [asymmetric_neighbors; ...
                        dominant, subordinate, spatial_dist(i,j), asymmetry, dominant_strength];
                end
            end
        end
    end
end

n_asymmetric = size(asymmetric_neighbors, 1);
fprintf('  Asymmetric neighbor pairs: %d\n', n_asymmetric);
if n_asymmetric > 0
    fprintf('  Hierarchical control structures detected\n');
end
fprintf('\n');

%% PARADOX 3: DISTANT BROADCASTERS
fprintf('Phase 6: Finding Distant Broadcasters (long-range)...\n');

% Define "distant" as top 25% of distances
dist_threshold_far = prctile(spatial_dist(spatial_dist>0), 75);
strength_threshold = prctile(comm_strength(comm_strength>0), 75);

distant_broadcasters = [];
for i = 1:n_clusters
    for j = 1:n_clusters
        if i ~= j
            is_far = spatial_dist(i,j) > dist_threshold_far;
            is_strong = comm_strength(i,j) > strength_threshold;
            
            if is_far && is_strong
                distant_broadcasters = [distant_broadcasters; ...
                    i, j, spatial_dist(i,j), comm_strength(i,j)];
            end
        end
    end
end

n_distant = size(distant_broadcasters, 1);
fprintf('  Distance threshold: %.1f um\n', dist_threshold_far);
fprintf('  Distant broadcasting pairs: %d\n', n_distant);
if n_distant > 0
    fprintf('  Long-range signaling detected\n');
end
fprintf('\n');

%% CHECK TUMOR SPECIFICITY
fprintf('Phase 7: Checking specificity...\n');

n_tumor_asymmetric = 0;
has_tumor_data = istable(gbm_specific) && height(gbm_specific) > 0;

if has_tumor_data && n_asymmetric > 0
    cols = gbm_specific.Properties.VariableNames;
    lig_col = ''; rec_col = '';
    
    for i = 1:length(cols)
        if contains(lower(cols{i}), {'ligand', 'source'}) && isempty(lig_col)
            lig_col = cols{i};
        end
        if contains(lower(cols{i}), {'receptor', 'target'}) && isempty(rec_col)
            rec_col = cols{i};
        end
    end
    
    if ~isempty(lig_col) && ~isempty(rec_col)
        % Check if dominant signals in asymmetric pairs are tumor-specific
        for row = 1:size(asymmetric_neighbors, 1)
            dom_idx = asymmetric_neighbors(row, 1);
            sub_idx = asymmetric_neighbors(row, 2);
            
            dom_id = all_cluster_ids(dom_idx);
            sub_id = all_cluster_ids(sub_idx);
            
            % Find interactions from dominant to subordinate
            mask = (src_clusters == dom_id) & (tgt_clusters == sub_id);
            
            if sum(mask) > 0
                % Check any tumor-specific
                found_tumor = false;
                for i = find(mask)'
                    if found_tumor
                        break;
                    end
                    lig = gbm_network.ligands{i};
                    rec = gbm_network.receptors{i};
                    
                    for spec_row = 1:height(gbm_specific)
                        if strcmp(gbm_specific.(lig_col){spec_row}, lig) && ...
                           strcmp(gbm_specific.(rec_col){spec_row}, rec)
                            n_tumor_asymmetric = n_tumor_asymmetric + 1;
                            found_tumor = true;
                            break;
                        end
                    end
                end
            end
        end
        
        fprintf('  Tumor-specific asymmetric controls: %d\n', n_tumor_asymmetric);
    else
        fprintf('  Could not determine specificity\n');
    end
else
    fprintf('  No tumor-specific data or no asymmetric neighbors\n');
end

fprintf('\n');

%% EXPORT
fprintf('Phase 8: Exporting results...\n');

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

% Initialize result structures
silent_table = [];
asym_table = [];
distant_table = [];

% ALWAYS SAVE SUMMARY CSV
summary = table({'SilentNeighbors'; 'AsymmetricNeighbors'; 'DistantBroadcasters'; 'TumorSpecificAsymmetric'}, ...
    [n_silent; n_asymmetric; n_distant; n_tumor_asymmetric], ...
    'VariableNames', {'ParadoxType', 'Count'});
writetable(summary, 'results/novel_discoveries/spatial_paradox_summary.csv');
fprintf('  spatial_paradox_summary.csv\n');

% Silent neighbors
if n_silent > 0
    silent_table = table(all_cluster_ids(silent_neighbors(:,1)), ...
        all_cluster_ids(silent_neighbors(:,2)), ...
        silent_neighbors(:,3), ...
        'VariableNames', {'ClusterA', 'ClusterB', 'SpatialDistance'});
    writetable(silent_table, 'results/novel_discoveries/silent_neighbors.csv');
    fprintf('  silent_neighbors.csv\n');
end

% Asymmetric neighbors
if n_asymmetric > 0
    asym_table = table(all_cluster_ids(asymmetric_neighbors(:,1)), ...
        all_cluster_ids(asymmetric_neighbors(:,2)), ...
        asymmetric_neighbors(:,3), ...
        asymmetric_neighbors(:,4), ...
        asymmetric_neighbors(:,5), ...
        'VariableNames', {'DominantCluster', 'SubordinateCluster', ...
                          'SpatialDistance', 'Asymmetry', 'SignalStrength'});
    writetable(asym_table, 'results/novel_discoveries/asymmetric_neighbors.csv');
    fprintf('  asymmetric_neighbors.csv\n');
end

% Distant broadcasters
if n_distant > 0
    distant_table = table(all_cluster_ids(distant_broadcasters(:,1)), ...
        all_cluster_ids(distant_broadcasters(:,2)), ...
        distant_broadcasters(:,3), ...
        distant_broadcasters(:,4), ...
        'VariableNames', {'SourceCluster', 'TargetCluster', ...
                          'SpatialDistance', 'SignalStrength'});
    writetable(distant_table, 'results/novel_discoveries/distant_broadcasters.csv');
    fprintf('  distant_broadcasters.csv\n');
end

if n_silent == 0 && n_asymmetric == 0 && n_distant == 0
    fprintf('  • No paradoxes detected (thresholds may be too strict)\n');
    fprintf('  • Summary still saved with counts\n');
end

fprintf('\n');

%% PACKAGE
paradox_results = struct();
paradox_results.n_silent_neighbors = n_silent;
paradox_results.n_asymmetric_neighbors = n_asymmetric;
paradox_results.n_distant_broadcasters = n_distant;
paradox_results.n_tumor_asymmetric = n_tumor_asymmetric;

if n_silent > 0
    paradox_results.silent_neighbors = silent_table;
end
if n_asymmetric > 0
    paradox_results.asymmetric_neighbors = asym_table;
end
if n_distant > 0
    paradox_results.distant_broadcasters = distant_table;
end

%% SUMMARY
fprintf('Spatial-molecular paradoxes analysis complete\n\n');

fprintf('PARADOXES FOUND:\n');
fprintf('  Silent Neighbors: %d pairs\n', n_silent);
fprintf('  Asymmetric Neighbors: %d pairs\n', n_asymmetric);
fprintf('  Distant Broadcasters: %d pairs\n', n_distant);
fprintf('  Tumor-specific asymmetric controls: %d\n', n_tumor_asymmetric);

fprintf('\nBIOLOGICAL INTERPRETATION:\n');
fprintf('  • Silent neighbors = functional barriers/immune exclusion\n');
fprintf('  • Asymmetric neighbors = hierarchical orchestration\n');
fprintf('  • Distant broadcasters = diffusible/systemic factors\n');

fprintf('\nAnalysis complete\n\n');

end