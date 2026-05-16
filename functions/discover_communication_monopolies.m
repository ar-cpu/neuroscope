function monopoly_results = discover_communication_monopolies(gbm_network, gbm_clusters, gbm_specific)
% DISCOVER_COMMUNICATION_MONOPOLIES  Characterize the redundancy structure
%   of the cell-cell communication network by identifying exclusive signaling
%   channels (monopolies), high-traffic cluster-pair bottlenecks, and
%   low-redundancy signals that are vulnerable to therapeutic disruption.
%
% Definitions
% -----------
%   Monopoly        — A ligand-receptor pair that is used by exactly one
%                     source→target cluster pair. These represent private
%                     communication channels; disrupting such a pair eliminates
%                     the signal entirely rather than merely reducing it.
%   Bottleneck      — A cluster pair (src, tgt) that carries the highest
%                     aggregate interaction count. Removing bottleneck cluster
%                     pairs would maximally fragment the network.
%   Low-redundancy  — L-R pairs used by ≤3 cluster pairs but more than one
%                     (i.e., they exist in the network but are not widely
%                     shared). These are intermediate-vulnerability targets.
%
% Condition-specificity of monopolized channels is determined by cross-
% referencing with the gbm_specific table (FDR-significant interactions).
%
% Usage:
%   monopoly_results = discover_communication_monopolies(network, clusters, condition_specific)
%
% Inputs:
%   gbm_network  - Communication network struct with fields:
%                    .source_clusters - [N×1] sending cluster IDs
%                    .target_clusters - [N×1] receiving cluster IDs
%                    .ligands         - {N×1} ligand gene symbols
%                    .receptors       - {N×1} receptor gene symbols
%   gbm_clusters - Cluster assignment table. Not used quantitatively here;
%                  retained for API consistency with other discovery functions.
%   gbm_specific - Table of condition-specific interactions. Used to annotate
%                  whether monopolized channels are condition-enriched.
%
% Outputs:
%   monopoly_results.n_lr_pairs              - Total unique L-R pairs
%   monopoly_results.n_monopolies            - Exclusive-channel L-R pairs
%   monopoly_results.monopoly_fraction       - Fraction of L-R pairs monopolized
%   monopoly_results.n_tumor_monopolies      - Condition-specific monopolies
%   monopoly_results.n_low_redundancy        - L-R pairs with 2–3 cluster users
%   monopoly_results.n_bottlenecks           - Top bottleneck cluster pairs (≤50)
%   monopoly_results.summary                 - Key-metric summary table
%   monopoly_results.monopoly_table          - Monopoly details (if any)
%   monopoly_results.bottleneck_table        - Bottleneck cluster-pair table
%   monopoly_results.vulnerable_table        - Low-redundancy signals (if any)
%   monopoly_results.mean_redundancy         - Mean cluster pairs per L-R pair
%   monopoly_results.max_bottleneck_strength - Highest cluster-pair edge weight

fprintf('Communication monopoly and bottleneck analysis\n\n');

%% Phase 1: Input validation

if ~isstruct(gbm_network)
    error('gbm_network must be struct');
end

required_fields = {'source_clusters', 'target_clusters', 'ligands', 'receptors'};
for i = 1:length(required_fields)
    if ~isfield(gbm_network, required_fields{i})
        error('Missing required field: %s', required_fields{i});
    end
end

src_clusters = gbm_network.source_clusters;
tgt_clusters = gbm_network.target_clusters;
ligands = gbm_network.ligands;
receptors = gbm_network.receptors;

n_interactions = length(src_clusters);

fprintf('Network structure validated\n');
fprintf('Total interactions: %d\n', n_interactions);
fprintf('Source clusters: %d unique\n', length(unique(src_clusters)));
fprintf('Target clusters: %d unique\n', length(unique(tgt_clusters)));
fprintf('\n');

%% Phase 2: Build L-R pair catalog and cluster-pair usage map
fprintf('Phase 2: Cataloging ligand-receptor pairs...\n');

% Create unique L-R pair identifiers
lr_pairs = cell(n_interactions, 1);
for i = 1:n_interactions
    lr_pairs{i} = [ligands{i}, '->', receptors{i}];
end

unique_lr_pairs = unique(lr_pairs);
n_lr_pairs = length(unique_lr_pairs);

fprintf('Unique L-R pairs: %d\n', n_lr_pairs);
fprintf('Building pair usage matrix...\n');

% For each L-R pair, find which cluster pairs use it
pair_users = cell(n_lr_pairs, 1);  % Cell array of cluster pairs using each L-R

for p = 1:n_lr_pairs
    lr_pair = unique_lr_pairs{p};
    
    % Find all interactions with this L-R pair
    mask = strcmp(lr_pairs, lr_pair);
    
    % Get unique cluster pairs using this L-R
    src_using = src_clusters(mask);
    tgt_using = tgt_clusters(mask);
    
    % Store as unique cluster pair identifiers
    cluster_pairs = unique([src_using, tgt_using], 'rows');
    pair_users{p} = cluster_pairs;
end

fprintf('Pair usage catalog complete\n\n');

%% Phase 3: Identify monopolized channels
fprintf('Phase 3: Identifying communication monopolies...\n');

% MONOPOLY = L-R pair used by EXACTLY ONE cluster pair
monopolies = [];
monopoly_lr_pairs = {};
monopoly_cluster_pairs = [];

for p = 1:n_lr_pairs
    n_users = size(pair_users{p}, 1);
    
    if n_users == 1
        % This is a MONOPOLY - exclusive channel
        monopolies = [monopolies; p];
        monopoly_lr_pairs{end+1} = unique_lr_pairs{p};
        monopoly_cluster_pairs = [monopoly_cluster_pairs; pair_users{p}];
    end
end

n_monopolies = length(monopolies);
monopoly_fraction = n_monopolies / n_lr_pairs;

fprintf('Monopolized L-R pairs: %d / %d (%.1f%%)\n', ...
    n_monopolies, n_lr_pairs, 100*monopoly_fraction);

if n_monopolies > 0
    fprintf('Exclusive communication channels detected\n');
else
    fprintf('No complete monopolies (all signals shared)\n');
end
fprintf('\n');

%% Phase 4: Redundancy analysis
fprintf('Phase 4: Analyzing communication redundancy...\n');

% For each L-R pair, count how many cluster pairs can use it
redundancy_counts = zeros(n_lr_pairs, 1);
for p = 1:n_lr_pairs
    redundancy_counts(p) = size(pair_users{p}, 1);
end

fprintf('Mean redundancy: %.1f cluster pairs per signal\n', mean(redundancy_counts));
fprintf('Max redundancy: %d cluster pairs\n', max(redundancy_counts));
fprintf('Low redundancy (≤3 pairs): %d signals (%.1f%%)\n', ...
    sum(redundancy_counts <= 3), 100*sum(redundancy_counts <= 3)/n_lr_pairs);

% Identify low-redundancy signals (vulnerable)
low_redundancy_idx = find(redundancy_counts <= 3 & redundancy_counts > 1);
n_low_redundancy = length(low_redundancy_idx);

fprintf('Vulnerable signals identified: %d\n', n_low_redundancy);
fprintf('\n');

%% Phase 5: Identify cluster-pair bottlenecks
fprintf('Phase 5: Identifying network bottlenecks...\n');

% Build cluster pair interaction strength matrix
all_cluster_ids = unique([src_clusters; tgt_clusters]);
n_clusters = length(all_cluster_ids);

% Map cluster IDs to indices
cluster_to_idx = containers.Map(all_cluster_ids, 1:n_clusters);

% Count interactions between each cluster pair
pair_strength = sparse(n_clusters, n_clusters);
for i = 1:n_interactions
    src_idx = cluster_to_idx(src_clusters(i));
    tgt_idx = cluster_to_idx(tgt_clusters(i));
    pair_strength(src_idx, tgt_idx) = pair_strength(src_idx, tgt_idx) + 1;
end

% BOTTLENECK = Cluster pair with high interaction count
% (removing it would eliminate many signals)
[src_idx, tgt_idx, strengths] = find(pair_strength);

% Sort by strength
[sorted_strengths, sort_order] = sort(strengths, 'descend');

% Top bottlenecks
n_bottlenecks = min(50, length(sorted_strengths));
bottleneck_pairs = zeros(n_bottlenecks, 3);

for i = 1:n_bottlenecks
    idx = sort_order(i);
    bottleneck_pairs(i, 1) = all_cluster_ids(src_idx(idx));
    bottleneck_pairs(i, 2) = all_cluster_ids(tgt_idx(idx));
    bottleneck_pairs(i, 3) = sorted_strengths(i);
end

fprintf('Top bottleneck pairs identified: %d\n', n_bottlenecks);
fprintf('Max interaction count: %d\n', max(sorted_strengths));
fprintf('Mean top bottleneck strength: %.1f\n', mean(sorted_strengths(1:n_bottlenecks)));
fprintf('\n');

%% Phase 6: Cross-reference monopolies with condition-specific interactions

n_tumor_monopolies = 0;
tumor_monopoly_idx = [];

if istable(gbm_specific) && height(gbm_specific) > 0
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
    
    if ~isempty(lig_col) && ~isempty(rec_col) && n_monopolies > 0
        for m = 1:n_monopolies
            monopoly_pair = monopoly_lr_pairs{m};
            parts = strsplit(monopoly_pair, '->');
            
            if length(parts) == 2
                lig = parts{1};
                rec = parts{2};
                
                for row = 1:height(gbm_specific)
                    if strcmp(gbm_specific.(lig_col){row}, lig) && ...
                       strcmp(gbm_specific.(rec_col){row}, rec)
                        n_tumor_monopolies = n_tumor_monopolies + 1;
                        tumor_monopoly_idx = [tumor_monopoly_idx; m];
                        break;
                    end
                end
            end
        end
        
        fprintf('Tumor-specific monopolies: %d / %d (%.1f%%)\n', ...
            n_tumor_monopolies, n_monopolies, 100*n_tumor_monopolies/max(1,n_monopolies));
    else
        fprintf('Could not determine tumor specificity\n');
    end
else
    fprintf('No tumor-specific data available\n');
end

fprintf('\n');

%% Phase 7: Assemble results tables

% Monopoly table
if n_monopolies > 0
    monopoly_table = table();
    monopoly_table.LR_Pair = monopoly_lr_pairs';
    monopoly_table.SourceCluster = monopoly_cluster_pairs(:, 1);
    monopoly_table.TargetCluster = monopoly_cluster_pairs(:, 2);
    
    % Count how many times this monopoly is used
    monopoly_counts = zeros(n_monopolies, 1);
    for m = 1:n_monopolies
        lr = monopoly_lr_pairs{m};
        src = monopoly_cluster_pairs(m, 1);
        tgt = monopoly_cluster_pairs(m, 2);
        
        mask = strcmp(lr_pairs, lr) & src_clusters == src & tgt_clusters == tgt;
        monopoly_counts(m) = sum(mask);
    end
    
    monopoly_table.InteractionCount = monopoly_counts;
    monopoly_table.TumorSpecific = false(n_monopolies, 1);
    
    if ~isempty(tumor_monopoly_idx)
        monopoly_table.TumorSpecific(tumor_monopoly_idx) = true;
    end
    
    monopoly_table = sortrows(monopoly_table, 'InteractionCount', 'descend');
else
    monopoly_table = table();
end

% Bottleneck table
bottleneck_table = table();
bottleneck_table.SourceCluster = bottleneck_pairs(:, 1);
bottleneck_table.TargetCluster = bottleneck_pairs(:, 2);
bottleneck_table.InteractionCount = bottleneck_pairs(:, 3);

% For each bottleneck, count unique L-R pairs
bottleneck_diversity = zeros(n_bottlenecks, 1);
for b = 1:n_bottlenecks
    src = bottleneck_pairs(b, 1);
    tgt = bottleneck_pairs(b, 2);
    
    mask = src_clusters == src & tgt_clusters == tgt;
    unique_pairs = unique(lr_pairs(mask));
    bottleneck_diversity(b) = length(unique_pairs);
end

bottleneck_table.UniqueLRPairs = bottleneck_diversity;

% Redundancy table (low redundancy signals)
if n_low_redundancy > 0
    vulnerable_table = table();
    vulnerable_table.LR_Pair = unique_lr_pairs(low_redundancy_idx);
    vulnerable_table.NumClusterPairs = redundancy_counts(low_redundancy_idx);
    vulnerable_table = sortrows(vulnerable_table, 'NumClusterPairs', 'ascend');
else
    vulnerable_table = table();
end

fprintf('Results tables created\n\n');

%% Phase 8: Export

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

% Summary
summary = table();
summary.Metric = {'TotalLRPairs'; 'Monopolies'; 'TumorSpecificMonopolies'; ...
                  'LowRedundancySignals'; 'TopBottlenecks'};
summary.Count = [n_lr_pairs; n_monopolies; n_tumor_monopolies; ...
                 n_low_redundancy; n_bottlenecks];
writetable(summary, 'results/novel_discoveries/monopoly_summary.csv');
fprintf('  monopoly_summary.csv\n');

if n_monopolies > 0
    writetable(monopoly_table, 'results/novel_discoveries/monopolies.csv');
    fprintf('monopolies.csv\n');
end

writetable(bottleneck_table, 'results/novel_discoveries/bottlenecks.csv');
fprintf('bottlenecks.csv\n');

if n_low_redundancy > 0
    writetable(vulnerable_table, 'results/novel_discoveries/vulnerable_signals.csv');
    fprintf('vulnerable_signals.csv\n');
end

fprintf('\n');

%% Package
monopoly_results = struct();
monopoly_results.n_lr_pairs = n_lr_pairs;
monopoly_results.n_monopolies = n_monopolies;
monopoly_results.monopoly_fraction = monopoly_fraction;
monopoly_results.n_tumor_monopolies = n_tumor_monopolies;
monopoly_results.n_low_redundancy = n_low_redundancy;
monopoly_results.n_bottlenecks = n_bottlenecks;
monopoly_results.summary = summary;

if n_monopolies > 0
    monopoly_results.monopoly_table = monopoly_table;
end
monopoly_results.bottleneck_table = bottleneck_table;
if n_low_redundancy > 0
    monopoly_results.vulnerable_table = vulnerable_table;
end

monopoly_results.mean_redundancy = mean(redundancy_counts);
monopoly_results.max_bottleneck_strength = max(sorted_strengths);

%% Summary

fprintf('Communication monopoly analysis complete\n\n');
fprintf('Results:\n');
fprintf('  L-R pairs analyzed:        %d\n', n_lr_pairs);
fprintf('  Monopolies (exclusive):    %d (%.1f%%)\n', n_monopolies, 100*monopoly_fraction);
fprintf('  Condition-specific mono:   %d\n', n_tumor_monopolies);
fprintf('  Low-redundancy signals:    %d\n', n_low_redundancy);
fprintf('  Bottleneck pairs:          %d\n\n', n_bottlenecks);

if n_monopolies > 0
    fprintf('Top 5 monopolized channels:\n');
    for i = 1:min(5, height(monopoly_table))
        fprintf('  %d. %s (Cluster %d->%d, %d interactions', ...
            i, monopoly_table.LR_Pair{i}, ...
            monopoly_table.SourceCluster(i), ...
            monopoly_table.TargetCluster(i), ...
            monopoly_table.InteractionCount(i));
        if monopoly_table.TumorSpecific(i)
            fprintf(', condition-specific');
        end
        fprintf(')\n');
    end
    fprintf('\n');
end

fprintf('Top 5 network bottlenecks:\n');
for i = 1:min(5, height(bottleneck_table))
    fprintf('  %d. Cluster %d->%d (%d interactions, %d unique L-R pairs)\n', ...
        i, bottleneck_table.SourceCluster(i), ...
        bottleneck_table.TargetCluster(i), ...
        bottleneck_table.InteractionCount(i), ...
        bottleneck_table.UniqueLRPairs(i));
end

fprintf('\n');

end