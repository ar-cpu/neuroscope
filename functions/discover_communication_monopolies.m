function monopoly_results = discover_communication_monopolies(gbm_network, gbm_clusters, gbm_specific)
% DISCOVER_COMMUNICATION_MONOPOLIES  Identify exclusive signaling channels
%   (monopolies), high-traffic cluster-pair bottlenecks, and low-redundancy
%   signals vulnerable to therapeutic disruption.
%
% Usage:
%   monopoly_results = discover_communication_monopolies(network, clusters, condition_specific)
%
% Inputs:
%   gbm_network  - Network struct with fields: source_clusters, target_clusters,
%                  ligands, receptors.
%   gbm_clusters - Cluster assignment table (retained for API consistency).
%   gbm_specific - Table of condition-specific interactions for annotation.
%
% Outputs:
%   monopoly_results.n_monopolies       - L-R pairs used by exactly one cluster pair
%   monopoly_results.monopoly_fraction  - Fraction of L-R pairs monopolized
%   monopoly_results.n_tumor_monopolies - Condition-specific monopolies
%   monopoly_results.n_low_redundancy   - L-R pairs with 2-3 cluster users
%   monopoly_results.monopoly_table, bottleneck_table, vulnerable_table

fprintf('Communication monopoly and bottleneck analysis\n\n');

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

fprintf('Total interactions: %d\n', n_interactions);
fprintf('Source clusters: %d unique\n', length(unique(src_clusters)));
fprintf('Target clusters: %d unique\n\n', length(unique(tgt_clusters)));

fprintf('Phase 2: Cataloging ligand-receptor pairs...\n');

lr_pairs = cell(n_interactions, 1);
for i = 1:n_interactions
    lr_pairs{i} = [ligands{i}, '->', receptors{i}];
end

unique_lr_pairs = unique(lr_pairs);
n_lr_pairs = length(unique_lr_pairs);

fprintf('Unique L-R pairs: %d\n', n_lr_pairs);

pair_users = cell(n_lr_pairs, 1);

for p = 1:n_lr_pairs
    lr_pair = unique_lr_pairs{p};
    mask = strcmp(lr_pairs, lr_pair);
    src_using = src_clusters(mask);
    tgt_using = tgt_clusters(mask);
    cluster_pairs = unique([src_using, tgt_using], 'rows');
    pair_users{p} = cluster_pairs;
end

fprintf('\nPhase 3: Identifying communication monopolies...\n');

monopolies = [];
monopoly_lr_pairs = {};
monopoly_cluster_pairs = [];

for p = 1:n_lr_pairs
    if size(pair_users{p}, 1) == 1
        monopolies = [monopolies; p];
        monopoly_lr_pairs{end+1} = unique_lr_pairs{p};
        monopoly_cluster_pairs = [monopoly_cluster_pairs; pair_users{p}];
    end
end

n_monopolies = length(monopolies);
monopoly_fraction = n_monopolies / n_lr_pairs;

fprintf('Monopolized L-R pairs: %d / %d (%.1f%%)\n', ...
    n_monopolies, n_lr_pairs, 100*monopoly_fraction);

fprintf('\nPhase 4: Analyzing communication redundancy...\n');

redundancy_counts = zeros(n_lr_pairs, 1);
for p = 1:n_lr_pairs
    redundancy_counts(p) = size(pair_users{p}, 1);
end

fprintf('Mean redundancy: %.1f cluster pairs per signal\n', mean(redundancy_counts));
fprintf('Low redundancy (<=3 pairs): %d signals (%.1f%%)\n', ...
    sum(redundancy_counts <= 3), 100*sum(redundancy_counts <= 3)/n_lr_pairs);

low_redundancy_idx = find(redundancy_counts <= 3 & redundancy_counts > 1);
n_low_redundancy = length(low_redundancy_idx);

fprintf('Vulnerable signals identified: %d\n\n', n_low_redundancy);

fprintf('Phase 5: Identifying network bottlenecks...\n');

all_cluster_ids = unique([src_clusters; tgt_clusters]);
n_clusters = length(all_cluster_ids);
cluster_to_idx = containers.Map(all_cluster_ids, 1:n_clusters);

pair_strength = sparse(n_clusters, n_clusters);
for i = 1:n_interactions
    src_idx = cluster_to_idx(src_clusters(i));
    tgt_idx = cluster_to_idx(tgt_clusters(i));
    pair_strength(src_idx, tgt_idx) = pair_strength(src_idx, tgt_idx) + 1;
end

[src_idx, tgt_idx, strengths] = find(pair_strength);
[sorted_strengths, sort_order] = sort(strengths, 'descend');

n_bottlenecks = min(50, length(sorted_strengths));
bottleneck_pairs = zeros(n_bottlenecks, 3);

for i = 1:n_bottlenecks
    idx = sort_order(i);
    bottleneck_pairs(i, 1) = all_cluster_ids(src_idx(idx));
    bottleneck_pairs(i, 2) = all_cluster_ids(tgt_idx(idx));
    bottleneck_pairs(i, 3) = sorted_strengths(i);
end

fprintf('Top bottleneck pairs: %d | Max interaction count: %d\n\n', ...
    n_bottlenecks, max(sorted_strengths));

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
            parts = strsplit(monopoly_lr_pairs{m}, '->');
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

if n_monopolies > 0
    monopoly_table = table();
    monopoly_table.LR_Pair = monopoly_lr_pairs';
    monopoly_table.SourceCluster = monopoly_cluster_pairs(:, 1);
    monopoly_table.TargetCluster = monopoly_cluster_pairs(:, 2);

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

bottleneck_table = table();
bottleneck_table.SourceCluster = bottleneck_pairs(:, 1);
bottleneck_table.TargetCluster = bottleneck_pairs(:, 2);
bottleneck_table.InteractionCount = bottleneck_pairs(:, 3);

bottleneck_diversity = zeros(n_bottlenecks, 1);
for b = 1:n_bottlenecks
    src = bottleneck_pairs(b, 1);
    tgt = bottleneck_pairs(b, 2);
    mask = src_clusters == src & tgt_clusters == tgt;
    unique_pairs = unique(lr_pairs(mask));
    bottleneck_diversity(b) = length(unique_pairs);
end
bottleneck_table.UniqueLRPairs = bottleneck_diversity;

if n_low_redundancy > 0
    vulnerable_table = table();
    vulnerable_table.LR_Pair = unique_lr_pairs(low_redundancy_idx);
    vulnerable_table.NumClusterPairs = redundancy_counts(low_redundancy_idx);
    vulnerable_table = sortrows(vulnerable_table, 'NumClusterPairs', 'ascend');
else
    vulnerable_table = table();
end

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

summary = table();
summary.Metric = {'TotalLRPairs'; 'Monopolies'; 'TumorSpecificMonopolies'; ...
                  'LowRedundancySignals'; 'TopBottlenecks'};
summary.Count = [n_lr_pairs; n_monopolies; n_tumor_monopolies; ...
                 n_low_redundancy; n_bottlenecks];
writetable(summary, 'results/novel_discoveries/monopoly_summary.csv');

if n_monopolies > 0
    writetable(monopoly_table, 'results/novel_discoveries/monopolies.csv');
end

writetable(bottleneck_table, 'results/novel_discoveries/bottlenecks.csv');

if n_low_redundancy > 0
    writetable(vulnerable_table, 'results/novel_discoveries/vulnerable_signals.csv');
end

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
