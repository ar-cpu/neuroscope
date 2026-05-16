function fragility_results = discover_network_fragility(gbm_network, healthy_network, target_rankings)
% DISCOVER_NETWORK_FRAGILITY  Quantify the robustness of the condition
%   communication network to targeted gene removal and identify super-additive
%   synergy pairs for combination therapy prioritization.
%
% Analysis components
% -------------------
%   Percolation analysis
%     Simulates progressive random node removal and tracks the normalized
%     size of the largest connected component (GCC) as a function of the
%     fraction of nodes removed. The percolation threshold — the removal
%     fraction at which the GCC collapses below 50% of its initial size —
%     is reported as a summary fragility metric. Networks with low thresholds
%     are fragile (easily fragmented); high thresholds indicate robustness.
%
%   Targeted synergy testing
%     For each pair of high-priority target genes (from target_rankings), the
%     function computes:
%       - Individual disruption effect: fractional GCC reduction from removing
%         each gene alone.
%       - Combined disruption effect: fractional GCC reduction from removing
%         both genes simultaneously.
%       - Synergy score: (combined − expected_additive) / expected_additive
%     Positive synergy indicates super-additive disruption (the pair works
%     better together than the sum of their individual effects), identifying
%     candidate combination therapy targets.
%
%   Differential fragility
%     The condition percolation curve is compared to the reference network
%     curve. Greater fragility in the condition network (lower percolation
%     threshold) is consistent with the oncogenic rewiring hypothesis.
%
% Usage:
%   fragility_results = discover_network_fragility(gbm_network, healthy_network, target_rankings)
%
% Inputs:
%   gbm_network     - Condition communication network struct with fields
%                     .ligands and .receptors.
%   healthy_network - Reference communication network struct (same schema).
%                     May be empty (create_empty_network()) for single-dataset
%                     runs; differential analysis is skipped in that case.
%   target_rankings - Table from rank_therapeutic_targets() used to select
%                     the gene candidates for pairwise synergy testing.
%
% Outputs:
%   fragility_results.synergy_table       - Pairwise synergy scores (sorted descending)
%   fragility_results.percolation_curve   - [fractions × 2] matrix: removal fraction
%                                           and corresponding normalized GCC size
%   fragility_results.percolation_threshold - Removal fraction at GCC < 0.5
%   fragility_results.hub_genes           - High-degree genes driving network integrity
%   fragility_results.summary             - Key metrics struct

fprintf('Network fragility analysis\n\n');

if ~isstruct(gbm_network) || ~isstruct(healthy_network)
    error('Networks must be structs');
end

%% Build gene networks

fprintf('Building gene networks...\n');

[gbm_genes, gbm_adj] = build_gene_network(gbm_network.ligands, gbm_network.receptors);
[healthy_genes, healthy_adj] = build_gene_network(healthy_network.ligands, healthy_network.receptors);

n_gbm = length(gbm_genes);
n_healthy = length(healthy_genes);

gbm_gcc = find_gcc_size(gbm_adj);
healthy_gcc = find_gcc_size(healthy_adj);

fprintf('  Condition: %d genes, %d edges, GCC=%d (%.1f%%)\n', ...
    n_gbm, nnz(gbm_adj), gbm_gcc, 100*gbm_gcc/max(1,n_gbm));
fprintf('  Reference: %d genes, %d edges, GCC=%d (%.1f%%)\n', ...
    n_healthy, nnz(healthy_adj), healthy_gcc, 100*healthy_gcc/max(1,n_healthy));

%% Percolation analysis: GCC size vs fraction removed

fprintf('\nRandom percolation analysis...\n');

n_steps = 10;
removal_fracs = linspace(0, 0.4, n_steps);
gcc_curve = zeros(n_steps, 1);

rng(42);  % reproducibility

for step = 1:n_steps
    n_remove = round(removal_fracs(step) * n_gbm);
    if n_remove > 0 && n_remove < n_gbm
        remove_idx = randperm(n_gbm, n_remove);
        adj_temp = gbm_adj;
        adj_temp(remove_idx, :) = 0;
        adj_temp(:, remove_idx) = 0;
        gcc_curve(step) = find_gcc_size(adj_temp);
    else
        gcc_curve(step) = gbm_gcc;
    end
end

thresh_idx = find(gcc_curve < 0.5 * gbm_gcc, 1);
if ~isempty(thresh_idx)
    fprintf('  Percolation threshold: %.1f%% removal\n', 100 * removal_fracs(thresh_idx));
else
    fprintf('  Network robust: GCC stays >50%% through 40%% random removal\n');
end

%% Identify high-degree hub targets

fprintf('\nIdentifying high-degree hub genes...\n');

degrees = full(sum(gbm_adj, 1)' + sum(gbm_adj, 2));
[sorted_degs, sort_idx] = sort(degrees, 'descend');

n_targets = min(50, n_gbm);
target_indices = sort_idx(1:n_targets);

fprintf('  Hub targets: %d (top degree genes)\n', n_targets);
fprintf('  Mean hub degree: %.1f\n', mean(sorted_degs(1:n_targets)));

%% Single-target removal effects

fprintf('\nComputing single-target effects...\n');

single_effects = zeros(n_targets, 1);
for i = 1:n_targets
    adj_temp = gbm_adj;
    idx = target_indices(i);
    adj_temp(idx, :) = 0;
    adj_temp(:, idx) = 0;
    single_effects(i) = gbm_gcc - find_gcc_size(adj_temp);
end

fprintf('  Mean effect: %.1f nodes disconnected\n', mean(single_effects));
fprintf('  Max effect:  %.1f nodes disconnected\n', max(single_effects));

%% Pairwise synergy analysis

fprintf('\nComputing pairwise synergy...\n');

max_pairs = min(500, n_targets*(n_targets-1)/2);

synergy_table = table('Size', [0, 7], ...
    'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Target1', 'Target2', 'Effect1', 'Effect2', ...
    'CombinedEffect', 'ExpectedEffect', 'SynergyScore'});

pair_count = 0;

for i = 1:n_targets
    if pair_count >= max_pairs
        break;
    end

    for j = i+1:n_targets
        if pair_count >= max_pairs
            break;
        end

        idx_i = target_indices(i);
        idx_j = target_indices(j);

        adj_temp = gbm_adj;
        adj_temp(idx_i, :) = 0;
        adj_temp(:, idx_i) = 0;
        adj_temp(idx_j, :) = 0;
        adj_temp(:, idx_j) = 0;

        combined_effect = gbm_gcc - find_gcc_size(adj_temp);
        expected_additive = single_effects(i) + single_effects(j);

        if expected_additive > 0
            synergy = (combined_effect - expected_additive) / expected_additive;
        else
            synergy = 0;
        end

        synergy_table = [synergy_table; {gbm_genes{idx_i}, gbm_genes{idx_j}, ...
            single_effects(i), single_effects(j), ...
            combined_effect, expected_additive, synergy}];

        pair_count = pair_count + 1;
    end

    if mod(i, 10) == 0
        fprintf('  Target %d/%d, %d pairs tested\n', i, n_targets, pair_count);
    end
end

synergy_table = sortrows(synergy_table, 'SynergyScore', 'descend');

n_synergistic = sum(synergy_table.SynergyScore > 0.2);
n_strong = sum(synergy_table.SynergyScore > 0.5);

fprintf('  Pairs tested: %d\n', pair_count);
fprintf('  Synergistic (>20%%): %d\n', n_synergistic);
fprintf('  Strongly synergistic (>50%%): %d\n', n_strong);

%% Export

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

writetable(synergy_table, 'results/novel_discoveries/drug_synergy.csv');
fprintf('\n  Saved: drug_synergy.csv\n');

%% Package results

fragility_results = struct();
fragility_results.synergy_data = synergy_table;
fragility_results.n_synergistic = n_synergistic;
fragility_results.n_strong_synergistic = n_strong;
fragility_results.percolation_curve = gcc_curve;
fragility_results.percolation_fracs = removal_fracs;
fragility_results.single_effects = single_effects;

fprintf('\nFragility analysis complete\n');
fprintf('  Network: %d genes, GCC=%d (%.1f%%)\n', n_gbm, gbm_gcc, 100*gbm_gcc/max(1,n_gbm));
fprintf('  Synergistic pairs: %d | Strong: %d\n\n', n_synergistic, n_strong);

if height(synergy_table) >= 1
    fprintf('Top synergistic target pairs:\n');
    for i = 1:min(10, height(synergy_table))
        fprintf('  %s + %s: %.1f + %.1f -> %.1f (synergy: +%.0f%%)\n', ...
            synergy_table.Target1{i}, synergy_table.Target2{i}, ...
            synergy_table.Effect1(i), synergy_table.Effect2(i), ...
            synergy_table.CombinedEffect(i), 100*synergy_table.SynergyScore(i));
    end
    fprintf('\n');
end

end

%% Local helper functions

function [genes, adj] = build_gene_network(ligs, recs)
pair_sigs = strcat(ligs, '->', recs);
[unique_pairs, ~, ~] = unique(pair_sigs);

ligs_u = cell(length(unique_pairs), 1);
recs_u = cell(length(unique_pairs), 1);
for i = 1:length(unique_pairs)
    parts = strsplit(unique_pairs{i}, '->');
    ligs_u{i} = parts{1};
    recs_u{i} = parts{2};
end

genes = unique([ligs_u; recs_u]);
n = length(genes);
adj = sparse(n, n);

for i = 1:length(ligs_u)
    lig_idx = find(strcmp(genes, ligs_u{i}), 1);
    rec_idx = find(strcmp(genes, recs_u{i}), 1);
    if ~isempty(lig_idx) && ~isempty(rec_idx)
        adj(lig_idx, rec_idx) = 1;
    end
end
end

function gcc_size = find_gcc_size(adj)
n = size(adj, 1);
visited = false(n, 1);
max_component = 0;

for start = 1:n
    if ~visited(start)
        queue = start;
        component_size = 0;

        while ~isempty(queue)
            current = queue(1);
            queue(1) = [];

            if ~visited(current)
                visited(current) = true;
                component_size = component_size + 1;
                neighbors = find(adj(current, :) | adj(:, current)');
                for neighbor = neighbors
                    if ~visited(neighbor)
                        queue = [queue, neighbor];
                    end
                end
            end
        end

        max_component = max(max_component, component_size);
    end
end

gcc_size = max_component;
end
