function example()
% EXAMPLE  Minimal end-to-end NEUROSCOPE demo with 500-cell synthetic data.
%   No real data needed — runs every pipeline step and prints a summary.
%
% Usage: example()
% Requires: addpath(genpath('path/to/neuroscope'))

fprintf('=== NEUROSCOPE Example (Synthetic Data) ===\n\n');

% synthetic expression matrices

rng(42);  % reproducibility

n_condition = 500;
n_reference = 300;
n_genes     = 200;

expr_condition = sparse(max(0, round(poissrnd(1, n_condition, n_genes))));
expr_reference = sparse(max(0, round(poissrnd(0.8, n_reference, n_genes))));

upregulated = 1:20;
expr_condition(:, upregulated) = expr_condition(:, upregulated) * 3;

fprintf('Synthetic data generated:\n');
fprintf('  Condition: %d cells × %d genes\n', n_condition, n_genes);
fprintf('  Reference: %d cells × %d genes\n\n', n_reference, n_genes);


gene_symbols = arrayfun(@(i) sprintf('GENE%d', i), 1:n_genes, ...
    'UniformOutput', false)';

ligand_genes  = {'VEGFA', 'EGF', 'TGFB1', 'IL6', 'CXCL12'};
receptor_genes= {'KDR',   'EGFR','TGFBR1','IL6R','CXCR4'};

for k = 1:length(ligand_genes)
    gene_symbols{k}                      = ligand_genes{k};
    gene_symbols{n_genes/2 + k}          = receptor_genes{k};
end

genes_condition = table(gene_symbols, 'VariableNames', {'gene_name'});
genes_reference = genes_condition;


cluster_ids_condition = repmat((1:5)', ceil(n_condition/5), 1);
cluster_ids_condition = cluster_ids_condition(1:n_condition);
clusters_condition = table((1:n_condition)', cluster_ids_condition, ...
    'VariableNames', {'Cell_ID', 'Cluster'});

cluster_ids_reference = repmat((1:5)', ceil(n_reference/5), 1);
cluster_ids_reference = cluster_ids_reference(1:n_reference);
clusters_reference = table((1:n_reference)', cluster_ids_reference, ...
    'VariableNames', {'Cell_ID', 'Cluster'});


lr_pairs = [ligand_genes', receptor_genes'];


fprintf('Step 1: Building communication network...\n');
network_condition = build_communication_network( ...
    expr_condition, genes_condition, clusters_condition, lr_pairs, 'condition');

fprintf('  Interactions found: %d\n\n', length(network_condition.weights));


fprintf('Step 2: Identifying condition-specific interactions...\n');
network_reference = build_communication_network( ...
    expr_reference, genes_reference, clusters_reference, lr_pairs, 'reference');

condition_specific = identify_gbm_specific_interactions( ...
    network_condition, network_reference, expr_condition, expr_reference, ...
    genes_condition, clusters_condition, 500);

fprintf('  Condition-specific pairs: %d\n\n', height(condition_specific));


fprintf('Step 3: Calculating interaction statistics...\n');
stats = calculate_interaction_statistics(network_condition, condition_specific);
fprintf('  Statistics computed for %d interactions\n\n', height(stats));


fprintf('Step 4: Ranking therapeutic targets...\n');
target_rankings = rank_therapeutic_targets(network_condition, condition_specific, stats);
fprintf('  Top target: %s → %s (Score=%.3f)\n\n', ...
    target_rankings.Ligand{1}, target_rankings.Receptor{1}, ...
    target_rankings.TherapeuticScore(1));


fprintf('Step 5: Identifying communication hubs...\n');
hubs = identify_communication_hubs(network_condition);
fprintf('  Ligand hubs: %d | Receptor hubs: %d\n\n', ...
    height(hubs.ligands), height(hubs.receptors));


fprintf('=== Example complete ===\n');
fprintf('For a full run on real data, launch NEUROSCOPE.m or open app1.mlapp.\n\n');

end
