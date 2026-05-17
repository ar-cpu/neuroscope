function motif_results = discover_network_motifs(gbm_network, healthy_network, condition_specific)
% DISCOVER_NETWORK_MOTIFS  3-node subgraph motif analysis (Chain/Loop/Fan)
%   with Maslov-Sneppen z-scores and edge essentiality via perturbation.
%
% motif_results = discover_network_motifs(gbm_network, healthy_network, condition_specific)
%
% Outputs:
%   .perturbation_table - Per-motif-edge: EssentialityScore, MotifTypeZScore
%   .motif_counts       - [Chain, Loop, Fan] observed counts
%   .motif_zscores      - [Chain, Loop, Fan] z-scores vs degree-preserving nulls

fprintf('Network motif analysis\n\n');

if ~isstruct(gbm_network) || ~isstruct(healthy_network)
    error('Networks must be structs');
end

%% Build unique-pair networks

fprintf('Building gene networks from unique L-R pairs...\n');

[gbm_genes, gbm_adj] = build_gene_network(gbm_network.ligands, gbm_network.receptors);
[~, healthy_adj_tmp] = build_gene_network(healthy_network.ligands, healthy_network.receptors);

n_gbm = length(gbm_genes);
n_healthy_genes = size(healthy_adj_tmp, 1);

fprintf('  Condition network: %d genes, %d edges (mean degree: %.2f)\n', ...
    n_gbm, nnz(gbm_adj), nnz(gbm_adj)/max(1, n_gbm));
fprintf('  Reference network: %d genes, %d edges\n', ...
    n_healthy_genes, nnz(healthy_adj_tmp));

if nnz(gbm_adj) / max(1, n_gbm) < 2
    fprintf('  Note: Sparse network (mean degree < 2).\n');
end

%% Find 3-node motifs

fprintf('\nSearching for 3-node motifs...\n');

max_search = min(300, n_gbm);
motifs = find_motifs(gbm_adj, max_search);

fprintf('  Motifs found: %d\n', length(motifs));

if isempty(motifs)
    warning('No motifs found - network too sparse for motif analysis');
    motif_results = make_empty_result();
    return;
end

%% Motif z-scores vs Maslov-Sneppen degree-preserving null model
%
% Tests whether each motif type (Chain, Loop, Fan) is statistically
% over- or under-represented relative to random networks with the same
% in/out degree sequence. Rewiring randomly swaps edge endpoints while
% preserving all node degrees, generating a configuration-model null.
%
% Reference: Milo et al. (2002) Science 298(5594):824-827.
%            Maslov & Sneppen (2002) Science 296(5569):910-913.

fprintf('\nComputing motif z-scores vs Maslov-Sneppen null model...\n');

motif_type_labels = {'Chain', 'Loop', 'Fan'};
n_types = length(motif_type_labels);

observed_counts = zeros(1, n_types);
for t = 1:n_types
    observed_counts(t) = sum(cellfun(@(m) strcmp(m.type, motif_type_labels{t}), motifs));
end
fprintf('  Observed: Chain=%d, Loop=%d, Fan=%d\n', ...
    observed_counts(1), observed_counts(2), observed_counts(3));

K_rand = 10;  % number of rewired null networks
rand_counts = zeros(K_rand, n_types);

for k = 1:K_rand
    gbm_adj_rand = maslov_sneppen_rewire(gbm_adj, 10 * nnz(gbm_adj));
    motifs_rand = find_motifs(gbm_adj_rand, max_search);
    for t = 1:n_types
        rand_counts(k, t) = sum(cellfun(@(m) strcmp(m.type, motif_type_labels{t}), motifs_rand));
    end
    if mod(k, 5) == 0
        fprintf('  Null network %d/%d\n', k, K_rand);
    end
end

motif_zscore_map = containers.Map();
for t = 1:n_types
    mu_rand = mean(rand_counts(:, t));
    sd_rand = std(rand_counts(:, t));
    if sd_rand > 0
        z = (observed_counts(t) - mu_rand) / sd_rand;
    else
        z = 0;
    end
    motif_zscore_map(motif_type_labels{t}) = z;
    sig_flag = '';
    if abs(z) > 1.96, sig_flag = ' *'; end
    fprintf('  %s: z=%.2f (obs=%d, rand=%.1f±%.1f)%s\n', ...
        motif_type_labels{t}, z, observed_counts(t), mu_rand, sd_rand, sig_flag);
end
fprintf('\n');

%% Perturbation analysis

fprintf('\nRunning edge perturbation analysis...\n');

n_test = min(100, length(motifs));

pert_table = table('Size', [0, 9], ...
    'VariableTypes', {'double', 'string', 'string', 'string', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'MotifID', 'MotifType', 'EdgeRemoved', 'Nodes', ...
    'BaselineReach', 'PerturbedReach', 'CascadeEffect', 'Essentiality', 'MotifTypeZScore'});

for m = 1:n_test
    motif_nodes = motifs{m}.nodes;
    motif_type = motifs{m}.type;

    % Collect all directed edges within the motif
    edges = [];
    for ii = 1:length(motif_nodes)
        for jj = 1:length(motif_nodes)
            if ii ~= jj && gbm_adj(motif_nodes(ii), motif_nodes(jj)) > 0
                edges = [edges; motif_nodes(ii), motif_nodes(jj)];
            end
        end
    end

    if isempty(edges), continue; end

    baseline = compute_reach_bfs(gbm_adj, motif_nodes);

    for e = 1:size(edges, 1)
        adj_pert = gbm_adj;
        adj_pert(edges(e,1), edges(e,2)) = 0;

        perturbed = compute_reach_bfs(adj_pert, motif_nodes);
        cascade = baseline - perturbed;
        essentiality = cascade / max(baseline, 1);

        edge_name = sprintf('%s->%s', gbm_genes{edges(e,1)}, gbm_genes{edges(e,2)});
        node_str = strjoin(gbm_genes(motif_nodes), ',');

        type_z = 0;
        if isKey(motif_zscore_map, motif_type)
            type_z = motif_zscore_map(motif_type);
        end

        pert_table = [pert_table; {m, motif_type, edge_name, node_str, ...
            baseline, perturbed, cascade, essentiality, type_z}];
    end

    if mod(m, 20) == 0
        fprintf('  Motif %d/%d\n', m, n_test);
    end
end

pert_table = sortrows(pert_table, 'Essentiality', 'descend');

n_critical = sum(pert_table.Essentiality > 0.3);
fprintf('  Total perturbations tested: %d\n', height(pert_table));
fprintf('  Critical edges (essentiality > 0.3): %d\n', n_critical);

%% Flag condition-specific edges

fprintf('\nFlagging condition-specific edges...\n');

pert_table.ConditionSpecific = flag_condition_specific_edges(pert_table, condition_specific);
n_tumor = sum(pert_table.ConditionSpecific & pert_table.Essentiality > 0.3);
fprintf('  Condition-specific critical edges: %d\n', n_tumor);

% Keep backward-compatible field name
pert_table.TumorSpecific = pert_table.ConditionSpecific;

%% Export

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

writetable(pert_table, 'results/novel_discoveries/motif_perturbation.csv');
fprintf('\n  Saved: motif_perturbation.csv\n');

%% Package results

motif_results = struct();
motif_results.perturbation_data = pert_table;
motif_results.n_critical = n_critical;
motif_results.n_tumor_critical = n_tumor;
motif_results.motif_zscores = motif_zscore_map;
motif_results.observed_motif_counts = observed_counts;

fprintf('\nMotif analysis complete\n');
fprintf('  Motifs: %d | Critical edges: %d | Condition-specific: %d\n\n', ...
    length(motifs), n_critical, n_tumor);

if height(pert_table) >= 1
    fprintf('Top critical edges:\n');
    for i = 1:min(5, height(pert_table))
        fprintf('  %s (%s, essentiality=%.3f)', ...
            pert_table.EdgeRemoved{i}, pert_table.MotifType{i}, pert_table.Essentiality(i));
        if pert_table.ConditionSpecific(i)
            fprintf(' [condition-specific]');
        end
        fprintf('\n');
    end
    fprintf('\n');
end

end

%% Local helper functions

function [genes, adj] = build_gene_network(ligs, recs)
% Build sparse adjacency from lists, taking only unique pairs
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

function motifs = find_motifs(adj, max_search)
% Find all connected 3-node subgraphs and classify them
n = min(size(adj, 1), max_search);
motifs = {};

for i = 1:n
    neighbors_i = find(adj(i,:) | adj(:,i)');

    for j_idx = 1:length(neighbors_i)
        j = neighbors_i(j_idx);
        if j <= i, continue; end

        neighbors_j = find(adj(j,:) | adj(:,j)');
        common = intersect(neighbors_i, neighbors_j);

        for k_idx = 1:length(common)
            k = common(k_idx);
            if k <= j, continue; end

            motif = struct();
            motif.nodes = [i, j, k];

            sub = full(adj(motif.nodes, motif.nodes));
            n_edges = nnz(sub);

            if n_edges >= 2
                if n_edges == 2
                    motif.type = 'Chain';
                elseif n_edges == 3
                    motif.type = 'Loop';
                else
                    motif.type = 'Fan';
                end
                motifs{end+1} = motif;
            end

            if length(motifs) >= 300
                return;
            end
        end
    end
end
end

function reach = compute_reach_bfs(adj, sources)
n = size(adj, 1);
visited = false(n, 1);
visited(sources) = true;

queue = sources(:)';
head = 1;

while head <= length(queue)
    curr = queue(head);
    head = head + 1;
    neighbors = find(adj(curr, :) & ~visited');
    if ~isempty(neighbors)
        visited(neighbors) = true;
        queue = [queue, neighbors];
    end
end

reach = sum(visited);
end

function is_specific = flag_condition_specific_edges(pert_table, condition_specific)
n = height(pert_table);
is_specific = false(n, 1);

if ~istable(condition_specific) || height(condition_specific) == 0
    return;
end

cols = condition_specific.Properties.VariableNames;
lig_col = '';
rec_col = '';
for i = 1:length(cols)
    if contains(lower(cols{i}), {'ligand', 'source'})
        lig_col = cols{i};
    end
    if contains(lower(cols{i}), {'receptor', 'target'})
        rec_col = cols{i};
    end
end

if isempty(lig_col) || isempty(rec_col)
    return;
end

for row = 1:n
    parts = strsplit(pert_table.EdgeRemoved{row}, '->');
    if length(parts) == 2
        for spec_row = 1:height(condition_specific)
            if strcmp(condition_specific.(lig_col){spec_row}, parts{1}) && ...
               strcmp(condition_specific.(rec_col){spec_row}, parts{2})
                is_specific(row) = true;
                break;
            end
        end
    end
end
end

function adj_rand = maslov_sneppen_rewire(adj, n_target_swaps)
% Maslov-Sneppen edge rewiring: generates a random network with the same
% in-degree and out-degree sequence as the input by repeatedly swapping
% edge endpoints. Each attempted swap picks two edges (a→b) and (c→d)
% and rewires to (a→d) and (c→b), accepting only if no self-loops or
% duplicate edges would be created.
%
% Reference: Maslov & Sneppen (2002) Science 296(5569):910-913.

[rows, cols] = find(adj);
n_edges = length(rows);

if n_edges < 2
    adj_rand = adj;
    return;
end

adj_rand = adj;
swaps_done = 0;
max_tries = 20 * n_target_swaps;
tries = 0;

while swaps_done < n_target_swaps && tries < max_tries
    tries = tries + 1;

    e1 = randi(n_edges);
    e2 = randi(n_edges);
    if e1 == e2, continue; end

    a = rows(e1);  b = cols(e1);
    c = rows(e2);  d = cols(e2);

    % Skip if any node appears twice (would create self-loop or same-node issue)
    if a == c || a == d || b == c || b == d, continue; end

    % New edges would be a→d and c→b — skip if either already exists
    if adj_rand(a, d) > 0 || adj_rand(c, b) > 0, continue; end

    % Perform swap
    adj_rand(a, d) = adj_rand(a, b);
    adj_rand(c, b) = adj_rand(c, d);
    adj_rand(a, b) = 0;
    adj_rand(c, d) = 0;

    rows(e1) = a;  cols(e1) = d;
    rows(e2) = c;  cols(e2) = b;

    swaps_done = swaps_done + 1;
end
end

function result = make_empty_result()
result = struct();
result.perturbation_data = table();
result.n_critical = 0;
result.n_tumor_critical = 0;
result.motif_zscores = containers.Map();
result.observed_motif_counts = zeros(1, 3);
end
