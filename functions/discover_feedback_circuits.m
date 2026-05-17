function circuit_results = discover_feedback_circuits(gbm_network, condition_specific, max_circuit_length)
% DISCOVER_FEEDBACK_CIRCUITS  Enumerate directed feedback cycles and classify
%   by loop gain G = ∏ w_i (Amplifying G>1.05 / Near-critical / Attenuating).
%
% circuit_results = discover_feedback_circuits(gbm_network, condition_specific)
% circuit_results = discover_feedback_circuits(gbm_network, condition_specific, max_circuit_length)
%
% Inputs:
%   gbm_network          - Network struct (.ligands, .receptors, .weights)
%   condition_specific   - Table of condition-specific interactions. Used to
%                          annotate whether each circuit's edges are enriched
%                          in the condition vs. reference.
%   max_circuit_length   - Maximum number of nodes in a detected cycle.
%                          Default: 4. Longer cycles grow exponentially in
%                          number; values above 5 are not recommended.
%
% Outputs:
%   circuit_results.circuits       - Table of detected cycles with columns:
%                                    CircuitID, GeneA, GeneB, GeneC, GeneD,
%                                    Length, Strength, LoopGain, CircuitType,
%                                    IsConditionSpecific
%   circuit_results.n_circuits     - Total circuits found
%   circuit_results.n_amplifying   - Amplifying circuits (gain > 1.05)
%   circuit_results.n_attenuating  - Attenuating circuits (gain < 0.95)

if nargin < 3
    max_circuit_length = 4;
end

fprintf('Feedback circuit analysis\n\n');

if ~isstruct(gbm_network)
    error('Network must be a struct');
end

ligs = gbm_network.ligands;
recs = gbm_network.receptors;
weights = gbm_network.weights;

fprintf('  Total interactions: %d\n', length(ligs));

%% Build weighted gene network from unique L-R pairs

fprintf('\nBuilding gene network from unique L-R pairs...\n');

pair_sigs = strcat(ligs, '->', recs);
[unique_pairs, ~, pair_idx] = unique(pair_sigs);

% For each unique pair, compute mean weight
ligs_unique = cell(length(unique_pairs), 1);
recs_unique = cell(length(unique_pairs), 1);
weights_unique = zeros(length(unique_pairs), 1);

for i = 1:length(unique_pairs)
    parts = strsplit(unique_pairs{i}, '->');
    ligs_unique{i} = parts{1};
    recs_unique{i} = parts{2};
    weights_unique(i) = mean(weights(pair_idx == i));
end

fprintf('  Unique L-R pairs: %d\n', length(unique_pairs));

all_genes = unique([ligs_unique; recs_unique]);
n_genes = length(all_genes);

% Build weighted adjacency matrix
adj = sparse(n_genes, n_genes);
for i = 1:length(ligs_unique)
    lig_idx = find(strcmp(all_genes, ligs_unique{i}), 1);
    rec_idx = find(strcmp(all_genes, recs_unique{i}), 1);
    if ~isempty(lig_idx) && ~isempty(rec_idx)
        % Keep max weight if pair appears multiple times (shouldn't after unique())
        adj(lig_idx, rec_idx) = max(adj(lig_idx, rec_idx), weights_unique(i));
    end
end

n_edges = nnz(adj);
mean_deg = n_edges / n_genes;
fprintf('  Genes: %d | Edges: %d | Mean degree: %.2f\n', n_genes, n_edges, mean_deg);

if mean_deg < 2
    fprintf('  Note: Sparse network (mean degree < 2). Fewer circuits expected.\n');
end

%% Adjacency list for DFS

adj_list = cell(n_genes, 1);
for i = 1:n_genes
    adj_list{i} = find(adj(i, :));
end

%% Search for directed cycles using DFS

fprintf('\nSearching for feedback circuits (max length: %d)...\n', max_circuit_length);

max_circuits = 200;
max_start_nodes = min(100, n_genes);

all_circuits = {};
n_found = 0;

for start = 1:max_start_nodes
    if n_found >= max_circuits
        break;
    end

    found = find_cycles_dfs(start, start, adj_list, [], max_circuit_length, max_circuits - n_found);

    for i = 1:length(found)
        n_found = n_found + 1;
        all_circuits{n_found} = found{i};
    end

    if mod(start, 20) == 0 && n_found > 0
        fprintf('  Searched %d/%d genes, %d circuits found\n', start, max_start_nodes, n_found);
    end
end

fprintf('  Total circuits found: %d\n\n', n_found);

if n_found == 0
    fprintf('  No feedback circuits found.\n');
    fprintf('  The network appears to be feedforward-dominant (acyclic or near-acyclic).\n');
    fprintf('  This is a meaningful biological finding - feedforward architectures\n');
    fprintf('  implement signal integration without runaway feedback amplification.\n\n');
    circuit_results = make_empty_result();
    return;
end

%% Analyze circuits

fprintf('Analyzing circuit properties...\n');

circuit_table = table('Size', [n_found, 8], ...
    'VariableTypes', {'double', 'double', 'string', 'double', 'double', 'double', 'string', 'logical'}, ...
    'VariableNames', {'CircuitID', 'Length', 'Path', 'Strength', 'LoopGain', 'Robustness', 'CircuitType', 'ConditionSpecific'});

for i = 1:n_found
    circuit_nodes = all_circuits{i};
    n_nodes_in_circuit = length(circuit_nodes);

    gene_name_list = cell(1, n_nodes_in_circuit);
    for j = 1:n_nodes_in_circuit
        gene_name_list{j} = all_genes{circuit_nodes(j)};
    end

    path_str = strjoin(gene_name_list, ' -> ');

    % Compute circuit strength as geometric mean of edge weights along the loop
    edge_weights = zeros(n_nodes_in_circuit, 1);
    for j = 1:n_nodes_in_circuit
        src_node = circuit_nodes(j);
        tgt_node = circuit_nodes(mod(j, n_nodes_in_circuit) + 1);
        w = adj(src_node, tgt_node);
        edge_weights(j) = w;
    end

    if all(edge_weights > 0)
        strength = exp(mean(log(edge_weights)));  % geometric mean of weights
    else
        strength = 0;
    end

    % Loop gain: product of all edge weights around the circuit.
    % Drawn from control-theory analysis of gene regulatory networks:
    % gain > 1 indicates the circuit can self-amplify perturbations
    % (positive feedback, bistability), gain < 1 indicates damped negative
    % feedback (homeostatic), and gain ≈ 1 is near the stability boundary.
    %
    % Reference: Thomas & d'Ari (1990) Biological Feedback;
    %            Cinquin & Demongeot (2002) J Theor Biol 216(2):229-241.
    if all(edge_weights > 0)
        loop_gain = prod(edge_weights);
    else
        loop_gain = 0;
    end

    if loop_gain > 1.05
        circuit_type = 'Amplifying';     % may sustain oscillations or bistability
    elseif loop_gain > 0 && loop_gain < 0.95
        circuit_type = 'Attenuating';    % negative feedback; dampens perturbations
    elseif loop_gain > 0
        circuit_type = 'Near-critical';  % marginal stability
    else
        circuit_type = 'Broken';         % zero-weight edge severs the loop
    end

    % Robustness: shorter circuits are more robust (fewer edges that can be broken)
    robustness = 1.0 / n_nodes_in_circuit;

    is_condition_specific = is_circuit_condition_specific(gene_name_list, condition_specific);

    circuit_table.CircuitID(i) = i;
    circuit_table.Length(i) = n_nodes_in_circuit;
    circuit_table.Path{i} = path_str;
    circuit_table.Strength(i) = strength;
    circuit_table.LoopGain(i) = loop_gain;
    circuit_table.Robustness(i) = robustness;
    circuit_table.CircuitType{i} = circuit_type;
    circuit_table.ConditionSpecific(i) = is_condition_specific;
end

circuit_table = sortrows(circuit_table, 'Strength', 'descend');

n_amplifying   = sum(strcmp(circuit_table.CircuitType, 'Amplifying'));
n_attenuating  = sum(strcmp(circuit_table.CircuitType, 'Attenuating'));
n_near_crit    = sum(strcmp(circuit_table.CircuitType, 'Near-critical'));

fprintf('  High-strength circuits (>0.1): %d\n', sum(circuit_table.Strength > 0.1));
fprintf('  Amplifying circuits (gain>1):   %d\n', n_amplifying);
fprintf('  Attenuating circuits (gain<1):  %d\n', n_attenuating);
fprintf('  Near-critical (gain≈1):         %d\n', n_near_crit);
fprintf('  Condition-specific circuits:    %d\n\n', sum(circuit_table.ConditionSpecific));

%% Export

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

writetable(circuit_table, 'results/novel_discoveries/feedback_circuits.csv');
fprintf('  Saved: feedback_circuits.csv\n\n');

%% Package results

circuit_results = struct();
circuit_results.circuits = circuit_table;
circuit_results.n_circuits = n_found;
circuit_results.n_strong = sum(circuit_table.Strength > 0.1);
circuit_results.n_amplifying = sum(strcmp(circuit_table.CircuitType, 'Amplifying'));
circuit_results.n_attenuating = sum(strcmp(circuit_table.CircuitType, 'Attenuating'));
circuit_results.n_condition_specific = sum(circuit_table.ConditionSpecific);

% Backward compatibility
circuit_results.n_robust = sum(circuit_table.Robustness > 0.5);
circuit_results.n_tumor_specific = circuit_results.n_condition_specific;

fprintf('Feedback circuit analysis complete\n');
fprintf('  Circuits: %d | Strong: %d | Condition-specific: %d\n\n', ...
    n_found, circuit_results.n_strong, circuit_results.n_condition_specific);

if height(circuit_table) >= 1
    fprintf('Top circuits by strength:\n');
    for i = 1:min(5, height(circuit_table))
        fprintf('  %d. %s\n', i, circuit_table.Path{i});
        fprintf('     Length: %d | Strength: %.4f | Gain: %.4f [%s]', ...
            circuit_table.Length(i), circuit_table.Strength(i), ...
            circuit_table.LoopGain(i), circuit_table.CircuitType{i});
        if circuit_table.ConditionSpecific(i)
            fprintf(' [condition-specific]');
        end
        fprintf('\n');
    end
    fprintf('\n');
end

end

%% Local helper functions

function circuits = find_cycles_dfs(current, target, adj_list, path, max_len, max_find)
circuits = {};

if length(path) >= max_len
    return;
end

if length(path) > 1 && current == target
    circuits{1} = path;
    return;
end

if any(path == current)
    return;
end

path = [path, current];
neighbors = adj_list{current};

for i = 1:length(neighbors)
    if length(circuits) >= max_find
        return;
    end
    found = find_cycles_dfs(neighbors(i), target, adj_list, path, max_len, max_find - length(circuits));
    circuits = [circuits, found];
end
end

function is_specific = is_circuit_condition_specific(gene_names, condition_specific)
is_specific = false;

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

n = length(gene_names);
for j = 1:n
    next_j = mod(j, n) + 1;
    for row = 1:height(condition_specific)
        if strcmp(condition_specific.(lig_col){row}, gene_names{j}) && ...
           strcmp(condition_specific.(rec_col){row}, gene_names{next_j})
            is_specific = true;
            return;
        end
    end
end
end

function result = make_empty_result()
result = struct();
result.circuits = table();
result.n_circuits = 0;
result.n_strong = 0;
result.n_robust = 0;
result.n_condition_specific = 0;
result.n_tumor_specific = 0;
end
