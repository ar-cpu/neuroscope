function role_results = discover_communication_role_asymmetry(gbm_network, gbm_clusters, gbm_specific)
% DISCOVER_COMMUNICATION_ROLE_ASYMMETRY  Classify cell clusters as signal
%   broadcasters, listeners, or balanced integrators based on their
%   out-degree versus in-degree in the directed ligand-receptor network.
%
% Conceptual basis
% ----------------
% In a directed cell-cell communication network each interaction edge flows
% from a source (ligand-expressing) cluster to a target (receptor-expressing)
% cluster. A cluster's communication role is captured by the signed asymmetry
% of its degree:
%
%   RoleIndex = (outgoing − incoming) / (outgoing + incoming)  ∈ [−1, +1]
%
% Clusters with RoleIndex > +0.3 are classified as "broadcasters" (primary
% signal senders); those with RoleIndex < −0.3 are "listeners" (primary
% signal receivers); the remainder are "balanced" integrators. The 0.3
% threshold is a standard moderate-effect heuristic, not statistically
% calibrated — treat it as exploratory.
%
% The function additionally characterizes the dominant signaling vocabulary
% of top broadcasters and the source diversity received by top listeners,
% then cross-references results against the condition-specific interaction
% table to flag biologically relevant channels.
%
% Usage:
%   role_results = discover_communication_role_asymmetry(gbm_network, gbm_clusters, gbm_specific)
%
% Inputs:
%   gbm_network  - Communication network struct with fields:
%                    .source_clusters - [N×1] cluster ID of the sending cell
%                    .target_clusters - [N×1] cluster ID of the receiving cell
%                    .ligands         - {N×1} ligand gene symbols
%                    .receptors       - {N×1} receptor gene symbols
%   gbm_clusters - Table with a column whose name contains 'cluster',
%                  used to enumerate all valid cluster IDs. If not a table,
%                  unique IDs are inferred directly from the network.
%   gbm_specific - Table of condition-specific interactions from
%                  identify_gbm_specific_interactions(). Used to flag
%                  broadcaster signals that are condition-enriched.
%
% Outputs:
%   role_results.role_table           - Per-cluster role summary table
%   role_results.broadcaster_signals  - Top broadcaster clusters with their
%                                       dominant L-R pair and signal counts
%   role_results.listener_integration - Top listener clusters with source
%                                       diversity and total incoming signals
%   role_results.n_broadcasters       - Number of broadcaster clusters
%   role_results.n_listeners          - Number of listener clusters
%   role_results.n_balanced           - Number of balanced clusters
%   role_results.n_tumor_broadcasters - Broadcasters with condition-specific
%                                       dominant signals

fprintf('Communication Role Asymmetry Analysis\n\n');

fprintf('ANALYSIS:\n');
fprintf('  Identifying "Broadcaster" vs "Listener" clusters\n');
fprintf('  in directed cell-cell communication networks\n\n');

%% Phase 1: Input validation

if ~isstruct(gbm_network)
    error('gbm_network must be struct');
end

if ~isfield(gbm_network, 'source_clusters') || ~isfield(gbm_network, 'target_clusters')
    error('Network missing cluster information');
end

src_clusters = gbm_network.source_clusters;
tgt_clusters = gbm_network.target_clusters;

fprintf('  Total interactions: %d\n', length(src_clusters));
fprintf('  Directed network: source -> target\n');

% Get cluster IDs with explicit column vector conversion
if istable(gbm_clusters)
    cluster_col = '';
    cols = gbm_clusters.Properties.VariableNames;
    % Prefer the standardised 'Cluster' column before partial-name search
    if any(strcmp(cols, 'Cluster'))
        cluster_col = 'Cluster';
    else
        for i = 1:length(cols)
            if contains(lower(cols{i}), 'cluster')
                cluster_col = cols{i};
                break;
            end
        end
    end
    if isempty(cluster_col)
        cluster_col = cols{1};
    end
    all_cluster_ids = unique(gbm_clusters.(cluster_col));
else
    all_cluster_ids = unique([src_clusters; tgt_clusters]);
end

% Ensure numeric column vector — fall back to network cluster IDs if the
% selected column is non-numeric (e.g., a string supercluster_term column).
if iscell(all_cluster_ids) || ischar(all_cluster_ids)
    all_cluster_ids = unique([src_clusters; tgt_clusters]);
end
all_cluster_ids = double(all_cluster_ids(:));

n_clusters = length(all_cluster_ids);
fprintf('  Total clusters: %d\n', n_clusters);
fprintf('\n');

%% Phase 2: Out-degree and in-degree per cluster

outgoing_counts = zeros(n_clusters, 1);
incoming_counts = zeros(n_clusters, 1);

for c = 1:n_clusters
    cluster_id = all_cluster_ids(c);
    
    % Outgoing: This cluster is SOURCE
    outgoing_counts(c) = sum(src_clusters == cluster_id);
    
    % Incoming: This cluster is TARGET
    incoming_counts(c) = sum(tgt_clusters == cluster_id);
end

fprintf('  Mean outgoing per cluster: %.1f\n', mean(outgoing_counts));
fprintf('  Mean incoming per cluster: %.1f\n', mean(incoming_counts));
fprintf('  Max outgoing: %d\n', max(outgoing_counts));
fprintf('  Max incoming: %d\n', max(incoming_counts));
fprintf('\n');

%% Phase 3: Role index and cluster classification

% Role Index = (out - in) / (out + in)
role_index = zeros(n_clusters, 1);
total_signals = outgoing_counts + incoming_counts;

for c = 1:n_clusters
    if total_signals(c) > 0
        role_index(c) = (outgoing_counts(c) - incoming_counts(c)) / total_signals(c);
    else
        role_index(c) = 0;
    end
end

fprintf('  Role Index range: [%.3f, %.3f]\n', min(role_index), max(role_index));
fprintf('  Mean: %.3f\n', mean(role_index));
fprintf('  Std: %.3f\n', std(role_index));

% Classify roles
is_broadcaster = role_index > 0.3;
is_listener = role_index < -0.3;
is_balanced = abs(role_index) <= 0.3;

n_broadcasters = sum(is_broadcaster);
n_listeners = sum(is_listener);
n_balanced = sum(is_balanced);

fprintf('  Broadcasters (>0.3): %d (%.1f%%)\n', n_broadcasters, 100*n_broadcasters/n_clusters);
fprintf('  Listeners (<-0.3): %d (%.1f%%)\n', n_listeners, 100*n_listeners/n_clusters);
fprintf('  Balanced (+/-0.3): %d (%.1f%%)\n', n_balanced, 100*n_balanced/n_clusters);
fprintf('\n');

%% Phase 4: Dominant signaling vocabulary of top broadcasters

% For top broadcasters, what are their main signals?
[~, sort_idx] = sort(role_index, 'descend');
top_n_broadcasters = min(10, n_broadcasters);

if top_n_broadcasters > 0
    top_broadcaster_ids = all_cluster_ids(sort_idx(1:top_n_broadcasters));
else
    top_broadcaster_ids = [];
end

broadcaster_signals = table('Size', [0, 5], ...
    'VariableTypes', {'double', 'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'ClusterID', 'TopSignal', 'SignalCount', 'RoleIndex', 'TotalOut'});

if ~isempty(top_broadcaster_ids)
    ligs = gbm_network.ligands;
    recs = gbm_network.receptors;
    
    for i = 1:length(top_broadcaster_ids)
        cluster_id = top_broadcaster_ids(i);
        
        % Find all signals FROM this cluster
        from_this = src_clusters == cluster_id;
        
        if sum(from_this) == 0
            continue;
        end
        
        % Get L-R pairs
        from_indices = find(from_this);
        pairs_from_this = cell(length(from_indices), 1);
        for j = 1:length(from_indices)
            idx = from_indices(j);
            pairs_from_this{j} = [ligs{idx}, '->', recs{idx}];
        end
        
        % Find most common signal
        [unique_pairs, ~, ic] = unique(pairs_from_this);
        counts = accumarray(ic, 1);
        [max_count, max_idx] = max(counts);
        top_signal = unique_pairs{max_idx};
        
        % Get role index for this cluster
        cluster_idx = find(all_cluster_ids == cluster_id, 1);
        r_idx = role_index(cluster_idx);
        total_out = outgoing_counts(cluster_idx);
        
        broadcaster_signals = [broadcaster_signals; ...
            {cluster_id, top_signal, max_count, r_idx, total_out}];
    end
end

fprintf('  Top broadcaster signals identified: %d\n\n', height(broadcaster_signals));

%% Phase 5: Cross-reference broadcaster signals with condition-specific interactions

if istable(gbm_specific) && height(gbm_specific) > 0
    cols = gbm_specific.Properties.VariableNames;
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
    
    if ~isempty(lig_col) && ~isempty(rec_col)
        is_tumor_specific = false(height(broadcaster_signals), 1);
        
        for row = 1:height(broadcaster_signals)
            signal = broadcaster_signals.TopSignal{row};
            parts = strsplit(signal, '->');
            
            if length(parts) == 2
                for spec_row = 1:height(gbm_specific)
                    if strcmp(gbm_specific.(lig_col){spec_row}, parts{1}) && ...
                       strcmp(gbm_specific.(rec_col){spec_row}, parts{2})
                        is_tumor_specific(row) = true;
                        break;
                    end
                end
            end
        end
        
        broadcaster_signals.TumorSpecific = is_tumor_specific;
        n_tumor_broadcasters = sum(is_tumor_specific);
        
        fprintf('  Tumor-specific broadcaster signals: %d\n', n_tumor_broadcasters);
    else
        broadcaster_signals.TumorSpecific = false(height(broadcaster_signals), 1);
        fprintf('  Could not determine tumor specificity\n');
    end
else
    broadcaster_signals.TumorSpecific = false(height(broadcaster_signals), 1);
    fprintf('  No tumor-specific data available\n');
end

fprintf('\n');

%% Phase 6: Source diversity received by top listeners

[~, sort_idx] = sort(role_index, 'ascend');
top_n_listeners = min(10, n_listeners);

if top_n_listeners > 0
    top_listener_ids = all_cluster_ids(sort_idx(1:top_n_listeners));
else
    top_listener_ids = [];
end

listener_integration = table('Size', [0, 4], ...
    'VariableTypes', {'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'ClusterID', 'NumSources', 'TotalSignals', 'RoleIndex'});

if ~isempty(top_listener_ids)
    for i = 1:length(top_listener_ids)
        cluster_id = top_listener_ids(i);
        
        % Find all signals TO this cluster
        to_this = tgt_clusters == cluster_id;
        
        if sum(to_this) == 0
            continue;
        end
        
        % How many different source clusters?
        unique_sources = length(unique(src_clusters(to_this)));
        total_sigs = sum(to_this);
        
        cluster_idx = find(all_cluster_ids == cluster_id, 1);
        r_idx = role_index(cluster_idx);
        
        listener_integration = [listener_integration; ...
            {cluster_id, unique_sources, total_sigs, r_idx}];
    end
end

fprintf('  Top listeners identified: %d\n', height(listener_integration));
if height(listener_integration) > 0
    fprintf('  Mean sources per listener: %.1f\n', mean(listener_integration.NumSources));
end
fprintf('\n');

%% Phase 7: Export

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

% All vectors share the same length (n_clusters) for table construction.
role_table = table(all_cluster_ids, outgoing_counts, incoming_counts, ...
    total_signals, role_index, is_broadcaster, is_listener, is_balanced, ...
    'VariableNames', {'ClusterID', 'Outgoing', 'Incoming', 'Total', ...
                      'RoleIndex', 'IsBroadcaster', 'IsListener', 'IsBalanced'});

role_table = sortrows(role_table, 'RoleIndex', 'descend');

writetable(role_table, 'results/novel_discoveries/communication_roles.csv');
fprintf('  communication_roles.csv\n');

if height(broadcaster_signals) > 0
    writetable(broadcaster_signals, 'results/novel_discoveries/broadcaster_signals.csv');
    fprintf('  broadcaster_signals.csv\n');
end

if height(listener_integration) > 0
    writetable(listener_integration, 'results/novel_discoveries/listener_integration.csv');
    fprintf('  listener_integration.csv\n');
end

fprintf('\n');

%% Package
role_results = struct();
role_results.role_table = role_table;
role_results.broadcaster_signals = broadcaster_signals;
role_results.listener_integration = listener_integration;
role_results.n_broadcasters = n_broadcasters;
role_results.n_listeners = n_listeners;
role_results.n_balanced = n_balanced;

if height(broadcaster_signals) > 0 && ismember('TumorSpecific', broadcaster_signals.Properties.VariableNames)
    role_results.n_tumor_broadcasters = sum(broadcaster_signals.TumorSpecific);
else
    role_results.n_tumor_broadcasters = 0;
end

%% Summary
fprintf('Communication role asymmetry analysis complete\n\n');

fprintf('RESULTS:\n');
fprintf('  Total clusters analyzed: %d\n', n_clusters);
fprintf('  Broadcaster clusters: %d (%.1f%%)\n', n_broadcasters, 100*n_broadcasters/n_clusters);
fprintf('  Listener clusters: %d (%.1f%%)\n', n_listeners, 100*n_listeners/n_clusters);
fprintf('  Balanced clusters: %d (%.1f%%)\n', n_balanced, 100*n_balanced/n_clusters);
fprintf('  Tumor-specific broadcaster signals: %d\n', role_results.n_tumor_broadcasters);

if height(broadcaster_signals) > 0
    fprintf('\nTOP 5 BROADCASTER CLUSTERS:\n');
    for i = 1:min(5, height(broadcaster_signals))
        fprintf('  %d. Cluster %d: %s (%d signals', ...
            i, broadcaster_signals.ClusterID(i), ...
            broadcaster_signals.TopSignal{i}, ...
            broadcaster_signals.SignalCount(i));
        fprintf(', Role: %.3f', broadcaster_signals.RoleIndex(i));
        if broadcaster_signals.TumorSpecific(i)
            fprintf(', TUMOR-SPECIFIC');
        end
        fprintf(')\n');
    end
end

if height(listener_integration) > 0
    fprintf('\nTOP 5 LISTENER CLUSTERS:\n');
    for i = 1:min(5, height(listener_integration))
        fprintf('  %d. Cluster %d: Integrates %d sources (%d total signals, Role: %.3f)\n', ...
            i, listener_integration.ClusterID(i), ...
            listener_integration.NumSources(i), ...
            listener_integration.TotalSignals(i), ...
            listener_integration.RoleIndex(i));
    end
end

fprintf('\nAnalysis complete\n\n');

end