function identity = compute_communication_identity(cell_comm_profiles, cell_metadata)
% COMPUTE_COMMUNICATION_IDENTITY  ARI/NMI between communication states and
%   transcriptomic clusters, plus per-cell divergence score (CDS).
%   CDS_i = 1 - P(comm_state=s_i | trans_cluster=c_i); high-CDS cells
%   communicate differently from their transcriptomic neighbours.
%
% identity = compute_communication_identity(cell_comm_profiles, cell_metadata)
%
% Outputs (struct):
%   .ari / .nmi              - clustering agreement metrics
%   .divergence              - [n_cells x 1] per-cell CDS
%   .n_divergent             - cells with CDS > 0.5
%   .interpretation          - 'Decoupled' / 'Partial' / 'Coupled'
%   .divergent_cells_table   - table sorted by CDS

fprintf('Computing communication identity...\n');

identity = struct('ari', NaN, 'nmi', NaN, 'divergence', [], ...
                  'n_divergent', 0, 'divergent_fraction', 0, ...
                  'contingency', [], 'interpretation', '', ...
                  'divergent_cells_table', table());

% Validate inputs
if ~isstruct(cell_comm_profiles) || ...
   ~isfield(cell_comm_profiles, 'comm_cluster_assignments') || ...
   isempty(cell_comm_profiles.comm_cluster_assignments)
    fprintf('  Communication identity skipped: no comm_cluster_assignments in profiles.\n\n');
    return;
end

comm_states    = cell_comm_profiles.comm_cluster_assignments(:);  % [n_cells × 1]
trans_clusters = cell_metadata.Cluster(:);                        % [n_cells × 1]
n_cells        = length(comm_states);

if length(trans_clusters) ~= n_cells
    fprintf('  Communication identity skipped: cell count mismatch (%d vs %d).\n\n', ...
        n_cells, length(trans_clusters));
    return;
end

fprintf('  %d cells, %d communication states, %d transcriptomic clusters\n', ...
    n_cells, length(unique(comm_states)), length(unique(trans_clusters)));

% ARI + NMI
identity.ari = local_ari(comm_states, trans_clusters);
identity.nmi = local_nmi(comm_states, trans_clusters);

if identity.ari < 0.10
    identity.interpretation = 'Decoupled';
    interp_note = 'Communication identity is ORTHOGONAL to transcriptomic clusters.';
elseif identity.ari < 0.40
    identity.interpretation = 'Partial';
    interp_note = 'Communication identity PARTIALLY overlaps transcriptomic clusters.';
else
    identity.interpretation = 'Coupled';
    interp_note = 'Communication identity largely follows transcriptomic clusters.';
end

fprintf('  ARI = %.4f | NMI = %.4f | %s\n', identity.ari, identity.nmi, identity.interpretation);
fprintf('  %s\n', interp_note);

% contingency table
unique_states    = unique(comm_states);
unique_trans     = unique(trans_clusters);
n_states         = length(unique_states);
n_trans          = length(unique_trans);

% Map labels to 1..n indices for fast array indexing
state_map  = containers.Map(num2cell(double(unique_states)),  num2cell(1:n_states));
trans_map  = containers.Map(num2cell(double(unique_trans)),   num2cell(1:n_trans));

si = cellfun(@(x) state_map(x), num2cell(double(comm_states)));
ci = cellfun(@(x) trans_map(x), num2cell(double(trans_clusters)));

contingency = accumarray([si, ci], 1, [n_states, n_trans]);

% per-cell divergence score: CDS_i = 1 - contingency[s_i,c_i]/cluster_size[c_i]
cluster_sizes = sum(contingency, 1);          % [1 × n_trans]
pair_counts   = contingency(sub2ind([n_states, n_trans], si, ci));  % count(s_i, c_i)
cluster_n_i   = cluster_sizes(ci);            % cluster size for each cell

divergence = 1 - pair_counts(:) ./ max(cluster_n_i(:), 1);
divergence = max(0, min(1, divergence));      % clamp to [0,1]

n_divergent          = sum(divergence > 0.5);
divergent_fraction   = n_divergent / n_cells;

fprintf('  Divergent cells (CDS > 0.5): %d / %d  (%.1f%%)\n', ...
    n_divergent, n_cells, 100 * divergent_fraction);

if divergent_fraction > 0.15
    fprintf('  High divergence: substantial population communicates differently\n');
    fprintf('  from their transcriptomic neighbours — strong novel finding.\n');
end

% per-cell output table
if ismember('CellID', cell_metadata.Properties.VariableNames)
    cell_ids = cell_metadata.CellID;
else
    cell_ids = (1:n_cells)';
end

divergent_cells_table = table(cell_ids, trans_clusters, comm_states, divergence, ...
    'VariableNames', {'CellID', 'TranscriptomicCluster', 'CommunicationState', 'DivergenceScore'});
divergent_cells_table = sortrows(divergent_cells_table, 'DivergenceScore', 'descend');

% Save
out_dir = 'results';
if ~exist(out_dir, 'dir'), out_dir = '.'; end
try
    writetable(divergent_cells_table, fullfile(out_dir, 'communication_identity.csv'));
    fprintf('  Saved: communication_identity.csv\n');
catch
end

fprintf('\n');

% assemble output
identity.divergence           = divergence;
identity.n_divergent          = n_divergent;
identity.divergent_fraction   = divergent_fraction;
identity.contingency          = contingency;
identity.comm_states          = unique_states;
identity.trans_clusters       = unique_trans;
identity.divergent_cells_table = divergent_cells_table;

end


function ari = local_ari(a, b)
% Adjusted Rand Index via contingency table.
n  = length(a);
ua = unique(a);  ub = unique(b);
na = length(ua); nb = length(ub);

% Build contingency
a_idx = arrayfun(@(x) find(ua == x, 1), a);
b_idx = arrayfun(@(x) find(ub == x, 1), b);
C     = accumarray([a_idx, b_idx], 1, [na, nb]);

sum_cij2 = sum(C(:) .* (C(:)-1)) / 2;        % sum C(i,j) choose 2
row_sums  = sum(C, 2);
col_sums  = sum(C, 1);
sum_ai2   = sum(row_sums .* (row_sums-1)) / 2;
sum_bj2   = sum(col_sums .* (col_sums-1)) / 2;
total_c2  = n*(n-1)/2;

expected  = sum_ai2 * sum_bj2 / max(total_c2, 1);
max_val   = (sum_ai2 + sum_bj2) / 2;

if max_val - expected < 1e-12
    ari = 1.0;
else
    ari = (sum_cij2 - expected) / (max_val - expected);
end
end


function nmi = local_nmi(a, b)
% Normalised Mutual Information.
n  = length(a);
ua = unique(a);  ub = unique(b);

H_a = local_entropy(a, ua, n);
H_b = local_entropy(b, ub, n);

H_ab = 0;
for i = 1:length(ua)
    for j = 1:length(ub)
        p = sum(a == ua(i) & b == ub(j)) / n;
        if p > 0, H_ab = H_ab - p * log2(p); end
    end
end

MI = H_a + H_b - H_ab;
denom = H_a + H_b;
nmi = 0;
if denom > 1e-12
    nmi = 2 * MI / denom;
end
end


function H = local_entropy(labels, unique_labels, n)
H = 0;
for i = 1:length(unique_labels)
    p = sum(labels == unique_labels(i)) / n;
    if p > 0, H = H - p * log2(p); end
end
end
