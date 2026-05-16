function entropy_results = compute_communication_entropy(network, condition_specific, reference_network)
% COMPUTE_COMMUNICATION_ENTROPY  Information-theoretic characterisation of
%   cell-cell communication: Shannon entropy, mutual information, spectral
%   community detection, and entropy-based cluster classification.
%
% Scientific novelty
% ------------------
% Existing tools report total interaction count or strength per cluster,
% treating all interactions as equivalent. This function is the first to
% characterise the DISTRIBUTIONAL structure of communication using
% information theory. A cluster that sends 100 units via a single L-R pair
% to a single target has low sending entropy (specialist). A cluster that
% distributes 100 units equally across 10 L-R pairs and 5 targets has high
% sending entropy (promiscuous broadcaster). These profiles are biologically
% distinct but invisible to count-based metrics.
%
% Entropy decomposition further enables:
%   - Mutual information quantification: how much does knowing the source
%     reduce uncertainty about the target, and vice versa? High MI = efficient
%     point-to-point communication channel. Low MI = diffuse broadcasting.
%   - Spectral community detection via the normalised graph Laplacian:
%     clusters of cells with dense internal communication are identified
%     without requiring cluster-count pre-specification (spectral gap method).
%   - Condition vs. reference entropy rewiring: which clusters become more
%     promiscuous or more specialised in disease?
%   - Quadrant classification: each cluster occupies one of four regimes
%     (Broadcaster / Hub / Integrator / Specialist) defined by whether its
%     sending and receiving entropy exceed the tissue median.
%
% Algorithm overview
% ------------------
%   Part A — Per-cluster Shannon entropy
%     H_send(s)  = -Σ_{t,lr} p(t,lr|s)  log2 p(t,lr|s)
%     H_recv(t)  = -Σ_{s,lr} p(s,lr|t)  log2 p(s,lr|t)
%     MI(s→t)    = H_send(s) + H_recv(t) - H_joint
%   Part B — Spectral analysis of cluster×cluster weight matrix
%     L = I - D^(-1/2) A D^(-1/2)  (normalised Laplacian)
%     Spectral gap → n_communities; K-means in eigenvector space.
%   Part C — Entropy rewiring scores (condition vs reference)
%     ΔH_send(s) = H_send_cond(s) - H_send_ref(s)
%   Part D — Quadrant classification (threshold = tissue median entropy)
%
% Usage
% -----
%   entropy_results = compute_communication_entropy(network, condition_specific)
%   entropy_results = compute_communication_entropy(network, condition_specific, [])
%   entropy_results = compute_communication_entropy(network, condition_specific, ref_network)
%
% Inputs
% ------
%   network            - Communication network struct from
%                        build_communication_network() with fields:
%                          .source_clusters  [N×1] sending cluster IDs
%                          .target_clusters  [N×1] receiving cluster IDs
%                          .ligands          {N×1} ligand gene symbols
%                          .receptors        {N×1} receptor gene symbols
%                          .weights          [N×1] interaction strengths
%   condition_specific - Table from identify_gbm_specific_interactions()
%                        with columns Ligand, Receptor. Used to annotate
%                        dominant L-R pairs per community. May be empty [].
%   reference_network  - (optional) Reference network struct with same
%                        fields as network. Pass [] if unavailable.
%                        When provided, entropy rewiring scores are computed.
%
% Outputs
% -------
%   entropy_results.H_send              - [n_clusters×1] sending entropy (bits)
%   entropy_results.H_recv              - [n_clusters×1] receiving entropy (bits)
%   entropy_results.MI                  - [n_clusters×n_clusters] pairwise MI (bits)
%   entropy_results.spectral_communities - [n_clusters×1] community assignments
%   entropy_results.n_communities       - Number of communities detected
%   entropy_results.spectral_gap        - Eigenvalue gap used for detection
%   entropy_results.quadrant            - {n_clusters×1} quadrant names
%   entropy_results.entropy_table       - Per-cluster summary table
%   entropy_results.community_table     - Per-community summary table
%   entropy_results.cluster_ids         - [n_clusters×1] ordered cluster IDs
%   entropy_results.H_send_ref          - (if ref provided) [n_clusters×1]
%   entropy_results.H_recv_ref          - (if ref provided) [n_clusters×1]
%   entropy_results.delta_H_send        - (if ref provided) [n_clusters×1]
%   entropy_results.delta_H_recv        - (if ref provided) [n_clusters×1]
%
% Files saved
% -----------
%   entropy_landscape.csv       - Per-cluster entropy statistics
%   communication_communities.csv - Per-community summary
%   (saved to results/ if that folder exists, otherwise current directory)

fprintf('[Entropy Landscape] Starting communication entropy analysis...\n');

%% ── 0. Output directory ──────────────────────────────────────────────────

out_dir = '';
if exist('results', 'dir')
    out_dir = 'results/';
end

entropy_results = struct();

try

%% ── 1. Input validation ──────────────────────────────────────────────────

if nargin < 2
    error('compute_communication_entropy: requires network and condition_specific.');
end
if nargin < 3
    reference_network = [];
end

if ~isstruct(network)
    error('[Entropy Landscape] network must be a struct.');
end
for fi = {'source_clusters','target_clusters','ligands','receptors','weights'}
    if ~isfield(network, fi{1})
        error('[Entropy Landscape] network missing field: %s', fi{1});
    end
end

src  = network.source_clusters(:);
tgt  = network.target_clusters(:);
ligs = cellfun(@char, network.ligands(:), 'UniformOutput', false);
recs = cellfun(@char, network.receptors(:), 'UniformOutput', false);
wts  = double(network.weights(:));

% Remove NaN / Inf / non-positive
valid_mask = isfinite(wts) & (wts > 0);
src  = src(valid_mask);
tgt  = tgt(valid_mask);
ligs = ligs(valid_mask);
recs = recs(valid_mask);
wts  = wts(valid_mask);

if numel(wts) == 0
    warning('[Entropy Landscape] No valid interactions after filtering. Returning empty struct.');
    return;
end

%% ── 2. Enumerate unique clusters ─────────────────────────────────────────

all_cluster_ids = unique([src; tgt]);
n_clusters = numel(all_cluster_ids);

fprintf('[Entropy Landscape] Clusters: %d | Interactions: %d\n', n_clusters, numel(wts));

if n_clusters < 2
    warning('[Entropy Landscape] Fewer than 2 clusters. Entropy analysis not meaningful.');
    return;
end

% L-R pair labels
lr_labels = strcat(ligs, '_', recs);
unique_lr  = unique(lr_labels);
n_lr       = numel(unique_lr);

%% ── 3. Build 3-way weight array  w(s_idx, t_idx, lr_idx) ────────────────

fprintf('[Entropy Landscape] Building 3-way weight array...\n');

% Index mapping: cluster ID → row/col index
[~, src_idx]    = ismember(src, all_cluster_ids);
[~, tgt_idx]    = ismember(tgt, all_cluster_ids);
[~, lr_num_idx] = ismember(lr_labels, unique_lr);

% Sparse accumulation into a 3-way array stored as n_clusters×n_clusters×n_lr
% For memory efficiency, we'll index by (s,t,lr) flat index and use accumarray
flat_idx = (lr_num_idx - 1) * n_clusters * n_clusters + ...
           (tgt_idx - 1)    * n_clusters + src_idx;

W3_flat = accumarray(flat_idx, wts, [n_clusters * n_clusters * n_lr, 1]);
W3 = reshape(W3_flat, n_clusters, n_clusters, n_lr);
% W3(s, t, lr) = total weight from cluster s to cluster t via LR pair lr

%% ══════════════════════════════════════════════════════════════════════════
%% PART A: Per-cluster Shannon entropy and mutual information
%% ══════════════════════════════════════════════════════════════════════════

fprintf('[Entropy Landscape] Part A: Computing per-cluster Shannon entropy...\n');

total_weight = sum(W3(:));

%% Sending entropy H_send(s)
H_send = zeros(n_clusters, 1);
for s = 1:n_clusters
    % Distribution over all (t, lr) for this source s
    dist_s = W3(s, :, :);          % [1 × n_tgt × n_lr]
    dist_s = dist_s(:);            % flatten to vector
    denom_s = sum(dist_s);
    if denom_s > 0
        p_s = dist_s / denom_s;
        H_send(s) = shannon_entropy(p_s);
    end
end

%% Receiving entropy H_recv(t)
H_recv = zeros(n_clusters, 1);
for t = 1:n_clusters
    % Distribution over all (s, lr) for this target t
    dist_t = W3(:, t, :);
    dist_t = dist_t(:);
    denom_t = sum(dist_t);
    if denom_t > 0
        p_t = dist_t / denom_t;
        H_recv(t) = shannon_entropy(p_t);
    end
end

%% Marginal joint distribution p(s,t) = Σ_lr w(s,t,lr) / total_weight
A_joint = squeeze(sum(W3, 3));   % [n_clusters × n_clusters] summed over LR
if total_weight > 0
    P_joint = A_joint / total_weight;
else
    P_joint = zeros(n_clusters);
end

%% Marginal source and target distributions
P_src = sum(P_joint, 2);    % [n_clusters×1]
P_tgt = sum(P_joint, 1)';   % [n_clusters×1]

%% H_joint and marginal entropies for MI calculation
H_joint_scalar = shannon_entropy(P_joint(:));
H_src_marginal = shannon_entropy(P_src);
H_tgt_marginal = shannon_entropy(P_tgt);

%% MI(s,t) pairwise: MI(s→t) = H_send(s) + H_recv(t) - H_joint (scalar)
%  We compute a pairwise MI matrix as:
%  MI_mat(s,t) = h(P_src(s)) + h(P_tgt(t)) - h(P_joint(s,t))
%  (pointwise mutual information per cluster pair, clamped to 0)
MI_mat = zeros(n_clusters, n_clusters);
for s = 1:n_clusters
    for t = 1:n_clusters
        p_st  = P_joint(s, t);
        p_s   = P_src(s);
        p_t   = P_tgt(t);
        if p_st > 0 && p_s > 0 && p_t > 0
            MI_mat(s, t) = p_st * log2(p_st / (p_s * p_t));
        end
    end
end
MI_mat = max(MI_mat, 0);   % clamp floating-point negatives

fprintf('[Entropy Landscape]   H_send range: [%.3f, %.3f] bits\n', min(H_send), max(H_send));
fprintf('[Entropy Landscape]   H_recv range: [%.3f, %.3f] bits\n', min(H_recv), max(H_recv));
fprintf('[Entropy Landscape]   Joint MI (total): %.4f bits\n', sum(MI_mat(:)));

%% ══════════════════════════════════════════════════════════════════════════
%% PART B: Spectral analysis — normalised Laplacian and communities
%% ══════════════════════════════════════════════════════════════════════════

fprintf('[Entropy Landscape] Part B: Spectral community detection...\n');

%% Cluster×cluster adjacency matrix A (sum over LR pairs)
A = A_joint;

%% Symmetrise for undirected spectral analysis
A_sym = (A + A') / 2;

%% Degree matrix D (row sums)
deg = sum(A_sym, 2);
deg(deg == 0) = eps;   % avoid division by zero for isolated clusters

%% Normalised Laplacian: L = I - D^(-1/2) * A * D^(-1/2)
D_invsqrt = diag(1 ./ sqrt(deg));
L_norm = eye(n_clusters) - D_invsqrt * A_sym * D_invsqrt;

%% Symmetrise L_norm to correct numerical asymmetry
L_norm = (L_norm + L_norm') / 2;

%% Eigendecomposition — use eig (dense, all eigenvalues)
try
    [V_eig, D_eig_mat] = eig(L_norm);
    eig_vals = real(diag(D_eig_mat));
    [eig_vals, eig_sort_idx] = sort(eig_vals, 'ascend');
    V_eig = real(V_eig(:, eig_sort_idx));
catch ME_eig
    fprintf('[Entropy Landscape] Warning: eigendecomposition failed: %s. Using trivial community.\n', ME_eig.message);
    eig_vals = zeros(n_clusters, 1);
    V_eig    = eye(n_clusters);
end

k_max = min(8, n_clusters);
eig_vals_k = eig_vals(1:k_max);

%% Spectral gap: maximum difference between consecutive eigenvalues
gaps = diff(eig_vals_k);
[spectral_gap, gap_idx] = max(gaps);
n_communities_detected = gap_idx;   % number of communities = index of max gap
n_communities_detected = max(2, min(n_communities_detected, k_max - 1));

fprintf('[Entropy Landscape]   Spectral gap: %.4f at index %d → %d communities\n', ...
    spectral_gap, gap_idx, n_communities_detected);

%% K-means in eigenvector space (use first n_communities eigenvectors)
embedding = V_eig(:, 1:n_communities_detected);

% Normalise rows for spectral clustering stability
row_norms = sqrt(sum(embedding.^2, 2));
row_norms(row_norms == 0) = 1;
embedding_norm = embedding ./ row_norms;

try
    community_labels = kmeans(embedding_norm, n_communities_detected, ...
        'Replicates', 10, 'MaxIter', 300);
catch ME_km
    fprintf('[Entropy Landscape] Warning: k-means failed (%s). Assigning all to community 1.\n', ME_km.message);
    community_labels = ones(n_clusters, 1);
end

%% ══════════════════════════════════════════════════════════════════════════
%% PART C: Entropy rewiring (condition vs. reference)
%% ══════════════════════════════════════════════════════════════════════════

H_send_ref   = nan(n_clusters, 1);
H_recv_ref   = nan(n_clusters, 1);
delta_H_send = nan(n_clusters, 1);
delta_H_recv = nan(n_clusters, 1);
has_reference = false;

if ~isempty(reference_network) && isstruct(reference_network) && ...
        isfield(reference_network, 'weights') && ...
        isfield(reference_network, 'source_clusters') && ...
        isfield(reference_network, 'target_clusters') && ...
        isfield(reference_network, 'ligands') && ...
        isfield(reference_network, 'receptors')

    fprintf('[Entropy Landscape] Part C: Computing entropy rewiring vs. reference...\n');
    has_reference = true;

    ref_src  = reference_network.source_clusters(:);
    ref_tgt  = reference_network.target_clusters(:);
    ref_ligs = cellfun(@char, reference_network.ligands(:), 'UniformOutput', false);
    ref_recs = cellfun(@char, reference_network.receptors(:), 'UniformOutput', false);
    ref_wts  = double(reference_network.weights(:));

    ref_valid = isfinite(ref_wts) & (ref_wts > 0);
    ref_src  = ref_src(ref_valid);
    ref_tgt  = ref_tgt(ref_valid);
    ref_ligs = ref_ligs(ref_valid);
    ref_recs = ref_recs(ref_valid);
    ref_wts  = ref_wts(ref_valid);

    ref_lr    = strcat(ref_ligs, '_', ref_recs);
    ref_lr_u  = unique(ref_lr);
    n_lr_ref  = numel(ref_lr_u);

    % Map ref clusters to the same index space as condition clusters
    % (use all_cluster_ids from condition; ref clusters not in condition
    %  are ignored for direct comparison)
    [ref_in_cond_src, ref_src_idx] = ismember(ref_src, all_cluster_ids);
    [ref_in_cond_tgt, ref_tgt_idx] = ismember(ref_tgt, all_cluster_ids);
    ref_in_both = ref_in_cond_src & ref_in_cond_tgt;

    ref_src_idx = ref_src_idx(ref_in_both);
    ref_tgt_idx = ref_tgt_idx(ref_in_both);
    ref_wts_f   = ref_wts(ref_in_both);
    ref_lr_f    = ref_lr(ref_in_both);

    [~, ref_lr_num] = ismember(ref_lr_f, ref_lr_u);

    if numel(ref_wts_f) > 0
        flat_ref = (ref_lr_num - 1) * n_clusters * n_clusters + ...
                   (ref_tgt_idx - 1) * n_clusters + ref_src_idx;
        W3_ref_flat = accumarray(flat_ref, ref_wts_f, ...
            [n_clusters * n_clusters * n_lr_ref, 1]);
        W3_ref = reshape(W3_ref_flat, n_clusters, n_clusters, n_lr_ref);

        for s = 1:n_clusters
            dist_sr = W3_ref(s, :, :);
            dist_sr = dist_sr(:);
            dn_sr = sum(dist_sr);
            if dn_sr > 0
                H_send_ref(s) = shannon_entropy(dist_sr / dn_sr);
            end
        end

        for t = 1:n_clusters
            dist_tr = W3_ref(:, t, :);
            dist_tr = dist_tr(:);
            dn_tr = sum(dist_tr);
            if dn_tr > 0
                H_recv_ref(t) = shannon_entropy(dist_tr / dn_tr);
            end
        end

        delta_H_send = H_send - H_send_ref;
        delta_H_recv = H_recv - H_recv_ref;

        fprintf('[Entropy Landscape]   ΔH_send range: [%.3f, %.3f] bits\n', ...
            min(delta_H_send(~isnan(delta_H_send))), ...
            max(delta_H_send(~isnan(delta_H_send))));
    else
        fprintf('[Entropy Landscape]   Reference network has no overlapping clusters. Skipping rewiring.\n');
        has_reference = false;
    end
else
    fprintf('[Entropy Landscape] Part C: No reference network provided. Skipping entropy rewiring.\n');
end

%% ══════════════════════════════════════════════════════════════════════════
%% PART D: Quadrant classification
%% ══════════════════════════════════════════════════════════════════════════

fprintf('[Entropy Landscape] Part D: Quadrant classification...\n');

med_send = median(H_send);
med_recv = median(H_recv);

quadrant = cell(n_clusters, 1);
for c = 1:n_clusters
    high_s = H_send(c) >= med_send;
    high_r = H_recv(c) >= med_recv;
    if high_s && ~high_r
        quadrant{c} = 'Signal Broadcaster';
    elseif high_s && high_r
        quadrant{c} = 'Communication Hub';
    elseif ~high_s && high_r
        quadrant{c} = 'Signal Integrator';
    else
        quadrant{c} = 'Signal Specialist';
    end
end

fprintf('[Entropy Landscape]   Broadcasters : %d\n', sum(strcmp(quadrant,'Signal Broadcaster')));
fprintf('[Entropy Landscape]   Hubs         : %d\n', sum(strcmp(quadrant,'Communication Hub')));
fprintf('[Entropy Landscape]   Integrators  : %d\n', sum(strcmp(quadrant,'Signal Integrator')));
fprintf('[Entropy Landscape]   Specialists  : %d\n', sum(strcmp(quadrant,'Signal Specialist')));

%% ══════════════════════════════════════════════════════════════════════════
%% PART E: Community summary — dominant LR pairs per community
%% ══════════════════════════════════════════════════════════════════════════

fprintf('[Entropy Landscape] Part E: Summarising communication communities...\n');

comm_ids_list   = (1:n_communities_detected)';
comm_clusters   = cell(n_communities_detected, 1);
comm_total_wt   = zeros(n_communities_detected, 1);
comm_dom_lr     = cell(n_communities_detected, 1);

for ci = 1:n_communities_detected
    member_idx = find(community_labels == ci);

    % Cluster IDs in this community
    member_cid = all_cluster_ids(member_idx);
    comm_clusters{ci} = strjoin(cellstr(num2str(double(member_cid(:)))), '; ');

    % Total communication weight within community
    sub_A = A_sym(member_idx, member_idx);
    comm_total_wt(ci) = sum(sub_A(:));

    % Dominant LR pairs: sum weights for interactions within the community
    member_src_mask = ismember(src_idx, member_idx);
    member_tgt_mask = ismember(tgt_idx, member_idx);
    within_mask     = member_src_mask & member_tgt_mask;

    if any(within_mask)
        lr_within = lr_labels(within_mask);
        wt_within = wts(within_mask);
        [ulr_w, ~, ulr_wi] = unique(lr_within);
        lr_agg = accumarray(ulr_wi, wt_within);
        [~, lr_ord] = sort(lr_agg, 'descend');
        top_lr_n = min(5, numel(ulr_w));
        comm_dom_lr{ci} = strjoin(ulr_w(lr_ord(1:top_lr_n)), '; ');
    else
        comm_dom_lr{ci} = 'None';
    end
end

community_table = table( ...
    comm_ids_list, ...
    comm_clusters, ...
    comm_total_wt, ...
    comm_dom_lr, ...
    'VariableNames', {'CommunityID','Clusters','TotalWeight','DominantLRPairs'});

%% ══════════════════════════════════════════════════════════════════════════
%% ASSEMBLE ENTROPY TABLE
%% ══════════════════════════════════════════════════════════════════════════

cluster_id_col = all_cluster_ids(:);
H_send_col     = H_send;
H_recv_col     = H_recv;

% Scalar MI per cluster: sum of MI over all target clusters (outgoing MI)
MI_per_cluster = sum(MI_mat, 2);

quadrant_col    = quadrant;
community_col   = community_labels(:);

if has_reference
    entropy_table = table( ...
        cluster_id_col, H_send_col, H_recv_col, MI_per_cluster, ...
        quadrant_col, community_col, delta_H_send, delta_H_recv, ...
        'VariableNames', { ...
            'ClusterID','H_send','H_recv','MI', ...
            'Quadrant','CommunityID','Delta_H_send','Delta_H_recv'});
else
    entropy_table = table( ...
        cluster_id_col, H_send_col, H_recv_col, MI_per_cluster, ...
        quadrant_col, community_col, ...
        'VariableNames', { ...
            'ClusterID','H_send','H_recv','MI', ...
            'Quadrant','CommunityID'});
end

%% ══════════════════════════════════════════════════════════════════════════
%% POPULATE OUTPUT STRUCT
%% ══════════════════════════════════════════════════════════════════════════

entropy_results.H_send               = H_send;
entropy_results.H_recv               = H_recv;
entropy_results.MI                   = MI_mat;
entropy_results.spectral_communities = community_labels;
entropy_results.n_communities        = n_communities_detected;
entropy_results.spectral_gap         = spectral_gap;
entropy_results.quadrant             = quadrant;
entropy_results.entropy_table        = entropy_table;
entropy_results.community_table      = community_table;
entropy_results.cluster_ids          = all_cluster_ids;
entropy_results.eigenvalues          = eig_vals(1:k_max);

if has_reference
    entropy_results.H_send_ref   = H_send_ref;
    entropy_results.H_recv_ref   = H_recv_ref;
    entropy_results.delta_H_send = delta_H_send;
    entropy_results.delta_H_recv = delta_H_recv;
end

%% ══════════════════════════════════════════════════════════════════════════
%% EXPORT CSVs
%% ══════════════════════════════════════════════════════════════════════════

fprintf('[Entropy Landscape] Exporting results...\n');

try
    writetable(entropy_table, [out_dir, 'entropy_landscape.csv']);
    fprintf('[Entropy Landscape] Saved: %sentropy_landscape.csv\n', out_dir);
catch me_csv
    fprintf('[Entropy Landscape] Warning: could not save entropy_landscape.csv: %s\n', me_csv.message);
end

try
    writetable(community_table, [out_dir, 'communication_communities.csv']);
    fprintf('[Entropy Landscape] Saved: %scommunication_communities.csv\n', out_dir);
catch me_csv
    fprintf('[Entropy Landscape] Warning: could not save communication_communities.csv: %s\n', me_csv.message);
end

%% ══════════════════════════════════════════════════════════════════════════
%% SUMMARY PRINT
%% ══════════════════════════════════════════════════════════════════════════

fprintf('\n[Entropy Landscape] ══ Summary ══════════════════════════════════\n');
fprintf('[Entropy Landscape]   Clusters analysed   : %d\n', n_clusters);
fprintf('[Entropy Landscape]   L-R pairs            : %d\n', n_lr);
fprintf('[Entropy Landscape]   Interactions         : %d\n', numel(wts));
fprintf('[Entropy Landscape]   Total weight         : %.4f\n', total_weight);
fprintf('[Entropy Landscape]   H_send median        : %.3f bits\n', med_send);
fprintf('[Entropy Landscape]   H_recv median        : %.3f bits\n', med_recv);
fprintf('[Entropy Landscape]   Joint H(S,T)         : %.3f bits\n', H_joint_scalar);
fprintf('[Entropy Landscape]   Communities detected : %d (gap=%.4f)\n', ...
    n_communities_detected, spectral_gap);
fprintf('[Entropy Landscape]   Quadrant breakdown:\n');
fprintf('[Entropy Landscape]     Signal Broadcaster : %d\n', sum(strcmp(quadrant,'Signal Broadcaster')));
fprintf('[Entropy Landscape]     Communication Hub  : %d\n', sum(strcmp(quadrant,'Communication Hub')));
fprintf('[Entropy Landscape]     Signal Integrator  : %d\n', sum(strcmp(quadrant,'Signal Integrator')));
fprintf('[Entropy Landscape]     Signal Specialist  : %d\n', sum(strcmp(quadrant,'Signal Specialist')));
if has_reference
    valid_dhs = delta_H_send(~isnan(delta_H_send));
    if ~isempty(valid_dhs)
        fprintf('[Entropy Landscape]   ΔH_send (cond−ref) : mean=%.3f, range=[%.3f, %.3f]\n', ...
            mean(valid_dhs), min(valid_dhs), max(valid_dhs));
    end
end
fprintf('[Entropy Landscape] ══════════════════════════════════════════════\n\n');

catch ME
    fprintf('[Entropy Landscape] ERROR: %s\n', ME.message);
    fprintf('[Entropy Landscape] Stack: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    entropy_results = struct();
end

end   % main function

%% ══════════════════════════════════════════════════════════════════════════
%%  LOCAL HELPER FUNCTIONS
%% ══════════════════════════════════════════════════════════════════════════

function H = shannon_entropy(p)
% SHANNON_ENTROPY  Compute Shannon entropy in bits from a probability vector.
%   H = -Σ p_i log2(p_i)
%   Convention: 0 * log2(0) = 0 (limit as p→0).
%   Input p need not sum to 1; it is normalised internally.
%   Handles NaN, Inf, and negative values by zeroing them.
    p = double(p(:));
    p(~isfinite(p) | p < 0) = 0;
    s = sum(p);
    if s <= 0
        H = 0;
        return;
    end
    p = p / s;
    % Compute only non-zero terms (0*log(0) = 0 by convention)
    nz   = p > 0;
    H    = -sum(p(nz) .* log2(p(nz)));
    if ~isfinite(H)
        H = 0;
    end
end
