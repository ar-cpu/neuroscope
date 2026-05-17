function entropy_results = compute_communication_entropy(network, condition_specific, reference_network)
% COMPUTE_COMMUNICATION_ENTROPY  Information-theoretic characterisation of
%   cell-cell communication: Shannon entropy, mutual information, spectral
%   community detection, and entropy-based cluster classification.
%
% Usage:
%   entropy_results = compute_communication_entropy(network, condition_specific)
%   entropy_results = compute_communication_entropy(network, condition_specific, ref_network)
%
% Inputs:
%   network            - struct from build_communication_network() with fields:
%                        .source_clusters, .target_clusters, .ligands,
%                        .receptors, .weights
%   condition_specific - table from identify_gbm_specific_interactions() with
%                        columns Ligand, Receptor; may be empty []
%   reference_network  - (optional) same-format struct; pass [] to skip rewiring
%
% Outputs:
%   .H_send / .H_recv           [n_clusters x 1] sending/receiving entropy (bits)
%   .MI                         [n_clusters x n_clusters] pairwise MI (bits)
%   .spectral_communities       [n_clusters x 1] community assignments
%   .n_communities / .spectral_gap
%   .quadrant                   {n_clusters x 1} Broadcaster/Hub/Integrator/Specialist
%   .entropy_table / .community_table
%   .H_send_ref / .H_recv_ref / .delta_H_send / .delta_H_recv  (if ref provided)
%
% Saves: entropy_landscape.csv, communication_communities.csv -> results/ or ./

fprintf('[Entropy Landscape] Starting communication entropy analysis...\n');

out_dir = '';
if exist('results', 'dir')
    out_dir = 'results/';
end

entropy_results = struct();

try

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
ligs = cellfun(@char, network.ligands(:),   'UniformOutput', false);
recs = cellfun(@char, network.receptors(:), 'UniformOutput', false);
wts  = double(network.weights(:));

valid_mask = isfinite(wts) & (wts > 0);
src  = src(valid_mask);  tgt  = tgt(valid_mask);
ligs = ligs(valid_mask); recs = recs(valid_mask);
wts  = wts(valid_mask);

if numel(wts) == 0
    warning('[Entropy Landscape] No valid interactions after filtering. Returning empty struct.');
    return;
end

all_cluster_ids = unique([src; tgt]);
n_clusters = numel(all_cluster_ids);

fprintf('[Entropy Landscape] Clusters: %d | Interactions: %d\n', n_clusters, numel(wts));

if n_clusters < 2
    warning('[Entropy Landscape] Fewer than 2 clusters. Entropy analysis not meaningful.');
    return;
end

lr_labels = strcat(ligs, '_', recs);
unique_lr  = unique(lr_labels);
n_lr       = numel(unique_lr);

% Build W3(s, t, lr) = total weight from cluster s to cluster t via LR pair lr
fprintf('[Entropy Landscape] Building 3-way weight array...\n');

[~, src_idx]    = ismember(src, all_cluster_ids);
[~, tgt_idx]    = ismember(tgt, all_cluster_ids);
[~, lr_num_idx] = ismember(lr_labels, unique_lr);

flat_idx = (lr_num_idx - 1) * n_clusters * n_clusters + ...
           (tgt_idx   - 1) * n_clusters + src_idx;
W3_flat = accumarray(flat_idx, wts, [n_clusters * n_clusters * n_lr, 1]);
W3 = reshape(W3_flat, n_clusters, n_clusters, n_lr);

% Part A: Per-cluster Shannon entropy and mutual information
fprintf('[Entropy Landscape] Part A: Shannon entropy...\n');

total_weight = sum(W3(:));

H_send = zeros(n_clusters, 1);
for s = 1:n_clusters
    dist_s = W3(s, :, :);
    dist_s = dist_s(:);
    denom_s = sum(dist_s);
    if denom_s > 0
        H_send(s) = shannon_entropy(dist_s / denom_s);
    end
end

H_recv = zeros(n_clusters, 1);
for t = 1:n_clusters
    dist_t = W3(:, t, :);
    dist_t = dist_t(:);
    denom_t = sum(dist_t);
    if denom_t > 0
        H_recv(t) = shannon_entropy(dist_t / denom_t);
    end
end

A_joint = squeeze(sum(W3, 3));
P_joint = A_joint / max(total_weight, eps);
P_src   = sum(P_joint, 2);
P_tgt   = sum(P_joint, 1)';

H_joint_scalar = shannon_entropy(P_joint(:));
H_src_marginal = shannon_entropy(P_src);  %#ok<NASGU>
H_tgt_marginal = shannon_entropy(P_tgt);  %#ok<NASGU>

% Pointwise MI per cluster pair, clamped to 0
MI_mat = zeros(n_clusters, n_clusters);
for s = 1:n_clusters
    for t = 1:n_clusters
        p_st = P_joint(s, t);
        p_s  = P_src(s);
        p_t  = P_tgt(t);
        if p_st > 0 && p_s > 0 && p_t > 0
            MI_mat(s, t) = p_st * log2(p_st / (p_s * p_t));
        end
    end
end
MI_mat = max(MI_mat, 0);

fprintf('[Entropy Landscape]   H_send range: [%.3f, %.3f] bits\n', min(H_send), max(H_send));
fprintf('[Entropy Landscape]   H_recv range: [%.3f, %.3f] bits\n', min(H_recv), max(H_recv));
fprintf('[Entropy Landscape]   Joint MI (total): %.4f bits\n', sum(MI_mat(:)));

% Part B: Spectral community detection via normalised Laplacian
fprintf('[Entropy Landscape] Part B: Spectral community detection...\n');

A_sym = (A_joint + A_joint') / 2;
deg   = sum(A_sym, 2);
deg(deg == 0) = eps;
D_invsqrt = diag(1 ./ sqrt(deg));
L_norm     = eye(n_clusters) - D_invsqrt * A_sym * D_invsqrt;
L_norm     = (L_norm + L_norm') / 2;

try
    [V_eig, D_eig_mat] = eig(L_norm);
    eig_vals = real(diag(D_eig_mat));
    [eig_vals, eig_sort_idx] = sort(eig_vals, 'ascend');
    V_eig = real(V_eig(:, eig_sort_idx));
catch ME_eig
    fprintf('[Entropy Landscape] Warning: eigendecomposition failed: %s\n', ME_eig.message);
    eig_vals = zeros(n_clusters, 1);
    V_eig    = eye(n_clusters);
end

k_max        = min(8, n_clusters);
eig_vals_k   = eig_vals(1:k_max);
gaps         = diff(eig_vals_k);
[spectral_gap, gap_idx] = max(gaps);
n_communities_detected  = max(2, min(gap_idx, k_max - 1));

fprintf('[Entropy Landscape]   Spectral gap: %.4f at index %d -> %d communities\n', ...
    spectral_gap, gap_idx, n_communities_detected);

embedding  = V_eig(:, 1:n_communities_detected);
row_norms  = sqrt(sum(embedding.^2, 2));
row_norms(row_norms == 0) = 1;
embedding_norm = embedding ./ row_norms;

try
    community_labels = kmeans(embedding_norm, n_communities_detected, ...
        'Replicates', 10, 'MaxIter', 300);
catch ME_km
    fprintf('[Entropy Landscape] Warning: k-means failed (%s). All -> community 1.\n', ME_km.message);
    community_labels = ones(n_clusters, 1);
end

% Part C: Entropy rewiring (condition vs. reference)
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

    fprintf('[Entropy Landscape] Part C: Entropy rewiring vs. reference...\n');
    has_reference = true;

    ref_src  = reference_network.source_clusters(:);
    ref_tgt  = reference_network.target_clusters(:);
    ref_ligs = cellfun(@char, reference_network.ligands(:),   'UniformOutput', false);
    ref_recs = cellfun(@char, reference_network.receptors(:), 'UniformOutput', false);
    ref_wts  = double(reference_network.weights(:));

    ref_valid = isfinite(ref_wts) & (ref_wts > 0);
    ref_src  = ref_src(ref_valid);  ref_tgt  = ref_tgt(ref_valid);
    ref_ligs = ref_ligs(ref_valid); ref_recs = ref_recs(ref_valid);
    ref_wts  = ref_wts(ref_valid);

    ref_lr   = strcat(ref_ligs, '_', ref_recs);
    ref_lr_u = unique(ref_lr);
    n_lr_ref = numel(ref_lr_u);

    % Only compare clusters present in both condition and reference
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
            dist_sr = W3_ref(s, :, :); dist_sr = dist_sr(:);
            dn_sr = sum(dist_sr);
            if dn_sr > 0, H_send_ref(s) = shannon_entropy(dist_sr / dn_sr); end
        end
        for t = 1:n_clusters
            dist_tr = W3_ref(:, t, :); dist_tr = dist_tr(:);
            dn_tr = sum(dist_tr);
            if dn_tr > 0, H_recv_ref(t) = shannon_entropy(dist_tr / dn_tr); end
        end

        delta_H_send = H_send - H_send_ref;
        delta_H_recv = H_recv - H_recv_ref;

        fprintf('[Entropy Landscape]   dH_send range: [%.3f, %.3f] bits\n', ...
            min(delta_H_send(~isnan(delta_H_send))), max(delta_H_send(~isnan(delta_H_send))));
    else
        fprintf('[Entropy Landscape]   No overlapping clusters. Skipping rewiring.\n');
        has_reference = false;
    end
else
    fprintf('[Entropy Landscape] Part C: No reference network. Skipping.\n');
end

% Part D: Quadrant classification (threshold = tissue median)
fprintf('[Entropy Landscape] Part D: Quadrant classification...\n');

med_send = median(H_send);
med_recv = median(H_recv);

quadrant = cell(n_clusters, 1);
for c = 1:n_clusters
    high_s = H_send(c) >= med_send;
    high_r = H_recv(c) >= med_recv;
    if     high_s && ~high_r,  quadrant{c} = 'Signal Broadcaster';
    elseif high_s &&  high_r,  quadrant{c} = 'Communication Hub';
    elseif ~high_s && high_r,  quadrant{c} = 'Signal Integrator';
    else,                       quadrant{c} = 'Signal Specialist';
    end
end

fprintf('[Entropy Landscape]   Broadcasters: %d | Hubs: %d | Integrators: %d | Specialists: %d\n', ...
    sum(strcmp(quadrant,'Signal Broadcaster')), sum(strcmp(quadrant,'Communication Hub')), ...
    sum(strcmp(quadrant,'Signal Integrator')),  sum(strcmp(quadrant,'Signal Specialist')));

% Part E: Community summary — dominant LR pairs per community
fprintf('[Entropy Landscape] Part E: Community summaries...\n');

comm_ids_list = (1:n_communities_detected)';
comm_clusters_str = cell(n_communities_detected, 1);
comm_total_wt     = zeros(n_communities_detected, 1);
comm_dom_lr       = cell(n_communities_detected, 1);

for ci = 1:n_communities_detected
    member_idx = find(community_labels == ci);
    member_cid = all_cluster_ids(member_idx);
    comm_clusters_str{ci} = strjoin(cellstr(num2str(double(member_cid(:)))), '; ');

    sub_A = A_sym(member_idx, member_idx);
    comm_total_wt(ci) = sum(sub_A(:));

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

community_table = table(comm_ids_list, comm_clusters_str, comm_total_wt, comm_dom_lr, ...
    'VariableNames', {'CommunityID','Clusters','TotalWeight','DominantLRPairs'});

% Assemble entropy table
MI_per_cluster = sum(MI_mat, 2);
cluster_id_col = all_cluster_ids(:);
community_col  = community_labels(:);

if has_reference
    entropy_table = table(cluster_id_col, H_send, H_recv, MI_per_cluster, ...
        quadrant, community_col, delta_H_send, delta_H_recv, ...
        'VariableNames', {'ClusterID','H_send','H_recv','MI', ...
                          'Quadrant','CommunityID','Delta_H_send','Delta_H_recv'});
else
    entropy_table = table(cluster_id_col, H_send, H_recv, MI_per_cluster, ...
        quadrant, community_col, ...
        'VariableNames', {'ClusterID','H_send','H_recv','MI','Quadrant','CommunityID'});
end

% Populate output struct
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

% Export CSVs
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

fprintf('\n[Entropy Landscape] Summary:\n');
fprintf('  Clusters: %d | L-R pairs: %d | Interactions: %d\n', n_clusters, n_lr, numel(wts));
fprintf('  H_send median: %.3f bits | H_recv median: %.3f bits\n', med_send, med_recv);
fprintf('  Joint H(S,T): %.3f bits | Communities: %d (gap=%.4f)\n', ...
    H_joint_scalar, n_communities_detected, spectral_gap);
if has_reference
    valid_dhs = delta_H_send(~isnan(delta_H_send));
    if ~isempty(valid_dhs)
        fprintf('  dH_send (cond-ref): mean=%.3f, range=[%.3f, %.3f]\n', ...
            mean(valid_dhs), min(valid_dhs), max(valid_dhs));
    end
end
fprintf('\n');

catch ME
    fprintf('[Entropy Landscape] ERROR: %s\n', ME.message);
    fprintf('[Entropy Landscape] Stack: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    entropy_results = struct();
end

end   % main function


function H = shannon_entropy(p)
% Shannon entropy in bits; normalises p internally; 0*log(0) = 0.
    p = double(p(:));
    p(~isfinite(p) | p < 0) = 0;
    s = sum(p);
    if s <= 0, H = 0; return; end
    p  = p / s;
    nz = p > 0;
    H  = -sum(p(nz) .* log2(p(nz)));
    if ~isfinite(H), H = 0; end
end
