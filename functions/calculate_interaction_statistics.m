function stats = calculate_interaction_statistics(gbm_network, gbm_specific)
% CALCULATE_INTERACTION_STATISTICS  Compute a statistical profile for each
%   condition-specific L-R interaction: descriptive stats, bootstrap CIs,
%   effect sizes, spatial heterogeneity, betweenness centrality, Estrada-
%   Hatano subgraph centrality, and a weighted composite score.
%
% Usage:
%   stats = calculate_interaction_statistics(gbm_network, gbm_specific)
%
% Inputs:
%   gbm_network  - Network struct from build_communication_network() with
%                  fields: .ligands, .receptors, .weights.
%   gbm_specific - Table of condition-specific interactions from
%                  identify_gbm_specific_interactions().
%
% Outputs:
%   stats - Table with columns: Ligand, Receptor, Frequency, MeanStrength,
%     MedianStrength, MaxStrength, StdStrength, CV, MAD, CI_Lower, CI_Upper,
%     EffectSize_Cohen, EffectSize_Cliff, SpatialHeterogeneity, Betweenness,
%     Communicability, ClusterDiversity, CompositeScore.

fprintf('Interaction Statistics - Bootstrap: 100 iters\n\n');
fprintf('Computing advanced statistics for %d interactions...\n\n', height(gbm_specific));

n_interactions = height(gbm_specific);

stats = table();
stats.Ligand = gbm_specific.Ligand;
stats.Receptor = gbm_specific.Receptor;
stats.Frequency = zeros(n_interactions, 1);
stats.MeanStrength = zeros(n_interactions, 1);
stats.MedianStrength = zeros(n_interactions, 1);
stats.MaxStrength = zeros(n_interactions, 1);
stats.StdStrength = zeros(n_interactions, 1);
stats.CV = zeros(n_interactions, 1);
stats.MAD = zeros(n_interactions, 1);
stats.CI_Lower = zeros(n_interactions, 1);
stats.CI_Upper = zeros(n_interactions, 1);
stats.EffectSize_Cohen = zeros(n_interactions, 1);
stats.EffectSize_Cliff = zeros(n_interactions, 1);
stats.SpatialHeterogeneity = zeros(n_interactions, 1);
stats.Betweenness = zeros(n_interactions, 1);
stats.Communicability = zeros(n_interactions, 1);
stats.ClusterDiversity = zeros(n_interactions, 1);

% Phase 1: Descriptive statistics
for i = 1:n_interactions
    ligand = stats.Ligand{i};
    receptor = stats.Receptor{i};

    lig_matches = strcmp(gbm_network.ligands, ligand);
    rec_matches = strcmp(gbm_network.receptors, receptor);
    interaction_matches = lig_matches & rec_matches;

    if ~any(interaction_matches)
        continue;
    end

    weights = gbm_network.weights(interaction_matches);

    stats.Frequency(i) = sum(interaction_matches);
    stats.MeanStrength(i) = mean(weights);
    stats.MedianStrength(i) = median(weights);
    stats.MaxStrength(i) = max(weights);
    stats.StdStrength(i) = std(weights);

    if stats.Frequency(i) > 1
        stats.CV(i) = stats.StdStrength(i) / (stats.MeanStrength(i) + 1e-10);
        stats.MAD(i) = mad(weights, 1);
    end
end

fprintf('  Basic statistics complete\n\n');

% Phase 2: Bootstrap confidence intervals (100 iterations for speed)
fprintf('Phase 2: Bootstrap confidence intervals (100 iters)...\n');

n_bootstrap = 100;

for i = 1:n_interactions
    ligand = stats.Ligand{i};
    receptor = stats.Receptor{i};

    lig_matches = strcmp(gbm_network.ligands, ligand);
    rec_matches = strcmp(gbm_network.receptors, receptor);
    interaction_matches = lig_matches & rec_matches;

    if ~any(interaction_matches)
        stats.CI_Lower(i) = stats.MeanStrength(i);
        stats.CI_Upper(i) = stats.MeanStrength(i);
        continue;
    end

    weights = gbm_network.weights(interaction_matches);
    n = length(weights);

    if n < 5
        stats.CI_Lower(i) = max(0, stats.MeanStrength(i) - 1.96 * stats.StdStrength(i) / sqrt(n));
        stats.CI_Upper(i) = stats.MeanStrength(i) + 1.96 * stats.StdStrength(i) / sqrt(n);
        continue;
    end

    boot_means = zeros(n_bootstrap, 1);
    for b = 1:n_bootstrap
        boot_sample = weights(randsample(n, n, true));
        boot_means(b) = mean(boot_sample);
    end

    % Percentile method (faster than BCa)
    stats.CI_Lower(i) = prctile(boot_means, 2.5);
    stats.CI_Upper(i) = prctile(boot_means, 97.5);

    if mod(i, 25) == 0
        fprintf('  Progress: %d/%d\n', i, n_interactions);
    end
end

fprintf('  Bootstrap complete\n\n');

% Phase 3: Effect sizes (Cohen's d and Cliff's delta vs. background)
fprintf('Phase 3: Computing effect sizes...\n');

all_weights = gbm_network.weights;
baseline_mean = mean(all_weights);
baseline_std = std(all_weights);

for i = 1:n_interactions
    ligand = stats.Ligand{i};
    receptor = stats.Receptor{i};

    lig_matches = strcmp(gbm_network.ligands, ligand);
    rec_matches = strcmp(gbm_network.receptors, receptor);
    interaction_matches = lig_matches & rec_matches;

    if ~any(interaction_matches)
        continue;
    end

    weights = gbm_network.weights(interaction_matches);

    pooled_std = sqrt((stats.StdStrength(i)^2 + baseline_std^2) / 2);
    if pooled_std > 0
        stats.EffectSize_Cohen(i) = (stats.MeanStrength(i) - baseline_mean) / pooled_std;
    end

    baseline_sample = randsample(all_weights, min(1000, length(all_weights)));
    n1 = length(weights);
    n2 = length(baseline_sample);

    dominance = 0;
    for j = 1:n1
        dominance = dominance + sum(weights(j) > baseline_sample);
        dominance = dominance - sum(weights(j) < baseline_sample);
    end

    stats.EffectSize_Cliff(i) = dominance / (n1 * n2);
end

fprintf('  Effect sizes computed\n\n');

% Phase 4: Spatial heterogeneity and cluster diversity (Shannon entropy)
fprintf('Phase 4: Spatial heterogeneity...\n');

for i = 1:n_interactions
    ligand = stats.Ligand{i};
    receptor = stats.Receptor{i};

    lig_matches = strcmp(gbm_network.ligands, ligand);
    rec_matches = strcmp(gbm_network.receptors, receptor);
    interaction_matches = lig_matches & rec_matches;

    if ~any(interaction_matches)
        continue;
    end

    src_clusters = gbm_network.source_clusters(interaction_matches);
    tgt_clusters = gbm_network.target_clusters(interaction_matches);
    weights = gbm_network.weights(interaction_matches);

    unique_pairs = unique([src_clusters, tgt_clusters], 'rows');
    n_unique_pairs = size(unique_pairs, 1);

    if n_unique_pairs > 1
        pair_counts = zeros(n_unique_pairs, 1);
        for p = 1:n_unique_pairs
            pair_counts(p) = sum(src_clusters == unique_pairs(p,1) & tgt_clusters == unique_pairs(p,2));
        end
        pair_probs = pair_counts / sum(pair_counts);
        pair_probs = pair_probs(pair_probs > 0);
        stats.ClusterDiversity(i) = -sum(pair_probs .* log2(pair_probs));
    end

    if length(weights) > 1
        stats.SpatialHeterogeneity(i) = std(weights) / (mean(weights) + 1e-10);
    end
end

fprintf('  Spatial heterogeneity computed\n\n');

% Phase 5: Brandes betweenness centrality on the gene interaction graph
fprintf('Phase 5: Network centrality...\n');

all_genes = unique([gbm_network.ligands; gbm_network.receptors]);
n_genes = length(all_genes);

adj_matrix = zeros(n_genes, n_genes);
for edge_idx = 1:length(gbm_network.ligands)
    lig_idx = find(strcmp(all_genes, gbm_network.ligands{edge_idx}), 1);
    rec_idx = find(strcmp(all_genes, gbm_network.receptors{edge_idx}), 1);
    if ~isempty(lig_idx) && ~isempty(rec_idx)
        adj_matrix(lig_idx, rec_idx) = gbm_network.weights(edge_idx);
    end
end

betweenness_scores = approximate_betweenness_fast(adj_matrix);

for i = 1:n_interactions
    ligand = stats.Ligand{i};
    receptor = stats.Receptor{i};

    lig_idx = find(strcmp(all_genes, ligand), 1);
    rec_idx = find(strcmp(all_genes, receptor), 1);

    if ~isempty(lig_idx) && ~isempty(rec_idx)
        stats.Betweenness(i) = (betweenness_scores(lig_idx) + betweenness_scores(rec_idx)) / 2;
    end
end

fprintf('  Centrality computed\n\n');

% Phase 5b: Estrada-Hatano subgraph centrality via order-5 matrix exponential.
% SC(v) = diag(e^A)_v, normalized by spectral radius for series convergence.
% Ref: Estrada & Hatano (2008) Phys Rev E 77:036111.
fprintf('Phase 5b: Subgraph centrality (Estrada-Hatano, order-5 series)...\n');

if n_genes <= 2000
    spec_rad = max(abs(eig(adj_matrix)));
    if spec_rad > 0
        A_norm = adj_matrix / spec_rad;
    else
        A_norm = adj_matrix;
    end

    subgraph_centrality = ones(n_genes, 1);   % k=0 term
    A_power = A_norm;
    fact_k = 1;
    for k = 1:5
        fact_k = fact_k * k;
        subgraph_centrality = subgraph_centrality + diag(A_power) / fact_k;
        A_power = A_power * A_norm;
    end

    for i = 1:n_interactions
        lig_idx = find(strcmp(all_genes, stats.Ligand{i}), 1);
        rec_idx = find(strcmp(all_genes, stats.Receptor{i}), 1);
        if ~isempty(lig_idx) && ~isempty(rec_idx)
            stats.Communicability(i) = (subgraph_centrality(lig_idx) + subgraph_centrality(rec_idx)) / 2;
        end
    end
    fprintf('  Subgraph centrality computed (%d genes)\n\n', n_genes);
else
    fprintf('  Network too large for dense eigen-decomp (n=%d) — skipping subgraph centrality\n\n', n_genes);
end

% Phase 6: Composite score (weights sum to 1.0)
fprintf('Phase 6: Composite scores...\n');

stats.FrequencyNorm = normalize_to_range(stats.Frequency);
stats.StrengthNorm = normalize_to_range(stats.MeanStrength);
stats.EffectSizeNorm = normalize_to_range(abs(stats.EffectSize_Cohen));
stats.BetweennessNorm = normalize_to_range(stats.Betweenness);
stats.HeterogeneityNorm = normalize_to_range(stats.SpatialHeterogeneity);
stats.CommunicabilityNorm = normalize_to_range(stats.Communicability);

weights_composite = [0.22, 0.22, 0.20, 0.16, 0.10, 0.10];
stats.CompositeScore = ...
    weights_composite(1) * stats.FrequencyNorm + ...
    weights_composite(2) * stats.StrengthNorm + ...
    weights_composite(3) * stats.EffectSizeNorm + ...
    weights_composite(4) * stats.BetweennessNorm + ...
    weights_composite(5) * (1 - stats.HeterogeneityNorm) + ...
    weights_composite(6) * stats.CommunicabilityNorm;

fprintf('  Composite scores computed\n\n');

fprintf('Statistical analysis complete\n\n');
fprintf('SUMMARY:\n');
fprintf('  Total interactions: %d\n', n_interactions);
fprintf('  Mean frequency: %.1f\n', mean(stats.Frequency));
fprintf('  Mean strength: %.4f\n', mean(stats.MeanStrength));
fprintf('  Mean betweenness: %.4f\n\n', mean(stats.Betweenness));

end

function betweenness = approximate_betweenness_fast(adj_matrix)
% Brandes (2001) betweenness on a random sample of source nodes.
% Ref: Brandes J. Math. Sociol. 25:163-177; approximation per Bader et al. (2007).

n = size(adj_matrix, 1);
betweenness = zeros(n, 1);

n_samples = min(100, n);
sources = randsample(n, n_samples, false);

for s_idx = 1:n_samples
    s = sources(s_idx);

    stack  = zeros(n, 1);   sp_top = 0;
    pred   = cell(n, 1);    for k = 1:n, pred{k} = zeros(1, 0); end
    sigma  = zeros(n, 1);
    dist   = -ones(n, 1);

    sigma(s) = 1;
    dist(s)  = 0;

    queue = zeros(n, 1);
    q_head = 1;  q_tail = 1;
    queue(q_tail) = s;  q_tail = q_tail + 1;

    while q_head < q_tail
        v = queue(q_head);  q_head = q_head + 1;
        sp_top = sp_top + 1;
        stack(sp_top) = v;

        nbrs = find(adj_matrix(v, :) > 0);
        for w = nbrs
            if dist(w) < 0
                queue(q_tail) = w;  q_tail = q_tail + 1;
                dist(w) = dist(v) + 1;
            end
            if dist(w) == dist(v) + 1
                sigma(w) = sigma(w) + sigma(v);
                pred{w}(end+1) = v;
            end
        end
    end

    % Reverse BFS accumulation: delta(w) = sum_t sigma(w)/sigma(t) * (1+delta(t))
    delta = zeros(n, 1);
    while sp_top > 0
        w = stack(sp_top);  sp_top = sp_top - 1;
        for v = pred{w}
            if sigma(w) > 0
                delta(v) = delta(v) + (sigma(v) / sigma(w)) * (1 + delta(w));
            end
        end
        if w ~= s
            betweenness(w) = betweenness(w) + delta(w);
        end
    end
end

% Scale by sampling fraction to keep estimates comparable to full betweenness
betweenness = betweenness * (n / n_samples);
end

function normalized = normalize_to_range(values)
min_val = min(values);
max_val = max(values);

if max_val == min_val
    normalized = ones(size(values));
else
    normalized = (values - min_val) / (max_val - min_val);
end
end
