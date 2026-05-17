function autocrine_results = discover_autocrine_signaling(gbm_network, clusters, expr_matrix, genes_table, gbm_specific)
% DISCOVER_AUTOCRINE_SIGNALING  Identify L-R pairs where the same cluster
%   both produces the ligand and expresses the receptor (self-stimulatory loops).
%
% Usage:
%   autocrine_results = discover_autocrine_signaling(network, clusters, expr, genes, condition_specific)
%
% Inputs:
%   gbm_network  - Communication network struct or table.
%   clusters     - [cells×1] cluster ID array or table with a cluster column.
%   expr_matrix  - [cells×genes] sparse or dense expression matrix.
%   genes_table  - Table or cell array of gene symbols for expr_matrix columns.
%   gbm_specific - Table from identify_gbm_specific_interactions().
%
% Outputs:
%   autocrine_results.autocrine_loops  - Table of co-expressed autocrine pairs
%   autocrine_results.n_loops          - Total autocrine loops found
%   autocrine_results.n_tumor_specific - Loops also in condition_specific table

fprintf('Autocrine Signaling Analysis\n\n');

[ligs, recs, weights, has_weights] = extract_network_data(gbm_network);
fprintf('  Interactions: %d\n', length(ligs));
fprintf('  Has weights: %s\n', string(has_weights));

has_cluster_info = false;
source_clusts = [];
target_clusts = [];

if istable(gbm_network)
    cols = gbm_network.Properties.VariableNames;

    for i = 1:length(cols)
        if contains(lower(cols{i}), {'source_cluster', 'sender_cluster', 'from_cluster'})
            source_clusts = gbm_network.(cols{i});
            has_cluster_info = true;
        end
        if contains(lower(cols{i}), {'target_cluster', 'receiver_cluster', 'to_cluster'})
            target_clusts = gbm_network.(cols{i});
        end
    end
end

fprintf('  Network has cluster info: %s\n\n', string(has_cluster_info));

if istable(genes_table)
    cols = genes_table.Properties.VariableNames;
    gene_col = cols{1};
    for i = 1:length(cols)
        if contains(lower(cols{i}), {'gene', 'symbol', 'name'})
            gene_col = cols{i};
            break;
        end
    end
    gene_names = genes_table.(gene_col);
elseif iscell(genes_table)
    gene_names = genes_table;
else
    gene_names = genes_table;
end

fprintf('  Genes in expression matrix: %d\n', length(gene_names));

if istable(clusters)
    cluster_cols = clusters.Properties.VariableNames;
    cluster_col = '';
    % Prefer exact 'Cluster' column before partial-name search (avoids hitting 'supercluster_term')
    if any(strcmp(cluster_cols, 'Cluster'))
        cluster_col = 'Cluster';
    else
        for i = 1:length(cluster_cols)
            if contains(lower(cluster_cols{i}), {'cluster', 'group', 'type'})
                cluster_col = cluster_cols{i};
                break;
            end
        end
    end
    if isempty(cluster_col)
        warning('Cannot find cluster column - using first column');
        cluster_col = cluster_cols{1};
    end
    cluster_ids = clusters.(cluster_col);
    fprintf('  Cells: %d\n  Clusters: %d unique\n', length(cluster_ids), length(unique(cluster_ids)));
elseif isnumeric(clusters) || iscell(clusters)
    cluster_ids = clusters;
    fprintf('  Cells: %d\n  Clusters: %d unique\n', length(cluster_ids), length(unique(cluster_ids)));
else
    error('Clusters must be table, numeric array, or cell array');
end

fprintf('\n');

% unique() with 'rows' does not work on cell arrays; deduplicate via table
if iscell(ligs)
    pair_table = table(ligs, recs, 'VariableNames', {'Ligand', 'Receptor'});
    unique_pair_table = unique(pair_table, 'rows');
    unique_pairs = [unique_pair_table.Ligand, unique_pair_table.Receptor];
else
    unique_pairs = unique([ligs, recs], 'rows');
end

n_pairs = size(unique_pairs, 1);
fprintf('  Testing %d unique L-R pairs\n', n_pairs);

autocrine_loops = table('Size', [0, 6], ...
    'VariableTypes', {'string', 'string', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Ligand', 'Receptor', 'Cluster', ...
                      'LigandExpr', 'ReceptorExpr', 'CoExprStrength'});

n_test_pairs = min(500, n_pairs);

for pair_idx = 1:n_test_pairs
    lig = unique_pairs{pair_idx, 1};
    rec = unique_pairs{pair_idx, 2};

    lig_idx = find(strcmp(gene_names, lig), 1);
    rec_idx = find(strcmp(gene_names, rec), 1);

    if isempty(lig_idx) || isempty(rec_idx)
        continue;
    end

    for clust = unique(cluster_ids)'
        cluster_mask = cluster_ids == clust;
        if sum(cluster_mask) < 5, continue; end

        lig_expr = mean(full(expr_matrix(cluster_mask, lig_idx)));
        rec_expr = mean(full(expr_matrix(cluster_mask, rec_idx)));

        if lig_expr > 0.1 && rec_expr > 0.1
            coexpr = sqrt(lig_expr * rec_expr);
            autocrine_loops = [autocrine_loops; {lig, rec, clust, lig_expr, rec_expr, coexpr}];
        end
    end

    if mod(pair_idx, 100) == 0
        fprintf('  Progress: %d/%d pairs\n', pair_idx, n_test_pairs);
    end
end

fprintf('  Autocrine candidates: %d\n\n', height(autocrine_loops));

if height(autocrine_loops) == 0
    warning('No autocrine loops detected');
    autocrine_results = struct();
    autocrine_results.autocrine_loops = table();
    autocrine_results.n_loops = 0;
    return;
end

fprintf('Phase 2: Computing ligand-receptor correlations...\n');

autocrine_loops.Correlation = zeros(height(autocrine_loops), 1);

for loop_idx = 1:height(autocrine_loops)
    lig = autocrine_loops.Ligand{loop_idx};
    rec = autocrine_loops.Receptor{loop_idx};
    clust = autocrine_loops.Cluster(loop_idx);

    lig_idx = find(strcmp(gene_names, lig), 1);
    rec_idx = find(strcmp(gene_names, rec), 1);

    if ~isempty(lig_idx) && ~isempty(rec_idx)
        cluster_mask = cluster_ids == clust;
        if sum(cluster_mask) >= 3
            lig_vals = full(expr_matrix(cluster_mask, lig_idx));
            rec_vals = full(expr_matrix(cluster_mask, rec_idx));
            if std(lig_vals) > 0 && std(rec_vals) > 0
                autocrine_loops.Correlation(loop_idx) = corr(lig_vals, rec_vals);
            end
        end
    end
end

fprintf('  Mean correlation: %.3f\n\n', mean(autocrine_loops.Correlation));

if has_cluster_info && ~isempty(source_clusts) && ~isempty(target_clusts)
    fprintf('Phase 3: Computing autocrine/paracrine ratios...\n');

    autocrine_loops.AutocrineScore = zeros(height(autocrine_loops), 1);
    autocrine_loops.ParacrineScore = zeros(height(autocrine_loops), 1);
    autocrine_loops.AutocrineRatio = zeros(height(autocrine_loops), 1);

    for loop_idx = 1:height(autocrine_loops)
        lig = autocrine_loops.Ligand{loop_idx};
        rec = autocrine_loops.Receptor{loop_idx};
        clust = autocrine_loops.Cluster(loop_idx);

        auto_mask = strcmp(ligs, lig) & strcmp(recs, rec) & ...
                    source_clusts == clust & target_clusts == clust;
        auto_score = sum(weights(auto_mask));

        para_mask = strcmp(ligs, lig) & strcmp(recs, rec) & ...
                    source_clusts ~= target_clusts;
        para_score = sum(weights(para_mask));

        autocrine_loops.AutocrineScore(loop_idx) = auto_score;
        autocrine_loops.ParacrineScore(loop_idx) = para_score;

        total = auto_score + para_score;
        if total > 0
            autocrine_loops.AutocrineRatio(loop_idx) = auto_score / total;
        end
    end

    enriched_mask = autocrine_loops.AutocrineRatio > 0.6;
    autocrine_loops_enriched = autocrine_loops(enriched_mask, :);
    fprintf('  Autocrine-enriched (ratio>0.6): %d\n\n', height(autocrine_loops_enriched));
else
    fprintf('Phase 3: Skipping autocrine/paracrine (no cluster info in network)\n\n');
    autocrine_loops_enriched = autocrine_loops;
    autocrine_loops_enriched.AutocrineScore = autocrine_loops.CoExprStrength;
    autocrine_loops_enriched.ParacrineScore = zeros(height(autocrine_loops), 1);
    autocrine_loops_enriched.AutocrineRatio = ones(height(autocrine_loops), 1);
end

if height(autocrine_loops_enriched) == 0
    warning('No enriched loops after filtering');
    autocrine_results = struct();
    autocrine_results.autocrine_loops = autocrine_loops;
    autocrine_results.n_loops = height(autocrine_loops);
    return;
end

fprintf('Phase 4: Computing relative autonomy scores...\n');

% Clamp Pearson correlation to [0,1]: negative correlation is inconsistent
% with autocrine signaling and would invert the autonomy score.
corr_vals = autocrine_loops_enriched.Correlation;
autocrine_loops_enriched.AutonomyScore = ...
    max(0, corr_vals) .* autocrine_loops_enriched.AutocrineRatio;

max_autonomy = max(autocrine_loops_enriched.AutonomyScore);
if max_autonomy > 0
    autocrine_loops_enriched.AutonomyScore = ...
        autocrine_loops_enriched.AutonomyScore / max_autonomy;
end

fprintf('  Mean autonomy score: %.3f\n', mean(autocrine_loops_enriched.AutonomyScore));
fprintf('  High autonomy (>0.7): %d loops\n\n', ...
    sum(autocrine_loops_enriched.AutonomyScore > 0.7));

fprintf('Phase 5: Checking specificity...\n');

autocrine_loops_enriched.TumorSpecific = false(height(autocrine_loops_enriched), 1);

if istable(gbm_specific) && height(gbm_specific) > 0
    spec_cols = gbm_specific.Properties.VariableNames;
    lig_col = '';
    rec_col = '';

    for i = 1:length(spec_cols)
        if contains(lower(spec_cols{i}), {'ligand', 'source'})
            lig_col = spec_cols{i};
        end
        if contains(lower(spec_cols{i}), {'receptor', 'target'})
            rec_col = spec_cols{i};
        end
    end

    if ~isempty(lig_col) && ~isempty(rec_col)
        fprintf('  Using: %s -> %s\n', lig_col, rec_col);
        for loop_idx = 1:height(autocrine_loops_enriched)
            lig = autocrine_loops_enriched.Ligand{loop_idx};
            rec = autocrine_loops_enriched.Receptor{loop_idx};
            is_spec = any(strcmp(gbm_specific.(lig_col), lig) & ...
                         strcmp(gbm_specific.(rec_col), rec));
            autocrine_loops_enriched.TumorSpecific(loop_idx) = is_spec;
        end
        fprintf('  Tumor-specific: %d\n', sum(autocrine_loops_enriched.TumorSpecific));
    else
        fprintf('  Warning: Cannot find ligand/receptor columns in gbm_specific\n');
    end
else
    fprintf('  No specificity data provided\n');
end

fprintf('\nPhase 6: Prioritizing by autonomy and specificity...\n');

autocrine_loops_enriched.PriorityScore = ...
    autocrine_loops_enriched.AutonomyScore .* ...
    (1 + 2*double(autocrine_loops_enriched.TumorSpecific));

autocrine_loops_enriched = sortrows(autocrine_loops_enriched, 'PriorityScore', 'descend');

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

writetable(autocrine_loops_enriched, 'results/novel_discoveries/autocrine_loops.csv');
fprintf('  autocrine_loops.csv\n\n');

autocrine_results = struct();
autocrine_results.autocrine_loops = autocrine_loops_enriched;
autocrine_results.n_loops = height(autocrine_loops_enriched);
autocrine_results.n_tumor_specific = sum(autocrine_loops_enriched.TumorSpecific);
autocrine_results.has_cluster_info = has_cluster_info;

fprintf('RESULTS:\n');
fprintf('  Autocrine candidates: %d\n', height(autocrine_loops));
fprintf('  Enriched loops: %d\n', height(autocrine_loops_enriched));
fprintf('  Tumor-specific: %d\n', sum(autocrine_loops_enriched.TumorSpecific));

if height(autocrine_loops_enriched) > 0
    fprintf('\nTOP 5 AUTOCRINE LOOPS:\n');
    for i = 1:min(5, height(autocrine_loops_enriched))
        fprintf('  %d. %s -> %s', i, ...
            autocrine_loops_enriched.Ligand{i}, ...
            autocrine_loops_enriched.Receptor{i});
        if autocrine_loops_enriched.TumorSpecific(i)
            fprintf(' [TUMOR-SPECIFIC]');
        end
        fprintf('\n');
        fprintf('     Cluster: %d | Co-expr: %.3f | Autonomy: %.3f\n', ...
            autocrine_loops_enriched.Cluster(i), ...
            autocrine_loops_enriched.CoExprStrength(i), ...
            autocrine_loops_enriched.AutonomyScore(i));
    end
end

fprintf('\nAnalysis complete\n\n');

end

function [ligs, recs, weights, has_weights] = extract_network_data(network)
% Extract ligands, receptors, and weights from a struct or table network

if istable(network)
    cols = network.Properties.VariableNames;

    lig_col = '';
    for i = 1:length(cols)
        if contains(lower(cols{i}), {'ligand', 'source', 'from'})
            lig_col = cols{i};
            break;
        end
    end

    rec_col = '';
    for i = 1:length(cols)
        if contains(lower(cols{i}), {'receptor', 'target', 'to'})
            rec_col = cols{i};
            break;
        end
    end

    if isempty(lig_col) || isempty(rec_col)
        error('Cannot find ligand/receptor columns');
    end

    ligs = network.(lig_col);
    recs = network.(rec_col);

    weight_col = '';
    for i = 1:length(cols)
        if contains(lower(cols{i}), {'weight', 'score', 'strength'})
            weight_col = cols{i};
            break;
        end
    end

    if ~isempty(weight_col)
        weights = network.(weight_col);
        has_weights = true;
    else
        weights = ones(height(network), 1);
        has_weights = false;
    end

elseif isstruct(network)
    fields = fieldnames(network);

    lig_field = '';
    for i = 1:length(fields)
        if contains(lower(fields{i}), {'ligand', 'source', 'from'})
            lig_field = fields{i};
            break;
        end
    end

    rec_field = '';
    for i = 1:length(fields)
        if contains(lower(fields{i}), {'receptor', 'target', 'to'})
            rec_field = fields{i};
            break;
        end
    end

    if isempty(lig_field) || isempty(rec_field)
        error('Cannot find ligand/receptor fields');
    end

    ligs = network.(lig_field);
    recs = network.(rec_field);

    if isfield(network, 'weights')
        weights = network.weights;
        has_weights = true;
    elseif isfield(network, 'weight')
        weights = network.weight;
        has_weights = true;
    else
        weights = ones(size(ligs));
        has_weights = false;
    end
else
    error('Network must be table or struct');
end
end
