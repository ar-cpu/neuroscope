function autocrine_results = discover_autocrine_signaling(gbm_network, clusters, expr_matrix, genes_table, gbm_specific)
% DISCOVER_AUTOCRINE_SIGNALING  Identify ligand-receptor pairs where the
%   same cell cluster both produces the ligand and expresses the receptor,
%   indicative of self-stimulatory (autocrine) signaling loops.
%
% Autocrine signaling is a hallmark of cancer cell autonomy: cells secrete
% factors that act on themselves rather than requiring extrinsic stimulation.
% This function quantifies autocrine activity at two levels:
%
%   Interaction-level — from the network, which L-R pairs route from a
%     cluster back to itself? These are directly observed autocrine edges.
%   Expression-level — for each cluster, which L-R pairs have both the
%     ligand and receptor co-expressed at appreciable levels (mean > threshold)?
%     This captures potential autocrine loops even if the network graph edge
%     was not inferred (e.g., because permutation p-value did not pass FDR).
%
% Autonomy scoring uses the formulation:
%   AutocrineScore = max(0, Correlation) × AutocrineRatio
% where Correlation is the Pearson correlation between ligand and receptor
% expression across cells in the cluster, and AutocrineRatio is the fraction
% of total L-R interactions that are self-directed.
%
% Limitations (explicit):
%   - Autocrine loop presence does not imply autonomous growth; tumor growth
%     requires additional permissive oncogenic events.
%   - The Pearson correlation is a proxy for co-regulation, not physical
%     interaction between secreted ligand and cell-surface receptor.
%   - Scores are relative rankings within this dataset, not calibrated
%     against experimental autocrine assays.
%
% Usage:
%   autocrine_results = discover_autocrine_signaling(network, clusters, expr, genes, condition_specific)
%
% Inputs:
%   gbm_network  - Communication network struct or table.
%   clusters     - [cells×1] cluster ID array or table with a cluster column.
%   expr_matrix  - [cells×genes] sparse or dense expression matrix.
%   genes_table  - Table or cell array of gene symbols for expr_matrix columns.
%   gbm_specific - Table of condition-specific interactions from
%                  identify_gbm_specific_interactions(). Used to annotate
%                  whether autocrine loops are condition-enriched.
%
% Outputs:
%   autocrine_results.autocrine_loops  - Table of co-expressed autocrine pairs
%   autocrine_results.n_loops          - Total number of autocrine loops found
%   autocrine_results.n_significant    - Loops with FDR < 0.05
%   autocrine_results.n_condition_specific - Loops also in condition_specific table

fprintf('Autocrine Signaling Analysis\n\n');

%% Extract network data (handles struct and table formats)

[ligs, recs, weights, has_weights] = extract_network_data(gbm_network);
fprintf('  Interactions: %d\n', length(ligs));
fprintf('  Has weights: %s\n', string(has_weights));

% Check for cluster assignments in network
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

fprintf('  Network has cluster info: %s\n', string(has_cluster_info));
fprintf('\n');

%% EXTRACT GENE NAMES
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

%% VALIDATE CLUSTERS
if istable(clusters)
    cluster_cols = clusters.Properties.VariableNames;
    cluster_col = '';

    % Prefer the standardised 'Cluster' column added by NEUROSCOPE before
    % falling back to a partial-name search (which can hit 'supercluster_term').
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
    fprintf('  Cells: %d\n', length(cluster_ids));
    fprintf('  Clusters: %d unique\n', length(unique(cluster_ids)));
elseif isnumeric(clusters) || iscell(clusters)
    cluster_ids = clusters;
    fprintf('  Cells: %d\n', length(cluster_ids));
    fprintf('  Clusters: %d unique\n', length(unique(cluster_ids)));
else
    error('Clusters must be table, numeric array, or cell array');
end

fprintf('\n');

%% Phase 1: Identify autocrine L-R candidates

% unique() with 'rows' does not operate on cell arrays; use string
% concatenation for deduplication of L-R pair signatures.
if iscell(ligs)
    % Convert to table for unique rows
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

n_clusters = length(unique(cluster_ids));
n_test_pairs = min(500, n_pairs);

for pair_idx = 1:n_test_pairs
    lig = unique_pairs{pair_idx, 1};
    rec = unique_pairs{pair_idx, 2};
    
    lig_idx = find(strcmp(gene_names, lig), 1);
    rec_idx = find(strcmp(gene_names, rec), 1);
    
    if isempty(lig_idx) || isempty(rec_idx)
        continue;
    end
    
    % Check each cluster
    for clust = unique(cluster_ids)'
        cluster_mask = cluster_ids == clust;
        
        if sum(cluster_mask) < 5
            continue;
        end
        
        lig_expr = mean(full(expr_matrix(cluster_mask, lig_idx)));
        rec_expr = mean(full(expr_matrix(cluster_mask, rec_idx)));
        
        % Both expressed
        if lig_expr > 0.1 && rec_expr > 0.1
            coexpr = sqrt(lig_expr * rec_expr);
            
            new_row = {lig, rec, clust, lig_expr, rec_expr, coexpr};
            autocrine_loops = [autocrine_loops; new_row];
        end
    end
    
    if mod(pair_idx, 100) == 0
        fprintf('  Progress: %d/%d pairs\n', pair_idx, n_test_pairs);
    end
end

fprintf('  Autocrine candidates: %d\n', height(autocrine_loops));
fprintf('\n');

if height(autocrine_loops) == 0
    warning('No autocrine loops detected');
    autocrine_results = struct();
    autocrine_results.autocrine_loops = table();
    autocrine_results.n_loops = 0;
    return;
end

%% Phase 2: Ligand-receptor co-expression correlation
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

fprintf('  Mean correlation: %.3f\n', mean(autocrine_loops.Correlation));
fprintf('\n');

%% Phase 3: Autocrine vs paracrine ratio (requires cluster info)
if has_cluster_info && ~isempty(source_clusts) && ~isempty(target_clusts)
    fprintf('Phase 3: Computing autocrine/paracrine ratios...\n');
    
    autocrine_loops.AutocrineScore = zeros(height(autocrine_loops), 1);
    autocrine_loops.ParacrineScore = zeros(height(autocrine_loops), 1);
    autocrine_loops.AutocrineRatio = zeros(height(autocrine_loops), 1);
    
    for loop_idx = 1:height(autocrine_loops)
        lig = autocrine_loops.Ligand{loop_idx};
        rec = autocrine_loops.Receptor{loop_idx};
        clust = autocrine_loops.Cluster(loop_idx);
        
        % Autocrine: same cluster
        auto_mask = strcmp(ligs, lig) & strcmp(recs, rec) & ...
                    source_clusts == clust & target_clusts == clust;
        auto_score = sum(weights(auto_mask));
        
        % Paracrine: different clusters
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
    
    % Filter for autocrine-enriched
    enriched_mask = autocrine_loops.AutocrineRatio > 0.6;
    autocrine_loops_enriched = autocrine_loops(enriched_mask, :);
    
    fprintf('  Autocrine-enriched (ratio>0.6): %d\n', height(autocrine_loops_enriched));
    fprintf('\n');
else
    fprintf('Phase 3: Skipping autocrine/paracrine (no cluster info in network)\n');
    fprintf('  Using all co-expressed pairs as candidates\n\n');
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

%% Phase 4: Relative autocrine autonomy scoring
fprintf('Phase 4: Computing relative autonomy scores...\n');

% Autonomy = within-cluster co-expression correlation × autocrine fraction.
% We use the Pearson correlation computed in Phase 2 (cells within the
% same cluster) as the co-expression component, clamped to [0, 1] because
% negative correlation means the genes are anti-correlated within a cell -
% inconsistent with autocrine signaling.  The autocrine fraction (ratio of
% self-loop signal weight to total signal weight) weights that correlation
% by how exclusively the pair operates within the same cluster.
corr_vals = autocrine_loops_enriched.Correlation;
autocrine_loops_enriched.AutonomyScore = ...
    max(0, corr_vals) .* autocrine_loops_enriched.AutocrineRatio;

% Normalize
max_autonomy = max(autocrine_loops_enriched.AutonomyScore);
if max_autonomy > 0
    autocrine_loops_enriched.AutonomyScore = ...
        autocrine_loops_enriched.AutonomyScore / max_autonomy;
end

fprintf('  Mean autonomy score: %.3f\n', mean(autocrine_loops_enriched.AutonomyScore));
fprintf('  High autonomy (>0.7): %d loops\n', ...
    sum(autocrine_loops_enriched.AutonomyScore > 0.7));
fprintf('\n');

%% Phase 5: Cross-reference with condition-specific interactions
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

fprintf('\n');

%% Phase 6: Prioritization and export
fprintf('Phase 6: Prioritizing by autonomy and specificity...\n');

% Priority = autonomy × (1 + tumor_specific)
autocrine_loops_enriched.PriorityScore = ...
    autocrine_loops_enriched.AutonomyScore .* ...
    (1 + 2*double(autocrine_loops_enriched.TumorSpecific));

autocrine_loops_enriched = sortrows(autocrine_loops_enriched, 'PriorityScore', 'descend');

fprintf('  Prioritization complete\n\n');

%% EXPORT
fprintf('Exporting results...\n');

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

writetable(autocrine_loops_enriched, 'results/novel_discoveries/autocrine_loops.csv');
fprintf('  autocrine_loops.csv\n\n');

%% PACKAGE RESULTS
autocrine_results = struct();
autocrine_results.autocrine_loops = autocrine_loops_enriched;
autocrine_results.n_loops = height(autocrine_loops_enriched);
autocrine_results.n_tumor_specific = sum(autocrine_loops_enriched.TumorSpecific);
autocrine_results.has_cluster_info = has_cluster_info;

%% SUMMARY
fprintf('Autocrine signaling analysis complete\n\n');

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

%% HELPER FUNCTIONS

function [ligs, recs, weights, has_weights] = extract_network_data(network)
% Extract ligands, receptors, and weights from network

if istable(network)
    cols = network.Properties.VariableNames;
    
    % Find ligand column
    lig_col = '';
    for i = 1:length(cols)
        if contains(lower(cols{i}), {'ligand', 'source', 'from'})
            lig_col = cols{i};
            break;
        end
    end
    
    % Find receptor column
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
    
    % Find weight column
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
    
    % Find ligand field
    lig_field = '';
    for i = 1:length(fields)
        if contains(lower(fields{i}), {'ligand', 'source', 'from'})
            lig_field = fields{i};
            break;
        end
    end
    
    % Find receptor field
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
    
    % Find weight field
    if isfield(network, 'weights') || isfield(network, 'weight')
        if isfield(network, 'weights')
            weights = network.weights;
        else
            weights = network.weight;
        end
        has_weights = true;
    else
        weights = ones(size(ligs));
        has_weights = false;
    end
else
    error('Network must be table or struct');
end
end