function niche_results = discover_spatial_communication_niches(gbm_network, clusters, expr_matrix, genes_table, niche_config)
% DISCOVER_SPATIAL_COMMUNICATION_NICHES  Assign cells to functional niches
%   based on marker gene expression profiles, then test which ligand-receptor
%   interactions are statistically enriched within each niche.
%
% Niche assignment
% ----------------
% Each cell is scored against each niche by averaging the expression of the
% niche's marker genes. Cells are assigned to their highest-scoring niche
% (argmax). A plasticity index — the normalized Shannon entropy over the
% niche score distribution — identifies cells that simultaneously express
% markers from multiple niches, indicating transitional or plastic states.
% Boundary cells (plasticity > 0.5 or second-best score > 0.3) are flagged
% for potential inter-niche transition analysis.
%
% Enrichment testing strategy
% ---------------------------
%   If the network contains cluster assignments:
%     Binomial test — the proportion of a given L-R pair's occurrences that
%     fall within the niche is compared to the expected proportion under the
%     null that interactions are distributed uniformly across clusters.
%   Otherwise:
%     Wilcoxon rank-sum test on ligand and receptor expression (niche vs.
%     rest), combined via Fisher's method (chi-squared with 4 df) for a
%     single per-pair p-value.
%
% Custom niche markers
% --------------------
% Pass a niche_config struct to override the default solid-tumor markers.
% niche_config.niche_markers must be a struct where each field is a niche
% name and its value is a cell array of marker gene symbols.
%
% Usage:
%   niche_results = discover_spatial_communication_niches(network, clusters, expr, genes)
%   niche_results = discover_spatial_communication_niches(..., niche_config)
%
% Inputs:
%   gbm_network  - Communication network struct or table.
%   clusters     - [cells×1] cluster ID array, or table with a column
%                  whose name contains 'cluster', 'group', or 'type'.
%   expr_matrix  - [cells×genes] sparse or dense expression matrix.
%   genes_table  - Table (first or 'gene'/'symbol'/'name' column) or cell
%                  array of gene symbols corresponding to expr_matrix columns.
%   niche_config - (optional) Struct with field niche_markers.(niche_name).
%                  Defaults to four solid-tumor functional states:
%                  Hypoxic, Proliferating, Mesenchymal, Immune.
%
% Outputs:
%   niche_results.niche_assignments - {cells×1} assigned niche per cell
%   niche_results.plasticity        - [cells×1] normalized entropy score
%   niche_results.is_boundary       - [cells×1] logical boundary-cell flag
%   niche_results.transition_counts - [niches×niches] boundary-cell matrix
%   niche_results.niche_enrichment  - Enriched L-R pairs per niche (FDR-corrected)
%   niche_results.niche_summary     - Per-niche cell counts and enrichment stats

fprintf('Spatial communication niche analysis\n\n');

%% Handle optional niche config

if nargin < 5 || isempty(niche_config)
    % Default markers cover major functional states in solid tumors.
    % These work well for GBM, breast, lung, and other solid tumors.
    % For brain-specific or non-cancer datasets, pass custom niche_config.
    niche_config.niche_markers.Hypoxic      = {'HIF1A', 'VEGFA', 'CA9', 'LDHA', 'SLC2A1', 'PDK1', 'BNIP3'};
    niche_config.niche_markers.Proliferating = {'MKI67', 'TOP2A', 'PCNA', 'CDK1', 'CCNB1', 'MCM2'};
    niche_config.niche_markers.Mesenchymal  = {'VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'TWIST1', 'MMP2', 'CDH2'};
    niche_config.niche_markers.Immune       = {'PTPRC', 'CD3D', 'CD68', 'CD14', 'MS4A1', 'NKG7'};
    fprintf('  Using default solid-tumor niche markers\n');
    fprintf('  (pass niche_config for custom markers)\n\n');
else
    fprintf('  Using custom niche markers\n\n');
end

%% Extract gene names

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
    gene_names = cellstr(string(genes_table));
end

%% Extract cluster assignments

if istable(clusters)
    cluster_cols = clusters.Properties.VariableNames;
    cluster_col = cluster_cols{1};
    for i = 1:length(cluster_cols)
        if contains(lower(cluster_cols{i}), {'cluster', 'group', 'type'})
            cluster_col = cluster_cols{i};
            break;
        end
    end
    cluster_ids = clusters.(cluster_col);
elseif isnumeric(clusters) || iscell(clusters)
    cluster_ids = clusters;
else
    error('Clusters must be a table or array');
end

n_cells = length(cluster_ids);
fprintf('  Cells: %d | Clusters: %d unique\n\n', n_cells, length(unique(cluster_ids)));

%% Phase 1: Score cells for each niche

fprintf('Phase 1: Scoring cells for niche membership...\n');

niche_types = fieldnames(niche_config.niche_markers);
n_niches = length(niche_types);
niche_scores = zeros(n_cells, n_niches);

for i = 1:n_niches
    niche_name = niche_types{i};
    markers = niche_config.niche_markers.(niche_name);

    marker_idx = [];
    for j = 1:length(markers)
        idx = find(strcmp(gene_names, markers{j}), 1);
        if ~isempty(idx)
            marker_idx = [marker_idx; idx];
        end
    end

    found = length(marker_idx);
    total = length(markers);
    fprintf('  %s: %d/%d markers found\n', niche_name, found, total);

    if found > 0
        niche_scores(:, i) = mean(full(expr_matrix(:, marker_idx)), 2);
    else
        fprintf('    Warning: no markers found for %s. Scores will be zero.\n', niche_name);
    end
end

% Assign each cell to its top-scoring niche
[max_scores, niche_assignments_idx] = max(niche_scores, [], 2);
niche_assignments = niche_types(niche_assignments_idx);

fprintf('\nNiche assignments:\n');
for i = 1:n_niches
    n_assigned = sum(niche_assignments_idx == i);
    fprintf('  %s: %d (%.1f%%)\n', niche_types{i}, n_assigned, 100*n_assigned/n_cells);
end

%% Phase 2: Cell plasticity (entropy over niche scores)

fprintf('\nPhase 2: Computing plasticity...\n');

row_sums = sum(niche_scores, 2) + 1e-10;
niche_scores_norm = niche_scores ./ row_sums;
plasticity = -sum(niche_scores_norm .* log(niche_scores_norm + 1e-10), 2);
plasticity = plasticity / log(n_niches);

fprintf('  Mean plasticity: %.3f (0=pure, 1=ambiguous)\n', mean(plasticity));
fprintf('  High-plasticity cells (>0.7): %d (%.1f%%)\n', ...
    sum(plasticity > 0.7), 100*mean(plasticity > 0.7));

%% Phase 3: Boundary cells

fprintf('\nPhase 3: Identifying boundary cells...\n');

second_best_scores = zeros(n_cells, 1);
for i = 1:n_cells
    sorted_scores = sort(niche_scores(i, :), 'descend');
    if length(sorted_scores) >= 2
        second_best_scores(i) = sorted_scores(2);
    end
end

is_boundary = (plasticity > 0.5) | (second_best_scores > 0.3);
fprintf('  Boundary cells: %d (%.1f%%)\n', sum(is_boundary), 100*mean(is_boundary));

%% Phase 4: Niche transition matrix

fprintf('\nPhase 4: Estimating niche transition potential...\n');

transition_counts = zeros(n_niches, n_niches);

for i = 1:n_cells
    if is_boundary(i)
        from_niche = niche_assignments_idx(i);
        [sorted_scores, sorted_idx] = sort(niche_scores(i, :), 'descend');
        if length(sorted_scores) >= 2 && sorted_scores(2) > 0.2
            to_niche = sorted_idx(2);
            transition_counts(from_niche, to_niche) = transition_counts(from_niche, to_niche) + 1;
        end
    end
end

fprintf('  Non-zero transitions:\n');
for i = 1:n_niches
    for j = 1:n_niches
        if i ~= j && transition_counts(i, j) > 0
            total_bidirectional = transition_counts(i,j) + transition_counts(j,i);
            fprintf('    %s <-> %s: %d boundary cells\n', ...
                niche_types{i}, niche_types{j}, total_bidirectional);
        end
    end
end

%% Phase 5: Niche-enriched L-R interactions

fprintf('\nPhase 5: Testing niche-enriched interactions...\n');

% Extract network data
[ligs, recs, source_clusters, target_clusters, has_cluster_info] = extract_network_data(gbm_network);

fprintf('  Network interactions: %d\n', length(ligs));

niche_enrichment = table('Size', [0, 7], ...
    'VariableTypes', {'string', 'string', 'string', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Ligand', 'Receptor', 'Niche', 'ObservedCount', 'ExpectedCount', 'FoldEnrichment', 'PValue'});

if has_cluster_info
    fprintf('  Method: cluster-based enrichment (binomial test)\n');
    niche_enrichment = test_cluster_enrichment(ligs, recs, source_clusters, target_clusters, ...
        niche_assignments, niche_assignments_idx, cluster_ids, niche_types, niche_enrichment);
else
    fprintf('  Method: expression-based enrichment (Wilcoxon + Fisher)\n');
    niche_enrichment = test_expression_enrichment(ligs, recs, niche_assignments, ...
        niche_assignments_idx, niche_types, expr_matrix, gene_names, niche_enrichment);
end

if height(niche_enrichment) > 0
    [~, ~, ~, fdr] = fdr_bh(niche_enrichment.PValue);
    niche_enrichment.FDR = fdr;
    niche_enrichment = sortrows(niche_enrichment, 'FoldEnrichment', 'descend');

    fprintf('\n  Niche-enriched pairs: %d\n', height(niche_enrichment));
    fprintf('  Significant (FDR<0.05): %d\n', sum(niche_enrichment.FDR < 0.05));
else
    fprintf('\n  No enrichments detected\n');
    niche_enrichment.FDR = zeros(0, 1);
end

%% Export

fprintf('\nExporting results...\n');

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

plasticity_table = table(niche_assignments, plasticity, is_boundary, ...
    'VariableNames', {'Niche', 'Plasticity', 'IsBoundary'});
writetable(plasticity_table, 'results/novel_discoveries/niche_plasticity.csv');

transition_table = array2table(transition_counts, ...
    'RowNames', niche_types, 'VariableNames', niche_types);
writetable(transition_table, 'results/novel_discoveries/niche_transitions.csv', ...
    'WriteRowNames', true);

if height(niche_enrichment) > 0
    writetable(niche_enrichment, 'results/novel_discoveries/spatial_niche_enrichment.csv');
end

niche_summary = table();
for i = 1:n_niches
    mask = strcmp(niche_assignments, niche_types{i});
    n_cells_niche = sum(mask);
    n_boundary_niche = sum(is_boundary & mask);
    mean_plast = mean(plasticity(mask));
    n_enriched = 0;
    if height(niche_enrichment) > 0
        n_enriched = sum(strcmp(niche_enrichment.Niche, niche_types{i}) & ...
            niche_enrichment.FDR < 0.05);
    end
    niche_summary = [niche_summary; {niche_types{i}, n_cells_niche, n_boundary_niche, mean_plast, n_enriched}];
end
niche_summary.Properties.VariableNames = {'Niche', 'NCells', 'NBoundary', 'MeanPlasticity', 'NEnriched'};
writetable(niche_summary, 'results/novel_discoveries/spatial_niche_summary.csv');

%% Package and summarize

niche_results = struct();
niche_results.niche_assignments = niche_assignments;
niche_results.plasticity = plasticity;
niche_results.is_boundary = is_boundary;
niche_results.transition_counts = transition_counts;
niche_results.niche_enrichment = niche_enrichment;
niche_results.niche_summary = niche_summary;

fprintf('\nSpatial niche analysis complete\n');
fprintf('  Niches: %d | Boundary cells: %d (%.1f%%)\n', ...
    n_niches, sum(is_boundary), 100*mean(is_boundary));
if height(niche_enrichment) > 0
    fprintf('  Enriched interactions: %d | Significant: %d\n', ...
        height(niche_enrichment), sum(niche_enrichment.FDR < 0.05));

    fprintf('\nTop 5 niche-enriched interactions:\n');
    for i = 1:min(5, height(niche_enrichment))
        fprintf('  %s -> %s (%s): FE=%.2f, FDR=%.2e\n', ...
            niche_enrichment.Ligand{i}, niche_enrichment.Receptor{i}, ...
            niche_enrichment.Niche{i}, ...
            niche_enrichment.FoldEnrichment(i), niche_enrichment.FDR(i));
    end
end
fprintf('\n');

end

%% Local helper functions

function [ligs, recs, source_clusters, target_clusters, has_cluster_info] = extract_network_data(network)
has_cluster_info = false;
source_clusters = [];
target_clusters = [];

if istable(network)
    cols = network.Properties.VariableNames;
    lig_col = find_col(cols, {'ligand', 'source', 'from'}, {'cluster'});
    rec_col = find_col(cols, {'receptor', 'target', 'to'}, {'cluster'});
    if isempty(lig_col) || isempty(rec_col)
        error('Cannot identify ligand/receptor columns in network table');
    end
    ligs = network.(lig_col);
    recs = network.(rec_col);

    src_col = find_col(cols, {'source_cluster', 'sender_cluster', 'from_cluster'}, {});
    tgt_col = find_col(cols, {'target_cluster', 'receiver_cluster', 'to_cluster'}, {});
    if ~isempty(src_col) && ~isempty(tgt_col)
        source_clusters = network.(src_col);
        target_clusters = network.(tgt_col);
        has_cluster_info = true;
    end

elseif isstruct(network)
    fields = fieldnames(network);
    lig_field = find_col(fields, {'ligand', 'source'}, {'cluster'});
    rec_field = find_col(fields, {'receptor', 'target'}, {'cluster'});
    if isempty(lig_field) || isempty(rec_field)
        error('Cannot identify ligand/receptor fields in network struct');
    end
    ligs = network.(lig_field);
    recs = network.(rec_field);

    if isfield(network, 'source_clusters') && isfield(network, 'target_clusters')
        source_clusters = network.source_clusters;
        target_clusters = network.target_clusters;
        has_cluster_info = true;
    end
else
    error('Network must be a table or struct');
end
end

function col = find_col(cols, include_terms, exclude_terms)
col = '';
for i = 1:length(cols)
    col_lower = lower(cols{i});
    matches_include = any(cellfun(@(t) contains(col_lower, t), include_terms));
    matches_exclude = ~isempty(exclude_terms) && any(cellfun(@(t) contains(col_lower, t), exclude_terms));
    if matches_include && ~matches_exclude
        col = cols{i};
        return;
    end
end
end

function niche_enrichment = test_cluster_enrichment(ligs, recs, source_clusters, target_clusters, ...
    niche_assignments, niche_assignments_idx, cluster_ids, niche_types, niche_enrichment)

if iscell(source_clusters)
    source_clusters = cellfun(@str2double, source_clusters);
end
if iscell(target_clusters)
    target_clusters = cellfun(@str2double, target_clusters);
end

for niche_idx = 1:length(niche_types)
    niche = niche_types{niche_idx};
    niche_cell_mask = strcmp(niche_assignments, niche);
    niche_cluster_ids = unique(cluster_ids(niche_cell_mask));

    if iscell(niche_cluster_ids)
        niche_cluster_ids = cellfun(@str2double, niche_cluster_ids);
    end

    niche_mask = false(length(ligs), 1);
    for cid = niche_cluster_ids'
        if ~isnan(cid)
            niche_mask = niche_mask | (source_clusters == cid) | (target_clusters == cid);
        end
    end

    if sum(niche_mask) == 0, continue; end

    pair_table = table(ligs(niche_mask), recs(niche_mask), 'VariableNames', {'Ligand', 'Receptor'});
    unique_pairs = unique(pair_table, 'rows');

    for pair_idx = 1:height(unique_pairs)
        lig = unique_pairs.Ligand{pair_idx};
        rec = unique_pairs.Receptor{pair_idx};

        observed = sum(strcmp(ligs, lig) & strcmp(recs, rec) & niche_mask);
        total_pairs = sum(strcmp(ligs, lig) & strcmp(recs, rec));
        niche_proportion = sum(niche_mask) / length(ligs);
        expected = total_pairs * niche_proportion;

        if expected > 0 && observed >= 3
            fold_enrichment = observed / expected;
            if fold_enrichment > 1.5
                p_val = 1 - binocdf(observed - 1, total_pairs, niche_proportion);
                niche_enrichment = [niche_enrichment; {lig, rec, niche, observed, expected, fold_enrichment, p_val}];
            end
        end
    end
end
end

function niche_enrichment = test_expression_enrichment(ligs, recs, niche_assignments, ...
    niche_assignments_idx, niche_types, expr_matrix, gene_names, niche_enrichment)

max_test = 500;

pair_table = table(ligs, recs, 'VariableNames', {'Ligand', 'Receptor'});
unique_pairs = unique(pair_table, 'rows');
n_test = min(max_test, height(unique_pairs));

for niche_idx = 1:length(niche_types)
    niche = niche_types{niche_idx};
    niche_cells = strcmp(niche_assignments, niche);

    if sum(niche_cells) < 10, continue; end

    for pair_idx = 1:n_test
        lig = unique_pairs.Ligand{pair_idx};
        rec = unique_pairs.Receptor{pair_idx};

        lig_idx = find(strcmp(gene_names, lig), 1);
        rec_idx = find(strcmp(gene_names, rec), 1);
        if isempty(lig_idx) || isempty(rec_idx), continue; end

        lig_expr_niche = mean(full(expr_matrix(niche_cells, lig_idx)));
        rec_expr_niche = mean(full(expr_matrix(niche_cells, rec_idx)));
        lig_expr_other = mean(full(expr_matrix(~niche_cells, lig_idx)));
        rec_expr_other = mean(full(expr_matrix(~niche_cells, rec_idx)));

        if lig_expr_niche > 0.1 && rec_expr_niche > 0.1
            lig_fc = (lig_expr_niche + 0.01) / (lig_expr_other + 0.01);
            rec_fc = (rec_expr_niche + 0.01) / (rec_expr_other + 0.01);
            fold_enrichment = sqrt(lig_fc * rec_fc);

            if fold_enrichment > 1.5
                lig_vals_niche = full(expr_matrix(niche_cells, lig_idx));
                lig_vals_other = full(expr_matrix(~niche_cells, lig_idx));

                try
                    p_lig = ranksum(lig_vals_niche, lig_vals_other);
                catch
                    p_lig = 1;
                end

                rec_vals_niche = full(expr_matrix(niche_cells, rec_idx));
                rec_vals_other = full(expr_matrix(~niche_cells, rec_idx));

                try
                    p_rec = ranksum(rec_vals_niche, rec_vals_other);
                catch
                    p_rec = 1;
                end

                % Fisher's combined p-value
                p_lig = max(p_lig, 1e-300);
                p_rec = max(p_rec, 1e-300);
                p_val = 1 - chi2cdf(-2 * (log(p_lig) + log(p_rec)), 4);

                observed = sum(strcmp(ligs, lig) & strcmp(recs, rec));
                expected = observed / fold_enrichment;
                niche_enrichment = [niche_enrichment; {lig, rec, niche, observed, expected, fold_enrichment, p_val}];
            end
        end
    end

    fprintf('  %s: tested %d pairs\n', niche, n_test);
end
end
