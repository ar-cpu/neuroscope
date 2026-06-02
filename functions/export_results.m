function export_results(gbm_network, gbm_specific, target_rankings, intervention_strategies, dataset_name)
% EXPORT_RESULTS  Write all NEUROSCOPE outputs to results/ (CSV, JSON, GraphML).
%
% export_results(gbm_network, gbm_specific, target_rankings, intervention_strategies)
% export_results(..., dataset_name)
%
% Inputs:
%   gbm_network             - Network struct (.ligands, .receptors, .weights, ...)
%   gbm_specific            - Condition-specific interactions table (Ligand, Receptor, FDR)
%   target_rankings         - Rankings table from rank_communication_channels()
%   intervention_strategies - Table from propose_interventions(), or []
%   dataset_name            - File prefix (default: 'condition')

if nargin < 5 || isempty(dataset_name)
    dataset_name = 'condition';
end

if ~exist('results', 'dir'), mkdir('results'); end
if ~exist(fullfile('results','networks'), 'dir'), mkdir(fullfile('results','networks')); end
if ~exist(fullfile('results','novel_discoveries'), 'dir'), mkdir(fullfile('results','novel_discoveries')); end

fprintf('Results Export\n\n');

% Phase 1: Export full network
fprintf('Phase 1: Exporting network data...\n');

network_export = table();
network_export.Ligand = gbm_network.ligands;
network_export.Receptor = gbm_network.receptors;
network_export.Weight = gbm_network.weights;
network_export.SourceCluster = gbm_network.source_clusters;
network_export.TargetCluster = gbm_network.target_clusters;

if isfield(gbm_network, 'mi_scores')
    network_export.MutualInformation = gbm_network.mi_scores;
end

if isfield(gbm_network, 'significance')
    network_export.PValue = gbm_network.significance;
end

if isfield(gbm_network, 'raw_strength')
    network_export.RawStrength = gbm_network.raw_strength;
end

writetable(network_export, fullfile('results', 'networks', [dataset_name '_full_network.csv']));
fprintf('  Saved: %s_full_network.csv (%d interactions)\n', dataset_name, height(network_export));

% Phase 2: Export condition-specific interactions
fprintf('\nPhase 2: Exporting condition-specific interactions...\n');

% Add rich annotations
gbm_specific_export = gbm_specific;

% Add descriptions
descriptions = cell(height(gbm_specific), 1);
for i = 1:height(gbm_specific)
    if any(strcmp(gbm_specific.Properties.VariableNames, 'FDR'))
        if gbm_specific.FDR(i) < 0.001
            sig_level = 'highly significant';
        elseif gbm_specific.FDR(i) < 0.01
            sig_level = 'very significant';
        elseif gbm_specific.FDR(i) < 0.05
            sig_level = 'significant';
        else
            sig_level = 'not significant';
        end
    else
        sig_level = 'unknown significance';
    end
    
    descriptions{i} = sprintf('%s signals to %s (%s, FDR=%.4f)', ...
        gbm_specific.Ligand{i}, gbm_specific.Receptor{i}, sig_level, ...
        gbm_specific.FDR(i));
end

gbm_specific_export.Description = descriptions;

writetable(gbm_specific_export, fullfile('results', [dataset_name '_specific_interactions.csv']));
fprintf('  Saved: %s_specific_interactions.csv (%d interactions)\n', dataset_name, height(gbm_specific_export));

% Export significant only
significant_only = gbm_specific_export(gbm_specific_export.FDR < 0.05, :);
writetable(significant_only, fullfile('results', [dataset_name '_specific_significant_only.csv']));
fprintf('  Saved: %s_specific_significant_only.csv (%d interactions)\n', dataset_name, height(significant_only));

% Phase 3: Export rankings
fprintf('\nPhase 3: Exporting therapeutic target rankings...\n');

% Add clinical context
target_rankings_export = target_rankings;

% Add druggability flags (column may not exist in science-mode rankings)
if ismember('DruggabilityScore', target_rankings.Properties.VariableNames)
    druggable_flags = target_rankings.DruggabilityScore > 0.8;
    target_rankings_export.IsDruggable = druggable_flags;
end

% Add priority categories using InformationScore (new) or TherapeuticScore (legacy)
if ismember('InformationScore', target_rankings.Properties.VariableNames)
    score_col = target_rankings.InformationScore;
elseif ismember('TherapeuticScore', target_rankings.Properties.VariableNames)
    score_col = target_rankings.TherapeuticScore;
else
    score_col = zeros(height(target_rankings), 1);
end
priority_categories = cell(height(target_rankings), 1);
for i = 1:height(target_rankings)
    if score_col(i) > 0.8
        priority_categories{i} = 'Top Priority';
    elseif score_col(i) > 0.6
        priority_categories{i} = 'High Priority';
    elseif score_col(i) > 0.4
        priority_categories{i} = 'Medium Priority';
    else
        priority_categories{i} = 'Low Priority';
    end
end
target_rankings_export.PriorityCategory = priority_categories;

writetable(target_rankings_export, fullfile('results', [dataset_name '_target_rankings.csv']));
fprintf('  Saved: %s_target_rankings.csv (%d targets)\n', dataset_name, height(target_rankings_export));

% Export top 20 only
top20 = target_rankings_export(1:min(20, height(target_rankings_export)), :);
writetable(top20, fullfile('results', [dataset_name '_target_rankings_top20.csv']));
fprintf('  Saved: %s_target_rankings_top20.csv (top 20 targets)\n', dataset_name);

% Phase 4: Export intervention strategies
fprintf('\nPhase 4: Exporting intervention strategies...\n');

if ~isempty(intervention_strategies) && istable(intervention_strategies)
    writetable(intervention_strategies, 'results/intervention_strategies.csv');
    fprintf('  Saved: intervention_strategies.csv (%d strategies)\n', height(intervention_strategies));
    
    % Export high priority only
    if any(strcmp(intervention_strategies.Properties.VariableNames, 'Priority'))
        high_priority = intervention_strategies(strcmp(intervention_strategies.Priority, 'High'), :);
        writetable(high_priority, 'results/intervention_strategies_high_priority.csv');
        fprintf('  Saved: intervention_strategies_high_priority.csv (%d strategies)\n', height(high_priority));
    end
end

% Phase 5: Summary statistics
fprintf('\nPhase 5: Generating summary statistics...\n');

summary_stats = struct();

% Network stats
summary_stats.network.total_interactions = length(gbm_network.weights);
summary_stats.network.unique_ligands = length(unique(gbm_network.ligands));
summary_stats.network.unique_receptors = length(unique(gbm_network.receptors));
summary_stats.network.mean_weight = mean(gbm_network.weights);
summary_stats.network.median_weight = median(gbm_network.weights);
summary_stats.network.std_weight = std(gbm_network.weights);

% Condition-specific interaction stats
summary_stats.tumor_specific.total = height(gbm_specific);
summary_stats.tumor_specific.significant_fdr_0_05 = sum(gbm_specific.FDR < 0.05);
summary_stats.tumor_specific.significant_fdr_0_01 = sum(gbm_specific.FDR < 0.01);
summary_stats.tumor_specific.significant_fdr_0_001 = sum(gbm_specific.FDR < 0.001);

if any(strcmp(gbm_specific.Properties.VariableNames, 'FoldChange'))
    summary_stats.tumor_specific.mean_fold_change = mean(gbm_specific.FoldChange);
    summary_stats.tumor_specific.median_fold_change = median(gbm_specific.FoldChange);
    summary_stats.tumor_specific.high_fold_change_count = sum(gbm_specific.FoldChange > 2);
end

% Therapeutic stats
summary_stats.therapeutic.total_targets = height(target_rankings);
% Get score column (InformationScore in science mode, TherapeuticScore legacy)
if ismember('InformationScore', target_rankings.Properties.VariableNames)
    score_vals = target_rankings.InformationScore;
    score_label = 'InformationScore';
elseif ismember('TherapeuticScore', target_rankings.Properties.VariableNames)
    score_vals = target_rankings.TherapeuticScore;
    score_label = 'TherapeuticScore';
else
    score_vals = zeros(height(target_rankings), 1);
    score_label = 'Score';
end
if ismember('DruggabilityScore', target_rankings.Properties.VariableNames)
    drugg_vals = target_rankings.DruggabilityScore;
else
    drugg_vals = zeros(height(target_rankings), 1);
end
summary_stats.therapeutic.top_score = score_vals(1);
summary_stats.therapeutic.mean_score = mean(score_vals);
summary_stats.therapeutic.druggable_count = sum(drugg_vals > 0.8);

% Export as JSON-like text
fid = fopen('results/summary_statistics.txt', 'w');
fprintf(fid, 'NEUROSCOPE ANALYSIS SUMMARY\n\n');

fprintf(fid, 'NETWORK STATISTICS:\n');
fprintf(fid, '  Total interactions: %d\n', summary_stats.network.total_interactions);
fprintf(fid, '  Unique ligands: %d\n', summary_stats.network.unique_ligands);
fprintf(fid, '  Unique receptors: %d\n', summary_stats.network.unique_receptors);
fprintf(fid, '  Mean weight: %.4f\n', summary_stats.network.mean_weight);
fprintf(fid, '  Median weight: %.4f\n', summary_stats.network.median_weight);
fprintf(fid, '  Std deviation: %.4f\n\n', summary_stats.network.std_weight);

fprintf(fid, 'CONDITION-SPECIFIC INTERACTIONS:\n');
fprintf(fid, '  Total: %d\n', summary_stats.tumor_specific.total);
fprintf(fid, '  Significant (FDR<0.05): %d (%.1f%%)\n', summary_stats.tumor_specific.significant_fdr_0_05, ...
    100*summary_stats.tumor_specific.significant_fdr_0_05/summary_stats.tumor_specific.total);
fprintf(fid, '  Significant (FDR<0.01): %d (%.1f%%)\n', summary_stats.tumor_specific.significant_fdr_0_01, ...
    100*summary_stats.tumor_specific.significant_fdr_0_01/summary_stats.tumor_specific.total);
fprintf(fid, '  Significant (FDR<0.001): %d (%.1f%%)\n\n', summary_stats.tumor_specific.significant_fdr_0_001, ...
    100*summary_stats.tumor_specific.significant_fdr_0_001/summary_stats.tumor_specific.total);

if isfield(summary_stats.tumor_specific, 'mean_fold_change')
    fprintf(fid, '  Mean fold-change: %.2f\n', summary_stats.tumor_specific.mean_fold_change);
    fprintf(fid, '  Median fold-change: %.2f\n', summary_stats.tumor_specific.median_fold_change);
    fprintf(fid, '  High fold-change (>2x): %d\n\n', summary_stats.tumor_specific.high_fold_change_count);
end

fprintf(fid, 'THERAPEUTIC TARGETS:\n');
fprintf(fid, '  Total ranked: %d\n', summary_stats.therapeutic.total_targets);
fprintf(fid, '  Top therapeutic score: %.4f\n', summary_stats.therapeutic.top_score);
fprintf(fid, '  Mean therapeutic score: %.4f\n', summary_stats.therapeutic.mean_score);
fprintf(fid, '  Druggable targets: %d (%.1f%%)\n\n', summary_stats.therapeutic.druggable_count, ...
    100*summary_stats.therapeutic.druggable_count/summary_stats.therapeutic.total_targets);

fprintf(fid, 'TOP 10 THERAPEUTIC TARGETS:\n');
for i = 1:min(10, height(target_rankings))
    fprintf(fid, '  %2d. %s -> %s (Score: %.4f)\n', i, ...
        target_rankings.Ligand{i}, target_rankings.Receptor{i}, ...
        score_vals(i));
end

fprintf(fid, '\n');
fprintf(fid, 'Analysis completed: %s\n', datestr(now));
fclose(fid);

fprintf('  Saved: summary_statistics.txt\n');

% Phase 6: Network graph export (GraphML)
fprintf('\nPhase 6: Exporting network graph (GraphML)...\n');

try
    % Create directed graph
    G = digraph(gbm_network.ligands, gbm_network.receptors, gbm_network.weights);
    
    % Add node attributes
    all_nodes = G.Nodes.Name;
    node_types = cell(height(G.Nodes), 1);
    for i = 1:height(G.Nodes)
        if any(strcmp(gbm_network.ligands, all_nodes{i}))
            node_types{i} = 'Ligand';
        else
            node_types{i} = 'Receptor';
        end
    end
    G.Nodes.Type = node_types;
    
    % Export (if function available)
    if exist('writegraph', 'file')
        writegraph(G, fullfile('results', 'networks', [dataset_name '_network.graphml']));
        fprintf('  Saved: %s_network.graphml\n', dataset_name);
    else
        fprintf('  Warning: writegraph not available, skipping GraphML export\n');
    end
catch ME
    fprintf('  Warning: GraphML export failed: %s\n', ME.message);
end

% Phase 7: README file
fprintf('\nPhase 7: Generating README...\n');

readme_file = 'results/README.txt';
fid = fopen(readme_file, 'w');

fprintf(fid, 'NEUROSCOPE ANALYSIS RESULTS\n');
fprintf(fid, 'Cell-Cell Communication Network Analysis\n\n');

fprintf(fid, 'CONTENTS:\n\n');

fprintf(fid, '1. NETWORKS/\n');
fprintf(fid, '   - %s_full_network.csv: Complete communication network\n', dataset_name);
fprintf(fid, '     Columns: Ligand, Receptor, Weight, SourceCluster, TargetCluster\n\n');

fprintf(fid, '2. CONDITION-SPECIFIC INTERACTIONS\n');
fprintf(fid, '   - %s_specific_interactions.csv: All condition-specific interactions\n', dataset_name);
fprintf(fid, '   - %s_specific_significant_only.csv: FDR<0.05 only\n', dataset_name);
fprintf(fid, '     Key columns: Ligand, Receptor, Weight, FDR, FoldChange\n\n');

fprintf(fid, '3. THERAPEUTIC TARGETS\n');
fprintf(fid, '   - %s_target_rankings.csv: All ranked therapeutic targets\n', dataset_name);
fprintf(fid, '   - %s_target_rankings_top20.csv: Top 20 targets only\n', dataset_name);
fprintf(fid, '     Key columns: Ligand, Receptor, TherapeuticScore, DruggabilityScore\n');
fprintf(fid, '     Component scores: ExpressionScore, FrequencyScore, CentralityScore,\n');
fprintf(fid, '                       DisruptionScore, SpecificityScore, PathwayScore\n\n');

fprintf(fid, '4. INTERVENTION STRATEGIES\n');
fprintf(fid, '   - intervention_strategies.csv: All proposed strategies\n');
fprintf(fid, '   - intervention_strategies_high_priority.csv: High priority only\n');
fprintf(fid, '     Columns: StrategyID, Type, Target, Mechanism, Drugs, Priority\n\n');

fprintf(fid, '5. FIGURES/\n');
fprintf(fid, '   - poster_figures/*.png: Publication-quality figures (300 DPI)\n');
fprintf(fid, '     * network_comparison.png\n');
fprintf(fid, '     * volcano_plot.png\n');
fprintf(fid, '     * therapeutic_targets.png\n');
fprintf(fid, '     * communication_hubs.png\n');
fprintf(fid, '     * pathway_enrichment.png\n\n');

fprintf(fid, '6. MACHINE LEARNING MODELS\n');
fprintf(fid, '   - ml_models/autoencoder_metrics.csv: Autoencoder performance\n');
fprintf(fid, '   - ml_models/classifier_performance.csv: Malignancy classifier\n');
fprintf(fid, '   - ml_models/pathway_predictions.csv: Pathway activity predictions\n\n');

fprintf(fid, '7. NOVEL DISCOVERIES\n');
fprintf(fid, '   - novel_discoveries/novel_niche_enrichment.csv: Spatial niches\n');
fprintf(fid, '   - novel_discoveries/novel_receptor_promiscuity.csv: Hub receptors\n');
fprintf(fid, '   - novel_discoveries/novel_autocrine_loops.csv: Self-signaling\n');
fprintf(fid, '   - novel_discoveries/novel_feedback_circuits.csv: Amplification circuits\n\n');

fprintf(fid, 'INTERPRETATION GUIDE:\n\n');

fprintf(fid, 'TherapeuticScore (0-1):\n');
fprintf(fid, '  >0.8: Top priority targets with strong evidence\n');
fprintf(fid, '  0.6-0.8: High priority, clinically relevant\n');
fprintf(fid, '  0.4-0.6: Medium priority, worth investigating\n');
fprintf(fid, '  <0.4: Lower priority\n\n');

fprintf(fid, 'FDR (False Discovery Rate):\n');
fprintf(fid, '  <0.001: Highly significant (***)\n');
fprintf(fid, '  <0.01: Very significant (**)\n');
fprintf(fid, '  <0.05: Significant (*)\n');
fprintf(fid, '  ≥0.05: Not significant (ns)\n\n');

fprintf(fid, 'DruggabilityScore (0-1):\n');
fprintf(fid, '  1.0: FDA-approved drugs available\n');
fprintf(fid, '  0.7: Clinical trial stage\n');
fprintf(fid, '  0.4: Preclinical/experimental\n');
fprintf(fid, '  0.3: Theoretical targetability\n\n');

fprintf(fid, 'CITATION:\n');
fprintf(fid, 'If you use these results, please cite:\n');
fprintf(fid, '  "NEUROSCOPE: Advanced Cell-Cell Communication Network Analysis"\n');
fprintf(fid, '  Computational Biology Pipeline v4.0\n\n');

fprintf(fid, 'CONTACT:\n');
fprintf(fid, '  For questions about methodology or interpretation,\n');
fprintf(fid, '  refer to the NEUROSCOPE documentation.\n\n');

fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Analysis runtime: See main output log\n');

fclose(fid);
fprintf('  Saved: README.txt\n');

% Phase 8: Analysis metadata (JSON)
fprintf('\nPhase 8: Exporting analysis metadata...\n');

metadata = struct();
metadata.analysis_date = datestr(now);
metadata.pipeline_version = 'NEUROSCOPE';
metadata.datasets = struct('condition', dataset_name, 'reference', 'Reference dataset');
metadata.methods = struct(...
    'clustering', 'Louvain/Leiden with SNN graph', ...
    'statistics', 'Permutation + BH-FDR + Bootstrap CI', ...
    'machine_learning', 'Autoencoder + Random Forest', ...
    'network_analysis', 'GP smoothing + MI + Graph diffusion');
metadata.results = summary_stats;

% Export as formatted text (JSON-like)
fid = fopen('results/analysis_metadata.json', 'w');
fprintf(fid, '{\n');
fprintf(fid, '  "analysis_date": "%s",\n', metadata.analysis_date);
fprintf(fid, '  "pipeline_version": "NEUROSCOPE",\n');
fprintf(fid, '  "network_interactions": %d,\n', summary_stats.network.total_interactions);
fprintf(fid, '  "tumor_specific_interactions": %d,\n', summary_stats.tumor_specific.total);
fprintf(fid, '  "significant_interactions_fdr_0.05": %d,\n', summary_stats.tumor_specific.significant_fdr_0_05);
fprintf(fid, '  "therapeutic_targets": %d,\n', summary_stats.therapeutic.total_targets);
fprintf(fid, '  "top_therapeutic_score": %.4f,\n', summary_stats.therapeutic.top_score);
fprintf(fid, '  "druggable_targets": %d\n', summary_stats.therapeutic.druggable_count);
fprintf(fid, '}\n');
fclose(fid);

fprintf('  Saved: analysis_metadata.json\n');

% Final summary
fprintf('\nExport complete\n\n');

fprintf('FILES GENERATED:\n');
fprintf('  NETWORKS:\n');
fprintf('    results/networks/%s_full_network.csv\n', dataset_name);
fprintf('  CONDITION-SPECIFIC:\n');
fprintf('    results/%s_specific_interactions.csv\n', dataset_name);
fprintf('    results/%s_specific_significant_only.csv\n', dataset_name);
fprintf('  THERAPEUTIC:\n');
fprintf('    results/%s_target_rankings.csv\n', dataset_name);
fprintf('    results/%s_target_rankings_top20.csv\n', dataset_name);
fprintf('  STRATEGIES:\n');
fprintf('    results/intervention_strategies.csv\n');
fprintf('  DOCUMENTATION:\n');
fprintf('    results/README.txt\n');
fprintf('    results/summary_statistics.txt\n');
fprintf('    results/analysis_metadata.json\n\n');

fprintf('KEY FINDINGS:\n');
fprintf('  %d total interactions in network\n', summary_stats.network.total_interactions);
fprintf('  %d condition-specific interactions (FDR<0.05: %d)\n', ...
    summary_stats.tumor_specific.total, summary_stats.tumor_specific.significant_fdr_0_05);
fprintf('  %d therapeutic targets ranked\n', summary_stats.therapeutic.total_targets);
fprintf('  %d druggable targets identified\n\n', summary_stats.therapeutic.druggable_count);

fprintf('TOP THERAPEUTIC TARGET:\n');
fprintf('  %s -> %s (Score: %.4f)\n\n', ...
    target_rankings.Ligand{1}, target_rankings.Receptor{1}, score_vals(1));

fprintf('All results exported.\n\n');

end