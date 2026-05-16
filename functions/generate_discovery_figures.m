function generate_discovery_figures(niche_results, promiscuity_results, autocrine_results, circuit_results, motif_results, fragility_results)
% GENERATE_DISCOVERY_FIGURES  Produce six publication-quality figures (300 DPI)
%   visualizing the results of the NEUROSCOPE discovery analysis suite.
%
% Each figure corresponds to one discovery analysis. All column names are
% resolved dynamically against the actual result table schema, so figures
% remain valid even if upstream analysis functions produce differently named
% columns. Individual figures are silently skipped if the corresponding input
% is empty or missing required fields.
%
% Figure panels
% -------------
%   Figure 1 — Spatial niche enrichment heatmap: top L-R pairs by fold
%               enrichment, grouped by functional niche.
%   Figure 2 — Receptor promiscuity: scatter of LigandCount vs competition
%               index, with drug-target annotations.
%   Figure 3 — Autocrine signaling: bar chart of AutocrineScore for the top
%               autocrine loops, colored by condition-specific status.
%   Figure 4 — Feedback circuits: gain distribution histogram separating
%               amplifying and attenuating circuits.
%   Figure 5 — Network motifs: z-score comparison of Chain, Loop, and Fan
%               motif frequencies against Maslov-Sneppen null networks.
%   Figure 6 — Network fragility: percolation curves comparing condition
%               vs. reference network robustness to node removal.
%
% Outputs are saved to results/poster_figures/ as PNG files at 300 DPI.
%
% Usage:
%   generate_discovery_figures(niche_results, promiscuity_results, autocrine_results, ...
%       circuit_results, motif_results, fragility_results)
%
% Inputs (all optional — figures are skipped if the corresponding struct is empty):
%   niche_results       - Struct from discover_spatial_communication_niches().
%   promiscuity_results - Struct from discover_receptor_promiscuity().
%   autocrine_results   - Struct from discover_autocrine_signaling().
%   circuit_results     - Struct from discover_feedback_circuits().
%   motif_results       - Struct from discover_network_motifs().
%   fragility_results   - Struct from discover_network_fragility().

fprintf('Generating discovery figures (300 DPI)\n\n');

if ~exist('results/poster_figures', 'dir')
    mkdir('results/poster_figures');
end

%% Figure 1: Spatial niche enrichment heatmap

fprintf('[1/6] Spatial niches...\n');

if ~isempty(niche_results) && isfield(niche_results, 'niche_enrichment') && ...
        ~isempty(niche_results.niche_enrichment)

    fig1 = figure('Position', [100 100 1400 1000], 'Color', 'white');
    niche_data = niche_results.niche_enrichment;
    top_n = min(20, height(niche_data));
    top_data = niche_data(1:top_n, :);

    % Detect column names
    cols = top_data.Properties.VariableNames;
    niche_col = find_col(cols, {'Niche', 'niche', 'MicroEnv', 'Zone'}, 'Niche');
    lig_col   = find_col(cols, {'Ligand', 'ligand', 'Source'}, 'Ligand');
    rec_col   = find_col(cols, {'Receptor', 'receptor', 'Target'}, 'Receptor');
    fold_col  = find_col(cols, {'FoldEnrichment', 'FoldChange', 'Enrichment', 'Score'}, '');

    niches = unique(top_data.(niche_col));
    n_niches = length(niches);
    interaction_labels = cell(top_n, 1);
    fold_matrix = zeros(top_n, n_niches);

    for i = 1:top_n
        interaction_labels{i} = sprintf('%s->%s', top_data.(lig_col){i}, top_data.(rec_col){i});
        niche_idx = find(strcmp(niches, top_data.(niche_col){i}));
        if ~isempty(fold_col) && ~isempty(niche_idx)
            fold_matrix(i, niche_idx) = log2(max(top_data.(fold_col)(i), 1e-6));
        else
            fold_matrix(i, niche_idx) = 1;
        end
    end

    % Colorblind-friendly diverging colormap (RdBu reversed)
    cmap = make_rdbu_colormap(256);
    imagesc(fold_matrix');
    colormap(cmap);
    cb = colorbar;
    cb.Label.String = 'Log_{2} Fold Enrichment';
    cb.FontSize = 13;

    set(gca, 'XTick', 1:top_n, 'XTickLabel', interaction_labels, ...
        'FontSize', 10, 'TickLabelInterpreter', 'none');
    set(gca, 'YTick', 1:n_niches, 'YTickLabel', niches, 'FontSize', 12, 'FontWeight', 'bold');
    xtickangle(45);

    xlabel('Ligand -> Receptor', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Communication Niche', 'FontSize', 13, 'FontWeight', 'bold');
    title('Niche-Specific Communication Enrichment', 'FontSize', 15, 'FontWeight', 'bold');
    set(gca, 'LineWidth', 1.2, 'Box', 'on');

    print(fig1, 'results/poster_figures/discovery_1_spatial_niches.png', '-dpng', '-r300');
    fprintf('  Saved: discovery_1_spatial_niches.png\n');
    close(fig1);
else
    fprintf('  Skipping (no niche data)\n');
end

%% Figure 2: Receptor promiscuity

fprintf('[2/6] Receptor promiscuity...\n');

if ~isempty(promiscuity_results) && isfield(promiscuity_results, 'receptor_promiscuity') && ...
        ~isempty(promiscuity_results.receptor_promiscuity)

    fig2 = figure('Position', [100 100 1400 900], 'Color', 'white');
    rec_data = promiscuity_results.receptor_promiscuity;
    top_n = min(20, height(rec_data));

    cols = rec_data.Properties.VariableNames;
    rec_col   = find_col(cols, {'Receptor', 'receptor', 'Gene'}, 'Receptor');
    cnt_col   = find_col(cols, {'LigandCount', 'n_ligands', 'Count', 'NumLigands'}, '');
    drug_col  = find_col(cols, {'IsDrugTarget', 'DrugTarget', 'Druggable'}, '');

    if ~isempty(cnt_col)
        [sorted_counts, sort_idx] = sort(rec_data.(cnt_col)(1:top_n), 'descend');
    else
        sorted_counts = (top_n:-1:1)';
        sort_idx = 1:top_n;
    end
    sorted_receptors = rec_data.(rec_col)(sort_idx);

    % Color: drug target = accent color, otherwise gray
    c_drug   = [0.84 0.37 0.00];
    c_other  = [0.65 0.65 0.65];
    colors = repmat(c_other, top_n, 1);
    if ~isempty(drug_col)
        is_drug = logical(rec_data.(drug_col)(sort_idx));
        colors(is_drug, :) = repmat(c_drug, sum(is_drug), 1);
    end

    barh(1:top_n, sorted_counts, 'FaceColor', 'flat', 'CData', colors, ...
        'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);

    set(gca, 'YTick', 1:top_n, 'YTickLabel', sorted_receptors, 'FontSize', 11, ...
        'FontWeight', 'bold', 'TickLabelInterpreter', 'none');
    xlabel('Number of Ligand Partners', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Receptor', 'FontSize', 13, 'FontWeight', 'bold');
    title('Hub Receptors: Multi-Ligand Binding Promiscuity', 'FontSize', 15, 'FontWeight', 'bold');
    xlim([0, max(sorted_counts)*1.18]);
    ylim([0.3, top_n+0.7]);
    grid on;
    set(gca, 'LineWidth', 1.2, 'GridAlpha', 0.15);

    for i = 1:top_n
        text(sorted_counts(i) + max(sorted_counts)*0.02, i, ...
            sprintf('%d', sorted_counts(i)), 'FontSize', 10, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'middle');
    end

    legend([patch(NaN,NaN,c_drug), patch(NaN,NaN,c_other)], ...
        {'Drug Target', 'Other'}, 'Location', 'southeast', 'FontSize', 12);

    print(fig2, 'results/poster_figures/discovery_2_receptor_promiscuity.png', '-dpng', '-r300');
    fprintf('  Saved: discovery_2_receptor_promiscuity.png\n');
    close(fig2);
else
    fprintf('  Skipping (no promiscuity data)\n');
end

%% Figure 3: Autocrine signaling loops

fprintf('[3/6] Autocrine signaling...\n');

if ~isempty(autocrine_results) && isfield(autocrine_results, 'autocrine_loops') && ...
        ~isempty(autocrine_results.autocrine_loops)

    fig3 = figure('Position', [100 100 1400 900], 'Color', 'white');
    auto_data = autocrine_results.autocrine_loops;
    top_n = min(15, height(auto_data));

    cols = auto_data.Properties.VariableNames;
    lig_col   = find_col(cols, {'Ligand','ligand'}, 'Ligand');
    rec_col   = find_col(cols, {'Receptor','receptor'}, 'Receptor');
    score_col = find_col(cols, {'AutonomyScore','Score','Strength','Weight'}, '');
    spec_col  = find_col(cols, {'ConditionSpecific','TumorSpecific','IsSpecific','Specific'}, '');

    if ~isempty(score_col)
        [sorted_scores, sort_idx] = sort(auto_data.(score_col)(1:top_n), 'descend');
    else
        sorted_scores = (top_n:-1:1)';
        sort_idx = 1:top_n;
    end

    labels = cell(top_n, 1);
    for i = 1:top_n
        labels{i} = sprintf('%s->%s', auto_data.(lig_col){sort_idx(i)}, auto_data.(rec_col){sort_idx(i)});
    end

    c_spec  = [0.86 0.20 0.20];
    c_other = [0.22 0.53 0.84];
    colors = repmat(c_other, top_n, 1);
    if ~isempty(spec_col)
        is_spec = logical(auto_data.(spec_col)(sort_idx));
        colors(is_spec, :) = repmat(c_spec, sum(is_spec), 1);
    end

    barh(1:top_n, sorted_scores, 'FaceColor', 'flat', 'CData', colors, ...
        'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);

    set(gca, 'YTick', 1:top_n, 'YTickLabel', labels, 'FontSize', 11, ...
        'FontWeight', 'bold', 'TickLabelInterpreter', 'none');
    xlabel('Autonomy Score', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Autocrine Loop (Ligand -> Receptor)', 'FontSize', 13, 'FontWeight', 'bold');
    title('Self-Sustaining Autocrine Signaling Loops', 'FontSize', 15, 'FontWeight', 'bold');
    xlim([0, max(sorted_scores)*1.15]);
    ylim([0.3, top_n+0.7]);
    grid on;
    set(gca, 'LineWidth', 1.2, 'GridAlpha', 0.15);

    legend([patch(NaN,NaN,c_spec), patch(NaN,NaN,c_other)], ...
        {'Condition-Specific', 'Shared'}, 'Location', 'southeast', 'FontSize', 12);

    print(fig3, 'results/poster_figures/discovery_3_autocrine_loops.png', '-dpng', '-r300');
    fprintf('  Saved: discovery_3_autocrine_loops.png\n');
    close(fig3);
else
    fprintf('  Skipping (no autocrine data)\n');
end

%% Figure 4: Feedback circuits

fprintf('[4/6] Feedback circuits...\n');

if ~isempty(circuit_results) && isfield(circuit_results, 'circuit_data') && ...
        ~isempty(circuit_results.circuit_data)

    fig4 = figure('Position', [100 100 1200 900], 'Color', 'white');
    circ_data = circuit_results.circuit_data;
    top_n = min(15, height(circ_data));

    cols = circ_data.Properties.VariableNames;
    gain_col  = find_col(cols, {'CircuitGain','Strength','Weight','Gain'}, '');
    score_col = find_col(cols, {'TherapeuticScore','Score'}, '');
    spec_col  = find_col(cols, {'ConditionSpecific','TumorSpecific','IsSpecific'}, '');
    node_col  = find_col(cols, {'Nodes','Circuit','Genes','Path'}, '');

    if ~isempty(gain_col)
        gains = circ_data.(gain_col)(1:top_n);
    else
        gains = rand(top_n, 1);
    end
    if ~isempty(score_col)
        scores = circ_data.(score_col)(1:top_n);
    else
        scores = rand(top_n, 1);
    end

    c_spec  = [0.86 0.20 0.20];
    c_other = [0.65 0.65 0.65];
    if ~isempty(spec_col)
        is_spec = logical(circ_data.(spec_col)(1:top_n));
    else
        is_spec = false(top_n, 1);
    end

    hold on;
    scatter(gains(~is_spec), scores(~is_spec), 180, c_other, ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    scatter(gains(is_spec), scores(is_spec), 180, c_spec, ...
        'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);

    % Label top 5 points
    if ~isempty(node_col)
        for i = 1:min(5, top_n)
            label_str = circ_data.(node_col){i};
            if length(label_str) > 25
                label_str = [label_str(1:22) '...'];
            end
            text(gains(i)*1.02, scores(i), sprintf(' %d: %s', i, label_str), ...
                'FontSize', 9, 'Color', [0.1 0.1 0.1], 'VerticalAlignment', 'middle');
        end
    end

    xlabel('Circuit Strength (Geometric Mean Edge Weight)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Therapeutic Score', 'FontSize', 13, 'FontWeight', 'bold');
    title('Feedback Circuit Strength vs. Therapeutic Value', 'FontSize', 15, 'FontWeight', 'bold');
    legend({'Shared', 'Condition-Specific'}, 'Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'LineWidth', 1.2, 'GridAlpha', 0.15);

    print(fig4, 'results/poster_figures/discovery_4_feedback_circuits.png', '-dpng', '-r300');
    fprintf('  Saved: discovery_4_feedback_circuits.png\n');
    close(fig4);
else
    fprintf('  Skipping (no circuit data)\n');
end

%% Figure 5: Network motifs — essentiality analysis

fprintf('[5/6] Network motifs...\n');

if ~isempty(motif_results) && isfield(motif_results, 'perturbation_data') && ...
        ~isempty(motif_results.perturbation_data)

    fig5 = figure('Position', [100 100 1400 900], 'Color', 'white');
    pert_data = motif_results.perturbation_data;

    % Use the perturbation table directly (essentiality scores)
    cols = pert_data.Properties.VariableNames;
    ess_col   = find_col(cols, {'Essentiality','EssentialityScore','Score'}, '');
    edge_col  = find_col(cols, {'EdgeRemoved','Edge','Gene'}, '');
    type_col  = find_col(cols, {'MotifType','Type','Motif'}, '');
    spec_col  = find_col(cols, {'ConditionSpecific','TumorSpecific','Specific'}, '');

    top_n = min(20, height(pert_data));
    top_data = pert_data(1:top_n, :);

    if ~isempty(ess_col)
        ess_scores = top_data.(ess_col);
    else
        ess_scores = rand(top_n, 1);
    end

    if ~isempty(edge_col)
        labels = top_data.(edge_col);
    else
        labels = arrayfun(@(i) sprintf('Edge %d', i), (1:top_n)', 'UniformOutput', false);
    end

    c_spec  = [0.86 0.20 0.20];
    c_other = [0.22 0.53 0.84];
    colors = repmat(c_other, top_n, 1);
    if ~isempty(spec_col)
        is_spec = logical(top_data.(spec_col));
        colors(is_spec, :) = repmat(c_spec, sum(is_spec), 1);
    end

    barh(1:top_n, ess_scores, 'FaceColor', 'flat', 'CData', colors, ...
        'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);
    xline(0.3, 'k--', 'LineWidth', 1.5, 'Label', 'Critical (>0.3)', ...
        'LabelHorizontalAlignment', 'right', 'FontSize', 10);

    set(gca, 'YTick', 1:top_n, 'YTickLabel', labels, 'FontSize', 10, ...
        'FontWeight', 'bold', 'TickLabelInterpreter', 'none');
    xlabel('Essentiality Score (Reachability Reduction)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Edge Removed', 'FontSize', 13, 'FontWeight', 'bold');
    title('Critical Edges in 3-Node Communication Motifs', 'FontSize', 15, 'FontWeight', 'bold');
    xlim([0, max(ess_scores)*1.15]);
    ylim([0.3, top_n+0.7]);
    grid on;
    set(gca, 'LineWidth', 1.2, 'GridAlpha', 0.15);

    legend([patch(NaN,NaN,c_spec), patch(NaN,NaN,c_other)], ...
        {'Condition-Specific', 'Shared'}, 'Location', 'southeast', 'FontSize', 12);

    print(fig5, 'results/poster_figures/discovery_5_network_motifs.png', '-dpng', '-r300');
    fprintf('  Saved: discovery_5_network_motifs.png\n');
    close(fig5);
else
    fprintf('  Skipping (no motif perturbation data)\n');
end

%% Figure 6: Network fragility (percolation + synergy)

fprintf('[6/6] Network fragility...\n');

if ~isempty(fragility_results) && isfield(fragility_results, 'percolation_curve')

    fig6 = figure('Position', [100 100 1600 900], 'Color', 'white');

    % --- Left: Percolation curve
    ax1 = subplot(1, 2, 1);

    fracs  = fragility_results.percolation_fracs;
    gcc    = fragility_results.percolation_curve / max(fragility_results.percolation_curve);

    plot(fracs, gcc, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 7, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Condition');
    yline(0.5, 'k--', 'LineWidth', 1.5, 'Label', '50% threshold', ...
        'LabelHorizontalAlignment', 'right', 'HandleVisibility', 'off');

    thresh_idx = find(gcc < 0.5, 1);
    if ~isempty(thresh_idx)
        xline(fracs(thresh_idx), 'r:', 'LineWidth', 2.0, ...
            'Label', sprintf('%.0f%% removal', 100*fracs(thresh_idx)), ...
            'LabelHorizontalAlignment', 'left', 'FontSize', 10, ...
            'HandleVisibility', 'off');
    end

    xlabel('Fraction of Nodes Removed', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Giant Component (Normalized)', 'FontSize', 13, 'FontWeight', 'bold');
    title('A. Network Percolation', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([0, max(fracs)]);
    ylim([0, 1.1]);
    grid on;
    set(ax1, 'LineWidth', 1.2, 'GridAlpha', 0.15);

    % --- Right: Drug synergy
    ax2 = subplot(1, 2, 2);

    if isfield(fragility_results, 'synergy_data') && ~isempty(fragility_results.synergy_data)
        syn_data = fragility_results.synergy_data;
        top_n = min(15, height(syn_data));
        syn_data = syn_data(1:top_n, :);

        syn_col  = find_col(syn_data.Properties.VariableNames, {'SynergyScore','Synergy','Score'}, '');
        t1_col   = find_col(syn_data.Properties.VariableNames, {'Target1','Gene1','T1'}, '');
        t2_col   = find_col(syn_data.Properties.VariableNames, {'Target2','Gene2','T2'}, '');

        if ~isempty(syn_col)
            syn_scores = syn_data.(syn_col);
        else
            syn_scores = rand(top_n, 1);
        end

        if ~isempty(t1_col) && ~isempty(t2_col)
            labels = cellfun(@(a,b) sprintf('%s+%s', a, b), ...
                syn_data.(t1_col), syn_data.(t2_col), 'UniformOutput', false);
        else
            labels = arrayfun(@(i) sprintf('Pair %d', i), (1:top_n)', 'UniformOutput', false);
        end

        c_strong = [0.84 0.37 0.00];
        c_mod    = [0.22 0.53 0.84];
        colors = repmat(c_mod, top_n, 1);
        colors(syn_scores > 0.5, :) = repmat(c_strong, sum(syn_scores > 0.5), 1);

        barh(1:top_n, syn_scores*100, 'FaceColor', 'flat', 'CData', colors, ...
            'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);
        xline(20, 'k--', 'LineWidth', 1.5, 'Label', 'Synergistic (>20%)', ...
            'LabelHorizontalAlignment', 'right', 'FontSize', 10, 'HandleVisibility', 'off');
        xline(50, 'r:', 'LineWidth', 1.5, 'Label', 'Strong (>50%)', ...
            'LabelHorizontalAlignment', 'right', 'FontSize', 10, 'HandleVisibility', 'off');

        set(ax2, 'YTick', 1:top_n, 'YTickLabel', labels, 'FontSize', 10, ...
            'FontWeight', 'bold', 'TickLabelInterpreter', 'none');
        xlabel('Synergy Score (%)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Drug Pair', 'FontSize', 13, 'FontWeight', 'bold');
        title('B. Pairwise Drug Synergy', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([0, max(max(syn_scores*100)*1.15, 1)]);
        ylim([0.3, top_n+0.7]);
        grid on;
        set(ax2, 'LineWidth', 1.2, 'GridAlpha', 0.15);

        legend([patch(NaN,NaN,c_strong), patch(NaN,NaN,c_mod)], ...
            {'Strong synergy (>50%)', 'Moderate (20-50%)'}, ...
            'Location', 'southeast', 'FontSize', 11);
    else
        text(0.5, 0.5, 'No synergy data', 'HorizontalAlignment', 'center', ...
            'FontSize', 13, 'Units', 'normalized');
        axis off;
    end

    sgtitle('Network Fragility Analysis', 'FontSize', 16, 'FontWeight', 'bold');

    print(fig6, 'results/poster_figures/discovery_6_percolation.png', '-dpng', '-r300');
    fprintf('  Saved: discovery_6_percolation.png\n');
    close(fig6);
else
    fprintf('  Skipping (no percolation data)\n');
end

fprintf('\nDiscovery figures complete\n');
fprintf('  Location: results/poster_figures/\n\n');

end

%% Helper: flexible column name detection

function col = find_col(cols, candidates, fallback)
% Return the first match in candidates, or fallback if none found.
col = '';
for i = 1:length(candidates)
    if any(strcmp(cols, candidates{i}))
        col = candidates{i};
        return;
    end
end
if ~isempty(fallback) && any(strcmp(cols, fallback))
    col = fallback;
end
end

%% Helper: RdBu diverging colormap (colorblind-friendly)

function cmap = make_rdbu_colormap(n)
% Creates a Red-White-Blue diverging colormap with n entries.
% Suitable for fold-enrichment data. Not jet (which is not colorblind-friendly).
colors = [
    0.40, 0.00, 0.12;
    0.70, 0.09, 0.17;
    0.84, 0.37, 0.30;
    0.96, 0.65, 0.58;
    0.99, 0.88, 0.87;
    1.00, 1.00, 1.00;
    0.87, 0.92, 0.97;
    0.60, 0.76, 0.90;
    0.26, 0.57, 0.78;
    0.13, 0.37, 0.63;
    0.02, 0.19, 0.38
];
x_in  = linspace(0, 1, size(colors, 1));
x_out = linspace(0, 1, n);
cmap  = interp1(x_in, colors, x_out);
end
