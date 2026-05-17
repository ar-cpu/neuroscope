function generate_poster_figures(gbm_network, allen_nn_network, allen_mge_network, ...
                                gbm_specific, target_rankings, communication_hubs)
% GENERATE_POSTER_FIGURES  Five publication-quality figures at 300 DPI.
%
% Usage:
%   generate_poster_figures(gbm_network, allen_nn_network, allen_mge_network, ...
%       gbm_specific, target_rankings, communication_hubs)
%
% Inputs:
%   gbm_network       - Condition communication network struct.
%   allen_nn_network  - Reference network 1 struct.
%   allen_mge_network - Reference network 2 struct.
%   gbm_specific      - Table of condition-specific interactions (FDR-corrected).
%   target_rankings   - Table from rank_therapeutic_targets().
%   communication_hubs - Struct from identify_communication_hubs().
%
% Outputs saved to results/poster_figures/ as PNG (300 DPI).

fprintf('      Figure generation (publication-quality with statistics)\n\n');

set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultLineLineWidth', 1.5);
set(groot, 'DefaultFigureRenderer', 'painters');

colors = struct();
colors.gbm      = [0.8500, 0.3250, 0.0980];
colors.healthy  = [0.0000, 0.4470, 0.7410];
colors.specific = [0.9290, 0.6940, 0.1250];
colors.significant = [0.4660, 0.6740, 0.1880];
colors.ns       = [0.5, 0.5, 0.5];
colors.bg       = [0.96, 0.96, 0.96];

fprintf('      [1/5] Network comparison figure...\n');

fig1 = figure('Position', [100, 100, 2000, 700], 'Color', 'white');

subplot(1, 3, 1);
hold on;

y_data = [length(gbm_network.weights); length(allen_nn_network.weights); height(gbm_specific)];
x_pos = 1:3;

bar_colors = [
    0.85, 0.33, 0.10;
    0.00, 0.45, 0.74;
    0.93, 0.69, 0.13
];

b = bar(x_pos, y_data, 'FaceColor', 'flat', 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5, 'BarWidth', 0.65);
b.CData = bar_colors;

for i = 1:3
    y_val = y_data(i);
    text(i, y_val + max(y_data)*0.08, sprintf('%s', format_number(y_val)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 14, ...
        'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.9], ...
        'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1, 'Margin', 3);
end

set(gca, 'XTick', 1:3, 'XTickLabel', {'GBM', 'Healthy', 'Tumor-Specific'}, ...
    'FontSize', 13, 'FontWeight', 'bold');
ylabel('Number of Interactions', 'FontWeight', 'bold', 'FontSize', 15);
title('A. Total Interactions', 'FontSize', 17, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
ylim([0, max(y_data)*1.25]);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15);

subplot(1, 3, 2);
hold on;

gbm_weights = gbm_network.weights;
if ~isempty(allen_nn_network.weights)
    healthy_weights = allen_nn_network.weights;
else
    healthy_weights = zeros(10, 1);
end

[f_gbm, xi_gbm] = ksdensity(log10(gbm_weights + 1));
f_gbm = f_gbm / max(f_gbm) * 0.4;
fill([1-f_gbm, 1+flip(f_gbm)], [xi_gbm, flip(xi_gbm)], bar_colors(1,:), ...
    'FaceAlpha', 0.7, 'EdgeColor', bar_colors(1,:)*0.6, 'LineWidth', 2.5);

[f_h, xi_h] = ksdensity(log10(healthy_weights + 1));
f_h = f_h / max(f_h) * 0.4;
fill([2-f_h, 2+flip(f_h)], [xi_h, flip(xi_h)], bar_colors(2,:), ...
    'FaceAlpha', 0.7, 'EdgeColor', bar_colors(2,:)*0.6, 'LineWidth', 2.5);

med_gbm = median(log10(gbm_weights+1));
med_h   = median(log10(healthy_weights+1));

% White outline then black line for median visibility
plot([0.65, 1.35], [med_gbm, med_gbm], 'Color', 'w', 'LineWidth', 5);
plot([1.65, 2.35], [med_h, med_h],     'Color', 'w', 'LineWidth', 5);
plot([0.65, 1.35], [med_gbm, med_gbm], 'Color', 'k', 'LineWidth', 3);
plot([1.65, 2.35], [med_h, med_h],     'Color', 'k', 'LineWidth', 3);

set(gca, 'XTick', [1, 2], 'XTickLabel', {'GBM', 'Healthy'}, ...
    'FontSize', 13, 'FontWeight', 'bold');
ylabel('Log₁₀(Interaction Weight + 1)', 'FontWeight', 'bold', 'FontSize', 15);
title('B. Weight Distributions', 'FontSize', 17, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
xlim([0.4, 2.6]);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15);

subplot(1, 3, 3);
hold on;

gbm_lig = length(unique(gbm_network.ligands));
gbm_rec = length(unique(gbm_network.receptors));
if ~isempty(allen_nn_network.ligands)
    healthy_lig = length(unique(allen_nn_network.ligands));
    healthy_rec = length(unique(allen_nn_network.receptors));
else
    healthy_lig = 0;
    healthy_rec = 0;
end

x = [1 2];
y = [gbm_lig, gbm_rec; healthy_lig, healthy_rec];
b = bar(x, y, 'grouped', 'BarWidth', 0.75, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);
b(1).FaceColor = bar_colors(1,:);
b(2).FaceColor = bar_colors(2,:);

for i = 1:2
    for j = 1:2
        x_pos = i + (j-1.5)*0.28;
        y_val = y(i,j);
        text(x_pos, y_val + max(y(:))*0.05, sprintf('%d', y_val), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12, ...
            'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.85], ...
            'EdgeColor', [0.4 0.4 0.4], 'Margin', 2);
    end
end

set(gca, 'XTick', 1:2, 'XTickLabel', {'Ligands', 'Receptors'}, ...
    'FontSize', 13, 'FontWeight', 'bold');
ylabel('Unique Genes', 'FontWeight', 'bold', 'FontSize', 15);
title('C. Communication Genes', 'FontSize', 17, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
legend({'GBM', 'Healthy'}, 'Location', 'northwest', 'Box', 'on', 'FontSize', 12, ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.5);
ylim([0, max(y(:))*1.25]);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15);

saveas(fig1, 'results/poster_figures/network_comparison.png');
print(fig1, 'results/poster_figures/network_comparison.png', '-dpng', '-r300');
fprintf('      Saved: network_comparison.png\n');

writetable(table(y_data, 'RowNames', {'GBM', 'Healthy', 'Tumor-Specific'}, ...
    'VariableNames', {'InteractionCount'}), 'results/poster_figures/network_comparison_data.csv');

close(fig1);

fprintf('      [2/5] Volcano plot...\n');

fig2 = figure('Position', [100, 100, 1400, 1000], 'Color', 'white');
hold on;

if any(strcmp(gbm_specific.Properties.VariableNames, 'LogFoldChange'))
    log_fc = gbm_specific.LogFoldChange;
elseif any(strcmp(gbm_specific.Properties.VariableNames, 'Log2FC'))
    log_fc = gbm_specific.Log2FC;
elseif any(strcmp(gbm_specific.Properties.VariableNames, 'FoldChange'))
    log_fc = log2(gbm_specific.FoldChange + 1);
elseif any(strcmp(gbm_specific.Properties.VariableNames, 'CondStrength'))
    log_fc = log2(gbm_specific.CondStrength + 1);
else
    log_fc = zeros(height(gbm_specific), 1);
end

if any(strcmp(gbm_specific.Properties.VariableNames, 'FDR'))
    neg_log_fdr = -log10(gbm_specific.FDR + 1e-300);
else
    neg_log_fdr = ones(height(gbm_specific), 1) * 0.1;
end

significant = neg_log_fdr > -log10(0.05);
high_fc     = abs(log_fc) > 1;
is_target   = high_fc & significant;

color_ns       = [0.7, 0.7, 0.7];
color_target   = [0.85, 0.33, 0.10];
color_sig_only = [0.00, 0.45, 0.74];

scatter(log_fc(~significant), neg_log_fdr(~significant), 60, color_ns, ...
    'filled', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');

sig_low_fc = significant & ~high_fc;
scatter(log_fc(sig_low_fc), neg_log_fdr(sig_low_fc), 80, color_sig_only, ...
    'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');

scatter(log_fc(is_target), neg_log_fdr(is_target), 140, color_target, ...
    'filled', 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.8);

xline(1,  '--', 'LineWidth', 2.5, 'Color', [0.3 0.3 0.3], 'Alpha', 0.5);
xline(-1, '--', 'LineWidth', 2.5, 'Color', [0.3 0.3 0.3], 'Alpha', 0.5);
yline(-log10(0.05), '--', 'LineWidth', 2.5, 'Color', [0.3 0.3 0.3], 'Alpha', 0.5);

top_n = min(6, sum(is_target));
if top_n > 0
    [~, top_idx] = maxk(log_fc(is_target), top_n);
    target_indices = find(is_target);

    label_positions = zeros(top_n, 2);

    for i = 1:top_n
        idx   = target_indices(top_idx(i));
        label = sprintf('%s->%s', gbm_specific.Ligand{idx}, gbm_specific.Receptor{idx});

        x_pos = log_fc(idx);
        y_pos = neg_log_fdr(idx);

        x_offset = 0.15 * sign(x_pos);
        y_offset = max(neg_log_fdr) * 0.03;

        % Push label up if it overlaps a previous one
        if i > 1
            for j = 1:i-1
                dist = sqrt((x_pos - label_positions(j,1))^2 + (y_pos - label_positions(j,2))^2);
                if dist < 0.5
                    y_offset = y_offset + max(neg_log_fdr) * 0.08;
                end
            end
        end

        label_positions(i,:) = [x_pos, y_pos];

        plot([x_pos, x_pos+x_offset], [y_pos, y_pos+y_offset], ...
            'Color', [0.2 0.2 0.2], 'LineWidth', 1.2);

        text(x_pos+x_offset, y_pos+y_offset, label, ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1], ...
            'BackgroundColor', [1 1 1 0.95], 'EdgeColor', [0.3 0.3 0.3], ...
            'LineWidth', 1.2, 'Margin', 3, 'HorizontalAlignment', 'center');
    end
end

xlabel('Log₂ Fold Change (GBM / Healthy)', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('-Log₁₀(FDR)', 'FontWeight', 'bold', 'FontSize', 16);
title('Tumor-Specific Interactions: Volcano Plot', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);

legend({'Non-significant', 'Significant (FDR<0.05)', 'High Impact Targets'}, ...
    'Location', 'northwest', 'Box', 'on', 'FontSize', 12, 'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.5);

grid on; box on;
set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15, 'FontSize', 12);

n_sig = sum(significant);
stats_str = sprintf('Total: %s\nFDR<0.05: %s (%.1f%%)\n|FC|>2: %s (%.1f%%)\nTargets: %s', ...
    format_number(height(gbm_specific)), ...
    format_number(n_sig), 100*n_sig/height(gbm_specific), ...
    format_number(sum(high_fc)), 100*sum(high_fc)/height(gbm_specific), ...
    format_number(sum(is_target)));

annotation('textbox', [0.15, 0.72, 0.20, 0.18], ...
    'String', stats_str, ...
    'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.95], ...
    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 2, 'FitBoxToText', 'on');

saveas(fig2, 'results/poster_figures/volcano_plot.png');
print(fig2, 'results/poster_figures/volcano_plot.png', '-dpng', '-r300');
fprintf('      Saved: volcano_plot.png\n');

volcano_data = table(log_fc, neg_log_fdr, gbm_specific.Ligand, gbm_specific.Receptor, ...
    'VariableNames', {'LogFC', 'NegLogFDR', 'Ligand', 'Receptor'});
writetable(volcano_data, 'results/poster_figures/volcano_plot_data.csv');

close(fig2);

fprintf('      [3/5] Therapeutic targets figure...\n');

fig3 = figure('Position', [100, 100, 1600, 1100], 'Color', 'white');

top_n = min(20, height(target_rankings));
top_targets = target_rankings(1:top_n, :);

labels = cell(top_n, 1);
for i = 1:top_n
    labels{i} = sprintf('%s -> %s', top_targets.Ligand{i}, top_targets.Receptor{i});
end

if ismember('InformationScore', top_targets.Properties.VariableNames)
    poster_score_col = top_targets.InformationScore;
elseif ismember('TherapeuticScore', top_targets.Properties.VariableNames)
    poster_score_col = top_targets.TherapeuticScore;
else
    poster_score_col = zeros(height(top_targets), 1);
end
[sorted_scores, sort_idx] = sort(poster_score_col, 'ascend');
sorted_labels = labels(sort_idx);

colors_grad = zeros(top_n, 3);
for i = 1:top_n
    t = (i-1)/(top_n-1);
    if t < 0.33
        colors_grad(i, :) = [1, t*3, 0];
    elseif t < 0.67
        tt = (t-0.33)/0.34;
        colors_grad(i, :) = [1, 0.5+0.5*tt, 0];
    else
        tt = (t-0.67)/0.33;
        colors_grad(i, :) = [1-tt*0.7, 1, tt*0.2];
    end
end

y_pos = 1:top_n;

barh(y_pos, sorted_scores, 'FaceColor', 'flat', 'CData', colors_grad, ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);

if any(strcmp(top_targets.Properties.VariableNames, 'CI_Lower'))
    sorted_ci_lower = top_targets.CI_Lower(sort_idx);
    sorted_ci_upper = top_targets.CI_Upper(sort_idx);

    hold on;
    for i = 1:top_n
        ci_low  = max(0, sorted_ci_lower(i));
        ci_high = min(1, sorted_ci_upper(i));

        errorbar(sorted_scores(i), y_pos(i), ...
            sorted_scores(i) - ci_low, ...
            ci_high - sorted_scores(i), ...
            'horizontal', 'LineStyle', 'none', 'Color', [0.2 0.2 0.2], 'LineWidth', 2.5, ...
            'CapSize', 10);
    end
end

set(gca, 'YTick', y_pos);
set(gca, 'YTickLabel', sorted_labels, 'FontSize', 11, 'FontWeight', 'normal');
xlabel('Therapeutic Impact Score', 'FontWeight', 'bold', 'FontSize', 16);
ylabel('Ligand -> Receptor Interaction', 'FontWeight', 'bold', 'FontSize', 16);
title('Top 20 Therapeutic Targets (Ranked by Impact)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
xlim([0, 1.05]);
ylim([0.3, top_n+0.7]);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15);

for i = 1:top_n
    x_label = min(sorted_scores(i) + 0.05, 0.98);
    text(x_label, y_pos(i), sprintf('%.3f', sorted_scores(i)), ...
        'FontSize', 10, 'VerticalAlignment', 'middle', 'FontWeight', 'bold', ...
        'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.9], ...
        'EdgeColor', [0.3 0.3 0.3], 'Margin', 2);
end

for i = 1:top_n
    text(-0.03, y_pos(i), sprintf('#%d', top_n-i+1), ...
        'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', ...
        'Color', [0.3 0.3 0.3]);
end

saveas(fig3, 'results/poster_figures/therapeutic_targets.png');
print(fig3, 'results/poster_figures/therapeutic_targets.png', '-dpng', '-r300');
fprintf('      Saved: therapeutic_targets.png\n');

target_export = table(sorted_labels, sorted_scores, 'VariableNames', {'Target', 'Score'});
writetable(target_export, 'results/poster_figures/therapeutic_targets_data.csv');

close(fig3);

fprintf('      [4/5] Communication hubs figure...\n');

fig4 = figure('Position', [100, 100, 2000, 900], 'Color', 'white');

color_ligand   = [0.85, 0.33, 0.10];
color_receptor = [0.00, 0.45, 0.74];

subplot(1, 2, 1);
top_lig  = min(15, height(communication_hubs.ligands));
lig_data = communication_hubs.ligands(1:top_lig, :);

[sorted_counts, sort_idx] = sort(lig_data.ConnectionCount, 'ascend');
sorted_labels = lig_data.Ligand(sort_idx);
y_pos = 1:top_lig;

colors_lig_grad = zeros(top_lig, 3);
for i = 1:top_lig
    t = (i-1)/(top_lig-1);
    colors_lig_grad(i, :) = color_ligand * (0.5 + 0.5*t);
end

barh(y_pos, sorted_counts, 'FaceColor', 'flat', 'CData', colors_lig_grad, ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);

set(gca, 'YTick', y_pos);
set(gca, 'YTickLabel', sorted_labels, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Number of Connections', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Ligand', 'FontWeight', 'bold', 'FontSize', 15);
title('A. Top Ligand Hubs', 'FontSize', 17, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
xlim([0, max(sorted_counts)*1.18]);
ylim([0.3, top_lig+0.7]);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15);

for i = 1:top_lig
    x_label = sorted_counts(i) + max(sorted_counts)*0.03;
    text(x_label, y_pos(i), sprintf('%d', sorted_counts(i)), ...
        'FontSize', 11, 'VerticalAlignment', 'middle', 'FontWeight', 'bold', ...
        'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.9], ...
        'EdgeColor', [0.3 0.3 0.3], 'Margin', 2);
end

for i = 1:top_lig
    text(-max(sorted_counts)*0.03, y_pos(i), sprintf('#%d', top_lig-i+1), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', ...
        'Color', [0.4 0.4 0.4]);
end

subplot(1, 2, 2);
top_rec  = min(15, height(communication_hubs.receptors));
rec_data = communication_hubs.receptors(1:top_rec, :);

[sorted_counts, sort_idx] = sort(rec_data.ConnectionCount, 'ascend');
sorted_labels = rec_data.Receptor(sort_idx);
y_pos = 1:top_rec;

colors_rec_grad = zeros(top_rec, 3);
for i = 1:top_rec
    t = (i-1)/(top_rec-1);
    colors_rec_grad(i, :) = color_receptor * (0.5 + 0.5*t);
end

barh(y_pos, sorted_counts, 'FaceColor', 'flat', 'CData', colors_rec_grad, ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5);

set(gca, 'YTick', y_pos);
set(gca, 'YTickLabel', sorted_labels, 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Number of Connections', 'FontWeight', 'bold', 'FontSize', 15);
ylabel('Receptor', 'FontWeight', 'bold', 'FontSize', 15);
title('B. Top Receptor Hubs', 'FontSize', 17, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
xlim([0, max(sorted_counts)*1.18]);
ylim([0.3, top_rec+0.7]);
grid on; box on;
set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15);

for i = 1:top_rec
    x_label = sorted_counts(i) + max(sorted_counts)*0.03;
    text(x_label, y_pos(i), sprintf('%d', sorted_counts(i)), ...
        'FontSize', 11, 'VerticalAlignment', 'middle', 'FontWeight', 'bold', ...
        'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.9], ...
        'EdgeColor', [0.3 0.3 0.3], 'Margin', 2);
end

for i = 1:top_rec
    text(-max(sorted_counts)*0.03, y_pos(i), sprintf('#%d', top_rec-i+1), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', ...
        'Color', [0.4 0.4 0.4]);
end

saveas(fig4, 'results/poster_figures/communication_hubs.png');
print(fig4, 'results/poster_figures/communication_hubs.png', '-dpng', '-r300');
fprintf('      Saved: communication_hubs.png\n');

writetable(lig_data, 'results/poster_figures/ligand_hubs_data.csv');
writetable(rec_data, 'results/poster_figures/receptor_hubs_data.csv');

close(fig4);

fprintf('      [5/5] Pathway enrichment figure...\n');

fig5 = figure('Position', [100, 100, 1600, 900], 'Color', 'white');

if any(strcmp(gbm_specific.Properties.VariableNames, 'Pathway'))
    pathways = gbm_specific.Pathway;

    all_pathways = {};
    for i = 1:length(pathways)
        if ~isempty(pathways{i}) && ~strcmp(pathways{i}, 'None') && ~strcmp(pathways{i}, 'none')
            splits = strsplit(pathways{i}, ';');
            splits = splits(~cellfun(@isempty, splits));
            all_pathways = [all_pathways, splits];
        end
    end

    if isempty(all_pathways)
        text(0.5, 0.5, 'No pathway annotations found', ...
            'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', ...
            'Color', [0.3 0.3 0.3]);
        axis off;
    else
        [unique_paths, ~, path_idx] = unique(all_pathways);
        path_counts = accumarray(path_idx(:), 1);

        [sorted_counts, sort_idx] = sort(path_counts, 'descend');
        sorted_paths = unique_paths(sort_idx);

        n_paths = min(8, length(sorted_paths));

        pathway_colors = zeros(n_paths, 3);
        for i = 1:n_paths
            pathway_lower = lower(sorted_paths{i});
            if contains(pathway_lower, {'angio', 'vessel', 'vascular'})
                pathway_colors(i, :) = [0.90, 0.20, 0.20];
            elseif contains(pathway_lower, {'immun', 'inflam', 'cytokine'})
                pathway_colors(i, :) = [0.20, 0.50, 0.90];
            elseif contains(pathway_lower, {'growth', 'prolif', 'mitogen'})
                pathway_colors(i, :) = [0.20, 0.75, 0.30];
            elseif contains(pathway_lower, {'stem', 'differ'})
                pathway_colors(i, :) = [0.90, 0.55, 0.10];
            elseif contains(pathway_lower, {'invasion', 'migrat', 'metasta'})
                pathway_colors(i, :) = [0.70, 0.20, 0.80];
            elseif contains(pathway_lower, {'metabol', 'energy'})
                pathway_colors(i, :) = [0.85, 0.75, 0.10];
            elseif contains(pathway_lower, {'survival', 'apopt', 'death'})
                pathway_colors(i, :) = [0.50, 0.20, 0.20];
            else
                pathway_colors(i, :) = [0.50, 0.50, 0.50];
            end
        end

        b = bar(sorted_counts(1:n_paths), 'FaceColor', 'flat', ...
            'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.8);
        b.CData = pathway_colors;

        display_paths = cell(n_paths, 1);
        for i = 1:n_paths
            temp = strrep(sorted_paths{i}, '_', ' ');
            if ~isempty(temp)
                temp(1) = upper(temp(1));
            end
            display_paths{i} = temp;
        end

        set(gca, 'XTick', 1:n_paths, 'XTickLabel', display_paths, 'FontSize', 12, 'FontWeight', 'bold');
        xtickangle(45);
        ylabel('Number of Interactions', 'FontWeight', 'bold', 'FontSize', 16);
        title('Pathway Enrichment in Tumor-Specific Interactions', ...
            'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.1]);
        ylim([0, max(sorted_counts(1:n_paths))*1.25]);
        grid on; box on;
        set(gca, 'LineWidth', 1.5, 'GridAlpha', 0.15);

        for i = 1:n_paths
            y_label = sorted_counts(i) + max(sorted_counts(1:n_paths))*0.05;
            text(i, y_label, sprintf('%d', sorted_counts(i)), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 13, ...
                'Color', [0.1 0.1 0.1], 'BackgroundColor', [1 1 1 0.9], ...
                'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1.2, 'Margin', 2);
        end

        for i = 1:n_paths
            pct   = 100 * sorted_counts(i) / height(gbm_specific);
            y_pct = sorted_counts(i) + max(sorted_counts(1:n_paths))*0.15;
            text(i, y_pct, sprintf('(%.1f%%)', pct), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', ...
                'FontSize', 10, 'Color', [0.3 0.3 0.3]);
        end

        pathway_names_col  = sorted_paths(1:n_paths);
        pathway_counts_col = sorted_counts(1:n_paths);
        pathway_table = table(pathway_names_col(:), pathway_counts_col(:), ...
            'VariableNames', {'Pathway', 'Count'});
        writetable(pathway_table, 'results/poster_figures/pathway_enrichment_data.csv');
    end
else
    text(0.5, 0.5, 'Pathway annotation not available', ...
        'HorizontalAlignment', 'center', 'FontSize', 18, 'FontWeight', 'bold', ...
        'Color', [0.3 0.3 0.3]);
    axis off;
end

saveas(fig5, 'results/poster_figures/pathway_enrichment.png');
print(fig5, 'results/poster_figures/pathway_enrichment.png', '-dpng', '-r300');
fprintf('      Saved: pathway_enrichment.png\n');
close(fig5);

fprintf('\n      Figure generation complete (5 figures, 300 DPI)\n\n');

fprintf('      Generated files:\n');
fprintf('      network_comparison.png + .csv\n');
fprintf('      volcano_plot.png + .csv\n');
fprintf('      therapeutic_targets.png + .csv\n');
fprintf('      communication_hubs.png + .csv (ligands & receptors)\n');
fprintf('      pathway_enrichment.png + .csv\n\n');

end

function str = format_number(num)
% Format integer with comma thousands separators.
if num < 1000
    str = sprintf('%d', num);
else
    str = sprintf('%d', num);
    str = regexprep(str, '(\d)(?=(\d{3})+(?!\d))', '$1,');
end
end
