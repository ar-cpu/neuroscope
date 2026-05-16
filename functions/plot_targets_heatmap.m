function plot_targets_heatmap(targets)
% PLOT_TARGETS_HEATMAP  Render a false-color heatmap of the top therapeutic
%                        target ligand-receptor pairs scored across four
%                        quantitative metrics.
%
% Each row represents a ligand-receptor pair (top 10 by TherapeuticScore).
% Each column represents a normalized metric (column-wise min-max to [0, 1]).
% The hot colormap maps higher values to brighter colors, allowing rapid
% visual identification of pairs with uniformly strong evidence across all
% dimensions.
%
% Renders into the currently active axes. Intended to be called from the
% NEUROSCOPE App Designer GUI (app1.mlapp).
%
% Usage:
%   plot_targets_heatmap(targets)
%
% Inputs:
%   targets - Table from rank_therapeutic_targets() with (at minimum) columns:
%               Ligand, Receptor, Weight, Frequency, MeanStrength, TherapeuticScore.
%             Rows should be pre-sorted by TherapeuticScore descending.

    data_matrix = [targets.Weight, targets.Frequency, targets.MeanStrength, targets.TherapeuticScore];
    data_norm = data_matrix ./ max(data_matrix, [], 1);
    
    imagesc(data_norm);
    colormap('hot');
    colorbar;
    
    set(gca, 'XTick', 1:4);
    set(gca, 'XTickLabel', {'Weight', 'Frequency', 'Mean Strength', 'Score'});
    set(gca, 'YTick', 1:min(10, height(targets)));
    
    y_labels = cell(min(10, height(targets)), 1);
    for label_idx = 1:length(y_labels)
        y_labels{label_idx} = [targets.Ligand{label_idx} '-' targets.Receptor{label_idx}];
    end
    set(gca, 'YTickLabel', y_labels);
    
    xlabel('Metrics');
    ylabel('Ligand-Receptor Pairs');
end