function plot_hub_barplot(hubs, title_str)
% PLOT_HUB_BARPLOT  Render a horizontal bar chart of communication hub
%                    connectivity counts from identify_communication_hubs().
%
% Each bar represents one gene (ligand or receptor) with bar length equal to
% the number of distinct interaction edges in which that gene participates.
% The function auto-detects whether the input table represents ligands or
% receptors by inspecting title_str for the substring 'Ligand'.
%
% Renders into the currently active axes. Intended to be called from the
% NEUROSCOPE App Designer GUI (app1.mlapp).
%
% Usage:
%   plot_hub_barplot(hubs, title_str)
%
% Inputs:
%   hubs      - Table with columns ConnectionCount and either Ligand or
%               Receptor (output of identify_communication_hubs()).
%               Should be pre-sorted by ConnectionCount descending.
%   title_str - Axes title string. Must contain 'Ligand' to select the
%               Ligand column for y-axis labels; otherwise uses Receptor.

    if height(hubs) == 0
        barh([0]);
        title([title_str ' (No data)']);
        xlabel('Number of Connections');
        return;
    end
    
    barh(hubs.ConnectionCount);
    set(gca, 'YTick', 1:height(hubs));
    
    if contains(title_str, 'Ligand')
        set(gca, 'YTickLabel', hubs.Ligand);
    else
        set(gca, 'YTickLabel', hubs.Receptor);
    end
    
    xlabel('Number of Connections');
    title(title_str);
end