function plot_network_summary(network, title_str)
% PLOT_NETWORK_SUMMARY  Render a bar chart summarizing the topology of a
%                        cell-cell communication network in terms of total
%                        edges, unique ligands, and unique receptors.
%
% Renders into the currently active axes. Intended to be called from the
% NEUROSCOPE App Designer GUI (app1.mlapp) within a UIAxes panel.
%
% Usage:
%   plot_network_summary(network, title_str)
%
% Inputs:
%   network    - Communication network struct with fields:
%                  .weights   - [N×1] edge weight vector
%                  .ligands   - {N×1} ligand gene symbols
%                  .receptors - {N×1} receptor gene symbols
%   title_str  - Character vector used as the axes title.

    if isempty(network.weights)
        bar([0, 0, 0]);
        set(gca, 'XTickLabel', {'Edges', 'Ligands', 'Receptors'});
        title([title_str ' (No data)']);
        ylabel('Count');
        return;
    end
    
    n_edges = length(network.weights);
    n_ligands = length(unique(network.ligands));
    n_receptors = length(unique(network.receptors));
    
    categories_list = {'Edges', 'Ligands', 'Receptors'};
    values = [n_edges, n_ligands, n_receptors];
    
    bar(values);
    set(gca, 'XTickLabel', categories_list);
    title(title_str);
    ylabel('Count');
    
    mean_weight = mean(network.weights);
    text(2, max(values)*0.8, sprintf('Mean Weight: %.3f', mean_weight), ...
         'HorizontalAlignment', 'center', 'FontSize', 10);
end