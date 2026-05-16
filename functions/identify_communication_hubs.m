function hubs = identify_communication_hubs(gbm_network, ~)
% IDENTIFY_COMMUNICATION_HUBS  Rank ligands and receptors by connectivity
%                               degree within the cell-cell communication network.
%
% A communication hub is a gene (ligand or receptor) that participates in an
% unusually large number of distinct signaling interactions. High-degree hubs
% exert broad influence over the communication landscape and are frequently
% prioritized as therapeutic targets because disrupting a hub interaction
% suppresses multiple downstream pathways simultaneously.
%
% Connectivity is measured as the raw degree — the number of edges in the
% network that involve a given gene as either the sending (ligand) or receiving
% (receptor) partner. Both ligand and receptor hub lists are returned separately
% because the two roles have distinct biological interpretations: secreted
% ligand hubs drive paracrine reach, whereas receptor hubs reflect a cell type's
% sensitivity to diverse incoming signals.
%
% Usage:
%   hubs = identify_communication_hubs(network)
%
% Inputs:
%   gbm_network - Communication network struct with fields:
%                   .ligands          - {N×1} cell array of ligand gene symbols
%                   .receptors        - {N×1} cell array of receptor gene symbols
%                 Additional fields (weights, clusters, etc.) are not used here.
%
% Outputs:
%   hubs.ligands   - Table [unique_ligand × 2]: Ligand | ConnectionCount
%                    sorted by ConnectionCount descending.
%   hubs.receptors - Table [unique_receptor × 2]: Receptor | ConnectionCount
%                    sorted by ConnectionCount descending.

    ligand_counts = countcats(categorical(gbm_network.ligands));
    receptor_counts = countcats(categorical(gbm_network.receptors));
    
    unique_ligands = categories(categorical(gbm_network.ligands));
    unique_receptors = categories(categorical(gbm_network.receptors));
    
    ligand_hubs = table(unique_ligands, ligand_counts, ...
                       'VariableNames', {'Ligand', 'ConnectionCount'});
    receptor_hubs = table(unique_receptors, receptor_counts, ...
                         'VariableNames', {'Receptor', 'ConnectionCount'});
    
    ligand_hubs = sortrows(ligand_hubs, 'ConnectionCount', 'descend');
    receptor_hubs = sortrows(receptor_hubs, 'ConnectionCount', 'descend');
    
    hubs.ligands = ligand_hubs;
    hubs.receptors = receptor_hubs;
    
    fprintf('     Identified %d ligand and %d receptor hubs\n', ...
            height(ligand_hubs), height(receptor_hubs));
end