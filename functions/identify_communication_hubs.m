function hubs = identify_communication_hubs(gbm_network, ~)
% IDENTIFY_COMMUNICATION_HUBS  Rank genes by degree in the communication network.
%
% hubs = identify_communication_hubs(gbm_network)
%
% Outputs:
%   hubs.ligands   - Table: Ligand | ConnectionCount (descending)
%   hubs.receptors - Table: Receptor | ConnectionCount (descending)

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