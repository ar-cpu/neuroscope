function network = create_empty_network()
% CREATE_EMPTY_NETWORK  Return a zero-interaction network struct that satisfies
%                        the NEUROSCOPE network schema without containing any
%                        actual interaction data.
%
% This function is used as a safe placeholder when a dataset (typically a
% reference/healthy control) has too few cells to support meaningful network
% inference, or when the user chooses to run the pipeline in single-dataset
% mode. Downstream functions that receive this struct will detect the empty
% state via isempty(network.weights) and skip or adapt their analyses.
%
% Usage:
%   network = create_empty_network()
%
% Outputs:
%   network - Struct with the following fields initialized to empty/zero values:
%               .ligands         - {} (empty cell array)
%               .receptors       - {} (empty cell array)
%               .weights         - [] (empty double)
%               .source_clusters - [] (empty double)
%               .target_clusters - [] (empty double)
%               .dataset         - 'empty'
%               .mi_scores       - [] (mutual information scores)
%               .significance    - [] (permutation p-values)
%               .raw_strength    - [] (pre-smoothing interaction strength)

network = struct();
network.ligands = {};
network.receptors = {};
network.weights = [];
network.source_clusters = [];
network.target_clusters = [];
network.dataset = 'empty';

% Optional fields that advanced network builder may include
network.mi_scores = [];
network.significance = [];
network.raw_strength = [];

end