function save_results_for_gui(gbm_network, allen_nn_network, allen_mge_network, ...
                             gbm_specific, target_rankings, communication_hubs, ...
                             intervention_strategies)
% SAVE_RESULTS_FOR_GUI  Persist all primary analysis outputs to a single .mat
%                        file for consumption by the NEUROSCOPE MATLAB App Designer
%                        interface (app1.mlapp).
%
% The GUI loads results/gbm_analysis_results.mat on startup to populate its
% interactive panels (network graph, volcano plot, target table, hub bar chart).
% Saving as a single .mat file rather than multiple CSVs allows the app to
% restore the full MATLAB workspace state — including sparse matrices and struct
% hierarchies — without parsing overhead.
%
% Usage:
%   save_results_for_gui(gbm_network, allen_nn_network, allen_mge_network, ...
%       gbm_specific, target_rankings, communication_hubs, intervention_strategies)
%
% Inputs:
%   gbm_network             - Condition communication network struct.
%   allen_nn_network        - Reference network 1 struct (e.g., Allen Brain Atlas neurons).
%   allen_mge_network       - Reference network 2 struct (e.g., Allen MGE interneurons).
%   gbm_specific            - Table of condition-specific interactions.
%   target_rankings         - Table of ranked therapeutic targets.
%   communication_hubs      - Struct with .ligands and .receptors hub tables.
%   intervention_strategies - Table of proposed intervention strategies.
%
% Outputs:
%   Writes results/gbm_analysis_results.mat containing all seven input
%   variables under their original names.

    save('results/gbm_analysis_results.mat', 'gbm_network', 'allen_nn_network', ...
         'allen_mge_network', 'gbm_specific', 'target_rankings', ...
         'communication_hubs', 'intervention_strategies');
    
    fprintf('     Saved results for GUI\n');
end