% Lupus SLE PBMC benchmark (Perez et al. 2022, GSE174188)
%
% 149 donors with systemic lupus erythematosus, 40,395 cells, 11 PBMC subtypes.
% This dataset is large and sparse — balanced sampling (30 cells/cluster/donor)
% keeps memory manageable. The two non-default thresholds below are needed
% because of that sparsity.
%
% Run from the repo root:
%   matlab -batch "run('runs/lupus.m')"
%
% Data prep (first time only):
%   python scripts/prep_lupus.py

[this_dir, ~, ~] = fileparts(mfilename('fullpath'));
prev_dir = cd(fullfile(this_dir, '..'));
cleanup  = onCleanup(@() cd(prev_dir));
addpath(genpath('functions'));

config = struct();
config.dataset_name    = 'Lupus';
config.condition_name  = 'SLE';
config.reference_name  = '';
config.expr_path           = 'converted/lupus/lupus_expression_sparse.mat';
config.metadata_path       = 'converted/lupus/lupus_cell_metadata.csv';
config.cellphonedb_dir     = 'converted';
config.reference_expr_path = '';
config.reference_meta_path = '';
config.output_dir          = 'results/lupus';
config.n_permutations      = 200;
config.fdr_threshold       = 0.05;
config.run_ml              = true;
config.random_seed         = 42;
config.min_cells_pct       = 0.001;  % relaxed for sparse balanced-sampled data
config.expression_threshold = 0.01; % GP smoothing dilutes means in sparse data

NEUROSCOPE(config);
