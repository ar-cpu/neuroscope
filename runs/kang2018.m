% Kang 2018 IFN-beta PBMC benchmark (GSE96583)
%
% 8 donors, 2 conditions (control + IFN-beta stimulated), 8 PBMC cell types.
% This is the core validation dataset — IFN-beta stimulation is a known
% intervention, so the ground truth is well-defined.
%
% Run from the repo root:
%   matlab -batch "run('runs/kang2018.m')"
%
% Data prep (first time only):
%   python scripts/prep_kang2018.py

% Navigate to repo root so relative paths work
[this_dir, ~, ~] = fileparts(mfilename('fullpath'));
prev_dir = cd(fullfile(this_dir, '..'));
cleanup  = onCleanup(@() cd(prev_dir));
addpath(genpath('functions'));

config = struct();
config.dataset_name    = 'Kang2018';
config.condition_name  = 'PBMC_IFNbeta';
config.reference_name  = '';
config.expr_path           = 'converted/kang2018/kang_expression_sparse.mat';
config.metadata_path       = 'converted/kang2018/kang_cell_metadata.csv';
config.cellphonedb_dir     = 'converted';
config.reference_expr_path = '';
config.reference_meta_path = '';
config.output_dir          = 'results/kang2018';
config.n_permutations = 200;
config.fdr_threshold  = 0.05;
config.run_ml         = false;
config.random_seed    = 42;

NEUROSCOPE(config);
