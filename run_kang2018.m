% run_kang2018.m  —  NEUROSCOPE v5.0 on Kang et al. 2018 PBMC
% 8 donors × 2 conditions (ctrl + IFN-beta) = 16 pseudobulk samples
% KSG MI active (>=10 samples, MIN_N=10).
% Run: matlab -batch "cd('c:/arjun-project/neuroscope'); run_kang2018"

addpath(genpath('functions'));

config = struct();

% Identity
config.dataset_name    = 'Kang2018';
config.condition_name  = 'PBMC_IFNbeta';
config.reference_name  = '';

% Paths
config.expr_path           = 'converted/kang2018/kang_expression_sparse.mat';
config.metadata_path       = 'converted/kang2018/kang_cell_metadata.csv';
config.cellphonedb_dir     = 'converted';
config.reference_expr_path = '';
config.reference_meta_path = '';
config.output_dir          = 'results/kang2018';

% Run parameters
config.n_permutations = 200;
config.fdr_threshold  = 0.05;
config.run_ml         = false;
config.random_seed    = 42;

fprintf('\n=== NEUROSCOPE v5.0 — Kang et al. 2018 PBMC ===\n\n');
fprintf('8 donors x 2 conditions = 16 pseudobulk samples\n');
fprintf('KSG MI active (MIN_N=10 threshold)\n');
fprintf('Cell types: CD4T, CD8T, B, NK, CD14+Mono, FCGR3A+Mono, DC, MK\n\n');

NEUROSCOPE(config);

fprintf('\nDone. Results in: %s\n', config.output_dir);
