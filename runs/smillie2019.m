% Smillie 2019 ulcerative colitis (UC) colon benchmark (CellxGene)
%
% 18 UC donors, 9,385 colon immune cells, 17 cell types.
% Only non-blood tissue in the benchmark suite — tests whether NEUROSCOPE
% works beyond PBMC. Short runtime (~6 min).
%
% Run from the repo root:
%   matlab -batch "run('runs/smillie2019.m')"
%
% Data prep (first time only):
%   python scripts/prep_smillie2019.py

[this_dir, ~, ~] = fileparts(mfilename('fullpath'));
prev_dir = cd(fullfile(this_dir, '..'));
cleanup  = onCleanup(@() cd(prev_dir));
addpath(genpath('functions'));

config = struct();
config.dataset_name    = 'SmillieIBD';
config.condition_name  = 'UlcerativeColitis';
config.reference_name  = '';
config.expr_path           = 'converted/smillie2019/smillie_expression_sparse.mat';
config.metadata_path       = 'converted/smillie2019/smillie_cell_metadata.csv';
config.cellphonedb_dir     = 'converted';
config.reference_expr_path = '';
config.reference_meta_path = '';
config.output_dir          = 'results/smillie2019';
config.n_permutations = 200;
config.fdr_threshold  = 0.05;
config.run_ml         = false;
config.random_seed    = 42;

NEUROSCOPE(config);
