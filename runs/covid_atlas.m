% COVID Blood Atlas benchmark (Stephenson et al. 2021)
%
% 124 donors across 3 conditions (COVID-19, influenza, healthy).
% 24,284 PBMC cells, 7 cell types.
% Largest donor count in the benchmark suite — gives narrow bootstrap CIs.
%
% Run from the repo root:
%   matlab -batch "run('runs/covid_atlas.m')"
%
% Data prep (first time only):
%   python scripts/prep_covid_atlas.py
%
% Note: this run takes ~21 hours due to n_permutations=200 x 124 donors x 127 pairs.
% Set run_ml=true below to also compute NMF programs and causal inference.

[this_dir, ~, ~] = fileparts(mfilename('fullpath'));
prev_dir = cd(fullfile(this_dir, '..'));
cleanup  = onCleanup(@() cd(prev_dir));
addpath(genpath('functions'));

config = struct();
config.dataset_name    = 'COVIDAtlas';
config.condition_name  = 'PBMC_COVID19';
config.reference_name  = '';
config.expr_path           = 'converted/covid_atlas/covid_expression_sparse.mat';
config.metadata_path       = 'converted/covid_atlas/covid_cell_metadata.csv';
config.cellphonedb_dir     = 'converted';
config.reference_expr_path = '';
config.reference_meta_path = '';
config.output_dir          = 'results/covid_atlas';
config.n_permutations = 200;
config.fdr_threshold  = 0.05;
config.run_ml         = false;
config.random_seed    = 42;

NEUROSCOPE(config);
