% run_lupus.m  —  NEUROSCOPE v5.2 on Perez et al. 2022 SLE PBMC
% 162 donors (88 SLE + 49 healthy), ~1.2M cells, 9 PBMC cell types.
% Run: matlab -batch "cd('c:/arjun-project/neuroscope'); run_lupus"

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
config.n_permutations = 200;
config.fdr_threshold  = 0.05;
config.run_ml         = false;
config.random_seed    = 42;
config.min_cells_pct       = 0.001;  % 0.1% = 40 cells; needed for large sparse dataset
config.expression_threshold = 0.01;  % lower threshold for GP-smoothed sparse data

fprintf('\n=== NEUROSCOPE v5.2 — Perez et al. 2022 SLE PBMC ===\n\n');
fprintf('162 donors (88 SLE + 49 healthy)\n');
fprintf('Expected: tight bootstrap CIs from large n\n');
fprintf('Expected: stronger 2SLS causal validation at n=162 donors\n\n');

NEUROSCOPE(config);

fprintf('\nDone. Results in: %s\n', config.output_dir);
