% run_smillie2019.m  —  NEUROSCOPE v5.2 on Smillie et al. 2019 UC colon
% ~18 UC donors, ~12 healthy donors. Epithelial + stromal + immune cells.
% Non-blood tissue — tests generalizability beyond PBMC.
% Run: matlab -batch "cd('c:/arjun-project/neuroscope'); run_smillie2019"

addpath(genpath('functions'));

config = struct();
config.dataset_name    = 'SmillieUC';
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

fprintf('\n=== NEUROSCOPE v5.2 — Smillie et al. 2019 UC Colon ===\n\n');
fprintf('Non-blood tissue: tests generalizability beyond PBMC\n');
fprintf('Cell types: Colonocyte, Goblet, Fibroblast, Endothelial, T, B, Myeloid\n');
fprintf('Expected L-R axes: TGFB1/IL-signaling, WNT, Notch in epithelial-immune crosstalk\n\n');

NEUROSCOPE(config);

fprintf('\nDone. Results in: %s\n', config.output_dir);
