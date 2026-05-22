% run_smillie2019.m  —  NEUROSCOPE v5.2 on colon immune (Crohn's disease)
% CellxGene colon immune: 152K cells, 34 donors (18 Crohn + 16 healthy)
% Prepared subset: 9,385 Crohn disease cells, 18 donors, 17 cell types
% Non-blood tissue — tests generalizability beyond PBMC.
% Run: matlab -batch "cd('c:/arjun-project/neuroscope'); run_smillie2019"

addpath(genpath('functions'));

config = struct();
config.dataset_name    = 'SmillieIBD';
config.condition_name  = 'CrohnDisease';
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

fprintf('\n=== NEUROSCOPE v5.2 — Colon Immune (Crohn''s Disease) ===\n\n');
fprintf('Non-blood tissue: tests generalizability beyond PBMC\n');
fprintf('18 Crohn disease donors, 17 immune cell types\n');
fprintf('Expected L-R axes: IL-6, TNF, TGF-b1, IL-1b in intestinal inflammation\n\n');

NEUROSCOPE(config);

fprintf('\nDone. Results in: %s\n', config.output_dir);
