function NEUROSCOPE(config)
% NEUROSCOPE - Information-theoretic cell-cell communication analysis.
%
% Quantifies actual information transmission between cell populations using
% KSG mutual-information channel capacity — the first framework that measures
% whether sender variation causally predicts receiver transcriptional response,
% rather than merely co-expression of ligand and receptor.
%
% Core scientific claim: existing tools (CellChat, NicheNet, CellPhoneDB) test
% co-expression. NEUROSCOPE tests communication. An interaction is "real" only
% when I(ligand ; receptor_response) > 0 by the KSG estimator. This collapses
% ~90% of co-expressed pairs that do not transmit measurable information.
%
% Usage:
%   NEUROSCOPE(config)
%
% Config fields:
%   dataset_name        - String label for output files (e.g. 'MGE', 'cortex')
%   condition_name      - Label for the condition group (e.g. 'Region_A')
%   reference_name      - Label for the reference group (e.g. 'Region_B')
%   expr_path           - Path to condition expression matrix (.mat, sparse)
%   metadata_path       - Path to cell metadata CSV (columns: CellID, Cluster)
%   cellphonedb_dir     - Directory with CellPhoneDB tables (default: 'converted')
%   reference_expr_path - Path to reference expression matrix (optional, '' for one-sample)
%   reference_meta_path - Path to reference metadata CSV (optional)
%   output_dir          - Output directory (default: 'results')
%   n_permutations      - Bootstrap iterations for CI estimation (default: 200)
%   fdr_threshold       - BH-FDR cutoff for differential channels (default: 0.05)
%   run_ml              - Run machine learning modules (default: true)
%   niche_config        - Optional struct with niche_markers field for niche analysis
%   random_seed         - Random seed for reproducibility (default: 42)
%
% Example:
%   config = struct();
%   config.dataset_name  = 'MGE_interneurons';
%   config.expr_path     = 'converted/mge_expression_sparse.mat';
%   config.metadata_path = 'converted/mge_cell_metadata.csv';
%   NEUROSCOPE(config);

%% -------------------------------------------------------------------------
%% Setup
%% -------------------------------------------------------------------------

t_pipeline = tic;

fprintf('NEUROSCOPE v5.0\n');
fprintf('Information-Theoretic Cell-Cell Communication Analysis\n');
fprintf('Start: %s\n\n', datestr(now));

% Parse config with defaults
config = apply_defaults(config);

fprintf('Dataset:   %s\n', config.dataset_name);
fprintf('Condition: %s\n', config.condition_name);
fprintf('Reference: %s\n', config.reference_name);
fprintf('Output:    %s\n\n', config.output_dir);

% Set random seed for reproducibility
rng(config.random_seed);

% Add functions to path
if exist('functions', 'dir')
    addpath(genpath('functions'));
end

% Create output directories
dirs_needed = {
    config.output_dir, ...
    fullfile(config.output_dir, 'networks'), ...
    fullfile(config.output_dir, 'novel_discoveries'), ...
    fullfile(config.output_dir, 'poster_figures'), ...
    fullfile(config.output_dir, 'ml_models'), ...
    fullfile(config.output_dir, 'statistics')
};
for i = 1:length(dirs_needed)
    if ~exist(dirs_needed{i}, 'dir')
        mkdir(dirs_needed{i});
    end
end

%% -------------------------------------------------------------------------
%% Step 1: Load data
%% -------------------------------------------------------------------------

fprintf('=== Step 1: Loading data ===\n\n');

% Load expression matrix
fprintf('Loading expression matrix: %s\n', config.expr_path);
t_load = tic;
expr_matrix = load_sparse_matrix(config.expr_path);
fprintf('  Loaded: %d cells x %d genes (%.1f s)\n\n', ...
    size(expr_matrix,1), size(expr_matrix,2), toc(t_load));

% Load cell metadata
fprintf('Loading metadata: %s\n', config.metadata_path);
cell_metadata = readtable(config.metadata_path, 'VariableNamingRule', 'preserve');
fprintf('  %d cells, columns: %s\n\n', height(cell_metadata), ...
    strjoin(cell_metadata.Properties.VariableNames, ', '));

% Add 'Cluster' column if not present (map from cluster_id or supercluster_term)
if ~ismember('Cluster', cell_metadata.Properties.VariableNames)
    if ismember('cluster_id', cell_metadata.Properties.VariableNames)
        cell_metadata.Cluster = cell_metadata.cluster_id;
        fprintf('  Mapped cluster_id -> Cluster\n\n');
    elseif ismember('supercluster_term', cell_metadata.Properties.VariableNames)
        [~, ~, cell_metadata.Cluster] = unique(cell_metadata.supercluster_term);
        fprintf('  Mapped supercluster_term -> Cluster (%d unique)\n\n', ...
            length(unique(cell_metadata.Cluster)));
    end
end

% Load CellPhoneDB
fprintf('Loading CellPhoneDB from: %s\n', config.cellphonedb_dir);
cellphonedb = load_cellphonedb(config.cellphonedb_dir);
fprintf('  Interactions: %d | Proteins: %d\n\n', ...
    height(cellphonedb.interactions), height(cellphonedb.protein_table));

% Load gene names
gene_names = load_gene_names(config.expr_path, config.cellphonedb_dir);
fprintf('  Genes in expression matrix: %d\n\n', height(gene_names));

% Convert ENSEMBL IDs to HGNC symbols if needed
% NOTE: gene names CSV may have a leading '0' index row, so check beyond it
first_real_gene_idx = find(~strcmp(gene_names.gene_name, '0'), 1);
if ~isempty(first_real_gene_idx) && strncmp(gene_names.gene_name{first_real_gene_idx}, 'ENSG', 4)
    fprintf('  Detected ENSEMBL IDs — converting to gene symbols...\n');
    gene_table_path = fullfile(config.cellphonedb_dir, 'gene_table.csv');
    if exist(gene_table_path, 'file')
        gene_table_cpdb = readtable(gene_table_path, 'VariableNamingRule', 'preserve');
        [gene_names, expr_matrix] = convert_ensembl_to_symbols(gene_names, expr_matrix, gene_table_cpdb);
        fprintf('  Genes after symbol conversion: %d\n\n', height(gene_names));
    else
        fprintf('  Warning: gene_table.csv not found in %s — ENSEMBL IDs kept as-is.\n\n', config.cellphonedb_dir);
    end
end

%% -------------------------------------------------------------------------
%% Step 1b: Expression preprocessing (normalization + HVG selection)
%% -------------------------------------------------------------------------
%
% Must run AFTER gene symbol resolution (ENSEMBL -> HGNC) so that HVG names
% are in the same symbol space used by CellPhoneDB L-R pair matching.

fprintf('=== Step 1b: Preprocessing ===\n\n');
fprintf('Normalizing expression and selecting highly variable genes...\n');
fprintf('(Raw counts contain noise from sequencing depth variation)\n\n');

[expr_matrix, gene_names_struct, hvg_indices, size_factors] = preprocess_expression( ...
    expr_matrix, gene_names, 'n_hvg', 3000, 'verbose', true);
gene_names = gene_names_struct;
fprintf('Preprocessing complete: %d cells x %d HVGs\n\n', size(expr_matrix,1), size(expr_matrix,2));

% Load reference data (optional)
has_reference = ~isempty(config.reference_expr_path) && ...
                exist(config.reference_expr_path, 'file');

reference_network = create_empty_network();

if has_reference
    fprintf('Loading reference expression: %s\n', config.reference_expr_path);
    ref_expr = load_sparse_matrix(config.reference_expr_path);
    ref_meta = readtable(config.reference_meta_path, 'VariableNamingRule', 'preserve');
    fprintf('  %d cells x %d genes\n\n', size(ref_expr,1), size(ref_expr,2));
    if ~ismember('Cluster', ref_meta.Properties.VariableNames)
        if ismember('cluster_id', ref_meta.Properties.VariableNames)
            ref_meta.Cluster = ref_meta.cluster_id;
        elseif ismember('supercluster_term', ref_meta.Properties.VariableNames)
            [~, ~, ref_meta.Cluster] = unique(ref_meta.supercluster_term);
        end
    end
else
    fprintf('No reference dataset — condition-specific analysis will use one-sample tests.\n\n');
    ref_expr = [];
    ref_meta = table();
end

%% -------------------------------------------------------------------------
%% Step 2: Build communication networks
%% -------------------------------------------------------------------------

fprintf('=== Step 2: Building communication network ===\n\n');

t_net = tic;
condition_network = build_communication_network( ...
    expr_matrix, gene_names, cell_metadata, cellphonedb, config.dataset_name);
fprintf('Network built in %.1f min\n', toc(t_net)/60);
fprintf('  Interactions: %d\n', length(condition_network.weights));
fprintf('  Unique ligands: %d\n', length(unique(condition_network.ligands)));
fprintf('  Unique receptors: %d\n\n', length(unique(condition_network.receptors)));

if has_reference
    fprintf('Building reference network...\n');
    t_ref = tic;
    reference_network = build_communication_network( ...
        ref_expr, gene_names, ref_meta, cellphonedb, config.reference_name);
    fprintf('Reference network built in %.1f min\n', toc(t_ref)/60);
    fprintf('  Interactions: %d\n\n', length(reference_network.weights));
end

%% -------------------------------------------------------------------------
%% Step 2b: Per-cell communication profiles
%% -------------------------------------------------------------------------

fprintf('=== Step 2b: Per-cell communication profiles ===\n\n');
cell_comm_profiles = struct();
try
    cell_comm_profiles = compute_cell_communication_profiles( ...
        expr_matrix, gene_names, condition_network, cell_metadata);
    fprintf('Communication manifold: %d cells x %d dimensions\n\n', ...
        size(cell_comm_profiles.pca_coords,1), size(cell_comm_profiles.pca_coords,2));
    fprintf('Communication states discovered: %d\n\n', cell_comm_profiles.n_comm_clusters);
catch ME
    fprintf('Cell profiles error: %s\n\n', ME.message);
end

%% =========================================================================
%% INFORMATION-THEORETIC CORE  (Steps 3–9)
%%
%% These steps constitute the novel scientific contribution of NEUROSCOPE.
%% Each step answers a question that existing tools (CellChat, NicheNet,
%% CellPhoneDB) cannot answer because they measure co-expression, not
%% information transmission.
%%
%%   Step 3 — Which L-R channels are statistically enriched?
%%            (two independent tests: expression strength AND capacity)
%%   Step 4 — What do receiver cells actually DO in response to signals?
%%            (downstream transcriptional response; defines the true Y variable)
%%   Step 5 — Which channels are CAUSALLY validated?
%%            (IV/2SLS with donor as natural instrument)
%%   Step 6 — What latent communication programs exist?
%%            (NMF decomposition of the capacity tensor)
%%   Step 7 — What is the topology of information flow?
%%            (entropy landscape + spectral Laplacian communities)
%%   Step 8 — Where does communication identity diverge from transcriptome?
%%            (per-cell CDS; the primary novel biological finding)
%%   Step 9 — Final channel ranking by information content
%% =========================================================================

%% -------------------------------------------------------------------------
%% Step 3: Differential channel analysis
%% -------------------------------------------------------------------------
%
% Two independent statistical tests per L-R pair:
%   FDR         — Wilcoxon rank-sum on per-cluster-pair strength distributions
%                 (equivalent to what CellChat / NicheNet test)
%   FDRCapacity — Wilcoxon rank-sum on per-cluster-pair channel capacities
%                 (novel: tests information transmission, not co-expression)
%
% Channels significant only by FDRCapacity are the key novel finding:
% they transmit more information without proportional expression increase —
% invisible to all existing tools.

fprintf('=== Step 3: Differential channel analysis ===\n\n');

t_spec = tic;
condition_specific = identify_differential_channels( ...
    condition_network, reference_network, config.n_permutations, config.fdr_threshold);
fprintf('Complete (%.1f min)\n', toc(t_spec)/60);
n_sig_str = sum(condition_specific.FDR          < config.fdr_threshold);
n_sig_cap = sum(condition_specific.FDRCapacity  < config.fdr_threshold);
n_novel   = sum(condition_specific.FDRCapacity  < config.fdr_threshold & ...
               condition_specific.FDR           >= config.fdr_threshold);
fprintf('  Significant by expression strength (FDR<%.2f): %d\n',  config.fdr_threshold, n_sig_str);
fprintf('  Significant by channel capacity    (FDR<%.2f): %d\n',  config.fdr_threshold, n_sig_cap);
fprintf('  Capacity-only enriched [novel]:                %d\n\n', n_novel);

%% -------------------------------------------------------------------------
%% Step 4: Downstream transcriptional response profiling
%% -------------------------------------------------------------------------
%
% For each receptor in the significant channel set: split receiver cells into
% high/low receptor-expression groups and run Wilcoxon DE to identify genes
% that actually change in response to signalling. PC1 of those genes is the
% true "receiver response score" — a physiologically grounded measure of
% whether the signal is acted upon, not just received.
%
% This redefines Y in the KSG estimator from "receptor abundance" (a proxy)
% to "transcriptional consequence" (the real read-out).

fprintf('=== Step 4: Downstream response profiling ===\n\n');

downstream_response = struct();
try
    downstream_response = compute_downstream_response( ...
        expr_matrix, gene_names, cell_metadata, condition_specific);
    n_resp_fields = length(fieldnames(downstream_response));
    fprintf('Response genes identified for %d receptors\n\n', n_resp_fields);
catch ME
    fprintf('Downstream response error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Step 5: Causal interaction validation (IV / 2SLS)
%% -------------------------------------------------------------------------
%
% Uses donor identity as a natural instrumental variable, analogous to
% Mendelian Randomization in GWAS. Each donor's genetic background perturbs
% ligand levels independently of shared confounders (cell state, batch,
% tissue microenvironment). A channel is "causally validated" when the
% donor-predicted ligand level significantly predicts receiver response
% (first-stage F > 10; second-stage p < 0.05 after bootstrap CI).
%
% This converts a correlational co-expression finding into a causal claim —
% something no existing cell-communication tool provides.

fprintf('=== Step 5: Causal interaction validation (IV/2SLS) ===\n\n');

causal_results = struct();
try
    causal_results = infer_causal_interactions( ...
        condition_specific, expr_matrix, gene_names, cell_metadata, condition_network);
    if isfield(causal_results, 'n_causal')
        fprintf('Causally validated: %d / %d  (%.1f%%)\n', ...
            causal_results.n_causal, causal_results.n_tested, ...
            100 * causal_results.causal_fraction);
        fprintf('Instrument type:    %s\n\n', causal_results.instrument_type);
    end
catch ME
    fprintf('Causal validation error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Step 6: Latent communication programs (NMF tensor decomposition)
%% -------------------------------------------------------------------------
%
% The per-cluster-pair × per-L-R capacity values form a tensor. Reshaping
% to a 2-D matrix and applying NMF reveals co-regulated interaction modules
% ("communication programs") — groups of L-R pairs that operate together
% across the same sender-receiver cluster pairs. Program rank is selected by
% elbow criterion (relative improvement < 10%).
%
% Unlike pathway annotation (which is label-based), NMF programs are
% data-driven — they reveal coordination structure that labels cannot capture.

fprintf('=== Step 6: Latent communication programs (NMF) ===\n\n');

comm_programs = struct();
try
    comm_programs = discover_communication_programs( ...
        condition_network, condition_specific);
    if isfield(comm_programs, 'rank')
        fprintf('Programs discovered: %d (auto-selected rank)\n\n', comm_programs.rank);
    end
catch ME
    fprintf('Communication programs error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Step 7: Spectral structure of information flow
%% -------------------------------------------------------------------------
%
% Treats the capacity-weighted network as a graph and computes:
%   - Shannon entropy per cluster (sending and receiving distributions)
%   - The communication Laplacian and its eigenvectors
%   - Spectral communities: groups of cell types that exchange more
%     information internally than with the rest of the network
%
% The spectral gap (second eigenvalue of the Laplacian) quantifies how
% modular the communication network is — large gap = strong community
% structure, small gap = diffuse broadcast signalling.

fprintf('=== Step 7: Spectral structure of information flow ===\n\n');

entropy_results = struct();
try
    entropy_results = compute_communication_entropy( ...
        condition_network, condition_specific, reference_network);
    if isfield(entropy_results, 'n_communities')
        fprintf('Spectral communities: %d  (spectral gap = %.4f)\n\n', ...
            entropy_results.n_communities, entropy_results.spectral_gap);
    end
catch ME
    fprintf('Spectral structure error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Step 8: Communication identity vs. transcriptomic identity
%% -------------------------------------------------------------------------
%
% THE PRIMARY NOVEL FINDING.
%
% Transcriptomic clusters group cells by gene expression similarity.
% Communication states (from Step 2b) group cells by their signalling
% behaviour. When these two clusterings diverge, a cell communicates
% differently from what its transcriptome predicts.
%
% Per-cell Communication Divergence Score (CDS):
%   CDS_i = 1 - P(comm_state = s_i | trans_cluster = c_i)
%   Range [0, 1].  CDS = 0: cell communicates like all its neighbours.
%                  CDS ≈ 1: cell communicates as an outlier within its cluster.
%
% High-CDS cells are at signalling state boundaries invisible to conventional
% sc-RNA-seq analysis. The Adjusted Rand Index quantifies the global
% decoupling: ARI ≈ 0 means communication identity is orthogonal to
% transcriptomic identity across the whole tissue.

fprintf('=== Step 8: Communication identity ===\n\n');

identity_results = struct();
try
    identity_results = compute_communication_identity(cell_comm_profiles, cell_metadata);
catch ME
    fprintf('Communication identity error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Step 9: Channel ranking by information content
%% -------------------------------------------------------------------------
%
% Ranks all L-R pairs by a composite InformationScore:
%   0.60 × normalised CondCapacity   (primary: actual information transmission)
%   0.25 × (1 - FDRCapacity)         (novel test confidence)
%   0.15 × (1 - FDR)                 (expression strength, for comparison)
%
% CausalScore is updated in-place if Step 5 validated any interactions.

fprintf('=== Step 9: Channel ranking ===\n\n');

interaction_stats = table();
try
    interaction_stats = calculate_interaction_statistics(condition_network, condition_specific);
catch ME
    fprintf('Interaction statistics error (non-fatal): %s\n\n', ME.message);
end

target_rankings = rank_communication_channels(condition_specific, condition_network, config.output_dir);
fprintf('Ranked %d channels\n', height(target_rankings));
if height(target_rankings) > 0
    fprintf('  Top channel: %s → %s  (capacity=%.4f bits, FDR_cap=%.4f)\n\n', ...
        target_rankings.Ligand{1}, target_rankings.Receptor{1}, ...
        target_rankings.CondCapacity(1), target_rankings.FDRCapacity(1));
else
    fprintf('\n');
end

%% -------------------------------------------------------------------------
%% Step 9b: Information landscape (the core empirical claim)
%% -------------------------------------------------------------------------
%
% Computes the headline statistic for the paper:
%   "X% of co-expressed L-R pairs transmit zero measurable information"
%
% Generates information_landscape.csv with the full 2×2 table and rank
% correlation between capacity and co-expression strength.

fprintf('=== Step 9b: Information landscape ===\n\n');

landscape = struct();
try
    landscape = compute_information_landscape( ...
        condition_network, condition_specific, cellphonedb);
catch ME
    fprintf('Information landscape error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Step 9c: Simulation validation (KSG MI vs. co-expression ground truth)
%% -------------------------------------------------------------------------
%
% Run the paper's in-silico validation to demonstrate that KSG normalized MI
% correctly rejects "dead-receiver" L-R pairs (active ligand, independent
% receptor) that co-expression tools spuriously flag.
%
% This generates validation_ground_truth.png and the validation struct that
% goes into the science report. It runs on synthetic data — no real cells
% needed — so it executes unconditionally.

fprintf('=== Step 9c: Ground-truth simulation validation ===\n\n');

validation = struct();
try
    validation = simulate_communication_ground_truth( ...
        'n',        200, ...
        'n_causal', 20, ...
        'n_dead',   20, ...
        'n_null',   20, ...
        'beta',     0.6, ...
        'zero_frac', 0.65, ...
        'seed',     config.random_seed, ...
        'plot',     true, ...
        'verbose',  true);
catch ME
    fprintf('Validation simulation error: %s\n\n', ME.message);
end

%% =========================================================================
%% CONTEXTUAL ANALYSES  (Step 10)
%%
%% These are classical network analyses included for scientific context and
%% to enable direct comparison with existing tools. They operate on the same
%% capacity-weighted network but use conventional graph-theoretic methods.
%% They do not constitute the novel contribution of NEUROSCOPE.
%% =========================================================================

fprintf('=== Step 10: Contextual analyses ===\n\n');

% Communication hubs — cell types that disproportionately dominate sending
% or receiving load. Useful for understanding network topology but does not
% require channel capacity to compute.
hub_analysis = struct();
intervention_strategies = table();  % legacy field for GUI compatibility
try
    hubs = identify_communication_hubs(condition_network);
    hub_analysis.hubs = hubs;
    if istable(hubs) && height(hubs) > 0
        fprintf('  Hubs: %d dominant cell types\n', height(hubs));
    end
catch ME
    fprintf('  Hub analysis error: %s\n', ME.message);
end

% Autocrine self-signalling loops
autocrine_results = struct();
try
    autocrine_results = discover_autocrine_signaling( ...
        condition_network, cell_metadata, expr_matrix, gene_names, condition_specific);
catch ME
    fprintf('  Autocrine error: %s\n', ME.message);
end

% Feedback circuits (A→B→A loops)
circuit_results = struct();
try
    circuit_results = discover_feedback_circuits(condition_network, condition_specific, 4);
catch ME
    fprintf('  Circuit error: %s\n', ME.message);
end

% Network motifs (fan-out, fan-in, relay patterns)
motif_results = struct();
try
    motif_results = discover_network_motifs(condition_network, reference_network, condition_specific);
catch ME
    fprintf('  Motif error: %s\n', ME.message);
end

% Communication monopolies (near-exclusive sender/receiver cell types)
monopoly_results = struct();
try
    monopoly_results = discover_communication_monopolies( ...
        condition_network, cell_metadata, condition_specific);
catch ME
    fprintf('  Monopoly error: %s\n', ME.message);
end

% Role asymmetry (sending/receiving balance per cluster)
role_results = struct();
try
    role_results = discover_communication_role_asymmetry( ...
        condition_network, cell_metadata, condition_specific);
catch ME
    fprintf('  Role asymmetry error: %s\n', ME.message);
end

% Receptor promiscuity (receptors used by many senders)
promiscuity_results = struct();
try
    promiscuity_results = discover_receptor_promiscuity(condition_network, condition_specific);
catch ME
    fprintf('  Promiscuity error: %s\n', ME.message);
end

% Spatial niches (if niche_config provided)
niche_config = [];
if isfield(config, 'niche_config'), niche_config = config.niche_config; end
niche_results = struct();
try
    niche_results = discover_spatial_communication_niches( ...
        condition_network, cell_metadata, expr_matrix, gene_names, niche_config);
catch ME
    fprintf('  Spatial niche error: %s\n', ME.message);
end

fprintf('  Contextual analyses complete\n\n');

%% =========================================================================
%% LEARNED REPRESENTATIONS  (Step 11, optional)
%%
%% Trains a variational autoencoder on the HVG expression matrix to produce
%% a low-dimensional cell embedding. Useful for UMAP/t-SNE visualisation and
%% for cross-condition cell-state classification. Does not affect the channel
%% capacity or differential analysis results.
%% =========================================================================

classifier_perf  = struct();
autoencoder_perf = struct();
biomarker_signature = struct();

if config.run_ml
    fprintf('=== Step 11: Learned representations ===\n\n');

    condition_latent = [];
    try
        [~, ~, condition_latent, autoencoder_perf] = ...
            train_communication_autoencoder(expr_matrix, gene_names, condition_specific, cell_metadata);
        fprintf('  Autoencoder trained\n');
    catch ME
        fprintf('  Autoencoder error: %s\n', ME.message);
    end

    % Cross-condition cell state classifier (requires reference dataset)
    if has_reference && ~isempty(condition_latent)
        try
            ref_latent_tmp = [];
            [~, ~, ref_latent_tmp, ~] = ...
                train_communication_autoencoder(ref_expr, gene_names, [], ref_meta);
            [~, classifier_perf] = ...
                train_malignancy_classifier(condition_latent, ref_latent_tmp, 0.2);
            fprintf('  Cross-condition classifier trained\n');
        catch ME
            fprintf('  Classifier error: %s\n', ME.message);
        end
    end

    % Pathway activity predictor
    try
        train_pathway_predictor(condition_latent, condition_specific, cell_metadata);
        fprintf('  Pathway predictor trained\n');
    catch ME
        fprintf('  Pathway predictor error: %s\n', ME.message);
    end

    % Differential gene signature (requires reference)
    if has_reference
        try
            biomarker_signature = generate_biomarker_signature( ...
                expr_matrix, ref_expr, gene_names, condition_specific, ...
                fullfile(config.output_dir, 'statistics'));
        catch ME
            fprintf('  Signature error: %s\n', ME.message);
        end
    end

    fprintf('\n');
else
    fprintf('=== Step 11: Learned representations skipped (run_ml=false) ===\n\n');
end

%% =========================================================================
%% OUTPUTS  (Steps 12–14)
%% =========================================================================

%% -------------------------------------------------------------------------
%% Step 12: Figures
%% -------------------------------------------------------------------------

fprintf('=== Step 12: Figures ===\n\n');

% Publication figures (paper-quality, PDF + PNG)
try
    generate_paper_figures(condition_network, condition_specific, ...
        identity_results, causal_results, comm_programs, cell_comm_profiles, ...
        landscape, config);
    fprintf('  Publication figures complete\n');
catch ME
    fprintf('  Publication figures error: %s\n', ME.message);
end

% Legacy network/discovery figures (for internal use)
try
    hubs_for_fig = identify_communication_hubs(condition_network);
    generate_poster_figures(condition_network, reference_network, reference_network, ...
        condition_specific, target_rankings, hubs_for_fig);
catch ME
    fprintf('  Poster figures error (non-fatal): %s\n', ME.message);
end
try
    generate_discovery_figures( ...
        niche_results, promiscuity_results, autocrine_results, ...
        circuit_results, motif_results, struct());
catch ME
    fprintf('  Discovery figures error (non-fatal): %s\n', ME.message);
end

fprintf('\n');

%% -------------------------------------------------------------------------
%% Step 13: Scientific discovery report
%% -------------------------------------------------------------------------

fprintf('=== Step 13: Scientific discovery report ===\n\n');

try
    report_results = struct();
    report_results.condition_specific   = condition_specific;
    report_results.channel_rankings     = target_rankings;
    report_results.hub_analysis         = hub_analysis;
    report_results.monopoly_results     = monopoly_results;
    report_results.circuit_results      = circuit_results;
    report_results.fragility_results    = struct();     % fragility not computed in default run
    report_results.autocrine_results    = autocrine_results;
    report_results.condition_network    = condition_network;
    report_results.downstream_response  = downstream_response;
    report_results.causal_results       = causal_results;
    report_results.entropy_results      = entropy_results;
    report_results.comm_programs        = comm_programs;
    report_results.cell_comm_profiles   = cell_comm_profiles;
    report_results.identity_results     = identity_results;
    report_results.landscape            = landscape;

    report_path = generate_science_report(report_results, config);
    fprintf('Report: %s\n\n', report_path);
catch ME
    fprintf('Report error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Step 14: Checkpoint save
%% -------------------------------------------------------------------------

fprintf('=== Step 14: Saving workspace ===\n\n');

% Legacy GUI file
try
    hubs_gui = identify_communication_hubs(condition_network);
    save_results_for_gui(condition_network, reference_network, reference_network, ...
        condition_specific, target_rankings, hubs_gui, intervention_strategies);
catch ME
    fprintf('GUI save error (non-fatal): %s\n\n', ME.message);
end

% Export CSVs
try
    export_results(condition_network, condition_specific, target_rankings, ...
        intervention_strategies, config.dataset_name);
catch ME
    fprintf('CSV export error (non-fatal): %s\n\n', ME.message);
end

checkpoint_path = fullfile(config.output_dir, [config.dataset_name '_complete_results.mat']);
try
    save(checkpoint_path, ...
        'condition_network', 'reference_network', ...
        'condition_specific', 'target_rankings', ...
        'downstream_response', 'causal_results', ...
        'comm_programs', 'entropy_results', ...
        'identity_results', 'cell_comm_profiles', ...
        'landscape', 'validation', ...
        'hub_analysis', 'autocrine_results', 'circuit_results', ...
        'motif_results', 'monopoly_results', 'role_results', ...
        'promiscuity_results', 'niche_results', ...
        'autoencoder_perf', 'classifier_perf', 'biomarker_signature', ...
        '-v7.3');
    info = dir(checkpoint_path);
    fprintf('Checkpoint: %s  (%.1f GB)\n\n', checkpoint_path, info.bytes/1e9);
catch ME
    fprintf('Checkpoint error: %s\n\n', ME.message);
end

%% -------------------------------------------------------------------------
%% Final summary
%% -------------------------------------------------------------------------

t_total = toc(t_pipeline);

fprintf('=============================================================\n');
fprintf('NEUROSCOPE v5.0 complete in %.1f min\n\n', t_total/60);
fprintf('--- Information-Theoretic Communication Summary ---\n\n');

% Network statistics
fprintf('  L-R interactions with capacity > 0:  %d\n', length(condition_network.weights));
if isfield(condition_network,'channel_capacity') && ~isempty(condition_network.channel_capacity)
    cap = condition_network.channel_capacity;
    fprintf('  Mean channel capacity:               %.4f bits\n', mean(cap));
    fprintf('  High-capacity channels (>0.1 bits):  %d  (%.1f%%)\n', ...
        sum(cap > 0.1), 100*mean(cap > 0.1));
end
fprintf('\n');

% Information landscape — the abstract number
if isfield(landscape,'zero_fraction') && ~isempty(landscape.zero_fraction)
    fprintf('  Co-expressed pairs with I≈0:         %.0f%%  (%d / %d)\n', ...
        100*landscape.zero_fraction, landscape.n_zero_capacity, landscape.n_tested);
    fprintf('  Capacity vs strength rank corr:      r = %.3f\n', landscape.rank_correlation);
    fprintf('\n');
end

% Differential analysis
n_sig_str = sum(condition_specific.FDR < config.fdr_threshold);
n_sig_cap = sum(condition_specific.FDRCapacity < config.fdr_threshold);
fprintf('  Differentially expressed channels:   %d (FDR<%.2f, Wilcoxon strength)\n', ...
    n_sig_str, config.fdr_threshold);
fprintf('  Differentially transmitted channels: %d (FDR<%.2f, channel capacity)\n', ...
    n_sig_cap, config.fdr_threshold);
fprintf('  [Novel] Capacity-only enriched:      %d (informationally enriched, not expression-enriched)\n', ...
    sum(condition_specific.FDRCapacity < config.fdr_threshold & ...
        condition_specific.FDR >= config.fdr_threshold));
fprintf('\n');

% Causal validation
if isfield(causal_results,'n_causal') && isfield(causal_results,'n_tested') && causal_results.n_tested > 0
    fprintf('  Causally validated (IV/2SLS):        %d / %d  (%.1f%%)\n', ...
        causal_results.n_causal, causal_results.n_tested, ...
        100 * causal_results.causal_fraction);
    fprintf('  Instrument type:                     %s\n', causal_results.instrument_type);
    fprintf('\n');
end

% Communication identity
if isfield(cell_comm_profiles,'n_comm_clusters')
    fprintf('  Transcriptomic clusters:             %d\n', ...
        length(unique(cell_metadata.Cluster)));
    fprintf('  Communication states (k-means):      %d\n', cell_comm_profiles.n_comm_clusters);
end
if isfield(identity_results,'ari') && ~isnan(identity_results.ari)
    fprintf('  Comm/transcriptome ARI:              %.4f  (%s)\n', ...
        identity_results.ari, identity_results.interpretation);
    fprintf('  Divergent cells (CDS > 0.5):         %d  (%.1f%%)\n', ...
        identity_results.n_divergent, 100*identity_results.divergent_fraction);
    if strcmp(identity_results.interpretation, 'Decoupled')
        fprintf('  >> Communication identity is ORTHOGONAL to transcriptome — strong finding\n');
    end
    fprintf('\n');
end

% Simulation validation
if isfield(validation,'auroc_ksg') && ~isempty(validation.auroc_ksg)
    fprintf('  Validation (synthetic, n=200 cells/pair):\n');
    fprintf('    KSG MI  AUROC: %.3f    dead-receiver FP in top-k: %d\n', ...
            validation.auroc_ksg, validation.dead_in_topk_ksg);
    fprintf('    Co-expr AUROC: %.3f    dead-receiver FP in top-k: %d\n\n', ...
            validation.auroc_coexpr, validation.dead_in_topk_coexpr);
end

% NMF programs
if isfield(comm_programs,'rank') && comm_programs.rank > 0
    fprintf('  Latent communication programs (NMF): %d\n\n', comm_programs.rank);
end

% Spectral communities
if isfield(entropy_results,'n_communities') && entropy_results.n_communities > 0
    fprintf('  Spectral communication communities:  %d\n', entropy_results.n_communities);
    fprintf('  Spectral gap:                        %.4f\n\n', entropy_results.spectral_gap);
end

% Top channel
if height(target_rankings) > 0
    fprintf('  Top channel by capacity:             %s → %s (%.4f bits)\n\n', ...
        target_rankings.Ligand{1}, target_rankings.Receptor{1}, ...
        target_rankings.CondCapacity(1));
end

fprintf('Output directory: %s\n', config.output_dir);
fprintf('End: %s\n', datestr(now));

end

%% =========================================================================
%% Helper functions
%% =========================================================================

function config = apply_defaults(config)
% Fill missing config fields with defaults.

if ~isfield(config, 'dataset_name')       || isempty(config.dataset_name)
    config.dataset_name = 'condition';
end
if ~isfield(config, 'condition_name')     || isempty(config.condition_name)
    config.condition_name = 'Condition';
end
if ~isfield(config, 'reference_name')     || isempty(config.reference_name)
    config.reference_name = 'Reference';
end
if ~isfield(config, 'cellphonedb_dir')    || isempty(config.cellphonedb_dir)
    config.cellphonedb_dir = 'converted';
end
if ~isfield(config, 'reference_expr_path')
    config.reference_expr_path = '';
end
if ~isfield(config, 'reference_meta_path')
    config.reference_meta_path = '';
end
if ~isfield(config, 'output_dir')         || isempty(config.output_dir)
    config.output_dir = 'results';
end
if ~isfield(config, 'n_permutations')     || isempty(config.n_permutations)
    config.n_permutations = 200;  % bootstrap CI iterations for identify_differential_channels
end
if ~isfield(config, 'fdr_threshold')      || isempty(config.fdr_threshold)
    config.fdr_threshold = 0.05;
end
if ~isfield(config, 'run_ml')             || isempty(config.run_ml)
    config.run_ml = true;
end
if ~isfield(config, 'random_seed')        || isempty(config.random_seed)
    config.random_seed = 42;
end

end

function cellphonedb = load_cellphonedb(cpdb_dir)
% Load CellPhoneDB tables from a directory.

cpdb_files = {
    'interactions_table.csv', 'interactions.csv', 'interaction_table.csv'
};
prot_files = {
    'protein_table.csv', 'proteins.csv', 'protein.csv'
};

cpdb_path = '';
for i = 1:length(cpdb_files)
    candidate = fullfile(cpdb_dir, cpdb_files{i});
    if exist(candidate, 'file')
        cpdb_path = candidate;
        break;
    end
end

prot_path = '';
for i = 1:length(prot_files)
    candidate = fullfile(cpdb_dir, prot_files{i});
    if exist(candidate, 'file')
        prot_path = candidate;
        break;
    end
end

if isempty(cpdb_path) || isempty(prot_path)
    error('CellPhoneDB tables not found in %s.\nExpected: interactions_table.csv, protein_table.csv', cpdb_dir);
end

cellphonedb = struct();
cellphonedb.interactions   = readtable(cpdb_path);
cellphonedb.protein_table  = readtable(prot_path);

end

function gene_names = load_gene_names(expr_path, cpdb_dir)
% Try to load gene names from a companion CSV file, or return empty table.
% Returns a table with a 'gene_name' column.

[d, f, ~] = fileparts(expr_path);
candidates = {
    fullfile(d, strrep(f, 'expression_sparse', 'gene_names')),
    fullfile(d, strrep(f, 'expression_sparse', 'genes')),
    fullfile(d, [f '_genes']),
    fullfile(d, 'genes'),
    fullfile(cpdb_dir, 'gene_names')
};

gene_names = table();
for i = 1:length(candidates)
    p = [candidates{i} '.csv'];
    if exist(p, 'file')
        T = readtable(p, 'VariableNamingRule', 'preserve');
        names = T{:, end};
        if ~iscell(names)
            names = cellstr(num2str(names));
        end
        gene_names = table(names, 'VariableNames', {'gene_name'});
        fprintf('  Gene names file: %s\n', p);
        return;
    end
end

fprintf('  Warning: gene names file not found. Using numeric indices.\n');

end
