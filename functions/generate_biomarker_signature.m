function signature = generate_biomarker_signature(expr_condition, expr_reference, ...
    gene_names, condition_specific, output_dir)
% GENERATE_BIOMARKER_SIGNATURE
%
% Identifies a minimal gene panel that discriminates the condition from the
% reference using LASSO-regularized logistic regression with 5-fold
% cross-validation. Only ligands and receptors from the condition-specific
% interaction table are considered as candidate genes, focusing the panel
% on genes driving the altered communication landscape.
%
% Model selection uses the 1-SE rule (Hastie, Tibshirani & Friedman 2009):
% the most parsimonious lambda whose CV deviance falls within one standard
% error of the minimum. This biases toward fewer genes while accepting a
% small performance penalty — important for clinical feasibility.
%
% A panel with fewer, more discriminative genes is more likely to succeed
% as a clinical assay (e.g., NanoString, qPCR, IHC) than one requiring
% hundreds of genes.
%
% All performance metrics (AUC, sensitivity, specificity, PPV, NPV) are
% estimated on a held-out 20% test set, NOT the training data. Reported
% numbers are therefore unbiased estimates of true out-of-sample performance.
%
% Inputs:
%   expr_condition    - [cells x genes] expression matrix, condition cells
%   expr_reference    - [cells x genes] expression matrix, reference cells
%                       (empty matrix [] if no reference available)
%   gene_names        - Cell array or struct with .gene_name field
%   condition_specific - Table from identify_gbm_specific_interactions
%   output_dir        - Directory for CSV output
%
% Outputs:
%   signature.genes          - Gene panel (cell array of strings)
%   signature.weights        - LASSO log-odds coefficients (condition = positive)
%   signature.n_genes        - Number of genes in the panel
%   signature.auc            - Held-out AUC-ROC
%   signature.sensitivity    - Sensitivity at 0.5 probability threshold
%   signature.specificity    - Specificity at 0.5 probability threshold
%   signature.ppv            - Positive predictive value
%   signature.npv            - Negative predictive value
%   signature.roc_fpr        - ROC false positive rate vector
%   signature.roc_tpr        - ROC true positive rate vector
%   signature.evidence_tier  - Tier A/B/C/D based on held-out AUC
%   signature.table          - Table of selected genes with weights

fprintf('Biomarker Signature Generation\n');
fprintf('  Strategy: LASSO logistic regression (5-fold CV, 1-SE rule)\n');
fprintf('  Candidate pool: communication genes from condition-specific interactions\n\n');

%% Resolve gene names

if istable(gene_names)
    if any(strcmp(gene_names.Properties.VariableNames, 'gene_name'))
        gnames = gene_names.gene_name;
    else
        gnames = gene_names{:, 1};
    end
elseif isstruct(gene_names) && isfield(gene_names, 'gene_name')
    gnames = gene_names.gene_name;
else
    gnames = gene_names;
end
if ~iscell(gnames)
    gnames = cellstr(gnames);
end

%% Candidate communication genes

if istable(condition_specific) && height(condition_specific) > 0 && ...
   any(strcmp(condition_specific.Properties.VariableNames, 'Ligand'))

    ligs = condition_specific.Ligand;
    recs = condition_specific.Receptor;
    if ~iscell(ligs), ligs = cellstr(ligs); end
    if ~iscell(recs), recs = cellstr(recs); end
    comm_genes = unique([ligs; recs]);
else
    % Fallback: most variable genes
    fprintf('  Warning: using top variable genes (no condition_specific table)\n');
    expr_var = var(full(expr_condition), 0, 1);
    [~, top_idx] = sort(expr_var, 'descend');
    n_top = min(300, length(gnames));
    comm_genes = gnames(top_idx(1:n_top));
end

% Map to column indices
gene_idx = [];
gene_panel_names = {};
for i = 1:length(comm_genes)
    idx = find(strcmp(gnames, comm_genes{i}), 1);
    if ~isempty(idx)
        gene_idx(end+1) = idx;                   %#ok<AGROW>
        gene_panel_names{end+1} = comm_genes{i}; %#ok<AGROW>
    end
end

n_candidates = length(gene_panel_names);
fprintf('  Communication gene candidates: %d (of %d unique L-R genes)\n', ...
    n_candidates, length(comm_genes));

if n_candidates < 3
    fprintf('  Too few communication genes mapped — using top variable genes.\n');
    expr_var = var(full(expr_condition), 0, 1);
    [~, top_idx] = sort(expr_var, 'descend');
    gene_idx = top_idx(1:min(200, length(gnames)));
    gene_panel_names = gnames(gene_idx);
end

%% Build labeled dataset

X_cond = full(expr_condition(:, gene_idx));
y_cond = ones(size(X_cond, 1), 1);

if ~isempty(expr_reference)
    X_ref = full(expr_reference(:, gene_idx));
    y_ref = zeros(size(X_ref, 1), 1);
else
    X_ref = zeros(0, length(gene_idx));
    y_ref = zeros(0, 1);
end

X_all = [X_cond; X_ref];
y_all = [y_cond; y_ref];
n_total = size(X_all, 1);
n_cond  = size(X_cond, 1);
n_ref   = size(X_ref, 1);

fprintf('  Dataset: %d condition cells, %d reference cells\n', n_cond, n_ref);

if n_total < 30 || n_ref == 0
    fprintf('  Insufficient data for LASSO (need reference cells).\n');
    signature = make_empty_signature();
    return;
end

%% Log-normalize (matches VAE preprocessing)
X_all = log1p(X_all);

%% Stratified 80/20 train/test split
rng(42);
cv_split = cvpartition(y_all, 'HoldOut', 0.20, 'Stratify', true);
X_train = X_all(training(cv_split), :);
y_train = y_all(training(cv_split));
X_test  = X_all(test(cv_split),     :);
y_test  = y_all(test(cv_split));

fprintf('  Split: %d train / %d test\n\n', length(y_train), length(y_test));

%% LASSO logistic regression with 5-fold cross-validation
%
% lassoglm with 'binomial' / 'logit' fits a L1-penalized logistic
% regression over a path of 50 lambda values. CV deviance is the
% binomial log-likelihood on held-out folds.

fprintf('  Fitting LASSO path (50 lambda values, 5-fold CV)...\n');

try
    [B, FitInfo] = lassoglm(X_train, logical(y_train), 'binomial', ...
        'Link',       'logit', ...
        'Alpha',      1, ...           % pure L1
        'CV',         5, ...
        'Standardize', true, ...
        'NumLambda',  50, ...
        'MaxIter',    1000);
catch ME
    fprintf('  lassoglm error: %s\n', ME.message);
    signature = make_empty_signature();
    return;
end

% Model selection: 1-SE rule
idx_1se = FitInfo.Index1SE;
if isempty(idx_1se)
    idx_1se = FitInfo.IndexMinDeviance;
end

beta_selected  = B(:, idx_1se);
intercept_1se  = FitInfo.Intercept(idx_1se);
nonzero_mask   = (beta_selected ~= 0);
n_selected     = sum(nonzero_mask);

fprintf('  Lambda (1-SE): %.5f | Genes selected: %d\n', FitInfo.Lambda(idx_1se), n_selected);

% If 1-SE is too aggressive, fall back to min-deviance
if n_selected == 0
    fprintf('  Zero genes at 1-SE lambda — using min-deviance lambda instead.\n');
    idx_1se        = FitInfo.IndexMinDeviance;
    beta_selected  = B(:, idx_1se);
    intercept_1se  = FitInfo.Intercept(idx_1se);
    nonzero_mask   = (beta_selected ~= 0);
    n_selected     = sum(nonzero_mask);
    fprintf('  Lambda (min-dev): %.5f | Genes selected: %d\n', FitInfo.Lambda(idx_1se), n_selected);
end

% Sort by absolute weight
weights_nonzero = beta_selected(nonzero_mask);
genes_nonzero   = gene_panel_names(nonzero_mask);
[~, sort_ord]   = sort(abs(weights_nonzero), 'descend');
selected_genes   = genes_nonzero(sort_ord);
selected_weights = weights_nonzero(sort_ord);

fprintf('\n  Gene panel (top 10 by |weight|):\n');
for i = 1:min(10, n_selected)
    arrow = '↑';
    if selected_weights(i) < 0, arrow = '↓'; end
    fprintf('    %2d. %-12s %s  (log-odds coef = %+.4f)\n', ...
        i, selected_genes{i}, arrow, selected_weights(i));
end
if n_selected > 10
    fprintf('    ... and %d more genes\n', n_selected - 10);
end

%% Evaluate on held-out test set

fprintf('\n  Evaluating on held-out 20%% test set...\n');

% Reconstruct beta in original column order for prediction
beta_full = zeros(length(gene_idx), 1);
beta_full(nonzero_mask) = beta_selected(nonzero_mask);

log_odds_test     = X_test * beta_full + intercept_1se;
prob_condition    = 1 ./ (1 + exp(-log_odds_test));
y_pred            = double(prob_condition > 0.5);

tp = sum(y_pred == 1 & y_test == 1);
tn = sum(y_pred == 0 & y_test == 0);
fp = sum(y_pred == 1 & y_test == 0);
fn = sum(y_pred == 0 & y_test == 1);

sensitivity = tp / max(tp + fn, 1);
specificity = tn / max(tn + fp, 1);
ppv         = tp / max(tp + fp, 1);
npv         = tn / max(tn + fn, 1);

[roc_fpr, roc_tpr, ~, auc] = perfcurve(y_test, prob_condition, 1);

fprintf('  AUC-ROC:     %.3f\n', auc);
fprintf('  Sensitivity: %.1f%%\n', sensitivity * 100);
fprintf('  Specificity: %.1f%%\n', specificity * 100);
fprintf('  PPV:         %.1f%%\n', ppv * 100);
fprintf('  NPV:         %.1f%%\n\n', npv * 100);

%% Evidence tier

if auc >= 0.90
    evidence_tier = 'A';
    tier_desc = 'Excellent (AUC ≥ 0.90) — suitable for prospective validation';
elseif auc >= 0.80
    evidence_tier = 'B';
    tier_desc = 'Good (AUC ≥ 0.80) — requires larger validation cohort';
elseif auc >= 0.70
    evidence_tier = 'C';
    tier_desc = 'Moderate (AUC ≥ 0.70) — proof-of-concept only';
else
    evidence_tier = 'D';
    tier_desc = 'Insufficient (AUC < 0.70) — not suitable for diagnostic use';
end

fprintf('  Biomarker signature evidence tier: %s\n  %s\n\n', evidence_tier, tier_desc);

%% Save

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

directions = cell(n_selected, 1);
for i = 1:n_selected
    if selected_weights(i) > 0
        directions{i} = 'Condition-enriched';
    else
        directions{i} = 'Reference-enriched';
    end
end

sig_table = table(selected_genes(:), selected_weights, directions, ...
    'VariableNames', {'Gene', 'LogOdds_Coefficient', 'Direction'});

sig_path = fullfile(output_dir, 'biomarker_signature.csv');
writetable(sig_table, sig_path);
fprintf('  Saved: %s\n\n', sig_path);

%% Package

signature = struct();
signature.genes         = selected_genes;
signature.weights       = selected_weights;
signature.n_genes       = n_selected;
signature.auc           = auc;
signature.sensitivity   = sensitivity;
signature.specificity   = specificity;
signature.ppv           = ppv;
signature.npv           = npv;
signature.roc_fpr       = roc_fpr;
signature.roc_tpr       = roc_tpr;
signature.evidence_tier = evidence_tier;
signature.tier_desc     = tier_desc;
signature.lambda        = FitInfo.Lambda(idx_1se);
signature.table         = sig_table;
signature.n_candidates  = n_candidates;

fprintf('Biomarker signature complete: %d genes selected from %d candidates | AUC=%.3f [Tier %s]\n\n', ...
    n_selected, n_candidates, auc, evidence_tier);

end

function s = make_empty_signature()
s = struct('genes', {{}}, 'weights', [], 'n_genes', 0, ...
    'auc', NaN, 'sensitivity', NaN, 'specificity', NaN, ...
    'ppv', NaN, 'npv', NaN, 'roc_fpr', [], 'roc_tpr', [], ...
    'evidence_tier', 'D', 'tier_desc', 'Insufficient data', ...
    'lambda', NaN, 'table', table(), 'n_candidates', 0);
end
