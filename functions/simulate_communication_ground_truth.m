function results = simulate_communication_ground_truth(varargin)
% SIMULATE_COMMUNICATION_GROUND_TRUTH  Validate the KSG normalized-MI estimator
%   against known ground truth on synthetic data mimicking scRNA-seq distributions.
%
% PURPOSE
% -------
% This script generates the validation data for the NEUROSCOPE paper. It
% demonstrates the central claim:
%
%   Co-expression-based scores (geometric mean of cluster means) cannot
%   distinguish truly communicating L-R pairs from pairs that are merely
%   co-expressed. The KSG normalized MI can.
%
% THREE SIMULATION TYPES
% ----------------------
%
%   Type 1: CAUSAL
%     Y = β·X + ε,  X ~ zero-inflated lognormal
%     Ground truth: I(X;Y) > 0 (Y is causally downstream of X)
%     Both co-expression and KSG MI should flag this pair.
%
%   Type 2: DEAD RECEIVER
%     X ~ zero-inflated lognormal (variable ligand)
%     Y ~ zero-inflated lognormal (constant, independent of X)
%     Ground truth: I(X;Y) = 0 (receptor varies but ignores ligand)
%     Co-expression will flag this pair (both genes are "on").
%     KSG MI should correctly identify I ≈ 0.
%     ** This is the key failure mode of existing tools. **
%
%   Type 3: NULL
%     X ~ zero-inflated lognormal
%     Y = constant + tiny noise
%     Ground truth: I(X;Y) = 0
%     Neither tool should flag this (Y has no variance to explain).
%
% EVALUATION METRICS
% ------------------
%   - Precision@k  (fraction of top-k flagged pairs that are truly causal)
%   - Recall@k     (fraction of causal pairs in top-k)
%   - AUROC        (area under receiver-operating-characteristic curve)
%   - Calibration  (estimated MI vs. theoretical MI under Gaussian approximation)
%
% USAGE
%   results = simulate_communication_ground_truth()
%   results = simulate_communication_ground_truth('n', 150, 'n_causal', 20, 'seed', 42)
%
% INPUTS (name-value)
%   'n'           - cells per cluster (default 200)
%   'n_causal'    - number of causal pairs (default 20)
%   'n_dead'      - number of dead-receiver pairs (default 20)
%   'n_null'      - number of null pairs (default 20)
%   'beta'        - causal coupling strength (default 0.6)
%   'noise_sd'    - noise std in causal model (default 0.8)
%   'zero_frac'   - fraction of zero counts (zero-inflation, default 0.6)
%   'seed'        - RNG seed (default 1)
%   'plot'        - show summary figures (default true)
%   'verbose'     - print detailed output (default true)

% ── Parse arguments ──────────────────────────────────────────────────────────
p = inputParser();
addParameter(p, 'n',        200);
addParameter(p, 'n_causal', 20);
addParameter(p, 'n_dead',   20);
addParameter(p, 'n_null',   20);
addParameter(p, 'beta',     0.6);
addParameter(p, 'noise_sd', 0.8);
addParameter(p, 'zero_frac', 0.6);
addParameter(p, 'seed',     1);
addParameter(p, 'plot',     true);
addParameter(p, 'verbose',  true);
parse(p, varargin{:});
opt = p.Results;

rng(opt.seed);
n          = opt.n;
n_causal   = opt.n_causal;
n_dead     = opt.n_dead;
n_null     = opt.n_null;
beta       = opt.beta;
noise_sd   = opt.noise_sd;
zf         = opt.zero_frac;
n_total    = n_causal + n_dead + n_null;

if opt.verbose
    fprintf('\n=== NEUROSCOPE ground-truth simulation ===\n\n');
    fprintf('  n cells per pair  : %d\n', n);
    fprintf('  causal pairs      : %d\n', n_causal);
    fprintf('  dead-receiver     : %d\n', n_dead);
    fprintf('  null pairs        : %d\n', n_null);
    fprintf('  coupling β        : %.2f\n', beta);
    fprintf('  noise σ           : %.2f\n', noise_sd);
    fprintf('  zero-inflation    : %.0f%%\n\n', zf * 100);
end

% ── Helper: zero-inflated lognormal sample ────────────────────────────────────
% Mimics scRNA-seq UMI counts after log1p normalisation.
zilnorm = @(sz) zero_inflate(lognrnd(0.5, 0.8, sz, 1), zf);

% ── Generate synthetic L-R pairs ─────────────────────────────────────────────
mi_ksg      = zeros(n_total, 1);    % our estimator
coexpr      = zeros(n_total, 1);    % competitor: geom mean of means
label       = zeros(n_total, 1);    % 1=causal, 0=not causal
pair_type   = cell(n_total, 1);

idx = 0;

% --- Causal pairs ---
% n observations simulate n donors. X(i) = ligand mean for donor i,
% Y(i) = receptor mean for donor i. This is the correct pairing contract.
for i = 1:n_causal
    idx = idx + 1;
    X = zilnorm(n);
    % Y is linearly driven by X — same donor, downstream effect
    Y_raw = beta .* X + noise_sd .* randn(n, 1);
    Y_raw = max(Y_raw, 0);           % non-negative expression
    Y = zero_inflate(Y_raw, zf * 0.5); % receiver has lower drop-out

    [mi_norm, ~, ~, ~] = estimate_channel_capacity(X, Y, 'n_perm', 0);
    mi_ksg(idx)   = mi_norm;
    coexpr(idx)   = geom_mean_score(X, Y);
    label(idx)    = 1;
    pair_type{idx} = 'causal';
end

% --- Dead-receiver pairs ---
% Key scenario: X varies (active ligand), Y is independent (receptor expressed
% but not listening). Co-expression tools flag this; KSG should not.
for i = 1:n_dead
    idx = idx + 1;
    X = zilnorm(n);
    Y = zilnorm(n);   % completely independent of X

    [mi_norm, ~, ~, ~] = estimate_channel_capacity(X, Y, 'n_perm', 0);
    mi_ksg(idx)   = mi_norm;
    coexpr(idx)   = geom_mean_score(X, Y);
    label(idx)    = 0;
    pair_type{idx} = 'dead_receiver';
end

% --- Null pairs ---
% Y is nearly constant — no variance to transmit into.
for i = 1:n_null
    idx = idx + 1;
    X = zilnorm(n);
    Y = 0.01 * abs(randn(n, 1));    % essentially constant

    [mi_norm, ~, ~, ~] = estimate_channel_capacity(X, Y, 'n_perm', 0);
    mi_ksg(idx)   = mi_norm;
    coexpr(idx)   = geom_mean_score(X, Y);
    label(idx)    = 0;
    pair_type{idx} = 'null';
end

% Replace NaN (small-n guard triggered) with 0
mi_ksg(isnan(mi_ksg))   = 0;
coexpr(isnan(coexpr))   = 0;

% ── Precision-Recall and AUROC ────────────────────────────────────────────────

% Sort by score descending
[~, ord_ksg]    = sort(mi_ksg,  'descend');
[~, ord_coexpr] = sort(coexpr,  'descend');

% AUROC via trapezoidal rule
auroc_ksg    = compute_auroc(mi_ksg,  label);
auroc_coexpr = compute_auroc(coexpr,  label);

% Precision @ k for k = n_causal
prec_ksg    = sum(label(ord_ksg(1:n_causal)))    / n_causal;
prec_coexpr = sum(label(ord_coexpr(1:n_causal))) / n_causal;

% Recall @ k
rec_ksg     = sum(label(ord_ksg(1:n_causal)))    / n_causal;
rec_coexpr  = sum(label(ord_coexpr(1:n_causal))) / n_causal;

% False discovery in dead-receiver set (the critical test)
% "How many dead-receiver pairs does each method put in its top-k?"
n_dead_in_top_ksg    = sum(strcmp(pair_type(ord_ksg(1:n_causal)),    'dead_receiver'));
n_dead_in_top_coexpr = sum(strcmp(pair_type(ord_coexpr(1:n_causal)), 'dead_receiver'));

if opt.verbose
    fprintf('=== Results ===\n\n');
    fprintf('  %-30s  %8s  %8s\n', 'Metric', 'KSG-MI', 'CoExpr');
    fprintf('  %-30s  %8.3f  %8.3f\n', 'AUROC',             auroc_ksg,    auroc_coexpr);
    fprintf('  %-30s  %8.3f  %8.3f\n', ['Precision@', num2str(n_causal)], prec_ksg, prec_coexpr);
    fprintf('  %-30s  %8.3f  %8.3f\n', ['Recall@',    num2str(n_causal)], rec_ksg,  rec_coexpr);
    fprintf('  %-30s  %8d  %8d\n', 'Dead-receiver in top-k',  ...
            n_dead_in_top_ksg, n_dead_in_top_coexpr);
    fprintf('\n');

    if n_dead_in_top_ksg < n_dead_in_top_coexpr
        fprintf('  INTERPRETATION: KSG-MI places fewer dead-receiver pairs in top-%d\n', n_causal);
        fprintf('  than co-expression (%d vs %d). The difference is the false-positive\n', ...
                n_dead_in_top_ksg, n_dead_in_top_coexpr);
        fprintf('  load that existing tools impose on experimental follow-up.\n\n');
    else
        fprintf('  NOTE: No advantage observed in this simulation configuration.\n');
        fprintf('  Try increasing beta (stronger causal signal) or zero_frac.\n\n');
    end
end

% ── Figure ────────────────────────────────────────────────────────────────────
if opt.plot
    c_bg    = [0.055, 0.067, 0.090];
    c_panel = [0.086, 0.106, 0.133];
    c_text  = [0.902, 0.929, 0.953];
    c_grid  = [0.176, 0.200, 0.224];
    c_caus  = [0.40, 0.76, 1.00];   % blue  — causal
    c_dead  = [1.00, 0.55, 0.40];   % coral — dead receiver
    c_null  = [0.60, 0.60, 0.60];   % grey  — null

    fig = figure('Name', 'Ground-Truth Validation', 'NumberTitle', 'off', ...
                 'Color', c_bg, 'Units', 'inches', 'Position', [0.5 0.5 13 8]);

    % Panel 1: KSG-MI score distribution by type
    ax1 = subplot(2, 3, 1);
    hold on;
    histogram(mi_ksg(strcmp(pair_type, 'causal')),        15, 'FaceColor', c_caus, 'FaceAlpha', 0.75);
    histogram(mi_ksg(strcmp(pair_type, 'dead_receiver')), 15, 'FaceColor', c_dead, 'FaceAlpha', 0.75);
    histogram(mi_ksg(strcmp(pair_type, 'null')),          15, 'FaceColor', c_null, 'FaceAlpha', 0.75);
    xlabel('KSG norm. MI', 'Color', c_text);
    ylabel('Count',        'Color', c_text);
    title('KSG MI by type', 'Color', c_text, 'FontWeight', 'bold');
    legend({'Causal','Dead receiver','Null'}, 'TextColor', c_text, 'Location', 'northeast');
    style_ax(ax1, c_panel, c_text, c_grid);

    % Panel 2: Co-expression score distribution by type
    ax2 = subplot(2, 3, 2);
    hold on;
    histogram(coexpr(strcmp(pair_type, 'causal')),        15, 'FaceColor', c_caus, 'FaceAlpha', 0.75);
    histogram(coexpr(strcmp(pair_type, 'dead_receiver')), 15, 'FaceColor', c_dead, 'FaceAlpha', 0.75);
    histogram(coexpr(strcmp(pair_type, 'null')),          15, 'FaceColor', c_null, 'FaceAlpha', 0.75);
    xlabel('Co-expression score', 'Color', c_text);
    title('Co-expression by type', 'Color', c_text, 'FontWeight', 'bold');
    style_ax(ax2, c_panel, c_text, c_grid);

    % Panel 3: Scatter KSG MI vs co-expression
    ax3 = subplot(2, 3, 3);
    hold on;
    scatter(coexpr(strcmp(pair_type, 'causal')),        mi_ksg(strcmp(pair_type, 'causal')),        30, c_caus, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter(coexpr(strcmp(pair_type, 'dead_receiver')), mi_ksg(strcmp(pair_type, 'dead_receiver')), 30, c_dead, 'filled', 'MarkerFaceAlpha', 0.8);
    scatter(coexpr(strcmp(pair_type, 'null')),          mi_ksg(strcmp(pair_type, 'null')),          20, c_null, 'filled', 'MarkerFaceAlpha', 0.6);
    xlabel('Co-expression score', 'Color', c_text);
    ylabel('KSG norm. MI',        'Color', c_text);
    title('Discordant pairs = dead receivers', 'Color', c_text, 'FontWeight', 'bold');
    style_ax(ax3, c_panel, c_text, c_grid);

    % Panel 4 & 5: ROC curves
    [tpr_ksg,    fpr_ksg]    = roc_curve(mi_ksg,  label);
    [tpr_coexpr, fpr_coexpr] = roc_curve(coexpr,  label);

    ax4 = subplot(2, 3, 4:5);
    hold on;
    plot(fpr_ksg,    tpr_ksg,    '-',  'Color', c_caus, 'LineWidth', 2.0, ...
         'DisplayName', sprintf('KSG MI  (AUROC = %.3f)', auroc_ksg));
    plot(fpr_coexpr, tpr_coexpr, '--', 'Color', c_dead, 'LineWidth', 1.8, ...
         'DisplayName', sprintf('Co-expression  (AUROC = %.3f)', auroc_coexpr));
    plot([0 1], [0 1], ':', 'Color', c_grid * 1.5, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    xlabel('False positive rate', 'Color', c_text, 'FontWeight', 'bold');
    ylabel('True positive rate',  'Color', c_text, 'FontWeight', 'bold');
    title('ROC curve: causal vs. non-causal L-R pairs', 'Color', c_text, 'FontWeight', 'bold');
    lh = legend('Location', 'southeast', 'FontSize', 9, 'TextColor', c_text);
    lh.Color = c_panel; lh.EdgeColor = c_grid;
    xlim([0 1]); ylim([0 1]);
    style_ax(ax4, c_panel, c_text, c_grid);

    % Panel 6: bar chart — dead-receivers in top-k
    ax6 = subplot(2, 3, 6);
    hold on;
    bar_data = [n_dead_in_top_coexpr, n_dead_in_top_ksg];
    b = bar([1 2], bar_data, 0.5);
    b.FaceColor = 'flat';
    b.CData = [c_dead; c_caus];
    set(ax6, 'XTick', [1 2], 'XTickLabel', {'Co-expression', 'KSG MI'});
    ylabel(sprintf('Dead-receiver pairs in top %d', n_causal), 'Color', c_text, 'FontWeight', 'bold');
    title('False positives: dead receivers promoted', 'Color', c_text, 'FontWeight', 'bold');
    yline(0, '-', 'Color', c_grid * 1.5, 'LineWidth', 1.0);
    style_ax(ax6, c_panel, c_text, c_grid);

    sgtitle({'NEUROSCOPE validation: KSG MI vs co-expression', ...
             sprintf('n=%d cells/pair, β=%.1f, zero-inflation=%.0f%%', n, beta, zf*100)}, ...
            'Color', c_text, 'FontSize', 11, 'FontWeight', 'bold');

    % Save alongside this script
    out_path = fullfile(fileparts(mfilename('fullpath')), '..', 'results', ...
                        'validation_ground_truth.png');
    try
        exportgraphics(fig, out_path, 'Resolution', 300, 'BackgroundColor', c_bg);
        if opt.verbose
            fprintf('  Figure saved: %s\n\n', out_path);
        end
    catch
        % results dir may not exist yet — not fatal
    end
end

% ── Pack results struct ───────────────────────────────────────────────────────
results.mi_ksg              = mi_ksg;
results.coexpr              = coexpr;
results.label               = label;
results.pair_type           = pair_type;
results.auroc_ksg           = auroc_ksg;
results.auroc_coexpr        = auroc_coexpr;
results.precision_ksg       = prec_ksg;
results.precision_coexpr    = prec_coexpr;
results.dead_in_topk_ksg    = n_dead_in_top_ksg;
results.dead_in_topk_coexpr = n_dead_in_top_coexpr;
results.params              = opt;

end  % main


% ─────────────────────────────────────────────────────────────────────────────
% Local helpers
% ─────────────────────────────────────────────────────────────────────────────

function v = zero_inflate(v, frac)
% Set a random fraction of elements to zero
mask = rand(size(v)) < frac;
v(mask) = 0;
end


function s = geom_mean_score(X, Y)
% Co-expression score: geometric mean of marginal means.
% This is what CellChat / NicheNet / CellPhoneDB use (cluster-level mean
% replaced here with sample mean, which is equivalent in single-cluster mode).
s = sqrt(mean(X) * mean(Y));
if isnan(s), s = 0; end
end


function auroc = compute_auroc(scores, labels)
% Compute area under the ROC curve via Mann-Whitney U / trapezoidal rule.
[~, ord] = sort(scores, 'descend');
labels_sorted = labels(ord);
n_pos = sum(labels);
n_neg = sum(~labels);
if n_pos == 0 || n_neg == 0
    auroc = 0.5;
    return;
end
tp = cumsum(labels_sorted);
fp = cumsum(~labels_sorted);
tpr = tp / n_pos;
fpr = fp / n_neg;
auroc = trapz([0; fpr], [0; tpr]);
end


function [tpr, fpr] = roc_curve(scores, labels)
[~, ord] = sort(scores, 'descend');
labels_s = labels(ord);
n_pos = sum(labels);
n_neg = sum(~labels);
tp  = [0; cumsum(labels_s)];
fp  = [0; cumsum(~labels_s)];
tpr = tp / max(n_pos, 1);
fpr = fp / max(n_neg, 1);
end


function style_ax(ax, c_panel, c_text, c_grid)
set(ax, 'Color', c_panel, 'XColor', c_text, 'YColor', c_text, ...
        'GridColor', c_grid, 'GridAlpha', 0.6, ...
        'XGrid', 'on', 'YGrid', 'on', 'TickDir', 'in', 'Box', 'on');
end
