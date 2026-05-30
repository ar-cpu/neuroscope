function [mi_norm, mi, pvalue, n_eff, ci_lower, ci_upper] = estimate_channel_capacity(X, Y, varargin)
% ESTIMATE_CHANNEL_CAPACITY  Normalized mutual information between a ligand/receptor
%   pair using the KSG k-nearest-neighbour MI estimator.
%
%   mi_norm = I(X ; Y) / H(X)   in [0, 1]
%
% USAGE
%   [mi_norm, mi, pvalue, n_eff]              = estimate_channel_capacity(X, Y)
%   [mi_norm, mi, pvalue, n_eff, ci_l, ci_u] = estimate_channel_capacity(X, Y, 'n_bootstrap', 500)
%   [mi_norm, mi, pvalue, n_eff]              = estimate_channel_capacity(X, Y, 'k', 5, 'n_perm', 199)
%
% INPUTS — PAIRING CONTRACT
%   X and Y MUST be paired: X(i) and Y(i) are from the same biological unit
%   (e.g. one per donor). Passing unpaired cell-level vectors destroys the
%   joint distribution and makes the MI estimate meaningless.
%
%   X           [n x 1]  ligand values, one per paired observation
%   Y           [n x 1]  receptor values, same pairing as X
%   'k'         KSG neighbourhood size (default 0 = adaptive, range 3-10)
%   'n_perm'    permutations for p-value (default 199; 0 to skip)
%   'max_n'     max observations before subsampling (default 300)
%   'seed'      RNG seed for reproducibility
%   'n_bootstrap' bootstrap resamples for 95% CI (default 500; 0 to skip)
%
% OUTPUTS
%   mi_norm   [0,1] normalized MI = I(X;Y)/H(X)
%   mi        mutual information in bits
%   pvalue    permutation p-value; NaN if n_perm=0
%   n_eff     effective sample size used
%   ci_lower  2.5th percentile of bootstrap distribution (NaN if not requested)
%   ci_upper  97.5th percentile of bootstrap distribution (NaN if not requested)

p = inputParser();
addParameter(p, 'k',           0);
addParameter(p, 'n_perm',      199);
addParameter(p, 'max_n',       300);
addParameter(p, 'seed',        []);
addParameter(p, 'n_bootstrap', 500);
parse(p, varargin{:});
k_arg       = p.Results.k;
n_perm      = p.Results.n_perm;
max_n       = p.Results.max_n;
n_bootstrap = p.Results.n_bootstrap;

X = real(double(full(X(:))));
Y = real(double(full(Y(:))));

ok = isfinite(X) & isfinite(Y);
X  = X(ok);
Y  = Y(ok);

MIN_N = 10;
if length(X) < MIN_N
    mi_norm = NaN;  mi = NaN;  pvalue = NaN;  n_eff = 0;
    ci_lower = NaN;  ci_upper = NaN;
    return;
end

% CRITICAL: X and Y are paired — subsample using a single shared index vector
% to preserve the joint distribution.
if ~isempty(p.Results.seed), rng(p.Results.seed); end
n = length(X);
if n > max_n
    idx = randsample(n, max_n);
    X   = X(idx);
    Y   = Y(idx);
    n   = max_n;
end

if k_arg > 0
    k = k_arg;
else
    k = max(3, min(10, floor(n^0.4)));
end

mi     = ksg_mi(X, Y, k);
H_x    = marginal_entropy_kl(X, k);
n_eff = n;
if H_x < 1e-10
    mi_norm = 0;  pvalue = 1;  ci_lower = 0;  ci_upper = 0;
    return;
end
mi_norm = max(0, min(1, mi / H_x));

if n_perm > 0
    null_mi = zeros(n_perm, 1);
    for b = 1:n_perm
        Y_perm    = Y(randperm(n));
        null_mi(b) = ksg_mi(X, Y_perm, k);
    end
    pvalue = (sum(null_mi >= mi) + 1) / (n_perm + 1);
else
    pvalue = NaN;
end

if n_bootstrap > 0 && nargout > 4
    boot_mi_norm = zeros(n_bootstrap, 1);
    for b = 1:n_bootstrap
        idx   = randsample(n, n, true);
        Xb    = X(idx);
        Yb    = Y(idx);
        k_b   = max(3, min(10, floor(n^0.4)));
        mi_b  = ksg_mi(Xb, Yb, k_b);
        H_b   = marginal_entropy_kl(Xb, k_b);
        if H_b < 1e-10
            boot_mi_norm(b) = 0;
        else
            boot_mi_norm(b) = max(0, min(1, mi_b / H_b));
        end
    end

    % BCa bias-correction (Efron 1987): z0 = median bias of bootstrap distribution
    prop_below = mean(boot_mi_norm < mi_norm);
    prop_below = max(1e-6, min(1-1e-6, prop_below));
    z0 = -sqrt(2) * erfcinv(2 * prop_below);

    % Acceleration constant a via jackknife (skewness of influence function)
    jack_est = zeros(n, 1);
    k_j = max(3, min(10, floor((n-1)^0.4)));
    for i = 1:n
        idx_j   = [1:i-1, i+1:n];
        mi_j    = ksg_mi(X(idx_j), Y(idx_j), k_j);
        H_j     = marginal_entropy_kl(X(idx_j), k_j);
        if H_j < 1e-10
            jack_est(i) = 0;
        else
            jack_est(i) = max(0, min(1, mi_j / H_j));
        end
    end
    jack_mean = mean(jack_est);
    diff3     = sum((jack_mean - jack_est).^3);
    diff2     = sum((jack_mean - jack_est).^2);
    if diff2 < 1e-14
        a = 0;
    else
        a = diff3 / (6 * diff2^1.5);
    end
    a = max(-0.5, min(0.5, a));

    % Adjusted percentile indices (Efron 1987, eq. 22)
    z_lo = -1.959964;
    z_hi =  1.959964;

    arg_lo = z0 + (z0 + z_lo) / (1 - a * (z0 + z_lo));
    arg_hi = z0 + (z0 + z_hi) / (1 - a * (z0 + z_hi));

    alpha_lo = 0.5 * erfc(-arg_lo / sqrt(2));
    alpha_hi = 0.5 * erfc(-arg_hi / sqrt(2));

    alpha_lo = max(0.001, min(0.499, alpha_lo));
    alpha_hi = max(0.501, min(0.999, alpha_hi));

    boot_sorted = sort(boot_mi_norm);
    lo_idx  = max(1,           round(alpha_lo * n_bootstrap));
    hi_idx  = min(n_bootstrap, round(alpha_hi * n_bootstrap));

    ci_lower = boot_sorted(lo_idx);
    ci_upper = boot_sorted(hi_idx);

    % Widen symmetrically if CI still doesn't bracket the point estimate
    if ci_lower > mi_norm
        ci_lower = max(0, mi_norm - (ci_upper - mi_norm));
    end
    if ci_upper < mi_norm
        ci_upper = min(1, mi_norm + (mi_norm - ci_lower));
    end
else
    ci_lower = NaN;
    ci_upper = NaN;
end

end  % main function


% KSG MI estimator — Algorithm 1 from Kraskov et al. 2004
function mi = ksg_mi(X, Y, k)
n  = length(X);
XY = [X(:), Y(:)];

% Robust standardization: IQR-scale rather than std, which is dominated by
% outliers in zero-inflated scRNA-seq data. Fallback to std when IQR ~ 0.
for col = 1:2
    v   = real(XY(:, col));
    iq  = iqr(v);
    if iq < 1e-12
        iq = std(v);
        if iq < 1e-12, iq = 1; end
    end
    XY(:, col) = (v - median(v)) / iq;
end
Xs = XY(:,1);
Ys = XY(:,2);

psi_sum = 0;
for i = 1:n
    d = max(abs(XY - XY(i,:)), [], 2);
    d(i) = Inf;
    sorted_d = sort(d);
    eps_i = sorted_d(k);

    nx = sum(abs(Xs - Xs(i)) < eps_i) - 1;
    ny = sum(abs(Ys - Ys(i)) < eps_i) - 1;
    nx = max(nx, 0);
    ny = max(ny, 0);

    psi_sum = psi_sum + psi(double(nx + 1)) + psi(double(ny + 1));
end

mi = max(0, (psi(k) + psi(n) - psi_sum/n) / log(2));
end


% Marginal entropy via Kozachenko-Leonenko estimator with zero-mass correction
function H = marginal_entropy_kl(X, k)
% Mixed-distribution decomposition for zero-inflated data (v5.2):
%   H(X) = H_bern(p0) + (1-p0) * H_continuous(X | X > 0)
% where p0 = P(X=0). Exact for any distribution with a point mass at 0.
n  = length(X);

zero_mask = (X == 0);
p0        = sum(zero_mask) / n;
n_pos     = sum(~zero_mask);

if p0 >= 1 - 1e-10
    H = 0;
    return;
end

if p0 > 1e-10
    H_bern = -p0 * log2(p0) - (1 - p0) * log2(1 - p0);

    if n_pos >= 2
        k_pos = max(2, min(k, floor(n_pos^0.4)));
        H_pos = kl_entropy_1d(X(~zero_mask), k_pos);
        H = H_bern + (1 - p0) * H_pos;
    else
        H = H_bern;
    end
else
    H = kl_entropy_1d(X, k);
end

H = max(0, H);
end


% Kozachenko-Leonenko 1-D entropy estimator on strictly continuous data
function H = kl_entropy_1d(X, k)
n  = length(X);
Xs = sort(X);
H  = 0;
for i = 1:n
    lo = max(1, i - k);
    hi = min(n, i + k);
    left_d  = Xs(i) - Xs(lo:i-1);
    right_d = Xs(i+1:hi) - Xs(i);
    all_d   = sort([left_d(:); right_d(:)]);
    if isempty(all_d)
        d = 1e-12;
    elseif length(all_d) >= k
        d = all_d(k);
    else
        d = all_d(end);
    end
    if d < 1e-12, d = 1e-12; end
    H = H + log2(2 * d);
end
H = H / n + log2(exp(1)) * (psi(n) - psi(k));
H = max(0, H);
end
