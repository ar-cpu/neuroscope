function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report)
%FDR_BH - Benjamini & Hochberg False Discovery Rate correction
%
% Executes the Benjamini & Hochberg (1995) and Benjamini & Yekutieli (2001)
% procedure for controlling the false discovery rate (FDR) of a family of
% hypothesis tests. FDR is the expected proportion of rejected hypotheses
% that are mistakenly rejected (i.e., the null hypothesis is actually true
% for those tests). FDR is a somewhat less conservative/more powerful method
% for correcting for multiple comparisons than procedures like Bonferroni
% correction that provide strong control of the family-wise error rate
% (i.e., the probability that one or more null hypotheses are mistakenly
% rejected).
%
% Usage:
%  [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report)
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - 'pdep' (default) or 'dep'. If 'pdep,' the original
%             Benjamini & Hochberg FDR procedure is used, which is guaranteed
%             to be accurate if the individual tests are independent or
%             positively dependent (e.g., Gaussian variables that are
%             positively correlated or independent). If 'dep,' the FDR
%             procedure described in Benjamini & Yekutieli (2001) is used,
%             which is guaranteed to be accurate for any test dependency
%             structure (e.g., if the tests are negatively dependent).
%   report  - If 'yes', a brief summary of the results is displayed. If 'no'
%             (default), no summary is displayed.
%
% Outputs:
%   h         - A binary vector or matrix of the same size as pvals. If
%               the ith element of h is 1, then the test that produced the
%               ith p-value in pvals is significant (i.e., the null hypothesis
%               of the test is rejected).
%   crit_p    - All uncorrected p-values less than or equal to crit_p are
%               significant (i.e., their null hypotheses are rejected). If
%               no p-values are significant, crit_p=0.
%   adj_ci_cvrg - The FCR-adjusted BH confidence interval coverage. This is
%                 equal to the desired false coverage rate minus the false
%                 discovery rate.
%   adj_p      - All adjusted p-values less than or equal to q are significant
%                (i.e., their null hypotheses are rejected). Note, adjusted
%                p-values can be greater than 1.
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
% Author: David Groppe
% Kutaslab, San Diego
% Modified for NeuroScope project

if nargin < 1
    error('You need to provide a vector or matrix of p-values.');
end

if nargin < 2 || isempty(q)
    q = 0.05;
end

if nargin < 3 || isempty(method)
    method = 'pdep';
end

if nargin < 4 || isempty(report)
    report = 'no';
end

% Check that q is valid
if q <= 0 || q >= 1
    error('q must be between 0 and 1.');
end

% Store original matrix dimensions
s = size(pvals);
if isvector(pvals)
    pvals = pvals(:);
end

% Make sure method is recognized
method = lower(method);
if ~strcmpi(method, 'pdep') && ~strcmpi(method, 'dep')
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end

% Get number of tests
m = numel(pvals);

% Sort p-values
[sorted_pvals, sort_ids] = sort(pvals(:));

% Create vector of p-value ranks
[~, ranks] = sort(sort_ids);

% Default to all null hypotheses being accepted (h=0)
h = zeros(m, 1);
adj_p = ones(m, 1);

if strcmpi(method, 'pdep')
    % BH procedure for independent or positively dependent tests
    thresh = (1:m)' * q / m;
    wtd_p = m * sorted_pvals ./ (1:m)';
    
elseif strcmpi(method, 'dep')
    % BH procedure for any dependency structure
    denom = sum(1 ./ (1:m));
    thresh = (1:m)' * q / (m * denom);
    wtd_p = m * denom * sorted_pvals ./ (1:m)';
end

% Find greatest p-value that is less than or equal to its critical value
rej = sorted_pvals <= thresh;
max_id = find(rej, 1, 'last');

if isempty(max_id)
    % No p-values are significant
    crit_p = 0;
    adj_ci_cvrg = NaN;
else
    % At least one p-value is significant
    h(sort_ids(1:max_id)) = 1;
    crit_p = sorted_pvals(max_id);
    
    % Calculate adjusted confidence interval coverage
    if nargout > 2
        adj_ci_cvrg = 1 - thresh(max_id);
    end
end

% Compute adjusted p-values
if nargout > 3
    % Compute adjusted p-values
    adj_p_sorted = zeros(m, 1);
    adj_p_sorted(m) = sorted_pvals(m);
    
    for i = m-1:-1:1
        adj_p_sorted(i) = min(adj_p_sorted(i+1), wtd_p(i));
    end
    
    % Put adjusted p-values back in original order
    adj_p(sort_ids) = adj_p_sorted;
    
    % Ensure adjusted p-values don't exceed 1
    adj_p = min(adj_p, 1);
end

% Reshape to original dimensions
h = reshape(h, s);
adj_p = reshape(adj_p, s);

% Report results if requested
if strcmpi(report, 'yes')
    n_sig = sum(h(:));
    fprintf('FDR Correction Results:\n');
    fprintf('  Total tests: %d\n', m);
    fprintf('  Significant tests: %d (%.1f%%)\n', n_sig, 100*n_sig/m);
    fprintf('  Critical p-value: %.6f\n', crit_p);
    fprintf('  Method: %s\n', method);
end

end