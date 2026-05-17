function condition_specific = identify_gbm_specific_interactions( ...
        condition_network, reference_network, ~, ~, ~, ~, n_permutations)
% IDENTIFY_GBM_SPECIFIC_INTERACTIONS  Backwards-compatible wrapper around
%   IDENTIFY_DIFFERENTIAL_CHANNELS. New code should call that function directly.
%   Extra arguments (expr_c, expr_r, genes, meta) are accepted but ignored.

if nargin < 7 || isempty(n_permutations)
    n_permutations = 200;
end

if isempty(reference_network)
    reference_network = create_empty_network();
end

condition_specific = identify_differential_channels( ...
    condition_network, reference_network, n_permutations);

end
