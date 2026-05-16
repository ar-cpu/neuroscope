function condition_specific = identify_gbm_specific_interactions( ...
        condition_network, reference_network, ~, ~, ~, ~, n_permutations)
% IDENTIFY_GBM_SPECIFIC_INTERACTIONS  (legacy name — tissue-agnostic wrapper)
%
% Identifies communication channels significantly enriched in the condition
% relative to a reference, using the KSG channel-capacity framework. This
% function is a backwards-compatible wrapper around IDENTIFY_DIFFERENTIAL_CHANNELS,
% which supersedes it. New code should call IDENTIFY_DIFFERENTIAL_CHANNELS directly.
%
% USAGE
%   condition_specific = identify_gbm_specific_interactions( ...
%       condition_network, reference_network)
%   condition_specific = identify_gbm_specific_interactions( ...
%       condition_network, reference_network, expr_c, expr_r, genes, meta, n_perm)
%
% INPUTS
%   condition_network  - Network struct from build_communication_network()
%                        for the condition of interest.
%   reference_network  - Network struct for the reference/control condition.
%                        Pass [] or create_empty_network() for one-sample mode.
%   expr_c, expr_r     - (ignored) Expression matrices; accepted for signature
%                        compatibility, not used.
%   genes, meta        - (ignored) Gene names / metadata; accepted for
%                        signature compatibility, not used.
%   n_permutations     - Number of bootstrap resamples for CI estimation.
%                        Default: 200.
%
% OUTPUT
%   condition_specific - Table with one row per tested L-R pair. Same schema
%                        as IDENTIFY_DIFFERENTIAL_CHANNELS: columns Ligand,
%                        Receptor, CondStrength, RefStrength, Log2FC,
%                        CondCapacity, RefCapacity, CondFrequency,
%                        RefFrequency, PValue, FDR, PValueCapacity,
%                        FDRCapacity, BootCI_lo, BootCI_hi, Direction,
%                        Pathway, ChannelType.

if nargin < 7 || isempty(n_permutations)
    n_permutations = 200;
end

if isempty(reference_network)
    reference_network = create_empty_network();
end

condition_specific = identify_differential_channels( ...
    condition_network, reference_network, n_permutations);

end
