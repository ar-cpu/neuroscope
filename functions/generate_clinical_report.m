function generate_clinical_report(network, condition_specific, target_rankings, config)
% GENERATE_CLINICAL_REPORT  (legacy name — wrapper for GENERATE_SCIENCE_REPORT)
%
% This function is superseded by GENERATE_SCIENCE_REPORT, which produces the
% same content without disease-specific framing. New code should call
% GENERATE_SCIENCE_REPORT directly. This wrapper is kept for backwards
% compatibility with existing scripts.
%
% USAGE
%   generate_clinical_report(network, condition_specific, target_rankings, config)
%
% INPUTS
%   network            - Communication network struct from
%                        build_communication_network().
%   condition_specific - Differential channel table from
%                        identify_differential_channels().
%   target_rankings    - Channel ranking table from
%                        rank_communication_channels() or
%                        rank_therapeutic_targets().
%   config             - Config struct with at minimum:
%                          .dataset_name  (string)
%                          .output_dir    (string path)
%                          .fdr_threshold (scalar, default 0.05)
%
% See also: GENERATE_SCIENCE_REPORT

if nargin < 4 || isempty(config)
    config = struct();
end
if ~isfield(config, 'dataset_name'),  config.dataset_name  = 'Dataset';  end
if ~isfield(config, 'output_dir'),    config.output_dir    = 'results';   end
if ~isfield(config, 'fdr_threshold'), config.fdr_threshold = 0.05;        end

generate_science_report(network, condition_specific, target_rankings, config);

end
