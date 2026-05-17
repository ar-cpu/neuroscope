function generate_clinical_report(network, condition_specific, target_rankings, config)
% GENERATE_CLINICAL_REPORT  Backwards-compatible wrapper for GENERATE_SCIENCE_REPORT.
%   New code should call generate_science_report directly.

if nargin < 4 || isempty(config)
    config = struct();
end
if ~isfield(config, 'dataset_name'),  config.dataset_name  = 'Dataset';  end
if ~isfield(config, 'output_dir'),    config.output_dir    = 'results';   end
if ~isfield(config, 'fdr_threshold'), config.fdr_threshold = 0.05;        end

generate_science_report(network, condition_specific, target_rankings, config);

end
