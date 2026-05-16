function ranked = rank_communication_channels(condition_specific, condition_network, out_dir)
% RANK_COMMUNICATION_CHANNELS  Rank L-R pairs by information content.
%
% Replaces the GBM-specific therapeutic target ranking with a science-first
% ordering based on channel capacity (information transmission), statistical
% confidence, and channel topology classification.
%
% Ranking criteria (in priority order):
%   1. CondCapacity     — KSG channel capacity in bits (primary: information transmission)
%   2. FDRCapacity      — BH-corrected p-value for capacity enrichment (novel test)
%   3. FDR              — BH-corrected p-value for strength enrichment (traditional)
%   4. BootCI_lo        — Lower bootstrap CI bound (excludes zero-capacity channels)
%
% OUTPUT columns (superset of the old target_rankings for backward compat):
%   Ligand, Receptor, CondCapacity, RefCapacity, CondStrength, RefStrength,
%   Log2FC, FDR, FDRCapacity, BootCI_lo, BootCI_hi, Direction,
%   Pathway, ChannelType, InformationScore, CausalScore
%
% The InformationScore combines capacity + FDR into a single [0,1] composite.
% CausalScore is a placeholder (0) updated by infer_causal_interactions().

fprintf('Ranking communication channels by information content...\n');

if isempty(condition_specific) || height(condition_specific) == 0
    ranked = table();
    fprintf('  No channels to rank.\n\n');
    return;
end

n = height(condition_specific);

%% Composite information score
% InformationScore = 0.60 * normalized_capacity
%                  + 0.25 * (1 - FDRCapacity)   [novel test]
%                  + 0.15 * (1 - FDR)            [strength]
%
% Each term in [0,1]; total in [0,1].

cap_norm = condition_specific.CondCapacity ./ (max(condition_specific.CondCapacity) + 1e-10);
cap_norm(isnan(cap_norm)) = 0;

fdr_cap_score = 1 - min(condition_specific.FDRCapacity, 1);
fdr_str_score = 1 - min(condition_specific.FDR,         1);

InformationScore = 0.60 .* cap_norm ...
                 + 0.25 .* fdr_cap_score ...
                 + 0.15 .* fdr_str_score;

% CausalScore: filled by infer_causal_interactions(); placeholder here
CausalScore = zeros(n, 1);

%% Assemble output table
ranked = condition_specific;
ranked.InformationScore = InformationScore;
ranked.CausalScore      = CausalScore;

%% Sort: descending InformationScore, then ascending FDRCapacity
[~, idx] = sortrows([InformationScore, fdr_cap_score, fdr_str_score], [-1, -2, -3]);
ranked = ranked(idx, :);

%% Print summary
fprintf('  Total channels: %d\n', n);
fprintf('  High-capacity (>0.1 bits): %d\n', sum(condition_specific.CondCapacity > 0.1));
if height(ranked) > 0
    fprintf('  Top channel:  %s → %s\n', ranked.Ligand{1}, ranked.Receptor{1});
    fprintf('    Capacity:   %.4f bits\n', ranked.CondCapacity(1));
    fprintf('    FDR (cap):  %.4f\n',      ranked.FDRCapacity(1));
    fprintf('    Channel:    %s\n',         ranked.ChannelType{1});
    fprintf('    Pathway:    %s\n',         ranked.Pathway{1});
end
fprintf('\n');

% Save — use caller-supplied out_dir (dataset-specific) if provided
if nargin < 3 || isempty(out_dir)
    out_dir = 'results';
end
if ~exist(out_dir, 'dir'), out_dir = '.'; end
try
    writetable(ranked, fullfile(out_dir, 'channel_rankings.csv'));
    fprintf('  Saved: %s/channel_rankings.csv\n\n', out_dir);
catch ME
    fprintf('  Warning: save failed: %s\n\n', ME.message);
end

end
