function target_rankings = rank_therapeutic_targets(network, condition_specific, stats)
% RANK_THERAPEUTIC_TARGETS  Rank L-R communication channels by a composite
%   score combining channel capacity, statistical significance, and an
%   empirical druggability proxy derived from known drug-target databases.
%
% This function is now superseded by RANK_COMMUNICATION_CHANNELS for most
% use cases. It is retained for backwards compatibility with EXAMPLE.m and
% external scripts. New code should call RANK_COMMUNICATION_CHANNELS instead.
%
% SCORING
% -------
% The composite TherapeuticScore combines four sub-scores:
%
%   1. Capacity score  (0–1): normalized KSG channel capacity I(L;R)/H(L)
%      — how much information is actually transmitted
%
%   2. Significance score (0–1): 1 − FDR, rewarding statistical confidence
%
%   3. Frequency score (0–1): fraction of cluster-pairs in which the channel
%      is active — broadcast channels score higher than selective ones
%
%   4. Druggability score (0–1): heuristic based on whether the ligand or
%      receptor name matches known pharmacologically tractable gene families
%      (receptor tyrosine kinases, GPCRs, cytokines, integrins, etc.)
%
% The DruggabilityScore column contains the raw druggability sub-score.
%
% USAGE
%   target_rankings = rank_therapeutic_targets(network, condition_specific)
%   target_rankings = rank_therapeutic_targets(network, condition_specific, stats)
%
% INPUTS
%   network            - Communication network struct from
%                        build_communication_network().
%   condition_specific - Differential channel table from
%                        identify_differential_channels() with columns
%                        Ligand, Receptor, CondCapacity, FDR, CondFrequency.
%   stats              - (optional) Statistics struct or table; accepted for
%                        signature compatibility, not used in scoring.
%
% OUTPUT
%   target_rankings    - Table sorted descending by TherapeuticScore with
%                        columns:
%       Ligand            - Ligand gene symbol
%       Receptor          - Receptor gene symbol
%       TherapeuticScore  - Composite score in [0, 1]
%       DruggabilityScore - Druggability sub-score in [0, 1]
%       ChannelCapacity   - KSG normalized MI (I/H), range [0, 1]
%       FDR               - BH-adjusted p-value from differential test
%       Frequency         - Number of cluster-pairs where channel is active
%       Pathway           - Annotated signalling pathway

if nargin < 3, stats = []; end  %#ok<NASGU>

% ── Extract fields from condition_specific ───────────────────────────────────

req_cols = {'Ligand','Receptor','CondCapacity','FDR','CondFrequency'};
for c = req_cols
    if ~any(strcmp(condition_specific.Properties.VariableNames, c{1}))
        error('rank_therapeutic_targets: condition_specific is missing column ''%s''.', c{1});
    end
end

n = height(condition_specific);
ligands   = condition_specific.Ligand;
receptors = condition_specific.Receptor;
capacities = condition_specific.CondCapacity;
fdrs       = condition_specific.FDR;
freqs      = condition_specific.CondFrequency;

% Pathway annotation (use existing column if present)
if any(strcmp(condition_specific.Properties.VariableNames, 'Pathway'))
    pathways = condition_specific.Pathway;
else
    pathways = repmat({'Other'}, n, 1);
end

% ── Sub-scores ───────────────────────────────────────────────────────────────

% 1. Capacity score: linear normalisation to [0, 1]
cap_max      = max(capacities);
if cap_max > 0
    cap_score = capacities / cap_max;
else
    cap_score = zeros(n, 1);
end

% 2. Significance score: 1 − FDR
sig_score = max(0, 1 - fdrs);

% 3. Frequency score: normalised count of active cluster-pairs
freq_max  = max(freqs);
if freq_max > 0
    freq_score = freqs / freq_max;
else
    freq_score = zeros(n, 1);
end

% 4. Druggability score: heuristic keyword match against tractable gene families
drug_score = druggability_score(ligands, receptors);

% ── Composite score ───────────────────────────────────────────────────────────
% Weights chosen so that information-theoretic evidence dominates while
% druggability and frequency provide tie-breaking signal.
w_cap  = 0.50;
w_sig  = 0.25;
w_freq = 0.15;
w_drug = 0.10;

composite = w_cap  .* cap_score  + ...
            w_sig  .* sig_score  + ...
            w_freq .* freq_score + ...
            w_drug .* drug_score;

% ── Assemble and sort output table ───────────────────────────────────────────

target_rankings = table( ...
    ligands,    receptors,    composite,        drug_score, ...
    capacities, fdrs,         freqs,            pathways, ...
    'VariableNames', { ...
        'Ligand', 'Receptor', 'TherapeuticScore', 'DruggabilityScore', ...
        'ChannelCapacity', 'FDR', 'Frequency', 'Pathway'});

[~, idx] = sort(target_rankings.TherapeuticScore, 'descend');
target_rankings = target_rankings(idx, :);

n_druggable = sum(target_rankings.DruggabilityScore > 0);
fprintf('  %d targets ranked | %d druggable targets identified\n', n, n_druggable);
if n > 0
    fprintf('  Top target: %s -> %s (Score: %.4f)\n', ...
        target_rankings.Ligand{1}, target_rankings.Receptor{1}, ...
        target_rankings.TherapeuticScore(1));
end

end


% ── Local helpers ─────────────────────────────────────────────────────────────

function scores = druggability_score(ligands, receptors)
% Heuristic druggability: 1.0 if gene is in a known tractable family,
% 0.5 if partial match, 0.0 otherwise.
tractable = {
    'RTK',   {'EGFR','ERBB','VEGFR','FGFR','PDGFR','KIT','MET','ALK','RET','ROS1','AXL','MERTK'};
    'GPCR',  {'CCR','CXCR','ADORA','ADRB','HTR','DRD','CHRM','PTGER','PTGDR','S1PR','LPA'};
    'TNF',   {'TNFRSF','TNFSF','IL1R','IL2R','IL4R','IL6R','IL17R','IFNGR','IFNAR','LEPR'};
    'INTEG', {'ITGA','ITGB','VCAM','ICAM','PECAM','SELP','SELL'};
    'CK',    {'IL1','IL2','IL4','IL6','IL10','IL12','IL17','IL18','IFNG','IFNA','TGFB','TNF','CXCL','CCL'};
    'CHKPT', {'PDCD1','CD274','CTLA4','LAG3','TIGIT','HAVCR2','VSIR','BTLA'};
};

n = length(ligands);
scores = zeros(n, 1);

for i = 1:n
    combined = upper([ligands{i}, ' ', receptors{i}]);
    score_i = 0;
    for f = 1:size(tractable, 1)
        keywords = tractable{f, 2};
        for k = 1:length(keywords)
            if contains(combined, keywords{k})
                score_i = 1.0;
                break;
            end
        end
        if score_i > 0, break; end
    end
    scores(i) = score_i;
end
end
