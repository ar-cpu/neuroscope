function target_rankings = rank_therapeutic_targets(network, condition_specific, stats)
% RANK_THERAPEUTIC_TARGETS  Rank L-R channels by composite therapeutic score.
%
% Legacy; use rank_communication_channels for new code.
%
% Composite TherapeuticScore (weights: 0.50 / 0.25 / 0.15 / 0.10):
%   1. Capacity score     (0-1): normalised KSG channel capacity I(L;R)/H(L)
%   2. Significance score (0-1): 1 - FDR
%   3. Frequency score    (0-1): fraction of cluster-pairs where channel is active
%   4. Druggability score (0-1): heuristic match against tractable gene families
%
%   target_rankings = rank_therapeutic_targets(network, condition_specific)
%   target_rankings = rank_therapeutic_targets(network, condition_specific, stats)
%
% Inputs
%   network            - Communication network struct from build_communication_network()
%   condition_specific - Differential channel table from identify_differential_channels();
%                        must have columns: Ligand, Receptor, CondCapacity, FDR, CondFrequency
%   stats              - (optional) Accepted for signature compatibility; not used
%
% Output
%   target_rankings - Table sorted descending by TherapeuticScore; columns:
%       Ligand, Receptor, TherapeuticScore, DruggabilityScore,
%       ChannelCapacity, FDR, Frequency, Pathway

if nargin < 3, stats = []; end  %#ok<NASGU>

req_cols = {'Ligand','Receptor','CondCapacity','FDR','CondFrequency'};
for c = req_cols
    if ~any(strcmp(condition_specific.Properties.VariableNames, c{1}))
        error('rank_therapeutic_targets: condition_specific is missing column ''%s''.', c{1});
    end
end

n = height(condition_specific);
ligands    = condition_specific.Ligand;
receptors  = condition_specific.Receptor;
capacities = condition_specific.CondCapacity;
fdrs       = condition_specific.FDR;
freqs      = condition_specific.CondFrequency;

if any(strcmp(condition_specific.Properties.VariableNames, 'Pathway'))
    pathways = condition_specific.Pathway;
else
    pathways = repmat({'Other'}, n, 1);
end

% 1. Capacity score
cap_max = max(capacities);
if cap_max > 0
    cap_score = capacities / cap_max;
else
    cap_score = zeros(n, 1);
end

% 2. Significance score
sig_score = max(0, 1 - fdrs);

% 3. Frequency score
freq_max = max(freqs);
if freq_max > 0
    freq_score = freqs / freq_max;
else
    freq_score = zeros(n, 1);
end

% 4. Druggability score
drug_score = druggability_score(ligands, receptors);

% Composite: capacity dominates; druggability/frequency break ties
w_cap  = 0.50;
w_sig  = 0.25;
w_freq = 0.15;
w_drug = 0.10;

composite = w_cap  .* cap_score  + ...
            w_sig  .* sig_score  + ...
            w_freq .* freq_score + ...
            w_drug .* drug_score;

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


% Druggability heuristic: 1.0 if ligand/receptor is in a known tractable family, 0 otherwise

function scores = druggability_score(ligands, receptors)
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
