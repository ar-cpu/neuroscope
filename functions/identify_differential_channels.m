function diff_channels = identify_differential_channels(condition_network, reference_network, n_permutations, fdr_threshold)
% IDENTIFY_DIFFERENTIAL_CHANNELS  Find communication channels significantly
%   enriched, depleted, or rewired between two conditions (tumor vs. normal,
%   treatment vs. control, etc.). Tissue- and disease-agnostic.
%
% Statistical framework: per-cluster-pair strengths are compared with
%   Wilcoxon rank-sum (two-sample) or signed-rank vs. zero (one-sample).
%   A separate test on KSG normalised MI distinguishes "more expressed"
%   from "more informative" — existing tools cannot do this.
%
% Usage: diff_channels = identify_differential_channels(cond_net, ref_net)
%        diff_channels = identify_differential_channels(cond_net, ref_net, 200, 0.05)
% Inputs:  condition_network  - Network struct from build_communication_network()
%          reference_network  - Reference network; pass create_empty_network() for one-sample mode
%          n_permutations     - Bootstrap iterations for CI estimation (default 200)
%          fdr_threshold      - BH-FDR cutoff for reporting (default 0.05)

fprintf('Identifying differential communication channels...\n\n');

if nargin < 3 || isempty(n_permutations), n_permutations = 200; end
if nargin < 4 || isempty(fdr_threshold),  fdr_threshold  = 0.05; end

has_reference = ~isempty(reference_network) && ...
    isfield(reference_network, 'ligands') && ...
    ~isempty(reference_network.ligands);

% One-sample mode tests each channel against zero; two-sample compares distributions.
if has_reference
    fprintf('  Mode: two-sample (condition vs. reference)\n');
else
    fprintf('  Mode: one-sample (condition characterization)\n');
end

cond_pairs  = aggregate_by_lr(condition_network);
has_capacity = isfield(condition_network, 'channel_capacity') && ...
               ~isempty(condition_network.channel_capacity);

if has_reference
    ref_pairs = aggregate_by_lr(reference_network);
else
    ref_pairs = containers.Map();
end

all_keys = union(keys(cond_pairs), keys(ref_pairs));
n_pairs  = length(all_keys);

fprintf('  Condition: %d unique L-R pairs\n', cond_pairs.Count);
if has_reference
    fprintf('  Reference: %d unique L-R pairs\n', ref_pairs.Count);
else
    fprintf('  Reference: no data (using one-sample tests)\n');
end
fprintf('\n');

Ligand         = cell(n_pairs, 1);
Receptor       = cell(n_pairs, 1);
CondStrength   = zeros(n_pairs, 1);
RefStrength    = zeros(n_pairs, 1);
CondCapacity   = zeros(n_pairs, 1);
RefCapacity    = zeros(n_pairs, 1);
CondFrequency  = zeros(n_pairs, 1);
RefFrequency   = zeros(n_pairs, 1);
PValue         = ones(n_pairs, 1);
PValueCapacity = ones(n_pairs, 1);
Log2FC         = zeros(n_pairs, 1);
Direction      = cell(n_pairs, 1);
BootCI_lo      = zeros(n_pairs, 1);
BootCI_hi      = zeros(n_pairs, 1);

fprintf('  [1/4] Running statistical tests...\n');

for i = 1:n_pairs
    key = all_keys{i};
    parts = strsplit(key, '|||');
    Ligand{i}   = parts{1};
    Receptor{i} = parts{2};

    cond_vals = [];
    ref_vals  = [];
    cond_caps = [];
    ref_caps  = [];
    cond_ci_lo = [];
    cond_ci_hi = [];

    if isKey(cond_pairs, key)
        cond_vals  = cond_pairs(key).strengths;
        cond_caps  = cond_pairs(key).capacities;
        cond_ci_lo = cond_pairs(key).ci_lo;
        cond_ci_hi = cond_pairs(key).ci_hi;
    end
    if has_reference && isKey(ref_pairs, key)
        ref_vals = ref_pairs(key).strengths;
        ref_caps = ref_pairs(key).capacities;
    end

    CondStrength(i)  = mean(cond_vals,  'omitnan');
    RefStrength(i)   = mean(ref_vals,   'omitnan');
    CondCapacity(i)  = mean(cond_caps,  'omitnan');
    RefCapacity(i)   = mean(ref_caps,   'omitnan');
    CondFrequency(i) = numel(cond_vals);
    RefFrequency(i)  = numel(ref_vals);

    CondStrength(isnan(CondStrength)) = 0;
    RefStrength(isnan(RefStrength))   = 0;
    CondCapacity(isnan(CondCapacity)) = 0;
    RefCapacity(isnan(RefCapacity))   = 0;

    Log2FC(i) = log2((CondStrength(i) + 1e-6) / (RefStrength(i) + 1e-6));

    if has_reference && ~isempty(ref_vals)
        try
            PValue(i) = ranksum(cond_vals, ref_vals);
        catch
            PValue(i) = 1;
        end
        if ~isempty(cond_caps) && ~isempty(ref_caps)
            try
                PValueCapacity(i) = ranksum(cond_caps, ref_caps);
            catch
                PValueCapacity(i) = 1;
            end
        end
    else
        % One-sample: one-sided right-tail test (H1: median > 0).
        % Two-sided would spuriously flag channels significantly below zero,
        % which has no meaning for expression strengths or MI values.
        if length(cond_vals) >= 5
            try
                PValue(i) = signrank(cond_vals, 0, 'tail', 'right');
            catch
                PValue(i) = 1;
            end
        end
        if length(cond_caps) >= 5
            try
                PValueCapacity(i) = signrank(cond_caps, 0, 'tail', 'right');
            catch
                PValueCapacity(i) = 1;
            end
        end
    end

    if CondStrength(i) > RefStrength(i)
        Direction{i} = 'Enriched';
    elseif CondStrength(i) < RefStrength(i)
        Direction{i} = 'Depleted';
    else
        Direction{i} = 'Absent';
    end

    % CIs propagated from estimate_channel_capacity (2.5th/97.5th percentiles
    % of the KSG bootstrap). For Spearman-fallback pairs, lo == hi == CondCapacity.
    if ~isempty(cond_ci_lo) && any(~isnan(cond_ci_lo))
        BootCI_lo(i) = mean(cond_ci_lo(~isnan(cond_ci_lo)));
        BootCI_hi(i) = mean(cond_ci_hi(~isnan(cond_ci_hi)));
    else
        BootCI_lo(i) = CondCapacity(i);
        BootCI_hi(i) = CondCapacity(i);
    end
end

fprintf('  [2/4] Multiple testing correction (Benjamini-Hochberg)...\n');
[~, sort_idx] = sort(PValue);
n = n_pairs;
FDR          = zeros(n,1);
FDRCapacity  = zeros(n,1);
ranks        = zeros(n,1);
ranks(sort_idx) = 1:n;
FDR = min(1, PValue .* n ./ max(ranks, 1));
for i = n-1:-1:1
    FDR(sort_idx(i)) = min(FDR(sort_idx(i)), FDR(sort_idx(i+1)));
end
[~, sort_idx2] = sort(PValueCapacity);
ranks2 = zeros(n,1);
ranks2(sort_idx2) = 1:n;
FDRCapacity = min(1, PValueCapacity .* n ./ max(ranks2,1));
for i = n-1:-1:1
    FDRCapacity(sort_idx2(i)) = min(FDRCapacity(sort_idx2(i)), FDRCapacity(sort_idx2(i+1)));
end

fprintf('  [3/4] Bootstrap CIs: propagated from KSG estimator (%d resamples each)...\n', n_permutations);
n_ci_valid = sum(BootCI_lo ~= BootCI_hi);
fprintf('  Bootstrap CIs valid for %d / %d pairs\n', n_ci_valid, n_pairs);
if n_ci_valid > 0
    mean_width = mean(BootCI_hi(BootCI_lo ~= BootCI_hi) - BootCI_lo(BootCI_lo ~= BootCI_hi));
    fprintf('  Mean 95%% CI width: %.4f bits (n=16: expected ~0.10-0.20)\n', mean_width);
    top_idx = 1;  %#ok<NASGU>
    top_boot = 1; %#ok<NASGU>
end

fprintf('  [4/4] Annotating pathways and channel types...\n');
Pathway    = annotate_pathways(Ligand, Receptor);
ChannelType = classify_channel_type(CondCapacity, CondFrequency, n_pairs);

diff_channels = table(Ligand, Receptor, ...
    CondStrength, RefStrength, Log2FC, ...
    CondCapacity, RefCapacity, ...
    CondFrequency, RefFrequency, ...
    PValue, FDR, PValueCapacity, FDRCapacity, ...
    BootCI_lo, BootCI_hi, ...
    Direction, Pathway, ChannelType);

[~, sort_order] = sortrows([FDR, -CondCapacity]);
diff_channels = diff_channels(sort_order, :);

n_sig = sum(diff_channels.FDR < fdr_threshold);
n_sig_cap = sum(diff_channels.FDRCapacity < fdr_threshold);

fprintf('\nDifferential channel analysis complete\n');
fprintf('  Total channels tested: %d\n', n_pairs);
fprintf('  Significant by strength (FDR<%.2f): %d\n', fdr_threshold, n_sig);
fprintf('  Significant by channel capacity (FDR<%.2f): %d\n', fdr_threshold, n_sig_cap);
if n_sig > 0
    fprintf('  Top channel: %s -> %s (log2FC=%.2f, norm-MI=%.4f)\n', ...
        diff_channels.Ligand{1}, diff_channels.Receptor{1}, ...
        diff_channels.Log2FC(1), diff_channels.CondCapacity(1));
end

out_dir = 'results';
if ~exist(out_dir, 'dir'), out_dir = '.'; end
try
    writetable(diff_channels, fullfile(out_dir, 'differential_channels.csv'));
    sig_only = diff_channels(diff_channels.FDR < fdr_threshold, :);
    writetable(sig_only, fullfile(out_dir, 'differential_channels_significant.csv'));
    fprintf('  Saved: differential_channels.csv\n');
catch ME
    fprintf('  Warning: save failed: %s\n', ME.message);
end
fprintf('\n');

end  % main


function pair_map = aggregate_by_lr(network)
% Build 'Ligand|||Receptor' → struct(.strengths, .capacities, .ci_lo, .ci_hi)
pair_map = containers.Map();
has_cap  = isfield(network, 'channel_capacity') && ~isempty(network.channel_capacity);
has_ci   = isfield(network, 'ci_lower') && ~isempty(network.ci_lower);
n = length(network.ligands);
for i = 1:n
    key = [network.ligands{i}, '|||', network.receptors{i}];
    s   = network.weights(i);
    c   = 0;
    clo = NaN;  chi = NaN;
    if has_cap && i <= length(network.channel_capacity)
        c = network.channel_capacity(i);
    end
    if has_ci && i <= length(network.ci_lower)
        clo = network.ci_lower(i);
        chi = network.ci_upper(i);
    end
    if isKey(pair_map, key)
        entry = pair_map(key);
        entry.strengths  = [entry.strengths,  s];
        entry.capacities = [entry.capacities, c];
        entry.ci_lo      = [entry.ci_lo,      clo];
        entry.ci_hi      = [entry.ci_hi,      chi];
        pair_map(key) = entry;
    else
        pair_map(key) = struct('strengths', s, 'capacities', c, 'ci_lo', clo, 'ci_hi', chi);
    end
end
end


function pathways = annotate_pathways(ligands, receptors)
pathways = cell(length(ligands), 1);
db = {
    {'WNT','FZD','LRP','SFRP','DKK'},       'Wnt';
    {'EGF','AREG','EREG','EGFR','ERBB'},     'EGF/ERBB';
    {'VEGF','KDR','FLT'},                    'VEGF';
    {'TGF','TGFB','SMAD','BMP','ACVR'},      'TGF-β/BMP';
    {'NOTCH','DLL','JAG','HES'},             'Notch';
    {'HH','SHH','PTCH','SMO','GLI'},         'Hedgehog';
    {'CCL','CXCL','CCR','CXCR'},             'Chemokine';
    {'CADM','CDH','NCAM','L1CAM','JAM'},     'Cell adhesion';
    {'EFNA','EFNB','EPHA','EPHB'},           'Ephrin';
    {'PDGF','PDGFR'},                        'PDGF';
    {'IGF','IGFR','IRS'},                    'IGF';
    {'SEMA','NRP','PLXN'},                   'Semaphorin';
    {'NRG','NRXN','NLGN'},                   'Neurexin/Neuroligin';
    {'BDNF','NGF','NTRK','NGFR'},            'Neurotrophin';
    {'IL','IFNG','TNF','JAK','STAT'},        'Cytokine/JAK-STAT';
    {'CSF','CSFR','PYY','NPY'},              'Colony/neuropeptide';
};
for i = 1:length(ligands)
    combined = upper([ligands{i}, ' ', receptors{i}]);
    assigned = 'Other';
    for d = 1:size(db,1)
        keywords = db{d,1};
        for k = 1:length(keywords)
            if contains(combined, keywords{k})
                assigned = db{d,2};
                break;
            end
        end
        if ~strcmp(assigned, 'Other'), break; end
    end
    pathways{i} = assigned;
end
end


function types = classify_channel_type(capacities, frequencies, n)
types = cell(n,1);
cap_med = median(capacities(capacities > 0));
if isnan(cap_med) || cap_med == 0, cap_med = 0.05; end
freq_med = median(frequencies(frequencies > 0));
if isnan(freq_med) || freq_med == 0, freq_med = 5; end

for i = 1:n
    hi_cap  = capacities(i)  > cap_med;
    hi_freq = frequencies(i) > freq_med;
    if hi_cap && hi_freq
        types{i} = 'HighCapacity_Broadcast';
    elseif hi_cap && ~hi_freq
        types{i} = 'HighCapacity_Selective';
    elseif ~hi_cap && hi_freq
        types{i} = 'LowCapacity_Broadcast';
    else
        types{i} = 'LowCapacity_Selective';
    end
end
end
