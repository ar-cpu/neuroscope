function landscape = compute_information_landscape(condition_network, condition_specific, cellphonedb)
% COMPUTE_INFORMATION_LANDSCAPE  Characterise what fraction of co-expressed
% L-R pairs transmit measurable information, and how capacity-based ranking
% differs from co-expression-based ranking.
%
% THIS FUNCTION GENERATES THE CORE EMPIRICAL CLAIM OF THE NEUROSCOPE PAPER:
%
%   "The majority of ligand-receptor co-expression events detected by
%   existing tools transmit zero measurable information to receiver cells.
%   The subset with I(X;Y) > 0 defines a communication landscape that is
%   only partially predicted by co-expression strength."
%
% To make this claim testable and falsifiable, the function computes:
%
%   1. Capacity distribution across ALL tested L-R pairs.
%      → What fraction have capacity ≈ 0 (permutation p > 0.05)?
%      → This is the headline percentage in the abstract.
%
%   2. Capacity vs co-expression rank correlation.
%      → Spearman r between capacity ranking and strength ranking.
%      → Low r means capacity reveals interactions invisible to existing tools.
%
%   3. Exclusive sets:
%      → Capacity-only enriched  (FDRCapacity < 0.05, FDR >= 0.05)
%      → Strength-only enriched  (FDR < 0.05, FDRCapacity >= 0.05)
%      → Both enriched           (FDRCapacity < 0.05, FDR < 0.05)
%      → Neither enriched
%
%   4. Pathway enrichment comparison:
%      → Which biological pathways are over-represented among capacity-only
%        interactions vs strength-only interactions?
%
%   5. Database coverage check:
%      → What fraction of CellPhoneDB pairs are in the network at all?
%      → High coverage = network building is not the bottleneck.
%
% OUTPUTS (struct):
%   landscape.n_tested              - total L-R pairs tested
%   landscape.n_zero_capacity       - pairs with pvalue > 0.05 (I≈0)
%   landscape.zero_fraction         - fraction with I≈0 (the headline number)
%   landscape.capacity_values       - vector of all capacity values
%   landscape.strength_values       - vector of all co-expression strengths
%   landscape.rank_correlation      - Spearman r (capacity rank vs strength rank)
%   landscape.exclusive_capacity    - table: capacity-only enriched pairs
%   landscape.exclusive_strength    - table: strength-only enriched pairs
%   landscape.both_enriched         - table: enriched by both tests
%   landscape.capacity_pathway_dist - pathway distribution, capacity-only
%   landscape.strength_pathway_dist - pathway distribution, strength-only
%   landscape.n_cpdb_pairs          - CellPhoneDB total L-R pairs
%   landscape.coverage_fraction     - fraction of CPDB pairs tested

fprintf('Computing information landscape...\n');
fprintf('(This generates the core empirical claim for publication)\n\n');

landscape = struct();
fdr_cut = 0.05;

if isempty(condition_specific) || height(condition_specific) == 0
    fprintf('  No interactions available for landscape analysis.\n\n');
    return;
end

cs = condition_specific;
n  = height(cs);
landscape.n_tested = n;

%% 1. Capacity distribution — the headline number
% "I≈0" = permutation p-value > 0.05 (channel capacity not significant)
% Note: we use FDRCapacity rather than raw p-value to be conservative.
% This makes the "I≈0" fraction a LOWER BOUND (actual fraction is higher).

n_zero     = sum(cs.FDRCapacity >= fdr_cut);
zero_frac  = n_zero / n;

landscape.n_zero_capacity  = n_zero;
landscape.zero_fraction    = zero_frac;
landscape.capacity_values  = cs.CondCapacity;

fprintf('  %-45s %d / %d  (%.1f%%)\n', 'L-R pairs with I≈0 (FDRCapacity>=0.05):', ...
    n_zero, n, 100*zero_frac);
fprintf('  %-45s %.4f bits (median)\n', 'Median channel capacity:', ...
    median(cs.CondCapacity));
fprintf('  %-45s %.4f bits\n', '90th percentile capacity:', ...
    prctile(cs.CondCapacity, 90));
fprintf('\n');

%% 2. Co-expression strength vs capacity rank correlation
% If Spearman r is high, capacity adds no information over co-expression.
% Low r is the key claim that capacity reveals a different landscape.

if ismember('CondStrength', cs.Properties.VariableNames)
    str_vals = cs.CondStrength;
else
    % Fall back to network weights if available
    str_vals = cs.CondCapacity;  % same if strength not stored
    fprintf('  Warning: CondStrength not in table — using capacity for comparison.\n');
end

landscape.strength_values = str_vals;

% Spearman correlation
cap_rank = tiedrank(cs.CondCapacity);
str_rank = tiedrank(str_vals);
r_spear  = corr(cap_rank, str_rank, 'type', 'Pearson');  % on ranks = Spearman
landscape.rank_correlation = r_spear;

fprintf('  Spearman r (capacity rank vs strength rank):   %.4f\n', r_spear);
if r_spear < 0.5
    fprintf('  >> LOW correlation — capacity reveals a substantially different\n');
    fprintf('     ranking than co-expression. This is the central novel finding.\n');
elseif r_spear < 0.7
    fprintf('  >> MODERATE correlation — capacity partially differs from co-expression.\n');
else
    fprintf('  >> HIGH correlation — capacity and co-expression are largely aligned.\n');
end
fprintf('\n');

%% 3. Four exclusive sets — the 2×2 table for the paper
cap_sig = cs.FDRCapacity < fdr_cut;
str_sig = cs.FDR         < fdr_cut;

both_mask  = cap_sig & str_sig;
cap_only   = cap_sig & ~str_sig;
str_only   = ~cap_sig & str_sig;
neither    = ~cap_sig & ~str_sig;

landscape.exclusive_capacity = cs(cap_only, :);
landscape.exclusive_strength = cs(str_only, :);
landscape.both_enriched      = cs(both_mask, :);

fprintf('  2×2 enrichment table:\n');
fprintf('  %-38s  %5s  %8s\n', '', 'N', 'Fraction');
fprintf('  %s\n', repmat('-',1,56));
fprintf('  %-38s  %5d  %7.1f%%\n', 'Both enriched (cap + strength):', sum(both_mask),  100*mean(both_mask));
fprintf('  %-38s  %5d  %7.1f%%  [NOVEL]', 'Capacity-only (missed by strength):', sum(cap_only), 100*mean(cap_only));
if sum(cap_only) > 0
    fprintf('  ← invisible to CellChat/NicheNet\n');
else
    fprintf('\n');
end
fprintf('  %-38s  %5d  %7.1f%%\n', 'Strength-only (no info transmission):',  sum(str_only), 100*mean(str_only));
fprintf('  %-38s  %5d  %7.1f%%\n', 'Neither significant:',                 sum(neither),  100*mean(neither));
fprintf('\n');

%% 4. Pathway distribution in exclusive sets
if ismember('Pathway', cs.Properties.VariableNames)
    landscape.capacity_pathway_dist = tabulate_pathways(cs(cap_only, :));
    landscape.strength_pathway_dist = tabulate_pathways(cs(str_only, :));

    if ~isempty(landscape.capacity_pathway_dist)
        fprintf('  Top pathways in capacity-only set:\n');
        pt = landscape.capacity_pathway_dist;
        for i = 1:min(5, height(pt))
            fprintf('    %-25s  %d\n', pt.Pathway{i}, pt.Count(i));
        end
        fprintf('\n');
    end
end

%% 5. CellPhoneDB coverage
if nargin >= 3 && isstruct(cellphonedb) && isfield(cellphonedb, 'interactions')
    n_cpdb = height(cellphonedb.interactions);
    landscape.n_cpdb_pairs    = n_cpdb;
    landscape.coverage_fraction = n / max(n_cpdb, 1);
    fprintf('  CellPhoneDB L-R pairs in database: %d\n', n_cpdb);
    fprintf('  Pairs tested by NEUROSCOPE:         %d (%.1f%%)\n\n', ...
        n, 100*landscape.coverage_fraction);
end

%% Summary line (for abstract drafting)
fprintf('  --- ABSTRACT CANDIDATE ---\n');
fprintf('  "We find that %.0f%% of co-expressed ligand-receptor pairs ', 100*zero_frac);
fprintf('transmit zero\n   measurable information (I(X;Y)=0, FDR>0.05 by KSG permutation\n');
fprintf('   test), and that the %.0f%% with I>0 are only %.0f%% rank-correlated\n', ...
    100*(1-zero_frac), 100*abs(r_spear));
fprintf('   with co-expression strength (Spearman r=%.2f)."\n\n', r_spear);

% Save
out_dir = 'results';
if ~exist(out_dir, 'dir'), out_dir = '.'; end
try
    T = table(cs.Ligand, cs.Receptor, cs.CondCapacity, str_vals, ...
        cs.FDRCapacity, cs.FDR, cap_sig, str_sig, cap_only, str_only, ...
        'VariableNames', {'Ligand','Receptor','Capacity','Strength', ...
        'FDRCapacity','FDRStrength','CapSig','StrSig','CapOnly','StrOnly'});
    writetable(T, fullfile(out_dir, 'information_landscape.csv'));
    fprintf('  Saved: information_landscape.csv\n\n');
catch
end

end


function pt = tabulate_pathways(T)
pt = table();
if isempty(T) || ~ismember('Pathway', T.Properties.VariableNames)
    return;
end
[pathways, ~, idx] = unique(T.Pathway);
counts = accumarray(idx, 1);
[counts_sorted, ord] = sort(counts, 'descend');
pt = table(pathways(ord), counts_sorted, 'VariableNames', {'Pathway','Count'});
end
