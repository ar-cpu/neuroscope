function interventions = propose_interventions(target_rankings, network, config)
% PROPOSE_INTERVENTIONS  Generate ranked intervention strategies from
%   communication channel rankings.
%
% For each high-scoring target, three intervention types are evaluated:
%
%   Inhibition   — block ligand secretion or receptor binding. Relevant when
%                  a channel transmits pro-disease information (e.g., a
%                  suppressive checkpoint or a tumour-promoting growth signal).
%
%   Activation   — augment signalling through an understimulated channel.
%                  Relevant for protective or anti-inflammatory pathways.
%
%   Combination  — simultaneous targeting of ligand and receptor. Scored as
%                  the geometric mean of the two individual strategy scores,
%                  with an additive synergy bonus when both sub-units score
%                  above 0.5.
%
% Intervention priority is a linear combination of:
%   - TherapeuticScore from rank_therapeutic_targets / rank_communication_channels
%   - DruggabilityScore
%   - A channel-capacity bonus (high MI → stronger evidence of actual transmission)
%
% USAGE
%   interventions = propose_interventions(target_rankings)
%   interventions = propose_interventions(target_rankings, network, config)
%
% INPUTS
%   target_rankings - Table from rank_therapeutic_targets() or
%                     rank_communication_channels(). Must contain columns:
%                     Ligand, Receptor, TherapeuticScore (or InformationScore),
%                     DruggabilityScore, ChannelCapacity (or CondCapacity).
%   network         - (optional) Communication network struct; used only for
%                     frequency-based filtering. Ignored when empty.
%   config          - (optional) Struct with fields:
%                       .n_top           Max targets to expand (default: 20)
%                       .min_score       Minimum TherapeuticScore to include
%                                        (default: 0.1)
%                       .output_dir      Directory for CSV export (default: 'results')
%                       .dataset_name    String prefix for output files
%
% OUTPUT
%   interventions - Table sorted by Priority with columns:
%       Ligand, Receptor, Strategy, Priority, TherapeuticScore,
%       DruggabilityScore, ChannelCapacity, Rationale

if nargin < 2, network = []; end
if nargin < 3, config  = struct(); end

if ~isfield(config, 'n_top'),        config.n_top        = 20;       end
if ~isfield(config, 'min_score'),    config.min_score    = 0.1;      end
if ~isfield(config, 'output_dir'),   config.output_dir   = 'results'; end
if ~isfield(config, 'dataset_name'), config.dataset_name = 'Dataset'; end

% ── Normalise column names ────────────────────────────────────────────────────

cols = target_rankings.Properties.VariableNames;
score_col = 'TherapeuticScore';
if ~any(strcmp(cols, score_col))
    if any(strcmp(cols, 'InformationScore'))
        score_col = 'InformationScore';
    else
        error('propose_interventions: target_rankings must contain TherapeuticScore or InformationScore.');
    end
end

cap_col = 'ChannelCapacity';
if ~any(strcmp(cols, cap_col))
    if any(strcmp(cols, 'CondCapacity')),   cap_col = 'CondCapacity';
    else,                                   cap_col = score_col;  end
end

drug_col = 'DruggabilityScore';
has_drug = any(strcmp(cols, drug_col));

% ── Select top targets ────────────────────────────────────────────────────────

scores   = target_rankings.(score_col);
keep_idx = scores >= config.min_score;
top_tbl  = target_rankings(keep_idx, :);
if height(top_tbl) > config.n_top
    top_tbl = top_tbl(1:config.n_top, :);
end

n_targets = height(top_tbl);
if n_targets == 0
    warning('propose_interventions: no targets meet the minimum score threshold (%.2f).', ...
        config.min_score);
    interventions = table();
    return;
end

% ── Build intervention rows ───────────────────────────────────────────────────

n_rows = n_targets * 3;  % Inhibition, Activation, Combination per target
Ligand          = cell(n_rows, 1);
Receptor        = cell(n_rows, 1);
Strategy        = cell(n_rows, 1);
Priority        = zeros(n_rows, 1);
TherapeuticScore= zeros(n_rows, 1);
DruggabilityScore = zeros(n_rows, 1);
ChannelCapacity = zeros(n_rows, 1);
Rationale       = cell(n_rows, 1);

row = 0;
for i = 1:n_targets
    lig  = top_tbl.Ligand{i};
    rec  = top_tbl.Receptor{i};
    ts   = top_tbl.(score_col)(i);
    cap  = top_tbl.(cap_col)(i);
    ds   = 0;
    if has_drug, ds = top_tbl.(drug_col)(i); end

    cap_bonus = min(1, cap * 2);  % bonus for high channel capacity

    % Strategy 1: Inhibition
    row = row + 1;
    p_inh = 0.6 * ts + 0.3 * ds + 0.1 * cap_bonus;
    Ligand{row}           = lig;
    Receptor{row}         = rec;
    Strategy{row}         = 'Inhibition';
    Priority(row)         = p_inh;
    TherapeuticScore(row) = ts;
    DruggabilityScore(row)= ds;
    ChannelCapacity(row)  = cap;
    Rationale{row}        = sprintf(['Block %s→%s axis. '...
        'Channel capacity %.3f bits; MI-confirmed information transmission '...
        'justifies targeting even at high per-cell sparsity.'], lig, rec, cap);

    % Strategy 2: Activation
    row = row + 1;
    p_act = 0.5 * ts + 0.2 * ds + 0.3 * cap_bonus;
    Ligand{row}           = lig;
    Receptor{row}         = rec;
    Strategy{row}         = 'Activation';
    Priority(row)         = p_act;
    TherapeuticScore(row) = ts;
    DruggabilityScore(row)= ds;
    ChannelCapacity(row)  = cap;
    Rationale{row}        = sprintf(['Augment %s→%s signalling. '...
        'High capacity (%.3f) indicates robust communication; '...
        'activation may amplify protective or regenerative effects.'], lig, rec, cap);

    % Strategy 3: Combination (ligand + receptor co-targeting)
    row = row + 1;
    synergy   = double(ts > 0.5 && ds > 0.5) * 0.05;
    p_combo   = sqrt(p_inh * p_act) + synergy;
    Ligand{row}           = lig;
    Receptor{row}         = rec;
    Strategy{row}         = 'Combination';
    Priority(row)         = p_combo;
    TherapeuticScore(row) = ts;
    DruggabilityScore(row)= ds;
    ChannelCapacity(row)  = cap;
    Rationale{row}        = sprintf(['Dual blockade of %s (ligand) and %s (receptor). '...
        'Combination reduces resistance risk when channel capacity (%.3f) '...
        'exceeds 0.3, indicating strong inter-donor co-variation.'], lig, rec, cap);
end

% ── Assemble, sort, export ────────────────────────────────────────────────────

interventions = table(Ligand, Receptor, Strategy, Priority, ...
    TherapeuticScore, DruggabilityScore, ChannelCapacity, Rationale);

[~, idx] = sort(interventions.Priority, 'descend');
interventions = interventions(idx, :);

% Save
out_dir = config.output_dir;
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
try
    writetable(interventions, fullfile(out_dir, 'intervention_strategies.csv'));
    high_priority = interventions(interventions.Priority >= 0.5, :);
    writetable(high_priority, fullfile(out_dir, 'intervention_strategies_high_priority.csv'));
    fprintf('  Saved: intervention_strategies.csv (%d strategies, %d high-priority)\n', ...
        height(interventions), height(high_priority));
catch ME
    fprintf('  Warning: could not save intervention strategies: %s\n', ME.message);
end

end
