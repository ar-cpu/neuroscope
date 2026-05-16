function promiscuity_results = discover_receptor_promiscuity(gbm_network, gbm_specific)
% DISCOVER_RECEPTOR_PROMISCUITY  Identify receptors that bind multiple
%   distinct ligands (hub receptors) and quantify the competitive dynamics
%   among those ligands for receptor occupancy.
%
% Biological motivation
% ---------------------
% Receptors that interact with many different ligands present a particular
% challenge for targeted therapy: blocking a single ligand may be insufficient
% if alternative ligands can substitute. This function characterizes receptor
% promiscuity from three complementary angles:
%
%   1. Ligand count and statistical enrichment — how many distinct ligands
%      bind each receptor, tested against a hypergeometric null (analytical
%      approximation) with BH-FDR correction.
%   2. Competition index — Shannon entropy of the ligand weight distribution,
%      normalized to [0, 1]. A score of 0 indicates a single dominant ligand
%      (monopoly); a score of 1 indicates all ligands contribute equally
%      (maximum competition).
%   3. Relative resistance — a heuristic composite of ligand count and
%      competition index that ranks receptors by difficulty of therapeutic
%      blockade. This is not a validated mechanistic model; treat it as an
%      exploratory prioritization tool.
%
% Accepts both table and struct network formats; column/field roles are
% inferred from names containing 'ligand'/'source' and 'receptor'/'target'.
%
% Usage:
%   promiscuity_results = discover_receptor_promiscuity(gbm_network, gbm_specific)
%
% Inputs:
%   gbm_network  - Communication network as a struct (with .ligands,
%                  .receptors, and optionally .weights fields) or a table
%                  with analogously named columns.
%   gbm_specific - Table of condition-specific interactions. Not used
%                  quantitatively in this function but passed for API
%                  consistency with other discovery functions.
%
% Outputs:
%   promiscuity_results.receptor_data   - Per-receptor table with columns:
%         Receptor, LigandCount, Ligands, TotalWeight, Expected, ZScore,
%         PValue, FDR, CompetitionIndex, DominantLigand, BackupLigand,
%         RelativeResistance, IsDrugTarget
%   promiscuity_results.n_hub_receptors - Receptors binding ≥ 3 ligands
%   promiscuity_results.n_significant   - Receptors with FDR < 0.05
%   promiscuity_results.n_drug_targets  - Hub receptors with known drugs

fprintf('Receptor Promiscuity Analysis\n\n');

%% Extract network fields (handles both struct and table input)
fprintf('Extracting network structure...\n');

% Identify format and extract columns
if istable(gbm_network)
    fprintf('  Network format: table\n');
    col_names = gbm_network.Properties.VariableNames;
    fprintf('  Columns: %s\n', strjoin(col_names, ', '));
    
    % Find ligand column
    lig_col = '';
    for i = 1:length(col_names)
        if contains(lower(col_names{i}), {'ligand', 'source', 'from'})
            lig_col = col_names{i};
            break;
        end
    end
    
    % Find receptor column
    rec_col = '';
    for i = 1:length(col_names)
        if contains(lower(col_names{i}), {'receptor', 'target', 'to'})
            rec_col = col_names{i};
            break;
        end
    end
    
    % Find weight column (optional)
    weight_col = '';
    for i = 1:length(col_names)
        if contains(lower(col_names{i}), {'weight', 'score', 'strength'})
            weight_col = col_names{i};
            break;
        end
    end
    
    if isempty(lig_col) || isempty(rec_col)
        error('Cannot identify ligand/receptor columns. Available: %s', ...
            strjoin(col_names, ', '));
    end
    
    fprintf('  Using: %s -> %s\n', lig_col, rec_col);
    
    ligands = gbm_network.(lig_col);
    receptors = gbm_network.(rec_col);
    
    if ~isempty(weight_col)
        fprintf('  Weights: %s\n', weight_col);
        weights = gbm_network.(weight_col);
    else
        fprintf('  No weight column - using counts\n');
        weights = ones(height(gbm_network), 1);
    end
    
elseif isstruct(gbm_network)
    fprintf('  Network format: struct\n');
    fields = fieldnames(gbm_network);
    fprintf('  Fields: %s\n', strjoin(fields, ', '));
    
    % Find ligand field
    lig_field = '';
    for i = 1:length(fields)
        if contains(lower(fields{i}), {'ligand', 'source', 'from'})
            lig_field = fields{i};
            break;
        end
    end
    
    % Find receptor field
    rec_field = '';
    for i = 1:length(fields)
        if contains(lower(fields{i}), {'receptor', 'target', 'to'})
            rec_field = fields{i};
            break;
        end
    end
    
    if isempty(lig_field) || isempty(rec_field)
        error('Cannot identify ligand/receptor fields. Available: %s', ...
            strjoin(fields, ', '));
    end
    
    ligands = gbm_network.(lig_field);
    receptors = gbm_network.(rec_field);
    
    % Look for weights
    if isfield(gbm_network, 'weights') || isfield(gbm_network, 'weight')
        if isfield(gbm_network, 'weights')
            weights = gbm_network.weights;
        else
            weights = gbm_network.weight;
        end
    else
        weights = ones(size(ligands));
    end
else
    error('Network must be table or struct, got: %s', class(gbm_network));
end

fprintf('  Total interactions: %d\n', length(ligands));
fprintf('  Unique ligands: %d\n', length(unique(ligands)));
fprintf('  Unique receptors: %d\n', length(unique(receptors)));
fprintf('\n');

%% Phase 1: Count distinct ligands per receptor
fprintf('Phase 1: Counting ligands per receptor...\n');

unique_receptors = unique(receptors);
n_receptors = length(unique_receptors);

receptor_data = table('Size', [0, 4], ...
    'VariableTypes', {'string', 'double', 'string', 'double'}, ...
    'VariableNames', {'Receptor', 'LigandCount', 'Ligands', 'TotalWeight'});

for rec_idx = 1:n_receptors
    receptor = unique_receptors{rec_idx};
    
    rec_mask = strcmp(receptors, receptor);
    rec_ligands = unique(ligands(rec_mask));
    n_ligands = length(rec_ligands);
    
    if n_ligands >= 3
        total_weight = sum(weights(rec_mask));
        new_row = {receptor, n_ligands, strjoin(rec_ligands, '; '), total_weight};
        receptor_data = [receptor_data; new_row];
    end
end

receptor_data = sortrows(receptor_data, 'LigandCount', 'descend');

fprintf('  Receptors analyzed: %d\n', n_receptors);
fprintf('  Promiscuous (>=3 ligands): %d\n', height(receptor_data));
fprintf('\n');

%% Phase 2: Statistical significance (analytical null approximation)
fprintf('Phase 2: Testing significance (analytical approximation)...\n');

n_interactions = length(ligands);
n_unique_ligands = length(unique(ligands));

receptor_data.Expected = zeros(height(receptor_data), 1);
receptor_data.ZScore = zeros(height(receptor_data), 1);
receptor_data.PValue = ones(height(receptor_data), 1);

for rec_idx = 1:height(receptor_data)
    receptor_name = receptor_data.Receptor{rec_idx};
    observed = receptor_data.LigandCount(rec_idx);
    
    % Count interactions for this receptor
    rec_mask = strcmp(receptors, receptor_name);
    n_rec_int = sum(rec_mask);
    
    if n_rec_int > 0
        % Expected unique ligands under random sampling
        % E[unique] = n_total * (1 - (1 - 1/n_total)^n_samples)
        p_single = 1.0 / n_unique_ligands;
        expected = n_unique_ligands * (1 - (1 - p_single)^n_rec_int);
        
        % Variance (approximation)
        var_unique = n_unique_ligands * (1 - (1-p_single)^n_rec_int) * ...
                     ((1-p_single)^n_rec_int);
        
        receptor_data.Expected(rec_idx) = expected;
        
        if var_unique > 0
            z = (observed - expected) / sqrt(var_unique);
            receptor_data.ZScore(rec_idx) = z;
            receptor_data.PValue(rec_idx) = 1 - normcdf(z);
        end
    end
end

fprintf('  Significance computed\n');
fprintf('  Note: Uses analytical approximation (not exact permutation)\n');
fprintf('\n');

%% Phase 3: Benjamini-Hochberg FDR correction
fprintf('Phase 3: Multiple testing correction...\n');

% Benjamini-Hochberg FDR
[sorted_p, sort_idx] = sort(receptor_data.PValue);
m = length(sorted_p);
bh_crit = ((1:m)' / m) * 0.05;

% Find adjusted p-values properly
fdr_vals = ones(m, 1);
for i = 1:m
    fdr_vals(i) = min(sorted_p(i) * m / i, 1);
end

% Enforce monotonicity
for i = (m-1):-1:1
    fdr_vals(i) = min(fdr_vals(i), fdr_vals(i+1));
end

receptor_data.FDR = zeros(height(receptor_data), 1);
receptor_data.FDR(sort_idx) = fdr_vals;

n_sig = sum(receptor_data.FDR < 0.05);
fprintf('  Significant (FDR < 0.05): %d\n', n_sig);
fprintf('\n');

%% Phase 4: Competition index (normalized Shannon entropy of ligand strengths)
fprintf('Phase 4: Computing competition index...\n');

receptor_data.CompetitionIndex = zeros(height(receptor_data), 1);
receptor_data.DominantLigand = strings(height(receptor_data), 1);
receptor_data.BackupLigand = strings(height(receptor_data), 1);

for rec_idx = 1:height(receptor_data)
    receptor_name = receptor_data.Receptor{rec_idx};
    rec_mask = strcmp(receptors, receptor_name);
    
    rec_ligs = ligands(rec_mask);
    rec_weights = weights(rec_mask);
    
    % Aggregate by ligand
    unique_ligs = unique(rec_ligs);
    lig_strengths = zeros(length(unique_ligs), 1);
    
    for lig_idx = 1:length(unique_ligs)
        lig_mask = strcmp(rec_ligs, unique_ligs{lig_idx});
        lig_strengths(lig_idx) = sum(rec_weights(lig_mask));
    end
    
    % Avoid division by zero
    total_strength = sum(lig_strengths);
    if total_strength > 0
        % Shannon entropy as competition measure
        probs = lig_strengths / total_strength;
        entropy = -sum(probs .* log2(probs + 1e-10));
        max_entropy = log2(length(probs));
        
        if max_entropy > 0
            receptor_data.CompetitionIndex(rec_idx) = entropy / max_entropy;
        end
        
        % Identify top 2 ligands
        [~, sort_idx] = sort(lig_strengths, 'descend');
        receptor_data.DominantLigand{rec_idx} = unique_ligs{sort_idx(1)};
        if length(sort_idx) >= 2
            receptor_data.BackupLigand{rec_idx} = unique_ligs{sort_idx(2)};
        end
    end
end

fprintf('  Mean competition: %.3f (0=monopoly, 1=balanced)\n', ...
    mean(receptor_data.CompetitionIndex));
fprintf('\n');

%% Phase 5: Relative resistance scoring (heuristic)
fprintf('Phase 5: Computing relative resistance scores...\n');

% Heuristic: more ligands + higher competition = harder to block
receptor_data.RelativeResistance = zeros(height(receptor_data), 1);

for rec_idx = 1:height(receptor_data)
    n_ligs = receptor_data.LigandCount(rec_idx);
    comp = receptor_data.CompetitionIndex(rec_idx);
    
    % Simple multiplicative model (heuristic, not mechanistic)
    receptor_data.RelativeResistance(rec_idx) = n_ligs * (1 + comp);
end

% Normalize to [0, 1]
max_res = max(receptor_data.RelativeResistance);
if max_res > 0
    receptor_data.RelativeResistance = receptor_data.RelativeResistance / max_res;
end

fprintf('  Mean resistance score: %.3f\n', mean(receptor_data.RelativeResistance));
fprintf('  Note: Heuristic score, not a validated mechanistic model\n');
fprintf('\n');

%% Phase 6: Known drug target annotation
fprintf('Phase 6: Drug target annotation...\n');

drug_targets = {'EGFR', 'PDGFRA', 'PDGFRB', 'MET', 'FGFR', 'KDR', 'FLT', 'TGFBR'};

receptor_data.IsDrugTarget = false(height(receptor_data), 1);

for rec_idx = 1:height(receptor_data)
    rec_name = receptor_data.Receptor{rec_idx};
    for drug_idx = 1:length(drug_targets)
        if contains(rec_name, drug_targets{drug_idx})
            receptor_data.IsDrugTarget(rec_idx) = true;
            break;
        end
    end
end

n_targets = sum(receptor_data.IsDrugTarget);
fprintf('  Drug targets: %d\n', n_targets);
fprintf('\n');

%% EXPORT
fprintf('Exporting results...\n');

if ~exist('results/novel_discoveries', 'dir')
    mkdir('results/novel_discoveries');
end

writetable(receptor_data, 'results/novel_discoveries/receptor_promiscuity.csv');
fprintf('  receptor_promiscuity.csv\n\n');

%% Package
promiscuity_results = struct();
promiscuity_results.receptor_data = receptor_data;
promiscuity_results.n_hub_receptors = height(receptor_data);
promiscuity_results.n_significant = n_sig;
promiscuity_results.n_drug_targets = n_targets;

%% SUMMARY
fprintf('Promiscuity analysis complete\n\n');

fprintf('RESULTS:\n');
fprintf('  Hub receptors (>=3 ligands): %d\n', height(receptor_data));
fprintf('  Significant (FDR<0.05): %d\n', n_sig);
fprintf('  Known drug targets: %d\n', n_targets);

if height(receptor_data) > 0
    fprintf('\nTOP 5 HUB RECEPTORS:\n');
    for i = 1:min(5, height(receptor_data))
        fprintf('  %d. %s: %d ligands', i, ...
            receptor_data.Receptor{i}, ...
            receptor_data.LigandCount(i));
        
        if receptor_data.FDR(i) < 0.05
            fprintf(' (FDR=%.4f)', receptor_data.FDR(i));
        end
        
        if receptor_data.IsDrugTarget(i)
            fprintf(' [DRUG TARGET]');
        end
        fprintf('\n');
        
        fprintf('     Competition: %.3f | Resistance: %.3f\n', ...
            receptor_data.CompetitionIndex(i), ...
            receptor_data.RelativeResistance(i));
        
        if ~isempty(receptor_data.DominantLigand{i})
            fprintf('     Dominant: %s', receptor_data.DominantLigand{i});
            if ~isempty(receptor_data.BackupLigand{i})
                fprintf(' | Backup: %s', receptor_data.BackupLigand{i});
            end
            fprintf('\n');
        end
    end
end

fprintf('\nAnalysis complete\n\n');

end