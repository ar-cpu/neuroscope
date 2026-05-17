function atlas = build_communication_atlas(networks, condition_specifics, dataset_names, varargin)
% BUILD_COMMUNICATION_ATLAS  Multi-dataset communication atlas: conserved and
%   specialized signaling programs across biological conditions or tissues.
%
% Usage:
%   atlas = build_communication_atlas(networks, condition_specifics, dataset_names)
%   atlas = build_communication_atlas(..., 'capacity_threshold', 0.05, 'verbose', true)
%
% Inputs:
%   networks            - {N×1} cell array of network structs (fields: .ligands, .receptors,
%                         optionally .channel_capacity, .weights)
%   condition_specifics - {N×1} cell array of tables with columns Ligand, Receptor, FDR.
%                         Pass {} or empty tables if unavailable.
%   dataset_names       - {N×1} cell array of dataset name strings
%
% Optional parameters:
%   'capacity_threshold'       - Min channel capacity to count as active (default 0.05)
%   'conservation_cosine_high' - Cosine threshold for conserved NMF programs (default 0.6)
%   'conservation_cosine_low'  - Below this = specialized program (default 0.3)
%   'fdr_threshold'            - FDR cutoff for condition-specific channels (default 0.05)
%   'verbose'                  - Print progress (default true)
%
% Output:
%   atlas.lr_table, .conserved_programs, .specialized_programs, .diversity_index,
%   .n_datasets, .dataset_names, .cross_dataset_matrix

p = inputParser();
p.addRequired('networks');
p.addRequired('condition_specifics');
p.addRequired('dataset_names');
p.addParameter('capacity_threshold',       0.05,  @isnumeric);
p.addParameter('conservation_cosine_high', 0.60,  @isnumeric);
p.addParameter('conservation_cosine_low',  0.30,  @isnumeric);
p.addParameter('fdr_threshold',            0.05,  @isnumeric);
p.addParameter('verbose',                  true,  @islogical);
p.parse(networks, condition_specifics, dataset_names, varargin{:});

cap_thresh    = p.Results.capacity_threshold;
cos_high      = p.Results.conservation_cosine_high;
cos_low       = p.Results.conservation_cosine_low;
fdr_thresh    = p.Results.fdr_threshold;
verbose       = p.Results.verbose;

if ~iscell(networks),           networks = {networks};                  end
if ~iscell(condition_specifics), condition_specifics = {condition_specifics}; end
if ~iscell(dataset_names),      dataset_names = {dataset_names};       end

N = length(networks);

if verbose
    fprintf('build_communication_atlas:\n');
    fprintf('  Datasets: %d\n', N);
    for i = 1:N
        fprintf('    [%d] %s\n', i, dataset_names{i});
    end
end

if N < 2
    fprintf('  Warning: only %d dataset provided. Atlas requires >= 2.\n', N);
end

% Step 1: Cross-dataset L-R pair alignment.
if verbose, fprintf('\n--- Step 1: Cross-dataset L-R pair alignment ---\n'); end

all_lr_pairs = {};
for d = 1:N
    net = networks{d};
    if isempty(net) || ~isfield(net, 'ligands'), continue; end

    ligs = net.ligands(:);
    recs = net.receptors(:);

    if ~iscell(ligs), ligs = cellstr(num2str(ligs)); end
    if ~iscell(recs), recs = cellstr(num2str(recs)); end

    pairs_d = strcat(ligs, '__', recs);
    all_lr_pairs = union(all_lr_pairs, pairs_d);
end

n_lr_total = length(all_lr_pairs);

if verbose
    fprintf('  Unique L-R pairs across all datasets: %d\n', n_lr_total);
end

lr_ligands   = cell(n_lr_total, 1);
lr_receptors = cell(n_lr_total, 1);
for i = 1:n_lr_total
    parts = strsplit(all_lr_pairs{i}, '__');
    lr_ligands{i}   = parts{1};
    lr_receptors{i} = parts{2};
end

cross_dataset_matrix = zeros(n_lr_total, N);

for d = 1:N
    net = networks{d};
    if isempty(net) || ~isfield(net, 'ligands'), continue; end

    ligs_d = net.ligands(:);
    recs_d = net.receptors(:);
    if ~iscell(ligs_d), ligs_d = cellstr(num2str(ligs_d)); end
    if ~iscell(recs_d), recs_d = cellstr(num2str(recs_d)); end

    pairs_d = strcat(ligs_d, '__', recs_d);

    if isfield(net, 'channel_capacity') && ~isempty(net.channel_capacity)
        cap_d = net.channel_capacity(:);
    elseif isfield(net, 'weights') && ~isempty(net.weights)
        cap_d = net.weights(:);
    else
        cap_d = ones(length(pairs_d), 1);
    end

    [~, ia, ib] = intersect(all_lr_pairs, pairs_d, 'stable');
    cross_dataset_matrix(ia, d) = cap_d(ib);
end

if verbose
    n_any = sum(any(cross_dataset_matrix > cap_thresh, 2));
    fprintf('  Active (capacity>%.3f) in at least one dataset: %d\n', cap_thresh, n_any);
    for d = 1:N
        n_active_d = sum(cross_dataset_matrix(:, d) > cap_thresh);
        fprintf('  %s: %d active channels\n', dataset_names{d}, n_active_d);
    end
end

% Step 2: Conservation scores.
if verbose, fprintf('\n--- Step 2: Conservation scoring ---\n'); end

is_active = cross_dataset_matrix > cap_thresh;
conservation_score = sum(is_active, 2) / N;
present_in_n = sum(is_active, 2);

is_conserved  = conservation_score == 1.0;
is_specialized = present_in_n == 1;

present_in_datasets_str = cell(n_lr_total, 1);
for i = 1:n_lr_total
    ds = dataset_names(is_active(i, :));
    present_in_datasets_str{i} = strjoin(ds, ';');
end

if verbose
    fprintf('  Conserved channels (all %d datasets): %d\n', N, sum(is_conserved));
    fprintf('  Specialized channels (1 dataset only): %d\n', sum(is_specialized));
    if N >= 2
        for d = 1:N
            n_ds_only = sum(is_active(:,d) & (present_in_n == 1));
            fprintf('  Dataset-specific to %s: %d channels\n', dataset_names{d}, n_ds_only);
        end
    end
end

% Step 3: NMF communication program matching via cosine similarity.
if verbose, fprintf('\n--- Step 3: NMF program matching ---\n'); end

all_programs      = {};
all_program_names = {};
all_program_datasets = [];

for d = 1:N
    net = networks{d};

    if ~isfield(net, 'programs') || isempty(net.programs)
        % No NMF programs: treat full capacity vector as a single surrogate program.
        ligs_d = net.ligands(:);
        recs_d = net.receptors(:);
        if ~iscell(ligs_d), ligs_d = cellstr(num2str(ligs_d)); end
        if ~iscell(recs_d), recs_d = cellstr(num2str(recs_d)); end

        if isfield(net, 'channel_capacity') && ~isempty(net.channel_capacity)
            cap_d = net.channel_capacity(:);
        elseif isfield(net, 'weights') && ~isempty(net.weights)
            cap_d = net.weights(:);
        else
            cap_d = ones(length(ligs_d), 1);
        end

        pairs_d = strcat(ligs_d, '__', recs_d);

        loading_vec = zeros(1, n_lr_total);
        [~, ia, ib] = intersect(all_lr_pairs, pairs_d, 'stable');
        loading_vec(ia) = cap_d(ib);

        if norm(loading_vec) > 0
            loading_vec = loading_vec / norm(loading_vec);
        end

        all_programs{end+1}        = loading_vec;
        all_program_names{end+1}   = sprintf('%s_P1', dataset_names{d});
        all_program_datasets(end+1) = d;
    else
        prog_struct = net.programs;
        if isfield(prog_struct, 'W')
            W = prog_struct.W;
            n_prog = size(W, 2);

            ligs_d = net.ligands(:);
            recs_d = net.receptors(:);
            if ~iscell(ligs_d), ligs_d = cellstr(num2str(ligs_d)); end
            if ~iscell(recs_d), recs_d = cellstr(num2str(recs_d)); end
            pairs_d = strcat(ligs_d, '__', recs_d);
            [~, ia, ib] = intersect(all_lr_pairs, pairs_d, 'stable');

            for prog_k = 1:n_prog
                loading_vec = zeros(1, n_lr_total);
                loading_vec(ia) = W(ib, prog_k);
                if norm(loading_vec) > 0
                    loading_vec = loading_vec / norm(loading_vec);
                end
                all_programs{end+1}        = loading_vec;
                all_program_names{end+1}   = sprintf('%s_P%d', dataset_names{d}, prog_k);
                all_program_datasets(end+1) = d;
            end
        end
    end
end

n_programs_total = length(all_programs);

if verbose
    fprintf('  Total programs across all datasets: %d\n', n_programs_total);
end

cosine_matrix = zeros(n_programs_total, n_programs_total);
for i = 1:n_programs_total
    for j = i:n_programs_total
        vi = all_programs{i};
        vj = all_programs{j};
        if norm(vi) > 0 && norm(vj) > 0
            c = dot(vi, vj);   % vectors already unit-normalized
        else
            c = 0;
        end
        cosine_matrix(i,j) = c;
        cosine_matrix(j,i) = c;
    end
end

conserved_program_rows = {};
for i = 1:n_programs_total
    for j = i+1:n_programs_total
        if all_program_datasets(i) ~= all_program_datasets(j)
            if cosine_matrix(i,j) >= cos_high
                conserved_program_rows{end+1} = {
                    all_program_names{i}, all_program_names{j}, ...
                    cosine_matrix(i,j), ...
                    dataset_names{all_program_datasets(i)}, ...
                    dataset_names{all_program_datasets(j)}
                };
            end
        end
    end
end

is_specialized_prog = false(n_programs_total, 1);
for i = 1:n_programs_total
    max_other_cosine = 0;
    for j = 1:n_programs_total
        if all_program_datasets(j) ~= all_program_datasets(i)
            max_other_cosine = max(max_other_cosine, cosine_matrix(i,j));
        end
    end
    is_specialized_prog(i) = (max_other_cosine < cos_low);
end

if verbose
    fprintf('  Conserved program pairs (cosine >= %.2f): %d\n', cos_high, length(conserved_program_rows));
    fprintf('  Specialized programs (max other-dataset cosine < %.2f): %d\n', ...
        cos_low, sum(is_specialized_prog));
end

if ~isempty(conserved_program_rows)
    cp_matrix = vertcat(conserved_program_rows{:});
    conserved_programs = table( ...
        cp_matrix(:,1), cp_matrix(:,2), cell2mat(cp_matrix(:,3)), ...
        cp_matrix(:,4), cp_matrix(:,5), ...
        'VariableNames', {'Program1', 'Program2', 'CosineSimilarity', 'Dataset1', 'Dataset2'});
else
    conserved_programs = table();
end

spec_prog_names    = all_program_names(is_specialized_prog)';
spec_prog_datasets = dataset_names(all_program_datasets(is_specialized_prog))';
if ~isempty(spec_prog_names)
    specialized_programs = table(spec_prog_names, spec_prog_datasets, ...
        'VariableNames', {'ProgramName', 'Dataset'});
else
    specialized_programs = table();
end

dominant_program_per_lr = cell(n_lr_total, 1);
for i = 1:n_lr_total
    best_prog = 1;
    best_load = -Inf;
    for p = 1:n_programs_total
        if all_programs{p}(i) > best_load
            best_load = all_programs{p}(i);
            best_prog = p;
        end
    end
    dominant_program_per_lr{i} = all_program_names{best_prog};
end

% Step 4: Shannon entropy of capacity distribution per dataset.
if verbose, fprintf('\n--- Step 4: Communication diversity index ---\n'); end

diversity_index = zeros(N, 1);
for d = 1:N
    cap_d = cross_dataset_matrix(:, d);
    total_cap = sum(cap_d);
    if total_cap > 0
        p = cap_d(cap_d > 0) / total_cap;
        diversity_index(d) = -sum(p .* log2(p));
    end
end

if verbose
    for d = 1:N
        fprintf('  %s: diversity index = %.3f bits\n', dataset_names{d}, diversity_index(d));
    end
end

% Step 5: Biological function annotation (keyword heuristic, not curated ontology).
if verbose, fprintf('\n--- Step 5: Biological function annotation ---\n'); end

bio_function = cell(n_lr_total, 1);

func_patterns = {
    {'BDNF','NGF','NT3','NT4','NTRK'},        'Neurotrophin signaling';
    {'CXCL','CCL','CXCR','CCR'},              'Chemokine signaling';
    {'WNT','FZD','LRP'},                       'Wnt signaling';
    {'EGF','EGFR','ERBB'},                     'EGF/ERBB signaling';
    {'IGF','IGFR','INSR'},                     'Insulin/IGF signaling';
    {'VEGF','FLT','KDR'},                      'VEGF/Angiogenesis';
    {'NOTCH','DLL','JAG'},                     'Notch signaling';
    {'TGF','TGFBR','BMP','SMAD'},             'TGFbeta/BMP signaling';
    {'FGF','FGFR'},                            'FGF signaling';
    {'PDGF','PDGFR'},                          'PDGF signaling';
    {'HGF','MET'},                             'HGF/MET signaling';
    {'SLIT','ROBO'},                           'SLIT/ROBO axon guidance';
    {'EPHRIN','EPH','EFNA','EFNB'},           'Ephrin signaling';
    {'SEMA','NRP','PLXN'},                     'Semaphorin signaling';
    {'NCAM','CDH','CNTN'},                     'Cell adhesion';
    {'GAP43','NRXN','NLGN'},                  'Synaptic organization';
    {'IL','TNFR','TNF','IFNG','IFNA'},        'Cytokine/Immune signaling';
    {'GABA','GABR'},                           'GABAergic signaling';
    {'GLUT','GRIN','AMPA','NMDA'},            'Glutamatergic signaling';
    {'DRD','HTR','ADRA','ADRB'},              'Monoamine signaling';
};

for i = 1:n_lr_total
    lr_str = upper([lr_ligands{i}, ' ', lr_receptors{i}]);
    assigned = false;
    for f = 1:length(func_patterns)
        keywords = func_patterns{f}{1};
        if any(cellfun(@(kw) contains(lr_str, kw), keywords))
            bio_function{i} = func_patterns{f}{2};
            assigned = true;
            break;
        end
    end
    if ~assigned
        bio_function{i} = 'Other/Unknown';
    end
end

[func_labels, ~, func_idx] = unique(bio_function);
func_counts = histcounts(func_idx, 1:length(func_labels)+1);
[func_counts_sorted, sort_ord] = sort(func_counts, 'descend');

if verbose
    fprintf('  Biological function distribution (top 10):\n');
    for k = 1:min(10, length(func_labels))
        fprintf('    %-30s %d pairs\n', func_labels{sort_ord(k)}, func_counts_sorted(k));
    end
end

% Step 6: Condition-specific annotation.
if verbose, fprintf('\n--- Step 6: Condition-specific annotation ---\n'); end

is_condition_specific = false(n_lr_total, 1);

for d = 1:N
    if d > length(condition_specifics) || isempty(condition_specifics{d})
        continue;
    end
    cs = condition_specifics{d};
    if ~istable(cs) || ~ismember('FDR', cs.Properties.VariableNames)
        continue;
    end
    if ~ismember('Ligand', cs.Properties.VariableNames) || ...
       ~ismember('Receptor', cs.Properties.VariableNames)
        continue;
    end

    sig_mask = cs.FDR < fdr_thresh;
    if ~any(sig_mask), continue; end

    sig_pairs = strcat(cs.Ligand(sig_mask), '__', cs.Receptor(sig_mask));
    [~, ia] = ismember(sig_pairs, all_lr_pairs);
    ia = ia(ia > 0);
    is_condition_specific(ia) = true;
end

if verbose
    fprintf('  Condition-specific pairs (FDR<%.2f in any dataset): %d\n', ...
        fdr_thresh, sum(is_condition_specific));
end

% Step 7: Assemble atlas table.
if verbose, fprintf('\n--- Step 7: Assembling atlas table ---\n'); end

mean_capacity = mean(cross_dataset_matrix, 2);
max_capacity  = max(cross_dataset_matrix, [], 2);

atlas_table = table( ...
    all_lr_pairs(:), ...
    lr_ligands(:), ...
    lr_receptors(:), ...
    conservation_score, ...
    present_in_n, ...
    present_in_datasets_str, ...
    dominant_program_per_lr, ...
    mean_capacity, ...
    max_capacity, ...
    is_condition_specific, ...
    bio_function, ...
    'VariableNames', { ...
        'LRPair', 'Ligand', 'Receptor', ...
        'ConservationScore', 'PresentInNDatasets', 'PresentInDatasets', ...
        'DominantProgram', 'MeanCapacity', 'MaxCapacity', ...
        'IsConditionSpecific', 'BiologicalFunction' ...
    });

[~, sort_order] = sortrows( ...
    [atlas_table.ConservationScore, atlas_table.MeanCapacity], [-1, -2]);
atlas_table = atlas_table(sort_order, :);

if verbose
    fprintf('  Atlas table: %d L-R pairs\n', height(atlas_table));
    fprintf('  Top 5 most conserved, high-capacity channels:\n');
    for k = 1:min(5, height(atlas_table))
        row = atlas_table(k,:);
        fprintf('    %s -> %s  conserved=%.2f  capacity=%.4f  [%s]\n', ...
            row.Ligand{1}, row.Receptor{1}, row.ConservationScore, ...
            row.MeanCapacity, row.BiologicalFunction{1});
    end
end

atlas = struct();
atlas.lr_table              = atlas_table;
atlas.conserved_programs    = conserved_programs;
atlas.specialized_programs  = specialized_programs;
atlas.diversity_index       = diversity_index;
atlas.n_datasets            = N;
atlas.dataset_names         = dataset_names;
atlas.cross_dataset_matrix  = cross_dataset_matrix;
atlas.all_lr_pairs          = all_lr_pairs;
atlas.cosine_matrix         = cosine_matrix;
atlas.program_names         = all_program_names(:);
atlas.conservation_score    = conservation_score;

if exist('results', 'dir')
    try
        writetable(atlas_table, fullfile('results', 'communication_atlas.csv'));
        if ~isempty(conserved_programs) && height(conserved_programs) > 0
            writetable(conserved_programs, fullfile('results', 'conserved_programs.csv'));
        end
        if ~isempty(specialized_programs) && height(specialized_programs) > 0
            writetable(specialized_programs, fullfile('results', 'specialized_programs.csv'));
        end
        div_table = table(dataset_names(:), diversity_index, ...
            'VariableNames', {'Dataset', 'DiversityIndex_bits'});
        writetable(div_table, fullfile('results', 'communication_diversity.csv'));
        if verbose, fprintf('  Atlas CSVs saved to results/\n'); end
    catch
    end
end

if verbose
    fprintf('\n=== build_communication_atlas summary ===\n');
    fprintf('  Datasets analyzed:         %d\n', N);
    fprintf('  Total L-R pairs in atlas:  %d\n', height(atlas_table));
    fprintf('  Fully conserved channels:  %d (score=1.0)\n', sum(conservation_score==1));
    fprintf('  Dataset-specific channels: %d (score=1/%d)\n', sum(present_in_n==1), N);
    fprintf('  Conserved program pairs:   %d (cosine>=%.2f)\n', ...
        height(conserved_programs), cos_high);
    fprintf('  Specialized programs:      %d (max cosine<%.2f)\n', ...
        height(specialized_programs), cos_low);
    fprintf('build_communication_atlas done.\n\n');
end

end
