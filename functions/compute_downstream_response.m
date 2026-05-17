function [response_matrix, gene_programs, response_scores] = compute_downstream_response( ...
        expr_matrix, gene_names, cell_metadata, condition_specific, varargin)
% COMPUTE_DOWNSTREAM_RESPONSE  Transcriptional response profiling in receiver cells.
%   Splits receivers into high/low receptor-expression groups, runs Wilcoxon
%   DE, takes PC1 of response genes as the response score for channel capacity.
%   Falls back to receptor expression directly when fewer than min_response_genes found.
%
% [response_matrix, gene_programs, response_scores] = ...
%     compute_downstream_response(expr_matrix, gene_names, cell_metadata, condition_specific,
%                                 'max_cells', 500, 'fdr_cutoff', 0.1)
%
% Inputs:
%   expr_matrix        - [cells × genes] sparse log-normalized expression
%   gene_names         - {genes×1} cell array of gene symbols
%   cell_metadata      - Table with 'Cluster' column
%   condition_specific - Table with Ligand, Receptor columns
%
% Outputs:
%   response_matrix - struct per receptor: .genes, .logFC, .fdr, .pc1_scores,
%                     .n_response_genes, .source ('downstream'|'receptor_fallback')
%   gene_programs   - struct of response programs shared across receptors
%   response_scores - struct of pairwise response score matrices

fprintf('\n=== Downstream Response Analysis ===\n\n');

p = inputParser();
addParameter(p, 'max_cells',          500);
addParameter(p, 'fdr_cutoff',         0.10);
addParameter(p, 'min_response_genes', 10);
addParameter(p, 'max_response_genes', 50);
parse(p, varargin{:});
opt = p.Results;

if istable(gene_names)
    if ismember('gene_name', gene_names.Properties.VariableNames)
        gene_names = gene_names.gene_name;
    else
        gene_names = gene_names{:,1};
    end
end
gene_names = cellstr(gene_names);

cluster_ids = cell_metadata.Cluster;
if iscategorical(cluster_ids), cluster_ids = double(cluster_ids); end

unique_receptors = unique(condition_specific.Receptor);
n_receptors      = length(unique_receptors);

response_matrix = struct();
response_scores = struct();
all_response_genes = {};

fprintf('Receptors to profile: %d\n\n', n_receptors);

for r_idx = 1:n_receptors
    receptor = unique_receptors{r_idx};
    rec_col  = gene_names_to_col(receptor, gene_names);

    if isnan(rec_col)
        fprintf('  [%d/%d] %s — gene not in matrix, skipping\n', r_idx, n_receptors, receptor);
        continue;
    end

    lr_rows     = strcmp(condition_specific.Receptor, receptor);
    source_clus = [];
    target_clus = [];
    if ismember('SourceCluster', condition_specific.Properties.VariableNames)
        source_clus = condition_specific.SourceCluster(lr_rows);
        target_clus = condition_specific.TargetCluster(lr_rows);
    end

    % Use all cells if no cluster restriction available
    if isempty(target_clus)
        receiver_mask = true(size(cluster_ids));
    else
        receiver_mask = ismember(cluster_ids, target_clus);
    end
    receiver_idx = find(receiver_mask);

    if length(receiver_idx) < 20
        fprintf('  [%d/%d] %s — too few receiver cells (%d), skipping\n', ...
            r_idx, n_receptors, receptor, length(receiver_idx));
        continue;
    end

    % Subsample
    if length(receiver_idx) > opt.max_cells
        receiver_idx = receiver_idx(randsample(length(receiver_idx), opt.max_cells));
    end

    rec_expr = full(expr_matrix(receiver_idx, rec_col));

    q33 = quantile(rec_expr, 1/3);
    q67 = quantile(rec_expr, 2/3);

    high_mask = rec_expr >= q67;
    low_mask  = rec_expr <= q33;

    n_high = sum(high_mask);
    n_low  = sum(low_mask);

    if n_high < 5 || n_low < 5
        fprintf('  [%d/%d] %s — insufficient variance in receiver expression, using fallback\n', ...
            r_idx, n_receptors, receptor);
        response_matrix.(matlab.lang.makeValidName(receptor)).source = 'receptor_fallback';
        response_matrix.(matlab.lang.makeValidName(receptor)).pc1_scores = rec_expr;
        continue;
    end

    high_idx = receiver_idx(high_mask);
    low_idx  = receiver_idx(low_mask);

    fprintf('  [%d/%d] %s — receivers: %d total, %d high-signal, %d low-signal\n', ...
        r_idx, n_receptors, receptor, length(receiver_idx), n_high, n_low);

    n_genes   = size(expr_matrix, 2);
    pvals     = ones(n_genes, 1);
    logFC     = zeros(n_genes, 1);

    % Batch-test genes expressed in at least 5% of receivers
    min_expr_frac = 0.05;
    expr_high = full(expr_matrix(high_idx, :));
    expr_low  = full(expr_matrix(low_idx,  :));

    pct_high = mean(expr_high > 0, 1)';
    pct_low  = mean(expr_low  > 0, 1)';
    test_mask = (pct_high > min_expr_frac | pct_low > min_expr_frac);
    test_idx  = find(test_mask);

    % Compute log2FC
    mean_high = mean(expr_high(:, test_mask), 1)';
    mean_low  = mean(expr_low(:,  test_mask), 1)';
    logFC(test_idx) = log2((mean_high + 0.1) ./ (mean_low + 0.1));

    % Wilcoxon (use fast vectorized ranksum via Matlab)
    for g = 1:length(test_idx)
        gi = test_idx(g);
        x  = expr_high(:, gi);
        y  = expr_low(:,  gi);
        if var(x) < 1e-12 && var(y) < 1e-12
            continue;
        end
        try
            pvals(gi) = ranksum(x, y);
        catch
        end
    end

    % BH correction
    [~, sort_order] = sort(pvals);
    n_tested = sum(test_mask);
    fdr      = zeros(n_genes, 1);
    ranks    = zeros(n_genes, 1);
    ranks(sort_order) = 1:n_genes;
    fdr = min(1, pvals .* n_tested ./ max(ranks, 1));
    % Monotone enforcement (Benjamini-Hochberg)
    for i = n_genes-1:-1:1
        fdr(sort_order(i)) = min(fdr(sort_order(i)), fdr(sort_order(i+1)));
    end

    % Select response genes
    sig_mask  = fdr < opt.fdr_cutoff & abs(logFC) > 0.5 & test_mask;
    % Remove receptor itself from response genes
    sig_mask(rec_col) = false;
    sig_idx   = find(sig_mask);

    % Sort by |logFC|, cap
    [~, sort_by_fc] = sort(abs(logFC(sig_idx)), 'descend');
    sig_idx = sig_idx(sort_by_fc(1:min(end, opt.max_response_genes)));

    n_response = length(sig_idx);
    fprintf('     Response genes found: %d  (FDR<%.2f, |log2FC|>0.5)\n', ...
        n_response, opt.fdr_cutoff);

    all_response_genes = [all_response_genes; gene_names(sig_idx)]; %#ok<AGROW>

    rfield = matlab.lang.makeValidName(receptor);

    if n_response >= opt.min_response_genes
            % PC1 of response gene matrix = response score
        all_receiver_idx = find(receiver_mask);
        if length(all_receiver_idx) > opt.max_cells
            all_receiver_idx = all_receiver_idx(randsample(length(all_receiver_idx), opt.max_cells));
        end
        resp_mat = full(expr_matrix(all_receiver_idx, sig_idx));  % n_recv × n_genes
        resp_mat = resp_mat - mean(resp_mat, 1);

        [scores, ~, ~] = svd(resp_mat, 'econ');
        pc1 = scores(:, 1);   % first left singular vector (U[:,1]) = PC1 scores
        % Orient so high-receptor cells have positive PC1
        rec_expr_all = full(expr_matrix(all_receiver_idx, rec_col));
        if corr(pc1, rec_expr_all, 'type', 'Spearman') < 0
            pc1 = -pc1;
        end

        response_matrix.(rfield).source            = 'downstream';
        response_matrix.(rfield).genes             = gene_names(sig_idx);
        response_matrix.(rfield).logFC             = logFC(sig_idx);
        response_matrix.(rfield).fdr               = fdr(sig_idx);
        response_matrix.(rfield).pc1_scores        = pc1;
        response_matrix.(rfield).n_response_genes  = n_response;
        response_matrix.(rfield).receiver_idx      = all_receiver_idx;
    else
        % Fallback: use receptor expression directly
        response_matrix.(rfield).source            = 'receptor_fallback';
        response_matrix.(rfield).genes             = {receptor};
        response_matrix.(rfield).logFC             = 0;
        response_matrix.(rfield).fdr               = 1;
        all_recv2 = find(receiver_mask);
        response_matrix.(rfield).pc1_scores        = full(expr_matrix(all_recv2, rec_col));
        response_matrix.(rfield).n_response_genes  = 0;
        response_matrix.(rfield).receiver_idx      = all_recv2;
    end
end

fprintf('\nIdentifying shared response gene programs across receptors...\n');
gene_programs = struct('genes', {{}}, 'receptors', {{}}, 'n_programs', 0);

unique_resp_genes = unique(all_response_genes);
if length(unique_resp_genes) > 10
    % Count how often each gene appears across receptors
    gene_freq = cellfun(@(g) sum(strcmp(all_response_genes, g)), unique_resp_genes);
    [sorted_freq, sf_idx] = sort(gene_freq, 'descend');
    top_genes = unique_resp_genes(sf_idx(sorted_freq >= 2));   % appear in ≥2 receptors
    n_programs = min(5, length(top_genes));
    gene_programs.genes      = top_genes(1:n_programs);
    gene_programs.n_programs = n_programs;
    fprintf('  Shared response genes (appear in ≥2 receptors): %d\n', length(top_genes));
end

fprintf('\nDownstream response analysis complete\n');
n_downstream = sum(structfun(@(s) isfield(s,'source') && strcmp(s.source,'downstream'), response_matrix));
n_fallback   = sum(structfun(@(s) isfield(s,'source') && strcmp(s.source,'receptor_fallback'), response_matrix));
fprintf('  Receptors with downstream response: %d\n', n_downstream);
fprintf('  Receptors using receptor fallback:  %d\n', n_fallback);

% Export response gene table
try
    out_rows = {};
    fnames = fieldnames(response_matrix);
    for fi = 1:length(fnames)
        rm = response_matrix.(fnames{fi});
        if isfield(rm, 'genes') && ~isempty(rm.genes)
            for gi = 1:length(rm.genes)
                fc_val = 0;
                if length(rm.logFC) >= gi, fc_val = rm.logFC(gi); end
                fdr_val = 1;
                if length(rm.fdr) >= gi, fdr_val = rm.fdr(gi); end
                out_rows{end+1} = {fnames{fi}, rm.genes{gi}, fc_val, fdr_val, rm.source}; %#ok<AGROW>
            end
        end
    end
    if ~isempty(out_rows)
        T = cell2table(vertcat(out_rows{:}), ...
            'VariableNames', {'Receptor','ResponseGene','log2FC','FDR','Source'});
        out_dir = 'results';
        if ~exist(out_dir,'dir'), out_dir = '.'; end
        writetable(T, fullfile(out_dir, 'downstream_response_genes.csv'));
        fprintf('  Saved: downstream_response_genes.csv (%d gene-receptor pairs)\n', height(T));
    end
catch ME
    fprintf('  Warning: could not save response genes: %s\n', ME.message);
end

end  % main


function col = gene_names_to_col(gene, gene_names)
idx = find(strcmpi(gene_names, gene), 1);
if isempty(idx)
    col = NaN;
else
    col = idx;
end
end
