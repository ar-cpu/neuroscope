function [norm_matrix, hvg_names, hvg_indices, size_factors] = preprocess_expression(expr_matrix, gene_names, varargin)
% PREPROCESS_EXPRESSION  Size-factor normalization and HVG selection for scRNA-seq data.
%
% [norm_matrix, hvg_names, hvg_indices, size_factors] = ...
%     preprocess_expression(expr_matrix, gene_names, 'n_hvg', 3000, 'min_cells', 0.01)
%
% Inputs:
%   expr_matrix  - [n_cells x n_genes] raw count matrix (sparse or dense)
%   gene_names   - cell array, table (with 'gene_name' col), or struct
%
% Outputs:
%   norm_matrix  - [n_cells x n_hvg] log1p-normalized dense double
%   hvg_names    - table with 'gene_name' column for selected HVGs
%   hvg_indices  - [n_hvg x 1] indices into original gene_names
%   size_factors - [n_cells x 1] per-cell normalization factors (median-ratio)

p = inputParser();
p.addRequired('expr_matrix');
p.addRequired('gene_names');
p.addParameter('n_hvg',     3000,  @(x) isnumeric(x) && x > 0);
p.addParameter('min_cells', 0.01,  @(x) isnumeric(x) && x >= 0 && x <= 1);
p.addParameter('verbose',   true,  @islogical);
p.parse(expr_matrix, gene_names, varargin{:});

n_hvg_target  = p.Results.n_hvg;
min_cell_frac = p.Results.min_cells;
verbose       = p.Results.verbose;

if istable(gene_names)
    if ismember('gene_name', gene_names.Properties.VariableNames)
        gene_names_cell = gene_names.gene_name;
    else
        gene_names_cell = gene_names{:,1};
    end
    if ~iscell(gene_names_cell)
        gene_names_cell = cellstr(num2str(gene_names_cell));
    end
elseif iscell(gene_names)
    gene_names_cell = gene_names;
elseif isstruct(gene_names) && isfield(gene_names, 'gene_name')
    gene_names_cell = gene_names.gene_name;
    if ~iscell(gene_names_cell)
        gene_names_cell = cellstr(num2str(gene_names_cell));
    end
else
    error('preprocess_expression: gene_names must be a table, cell array, or struct with gene_name field.');
end

gene_names_cell = gene_names_cell(:);

[n_cells, n_genes_orig] = size(expr_matrix);

if length(gene_names_cell) ~= n_genes_orig
    error('preprocess_expression: gene_names length (%d) does not match expr_matrix columns (%d).', ...
        length(gene_names_cell), n_genes_orig);
end

if verbose
    fprintf('preprocess_expression: %d cells x %d genes (input)\n', n_cells, n_genes_orig);
end

% Step 1: Size factor normalization
if verbose
    fprintf('\n--- Step 1: Size factor normalization ---\n');
end

t_norm = tic;

if issparse(expr_matrix)
    lib_sizes = full(sum(expr_matrix, 2));
else
    lib_sizes = sum(expr_matrix, 2);
end

zero_lib = lib_sizes == 0;
n_zero = sum(zero_lib);
if n_zero > 0
    if verbose
        fprintf('  Warning: %d cells have zero total counts. Setting size factor to 1.\n', n_zero);
    end
    lib_sizes(zero_lib) = 1;
end

lib_median = median(lib_sizes(~zero_lib));
lib_mean   = mean(lib_sizes(~zero_lib));
lib_q25    = quantile(lib_sizes(~zero_lib), 0.25);
lib_q75    = quantile(lib_sizes(~zero_lib), 0.75);

if verbose
    fprintf('  Library size distribution (non-zero cells):\n');
    fprintf('    Median: %.0f counts\n', lib_median);
    fprintf('    Mean:   %.0f counts\n', lib_mean);
    fprintf('    IQR:    %.0f - %.0f counts\n', lib_q25, lib_q75);
    fprintf('    Min:    %.0f  Max: %.0f counts\n', min(lib_sizes(~zero_lib)), max(lib_sizes(~zero_lib)));
end

size_factors = lib_sizes / lib_median;

BATCH_SIZE = 5000;

if verbose
    fprintf('  Applying size factor normalization (batch size %d)...\n', BATCH_SIZE);
end

norm_matrix_full = zeros(n_cells, n_genes_orig, 'single');

n_batches = ceil(n_cells / BATCH_SIZE);
for b = 1:n_batches
    idx_start = (b-1)*BATCH_SIZE + 1;
    idx_end   = min(b*BATCH_SIZE, n_cells);
    batch_idx = idx_start:idx_end;

    if issparse(expr_matrix)
        batch_counts = full(expr_matrix(batch_idx, :));
    else
        batch_counts = expr_matrix(batch_idx, :);
    end
    batch_counts = max(0, real(double(batch_counts)));  % clip signed-int / complex storage artifacts

    sf_batch   = size_factors(batch_idx);
    sf_batch   = max(sf_batch, 1e-6);  % guard against zero/negative size factors
    batch_norm = bsxfun(@rdivide, batch_counts, sf_batch);
    norm_matrix_full(batch_idx, :) = single(log1p(batch_norm));
end

if verbose
    fprintf('  Normalization complete (%.1f s)\n', toc(t_norm));
    n_zero_entries = sum(sum(norm_matrix_full == 0));
    sparsity = n_zero_entries / numel(norm_matrix_full) * 100;
    fprintf('  Normalized matrix sparsity: %.1f%%\n', sparsity);
end

% Step 2: Highly variable gene (HVG) selection via log-log polynomial mean-variance trend
if verbose
    fprintf('\n--- Step 2: Highly variable gene selection ---\n');
end

t_hvg = tic;

n_cells_expressing = sum(norm_matrix_full > 0, 1);
min_cells_abs       = ceil(min_cell_frac * n_cells);
passes_filter       = n_cells_expressing >= min_cells_abs;
n_pass              = sum(passes_filter);

if verbose
    fprintf('  Gene expression filter (min %.1f%% cells = %d cells):\n', ...
        min_cell_frac * 100, min_cells_abs);
    fprintf('    Genes passing filter: %d / %d (%.1f%%)\n', ...
        n_pass, n_genes_orig, 100 * n_pass / n_genes_orig);
end

if n_pass < n_hvg_target
    if verbose
        fprintf('  Warning: fewer genes pass filter (%d) than n_hvg (%d). Using all passing genes.\n', ...
            n_pass, n_hvg_target);
    end
    n_hvg_actual = n_pass;
else
    n_hvg_actual = n_hvg_target;
end

filtered_cols = find(passes_filter);
norm_filtered = double(norm_matrix_full(:, filtered_cols));

gene_means = mean(norm_filtered, 1);
gene_vars  = var(norm_filtered, 0, 1);

valid_mv = (gene_means > 0) & (gene_vars > 0);
n_valid  = sum(valid_mv);

if verbose
    fprintf('  Genes with mean>0 and var>0: %d / %d\n', n_valid, n_pass);
end

log_mean_valid = log(gene_means(valid_mv));
log_var_valid  = log(gene_vars(valid_mv));

X_design  = [ones(sum(valid_mv), 1), log_mean_valid(:), log_mean_valid(:).^2];
[Q, R]    = qr(X_design, 0);
poly_coef = R \ (Q' * log_var_valid(:));

log_mean_all   = log(max(gene_means, 1e-10));
fitted_log_var = poly_coef(1) + poly_coef(2)*log_mean_all + poly_coef(3)*log_mean_all.^2;

residual_var = zeros(1, n_pass);
for g = 1:n_pass
    if gene_vars(g) > 0
        residual_var(g) = gene_vars(g) / exp(fitted_log_var(g));
    else
        residual_var(g) = 0;
    end
end

residual_var(~isfinite(residual_var)) = 0;

[~, sort_order]      = sort(residual_var, 'descend');
selected_in_filtered = sort_order(1:n_hvg_actual);

hvg_indices = filtered_cols(selected_in_filtered)';
hvg_names   = gene_names_cell(hvg_indices);

norm_matrix = double(norm_matrix_full(:, hvg_indices));

clear norm_matrix_full norm_filtered;

if verbose
    fprintf('  Mean-variance trend polynomial coefficients: [%.3f, %.3f, %.3f]\n', ...
        poly_coef(1), poly_coef(2), poly_coef(3));
    fprintf('  Residual variance range: %.2f - %.2f\n', ...
        min(residual_var(residual_var>0)), max(residual_var));

    top10_residuals = residual_var(selected_in_filtered(1:min(10,n_hvg_actual)));
    top10_names     = hvg_names(1:min(10,n_hvg_actual));
    fprintf('  Top 10 HVGs by residual variance:\n');
    for k = 1:length(top10_names)
        fprintf('    %2d. %-12s  residual_var = %.3f\n', k, top10_names{k}, top10_residuals(k));
    end

    fprintf('  HVG selection complete (%.1f s)\n', toc(t_hvg));
end

% Step 3: Return filtered matrix
if verbose
    fprintf('\n--- Step 3: Final matrix ---\n');
    fprintf('  Output: %d cells x %d HVGs\n', size(norm_matrix, 1), size(norm_matrix, 2));
    fprintf('  Size factors: median=%.3f, IQR=[%.3f, %.3f]\n', ...
        median(size_factors), quantile(size_factors, 0.25), quantile(size_factors, 0.75));
    fprintf('  norm_matrix: %.1f MB (dense double)\n', ...
        numel(norm_matrix) * 8 / 1e6);
    fprintf('\npreprocess_expression done.\n\n');
end

if exist('results', 'dir')
    hvg_table = table(hvg_names, hvg_indices, residual_var(selected_in_filtered)', ...
        gene_means(selected_in_filtered)', gene_vars(selected_in_filtered)', ...
        'VariableNames', {'gene_name', 'original_index', 'residual_variance', 'mean_expr', 'variance'});
    out_path = fullfile('results', 'hvg_selection.csv');
    try
        writetable(hvg_table, out_path);
        if verbose
            fprintf('  HVG table saved: %s\n', out_path);
        end
    catch
    end

    sf_table = table((1:n_cells)', size_factors, lib_sizes, ...
        'VariableNames', {'cell_index', 'size_factor', 'library_size'});
    sf_path = fullfile('results', 'size_factors.csv');
    try
        writetable(sf_table, sf_path);
    catch
    end
end

hvg_names_table = table(hvg_names, 'VariableNames', {'gene_name'});
hvg_names = hvg_names_table;

end
