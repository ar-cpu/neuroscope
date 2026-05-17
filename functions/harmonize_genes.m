function [genes_harmonized, expr_harmonized] = harmonize_genes(genes_table, expr_matrix, gene_mapping)
% HARMONIZE_GENES  Two-pass gene ID mapping (exact lookup then heuristic passthrough)
%   with automatic dimension-mismatch correction for off-by-one index artifacts.
%
% [genes_harmonized, expr_harmonized] = harmonize_genes(genes_table, expr_matrix, gene_mapping)
%
% Inputs:
%   genes_table  - Table with gene IDs in first column
%   expr_matrix  - [cells x genes] sparse or dense matrix
%   gene_mapping - containers.Map from ID to canonical symbol (create_gene_mapping)
%
% Outputs:
%   genes_harmonized - Table with 'gene_name' column (canonical symbols)
%   expr_harmonized  - Sparse matrix restricted to retained genes

    col_names = genes_table.Properties.VariableNames;
    gene_col = col_names{1};
    
    gene_list = genes_table.(gene_col);
    if ~iscell(gene_list)
        gene_list = cellstr(gene_list);
    end
    
    % Resolve off-by-one dimension mismatches introduced by index columns
    % or header rows present in some export formats.
    n_genes_csv = length(gene_list);
    n_genes_matrix = size(expr_matrix, 2);

    if n_genes_matrix == n_genes_csv + 1
        % Extra leading column in the matrix is a 0-based integer index artifact.
        fprintf('     Detected index column in matrix — dropping column 1\n');
        expr_matrix = expr_matrix(:, 2:end);
    elseif n_genes_matrix == n_genes_csv - 1
        % Extra leading entry in the gene list is a row-index placeholder (e.g., "0").
        fprintf('     Detected index entry in gene list — dropping row 1\n');
        gene_list = gene_list(2:end);
    elseif n_genes_matrix ~= n_genes_csv
        fprintf('     Dimension mismatch: matrix has %d columns but gene list has %d entries\n', ...
                n_genes_matrix, n_genes_csv);
        genes_harmonized = table(cell(0,1), 'VariableNames', {'gene_name'});
        expr_harmonized = sparse(size(expr_matrix, 1), 0);
        return;
    end
    
    keep_genes = false(length(gene_list), 1);
    mapped_genes = cell(length(gene_list), 1);
    
    % Pass 1: exact lookup in the reference gene map.
    for gene_idx = 1:length(gene_list)
        gene_str = char(gene_list{gene_idx});
        if ~isempty(gene_str) && isKey(gene_mapping, gene_str)
            mapped_genes{gene_idx} = gene_mapping(gene_str);
            keep_genes(gene_idx) = true;
        end
    end
    
    fprintf('     First pass: %d/%d genes\n', sum(keep_genes), length(gene_list));
    
    % Pass 2: heuristic acceptance of valid-looking identifiers not found in
    %         the reference map (ENSEMBL IDs, short alphanumeric symbols).
    for gene_idx = 1:length(gene_list)
        if ~keep_genes(gene_idx)
            gene_str = char(gene_list{gene_idx});
            if ~isempty(gene_str) && ~strcmp(gene_str, '0')
                if contains(gene_str, 'ENSG') || ...
                   (length(gene_str) <= 20 && length(gene_str) >= 2 && ...
                    all(isstrprop(gene_str, 'alphanum') | gene_str == '_' | gene_str == '-' | gene_str == '.') && ...
                    any(isstrprop(gene_str, 'alpha')))
                    mapped_genes{gene_idx} = gene_str;
                    keep_genes(gene_idx) = true;
                end
            end
        end
    end
    
    if sum(keep_genes) == 0
        fprintf('     No genes could be mapped to the reference namespace\n');
        genes_harmonized = table(cell(0,1), 'VariableNames', {'gene_name'});
        expr_harmonized = sparse(size(expr_matrix, 1), 0);
        return;
    end

    genes_harmonized = table(mapped_genes(keep_genes), 'VariableNames', {'gene_name'});
    expr_harmonized = expr_matrix(:, keep_genes);
    
    fprintf('     Final: %d/%d genes (%.1f%%)\n', sum(keep_genes), length(gene_list), ...
            100*sum(keep_genes)/length(gene_list));
end