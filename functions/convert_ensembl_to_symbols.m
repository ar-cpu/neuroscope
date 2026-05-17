function [genes_out, expr_out] = convert_ensembl_to_symbols(genes_table, expr_matrix, gene_table)
% CONVERT_ENSEMBL_TO_SYMBOLS  Map ENSEMBL IDs to HGNC gene symbols using
%   the CellPhoneDB gene reference table, then filter the expression matrix
%   to the successfully converted, deduplicated genes.
%
% Usage:
%   [genes_out, expr_out] = convert_ensembl_to_symbols(genes_table, expr_matrix, gene_table)
%
% Inputs:
%   genes_table  - Single-column table with variable 'gene_name' containing
%                  ENSEMBL IDs. A leading '0' index entry is skipped.
%   expr_matrix  - [cells × genes] sparse or dense expression matrix.
%   gene_table   - CellPhoneDB reference table with columns 'ensembl' and
%                  'gene_name'.
%
% Outputs:
%   genes_out - Table with variable 'gene_name' holding HGNC symbols.
%   expr_out  - Sparse expression matrix restricted to converted genes.

fprintf('  Converting ENSEMBL IDs to gene symbols...\n');

ensembl_ids = genes_table.gene_name;

ensembl_map = containers.Map();
for i = 1:height(gene_table)
    if ~isempty(gene_table.ensembl{i}) && ~isempty(gene_table.gene_name{i})
        ensembl_map(gene_table.ensembl{i}) = gene_table.gene_name{i};
    end
end

fprintf('    Created ENSEMBL mapping: %d entries\n', length(ensembl_map));

gene_symbols = cell(length(ensembl_ids), 1);
keep_idx = false(length(ensembl_ids), 1);

for i = 1:length(ensembl_ids)
    ensembl_id = ensembl_ids{i};
    if strcmp(ensembl_id, '0')
        continue;
    end
    if isKey(ensembl_map, ensembl_id)
        gene_symbols{i} = ensembl_map(ensembl_id);
        keep_idx(i) = true;
    end
end

% Deduplicate: keep first occurrence for each symbol
valid_symbols = gene_symbols(keep_idx);
[unique_symbols, unique_idx] = unique(valid_symbols, 'stable');

keep_positions = find(keep_idx);
final_idx = keep_positions(unique_idx);

expr_out = expr_matrix(:, final_idx);
genes_out = table(unique_symbols, 'VariableNames', {'gene_name'});

fprintf('    Converted: %d ENSEMBL IDs -> %d unique gene symbols\n', ...
    length(ensembl_ids), length(unique_symbols));

end
