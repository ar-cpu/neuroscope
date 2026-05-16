function [genes_out, expr_out] = convert_ensembl_to_symbols(genes_table, expr_matrix, gene_table)
% CONVERT_ENSEMBL_TO_SYMBOLS  Map ENSEMBL gene identifiers (ENSG…) to
%                              HGNC gene symbols using the CellPhoneDB gene
%                              reference table, then filter the expression
%                              matrix to the successfully converted genes.
%
% Several publicly available single-cell datasets (e.g., Allen Brain Cell
% Atlas, GTEx, Human Cell Atlas) distribute count matrices indexed by ENSEMBL
% gene IDs rather than gene symbols. This function resolves those IDs against
% the CellPhoneDB gene_table to produce a symbol-indexed matrix compatible
% with the rest of the NEUROSCOPE pipeline, which operates on HGNC symbols
% for ligand-receptor interaction lookup.
%
% Duplicate symbols arising from one-to-many ENSEMBL→symbol mappings are
% resolved by retaining the first occurrence (stable deduplication).
%
% Usage:
%   [genes_out, expr_out] = convert_ensembl_to_symbols(genes_table, expr_matrix, gene_table)
%
% Inputs:
%   genes_table  - Single-column table with variable name 'gene_name'
%                  containing ENSEMBL IDs (e.g., 'ENSG00000141510').
%                  A leading '0' index entry is skipped automatically.
%   expr_matrix  - [cells × genes] sparse or dense expression matrix.
%                  Column count must match height(genes_table).
%   gene_table   - CellPhoneDB gene reference table with (at minimum) columns
%                  'ensembl' and 'gene_name'. Loaded from gene_table.csv.
%
% Outputs:
%   genes_out - Single-column table with variable name 'gene_name' containing
%               the HGNC symbol for each successfully converted gene.
%   expr_out  - Sparse expression matrix with columns restricted to converted,
%               deduplicated genes, in the same cell order as expr_matrix.

fprintf('  Converting ENSEMBL IDs to gene symbols...\n');

% Extract ENSEMBL IDs
ensembl_ids = genes_table.gene_name;

% Create ENSEMBL to gene symbol mapping from CellPhoneDB
ensembl_map = containers.Map();

for i = 1:height(gene_table)
    if ~isempty(gene_table.ensembl{i}) && ~isempty(gene_table.gene_name{i})
        ensembl_id = gene_table.ensembl{i};
        gene_symbol = gene_table.gene_name{i};
        
        % Store mapping
        ensembl_map(ensembl_id) = gene_symbol;
    end
end

fprintf('    Created ENSEMBL mapping: %d entries\n', length(ensembl_map));

% Convert IDs to symbols
gene_symbols = cell(length(ensembl_ids), 1);
keep_idx = false(length(ensembl_ids), 1);

for i = 1:length(ensembl_ids)
    ensembl_id = ensembl_ids{i};
    
    % Skip the index row if present
    if strcmp(ensembl_id, '0')
        continue;
    end
    
    % Try to find in mapping
    if isKey(ensembl_map, ensembl_id)
        gene_symbols{i} = ensembl_map(ensembl_id);
        keep_idx(i) = true;
    end
end

% Remove duplicates (keep first occurrence)
valid_symbols = gene_symbols(keep_idx);
[unique_symbols, unique_idx] = unique(valid_symbols, 'stable');

% Get corresponding indices in original array
keep_positions = find(keep_idx);
final_idx = keep_positions(unique_idx);

% Filter expression matrix and create output
expr_out = expr_matrix(:, final_idx);
genes_out = table(unique_symbols, 'VariableNames', {'gene_name'});

fprintf('    Converted: %d ENSEMBL IDs -> %d unique gene symbols\n', ...
    length(ensembl_ids), length(unique_symbols));

end