function gene_mapping = create_gene_mapping(gene_synonym_table, gene_table)
% CREATE_GENE_MAPPING  Build a case-normalized gene identifier lookup table
%                      that maps synonyms, aliases, and ENSEMBL IDs to
%                      canonical HGNC gene symbols.
%
% The resulting containers.Map is consumed by harmonize_genes() to unify gene
% identifiers across datasets that may use different naming conventions
% (e.g., one dataset uses ENSG00000141510 while another uses TP53).
%
% Mapping priority (earlier entries are not overwritten by later ones):
%   1. Self-mapping of canonical HGNC symbols from gene_table.
%   2. Synonym → canonical mappings from gene_synonym_table.
%   3. ENSEMBL ID → canonical gene symbol from gene_table.ensembl_gene_id.
%   4. Hardcoded aliases for well-known immune checkpoint and signaling genes.
%
% All keys are stored in upper case to enable case-insensitive lookups.
%
% Usage:
%   gene_mapping = create_gene_mapping(gene_synonym_table, gene_table)
%
% Inputs:
%   gene_synonym_table - Table of synonym pairs. The function infers column
%                        roles by name (looks for 'synonym'/'alias' and
%                        'gene_name'/'official'/'symbol'). Falls back to
%                        columns 1 and 2 if names are not recognized.
%   gene_table         - Reference gene table. Must contain at least one
%                        column named 'gene_name' or 'hgnc_symbol'. Optionally
%                        contains 'ensembl_gene_id' for ENSG→symbol mapping.
%
% Outputs:
%   gene_mapping - containers.Map with KeyType='char', ValueType='char'.
%                  Keys are upper-cased gene identifiers; values are canonical
%                  mixed-case gene symbols.

fprintf('Building gene mapping dictionary...\n');

% Initialize map
gene_mapping = containers.Map('KeyType', 'char', 'ValueType', 'char');

%% Phase 1: Self-map canonical HGNC symbols
fprintf('  Processing canonical genes...\n');

if any(strcmp(gene_table.Properties.VariableNames, 'gene_name'))
    canonical_genes = gene_table.gene_name;
elseif any(strcmp(gene_table.Properties.VariableNames, 'hgnc_symbol'))
    canonical_genes = gene_table.hgnc_symbol;
else
    % Use first column
    canonical_genes = gene_table{:, 1};
end

% Ensure cell array
if ~iscell(canonical_genes)
    canonical_genes = cellstr(canonical_genes);
end

n_canonical = 0;
for i = 1:length(canonical_genes)
    gene = char(canonical_genes{i});
    
    if isempty(gene) || strcmp(gene, '0') || strcmp(gene, '')
        continue;
    end
    
    % Self-map canonical genes
    gene_upper = upper(strtrim(gene));
    if ~isKey(gene_mapping, gene_upper)
        gene_mapping(gene_upper) = gene;
        n_canonical = n_canonical + 1;
    end
end

fprintf('    Mapped %d canonical genes\n', n_canonical);

%% Phase 2: Synonym → canonical mappings from synonym table
fprintf('  Processing gene synonyms...\n');

n_synonyms = 0;

if height(gene_synonym_table) > 0
    % Get column names
    synonym_cols = gene_synonym_table.Properties.VariableNames;
    
    % Find synonym and target columns
    synonym_col = '';
    target_col = '';
    
    for col_idx = 1:length(synonym_cols)
        col_name = lower(synonym_cols{col_idx});
        if contains(col_name, 'synonym') || contains(col_name, 'alias')
            synonym_col = synonym_cols{col_idx};
        end
        if contains(col_name, 'gene_name') || contains(col_name, 'official') || contains(col_name, 'symbol')
            target_col = synonym_cols{col_idx};
        end
    end
    
    % If not found, use first two columns
    if isempty(synonym_col)
        synonym_col = synonym_cols{1};
    end
    if isempty(target_col)
        if length(synonym_cols) >= 2
            target_col = synonym_cols{2};
        else
            target_col = synonym_col; % Self-map
        end
    end
    
    % Build synonym mapping
    for i = 1:height(gene_synonym_table)
        synonym = char(gene_synonym_table{i, synonym_col});
        target = char(gene_synonym_table{i, target_col});
        
        if isempty(synonym) || isempty(target) || strcmp(synonym, '0') || strcmp(target, '0')
            continue;
        end
        
        synonym_upper = upper(strtrim(synonym));
        target_clean = strtrim(target);
        
        % Only add if not already mapped (canonical has priority)
        if ~isKey(gene_mapping, synonym_upper)
            gene_mapping(synonym_upper) = target_clean;
            n_synonyms = n_synonyms + 1;
        end
    end
end

fprintf('    Mapped %d synonyms\n', n_synonyms);

%% Phase 3: ENSEMBL ID → gene symbol mappings
fprintf('  Processing Ensembl IDs...\n');

n_ensembl = 0;

if any(strcmp(gene_table.Properties.VariableNames, 'ensembl_gene_id')) && ...
   any(strcmp(gene_table.Properties.VariableNames, 'gene_name'))
    
    ensembl_ids = gene_table.ensembl_gene_id;
    gene_names = gene_table.gene_name;
    
    if ~iscell(ensembl_ids)
        ensembl_ids = cellstr(ensembl_ids);
    end
    if ~iscell(gene_names)
        gene_names = cellstr(gene_names);
    end
    
    for i = 1:length(ensembl_ids)
        ensembl = char(ensembl_ids{i});
        name = char(gene_names{i});
        
        if isempty(ensembl) || isempty(name) || strcmp(ensembl, '0') || strcmp(name, '0')
            continue;
        end
        
        if startsWith(ensembl, 'ENSG')
            ensembl_upper = upper(strtrim(ensembl));
            name_clean = strtrim(name);
            
            if ~isKey(gene_mapping, ensembl_upper)
                gene_mapping(ensembl_upper) = name_clean;
                n_ensembl = n_ensembl + 1;
            end
        end
    end
end

fprintf('    Mapped %d Ensembl IDs\n', n_ensembl);

%% Phase 4: Hardcoded aliases for immune checkpoint and signaling genes
fprintf('  Adding common gene aliases...\n');

common_aliases = {
    'CD274', 'PD-L1';
    'PDCD1', 'PD-1';
    'CTLA4', 'CD152';
    'HAVCR2', 'TIM3';
    'LAG3', 'CD223';
    'TNFRSF18', 'GITR';
    'TNFRSF9', '4-1BB';
    'TNFRSF4', 'OX40';
    'TIGIT', 'VSTM3';
    'CD47', 'IAP';
};

n_aliases = 0;
for i = 1:size(common_aliases, 1)
    alias = upper(common_aliases{i, 1});
    target = common_aliases{i, 2};
    
    if ~isKey(gene_mapping, alias)
        gene_mapping(alias) = target;
        n_aliases = n_aliases + 1;
    end
end

fprintf('    Added %d common aliases\n', n_aliases);

%% Phase 5: Validation summary
fprintf('  Validating mapping...\n');

total_mappings = length(gene_mapping);
unique_targets = length(unique(values(gene_mapping)));

fprintf('    Total mappings: %d\n', total_mappings);
fprintf('    Unique targets: %d\n', unique_targets);
fprintf('    Compression ratio: %.2f:1\n', total_mappings / unique_targets);

%% FINAL SUMMARY
fprintf('  Gene mapping complete\n');
fprintf('    Canonical: %d\n', n_canonical);
fprintf('    Synonyms: %d\n', n_synonyms);
fprintf('    Ensembl: %d\n', n_ensembl);
fprintf('    Aliases: %d\n', n_aliases);
fprintf('    Total: %d unique mappings\n\n', total_mappings);

end