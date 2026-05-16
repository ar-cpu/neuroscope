function lr_pairs = get_cellphonedb_interactions(cellphonedb, gene_names)
% GET_CELLPHONEDB_INTERACTIONS  Extract validated ligand-receptor pairs from
%                                a parsed CellPhoneDB v2/v3 database, filtered
%                                to genes present in the expression dataset.
%
% CellPhoneDB stores interactions via an indirect mapping through a multidata
% identifier table. The resolution chain is:
%
%   interactions.multidata_1_id  ──►  multidata_table.id
%       └──► protein_table.protein_multidata_id  ──►  gene symbol
%
% The multidata table must be pre-exported to converted/multidata_table.csv
% by the CellPhoneDB Python package before running this function.
%
% If CellPhoneDB extraction yields fewer than 50 pairs (typically because the
% dataset covers a narrow gene panel), the result is supplemented with a
% manually curated set of high-confidence L-R pairs spanning canonical signaling
% pathways (VEGF, EGFR, PDGF, FGF, TGF-β, Notch, Wnt, BMP, chemokines,
% immune checkpoint, Ephrin, and Semaphorin families).
%
% Usage:
%   lr_pairs = get_cellphonedb_interactions(cellphonedb, gene_names)
%
% Inputs:
%   cellphonedb - Struct with fields:
%                   .interactions  - Table loaded from interactions.csv
%                   .protein_table - Table loaded from protein_table.csv
%   gene_names  - {G×1} cell array of gene symbols present in the expression
%                 matrix. Only pairs where both partners appear in this list
%                 are retained.
%
% Outputs:
%   lr_pairs - {N×2} cell array where column 1 contains ligand symbols and
%              column 2 contains receptor symbols. Pairs are deduplicated.
%              Raises an error if no pairs are available after all fallbacks.
%
% References:
%   Efremova et al. (2020) CellPhoneDB: inferring cell-cell communication
%   from combined expression of multi-subunit ligand-receptor complexes.
%   Nature Protocols 15:1484-1506.

fprintf('Extracting L-R pairs from CellPhoneDB...\n');

lr_pairs = {};

try
    %% Load multidata table

    if exist('converted/multidata_table.csv', 'file')
        multidata = readtable('converted/multidata_table.csv');
        fprintf('  Multidata table: %d entries\n', height(multidata));
    else
        error('Multidata table not found: converted/multidata_table.csv');
    end

    interactions = cellphonedb.interactions;
    protein_table = cellphonedb.protein_table;

    %% Load gene_table for protein_name -> hgnc_symbol mapping
    % This is the key fix: CellPhoneDB protein names (e.g. ST1A1_HUMAN) do NOT
    % always match HGNC symbols (e.g. SULT1A1). gene_table.csv has the correct
    % protein_name -> hgnc_symbol mapping.
    pn_to_hgnc = containers.Map('KeyType', 'char', 'ValueType', 'char');
    if exist('converted/gene_table.csv', 'file')
        gt = readtable('converted/gene_table.csv');
        for i = 1:height(gt)
            if ismember('protein_name', gt.Properties.VariableNames) && ...
               ismember('hgnc_symbol', gt.Properties.VariableNames)
                pn = upper(strtrim(gt.protein_name{i}));
                hs = upper(strtrim(gt.hgnc_symbol{i}));
                if ~isempty(pn) && ~isempty(hs)
                    pn_to_hgnc(pn) = hs;
                end
            end
        end
        fprintf('  gene_table protein->hgnc mappings: %d\n', length(pn_to_hgnc));
    end

    %% Load complex_composition_table for complex -> subunit resolution
    complex_to_proteins = containers.Map('KeyType', 'double', 'ValueType', 'any');
    if exist('converted/complex_composition_table.csv', 'file')
        comp_comp = readtable('converted/complex_composition_table.csv');
        for i = 1:height(comp_comp)
            cmid = comp_comp.complex_multidata_id(i);
            pmid = comp_comp.protein_multidata_id(i);
            if isKey(complex_to_proteins, cmid)
                complex_to_proteins(cmid) = [complex_to_proteins(cmid), pmid];
            else
                complex_to_proteins(cmid) = pmid;
            end
        end
        fprintf('  Complex compositions loaded: %d complexes\n', length(complex_to_proteins));
    end

    %% Build mapping: multidata_id -> gene symbol

    fprintf('  Building protein-gene mapping...\n');

    multidata_to_gene = containers.Map('KeyType', 'double', 'ValueType', 'any');

    % Map individual proteins using gene_table (protein_name -> hgnc_symbol)
    prot_mid_to_gene = containers.Map('KeyType', 'double', 'ValueType', 'char');
    for i = 1:height(protein_table)
        multidata_id = protein_table.protein_multidata_id(i);
        protein_name = upper(strtrim(protein_table.protein_name{i}));
        if isKey(pn_to_hgnc, protein_name)
            gene_symbol = pn_to_hgnc(protein_name);
        else
            gene_symbol = extract_gene_symbol(protein_table.protein_name{i});
        end
        multidata_to_gene(multidata_id) = gene_symbol;
        prot_mid_to_gene(multidata_id) = gene_symbol;
    end

    % Map complex multidata_ids to their best representative subunit gene
    % "best" = first subunit whose gene symbol is present in the expression matrix
    comp_keys = keys(complex_to_proteins);
    for ci = 1:length(comp_keys)
        cmid = comp_keys{ci};
        pmids = complex_to_proteins(cmid);
        chosen = '';
        for pi = 1:length(pmids)
            if isKey(prot_mid_to_gene, pmids(pi))
                g = prot_mid_to_gene(pmids(pi));
                if any(strcmp(gene_names, g))
                    chosen = g;
                    break;
                elseif isempty(chosen)
                    chosen = g;  % fallback if nothing in matrix
                end
            end
        end
        if ~isempty(chosen)
            multidata_to_gene(cmid) = chosen;
        end
    end

    fprintf('  Mapped %d multidata entries to gene symbols (proteins + complexes)\n', length(multidata_to_gene));

    %% Extract interaction pairs

    fprintf('  Processing %d interactions...\n', height(interactions));

    found_count = 0;
    skipped_count = 0;

    for int_idx = 1:height(interactions)
        try
            multidata_1_id = interactions.multidata_1_id(int_idx);
            multidata_2_id = interactions.multidata_2_id(int_idx);

            if isKey(multidata_to_gene, multidata_1_id) && isKey(multidata_to_gene, multidata_2_id)
                gene_1 = multidata_to_gene(multidata_1_id);
                gene_2 = multidata_to_gene(multidata_2_id);

                if any(strcmp(gene_names, gene_1)) && any(strcmp(gene_names, gene_2))
                    lr_pairs{end+1, 1} = gene_1;
                    lr_pairs{end, 2} = gene_2;
                    found_count = found_count + 1;
                else
                    skipped_count = skipped_count + 1;
                end
            else
                skipped_count = skipped_count + 1;
            end
        catch
            skipped_count = skipped_count + 1;
            continue;
        end

        if mod(int_idx, 500) == 0
            fprintf('  Processed %d/%d (found: %d)\n', int_idx, height(interactions), found_count);
        end
    end

    fprintf('  Found: %d valid pairs | Skipped: %d\n', found_count, skipped_count);

    %% Deduplicate

    if ~isempty(lr_pairs)
        n_before = size(lr_pairs, 1);
        [~, unique_idx] = unique(strcat(lr_pairs(:,1), '_', lr_pairs(:,2)), 'stable');
        lr_pairs = lr_pairs(unique_idx, :);
        n_removed = n_before - size(lr_pairs, 1);
        if n_removed > 0
            fprintf('  Removed %d duplicates\n', n_removed);
        end
        fprintf('  Unique pairs: %d\n', size(lr_pairs, 1));
    end

    %% Supplement with curated pairs if needed

    min_pairs_threshold = 50;

    if isempty(lr_pairs) || size(lr_pairs, 1) < min_pairs_threshold
        fprintf('  Only %d pairs from CellPhoneDB - supplementing with curated list\n', ...
            size(lr_pairs, 1));

        known_pairs = get_known_ligand_receptor_pairs();

        valid_known = false(size(known_pairs, 1), 1);
        for i = 1:size(known_pairs, 1)
            if any(strcmp(gene_names, known_pairs{i, 1})) && ...
               any(strcmp(gene_names, known_pairs{i, 2}))
                valid_known(i) = true;
            end
        end

        known_pairs_valid = known_pairs(valid_known, :);
        fprintf('  Adding %d curated pairs\n', sum(valid_known));

        if isempty(lr_pairs)
            lr_pairs = known_pairs_valid;
        else
            lr_pairs = [lr_pairs; known_pairs_valid];
        end

        if ~isempty(lr_pairs)
            [~, unique_idx] = unique(strcat(lr_pairs(:,1), '_', lr_pairs(:,2)), 'stable');
            lr_pairs = lr_pairs(unique_idx, :);
        end

        fprintf('  Combined total: %d unique pairs\n', size(lr_pairs, 1));
    end

catch ME
    fprintf('  Error parsing CellPhoneDB: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('  at %s line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
    fprintf('  Falling back to curated list...\n');
    lr_pairs = get_known_ligand_receptor_pairs();
    fprintf('  Using %d curated pairs\n', size(lr_pairs, 1));
end

if isempty(lr_pairs)
    error('No L-R pairs available. Check CellPhoneDB data files.');
end

fprintf('\n  L-R pairs extracted: %d\n\n', size(lr_pairs, 1));

end

%% Local helper functions

function gene_symbol = extract_gene_symbol(protein_name)
% Extract gene symbol from CellPhoneDB protein names (e.g., ESR1_HUMAN -> ESR1)

    gene_symbol = protein_name;

    gene_symbol = strrep(gene_symbol, '_HUMAN', '');
    gene_symbol = strrep(gene_symbol, '_MOUSE', '');

    % Take first part for complex entries (e.g., GENE1:GENE2)
    if contains(gene_symbol, ':')
        parts = strsplit(gene_symbol, ':');
        gene_symbol = parts{1};
    end

    gene_symbol = regexprep(gene_symbol, '[.-]\d+$', '');
    gene_symbol = upper(strtrim(gene_symbol));
end

function pairs = get_known_ligand_receptor_pairs()
% Curated high-confidence ligand-receptor pairs from literature and databases.
% Used as a fallback when CellPhoneDB extraction yields insufficient pairs.
% Pairs are drawn from well-characterized signaling pathways across multiple
% cell types and tissue contexts.

    pairs = {
        % Angiogenesis
        'VEGFA', 'FLT1'; 'VEGFA', 'KDR'; 'VEGFB', 'FLT1'; 'VEGFC', 'FLT4';
        'ANGPT1', 'TEK'; 'ANGPT2', 'TEK';

        % EGFR pathway
        'EGF', 'EGFR'; 'TGFA', 'EGFR'; 'HBEGF', 'EGFR';
        'AREG', 'EGFR'; 'EREG', 'EGFR'; 'EPGN', 'EGFR';

        % PDGF pathway
        'PDGFA', 'PDGFRA'; 'PDGFB', 'PDGFRB';
        'PDGFC', 'PDGFRA'; 'PDGFD', 'PDGFRB';

        % FGF pathway
        'FGF1', 'FGFR1'; 'FGF2', 'FGFR1'; 'FGF7', 'FGFR2';
        'FGF10', 'FGFR2'; 'FGF19', 'FGFR4';

        % HGF/MET
        'HGF', 'MET';

        % Cytokines and immune signaling
        'TNF', 'TNFRSF1A'; 'TNF', 'TNFRSF1B'; 'LTA', 'TNFRSF1A';
        'IL6', 'IL6R'; 'IL1B', 'IL1R1'; 'IL1A', 'IL1R1';
        'IL8', 'CXCR1'; 'IL10', 'IL10RA';

        % TGF-beta
        'TGFB1', 'TGFBR1'; 'TGFB1', 'TGFBR2';
        'TGFB2', 'TGFBR1'; 'TGFB3', 'TGFBR2';

        % Immune checkpoint
        'CD274', 'PDCD1'; 'PDCD1LG2', 'PDCD1';
        'CTLA4', 'CD80'; 'CD28', 'CD86';

        % Notch pathway
        'DLL1', 'NOTCH1'; 'DLL4', 'NOTCH1';
        'JAG1', 'NOTCH1'; 'JAG2', 'NOTCH2';
        'JAG1', 'NOTCH2'; 'DLL1', 'NOTCH2';

        % Wnt signaling
        'WNT1', 'FZD1'; 'WNT3A', 'FZD1'; 'WNT5A', 'FZD5'; 'WNT7A', 'FZD5';

        % BMP signaling
        'BMP2', 'BMPR1A'; 'BMP4', 'BMPR1A'; 'BMP7', 'BMPR2';

        % Chemokines
        'CXCL12', 'CXCR4'; 'CCL2', 'CCR2'; 'CCL5', 'CCR5';
        'CXCL10', 'CXCR3'; 'CCL3', 'CCR1'; 'CCL4', 'CCR5';

        % Macrophage/myeloid signaling
        'CSF1', 'CSF1R';

        % Ephrin signaling
        'EFNA1', 'EPHA2'; 'EFNB2', 'EPHB4';

        % Semaphorins
        'SEMA3A', 'NRP1'; 'SEMA4D', 'PLXNB1';
    };
end
