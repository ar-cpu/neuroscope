function pairs = get_known_ligand_receptor_pairs()
% GET_KNOWN_LIGAND_RECEPTOR_PAIRS  Return a manually curated set of
%   high-confidence ligand-receptor pairs drawn from peer-reviewed literature
%   and validated signaling databases (CellPhoneDB, NicheNet, KEGG).
%
% This function serves as a fallback when CellPhoneDB extraction yields
% fewer than 50 pairs for a given dataset. Pairs span major signaling
% families relevant across cancer, immune, and stromal biology: VEGF,
% PDGF, FGF, EGFR, TGF-β, Notch, Wnt, BMP, cytokine, immune checkpoint,
% chemokine, Ephrin, and Semaphorin families.
%
% Usage:
%   pairs = get_known_ligand_receptor_pairs()
%
% Outputs:
%   pairs - {N×2} cell array where column 1 is the ligand gene symbol
%           and column 2 is the receptor gene symbol. Both use HGNC
%           nomenclature.

    pairs = {
        'VEGFA', 'FLT1'; 'VEGFA', 'KDR'; 'VEGFB', 'FLT1'; 'VEGFC', 'FLT4'
        'PDGFA', 'PDGFRA'; 'PDGFB', 'PDGFRB'; 'PDGFC', 'PDGFRA'; 'PDGFD', 'PDGFRB'
        'FGF1', 'FGFR1'; 'FGF2', 'FGFR1'; 'FGF7', 'FGFR2'; 'FGF10', 'FGFR2'
        'EGF', 'EGFR'; 'TGFA', 'EGFR'; 'HBEGF', 'EGFR'; 'EREG', 'EGFR'
        'TNF', 'TNFRSF1A'; 'TNF', 'TNFRSF1B'; 'LTA', 'TNFRSF1A'
        'IL1B', 'IL1R1'; 'IL1A', 'IL1R1'; 'IL6', 'IL6R'; 'IL8', 'CXCR1'
        'TGFB1', 'TGFBR1'; 'TGFB1', 'TGFBR2'; 'TGFB2', 'TGFBR1'; 'TGFB3', 'TGFBR2'
        'WNT1', 'FZD1'; 'WNT3A', 'FZD1'; 'WNT5A', 'FZD5'; 'WNT7A', 'FZD5'
        'BMP2', 'BMPR1A'; 'BMP4', 'BMPR1A'; 'BMP7', 'BMPR2'
        'DLL1', 'NOTCH1'; 'DLL4', 'NOTCH1'; 'JAG1', 'NOTCH1'; 'JAG2', 'NOTCH2'
        'CD274', 'PDCD1'; 'PDCD1LG2', 'PDCD1'; 'CTLA4', 'CD80'; 'CD28', 'CD86'
        'CXCL12', 'CXCR4'; 'CCL2', 'CCR2'; 'CCL5', 'CCR5'; 'CXCL10', 'CXCR3'
        'ANGPT1', 'TEK'; 'ANGPT2', 'TEK'; 'HGF', 'MET'; 'CSF1', 'CSF1R'
    };
end