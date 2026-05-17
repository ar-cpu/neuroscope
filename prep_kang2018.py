"""
prep_kang2018.py  —  Convert Kang et al. 2018 (GSE96583) PBMC data
to NEUROSCOPE format.

Dataset: Kang et al. 2018 Nature Biotechnology, 8 donors × 2 conditions (ctrl + IFN-beta stim),
  ~25K PBMCs after QC, cell types: CD4 T, CD8 T, B, NK, CD14+ Mono, FCGR3A+ Mono, DC, MK
Outputs:  converted/kang2018/kang_expression_sparse.mat
          converted/kang2018/kang_cell_metadata.csv    (Cluster, donor_id)
          converted/kang2018/kang_gene_names.csv
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import scipy.io as sio
import gzip, os, sys
import warnings
warnings.filterwarnings('ignore')

BASE = 'c:/arjun-project/neuroscope/converted/kang2018'
OUT_DIR = BASE

print("=== Kang 2018 PBMC prep ===\n")

print("Loading metadata...")
meta = pd.read_csv(f'{BASE}/batch2_tsne.tsv.gz', sep='\t', index_col=0)
meta = meta[meta['multiplets'] == 'singlet'].copy()
print(f"  Singlet cells: {len(meta)}")
print(f"  Donors: {meta['ind'].nunique()}  = {sorted(meta['ind'].unique())}")
print(f"  Conditions: {meta['stim'].unique()}")
print(f"  Cell types: {meta['cell'].value_counts().to_dict()}")
print()

print("Loading gene names...")
genes_df = pd.read_csv(f'{BASE}/batch2_genes.tsv.gz', sep='\t', header=None,
                        names=['ensembl', 'symbol'])
gene_symbols = genes_df['symbol'].str.upper().values
n_genes = len(gene_symbols)
print(f"  Genes: {n_genes}")

print("Loading MTX count matrices (batch 2.1 and 2.2)...")

def load_mtx_gz(mtx_path, barcodes_path):
    """Load a gzipped MTX file and return (csc_matrix, barcode_list)."""
    with gzip.open(mtx_path, 'rt') as f:
        while True:
            line = f.readline()
            if not line.startswith('%'):
                break
        n_genes, n_cells, n_entries = [int(x) for x in line.strip().split()]
        print(f"    MTX: {n_genes} genes × {n_cells} cells, {n_entries} entries")
        rows, cols, data = [], [], []
        for line in f:
            parts = line.strip().split()
            rows.append(int(parts[0])-1)
            cols.append(int(parts[1])-1)
            data.append(float(parts[2]))
    mat = sp.csc_matrix((data, (rows, cols)), shape=(n_genes, n_cells))
    with gzip.open(barcodes_path, 'rt') as f:
        barcodes = [l.strip() for l in f]
    return mat, barcodes

mat_21, bc_21 = load_mtx_gz(f'{BASE}/GSM2560248_2.1.mtx.gz', f'{BASE}/GSM2560248_barcodes.tsv.gz')
mat_22, bc_22 = load_mtx_gz(f'{BASE}/GSM2560249_2.2.mtx.gz', f'{BASE}/GSM2560249_barcodes.tsv.gz')

print(f"  Batch 2.1: {mat_21.shape}")
print(f"  Batch 2.2: {mat_22.shape}")

bc_21_arr = np.array(bc_21)
bc_22_arr = np.array(bc_22)

meta_barcodes = set(meta.index.tolist())
in_21 = np.array([b in meta_barcodes for b in bc_21])
in_22 = np.array([b in meta_barcodes for b in bc_22])
print(f"\n  Barcodes in metadata: {in_21.sum()} from 2.1, {in_22.sum()} from 2.2")

idx_21 = np.where(in_21)[0]
idx_22 = np.where(in_22)[0]

X_21 = mat_21[:, idx_21].T  # cells × genes
X_22 = mat_22[:, idx_22].T

print(f"  Retained: {X_21.shape[0]} from batch 2.1, {X_22.shape[0]} from batch 2.2")

bc_21_kept = bc_21_arr[idx_21]
bc_22_kept = bc_22_arr[idx_22]
all_barcodes = np.concatenate([bc_21_kept, bc_22_kept])

X_all = sp.vstack([X_21, X_22])
print(f"  Combined matrix: {X_all.shape}  (cells × genes)")

meta_21 = meta.loc[bc_21_kept]
meta_22 = meta.loc[bc_22_kept]
meta_combined = pd.concat([meta_21, meta_22])
print(f"\n  Combined metadata: {meta_combined.shape}")
print(f"  Donors: {meta_combined['ind'].nunique()}")
print(f"  Conditions: {meta_combined['stim'].value_counts().to_dict()}")

print("\nNormalizing...")
X_dense = X_all.toarray().astype(np.float32)
lib_sizes = X_dense.sum(axis=1, keepdims=True)
lib_sizes = np.maximum(lib_sizes, 1)
X_norm = np.log1p(X_dense / lib_sizes * 1e4)  # log1p(TPM/10)
print(f"  Normalized matrix: {X_norm.shape}  max={X_norm.max():.2f}")

print("\nFiltering genes...")
lr_genes_set = set()
try:
    gt = pd.read_csv('c:/arjun-project/neuroscope/converted/gene_table.csv')
    lr_genes_set = set(gt['hgnc_symbol'].str.upper().dropna().values)
    lr_idx = np.where(np.isin(gene_symbols, list(lr_genes_set)))[0]
    print(f"  L-R genes in matrix: {len(lr_idx)}")
except Exception as e:
    print(f"  Warning: {e}")
    lr_idx = np.array([], dtype=int)

gene_var = np.var(X_norm, axis=0)
hvg_idx = np.argsort(gene_var)[::-1][:5000]

keep_cols = np.sort(np.union1d(hvg_idx, lr_idx))
print(f"  Keeping {len(keep_cols)} genes (HVG + L-R)")

X_filt = X_norm[:, keep_cols]
gene_names_filt = gene_symbols[keep_cols]

# donor_id = donor × condition (so each sample = distinct pseudobulk unit)
meta_combined = meta_combined.copy()
meta_combined['donor_sample'] = (meta_combined['ind'].astype(str) + '_' +
                                  meta_combined['stim'].astype(str))
print(f"\n  Pseudobulk donor-samples: {meta_combined['donor_sample'].nunique()}")
print(f"  Sample: {sorted(meta_combined['donor_sample'].unique())[:8]}")

ct_cats = pd.Categorical(meta_combined['cell'])
meta_out = pd.DataFrame({
    'Cluster': ct_cats.codes.astype(int),
    'CellTypeName': ct_cats,
    'donor_id': meta_combined['donor_sample'].values
})
meta_out = meta_out.reset_index(drop=True)

print(f"\n  Cluster mapping:")
print(meta_out.groupby(['Cluster','CellTypeName'])['donor_id'].nunique().rename('n_donors'))

print("\nSaving...")

meta_out.to_csv(f'{OUT_DIR}/kang_cell_metadata.csv', index=False)
print(f"  Metadata: {meta_out.shape}")

gene_df = pd.DataFrame({'gene_name': gene_names_filt})
gene_df.to_csv(f'{OUT_DIR}/kang_gene_names.csv', index=True)
print(f"  Genes: {len(gene_names_filt)}")

X_sparse = sp.csc_matrix(X_filt)
sio.savemat(f'{OUT_DIR}/kang_expression_sparse.mat',
            {'X': X_sparse}, do_compression=True)
sz = os.path.getsize(f'{OUT_DIR}/kang_expression_sparse.mat')/1e6
print(f"  Expression matrix: {sz:.1f} MB")

print("\n=== Prep complete ===")
print(f"  Cells: {X_filt.shape[0]}")
print(f"  Genes: {X_filt.shape[1]}")
print(f"  Donors (pseudobulk): {meta_out['donor_id'].nunique()}  (8 donors × 2 conditions)")
print(f"  Cell types: {meta_out['Cluster'].nunique()}")
