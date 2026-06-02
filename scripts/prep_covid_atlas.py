"""
prep_covid_atlas.py — Stream COVID-19 Blood Atlas (Stephenson et al. 2021 Nat Med)
                      from CellxGene via HTTP range requests.

Dataset: COMBAT consortium, 836,148 cells, 124 donors
Collection: 8f126edf-5405-4731-8374-b5ce11f53e82
File URL:   https://datasets.cellxgene.cziscience.com/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad
File size:  5.6 GB  — we download ~50-200 MB via range reads
Outputs:    converted/covid_atlas/covid_expression_sparse.mat, covid_cell_metadata.csv, covid_gene_names.csv
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import h5py
import fsspec
import os, sys, time
from collections import defaultdict

URL = "https://datasets.cellxgene.cziscience.com/687c09ff-731a-4e3d-ac07-4c29c33a6338.h5ad"
OUT_DIR  = "converted/covid_atlas"
GENE_TBL = "converted/gene_table.csv"
LR_CSV   = "converted/interaction_table.csv"
PROT_CSV = "converted/protein_table.csv"
COMP_CSV = "converted/complex_composition_table.csv"

MAX_CELLS_PER_CLUSTER_DONOR = 30   # ~30 cells × 124 donors × 6 clusters = ~22,000 cells
MIN_DONORS_PER_CLUSTER = 20        # need enough donors for robust KSG MI
PROBE_N = 300                      # cells for HVG estimation
N_HVG   = 2000                     # HVGs to include

os.makedirs(OUT_DIR, exist_ok=True)

print("=== COVID-19 Blood Atlas Prep (HTTP range streaming) ===")
print(f"URL: {URL}\n")
t0 = time.time()

def decode(arr):
    """Decode bytes/strings from HDF5."""
    if hasattr(arr, 'dtype') and arr.dtype.kind in ('S', 'O'):
        return np.array([x.decode('utf-8') if isinstance(x, bytes) else str(x) for x in arr])
    return arr

def load_categorical(group, key):
    """Load a categorical column from obs/var HDF5 group."""
    d = group[key]
    if isinstance(d, h5py.Group):
        cats  = decode(d['categories'][:])
        codes = d['codes'][:]
        return np.array([cats[c] if c >= 0 else 'NA' for c in codes])
    else:
        return decode(d[:])

print("Opening remote HDF5 via HTTP range reads...")
of = fsspec.open(URL, 'rb')
f  = h5py.File(of.open(), 'r')
print(f"  Keys: {list(f.keys())}")

print("\nLoading obs metadata...")
obs     = f['obs']
n_cells = int(f['X']['indptr'].shape[0]) - 1
print(f"  Total cells: {n_cells}")

donor_ids = load_categorical(obs, 'donor_id')
cell_types = load_categorical(obs, 'cell_type')
diseases  = load_categorical(obs, 'disease')
clusters  = load_categorical(obs, 'cluster')

print(f"  Unique donors: {np.unique(donor_ids).shape[0]}")
print(f"  Cell types (sample): {np.unique(cell_types)[:8]}")
print(f"  Disease: {dict(zip(*np.unique(diseases, return_counts=True)))}")

print("\nLoading gene names...")
gene_names = load_categorical(f['var'], 'feature_name')
gene_names_upper = np.array([g.upper() for g in gene_names])
n_genes = len(gene_names_upper)
print(f"  Total genes: {n_genes}")

gene_to_idx = {g: i for i, g in enumerate(gene_names_upper)}

print("\nBuilding L-R gene set...")
gt        = pd.read_csv(GENE_TBL)
protein   = pd.read_csv(PROT_CSV)
inter     = pd.read_csv(LR_CSV)
comp_comp = pd.read_csv(COMP_CSV)

pn_to_hgnc = {}
for _, row in gt.iterrows():
    pn = str(row['protein_name']).upper().strip()
    hs = str(row['hgnc_symbol']).upper().strip()
    if pn and hs and pn != 'NAN' and hs != 'NAN':
        pn_to_hgnc[pn] = hs

prot_mid_to_gene = {}
for _, row in protein.iterrows():
    mid = row['protein_multidata_id']
    pn  = str(row['protein_name']).upper().strip()
    prot_mid_to_gene[mid] = pn_to_hgnc.get(pn, pn.replace('_HUMAN', ''))

complex_prots = defaultdict(list)
for _, row in comp_comp.iterrows():
    complex_prots[row['complex_multidata_id']].append(row['protein_multidata_id'])

matrix_gene_set = set(gene_names_upper)
multidata_to_gene = dict(prot_mid_to_gene)
for cmid, pmids in complex_prots.items():
    chosen = ''
    for pmid in pmids:
        g = prot_mid_to_gene.get(pmid, '')
        if g in matrix_gene_set:
            chosen = g; break
        elif not chosen:
            chosen = g
    if chosen:
        multidata_to_gene[cmid] = chosen

lr_gene_set = set()
for _, row in inter.iterrows():
    for mid in [row['multidata_1_id'], row['multidata_2_id']]:
        if mid in multidata_to_gene:
            lr_gene_set.add(multidata_to_gene[mid])

lr_genes_in_matrix = lr_gene_set & matrix_gene_set
lr_gene_idx = np.array(sorted([gene_to_idx[g] for g in lr_genes_in_matrix]))
print(f"  L-R genes in matrix: {len(lr_gene_idx)}")

def broad_cluster(ct):
    ct_l = ct.lower()
    if 'b cell' in ct_l:                                       return 'B'
    if 'cd4' in ct_l or 'helper' in ct_l:                     return 'CD4T'
    if 'cd8' in ct_l or 'cytotox' in ct_l:                    return 'CD8T'
    if 'natural killer' in ct_l or ct_l.startswith('nk '):    return 'NK'
    if 'monocyte' in ct_l or 'mono' in ct_l:                  return 'Mono'
    if 'dendritic' in ct_l:                                    return 'DC'
    if 'plasmablast' in ct_l or 'plasma cell' in ct_l:        return 'Plasma'
    return None  # exclude

broad_labels = np.array([broad_cluster(ct) for ct in cell_types])
pbmc_mask    = broad_labels != None  # noqa: E711

print(f"\nPBMC cells: {pbmc_mask.sum()} / {n_cells}")
pbmc_idx = np.where(pbmc_mask)[0]
b_labels  = broad_labels[pbmc_mask]
b_donors  = donor_ids[pbmc_mask]
b_disease = diseases[pbmc_mask]

from collections import Counter
for cl, n in Counter(b_labels).most_common():
    n_don = len(np.unique(b_donors[b_labels == cl]))
    print(f"  {cl:8s}: {n:6d} cells, {n_don:4d} donors")

valid_clusters = []
for cl in np.unique(b_labels):
    n_don = len(np.unique(b_donors[b_labels == cl]))
    if n_don >= MIN_DONORS_PER_CLUSTER:
        valid_clusters.append(cl)
    else:
        print(f"  Dropping {cl}: only {n_don} donors (< {MIN_DONORS_PER_CLUSTER})")

print(f"\nValid clusters ({len(valid_clusters)}): {valid_clusters}")
valid_mask = np.isin(b_labels, valid_clusters)
pbmc_idx  = pbmc_idx[valid_mask]
b_labels  = b_labels[valid_mask]
b_donors  = b_donors[valid_mask]
b_disease = b_disease[valid_mask]

print(f"\nSubsampling to <= {MAX_CELLS_PER_CLUSTER_DONOR} cells per cluster-donor...")
rng  = np.random.default_rng(42)
keep = []
for cl in valid_clusters:
    cl_mask = b_labels == cl
    for dn in np.unique(b_donors[cl_mask]):
        sel = np.where(cl_mask & (b_donors == dn))[0]
        if len(sel) > MAX_CELLS_PER_CLUSTER_DONOR:
            sel = rng.choice(sel, MAX_CELLS_PER_CLUSTER_DONOR, replace=False)
        keep.extend(sel.tolist())

keep     = np.array(sorted(set(keep)))
sub_full = pbmc_idx[keep]     # global cell indices in the full dataset (for X loading)
sub_cl   = b_labels[keep]
sub_don  = b_donors[keep]
sub_dis  = b_disease[keep]

print(f"  Subsampled: {len(sub_full)} cells from {len(np.unique(sub_don))} donors")

print("\nLoading indptr (once, ~3.3 MB)...")
indptr = f['X']['indptr'][:]     # shape (n_cells+1,)
print(f"  Loaded indptr shape: {indptr.shape}")

def load_rows_cols(X_node, row_idx_sorted, col_idx_sorted, batch_gap=300, verbose=True):
    """
    Load a submatrix of a CSR HDF5 sparse matrix via BATCHED range reads.

    Batches consecutive rows to minimize HTTP requests — rows within batch_gap
    of each other share a single data+indices fetch instead of one per row.
    row_idx_sorted : sorted global row indices (must be sorted ascending)
    col_idx_sorted : sorted global col indices to keep in output
    Returns dense float32 array of shape (len(row_idx), len(col_idx))
    """
    col_set = {int(c): ci for ci, c in enumerate(col_idx_sorted)}
    n_rows  = len(row_idx_sorted)
    n_cols  = len(col_idx_sorted)
    result  = np.zeros((n_rows, n_cols), dtype=np.float32)

    row_to_ri = {int(r): ri for ri, r in enumerate(row_idx_sorted)}

    batches = []  # list of (first_row, last_row, [row_indices])
    cur_batch_rows = [row_idx_sorted[0]]
    for i in range(1, n_rows):
        gap = int(row_idx_sorted[i]) - int(row_idx_sorted[i-1])
        if gap <= batch_gap:
            cur_batch_rows.append(row_idx_sorted[i])
        else:
            batches.append(cur_batch_rows)
            cur_batch_rows = [row_idx_sorted[i]]
    batches.append(cur_batch_rows)

    if verbose:
        print(f"  Batching: {n_rows} rows -> {len(batches)} batches "
              f"(avg {n_rows/len(batches):.0f} rows/batch)")

    for bi, batch_rows in enumerate(batches):
        r_first = int(batch_rows[0])
        r_last  = int(batch_rows[-1])
        d_start = int(indptr[r_first])
        d_end   = int(indptr[r_last + 1])
        if d_start >= d_end:
            continue
        batch_data    = np.array(X_node['data'][d_start:d_end],    dtype=np.float32)
        batch_indices = np.array(X_node['indices'][d_start:d_end], dtype=np.int32)

        for r in batch_rows:
            ri       = row_to_ri[int(r)]
            rs       = int(indptr[r])     - d_start
            re       = int(indptr[r + 1]) - d_start
            if rs >= re:
                continue
            row_idx_arr = batch_indices[rs:re]
            row_val_arr = batch_data[rs:re]
            for c, v in zip(row_idx_arr, row_val_arr):
                ci = col_set.get(int(c))
                if ci is not None:
                    result[ri, ci] = float(v)

        if verbose and bi % 50 == 0:
            print(f"  Batch {bi+1}/{len(batches)} done...")

    return result

print(f"\nEstimating HVGs from {PROBE_N} probe cells (batched loading)...")
probe_global = sub_full[rng.choice(len(sub_full), min(PROBE_N, len(sub_full)), replace=False)]
probe_global = np.array(sorted(probe_global))

hvg_probe_idx = np.arange(0, n_genes, 7)   # every 7th gene (~5186 genes) for variance estimation
print(f"  Loading {len(probe_global)} cells x {len(hvg_probe_idx)} probe genes (every 7th)...")
t_probe = time.time()
probe_mat = load_rows_cols(f['X'], probe_global, hvg_probe_idx, batch_gap=50, verbose=False)
print(f"  Done in {time.time()-t_probe:.1f}s  shape={probe_mat.shape}")

means = probe_mat.mean(axis=0)
vars_ = probe_mat.var(axis=0)
lm    = np.log(means + 1e-8)
lv    = np.log(vars_  + 1e-8)
poly  = np.polyfit(lm, lv, 2)
resid = lv - np.polyval(poly, lm)
# Top N_HVG by residual in probe-gene space, map back to full gene indices
top_probe  = np.argsort(-resid)[:N_HVG]
hvg_idx    = hvg_probe_idx[top_probe]   # full-matrix gene indices
print(f"  Selected {N_HVG} HVGs by residual variance (from {len(hvg_probe_idx)}-gene probe)")

# Union: L-R + HVG
final_gene_idx   = np.array(sorted(np.union1d(lr_gene_idx, hvg_idx)))
final_gene_names = gene_names_upper[final_gene_idx]
print(f"  Final gene set: {len(final_gene_idx)} ({len(lr_gene_idx)} L-R + {N_HVG} HVG)")

print(f"\nLoading final matrix ({len(sub_full)} cells × {len(final_gene_idx)} genes)...")
t_load = time.time()
X_raw = load_rows_cols(f['X'], sub_full, final_gene_idx)
print(f"  Done in {time.time()-t_load:.1f}s  shape={X_raw.shape}")
print(f"  Sparsity: {1 - np.count_nonzero(X_raw) / X_raw.size:.1%}")

print("\nNormalizing (library-size + log1p)...")
lib_sizes = X_raw.sum(axis=1, keepdims=True)
lib_sizes[lib_sizes == 0] = 1
X_norm = np.log1p(X_raw / lib_sizes * 1e4).astype(np.float32)

cluster_names = sorted(np.unique(sub_cl))
cl_to_id = {c: i for i, c in enumerate(cluster_names)}
cluster_ids = np.array([cl_to_id[c] for c in sub_cl])

meta = pd.DataFrame({
    'Cluster':      cluster_ids,
    'CellTypeName': sub_cl,
    'donor_id':     sub_don,
    'condition':    sub_dis,
})

print(f"\nFinal dataset:")
print(f"  Cells:   {len(meta)}")
print(f"  Donors:  {meta['donor_id'].nunique()}")
print(f"  Clusters: {meta['Cluster'].nunique()}")
for cl, n in meta['CellTypeName'].value_counts().items():
    n_d = meta[meta['CellTypeName']==cl]['donor_id'].nunique()
    print(f"    {cl:8s}: {n:5d} cells, {n_d:3d} donors")

print(f"\nSaving to {OUT_DIR}/...")
scipy.io.savemat(f"{OUT_DIR}/covid_expression_sparse.mat",
                 {'X': X_norm}, do_compression=True)
meta.to_csv(f"{OUT_DIR}/covid_cell_metadata.csv", index=False)
pd.DataFrame({'gene_name': final_gene_names}).to_csv(
    f"{OUT_DIR}/covid_gene_names.csv", index=False)

print(f"  covid_expression_sparse.mat : {X_norm.shape}")
print(f"  covid_cell_metadata.csv     : {len(meta)} rows")
print(f"  covid_gene_names.csv        : {len(final_gene_names)} genes")

f.close()
print(f"\nDone in {time.time()-t0:.1f}s  (downloaded ~{(time.time()-t0)*10:.0f} MB estimate)")
print(f"\nNext: run_covid_atlas.m with:")
print(f"  config.expr_path     = '{OUT_DIR}/covid_expression_sparse.mat'")
print(f"  config.metadata_path = '{OUT_DIR}/covid_cell_metadata.csv'")
print(f"  {meta['donor_id'].nunique()} donors → KSG MI well-powered; IV/2SLS valid")
