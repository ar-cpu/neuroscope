"""
prep_lupus.py  —  Perez et al. 2022 SLE PBMC (GSE174188 / CellxGene)
  Single-cell eQTL + cell-type analysis of systemic lupus erythematosus.
  162 donors (88 SLE, 49 healthy, 25 other), ~1.2M cells.
  Streams H5AD from CellxGene via HTTP range reads — no full download needed.

Outputs:
  converted/lupus/lupus_expression_sparse.mat
  converted/lupus/lupus_cell_metadata.csv   (CellID, Cluster, donor_id, condition)
  converted/lupus/lupus_gene_names.csv

Run: python prep_lupus.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import h5py
import fsspec
import os, sys, time
from collections import defaultdict

# CellxGene: Perez et al. 2022 Science — SLE PBMC
# Collection: 436154da-bcf1-4130-9c8b-120ff9a888f2
URL = "https://datasets.cellxgene.cziscience.com/97a17473-e2b1-4f31-a544-44a60773e2dd.h5ad"
OUT_DIR  = "converted/lupus"
GENE_TBL = "converted/gene_table.csv"
LR_CSV   = "converted/interaction_table.csv"

MAX_CELLS_PER_CLUSTER_DONOR = 30
MIN_DONORS_PER_CLUSTER      = 20
PROBE_N = 300

print("=== Lupus PBMC prep (Perez et al. 2022) ===\n")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Open remote H5AD ──────────────────────────────────────────────────────────
print(f"Opening remote H5AD: {URL}")
t0 = time.time()
fs = fsspec.filesystem("http")
f  = h5py.File(fs.open(URL, "rb"), "r")
print(f"  Opened in {time.time()-t0:.1f}s")

# ── Cell metadata ─────────────────────────────────────────────────────────────
print("\nReading cell metadata...")

def decode_col(f, path):
    if path not in f:
        return None
    d = f[path]
    if "categories" in d.attrs:
        cats = [c.decode() if isinstance(c, bytes) else str(c)
                for c in d.attrs["categories"]]
        codes = d[()]
        return [cats[c] for c in codes]
    raw = d[()]
    if raw.dtype.kind in ("S", "O"):
        return [r.decode() if isinstance(r, bytes) else str(r) for r in raw]
    return raw.tolist()

obs_keys = list(f["obs"].keys())
print(f"  obs columns: {obs_keys}")

cell_type  = decode_col(f, "obs/cell_type")
donor_id   = decode_col(f, "obs/donor_id")
if donor_id is None:
    donor_id = decode_col(f, "obs/individual")
disease    = decode_col(f, "obs/disease")
if disease is None:
    disease = decode_col(f, "obs/Status")

n_cells_total = len(cell_type) if cell_type else f["X/data"].shape[0]
print(f"  Total cells: {n_cells_total:,}")
if donor_id:
    print(f"  Unique donors: {len(set(donor_id))}")
if disease:
    from collections import Counter
    print(f"  Conditions: {dict(Counter(disease))}")

# ── Identify SLE condition ────────────────────────────────────────────────────
# Keep SLE-affected donors for condition, healthy for reference
if disease:
    unique_diseases = list(set(disease))
    print(f"\n  Disease values: {unique_diseases}")
    # CellxGene uses ontology terms; look for SLE/lupus or healthy
    sle_terms    = [d for d in unique_diseases if any(k in d.lower() for k in
                    ['lupus', 'sle', 'systemic lupus', 'erythematosus'])]
    health_terms = [d for d in unique_diseases if any(k in d.lower() for k in
                    ['normal', 'healthy', 'control', 'na'])]
    if not sle_terms:
        sle_terms = [d for d in unique_diseases if 'disease' in d.lower() or 'abnormal' in d.lower()]
    if not health_terms:
        health_terms = [d for d in unique_diseases if d not in sle_terms]
    print(f"  SLE terms: {sle_terms}")
    print(f"  Healthy terms: {health_terms}")
else:
    sle_terms = []
    health_terms = []

# ── Gene list ─────────────────────────────────────────────────────────────────
print("\nReading gene list...")
var_keys = list(f["var"].keys())
print(f"  var columns: {var_keys}")

gene_symbols = None
for key in ["feature_name", "gene_symbols", "gene_name", "_index"]:
    if key in f["var"]:
        gene_symbols = decode_col(f, f"var/{key}")
        if gene_symbols:
            print(f"  Gene symbols from var/{key}: {len(gene_symbols)} genes")
            break

if gene_symbols is None:
    gene_symbols = [str(i) for i in range(f["X"].shape[1] if "X" in f else
                    f["X/data"].shape[0])]

gene_symbols = [g.upper() for g in gene_symbols]
n_genes_total = len(gene_symbols)

# ── Identify L-R genes to keep ────────────────────────────────────────────────
print("\nIdentifying L-R genes in expression matrix...")
lr_genes = set()
if os.path.exists(LR_CSV):
    lr_df = pd.read_csv(LR_CSV)
    for col in lr_df.columns:
        vals = lr_df[col].dropna().astype(str).str.upper().tolist()
        lr_genes.update(vals)

gene_idx_map = {g: i for i, g in enumerate(gene_symbols)}
lr_keep_idx  = sorted([gene_idx_map[g] for g in lr_genes if g in gene_idx_map])
print(f"  L-R genes found: {len(lr_keep_idx)} / {len(lr_genes)}")

# Probe for HVGs from a sample
print(f"\nProbing {PROBE_N} cells for HVG estimation...")
probe_idx = np.linspace(0, n_cells_total - 1, PROBE_N, dtype=int)
X_key = "X" if "X" in f and hasattr(f["X"], "shape") else None
is_sparse_group = X_key is None and "X/data" in f

if is_sparse_group:
    data_all  = f["X/data"][()]
    indices_a = f["X/indices"][()]
    indptr_a  = f["X/indptr"][()]

    probe_rows = []
    for ci in probe_idx:
        start, end = indptr_a[ci], indptr_a[ci+1]
        row = np.zeros(n_genes_total, dtype=np.float32)
        row[indices_a[start:end]] = data_all[start:end]
        probe_rows.append(row)
    probe_mat = np.vstack(probe_rows)
else:
    probe_mat = f["X"][probe_idx, :]

probe_means = probe_mat.mean(axis=0)
probe_vars  = probe_mat.var(axis=0)
eps = 1e-8

import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    log_m = np.log(probe_means + eps)
    log_v = np.log(probe_vars  + eps)
    coeffs = np.polyfit(log_m, log_v, 2)
    fitted = np.polyval(coeffs, log_m)
    residuals = log_v - fitted

top_hvg_idx = np.argsort(residuals)[::-1][:3000]
keep_idx = np.union1d(top_hvg_idx, lr_keep_idx).astype(int)
keep_idx.sort()
print(f"  Keeping {len(keep_idx)} genes (3000 HVG + {len(lr_keep_idx)} L-R)")

# ── Sample cells ──────────────────────────────────────────────────────────────
print("\nSampling cells (balanced by cluster × donor)...")

# Use SLE condition for the main condition
if sle_terms and disease:
    condition_mask = [d in sle_terms for d in disease]
else:
    condition_mask = [True] * n_cells_total

if cell_type and donor_id:
    df_meta = pd.DataFrame({
        "cell_type": cell_type,
        "donor_id":  donor_id,
        "condition": disease if disease else ["unknown"] * n_cells_total,
    })
    df_meta["is_condition"] = condition_mask
    df_meta.index.name = "cell_idx"

    # Encode cluster
    cluster_enc, cluster_cats = pd.factorize(df_meta["cell_type"])
    df_meta["Cluster"] = cluster_enc + 1

    # Filter: condition cells only, clusters with ≥ MIN_DONORS_PER_CLUSTER donors
    cond_df = df_meta[df_meta["is_condition"]].copy()
    cluster_donor_counts = cond_df.groupby("cell_type")["donor_id"].nunique()
    valid_clusters = cluster_donor_counts[cluster_donor_counts >= MIN_DONORS_PER_CLUSTER].index.tolist()
    print(f"  Clusters with ≥{MIN_DONORS_PER_CLUSTER} donors: {len(valid_clusters)}")

    keep_cells = []
    for cl in valid_clusters:
        cl_cells = cond_df[cond_df["cell_type"] == cl]
        for don, grp in cl_cells.groupby("donor_id"):
            sampled = grp.sample(min(MAX_CELLS_PER_CLUSTER_DONOR, len(grp)),
                                 random_state=42)
            keep_cells.extend(sampled.index.tolist())

    keep_cells = sorted(set(keep_cells))
    print(f"  Selected {len(keep_cells):,} cells from {cond_df.loc[keep_cells, 'donor_id'].nunique()} donors")

    meta_out = df_meta.loc[keep_cells, ["Cluster", "donor_id", "condition"]].copy()
    meta_out.insert(0, "CellID", np.arange(1, len(keep_cells)+1))
    meta_out.columns = ["CellID", "Cluster", "donor_id", "condition"]
else:
    keep_cells = list(range(min(n_cells_total, 30000)))
    meta_out = pd.DataFrame({
        "CellID": np.arange(1, len(keep_cells)+1),
        "Cluster": np.ones(len(keep_cells), dtype=int),
        "donor_id": ["unknown"] * len(keep_cells),
    })

# ── Load expression for selected cells ───────────────────────────────────────
print(f"\nLoading expression for {len(keep_cells):,} cells × {len(keep_idx)} genes...")
t1 = time.time()
keep_cells = np.array(keep_cells)

rows_list, cols_list, vals_list = [], [], []
chunk = 500

for i in range(0, len(keep_cells), chunk):
    batch_ci = keep_cells[i:i+chunk]
    if i % 5000 == 0:
        pct = 100 * i / len(keep_cells)
        print(f"  {i:,}/{len(keep_cells):,} ({pct:.0f}%)  {time.time()-t1:.0f}s", flush=True)

    if is_sparse_group:
        for local_j, ci in enumerate(batch_ci):
            start, end = indptr_a[ci], indptr_a[ci+1]
            raw_idx  = indices_a[start:end]
            raw_vals = data_all[start:end]
            # keep only genes in keep_idx
            mask = np.isin(raw_idx, keep_idx)
            if mask.any():
                mapped = np.searchsorted(keep_idx, raw_idx[mask])
                rows_list.append(np.full(mask.sum(), i + local_j, dtype=np.int32))
                cols_list.append(mapped.astype(np.int32))
                vals_list.append(raw_vals[mask].astype(np.float32))
    else:
        block = f["X"][batch_ci.tolist(), :]
        block = np.asarray(block)[:, keep_idx]
        r, c = block.nonzero()
        rows_list.append((r + i).astype(np.int32))
        cols_list.append(c.astype(np.int32))
        vals_list.append(block[r, c].astype(np.float32))

rows_a = np.concatenate(rows_list)
cols_a = np.concatenate(cols_list)
vals_a = np.concatenate(vals_list)

n_out_cells = len(keep_cells)
n_out_genes = len(keep_idx)
X_sparse = sp.csc_matrix((vals_a, (rows_a, cols_a)),
                          shape=(n_out_cells, n_out_genes))
print(f"  Expression matrix: {n_out_cells} × {n_out_genes}  (nnz={X_sparse.nnz:,})")
print(f"  Load time: {time.time()-t1:.1f}s")

f.close()

# ── Save ──────────────────────────────────────────────────────────────────────
print("\nSaving outputs...")
mat_path  = os.path.join(OUT_DIR, "lupus_expression_sparse.mat")
meta_path = os.path.join(OUT_DIR, "lupus_cell_metadata.csv")
gene_path = os.path.join(OUT_DIR, "lupus_gene_names.csv")

scipy.io.savemat(mat_path, {"expr_matrix": X_sparse}, do_compression=True)
meta_out.to_csv(meta_path, index=False)
gene_names_out = [gene_symbols[i] for i in keep_idx]
pd.DataFrame({"gene_name": gene_names_out}).to_csv(gene_path, index=False)

print(f"  {mat_path}")
print(f"  {meta_path}")
print(f"  {gene_path}")
print(f"\nDone. {n_out_cells} cells, {n_out_genes} genes, "
      f"{meta_out['donor_id'].nunique()} donors")
