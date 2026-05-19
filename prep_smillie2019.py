"""
prep_smillie2019.py  —  Smillie et al. 2019 Ulcerative Colitis colon (GSE134809)
  Non-blood tissue: shows NEUROSCOPE generalizes beyond PBMC.
  ~368K cells, 12 healthy + 18 UC donors, epithelial + stromal + immune.
  Streams from CellxGene CELLxGENE collection.

Dataset: "Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis"
  Cell 178(3):714-730.e22 (Smillie et al. 2019)
  CellxGene: https://cellxgene.cziscience.com/collections/

Outputs:
  converted/smillie2019/smillie_expression_sparse.mat
  converted/smillie2019/smillie_cell_metadata.csv
  converted/smillie2019/smillie_gene_names.csv

Run: python prep_smillie2019.py
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import h5py
import fsspec
import os, sys, time
from collections import Counter

# CellxGene Smillie 2019 — Ulcerative Colitis colon
URL = "https://datasets.cellxgene.cziscience.com/444b5b74-4b23-4155-abe3-a61e1b74e722.h5ad"
OUT_DIR  = "converted/smillie2019"
GENE_TBL = "converted/gene_table.csv"
LR_CSV   = "converted/interaction_table.csv"

MAX_CELLS_PER_CLUSTER_DONOR = 40
MIN_DONORS_PER_CLUSTER      = 8   # fewer donors than PBMC datasets — colon biopsy
PROBE_N = 300

print("=== Smillie 2019 Ulcerative Colitis Colon prep ===\n")
os.makedirs(OUT_DIR, exist_ok=True)

print(f"Opening remote H5AD: {URL}")
t0 = time.time()
fs = fsspec.filesystem("http")
f  = h5py.File(fs.open(URL, "rb"), "r")
print(f"  Opened in {time.time()-t0:.1f}s")

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

print("\nReading cell metadata...")
obs_keys = list(f["obs"].keys())
print(f"  obs columns: {obs_keys}")

cell_type = decode_col(f, "obs/cell_type")
donor_id  = decode_col(f, "obs/donor_id")
if donor_id is None:
    donor_id = decode_col(f, "obs/Subject")
disease   = decode_col(f, "obs/disease")
if disease is None:
    disease = decode_col(f, "obs/Health")

n_cells_total = len(cell_type) if cell_type else 0
print(f"  Total cells: {n_cells_total:,}")
if donor_id:
    print(f"  Unique donors: {len(set(donor_id))}")
if disease:
    print(f"  Conditions: {dict(Counter(disease))}")

# ── Gene list ─────────────────────────────────────────────────────────────────
print("\nReading gene list...")
var_keys = list(f["var"].keys())
gene_symbols = None
for key in ["feature_name", "gene_symbols", "gene_name", "_index"]:
    if key in f["var"]:
        gene_symbols = decode_col(f, f"var/{key}")
        if gene_symbols:
            print(f"  Gene symbols from var/{key}: {len(gene_symbols)}")
            break
if gene_symbols is None:
    gene_symbols = [str(i) for i in range(n_cells_total)]
gene_symbols = [g.upper() for g in gene_symbols]
n_genes_total = len(gene_symbols)

# ── L-R genes ─────────────────────────────────────────────────────────────────
lr_genes = set()
if os.path.exists(LR_CSV):
    lr_df = pd.read_csv(LR_CSV)
    for col in lr_df.columns:
        lr_genes.update(lr_df[col].dropna().astype(str).str.upper().tolist())
gene_idx_map = {g: i for i, g in enumerate(gene_symbols)}
lr_keep_idx  = sorted([gene_idx_map[g] for g in lr_genes if g in gene_idx_map])
print(f"  L-R genes found: {len(lr_keep_idx)}")

# ── HVG probe ─────────────────────────────────────────────────────────────────
print(f"\nProbing {PROBE_N} cells for HVG...")
probe_idx = np.linspace(0, n_cells_total - 1, PROBE_N, dtype=int)
is_sparse_group = "X/data" in f

if is_sparse_group:
    data_all  = f["X/data"][()]
    indices_a = f["X/indices"][()]
    indptr_a  = f["X/indptr"][()]
    probe_rows = []
    for ci in probe_idx:
        s, e = indptr_a[ci], indptr_a[ci+1]
        row = np.zeros(n_genes_total, dtype=np.float32)
        row[indices_a[s:e]] = data_all[s:e]
        probe_rows.append(row)
    probe_mat = np.vstack(probe_rows)
else:
    probe_mat = np.asarray(f["X"][probe_idx, :])

probe_means = probe_mat.mean(0); probe_vars = probe_mat.var(0)
eps = 1e-8
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    lm = np.log(probe_means + eps); lv = np.log(probe_vars + eps)
    coeffs = np.polyfit(lm, lv, 2)
    residuals = lv - np.polyval(coeffs, lm)
top_hvg_idx = np.argsort(residuals)[::-1][:3000]
keep_idx = np.union1d(top_hvg_idx, lr_keep_idx).astype(int)
keep_idx.sort()
print(f"  Keeping {len(keep_idx)} genes")

# ── Sample cells ──────────────────────────────────────────────────────────────
print("\nSampling cells (balanced by cluster × donor)...")

# UC condition cells for condition, healthy for reference
if disease:
    uc_terms  = [d for d in set(disease) if any(k in d.lower() for k in
                 ['ulcerative', 'uc ', 'colitis', 'inflammatory', 'ibd'])]
    if not uc_terms:
        # fallback: non-healthy
        health_terms = [d for d in set(disease) if any(k in d.lower() for k in
                        ['normal','healthy','control'])]
        uc_terms = [d for d in set(disease) if d not in health_terms]
    print(f"  UC terms: {uc_terms}")
else:
    uc_terms = []

if cell_type and donor_id:
    df_meta = pd.DataFrame({
        "cell_type": cell_type,
        "donor_id":  donor_id,
        "condition": disease if disease else ["unknown"]*n_cells_total,
    })
    cluster_enc, _ = pd.factorize(df_meta["cell_type"])
    df_meta["Cluster"] = cluster_enc + 1

    if uc_terms:
        cond_df = df_meta[df_meta["condition"].isin(uc_terms)].copy()
    else:
        cond_df = df_meta.copy()

    cluster_donor_counts = cond_df.groupby("cell_type")["donor_id"].nunique()
    valid_clusters = cluster_donor_counts[
        cluster_donor_counts >= MIN_DONORS_PER_CLUSTER].index.tolist()
    print(f"  Clusters with ≥{MIN_DONORS_PER_CLUSTER} donors: {len(valid_clusters)}")

    keep_cells = []
    for cl in valid_clusters:
        for don, grp in cond_df[cond_df["cell_type"]==cl].groupby("donor_id"):
            sampled = grp.sample(min(MAX_CELLS_PER_CLUSTER_DONOR, len(grp)), random_state=42)
            keep_cells.extend(sampled.index.tolist())

    keep_cells = sorted(set(keep_cells))
    print(f"  Selected {len(keep_cells):,} cells from {df_meta.loc[keep_cells,'donor_id'].nunique()} donors")
    meta_out = df_meta.loc[keep_cells, ["Cluster","donor_id","condition"]].copy()
    meta_out.insert(0, "CellID", np.arange(1, len(keep_cells)+1))
    meta_out.columns = ["CellID","Cluster","donor_id","condition"]
else:
    keep_cells = list(range(min(n_cells_total, 20000)))
    meta_out = pd.DataFrame({"CellID": np.arange(1, len(keep_cells)+1),
                              "Cluster": np.ones(len(keep_cells), int),
                              "donor_id": ["unknown"]*len(keep_cells)})

# ── Load expression ───────────────────────────────────────────────────────────
print(f"\nLoading {len(keep_cells):,} cells × {len(keep_idx)} genes...")
t1 = time.time()
keep_cells = np.array(keep_cells)
rows_list, cols_list, vals_list = [], [], []
chunk = 500

for i in range(0, len(keep_cells), chunk):
    batch_ci = keep_cells[i:i+chunk]
    if i % 5000 == 0:
        print(f"  {i:,}/{len(keep_cells):,} ({100*i/len(keep_cells):.0f}%)", flush=True)
    if is_sparse_group:
        for local_j, ci in enumerate(batch_ci):
            s, e = indptr_a[ci], indptr_a[ci+1]
            ri = indices_a[s:e]; rv = data_all[s:e]
            mask = np.isin(ri, keep_idx)
            if mask.any():
                mapped = np.searchsorted(keep_idx, ri[mask])
                rows_list.append(np.full(mask.sum(), i+local_j, dtype=np.int32))
                cols_list.append(mapped.astype(np.int32))
                vals_list.append(rv[mask].astype(np.float32))
    else:
        block = np.asarray(f["X"][batch_ci.tolist(), :])[:, keep_idx]
        r, c = block.nonzero()
        rows_list.append((r+i).astype(np.int32))
        cols_list.append(c.astype(np.int32))
        vals_list.append(block[r,c].astype(np.float32))

X_sparse = sp.csc_matrix(
    (np.concatenate(vals_list), (np.concatenate(rows_list), np.concatenate(cols_list))),
    shape=(len(keep_cells), len(keep_idx)))
print(f"  Done in {time.time()-t1:.1f}s  nnz={X_sparse.nnz:,}")
f.close()

# ── Save ──────────────────────────────────────────────────────────────────────
print("\nSaving...")
scipy.io.savemat(os.path.join(OUT_DIR, "smillie_expression_sparse.mat"),
                 {"expr_matrix": X_sparse}, do_compression=True)
meta_out.to_csv(os.path.join(OUT_DIR, "smillie_cell_metadata.csv"), index=False)
pd.DataFrame({"gene_name": [gene_symbols[i] for i in keep_idx]}).to_csv(
    os.path.join(OUT_DIR, "smillie_gene_names.csv"), index=False)
print(f"Done. {len(keep_cells):,} cells, {len(keep_idx)} genes, "
      f"{meta_out['donor_id'].nunique()} donors")
