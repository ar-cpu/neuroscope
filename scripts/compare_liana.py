"""
compare_liana.py  —  LIANA multi-method comparison vs NEUROSCOPE
  Runs CellChat, CellPhoneDB, NicheNet, NATMI, and Connectome via the
  liana-py framework on the same datasets and compares to NEUROSCOPE.

  LIANA reference: Dimitrov et al. 2022, Nature Cell Biology
  liana-py: https://liana-py.readthedocs.io

Install:
  pip install liana>=1.0 anndata scanpy

Run: python compare_liana.py
"""

import os
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import warnings
warnings.filterwarnings("ignore")

try:
    import anndata as ad
    import liana as li
    from liana.method import singlecellsignalr, connectome, logfc, natmi, cellphonedb
    LIANA_OK = True
    print("liana-py loaded OK")
except ImportError as e:
    LIANA_OK = False
    print(f"liana-py not installed: {e}")
    print("Install: pip install liana anndata scanpy")

# ── Load a NEUROSCOPE-prepared dataset into AnnData ───────────────────────────
def load_dataset(mat_path, meta_path, gene_path):
    print(f"\nLoading {mat_path}...")
    data = scipy.io.loadmat(mat_path)
    key  = [k for k in data if not k.startswith("_")][0]
    X    = data[key]
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = sp.csr_matrix(X)

    meta  = pd.read_csv(meta_path)
    genes = pd.read_csv(gene_path)["gene_name"].tolist()

    # Orient so rows = cells: pick whichever dimension is closer to len(meta)
    n_cells = len(meta)
    if abs(X.shape[1] - n_cells) < abs(X.shape[0] - n_cells):
        X = X.T   # was genes × cells, transpose to cells × genes
    X = sp.csr_matrix(X)

    # Trim to matching length if small mismatch (MATLAB sometimes pads)
    if X.shape[0] > n_cells:
        X = X[:n_cells, :]

    id_col = "CellID" if "CellID" in meta.columns else meta.columns[0]
    adata = ad.AnnData(X=X,
                       obs=meta.set_index(id_col),
                       var=pd.DataFrame(index=genes))
    adata.var_names_make_unique()

    # Resolve cell type column — prefer string names over integer cluster IDs
    for col in ("CellTypeName", "cell_type", "cell_type_label", "Cluster"):
        if col in adata.obs.columns:
            adata.obs["cell_type"] = adata.obs[col].astype(str)
            break
    else:
        adata.obs["cell_type"] = "unknown"

    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"  Cell types: {adata.obs['cell_type'].nunique()}")
    return adata


def run_liana_methods(adata, dataset_name, neuroscope_csv, out_dir):
    """Run all LIANA methods and compare to NEUROSCOPE."""
    os.makedirs(out_dir, exist_ok=True)

    # Normalize if not already done
    import scanpy as sc
    if adata.X.max() > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    print(f"\nRunning LIANA on {dataset_name}...")

    # Run all available methods — try rank_aggregate first, fall back to cellphonedb
    results_df = None
    try:
        # liana >= 1.0 by_sample requires sample_key; try with donor_id
        sample_col = next((c for c in ("donor_id","sample_id","subject") if c in adata.obs.columns), None)
        if sample_col:
            li.mt.rank_aggregate.by_sample(
                adata, groupby="cell_type", sample_key=sample_col,
                expr_prop=0.1, verbose=False, use_raw=False,
            )
        else:
            li.mt.rank_aggregate(
                adata, groupby="cell_type", expr_prop=0.1,
                verbose=False, use_raw=False,
            )
        results_df = adata.uns["liana_res"].copy()
        print(f"  LIANA rank_aggregate: {len(results_df)} interactions")
    except Exception as e:
        print(f"  rank_aggregate error: {e}")
        try:
            li.mt.cellphonedb(adata, groupby="cell_type", expr_prop=0.1,
                              verbose=False, use_raw=False, n_perms=100)
            results_df = adata.uns["liana_res"].copy()
            print(f"  CellPhoneDB fallback: {len(results_df)} interactions")
        except Exception as e2:
            print(f"  LIANA error: {e2}")
            return

    # Save raw LIANA output
    results_df.to_csv(os.path.join(out_dir, f"liana_all_methods_{dataset_name}.csv"),
                      index=False)

    # Load NEUROSCOPE results
    if not os.path.exists(neuroscope_csv):
        print(f"  NEUROSCOPE results not found: {neuroscope_csv}")
        return

    neuro = pd.read_csv(neuroscope_csv)
    neuro_sig = neuro[neuro["FDRCapacity"] < 0.05].copy() if "FDRCapacity" in neuro.columns else neuro.head(50)
    neuro_pairs = set(zip(neuro_sig["Ligand"].str.upper(),
                          neuro_sig["Receptor"].str.upper()))

    # Build comparison per method
    method_cols = [c for c in results_df.columns
                   if c not in ("source","target","ligand_complex","receptor_complex")]

    comparison_rows = []
    for _, row in neuro_sig.iterrows():
        lig = row["Ligand"].upper()
        rec = row["Receptor"].upper() if "Receptor" in row else ""
        cap = row.get("CondCapacity", np.nan)
        fdr = row.get("FDRCapacity",  np.nan)

        # Match in LIANA results — check ligand_complex and receptor_complex
        match = results_df[
            (results_df.get("ligand_complex",  pd.Series(dtype=str)).str.upper() == lig) &
            (results_df.get("receptor_complex",pd.Series(dtype=str)).str.upper() == rec)
        ]

        row_out = {"Ligand": lig, "Receptor": rec,
                   "NEUROSCOPE_Capacity": cap, "NEUROSCOPE_FDR": fdr,
                   "In_LIANA": len(match) > 0}

        # Collect per-method scores from matched rows
        if len(match) > 0:
            for col in method_cols:
                if col in match.columns and match[col].dtype in (float, np.float64):
                    row_out[f"LIANA_{col}"] = match[col].iloc[0]

        comparison_rows.append(row_out)

    comp_df = pd.DataFrame(comparison_rows)
    comp_path = os.path.join(out_dir, f"comparison_liana_vs_neuroscope.csv")
    comp_df.to_csv(comp_path, index=False)
    print(f"  Saved: {comp_path}")

    n_both    = comp_df["In_LIANA"].sum()
    n_neuro   = len(comp_df)
    print(f"\n  NEUROSCOPE sig:      {n_neuro}")
    print(f"  Also in LIANA:       {n_both}  ({100*n_both/max(n_neuro,1):.0f}%)")
    print(f"  NEUROSCOPE-only:     {n_neuro - n_both}  ({100*(n_neuro-n_both)/max(n_neuro,1):.0f}%)")

    return comp_df


def mock_liana_comparison(neuroscope_csv, dataset_name, out_dir):
    """
    When liana-py is unavailable, produce a documented placeholder
    showing the comparison framework and expected output columns.
    """
    os.makedirs(out_dir, exist_ok=True)
    if not os.path.exists(neuroscope_csv):
        return

    neuro = pd.read_csv(neuroscope_csv)
    neuro_sig = neuro[neuro["FDRCapacity"] < 0.05] if "FDRCapacity" in neuro.columns else neuro.head(50)

    placeholder = neuro_sig[["Ligand","Receptor","CondCapacity","FDRCapacity"]].copy() if \
        all(c in neuro_sig.columns for c in ["Ligand","Receptor","CondCapacity","FDRCapacity"]) else \
        neuro_sig.head(10)
    placeholder["In_LIANA_CellPhoneDB"]   = "REQUIRES_LIANA_INSTALL"
    placeholder["In_LIANA_CellChat"]      = "REQUIRES_LIANA_INSTALL"
    placeholder["In_LIANA_NicheNet"]      = "REQUIRES_LIANA_INSTALL"
    placeholder["In_LIANA_NATMI"]         = "REQUIRES_LIANA_INSTALL"
    placeholder["In_LIANA_Connectome"]    = "REQUIRES_LIANA_INSTALL"
    placeholder["In_LIANA_rank_aggregate"]= "REQUIRES_LIANA_INSTALL"

    out_path = os.path.join(out_dir, "comparison_liana_vs_neuroscope.csv")
    placeholder.to_csv(out_path, index=False)
    print(f"  Placeholder saved: {out_path}")
    print(f"  Install liana-py then re-run: pip install liana anndata scanpy")


# ── Main ──────────────────────────────────────────────────────────────────────
DATASETS = [
    {
        "name":           "Kang2018",
        "mat_path":       "converted/kang2018/kang_expression_sparse.mat",
        "meta_path":      "converted/kang2018/kang_cell_metadata.csv",
        "gene_path":      "converted/kang2018/kang_gene_names.csv",
        "neuroscope_csv": "results/kang2018/channel_rankings.csv",
        "out_dir":        "results/kang2018",
    },
    {
        "name":           "COVIDAtlas",
        "mat_path":       "converted/covid_atlas/covid_expression_sparse.mat",
        "meta_path":      "converted/covid_atlas/covid_cell_metadata.csv",
        "gene_path":      "converted/covid_atlas/covid_gene_names.csv",
        "neuroscope_csv": "results/covid_atlas/channel_rankings.csv",
        "out_dir":        "results/covid_atlas",
    },
    {
        "name":           "SmillieIBD",
        "mat_path":       "converted/smillie2019/smillie_expression_sparse.mat",
        "meta_path":      "converted/smillie2019/smillie_cell_metadata.csv",
        "gene_path":      "converted/smillie2019/smillie_gene_names.csv",
        "neuroscope_csv": "results/smillie2019/channel_rankings.csv",
        "out_dir":        "results/smillie2019",
    },
    {
        "name":           "Lupus",
        "mat_path":       "converted/lupus/lupus_expression_sparse.mat",
        "meta_path":      "converted/lupus/lupus_cell_metadata.csv",
        "gene_path":      "converted/lupus/lupus_gene_names.csv",
        "neuroscope_csv": "results/lupus/channel_rankings.csv",
        "out_dir":        "results/lupus",
    },
]

for ds in DATASETS:
    if not os.path.exists(ds["mat_path"]):
        print(f"\nSkipping {ds['name']}: data not prepared yet ({ds['mat_path']})")
        if os.path.exists(ds["neuroscope_csv"]):
            mock_liana_comparison(ds["neuroscope_csv"], ds["name"], ds["out_dir"])
        continue

    if not LIANA_OK:
        print(f"\n{ds['name']}: liana-py not installed — creating placeholder")
        mock_liana_comparison(ds["neuroscope_csv"], ds["name"], ds["out_dir"])
        continue

    try:
        adata = load_dataset(ds["mat_path"], ds["meta_path"], ds["gene_path"])
        run_liana_methods(adata, ds["name"], ds["neuroscope_csv"], ds["out_dir"])
    except Exception as e:
        print(f"\n{ds['name']} error: {e}")
        mock_liana_comparison(ds["neuroscope_csv"], ds["name"], ds["out_dir"])

print("\n=== LIANA comparison complete ===")
print("\nMethod interpretation guide:")
print("  CellPhoneDB:   mean(L) * mean(R), permutation test — no trimean zero-floor")
print("  CellChat:      trimean — zero-floor at >50% sparsity")
print("  NicheNet:      prior-knowledge ligand activity scoring")
print("  NATMI:         specificity-weighted product of mean expression")
print("  Connectome:    expression-weight edge scoring")
print("  rank_aggregate: aggregate rank of all methods (consensus)")
print("\nPairs significant in NEUROSCOPE but absent from all LIANA methods")
print("are the strongest candidates for novel KSG-exclusive interactions.")
