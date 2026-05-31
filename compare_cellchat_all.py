"""
compare_cellchat_all.py — CellChat comparison for all 4 benchmarks using
LIANA all-methods CSVs (scaled_weight = CellChat trimean score).

CellChat zero-floor: trimean=0 when Q1=Q2=0 (>50% sparsity), manifests as
NaN scaled_weight in LIANA output.
"""
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import os

DATASETS = [
    {
        "name":           "Kang2018",
        "neuroscope_csv": "results/kang2018/channel_rankings.csv",
        "liana_all_csv":  "results/kang2018/liana_all_methods_Kang2018.csv",
        "out_dir":        "results/kang2018",
    },
    {
        "name":           "COVIDAtlas",
        "neuroscope_csv": "results/covid_atlas/channel_rankings.csv",
        "liana_all_csv":  "results/covid_atlas/liana_all_methods_COVIDAtlas.csv",
        "out_dir":        "results/covid_atlas",
    },
    {
        "name":           "Lupus",
        "neuroscope_csv": "results/lupus/channel_rankings.csv",
        "liana_all_csv":  "results/lupus/liana_all_methods_Lupus.csv",
        "out_dir":        "results/lupus",
    },
    {
        "name":           "SmillieIBD",
        "neuroscope_csv": "results/smillie2019/channel_rankings.csv",
        "liana_all_csv":  "results/smillie2019/liana_all_methods_SmillieIBD.csv",
        "out_dir":        "results/smillie2019",
    },
]

print("=== CellChat comparison (all 4 benchmarks) ===")
print("Method: LIANA scaled_weight column (CellChat trimean score)")
print("NaN scaled_weight = zero-floored by CellChat (>50% sparsity)\n")

summary_rows = []

for ds in DATASETS:
    name = ds["name"]
    print(f"\n{'='*55}")
    print(f"Dataset: {name}")

    if not os.path.exists(ds["neuroscope_csv"]):
        print(f"  Missing: {ds['neuroscope_csv']}")
        continue
    if not os.path.exists(ds["liana_all_csv"]):
        print(f"  Missing: {ds['liana_all_csv']}")
        continue

    neuro = pd.read_csv(ds["neuroscope_csv"])
    liana = pd.read_csv(ds["liana_all_csv"], low_memory=False)

    # NEUROSCOPE sig pairs
    fdr_col = "FDRCapacity" if "FDRCapacity" in neuro.columns else "FDR"
    neuro_sig = neuro[pd.to_numeric(neuro[fdr_col], errors="coerce") < 0.05].copy()
    neuro_all = neuro.copy()
    neuro_sig["pair"] = neuro_sig["Ligand"].str.upper() + "_" + neuro_sig["Receptor"].str.upper()
    neuro_all["pair"] = neuro_all["Ligand"].str.upper() + "_" + neuro_all["Receptor"].str.upper()

    # CellChat: scaled_weight per (ligand, receptor) across all cluster-pairs
    liana["pair"] = liana["ligand_complex"].str.upper() + "_" + liana["receptor_complex"].str.upper()
    liana["sw"] = pd.to_numeric(liana["scaled_weight"], errors="coerce")

    # A pair is CellChat-detectable if ANY cluster-pair row has non-NaN scaled_weight
    cc_detected = set(liana[liana["sw"].notna()]["pair"])

    neuro_sig_pairs  = set(neuro_sig["pair"])
    neuro_all_pairs  = set(neuro_all["pair"])
    neuro_not_sig    = neuro_all_pairs - neuro_sig_pairs

    in_cc    = neuro_sig_pairs & cc_detected
    not_cc   = neuro_sig_pairs - cc_detected
    # Zero-floor: pair is in some LIANA rows but all scaled_weight NaN
    all_liana_pairs = set(liana["pair"])
    zero_floored = (neuro_sig_pairs - cc_detected) & all_liana_pairs

    print(f"  NEUROSCOPE FDR-sig: {len(neuro_sig_pairs)}")
    print(f"  In CellChat:        {len(in_cc)}  ({100*len(in_cc)/max(1,len(neuro_sig_pairs)):.0f}%)")
    print(f"  NEUROSCOPE-only:    {len(not_cc)}  ({100*len(not_cc)/max(1,len(neuro_sig_pairs)):.0f}%)")
    print(f"    (zero-floored):   {len(zero_floored)}")

    # Fisher exact: are NEUROSCOPE-sig pairs enriched for CellChat-detectable?
    a = len(neuro_sig_pairs & cc_detected)
    b = len(neuro_sig_pairs - cc_detected)
    c = len(neuro_not_sig & cc_detected)
    d = len(neuro_not_sig - cc_detected)
    if a + b + c + d > 0:
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="less")
        print(f"  Fisher OR (CC vs NEUROSCOPE-sig): {odds:.3f}  p={pval:.4f}")
    else:
        odds, pval = float("nan"), float("nan")

    # Per-pair output
    rows = []
    for _, row in neuro_sig.iterrows():
        pair = row["pair"]
        cap  = pd.to_numeric(row.get("CondCapacity", row.get("capacity", float("nan"))), errors="coerce")
        fdr  = pd.to_numeric(row.get(fdr_col, float("nan")), errors="coerce")
        in_cellchat   = pair in cc_detected
        zero_floor_fl = pair in zero_floored
        rows.append({
            "Ligand":            row["Ligand"],
            "Receptor":          row["Receptor"],
            "NEUROSCOPE_Capacity": cap,
            "NEUROSCOPE_FDR":    fdr,
            "In_CellChat":       in_cellchat,
            "CellChat_ZeroFloor": zero_floor_fl,
        })
    out_df = pd.DataFrame(rows)
    out_path = os.path.join(ds["out_dir"], "comparison_cellchat_vs_neuroscope.csv")
    os.makedirs(ds["out_dir"], exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")

    summary_rows.append({
        "Dataset":         name,
        "NEUROSCOPE_sig":  len(neuro_sig_pairs),
        "In_CellChat":     len(in_cc),
        "NEUROSCOPE_only": len(not_cc),
        "Exclusive_pct":   f"{100*len(not_cc)/max(1,len(neuro_sig_pairs)):.0f}%",
        "Zero_floored":    len(zero_floored),
        "Fisher_OR":       f"{odds:.3f}" if not (isinstance(odds, float) and np.isnan(odds)) else "N/A",
        "Fisher_p":        f"{pval:.4f}" if not (isinstance(pval, float) and np.isnan(pval)) else "N/A",
    })

print("\n\n=== CROSS-DATASET CELLCHAT SUMMARY ===")
summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False))

out_summary = "results/cellchat_comparison_summary.csv"
summary_df.to_csv(out_summary, index=False)
print(f"\nSaved: {out_summary}")
