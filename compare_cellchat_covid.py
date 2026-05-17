"""
compare_cellchat_covid.py  —  CellChat-equivalent scoring vs NEUROSCOPE channel capacity
                              on Stephenson et al. 2021 COVID Blood Atlas.
CellChat scoring (Jin et al. 2021 Nat Commun eq. 1): trimean = 0.25*Q1 + 0.5*Q2 + 0.25*Q3
Structural zero-floor: trimean = 0 whenever Q1 = Q2 = 0, i.e., >50% of cells are silent.
Outputs: results/covid_atlas/cellchat_results_full.csv, comparison_cellchat_vs_neuroscope.csv, comparison_summary.txt
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.io
import scipy.sparse as sp
import os, sys, time
from collections import defaultdict

EXPR_PATH  = "converted/covid_atlas/covid_expression_sparse.mat"
META_PATH  = "converted/covid_atlas/covid_cell_metadata.csv"
GENES_PATH = "converted/covid_atlas/covid_gene_names.csv"
LR_CSV     = "converted/interaction_table.csv"
PROT_CSV   = "converted/protein_table.csv"
GENE_TBL   = "converted/gene_table.csv"
COMP_CSV   = "converted/complex_composition_table.csv"
NEURO_RANK = "results/covid_atlas/channel_rankings.csv"
OUT_DIR    = "results/covid_atlas"

N_PERM_CC  = 1000
MIN_FRAC   = 0.10
ALPHA      = 0.05

print("=== CellChat vs NEUROSCOPE Comparison (COVID Blood Atlas) ===\n")
print("Reference: Jin et al. 2021 Nat Commun (CellChat); Kraskov et al. 2004 (KSG MI)\n")

print("Loading expression matrix...")
mat = scipy.io.loadmat(EXPR_PATH)
X = mat['X']
if sp.issparse(X):
    X = X.toarray()
X = X.astype(np.float32)
print(f"  {X.shape}  (cells x genes)")

meta = pd.read_csv(META_PATH)
genes_df = pd.read_csv(GENES_PATH)
gene_names = genes_df['gene_name'].str.upper().values
print(f"  Metadata: {meta.shape}")
print(f"  Donors: {meta['donor_id'].nunique()}")
print(f"  Clusters: {meta['Cluster'].nunique()}")

if X.shape[0] != len(meta):
    print(f"  WARNING: X has {X.shape[0]} rows but metadata has {len(meta)} rows — truncating X")
    X = X[:len(meta), :]
    print(f"  Aligned X to {X.shape}")
print()

gene_to_idx = {g: i for i, g in enumerate(gene_names)}

print("Building L-R pair list from CellPhoneDB...")

protein    = pd.read_csv(PROT_CSV)
gt         = pd.read_csv(GENE_TBL)
inter      = pd.read_csv(LR_CSV)
comp_comp  = pd.read_csv(COMP_CSV)

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
    prot_mid_to_gene[mid] = pn_to_hgnc.get(pn, pn.replace('_HUMAN','').replace('_MOUSE',''))

complex_proteins = defaultdict(list)
for _, row in comp_comp.iterrows():
    complex_proteins[row['complex_multidata_id']].append(row['protein_multidata_id'])

matrix_genes       = set(gene_names)
multidata_to_gene  = dict(prot_mid_to_gene)
for cmid, pmids in complex_proteins.items():
    chosen = ''
    for pmid in pmids:
        g = prot_mid_to_gene.get(pmid, '')
        if g in matrix_genes:
            chosen = g; break
        elif not chosen:
            chosen = g
    if chosen:
        multidata_to_gene[cmid] = chosen

lr_pairs = []
for _, row in inter.iterrows():
    id1, id2 = row['multidata_1_id'], row['multidata_2_id']
    if id1 in multidata_to_gene and id2 in multidata_to_gene:
        g1, g2 = multidata_to_gene[id1], multidata_to_gene[id2]
        if g1 in matrix_genes and g2 in matrix_genes:
            lr_pairs.append((g1, g2))

lr_pairs = list(set(lr_pairs))
print(f"  Valid L-R pairs: {len(lr_pairs)}")

lr_pairs_filt = [(L, R) for L, R in lr_pairs
                 if L in gene_to_idx and R in gene_to_idx]
print(f"  In expression matrix: {len(lr_pairs_filt)}")
print()

clusters        = meta['Cluster'].values.astype(int)
donor_ids       = meta['donor_id'].values
unique_clusters = sorted(np.unique(clusters))
unique_donors   = sorted(np.unique(donor_ids))
n_clusters      = len(unique_clusters)
n_donors        = len(unique_donors)

ct_map = (meta.groupby('Cluster')['CellTypeName'].first().to_dict()
          if 'CellTypeName' in meta.columns
          else {c: str(c) for c in unique_clusters})

def trimean(x):
    q1, q2, q3 = np.percentile(x, [25, 50, 75])
    return 0.25*q1 + 0.5*q2 + 0.25*q3

def expr_frac(x):
    return np.mean(x > 0)

print("Precomputing per-cluster gene statistics...")
t0 = time.time()
cluster_trimean = {}
cluster_frac    = {}
cluster_mean    = {}
for c in unique_clusters:
    mask            = clusters == c
    Xc              = X[mask, :]
    cluster_trimean[c] = np.array([trimean(Xc[:, g]) for g in range(Xc.shape[1])])
    cluster_frac[c]    = np.mean(Xc > 0, axis=0)
    cluster_mean[c]    = np.mean(Xc, axis=0)
print(f"  Done in {time.time()-t0:.1f}s")

rng     = np.random.default_rng(42)
results = []
tested  = 0

print(f"\nRunning CellChat significance test ({N_PERM_CC} permutations, Jin et al. 2021)...")
print(f"  Testing {len(lr_pairs_filt)} L-R pairs across {n_clusters*(n_clusters-1)} directed cluster pairs...")

for L, R in lr_pairs_filt:
    li = gene_to_idx[L]
    ri = gene_to_idx[R]

    for src in unique_clusters:
        l_frac = cluster_frac[src][li]
        if l_frac < MIN_FRAC:
            continue
        l_tm = cluster_trimean[src][li]
        if l_tm <= 0:
            continue

        for tgt in unique_clusters:
            if src == tgt:
                continue
            r_frac = cluster_frac[tgt][ri]
            if r_frac < MIN_FRAC:
                continue
            r_tm = cluster_trimean[tgt][ri]
            if r_tm <= 0:
                continue

            score = l_tm * r_tm

            null_scores = np.zeros(N_PERM_CC)
            for b in range(N_PERM_CC):
                perm        = rng.permutation(len(clusters))
                perm_cl     = clusters[perm]
                src_mask    = perm_cl == src
                tgt_mask    = perm_cl == tgt
                if src_mask.sum() < 5 or tgt_mask.sum() < 5:
                    null_scores[b] = 0.0
                    continue
                null_scores[b] = trimean(X[src_mask, li]) * trimean(X[tgt_mask, ri])

            pval = np.mean(null_scores >= score)
            results.append({
                'Ligand': L, 'Receptor': R,
                'SenderCluster': src, 'ReceiverCluster': tgt,
                'SenderType':   ct_map.get(src, str(src)),
                'ReceiverType': ct_map.get(tgt, str(tgt)),
                'CellChat_Score':   score,
                'CellChat_Pval':    pval,
                'L_Trimean': l_tm, 'R_Trimean': r_tm,
                'L_Frac': l_frac,  'R_Frac':    r_frac,
            })
            tested += 1

print(f"  Tested {tested} combinations in {time.time()-t0:.1f}s")

cc_df = pd.DataFrame(results)
print(f"  CellChat results: {len(cc_df)} combinations")

if len(cc_df) > 0:
    pvals = cc_df['CellChat_Pval'].values
    n     = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(n); ranks[order] = np.arange(1, n+1)
    fdr   = np.minimum(1, pvals * n / ranks)
    sorted_fdr = fdr[order]
    for i in range(n-2, -1, -1):
        sorted_fdr[i] = min(sorted_fdr[i], sorted_fdr[i+1])
    fdr = sorted_fdr[np.argsort(order)]
    cc_df['CellChat_FDR']         = fdr
    cc_df['CellChat_Significant'] = fdr < ALPHA

    cc_sig = cc_df[cc_df['CellChat_Significant']]
    print(f"  CellChat significant (FDR<{ALPHA}): {len(cc_sig)} combinations, "
          f"{cc_sig[['Ligand','Receptor']].drop_duplicates().shape[0]} unique L-R pairs")

print("\nLoading NEUROSCOPE results...")
if not os.path.exists(NEURO_RANK):
    print(f"  ERROR: {NEURO_RANK} not found — run run_covid_atlas.m first")
    sys.exit(1)

neuro_df = pd.read_csv(NEURO_RANK)
print(f"  Loaded: {len(neuro_df)} rows")

neuro_sig_col = None
for col in ['FDRCapacity', 'FDR_cap', 'FDR', 'fdr']:
    if col in neuro_df.columns:
        neuro_sig_col = col; break

neuro_lr_sig = set()
if neuro_sig_col:
    neuro_sig = neuro_df[neuro_df[neuro_sig_col] < ALPHA]
    if 'Ligand' in neuro_df.columns and 'Receptor' in neuro_df.columns:
        neuro_lr_sig = set(zip(neuro_sig['Ligand'].str.upper(), neuro_sig['Receptor'].str.upper()))
    print(f"  NEUROSCOPE significant pairs (FDR<{ALPHA}): {len(neuro_lr_sig)}")

cc_lr_sig = set()
if 'CellChat_Significant' in cc_df.columns:
    cc_lr_sig = set(
        zip(cc_df[cc_df['CellChat_Significant']]['Ligand'],
            cc_df[cc_df['CellChat_Significant']]['Receptor'])
    )

all_lr = neuro_lr_sig | cc_lr_sig
print(f"\nTotal unique L-R pairs in either method: {len(all_lr)}")

def max_trimean_for_lr(L, R):
    """Return max trimean score and fraction-expressing for ligand + receptor."""
    if L not in gene_to_idx or R not in gene_to_idx:
        return 0.0, 0.0, 0.0, 0.0
    li, ri = gene_to_idx[L], gene_to_idx[R]
    best_score  = 0.0
    best_l_frac = 0.0
    best_r_frac = 0.0
    best_l_tm   = 0.0
    for src in unique_clusters:
        l_tm   = cluster_trimean[src][li]
        l_frac = cluster_frac[src][li]
        for tgt in unique_clusters:
            if src == tgt: continue
            r_tm   = cluster_trimean[tgt][ri]
            r_frac = cluster_frac[tgt][ri]
            score  = l_tm * r_tm
            if score > best_score:
                best_score  = score
                best_l_frac = l_frac
                best_r_frac = r_frac
                best_l_tm   = l_tm
    return best_score, best_l_frac, best_r_frac, best_l_tm

print("\nBuilding comparison table...")
comparison = []
for (L, R) in all_lr:
    in_neuro = (L, R) in neuro_lr_sig
    in_cc    = (L, R) in cc_lr_sig

    if in_neuro and in_cc:
        category = 'Both'
    elif in_neuro and not in_cc:
        category = 'Capacity_Only'
    elif in_cc and not in_neuro:
        category = 'Strength_Only'
    else:
        category = 'Neither'

    cap     = np.nan
    ci_lo   = np.nan
    ci_hi   = np.nan
    fdr_cap = np.nan
    if 'Ligand' in neuro_df.columns:
        row = neuro_df[(neuro_df['Ligand'].str.upper()==L) & (neuro_df['Receptor'].str.upper()==R)]
        if len(row) > 0:
            for cap_col in ['CondCapacity', 'Capacity', 'capacity']:
                if cap_col in neuro_df.columns:
                    cap = row[cap_col].values[0]; break
            if 'BootCI_lo' in neuro_df.columns:
                ci_lo = row['BootCI_lo'].values[0]
                ci_hi = row['BootCI_hi'].values[0]
            if neuro_sig_col in neuro_df.columns:
                fdr_cap = row[neuro_sig_col].values[0]

    cc_row      = cc_df[(cc_df['Ligand']==L) & (cc_df['Receptor']==R)]
    cc_best     = cc_row['CellChat_Score'].max() if len(cc_row) > 0 else np.nan
    cc_best_fdr = (cc_row['CellChat_FDR'].min()
                   if ('CellChat_FDR' in cc_row.columns and len(cc_row) > 0) else np.nan)

    max_tm_score, l_frac, r_frac, l_tm = max_trimean_for_lr(L, R)
    # A pair hits the CellChat zero-floor when max trimean product = 0
    # (one or both genes have trimean=0 in all cluster combinations tested)
    zero_floor = (max_tm_score == 0.0) and in_neuro and not in_cc

    comparison.append({
        'Ligand':              L,
        'Receptor':            R,
        'Category':            category,
        'NEUROSCOPE_Capacity': cap,
        'NEUROSCOPE_CI_lo':    ci_lo,
        'NEUROSCOPE_CI_hi':    ci_hi,
        'NEUROSCOPE_FDR':      fdr_cap,
        'CellChat_BestScore':  cc_best,
        'CellChat_BestFDR':    cc_best_fdr,
        'MaxTrimeanProduct':   max_tm_score,
        'In_NEUROSCOPE':       in_neuro,
        'In_CellChat':         in_cc,
        'CellChat_ZeroFloor':  zero_floor,
    })

comp_df = pd.DataFrame(comparison).sort_values('Category')

n_both     = (comp_df['Category'] == 'Both').sum()
n_cap_only = (comp_df['Category'] == 'Capacity_Only').sum()
n_str_only = (comp_df['Category'] == 'Strength_Only').sum()
n_neither  = (comp_df['Category'] == 'Neither').sum()

contingency = np.array([[n_both, n_cap_only],
                         [n_str_only, max(1, n_neither)]])
odds_ratio, fisher_pval = stats.fisher_exact(contingency, alternative='two-sided')

print(f"\n{'='*65}")
print("COMPARISON RESULTS — 2x2 CONTINGENCY (COVID Atlas)")
print(f"{'='*65}")
print(f"                         CellChat sig   CellChat not sig")
print(f"  NEUROSCOPE sig (cap):  {n_both:12d}   {n_cap_only:12d}")
print(f"  NEUROSCOPE ns:         {n_str_only:12d}   {n_neither:12d}")
print(f"{'='*65}")
print(f"  Fisher exact test: OR = {odds_ratio:.3f},  p = {fisher_pval:.4f}")
print()
print(f"  Both methods:     {n_both:4d}  (validated by both)")
print(f"  Capacity-only:    {n_cap_only:4d}  [NOVEL — missed by CellChat]")
print(f"  Strength-only:    {n_str_only:4d}  (CellChat finds, NEUROSCOPE misses)")
print(f"  Neither:          {n_neither:4d}")
print()

n_zero_floor = comp_df['CellChat_ZeroFloor'].sum()
print(f"  CellChat zero-floor pairs: {n_zero_floor}/{n_cap_only} capacity-only pairs")
print(f"  (trimean = 0 in all cluster combinations — structural, not threshold-dependent)")

if n_cap_only > 0:
    cap_only_df = comp_df[comp_df['Category']=='Capacity_Only'].sort_values(
        'NEUROSCOPE_Capacity', ascending=False)
    print("\nTop capacity-only pairs (missed by CellChat, detected by NEUROSCOPE KSG MI):")
    has_ci = not cap_only_df['NEUROSCOPE_CI_lo'].isna().all()
    for _, row in cap_only_df.head(10).iterrows():
        ci_str = (f"  95%CI=[{row['NEUROSCOPE_CI_lo']:.4f},{row['NEUROSCOPE_CI_hi']:.4f}]"
                  if has_ci and not np.isnan(row['NEUROSCOPE_CI_lo'])
                  else "")
        floor_str = "  [ZERO-FLOOR]" if row['CellChat_ZeroFloor'] else f"  trimean={row['MaxTrimeanProduct']:.4f}"
        print(f"  {row['Ligand']:12s} -> {row['Receptor']:12s}  "
              f"cap={row['NEUROSCOPE_Capacity']:.4f}{ci_str}{floor_str}")

os.makedirs(OUT_DIR, exist_ok=True)
comp_df.to_csv(f"{OUT_DIR}/comparison_cellchat_vs_neuroscope.csv", index=False)
cc_df.to_csv(f"{OUT_DIR}/cellchat_results_full.csv",               index=False)

with open(f"{OUT_DIR}/comparison_summary.txt", 'w') as f:
    f.write("CellChat vs NEUROSCOPE Comparison — Stephenson et al. 2021 COVID Blood Atlas\n")
    f.write("=" * 70 + "\n\n")
    f.write(f"Dataset: {len(meta)} cells, {n_donors} donors, {n_clusters} cell types\n")
    f.write(f"L-R pairs tested: {len(lr_pairs_filt)}\n")
    f.write(f"CellChat permutations: {N_PERM_CC} (Jin et al. 2021 Nat Commun)\n")
    f.write(f"NEUROSCOPE: KSG MI, n={n_donors} pseudobulk samples, BCa bootstrap CIs\n\n")

    f.write("2x2 CONTINGENCY TABLE:\n")
    f.write("                         CellChat sig   CellChat not sig\n")
    f.write(f"  NEUROSCOPE sig (cap):  {n_both:12d}   {n_cap_only:12d}\n")
    f.write(f"  NEUROSCOPE ns:         {n_str_only:12d}   {n_neither:12d}\n\n")
    f.write(f"Fisher exact test: OR = {odds_ratio:.3f},  p = {fisher_pval:.4f}\n\n")
    f.write(f"  Both methods:          {n_both}\n")
    f.write(f"  Capacity-only (NOVEL): {n_cap_only}\n")
    f.write(f"  Strength-only:         {n_str_only}\n")
    f.write(f"  Neither:               {n_neither}\n\n")

    f.write("MECHANISTIC EXPLANATION — CellChat ZERO-FLOOR:\n")
    f.write(f"  {n_zero_floor}/{n_cap_only} capacity-only pairs have CellChat trimean = 0.0000\n")
    f.write("  for ALL cluster combinations tested.\n\n")
    f.write("  Structural property of trimean (Jin et al. 2021 eq. 1):\n")
    f.write("    trimean(x) = 0.25*Q1 + 0.5*Q2 + 0.25*Q3\n")
    f.write("    => trimean = 0 whenever Q1 = Q2 = 0 (i.e., >50% sparsity)\n\n")

    if n_cap_only > 0:
        cap_only_df_s = comp_df[comp_df['Category']=='Capacity_Only'].sort_values(
            'NEUROSCOPE_Capacity', ascending=False)
        f.write("TOP CAPACITY-ONLY PAIRS:\n")
        f.write(f"  {'Ligand':12s} {'Receptor':12s} {'Capacity':>10s} "
                f"{'95%CI_lo':>10s} {'95%CI_hi':>10s} {'ZeroFloor':>10s}\n")
        f.write("  " + "-"*70 + "\n")
        for _, row in cap_only_df_s.iterrows():
            ci_lo = f"{row['NEUROSCOPE_CI_lo']:.4f}" if not np.isnan(row['NEUROSCOPE_CI_lo']) else "  N/A"
            ci_hi = f"{row['NEUROSCOPE_CI_hi']:.4f}" if not np.isnan(row['NEUROSCOPE_CI_hi']) else "  N/A"
            f.write(f"  {row['Ligand']:12s} {row['Receptor']:12s} "
                    f"{row['NEUROSCOPE_Capacity']:10.4f} {ci_lo:>10s} {ci_hi:>10s} "
                    f"{'YES' if row['CellChat_ZeroFloor'] else 'no':>10s}\n")

    f.write(f"\nSaved: {OUT_DIR}/comparison_cellchat_vs_neuroscope.csv\n")
    f.write(f"Saved: {OUT_DIR}/cellchat_results_full.csv\n")

print(f"\nSaved: {OUT_DIR}/comparison_cellchat_vs_neuroscope.csv")
print(f"Saved: {OUT_DIR}/comparison_summary.txt")
print(f"Saved: {OUT_DIR}/cellchat_results_full.csv")

print("\n\n=== CellPhoneDB-style Comparison (COVID Atlas) ===")
print("Reference: Efremova et al. 2020 Nat Protoc (CellPhoneDB v3)")
print("Score: mean(L in A) * mean(R in B), minimum 10% cells expressing\n")

MIN_FRAC_CPD = 0.10
N_PERM_CPD   = 1000

t1 = time.time()
cpd_results = []
tested_cpd  = 0

print(f"Running CellPhoneDB permutation test ({N_PERM_CPD} permutations)...")
for L, R in lr_pairs_filt:
    li = gene_to_idx[L]
    ri = gene_to_idx[R]

    for src in unique_clusters:
        l_frac = cluster_frac[src][li]
        if l_frac < MIN_FRAC_CPD:
            continue
        l_mean = cluster_mean[src][li]

        for tgt in unique_clusters:
            if src == tgt:
                continue
            r_frac = cluster_frac[tgt][ri]
            if r_frac < MIN_FRAC_CPD:
                continue
            r_mean = cluster_mean[tgt][ri]
            score  = l_mean * r_mean
            if score <= 0:
                continue

            null_scores = np.zeros(N_PERM_CPD)
            for b in range(N_PERM_CPD):
                perm     = rng.permutation(len(clusters))
                perm_cl  = clusters[perm]
                src_mask = perm_cl == src
                tgt_mask = perm_cl == tgt
                if src_mask.sum() < 5 or tgt_mask.sum() < 5:
                    null_scores[b] = 0.0
                    continue
                null_scores[b] = np.mean(X[src_mask, li]) * np.mean(X[tgt_mask, ri])

            pval = np.mean(null_scores >= score)
            cpd_results.append({
                'Ligand': L, 'Receptor': R,
                'SenderCluster': src, 'ReceiverCluster': tgt,
                'CPD_Score': score, 'CPD_Pval':  pval,
                'L_Mean': l_mean,   'R_Mean': r_mean,
                'L_Frac': l_frac,   'R_Frac': r_frac,
            })
            tested_cpd += 1

print(f"  Tested {tested_cpd} combinations in {time.time()-t1:.1f}s")

cpd_df = pd.DataFrame(cpd_results)
cpd_sig_lr = set()
if len(cpd_df) > 0:
    pvals_c = cpd_df['CPD_Pval'].values
    n_c     = len(pvals_c)
    order_c = np.argsort(pvals_c)
    ranks_c = np.empty(n_c); ranks_c[order_c] = np.arange(1, n_c+1)
    fdr_c   = np.minimum(1, pvals_c * n_c / ranks_c)
    sf = fdr_c[order_c]
    for i in range(n_c-2, -1, -1):
        sf[i] = min(sf[i], sf[i+1])
    fdr_c = sf[np.argsort(order_c)]
    cpd_df['CPD_FDR'] = fdr_c
    cpd_df['CPD_Sig'] = fdr_c < ALPHA
    cpd_sig_lr = set(zip(cpd_df[cpd_df['CPD_Sig']]['Ligand'],
                         cpd_df[cpd_df['CPD_Sig']]['Receptor']))
    print(f"  CellPhoneDB significant (FDR<{ALPHA}): "
          f"{len(cpd_df[cpd_df['CPD_Sig']])} combinations, "
          f"{len(cpd_sig_lr)} unique L-R pairs")

n_both_cpd     = len(neuro_lr_sig & cpd_sig_lr)
n_cap_only_cpd = len(neuro_lr_sig - cpd_sig_lr)
n_str_only_cpd = len(cpd_sig_lr - neuro_lr_sig)

print(f"\n{'='*65}")
print("CELLPHONEDB COMPARISON RESULTS (COVID Atlas)")
print(f"{'='*65}")
print(f"  Both methods:     {n_both_cpd:4d}")
print(f"  Capacity-only:    {n_cap_only_cpd:4d}  [NOVEL — missed by CellPhoneDB]")
print(f"  Strength-only:    {n_str_only_cpd:4d}  (CellPhoneDB finds, NEUROSCOPE misses)")
print()

if n_cap_only_cpd > 0:
    cpd_only_pairs = neuro_lr_sig - cpd_sig_lr
    print("Capacity-only pairs vs CellPhoneDB:")
    for L, R in sorted(cpd_only_pairs):
        if L not in gene_to_idx or R not in gene_to_idx:
            print(f"  {L:12s} -> {R:12s}  [not in expression matrix]")
            continue
        li, ri = gene_to_idx[L], gene_to_idx[R]
        max_mean_score = max(
            (cluster_mean[s][li] * cluster_mean[t][ri])
            for s in unique_clusters for t in unique_clusters if s != t
        )
        is_zero = (max_mean_score == 0.0)
        row = neuro_df[(neuro_df['Ligand'].str.upper()==L) & (neuro_df['Receptor'].str.upper()==R)]
        cap = row['CondCapacity'].values[0] if len(row)>0 and 'CondCapacity' in neuro_df.columns else np.nan
        print(f"  {L:12s} -> {R:12s}  cap={cap:.4f}  "
              f"{'[MEAN-ZERO]' if is_zero else f'max_mean={max_mean_score:.5f}'}")

cpd_df.to_csv(f"{OUT_DIR}/cellphonedb_results_full.csv", index=False)

with open(f"{OUT_DIR}/comparison_summary.txt", 'a') as f:
    f.write("\n" + "=" * 70 + "\n")
    f.write("CELLPHONEDB COMPARISON\n")
    f.write("=" * 70 + "\n")
    f.write("Method: mean(L in A) * mean(R in B), min 10% cells, 1000 permutations\n")
    f.write("Reference: Efremova et al. 2020 Nat Protoc\n\n")
    f.write(f"  Both methods:          {n_both_cpd}\n")
    f.write(f"  Capacity-only (NOVEL): {n_cap_only_cpd}  (missed by CellPhoneDB)\n")
    f.write(f"  Strength-only:         {n_str_only_cpd}\n\n")
    f.write(f"\nSaved: {OUT_DIR}/cellphonedb_results_full.csv\n")

print(f"\nSaved: {OUT_DIR}/cellphonedb_results_full.csv")
print(f"\nTotal runtime: {time.time()-t0:.1f}s")
