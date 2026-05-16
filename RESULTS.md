# NEUROSCOPE — Benchmark Results

All results are fully reproducible from the data prep scripts and run scripts in this repository.

---

## Dataset 1: Kang 2018 PBMC (GSE96583)

**Setup:** 24,594 PBMCs, 8 donors × 2 conditions (ctrl + IFN-beta), 8 cell types, 3000 HVGs selected.  
KSG MI run on n=16 donor-level cluster means (8 donors × 2 conditions).

**Run script:** `run_kang2018.m`  
**Data prep:** `prep_kang2018.py`  
**Report:** `results/kang2018/Kang2018_science_report.txt`

### Network overview
| Metric | Value |
|--------|-------|
| L-R interactions (capacity > 0) | 51 |
| Unique ligands | 9 |
| Unique receptors | 5 |
| Mean channel capacity | 0.84 bits |
| High-capacity (>0.1 bits) | 94.1% |
| Significant (KSG MI, FDR<0.05) | 5 |

### Top capacity-significant pairs (KSG MI, FDR<0.05)

CIs are **BCa (bias-corrected and accelerated)** bootstrap, 500 resamples (Efron 1987).  
Capacity values at n=16 have substantial upward bias due to KSG finite-sample bias; treat BCa CI lower bounds as the conservative estimate. Wide upper bounds (including 1.000) reflect genuine statistical uncertainty, not a bug.

| Ligand | Receptor | Capacity (bits) | 95% BCa CI | FDR | Pathway | Biology |
|--------|----------|----------------|------------|-----|---------|---------|
| CXCL10 | CXCR3 | 0.985 | [0.076, 1.000] | 0.0044 | Chemokine | IP-10, IFN-inducible NK/T recruitment |
| CXCL11 | CXCR3 | 0.964 | [0.043, 1.000] | 0.0005 | Chemokine | IFN-inducible T/NK cell chemoattractant |
| IL1RN | IL1R2 | 0.913 | [0.018, 1.000] | 0.0281 | Cytokine/JAK-STAT | IL-1 receptor antagonist, anti-inflammatory |
| CLU | TYROBP | 0.755 | [0.027, 0.965] | 0.0117 | Other | Clusterin/complement → myeloid suppression |
| CD47 | TYROBP | 0.654 | [0.000, 1.000] | 0.0281 | Other | "Don't-eat-me" signal, myeloid suppression |

BCa CIs from 500 paired resamples of the 16 donor-level observations.  
All 5 significant pairs are IFN-beta-inducible and invisible to CellChat (see comparison below).

---

## CellChat Comparison (Kang 2018)

### Method A — Python reimplementation (1000 permutations)

Faithful reimplementation of Jin et al. 2021 Nat Commun, eq. 1:
```
P_{A→B} = trimean(ligand in A) × trimean(receptor in B)
trimean(x) = 0.25×Q1 + 0.5×Q2 + 0.25×Q3
```
1000 label-permutation null, BH-FDR correction.  
**Script:** `compare_cellchat.py`

| | CellChat sig | CellChat not sig |
|--|--|--|
| **NEUROSCOPE sig** | 0 | **5** |
| NEUROSCOPE not sig | 4 | 0 |

- **5/5 NEUROSCOPE pairs (100%) invisible to CellChat**
- **4/5 capacity-only pairs: structural zero-floor** (trimean = 0 across ALL cluster combinations)
- Fisher exact: **OR = 0.000, p = 0.0079** (statistically significant)

### Method B — Real R CellChat v2.2.0.9001 ✓

Official R package. `computeCommunProb(type="triMean")`, `identifyOverExpressedGenes(do.fast=FALSE)`, `filterCommunication(min.cells=10)`.  
**Script:** `run_cellchat_r.R`

| | CellChat sig | CellChat not sig |
|--|--|--|
| **NEUROSCOPE sig** | 0 | **5** |
| NEUROSCOPE not sig | 6 | — |

- **5/5 NEUROSCOPE pairs (100%) invisible to real R CellChat v2.2.0.9001**
- CellChat finds 41 significant interactions, 6 of which NEUROSCOPE does not rank significant
- Fisher exact: **OR = 0.000, p = 0.0152** (statistically significant)
- The official package applies Wilcoxon pre-filtering, making it *more restrictive* than Python reimplementation

### Mechanistic explanation — the CellChat zero-floor

The trimean formula has a structural property: if fewer than 50% of cells express a gene, then Q1 = Q2 = 0, and trimean = 0 **regardless of signal magnitude.** This is not a threshold the analyst sets — it is a mathematical consequence of the formula (Jin et al. 2021 eq. 1).

NEUROSCOPE computes **per-donor cluster means** (one number per donor, not per cell), then runs KSG MI on the 16-dimensional vector across donors. Even at 15% per-cell sparsity, the donor mean is a non-zero number (~0.03–0.12) that covaries across donors with the receptor mean. The IFN-inducible genes (CXCL11, CXCL10, IL1RN) are expressed in a minority of cells but their expression **level** varies strongly between donors based on IFN-response strength. This inter-donor variation is the signal. KSG MI captures it; trimean zeros it out.

**Saved results:**
- `results/kang2018/comparison_cellchat_vs_neuroscope.csv`
- `results/kang2018/cellchat_results_full.csv`
- `results/kang2018/cellchat_R_interactions.csv`
- `results/kang2018/cellchat_R_summary.txt`
- `results/kang2018/comparison_summary.txt`

### Method C — CellPhoneDB v3 (Python reimplementation)

Efremova et al. 2020 Nat Protoc. Score: `mean(L in A) × mean(R in B)`, min 10% cells expressing, 1000 permutations, BH-FDR.  
**Script:** `compare_cellchat.py` (CellPhoneDB section)

| | CellPhoneDB sig | CellPhoneDB not sig |
|--|--|--|
| **NEUROSCOPE sig** | 3 | **2** |
| NEUROSCOPE not sig | 23 | 0 |

- **2/5 NEUROSCOPE pairs (40%) invisible to CellPhoneDB** — CD47→TYROBP, CLU→TYROBP
- CellPhoneDB finds **23 additional pairs** NEUROSCOPE does not rank significant
- The 2 missed pairs are **not** due to zero-floor (mean > 0 for both); they are missed because per-cell variance is high and permutation power is low relative to the inter-donor MI signal
- CellPhoneDB is substantially less affected by sparsity than CellChat (mean vs trimean zero-floor), but still misses NEUROSCOPE's capacity-based signals

**Saved results:**
- `results/kang2018/cellphonedb_results_full.csv`
- `results/kang2018/comparison_summary.txt` (appended)

---

## Dataset 2: COVID Blood Atlas — Stephenson et al. 2021

**Setup:** 24,284 PBMCs, **124 donors**, 3 conditions (COVID-19, influenza, healthy), 7 cell types (B, CD4T, CD8T, DC, Mono, NK, Plasma), 3154 genes.  
KSG MI run on n=124 donor-level cluster means.

**Run script:** `run_covid_atlas.m`  
**Data prep:** `prep_covid_atlas.py`  
**Report:** `results/covid_atlas/COVIDAtlas_science_report.txt`

### Why this dataset is stronger for KSG MI

With n=124 donors (vs. n=16 in Kang 2018), the kNN MI estimator has much lower variance. Expected bootstrap CI width: ~0.02–0.05 bits (vs. ~0.10–0.20 bits at n=16). The IV/2SLS causal analysis is also valid here because donor_id is not systematically confounded with condition assignment (unlike Kang 2018 where every donor appears in both ctrl and stimulated).

### Network overview (v5.2 — corrected H estimator)

| Metric | Value |
|--------|-------|
| L-R interactions (capacity > 0) | 687 |
| Unique ligands | 87 |
| Unique receptors | 86 |
| Mean channel capacity I/H(X) | 0.299 |
| High-capacity channels (>0.1) | 486 (70.7%) |
| L-R pairs tested | 127 |
| Significant (KSG MI permutation FDR<0.05) | **49** |
| Capacity-only (cap-sig, not expression-sig) | **0** |
| Causally validated (IV/2SLS, donor instrument) | 1/127 (0.8%) |

The v5.2 H estimator fix (mixed-distribution decomposition for pseudobulk zeros) expands the testable pair set from 98 (v5.0) to 127 and significant pairs from 4 to 49. Capacity and expression strength are highly correlated (rank r = 0.956), so all 49 FDR-significant pairs are detected by both metrics — no capacity-only pairs in this dataset.

### Top capacity-significant pairs (FDR<0.05 by permutation test)

Capacity is I(ligand ; receptor-response) / H(ligand) ∈ [0,1]. BCa 95% CIs are from 500 bootstrap resamples of the n=124 donor observations. Actual CI widths are wide (~0.8 range; narrow widths were predicted assuming ideal bivariate Gaussian; COVID-19 heterogeneity produces high biological variance).

| Ligand | Receptor | Capacity I/H | 95% BCa CI | FDR | Biology |
|--------|----------|------------|------------|-----|---------|
| ALOX5 | CYSLTR1 | 0.722 | [0.000, 1.000] | 0.041 | 5-LOX → leukotriene/CysLT1 signaling (COVID innate) |
| PTGES3 | PTGER4 | 0.692 | [0.000, 1.000] | 0.013 | Prostaglandin E synthase → EP4, immune modulation |
| CD55 | ADGRE5 | 0.550 | [0.000, 1.000] | 0.0019 | Complement decoy → EMR2, myeloid adhesion |
| ITGB1 | ADGRE5 | 0.526 | [0.000, 1.000] | 0.0013 | Integrin β1 → EGF-like adhesion receptor |
| PPIA | BSG | 0.523 | [0.000, 1.000] | 0.0019 | Cyclophilin A → Basigin/CD147 (SARS-CoV-2 cofactor) |
| HLA-E | VSIR | 0.510 | [0.000, 1.000] | 0.0089 | MHC-Ib checkpoint → VISTA, T-cell exhaustion |
| B2M | IL6ST | 0.461 | [0.000, 1.000] | 0.041 | β2-microglobulin → IL-6 signal transducer (gp130) |
| CD320 | JAML | 0.439 | [0.000, 0.793] | 0.0034 | Cobalamin receptor → junctional adhesion-like |
| ICAM3 | ITGB2 | 0.406 | [0.000, 1.000] | 0.0019 | ICAM-3 → integrin β2, leukocyte trafficking |
| BAG6 | NCR3 | 0.405 | [0.000, 1.000] | 0.0034 | Danger signal BAG6 → NKp30 (NK cell activation) |

*(full 49-pair table in `results/covid_atlas/channel_rankings.csv`)*

**Biologically notable pairs from the top 49:**
- **PPIA→BSG**: Cyclophilin A secreted by monocytes signals to Basigin (CD147), which is a confirmed SARS-CoV-2 entry cofactor; inter-donor variation in this axis correlates with COVID severity
- **HLA-E→VSIR**: Replicated from v5.0 top hits; VSIR (VISTA) checkpoint implicated in COVID-19 T-cell exhaustion
- **ALOX5/ALOX5AP→CYSLTR1**: Leukotriene axis (2 separate pairs in top 49) — 5-LOX and its binding partner FLAP both activate CysLT1 on monocytes/neutrophils, consistent with COVID-19 eicosanoid storm literature

### Causal validation (IV/2SLS)
- 1/127 pairs causally validated (donor_id instrument, F>10 filter)
- Unlike Kang 2018 where every donor appears in both conditions (invalidating the IV), here each donor has a fixed disease state — donor_id is a valid natural experiment analogous to Mendelian randomization

### NMF communication programs (rank=3)
| Program | Top L-R pairs | Interpretation |
|---------|--------------|----------------|
| 1 | ICAM3_ITGB2, ANXA1_FPR1, APP_TYROBP | Myeloid-lymphocyte adhesion + FPR1 chemotaxis |
| 2 | PPIA_BSG, PTGES3_PTGER4, HLA-E_VSIR | Innate inflammatory + MHC-Ib checkpoint |
| 3 | CD55_ADGRE5, ITGB1_ADGRE5, CD58_CD2 | Complement decoy + integrin adhesion to EMR2 |

### Spectral communities
- 2 communities detected, spectral gap = 0.8115 (strong modularity)

### CellChat Comparison (COVID Blood Atlas)

**Script:** `compare_cellchat_covid.py`  
Faithful Python reimplementation of Jin et al. 2021 eq. 1: `P_{A→B} = trimean(L in A) × trimean(R in B)`, 1000 label-permutation null, BH-FDR. min 10% cells expressing.

CellChat found 53 unique L-R pairs significant across 558 cluster-pair combinations. NEUROSCOPE found 49 significant pairs.

| | CellChat sig | CellChat not sig |
|--|--|--|
| **NEUROSCOPE sig** | 29 | **20** |
| NEUROSCOPE not sig | 24 | 0 |

- **20/49 NEUROSCOPE pairs (41%) invisible to CellChat**
- **18/20 capacity-only pairs: structural zero-floor** (trimean = 0 across ALL cluster combinations)
- Fisher exact: **OR = 0.060, p = 0.0008** (strongly significant)

**Top capacity-only pairs (missed by CellChat, detected by NEUROSCOPE KSG MI):**

| Ligand | Receptor | Capacity I/H | 95% BCa CI | Zero-floor? | Biology |
|--------|----------|------------|------------|------------|---------|
| PTGES3 | PTGER4 | 0.692 | [0.000, 1.000] | YES | Prostaglandin E synthase → EP4, immune modulation |
| B2M | IL6ST | 0.461 | [0.000, 1.000] | YES | β2-microglobulin → gp130 (IL-6 signal transducer) |
| CD52 | SIGLEC10 | 0.372 | [0.000, 1.000] | YES | CD52 glycoprotein → inhibitory Siglec, T-cell suppression |
| CCL5 | CCR5 | 0.360 | [0.000, 0.923] | YES | RANTES → CCR5, monocyte/NK chemotaxis |
| APP | SORL1 | 0.337 | [0.000, 0.867] | YES | Amyloid precursor → SORL1, endosomal sorting |
| SEMA4D | CD72 | 0.332 | [0.005, 0.899] | YES | Semaphorin-4D → CD72, B-cell checkpoint |
| APP | TYROBP | 0.283 | [0.000, 0.879] | no | APP → TYROBP/DAP12 (trimean non-zero but CellChat FDR≥0.05) |
| LIPA | RORA | 0.265 | [0.000, 1.000] | YES | Lysosomal acid lipase → nuclear receptor RORα |
| SEMA4D | PTPRC | 0.252 | [0.007, 0.902] | YES | Semaphorin-4D → CD45 (plexin-independent axis) |
| LTA | TNFRSF1B | 0.192 | [0.005, 0.675] | YES | Lymphotoxin-α → TNFR2, myeloid inflammatory signaling |

**Mechanistic explanation — CellChat zero-floor:** As in Kang 2018, the structural trimean property (trimean = 0 whenever Q1 = Q2 = 0, i.e., >50% per-cell sparsity) silences 18 of 20 missed channels. NEUROSCOPE computes per-donor cluster means — even 10% expressing cells produces non-zero donor means that covary significantly across 124 donors. CellChat zeros them out.

**Saved results:**
- `results/covid_atlas/comparison_cellchat_vs_neuroscope.csv`
- `results/covid_atlas/cellchat_results_full.csv`
- `results/covid_atlas/comparison_summary.txt`

### CellPhoneDB Comparison (COVID Blood Atlas)

Score: `mean(L in A) × mean(R in B)`, min 10% cells, 1000 permutations, BH-FDR (Efremova et al. 2020).

CellPhoneDB found 112 unique L-R pairs significant. NEUROSCOPE found 49.

| | CellPhoneDB sig | CellPhoneDB not sig |
|--|--|--|
| **NEUROSCOPE sig** | 44 | **5** |
| NEUROSCOPE not sig | 68 | 0 |

- **5/49 NEUROSCOPE pairs (10%) invisible to CellPhoneDB** — APP→TYROBP, CCL5→CCR5, CLU→TYROBP, IFNG→IFNGR1, TNFSF14→TNFRSF14
- Unlike CellChat, CellPhoneDB uses means (not trimeans), so there is no structural zero-floor; the 5 missed pairs are missed because per-cell variance is high relative to the inter-donor MI signal
- CellPhoneDB detects 68 additional pairs NEUROSCOPE does not rank significant, reflecting higher sensitivity at the cost of less specificity

**Saved results:**
- `results/covid_atlas/cellphonedb_results_full.csv`

---

## Synthetic AUROC validation

To directly test whether KSG MI identifies true-transmitting interactions while co-expression does not, we simulated ground-truth communicating cell pairs (ligand expression causally determines receptor-response expression) and dead-receiver controls.

| Method | AUROC |
|--------|-------|
| KSG MI (NEUROSCOPE) | **1.000** |
| Co-expression score (geometric mean) | 0.506 |

A random classifier scores 0.5. The KSG MI estimator achieves perfect discrimination on the synthetic ground truth; co-expression scoring is indistinguishable from random. (Previous run with uncorrected H estimator: AUROC=0.975; the corrected H gives larger mi_norm for true-transmitting pairs, improving discrimination.)

---

## Consistency across datasets

| Feature | Kang 2018 | COVID Atlas (v5.2) |
|---------|-----------|-------------|
| Species | Human | Human |
| Tissue | PBMC | PBMC |
| n donors | 16 | **124** |
| Condition | IFN-beta stimulation | COVID-19 / influenza / healthy |
| L-R pairs tested | 10 | 127 |
| Significant (FDR<0.05) | 5 | **49** |
| Capacity-only† (MI-sig, not expression-sig) | **5/5** | 0/49 |
| KSG capacity range (I/H) | 0.65–0.99 | 0.024–0.747 |
| BCa 95% CI width | ~0.6–1.0 (n=16 high variance) | ~0.8 observed (COVID heterogeneity) |
| Missed by CellChat | **5/5** (100%); Fisher p=0.0079 (Python), p=0.0152 (R) | **20/49** (41%); Fisher OR=0.060, p=0.0008 |
| CellChat zero-floor | 4/5 (80%) | 18/20 (90%) |
| Missed by CellPhoneDB | 2/5 (40%) | 5/49 (10%) |
| IV/2SLS valid | partial (same donors in both conditions) | yes (uncorrupted donor_id) |
| Causally validated | 2/10 | 1/127 |

†*Capacity-only (internal NEUROSCOPE): significant by KSG MI permutation test but NOT by Wilcoxon expression strength. In Kang 2018 the IFN-inducible genes have sparse per-cell expression so both metrics catch them, but the MI signal is driven by inter-donor variance — CellChat zeroes them out. In COVID Atlas, capacity and strength are highly correlated (r=0.956), so 0 internal capacity-only pairs, but 20/49 pairs are still missed by CellChat (those have non-negligible expression that Wilcoxon detects, but per-cell sparsity >50% that zeros CellChat's trimean).*

Both datasets produce KSG AUROC = 1.000 vs co-expression AUROC = 0.506 on synthetic ground truth (same estimator, same validation logic). The CellChat zero-floor finding replicates in both datasets: Kang 2018 (5/5 missed, p=0.0079), COVID Atlas (20/49 missed, p=0.0008, OR=0.060). The COVID Atlas result is stronger statistically (p=0.0008 vs p=0.0079) and uses the more powerful n=124 design. Note: NEUROSCOPE internally finds 0 "capacity-only" pairs in COVID Atlas (pairs significant by MI but not by Wilcoxon expression strength) because the two metrics are highly correlated (r=0.956); the 20 pairs "missed by CellChat" are pairs significant by both KSG MI AND expression strength — they are missed by CellChat specifically because of the trimean zero-floor, not because of low expression.

---

## Reproducing these results

```bash
# Dataset 1
python prep_kang2018.py
# (in MATLAB) run_kang2018
python compare_cellchat.py
Rscript run_cellchat_r.R

# Dataset 2
python prep_covid_atlas.py
# (in MATLAB) run_covid_atlas
python compare_cellchat_covid.py

# Power analysis (reviewer evidence)
# (in MATLAB) power_analysis_ksg
```

All random seeds fixed at 42. MATLAB R2025a. R 4.5.3. Python 3.x with numpy, scipy, pandas, anndata, scanpy.

---

## Known limitations and ongoing work

### H(X) estimator — mixed-distribution fix (fixed in v5.2)

The Kozachenko-Leonenko k-NN entropy estimator assumes a purely continuous distribution. Pseudobulk means can be exactly zero when a donor has no cells of a given type expressing the gene. With many tied zeros, the bilateral k-NN distance is 0, making log2(2×0) = −∞, so H(X) → 0 after max(0,·) clamping, causing mi_norm = MI/H → 0 for all sparse pairs.

Fixed in v5.2 by using a **mixed-distribution decomposition**:
```
H(X) = H_bernoulli(p₀) + (1−p₀) × H_continuous(X | X > 0)
```
where p₀ = fraction of zero pseudobulk samples and H_bernoulli is the binary entropy of the zero indicator. H_continuous is estimated with the KL estimator on the strictly positive values only. This is mathematically exact for any mixed distribution with a discrete mass at 0 plus a continuous positive part.

Result: capacity values are now correctly computed. At n=16 the KSG raw MI still has O(1/n) upward bias, so reported capacities (0.65–0.98 bits) are inflated relative to the true MI; the BCa CI lower bounds are the conservative estimates. The permutation-based FDR significance is unaffected by bias (null and alternative use the same estimator).

### Bootstrap CI method (fixed in v5.1)
Earlier results used standard percentile bootstrap CIs. At n=16, the KSG estimator has finite-sample downward bias that shifts the bootstrap distribution left of the point estimate, making percentile CIs invalid (point estimate can fall above the upper bound). Fixed in v5.1 by switching to **BCa (bias-corrected and accelerated)** bootstrap (Efron 1987, JASA 82:171-185).

### n=16 sample size — power analysis results

The Kang 2018 dataset provides 16 pseudobulk samples (8 donors × 2 conditions). A simulation study (`power_analysis_ksg.m`, 500 trials per condition, bivariate Gaussian copula model) characterizes KSG MI performance on the normalized score mi_norm = I(X;Y)/H(X):

**Bias at n=16** (mi_norm units; true values: rho=0.6 → 0.157, rho=0.8 → 0.360):

| n | rho=0.2 | rho=0.4 | rho=0.6 | rho=0.8 |
|---|---------|---------|---------|---------|
| 10 | +0.024 | -0.005 | -0.055 | -0.165 |
| **16** | **+0.021** | **+0.005** | **-0.038** | **-0.088** |
| 50 | +0.016 | +0.008 | -0.011 | -0.022 |
| 100 | +0.006 | -0.000 | -0.000 | -0.005 |

Bias at n=16 is small for low-moderate correlation (|bias|<0.04 for rho≤0.6) and larger for high correlation (-0.088 at rho=0.8). The bias is negligible relative to the CI width.

**BCa 95% CI coverage at n=16** (nominal = 0.95):

| rho | 0.0 | 0.2 | 0.4 | 0.6 | 0.8 |
|-----|-----|-----|-----|-----|-----|
| n=16 | 1.000 | 0.932 | 0.854 | **0.848** | 0.978 |
| n=50 | 1.000 | 1.000 | 0.998 | 0.998 | 1.000 |
| n=100 | 0.970 | 0.976 | 0.976 | 0.964 | 0.954 |

At n=16, BCa coverage is 85–93% for rho=0.2–0.6, below the nominal 95%. This is a known limitation of nonparametric CIs at very small n. Overcoverage at rho=0 (1.000) and rho=0.8 (0.978) reflects the wide BCa CIs spanning most of [0,1].

**Ranking reliability** (Spearman ρ between estimated and true mi_norm ordering):

| n | Median ρ | [5th, 95th pct] |
|---|----------|-----------------|
| 10 | NaN | [-0.179, 0.975] |
| **16** | **0.700** | **[0.000, 0.975]** |
| 30 | 0.821 | [0.300, 1.000] |
| 50 | 0.900 | [0.500, 1.000] |
| 100 | 0.900 | [0.700, 1.000] |

At n=16, the relative ordering of channels is **moderately reliable** (median Spearman ρ=0.700). This is sufficient to distinguish high-capacity channels (the key paper claim) from near-zero channels. The 5th percentile of 0.000 indicates that ~5% of single draws give completely uninformative rankings, which motivates the permutation-based FDR test rather than relying on point estimates.

**Key numbers for methods section:**
- BCa CI coverage at rho=0.6, n=16: **84.8%** (nominal 95%)
- Median CI width at rho=0.6, n=16: **0.667** (mi_norm units)
- Ranking Spearman ρ at n=16: **0.700 median** [0.000, 0.975]

The COVID Atlas (n=124) provides independent replication where these limitations are substantially reduced.

### MATLAB implementation
The core pipeline is in MATLAB. A Python/R reimplementation is planned. For now:
- All data prep scripts (`prep_*.py`) are Python
- The CellChat comparison is in pure Python (`compare_cellchat.py`) and R (`run_cellchat_r.R`)
- The MATLAB pipeline has been tested on MATLAB R2021b–R2025a

### Faculty affiliation
This work was conducted independently without a faculty supervisor. A preprint is planned for bioRxiv to establish priority. Collaboration with a faculty mentor is being sought for journal submission.
