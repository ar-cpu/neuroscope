# NEUROSCOPE

**An information-theoretic pipeline for measuring cell-cell communication from single-cell RNA-seq data.**

NEUROSCOPE uses KSG mutual information channel capacity to detect which ligand-receptor pairs actually transmit information between cell populations — not merely which pairs are co-expressed. This collapses the majority of co-expressed pairs that carry no measurable signal, leaving only interactions with channel capacity > 0.

> **Core claim:** Every existing tool (CellChat, NicheNet, CellPhoneDB) scores co-expression. NEUROSCOPE scores information transmission. These are not the same thing.

---

## Table of Contents

- [Why NEUROSCOPE](#why-neuroscope)
- [Key Algorithm](#key-algorithm)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Data Format](#input-data-format)
- [Pipeline Modules](#pipeline-modules)
- [Discovery Analyses](#discovery-analyses)
- [Output Files](#output-files)
- [Reproducing Results](#reproducing-results)
- [Citation](#citation)
- [License](#license)

---

## Why NEUROSCOPE

CellChat and similar tools score L-R co-expression using the trimean statistic. For any gene expressed in fewer than 50% of cells, trimean = 0 regardless of expression level. This is a structural property of the trimean formula (Jin et al. 2021, eq. 1): when Q1 = Q2 = 0, the result is always zero, silencing biologically active interactions at realistic single-cell sparsity.

NEUROSCOPE replaces the trimean with **KSG mutual information channel capacity**: I(ligand_distribution; receptor_distribution) estimated via the Kraskov-Stögbauer-Grassberger k-NN estimator on donor-level cluster means. This captures inter-donor variation even at 15% per-cell sparsity.

**Benchmark 1 — Kang 2018 PBMC, GSE96583 (n = 18 donors, IFN-β stimulation):**
- 5/5 NEUROSCOPE top capacity-significant pairs were invisible to CellChat (CellChat trimean = 0 in all cluster combinations)
- Synthetic ground-truth AUROC: KSG = **1.000** vs. co-expression scoring = 0.506
- Top capacity-only pairs: CXCL11→CXCR3, IL1RN→IL1R2, CXCL10→CXCR3 (IFN-inducible chemokines and cytokine regulators — biologically coherent)

**Benchmark 2 — COVID Blood Atlas, Stephenson 2021 (n = 124 donors, 3 conditions):**
- 49 FDR-significant pairs from 127 tested (KSG channel capacity)
- CellChat missed 20/49 (41%) — 18 of those 20 attributable to the trimean zero-floor
- CellPhoneDB missed 5/49 (10%)
- Top pairs: ALOX5→CYSLTR1, PTGES3→PTGER4, CD55→ADGRE5 (eicosanoid and complement pathways)

---

## Key Algorithm

**Channel capacity** (implemented in `functions/estimate_channel_capacity.m`) is computed as:

1. Aggregate single-cell expression to donor-level cluster means (reduces noise, enables kNN)
2. Run KSG k-NN mutual information estimator on paired (ligand_mean, receptor_mean) vectors across donors
3. Normalize by sender entropy: capacity = I(X;Y) / H(X) ∈ [0, 1]
4. Bootstrap 500 resamplings of donor observations to produce 95% BCa confidence intervals

Entropy of donor-level means uses a mixed-distribution decomposition that handles zero-inflation:
H(X) = H_Bernoulli(p₀) + (1 − p₀) × H_continuous(X | X > 0)

Pairs with capacity ≤ 0 after permutation testing are not called as real interactions.

The final `InformationScore` reported in `channel_rankings.csv` combines:

```
InformationScore = 0.60 × CondCapacity
                 + 0.25 × (1 − FDR)
                 + 0.15 × CondFrequency_normalized
```

---

## Installation

**Requirements:**
- MATLAB R2021b or later
- Statistics and Machine Learning Toolbox
- Deep Learning Toolbox (for NMF and autoencoder modules; pipeline falls back to PCA when absent)

**Python requirements** (for data preparation and CellChat comparison):
```bash
pip install anndata scanpy scipy numpy pandas h5py
```

**Setup:**
```matlab
addpath(genpath('functions'));
```

CellPhoneDB interaction tables are bundled in `converted/` (CSV files only — see `.gitignore` for what is excluded).

---

## Quick Start

### Run on the Kang 2018 benchmark dataset

```matlab
% Prepare data first (Python):
% python prep_kang2018.py

% Then run NEUROSCOPE:
run_kang2018
```

### Run on the COVID Blood Atlas (Stephenson 2021, 124 donors)

```matlab
% python prep_covid_atlas.py
run_covid_atlas
```

### Run on your own dataset

```matlab
config = struct();
config.dataset_name        = 'MyDataset';
config.condition_name      = 'Condition';
config.reference_name      = 'Reference';   % leave '' if none
config.expr_path           = 'path/to/expression_sparse.mat';
config.metadata_path       = 'path/to/cell_metadata.csv';
config.output_dir          = 'results/';
config.n_permutations      = 50;
config.fdr_threshold       = 0.05;

NEUROSCOPE(config);
```

### Compare to CellChat (requires R 4.x + CellChat)

```bash
# Install CellChat first:
Rscript install_r_cellchat.R

# Run comparison:
python compare_cellchat.py
```

---

## Input Data Format

### Expression matrix
Sparse MATLAB matrix (cells × genes) saved as `.mat`. Variable name: `expr_matrix`.

### Cell metadata
CSV with one row per cell:

| Column | Required | Description |
|--------|----------|-------------|
| CellID | Yes | Unique cell identifier |
| Cluster | Yes | Integer cluster assignment |
| donor_id | Yes (for causal analysis) | Donor/patient identifier |
| x_coord, y_coord | No | Spatial coordinates |

### CellPhoneDB tables (bundled)
The `converted/` directory contains pre-processed CellPhoneDB v3 CSV tables:
- `interaction_table.csv`
- `protein_table.csv`
- `multidata_table.csv`
- `complex_table.csv`

---

## Pipeline Modules

| Step | Function | What it does |
|------|----------|--------------|
| 1b | `preprocess_expression` | scran-style normalization, HVG selection (3000 genes) |
| 2 | `build_communication_network` | KSG channel capacity scoring of all L-R pairs |
| 2b | `compute_cell_communication_profiles` | Per-cell sending/receiving scores → communication states |
| 3 | `identify_differential_channels` | Wilcoxon on channel strengths + capacity-based permutation test |
| 4b | `rank_communication_channels` | Rank by InformationScore (capacity + significance + frequency) |
| 7 | `identify_communication_hubs` | Hub cluster analysis (betweenness centrality) |
| 8b | `compute_downstream_response` | Transcriptional response in receiver cells |
| 8c | `infer_causal_interactions` | IV/2SLS causal inference using donor_id as instrument |
| 8d | `compute_communication_entropy` | Spectral communities from communication Laplacian |
| 8e | `discover_communication_programs` | NMF tensor decomposition |
| 11 | `generate_science_report` | Plain-text discovery report with bootstrap CIs |

---

## Discovery Analyses

All discovery functions are in `functions/` and can be run independently.

| Function | What it finds |
|----------|---------------|
| `discover_spatial_communication_niches` | Microenvironmental communication domains (Fisher's exact + BH-FDR) |
| `discover_receptor_promiscuity` | Hub receptors with high ligand Shannon entropy |
| `discover_autocrine_signaling` | Self-sustaining co-expressed L-R loops within clusters |
| `discover_feedback_circuits` | Directed cycles (geometric mean circuit strength) |
| `discover_network_motifs` | Enriched 3-node subgraph patterns vs. random expectation |
| `discover_network_fragility` | Percolation curves + pairwise drug synergy scores |
| `discover_communication_monopolies` | Exclusive cluster-to-cluster signaling channels |
| `discover_communication_role_asymmetry` | Broadcaster vs. Listener cluster classification |

---

## Output Files

```
results/{dataset}/
├── {dataset}_science_report.txt            — Top interactions with capacity + 95% BCa CIs
├── channel_rankings.csv                    — All channels ranked by InformationScore
├── {dataset}_specific_interactions.csv     — All differential channels
├── {dataset}_specific_significant_only.csv — FDR < 0.05 subset
├── networks/
│   └── {dataset}_full_network.csv          — Complete scored network
├── novel_discoveries/
│   ├── communication_roles.csv             — Broadcaster/Listener classification
│   ├── downstream_response_genes.csv       — Transcriptional response targets
│   ├── causal_interactions.csv             — IV/2SLS causal estimates
│   ├── entropy_landscape.csv              — Spectral community assignments
│   └── communication_programs.csv          — NMF program loadings
├── statistics/
│   ├── differential_channels.csv           — Wilcoxon + capacity test results
│   └── differential_channels_significant.csv
├── paper_figures/                           — Publication figures
└── poster_figures/
```

### Key columns in `channel_rankings.csv`

| Column | Description |
|--------|-------------|
| Ligand / Receptor | Gene symbols |
| CondCapacity | KSG channel capacity I(L;R)/H(L) in condition, range [0, 1] |
| RefCapacity | Capacity in reference condition |
| FDR / FDRCapacity | BH-adjusted p-values (expression and capacity tests) |
| BootCI_lo / BootCI_hi | 95% BCa bootstrap CI (500 donor resamples) |
| InformationScore | Composite rank score (capacity + FDR + frequency) |
| CausalScore | IV/2SLS causal estimate; NaN when instrument weak |
| Direction | 'Up', 'Down', or 'Both' |
| Pathway | Annotated signalling pathway |

---

## Reproducing Results

### Kang 2018 PBMC (GSE96583)

```bash
# 1. Download and prepare data:
python prep_kang2018.py

# 2. Run NEUROSCOPE (in MATLAB):
run_kang2018

# 3. Run CellChat comparison:
python compare_cellchat.py
```

Expected output: `results/kang2018/Kang2018_science_report.txt`

Key result: 5/5 top capacity-significant pairs have CellChat trimean = 0 in all cluster combinations. AUROC 1.000 vs. 0.506.

### COVID Blood Atlas (Stephenson 2021)

```bash
python prep_covid_atlas.py
# Then in MATLAB:
run_covid_atlas
# Then:
python compare_cellchat_covid.py
```

Expected output: `results/covid_atlas/COVIDAtlas_science_report.txt`

Key result: 49 FDR-significant pairs; CellChat misses 20/49 (41%), 18 attributable to zero-floor. Runtime ~135 min on a 2024 laptop.

---

## Citation

If you use NEUROSCOPE in your research, please cite:

> Ananda A. *NEUROSCOPE: Information-Theoretic Cell-Cell Communication Analysis via KSG Channel Capacity.* 2026.

---

## License

MIT License. See `LICENSE` for details.

---

## Contributing

Issues and pull requests welcome. Please ensure:
- All `fprintf` strings use plain ASCII (no Unicode box-drawing)
- No dataset-specific hardcoding — use `config.dataset_name` throughout
- Statistical methods are documented in function headers
- New functions follow the same struct output convention as existing ones

---

*NEUROSCOPE v5.2 — Information-theoretic cell communication, any tissue, any condition.*
