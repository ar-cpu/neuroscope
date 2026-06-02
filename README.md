# NEUROSCOPE

**Information-theoretic cell-cell communication analysis for single-cell RNA-seq data.**

NEUROSCOPE measures ligand-receptor communication using KSG mutual information channel capacity rather than co-expression scoring. The distinction matters: co-expression methods score whether a ligand and receptor are expressed in neighboring clusters; channel capacity scores whether ligand variation across donors actually predicts receptor variation — a much stricter test that filters out interactions with no measurable signal.

---

## Key Results

Benchmarked across four independent datasets spanning three diseases and two tissue types:

| Dataset | Donors | FDR-sig pairs | Not in LIANA | Not in CellChat | AUROC |
|---------|--------|---------------|--------------|-----------------|-------|
| Kang 2018 — IFN-β PBMC | 16 | 5 / 9 | 3/5 (60%) | 3/5 (60%) | **1.000** vs 0.506 |
| COVID Blood Atlas — PBMC | 124 | 49 / 127 | 32/49 (65%) | 20/49 (41%) | **1.000** vs 0.506 |
| Perez 2022 — SLE PBMC | 149 | 33 / 40 | 19/33 (58%) | 19/33 (58%) | **1.000** vs 0.506 |
| Smillie 2019 — Crohn colon | 18 | 17 / 34 | 12/17 (71%) | 12/17 (70%) | **1.000** vs 0.506 |

AUROC is from a synthetic ground-truth experiment (true-signal vs dead-receiver pairs, seed 42). Co-expression scoring achieves 0.506 ≈ chance across all four datasets.

Full results: [`results/NEUROSCOPE_v52_figures.pdf`](results/NEUROSCOPE_v52_figures.pdf) and [`results/RUN_LOG.md`](results/RUN_LOG.md).

---

## How It Works

Standard tools (CellChat, CellPhoneDB, LIANA) score `mean(ligand) × mean(receptor)` or its variants within a sample. When a gene is expressed in fewer than 50% of cells, CellChat's trimean collapses to zero regardless of expression level — a structural property of the formula that silences biologically active interactions in sparse data.

NEUROSCOPE instead computes per-donor cluster means (one value per donor per cluster), then estimates the KSG k-NN mutual information between paired ligand and receptor vectors across donors. This preserves inter-donor variation that is real biological signal, even when within-sample sparsity is high.

Channel capacity = I(ligand; receptor) / H(ligand), normalized to [0, 1]. Significance is assessed by permutation test (199 shuffles) with BH-FDR correction. Confidence intervals use BCa bootstrap (500 donor resamples).

---

## Installation

**MATLAB** (R2021b or later):
- Statistics and Machine Learning Toolbox
- Deep Learning Toolbox (optional — used for NMF; falls back to PCA if absent)

**Python** (data preparation and comparisons):
```bash
pip install anndata scanpy scipy numpy pandas h5py fsspec
```

**R** (NicheNet comparison only):
```bash
Rscript scripts/install_nichenetr.R
```

No other setup needed. CellPhoneDB reference tables are bundled in `converted/`.

---

## Quick Start

**Run on one of the benchmark datasets:**
```bash
# 1. Download and prepare expression data
python scripts/prep_kang2018.py

# 2. Run the pipeline (from repo root)
matlab -batch "run('runs/kang2018.m')"
```

Results appear in `results/kang2018/`.

**Run on your own data:**
```matlab
config = struct();
config.dataset_name  = 'MyStudy';
config.expr_path     = 'path/to/expression_sparse.mat';   % cells × genes
config.metadata_path = 'path/to/cell_metadata.csv';        % CellID, Cluster, donor_id
config.cellphonedb_dir = 'converted';
config.output_dir    = 'results/mystudy';
config.n_permutations = 200;
config.fdr_threshold  = 0.05;
config.run_ml         = false;   % set true for NMF + causal inference
config.random_seed    = 42;

NEUROSCOPE(config);
```

---

## Input Format

**Expression matrix** — MATLAB sparse matrix (cells × genes), saved as `.mat` with variable name `expr_matrix`. Prep scripts generate this from public H5AD files.

**Cell metadata** — CSV with one row per cell:

| Column | Required | Description |
|--------|----------|-------------|
| `CellID` | Yes | Unique integer cell index |
| `Cluster` | Yes | Integer cluster assignment (1-based) |
| `donor_id` | Yes | Donor / patient identifier (enables KSG cross-sample mode) |

---

## Reproducing the Benchmarks

All four benchmarks can be reproduced from scratch. Data is downloaded from public repositories (GEO, CellxGene).

```bash
# Kang 2018 IFN-beta PBMC
python scripts/prep_kang2018.py
matlab -batch "run('runs/kang2018.m')"

# COVID Blood Atlas (Stephenson 2021)
python scripts/prep_covid_atlas.py
matlab -batch "run('runs/covid_atlas.m')"

# Lupus SLE PBMC (Perez 2022)
python scripts/prep_lupus.py
matlab -batch "run('runs/lupus.m')"

# Crohn's disease colon (Smillie 2019)
python scripts/prep_smillie2019.py
matlab -batch "run('runs/smillie2019.m')"
```

**Comparison against existing tools:**
```bash
# LIANA multi-method (CellPhoneDB, CellChat, NATMI, Connectome)
pip install liana
python scripts/compare_liana.py

# NicheNet LR-network comparison
Rscript scripts/compare_nichenet.R
```

Expected runtimes on a standard laptop: Kang 3 min, Smillie 6 min, Lupus ~17 hr (149 donors × 40 pairs × 200 permutations), COVID Atlas ~21 hr.

---

## Repository Structure

```
NEUROSCOPE.m          # Main pipeline entry point
functions/            # All MATLAB pipeline modules
runs/                 # Per-dataset configs (edit these to change parameters)
  kang2018.m
  covid_atlas.m
  lupus.m
  smillie2019.m
  power_analysis_ksg.m
scripts/              # Data prep and comparison scripts
  prep_*.py           # Download + convert each dataset
  compare_*.py / .R   # Benchmark comparisons
  make_figure_doc.py  # Regenerate results/NEUROSCOPE_v52_figures.pdf
converted/            # CellPhoneDB reference tables (bundled)
results/              # Benchmark outputs — one subfolder per dataset
docs/                 # Cleanup notes and file inventory
```

---

## Output Files

Each dataset produces a subfolder under `results/` containing:

| File | Description |
|------|-------------|
| `channel_rankings.csv` | All L-R pairs ranked by channel capacity |
| `*_specific_interactions.csv` | Differential channels (condition vs reference) |
| `*_specific_significant_only.csv` | FDR < 0.05 subset |
| `*_target_rankings.csv` | Therapeutic target rankings |
| `*_science_report.txt` | Narrative report with top interactions and CIs |
| `comparison_liana_vs_neuroscope.csv` | LIANA overlap analysis |
| `comparison_cellchat_vs_neuroscope.csv` | CellChat overlap analysis |
| `comparison_nichenet_vs_neuroscope.csv` | NicheNet LR-network membership |

Key columns in `channel_rankings.csv`:

| Column | Description |
|--------|-------------|
| `CondCapacity` | KSG channel capacity I(L;R)/H(L) ∈ [0, 1] |
| `FDRCapacity` | BH-adjusted permutation p-value |
| `BootCI_lo/hi` | 95% BCa bootstrap CI (500 donor resamples) |
| `InformationScore` | Composite rank score |
| `CausalScore` | IV/2SLS causal estimate (NaN when instrument weak) |

---

## Citation

> Ananda A. *NEUROSCOPE: Information-Theoretic Cell-Cell Communication Analysis via KSG Channel Capacity.* 2026.

---

## License

MIT. See [`LICENSE`](LICENSE).
