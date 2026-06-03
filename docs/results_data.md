# `results/` Cleanup Manifest

**Date:** 2026-05-31
Files below were **untracked** from git (`git rm --cached`) to keep the public repo
small. **They still exist on your local disk and in git history.** See the top-level
[CLEANUP_MANIFEST.md](../CLEANUP_MANIFEST.md) for restore instructions.

Two categories were removed: (1) run **logs** (verbose stdout, not publication
material), and (2) **`liana_all_methods_*.csv`** raw per-cluster-pair dumps (the
summarized `comparison_*_vs_neuroscope.csv` files are the real, kept results).

---

## Files UNTRACKED from git (kept on disk)

### Raw LIANA intermediate dumps (superseded by comparison summaries)

| File | Size | What it is | Regenerate |
|------|------|-----------|-----------|
| `lupus/liana_all_methods_Lupus.csv` | **84 MB** | Per-donor × per-cluster-pair LIANA output, all 5 methods, every L-R pair, Lupus. | `python scripts/compare_liana.py --dataset Lupus` |
| `smillie2019/liana_all_methods_SmillieIBD.csv` | 5.6 MB | Same, Smillie IBD. | `python scripts/compare_liana.py --dataset SmillieIBD` |
| `kang2018/liana_all_methods_Kang2018.csv` | 5.3 MB | Same, Kang 2018. | `python scripts/compare_liana.py --dataset Kang2018` |
| `covid_atlas/liana_all_methods_COVIDAtlas.csv` | 200 KB | Same, COVID Atlas. | `python scripts/compare_liana.py --dataset COVIDAtlas` |

> The **kept** files `comparison_liana_vs_neuroscope.csv`,
> `comparison_cellchat_vs_neuroscope.csv`, and `comparison_nichenet_vs_neuroscope.csv`
> in each dataset folder are the per-pair summaries derived from these dumps. CellChat
> scores were extracted from the `scaled_weight` column of these dumps; that data is
> preserved in the kept `comparison_cellchat_vs_neuroscope.csv` files.

### Run logs (verbose stdout — not publication material)

| File | Size | What it is |
|------|------|-----------|
| `covid_atlas_run.log` | 11 MB | Full MATLAB stdout, COVID Atlas benchmark run |
| `lupus_run2.log` | 3.2 MB | Full MATLAB stdout, Lupus run with ML (run_ml=true) |
| `smillie2019_run.log` | 2.8 MB | Full MATLAB stdout, Smillie IBD benchmark run |
| `lupus_run.log` | 36 KB | MATLAB stdout, Lupus first run (run_ml=false) |
| `kang2018_run.log` | 28 KB | MATLAB stdout, Kang 2018 benchmark run |
| `nichenet_run.log` | 8 KB | Rscript stdout, NicheNet comparison |
| `liana_comparison.log` | 4 KB | LIANA comparison stdout (early run) |
| `liana_lupus.log` | 4 KB | LIANA Lupus stdout (run 1) |
| `liana_lupus2.log` | 4 KB | LIANA Lupus stdout (run 2) |
| `prep_lupus.log` | 4 KB | prep_lupus.py stdout (first attempt) |
| `prep_lupus_new.log` | 4 KB | prep_lupus.py stdout (final, 40,395 cells) |
| `prep_smillie2019.log` | 4 KB | prep_smillie2019.py stdout |
| `nichenetr_install.log` | 0 B | Empty install log |

**Total untracked: ~102 MB**

> The key information from these logs (final stats, AUROC, FDR-sig counts, runtimes)
> is preserved in human-readable form in **`RUN_LOG.md`** (kept) and in the Results
> section of **`NEUROSCOPE_figures.pdf`** (kept).

---

## Files KEPT in git (publishable results)

**Per-dataset folders** (`kang2018/`, `covid_atlas/`, `lupus/`, `smillie2019/`):
- `channel_rankings.csv` — core NEUROSCOPE output, all L-R pairs ranked by capacity
- `comparison_liana_vs_neuroscope.csv` — LIANA rank_aggregate comparison
- `comparison_cellchat_vs_neuroscope.csv` — CellChat comparison
- `comparison_nichenet_vs_neuroscope.csv` — NicheNet LR-network comparison
- `*_science_report.txt` — narrative discovery report
- `comparison_summary.txt` — text summary (where present)
- `cellchat_results_full.csv`, `cellphonedb_results_full.csv` — full CellChat/CPDB
  permutation-test outputs (Kang, COVID; small, kept for reproducibility)
- `cellchat_R_*.csv/.txt` — native R CellChat outputs (Kang)

**`results/` root:**
- `NEUROSCOPE_figures.pdf` — **the figure + Results deliverable** (1.9 MB)
- `RUN_LOG.md` — master human-readable benchmark log (all 4 datasets)
- `cellchat_comparison_summary.csv`, `nichenet_comparison_summary.csv` — cross-dataset summaries
- `Lupus_*` / `SmillieIBD_*` `specific_interactions` / `target_rankings` CSVs — result tables
- `communication_identity.csv`, `pathway_activity_scores.csv`,
  `pathway_predictor_weights.csv` — analysis outputs

> Per-run generated directories (`paper_figures/`, `ml_models/`, `statistics/`,
> `networks/`, `novel_discoveries/`, `poster_figures/`) were already gitignored before
> this cleanup. The figures they contain are all embedded in `NEUROSCOPE_figures.pdf`.
