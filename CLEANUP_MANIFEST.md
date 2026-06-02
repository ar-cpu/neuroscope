# NEUROSCOPE Repository Cleanup Manifest

**Date:** 2026-05-31  
**Purpose:** Prepare the repository for public GitHub publication.  
**Method:** Files were untracked from git (`git rm --cached`), NOT deleted from disk. Scripts were moved to a `scripts/` subfolder. Per-folder manifests have full details.

---

## Nothing was permanently destroyed

Every removed file:
- **still exists on your local disk** (working copy is untouched)
- **remains in git history** (recoverable from any commit before `3961d61`)

To restore a file to git: `git add -f <path>`  
To recover from history: `git checkout 3961d61^ -- <path>`

---

## What was done

### Commit `3961d61` — Untrack large data/log files (~460 MB total)

| Category | Files | Size |
|----------|-------|------|
| NicheNet prior matrices | `converted/nichenet/*.rds` | **262 MB** (GitHub hard-blocks >100 MB) |
| Lupus expression data | `converted/lupus/*.mat` + metadata | 88 MB |
| Smillie expression data | `converted/smillie2019/*.mat` + metadata | 5 MB |
| LIANA intermediate dumps | `results/**/liana_all_methods_*.csv` | 95 MB |
| Run logs | `results/**/*run*.log`, `*.log` | 20 MB |

See [converted/CLEANUP_MANIFEST.md](converted/CLEANUP_MANIFEST.md) and [results/CLEANUP_MANIFEST.md](results/CLEANUP_MANIFEST.md) for per-file regeneration commands.

### This commit — Move scripts + fix portable paths

Moved data-prep and comparison scripts to `scripts/`. Fixed hardcoded `C:/Users/arjun/...` absolute paths in 5 files (replaced with portable `Sys.getenv("R_LIBS_USER")` / relative paths so any user can run them).

---

## Repo root after cleanup

| Item | Purpose |
|------|---------|
| `NEUROSCOPE.m` | Main pipeline entry point |
| `functions/` | All MATLAB pipeline modules (untouched) |
| `run_*.m` | Per-dataset launcher configs |
| `make_figure_doc.py` | Figure PDF generator |
| `power_analysis_ksg.m` | KSG MI power analysis |
| `converted/` | CellPhoneDB reference tables (small, required to run) |
| `results/` | All benchmark outputs, comparison CSVs, figures PDF |
| `scripts/` | Data prep + comparison scripts (moved from root) |
| `README.md`, `RESULTS.md`, `LICENSE` | Documentation |

---

## `scripts/` folder contents

| Script | Category | What it does |
|--------|----------|-------------|
| `prep_kang2018.py` | PREP | Download + convert Kang 2018 IFN-beta PBMC from GEO |
| `prep_covid_atlas.py` | PREP | Stream COVID Blood Atlas from CellxGene |
| `prep_lupus.py` | PREP | Stream Lupus PBMC from CellxGene (Perez 2022) |
| `prep_smillie2019.py` | PREP | Download Smillie 2019 colon from CellxGene |
| `compare_liana.py` | COMPARE | Run LIANA rank_aggregate vs NEUROSCOPE (4 datasets) |
| `compare_cellchat.py` | COMPARE | Full CellChat permutation test (Kang 2018) |
| `compare_cellchat_all.py` | COMPARE | Extract CellChat scores from LIANA output (all 4) |
| `compare_cellchat_covid.py` | COMPARE | CellChat for COVID Atlas |
| `compare_nichenet.R` | COMPARE | NicheNet LR-network membership (all 4 datasets) |
| `run_cellchat_r.R` | COMPARE | Native R CellChat (Kang 2018 validation) |
| `install_nichenetr.R` | INSTALL | Install nichenetr from GitHub (saeyslab/nichenetr) |
| `install_r_cellchat.R` | INSTALL | Install R CellChat dependencies |

All scripts run from the **repository root**: `python scripts/prep_kang2018.py`
