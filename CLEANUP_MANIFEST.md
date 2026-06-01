# NEUROSCOPE Repository Cleanup Manifest

**Date:** 2026-05-31
**Purpose:** Prepare the repository for public GitHub publication.

**Method:** Files listed here were **untracked** (`git rm --cached`), NOT deleted.
Every one of them:
- **still exists on your local disk** (your working copy is untouched), and
- **remains in git history** (recoverable from any prior commit).

Nothing was permanently destroyed. All source code was kept. This document plus the
per-folder manifests record exactly what was removed from version control and how to
bring it back.

---

## Why this cleanup was needed

Large data/intermediate/log files had been committed to git because the `.gitignore`
patterns did not match their actual filenames:

1. **`converted/nichenet/ligand_target_matrix.rds` is 262 MB.** GitHub **hard-blocks
   any file larger than 100 MB** — a `git push` would have failed outright.
2. `.gitignore` ignored `converted/kang2018/` and `converted/covid_atlas/` but the
   newer `converted/lupus/`, `converted/smillie2019/`, and `converted/nichenet/`
   directories were never added, so their large matrices got committed.
3. `.gitignore` ignored `results/**/*_log.txt`, but the actual run logs are named
   `*_run.log` / `*.log`, so multi-megabyte logs slipped through.
4. `results/**/liana_all_methods_*.csv` are raw per-cluster-pair dumps (up to 84 MB) —
   the **summarized** `comparison_*_vs_neuroscope.csv` files are the real results and
   were kept.

Total removed from version control: **~460 MB** across 25 files.

---

## What was KEPT (all source code + all summarized results)

- **All source code:** `NEUROSCOPE.m`, everything in `functions/`, all `run_*.m`,
  all `prep_*.py`, all `compare_*.py` / `compare_*.R`, `make_figure_doc.py`,
  `power_analysis_ksg.m`, and the install helpers. (User decision: keep all scripts
  for reproducibility.)
- **The figure deliverable:** `results/NEUROSCOPE_v52_figures.pdf` (1.9 MB) — contains
  every benchmark figure and the full Results section.
- **All summarized results:** `channel_rankings.csv`, `comparison_*_vs_neuroscope.csv`,
  `*_science_report.txt`, `comparison_summary.txt`, `*_specific_interactions.csv`,
  `*_target_rankings.csv`, `cellchat_comparison_summary.csv`,
  `nichenet_comparison_summary.csv`, `RUN_LOG.md`.
- **CellPhoneDB reference tables** in `converted/` top level (small, required to run):
  `gene_table.csv`, `protein_table.csv`, `interaction_table.csv`, `complex_table.csv`,
  `complex_composition_table.csv`, `multidata_table.csv`,
  `gene_synonym_to_gene_name.csv`, `receptor_to_transcription_factor.csv`.

---

## How to restore any untracked file

The file is already on your disk — nothing is needed to keep using it locally. To put
it back under version control:

```bash
git add -f <path>
```

To recover a copy from git history if the local file is ever lost:

```bash
git checkout c5bea60 -- <path>     # c5bea60 = last commit before cleanup
```

To regenerate from scratch, see the per-folder manifests for the exact command.

---

## Detailed per-folder manifests

- [converted/CLEANUP_MANIFEST.md](converted/CLEANUP_MANIFEST.md) — 8 data files (~358 MB)
- [results/CLEANUP_MANIFEST.md](results/CLEANUP_MANIFEST.md) — 17 log/intermediate files (~102 MB)

**Last commit before cleanup:** `c5bea60` ("Add comprehensive Results section to PDF").
Every file is in its fully-tracked state in that commit and earlier.
