# `converted/` Cleanup Manifest

**Date:** 2026-05-31
Files below were **untracked** from git (`git rm --cached`) to keep the public repo
small. **They still exist on your local disk and in git history.** See the top-level
[CLEANUP_MANIFEST.md](../CLEANUP_MANIFEST.md) for restore instructions.

## What `converted/` contains

1. **CellPhoneDB reference tables** (top level) — the ligand-receptor interaction
   database the method needs to run. **KEPT in git** (~580 KB total, required).
2. **Per-dataset expression matrices + metadata** (subfolders) — produced by the
   `prep_*.py` scripts from public single-cell datasets. **UNTRACKED** (large,
   regenerable). Consistent with `.gitignore`, which already excluded
   `converted/kang2018/` and `converted/covid_atlas/`.

---

## Files UNTRACKED from git (kept on disk)

| File | Size | What it is | How to regenerate |
|------|------|-----------|-------------------|
| `nichenet/ligand_target_matrix.rds` | **262 MB** | NicheNet prior ligand→target regulatory-potential matrix (33,354 × 1,226), from Zenodo 7074291. **Exceeds GitHub's 100 MB hard limit.** | Auto-downloaded by `compare_nichenet.R` (cached here). Direct: Zenodo 7074291 `ligand_target_matrix_nsga2r_final.rds`. |
| `nichenet/lr_network.rds` | 20 KB | NicheNet curated human L-R network (4,986 interactions). | Auto-downloaded by `compare_nichenet.R`. Direct: Zenodo 7074291 `lr_network_human_21122021.rds`. |
| `lupus/lupus_expression_sparse.mat` | 87 MB | Sparse expression, 40,395 cells × 4,181 genes (Perez et al. 2022 SLE PBMC). | `python prep_lupus.py` |
| `lupus/lupus_cell_metadata.csv` | 1.7 MB | Cell metadata: CellID, Cluster, donor_id, condition (149 donors). | `python prep_lupus.py` |
| `lupus/lupus_gene_names.csv` | 36 KB | Gene symbols for lupus matrix columns. | `python prep_lupus.py` |
| `smillie2019/smillie_expression_sparse.mat` | 4.8 MB | Sparse expression, 9,385 cells × 3,000 genes (Smillie et al. 2019 colon). | `python prep_smillie2019.py` |
| `smillie2019/smillie_cell_metadata.csv` | 276 KB | Cell metadata for Smillie colon dataset (18 donors). | `python prep_smillie2019.py` |
| `smillie2019/smillie_gene_names.csv` | 24 KB | Gene symbols for smillie matrix columns. | `python prep_smillie2019.py` |

**Total untracked: ~358 MB**

---

## Files KEPT in git (CellPhoneDB reference DB — required to run)

| File | Size | What it is |
|------|------|-----------|
| `gene_table.csv` | 74 KB | Gene → HGNC symbol → protein mappings |
| `protein_table.csv` | 43 KB | Protein multidata entries |
| `interaction_table.csv` | 327 KB | Ligand-receptor interaction definitions |
| `complex_table.csv` | 15 KB | Protein complex definitions |
| `complex_composition_table.csv` | 11 KB | Complex → component protein mappings |
| `multidata_table.csv` | 118 KB | Unified protein + complex multidata index |
| `gene_synonym_to_gene_name.csv` | 24 KB | Gene synonym resolution |
| `receptor_to_transcription_factor.csv` | 3.6 KB | Receptor → downstream TF mappings |

These are the CellPhoneDB v3 reference tables (Efremova et al. 2020). Small, required
by `NEUROSCOPE.m` to identify L-R pairs, kept public.

> `converted/kang2018/` and `converted/covid_atlas/` expression data were already
> gitignored before this cleanup. Regenerate with `python prep_kang2018.py` and
> `python prep_covid_atlas.py`.
