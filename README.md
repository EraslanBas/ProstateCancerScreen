# Prostate cancer combinatorial Perturb-seq screen

Analysis code for a combinatorial CRISPR knockout screen in neuroendocrine prostate
cancer (NEPC). Thirteen transcription factors plus non-targeting controls (NTC) were
knocked out singly and in pairs and assayed by scRNA-seq at two time points
(`day04`, `day10`). After removing cells carrying three or more guide assignments,
291,301 cells across 92 perturbation groups enter the analysis.

The repository produces the panels of **Figure 4**, **Figure 5F** and
**Supplementary Figures S01 to S17**: 28 vector PDFs written to
`Src/PaperFigures/figures/`.

## External inputs

The pipeline starts from two AnnData objects that are **produced outside this
repository**:

| file | contents |
| --- | --- |
| `Data/NEPC_sec_screen_RNA_scanpy_raw_counts.h5ad` | raw UMI counts, all cells and genes |
| `Data/NEPC_merged_sec_screen_RNA_scanpy_final_singlets.HVG_downstream.h5ad` | the same cells after upstream QC, carrying the 35 `obs` columns everything downstream relies on |

Those 35 columns include the 14 guide one-hots (`ASCL1` ... `ZNF776`, `NTC`),
`N_genes_targeted`, `time_point`, the transcriptional state label `final_label`,
the six state signature scores plus `Doxo-1-differentiated_score`, cell-cycle
`S_score` / `G2M_score` / `phase`, and the QC metrics `UMI_counts` and
`gene_number`.

**The code that produces these two files is not in this repository.** That covers
demultiplexing, cell calling, ambient and doublet filtering, guide assignment,
the clustering that yields `final_label`, cell-cycle scoring and the state
signature scoring. Everything downstream of these two objects is reproducible
here; nothing upstream of them is.

Two further inputs have no code that regenerates them:

- **`Src/ComboScreen/Doxo_1_differentiated.DEGs.csv`** (committed, 214 genes) is a
  scanpy `rank_genes_groups` export that defines the Doxo1 differentiation score.
  That score is the backbone of Figures 4B, 4F, 4J, 4K and 5F. It is versioned here
  as a fixed input; the comparison that produced it is not recorded.
- **`Data/combined_sgRNA_assignment_df.final.csv`** is read by
  `Src/ComboScreen/InvestigateGuideAssignment.ipynb`, which is exploratory and not
  on the figure path.

## Pipeline

Run in this order. Steps 1 to 3 are one-off and expensive; step 4 onwards is what
you repeat when a figure changes.

| # | code | reads | writes |
| --- | --- | --- | --- |
| 1 | `Src/ComboScreen/00_Build_ComboScreen.py` | the two external h5ads | `Data/ComboScreen.h5ad` (449,267 x 23,423) |
| 2 | `Src/ComboScreen/ComputeDEPdex-Copy1.ipynb`, `-Copy2.ipynb` | `ComboScreen.h5ad` | `Day04DEGs.csv`, `Day10DEGs.csv` (per-perturbation DE vs NTC, via `pdex`) |
| 3 | `Src/PaperFigures/_add_fdr8000.py` | those two | `Day{04,10}DEGs_fdr8000.csv` (BH recomputed over an 8,000-gene universe) |
| 4 | `Src/PaperFigures/00_Prepare_data.ipynb` | `ComboScreen.h5ad`, the Doxo1 signature | `Data/ComboScreen_processed.h5ad` |
| 5 | the four figure notebooks | `ComboScreen_processed.h5ad`, the fdr8000 tables | 28 PDFs in `Src/PaperFigures/figures/` |

Step 4 normalises to 20,000 UMIs and log1p-transforms, keeps raw counts in
`.layers['counts']`, scores the Doxo1 signature, computes PCA and UMAP on a
transient scaled copy, and runs diffusion pseudotime rooted at the diffusion-space
medoid of the lowest 1% of Doxo1 score, oriented NEPC-N to NEPC-A.

### Figure notebooks

Each notebook is named after the panels it writes.

| notebook | main panels | supplementary |
|---|---|---|
| `00_Prepare_data.ipynb` | none | none |
| `Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb` | 4B, 4C, 4F, 4G, 4H, 4J | S07, S08, S10, S11, S16, S17 |
| `Fig4DE_S01-S06_S09__Overview_UMAP.ipynb` | 4D, 4E | S01 to S06, S09 |
| `Fig4K_S12-S14__Interaction_Doxo1.ipynb` | 4K | S12, S13, S14 |
| `Fig4L_Fig5F_S15__Perturbation_Gene_Effects.ipynb` | 4L, 5F | S15 |

The Pseudotime notebook must run before the Interaction and Perturbation
notebooks: those two read `figures/Fig2_state_enrichment_{day}.csv` to define
their perturbation set. Nothing enforces this in code.

Panel legends and the statistics as actually implemented are in
`Src/PaperFigures/Figure4_legends.md`, `Figure5_legends.md` and
`Supplementary_figures.md`.

### Panel naming

Notebooks label their own panels with an internal key (`pseudotime_a2_umap_doxo1`).
`PAPER_PANELS` in `Src/PaperFigures/_figutils.py` maps each key to its published
name (`Fig4B_umap_doxo1`), and `savefig` writes the published name. A key mapped to
`None` is deliberately not published.

This indirection means renumbering a figure is a one-line edit in `PAPER_PANELS`.
Run `python Src/PaperFigures/_panel_names.py` afterwards: it re-derives what each
notebook writes and reports any filename that no longer matches, and `--fix`
renames them.

## Environment

The figure notebooks run under the conda `base` environment (scanpy 1.9.1,
anndata 0.9.2, pandas 1.4.3, matplotlib 3.5.2, seaborn 0.11.2, statsmodels 0.13.2).
`Src/ComboScreen/Install_python_req.sh` creates a minimal venv. Step 2 additionally
needs the `pdex` package and runs in its own environment.

Execute a notebook headlessly with:

```
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=-1 <notebook>
```

## Data files

No data is tracked. `.gitignore` is an allowlist: it excludes everything, then
re-includes source and documentation extensions. The only data exceptions are three
small input tables that no code regenerates (the Doxo1 signature, `Programs_K14.csv`
and `Human_TFs.csv`); every derived table, h5ad and PDF stays out.

The h5ad objects are large (`ComboScreen.h5ad` 13 GB, `ComboScreen_processed.h5ad`
19 GB) and the per-perturbation DE tables are 300 to 400 MB each.

## Other directories

Not on the figure path, kept for reference:

- `Src/*.ipynb`, `Src/RunCNMF.py`: earlier pilot analysis (QC, optimal-transport
  distances between knockout and target cells, cNMF programs).
- `Src/PilotData/`: pilot dataset notebooks.
- `Src/ComboScreen/`: the wider ComboScreen analysis, including
  `01_Data_outlook.ipynb` (the original exploratory pass) and several exploratory
  notebooks. Only `00_Build_ComboScreen.py` and the two `ComputeDEPdex` notebooks
  are on the figure path.
- `Src/Src/`: an older duplicate of `Src/ComboScreen/` retained from before the
  reorganisation.
