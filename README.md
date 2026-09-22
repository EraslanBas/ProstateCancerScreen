# Prostate cancer combinatorial Perturb-seq screen

Analysis code for a combinatorial CRISPR knockout screen in neuroendocrine prostate
cancer (NEPC). Thirteen transcription factors plus non-targeting controls (NTC) were
knocked out singly and in pairs and assayed by scRNA-seq at two time points
(`day04`, `day10`). After removing cells carrying three or more guide assignments,
291,301 cells across 92 perturbation groups enter the analysis.

The screen was assayed as a **multiome** experiment, so this repository produces the
panels of **Figure 4**, **Figure 5** and the **supplementary figures** from two
modalities whose figure code shares `Src/PaperFigures/`:

- **RNA** (Perturb-seq) — Figure 4, Figure 5F and Supplementary Figures S01–S17.
- **ATAC** (chromatin accessibility) — Figure 5A–E and its supplements; the ATAC
  pipeline lives under `Src/ComboScreen/ATAC/` (see
  [ATAC (chromatin) analysis](#atac-chromatin-analysis)).

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

One further input has no code that regenerates it:

- **`Src/ComboScreen/Doxo_1_differentiated.DEGs.csv`** (committed, 214 genes) is a
  scanpy `rank_genes_groups` export that defines the Doxo1 differentiation score.
  That score is the backbone of Figures 4B, 4F, 4J, 4K and 5F. It is versioned here
  as a fixed input; the comparison that produced it is not recorded.

## Pipeline

Run in this order. Steps 1 to 3 are one-off and expensive; step 4 onwards is what
you repeat when a figure changes.

| # | code | reads | writes |
| --- | --- | --- | --- |
| 1 | `Src/ComboScreen/00_Build_ComboScreen.py` | the two external h5ads | `Data/ComboScreen.h5ad` (449,267 x 23,423) |
| 2 | `Src/ComboScreen/ComputeDE_day04.ipynb`, `ComputeDE_day10.ipynb` | `ComboScreen.h5ad` | `Day04DEGs.csv`, `Day10DEGs.csv` (per-perturbation DE vs NTC, via `pdex`; the two differ only in the reference group, `NTC_day04` and `NTC_day10`) |
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

The **ATAC** panels (Figure 5A–E and the S5 supplements) are produced by the ATAC
notebooks in the same folder (`Fig5A_*`–`Fig5E_*`, `FigS5*`), generated from
`_build_notebooks.py` / `_build_qc_notebook.py`; see the ATAC section below.

Panel legends and the statistics as actually implemented are in
`Src/PaperFigures/Figure4_legends.md`, `Figure5_legends.md` (all of Figure 5, RNA + ATAC)
and `Supplementary_figures.md` (all supplements, RNA + ATAC).

### Panel naming

Notebooks label their own panels with an internal key (`pseudotime_a2_umap_doxo1`).
`PAPER_PANELS` in `Src/PaperFigures/_figutils.py` maps each key to its published
name (`Fig4B_umap_doxo1`), and `savefig` writes the published name. A key mapped to
`None` is deliberately not published.

This indirection means renumbering a figure is a one-line edit in `PAPER_PANELS`.
Run `python Src/PaperFigures/_panel_names.py` afterwards: it re-derives what each
notebook writes and reports any filename that no longer matches, and `--fix`
renames them. This applies to the RNA notebooks; the ATAC notebooks name their panels
directly through `savepanel` in `paperfig_style.py`.

## Environment

The RNA figure notebooks run under the conda `base` environment (scanpy 1.9.1,
anndata 0.9.2, pandas 1.4.3, matplotlib 3.5.2, seaborn 0.11.2, statsmodels 0.13.2).
`Src/ComboScreen/Install_python_req.sh` creates a minimal venv. Step 2 additionally
needs the `pdex` package and runs in its own environment. The ATAC figure notebooks
run under a separate `scanpy_env` kernel (scanpy + leidenalg); the ATAC pipeline needs
R + ArchR.

Execute a notebook headlessly with:

```
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=-1 <notebook>
```

## Data files

`.gitignore` is an allowlist: it excludes everything, then re-includes source and
documentation extensions. Large data (h5ad objects, ArchR arrows, the big derived
matrices) and all PDFs stay out. Two sets of small result tables are committed as
exceptions, because the figures read them and no lightweight code regenerates them:
`Src/ComboScreen/Doxo_1_differentiated.DEGs.csv` (an RNA input) and the ATAC result
CSVs under `Src/ComboScreen/ATAC/CSV_Files/` (allowlisted so the chromatin panels
reproduce without the ArchR data; the >100 MB matrices there are still excluded).

The h5ad objects are large (`ComboScreen.h5ad` 13 GB, `ComboScreen_processed.h5ad`
19 GB) and the per-perturbation DE tables are 300 to 400 MB each.

## Repository contents

Every tracked file is on the figure path. `Src/ComboScreen/` holds the RNA steps that
build the inputs (`00_Build_ComboScreen.py`, `ComputeDE_day04.ipynb`,
`ComputeDE_day10.ipynb`), their shared imports (`libraries.py`), the Doxo1 signature
and the environment script, plus the **`ATAC/`** subtree (chromatin pipeline + result
CSVs). `Src/PaperFigures/` holds the shared figure code for both modalities: the RNA
preparation and figure notebooks with `_figutils.py` / `_panel_names.py`, the ATAC
notebooks with `paperfig_style.py` / `_build_notebooks.py`, and the panel legends.

Exploratory work and earlier analyses were removed after the first release:
the ComboScreen exploratory notebooks (`01_Data_outlook`, `02_TestDETFs`,
`03_PlotBetaMatrices`, `InteractionModel`, `InvestigateGuideAssignment`,
`TestDoxo1Signature`) and the pilot and single-knockout analyses, which ran on
different datasets (`adataALL.h5ad`, `adata_20K.h5ad`, `AbbasAnndata.h5ad`).
All of it remains in the git history.

## ATAC (chromatin) analysis

`Src/ComboScreen/ATAC/` holds the single-cell **ATAC** arm of the same multiome experiment:
reproducible peak calling, motif-based **TF motif accessibility**, and the chromatin
panels of **Figure 5A–E** plus supplements. Built on ArchR (hg38) + HOCOMOCO v11.

**Pipeline (`Src/ComboScreen/ATAC/`).** ArrowFiles + reproducible peaks (`CreateArrowFiles.r`,
`parallel_arrows/`, `01_*`), motif matrix + TF grouping (`02_*`, `06_*`), then
differential peaks and per-peak-Log2FC → motif regression giving TF motif
accessibility per cell state, pseudotime segment and perturbation (`03_*`–`14_*`).
These steps read the ArchR arrows/projects (large, external, on local storage) and
are **not** reproducible from a clone alone.

**Figures (`Src/PaperFigures/`, shared with the RNA figures).** The ATAC figure code
(`_build_notebooks.py`, `_build_qc_notebook.py`, `paperfig_style.py`, `_export_*`) lives
alongside the RNA notebooks; it generates one notebook per panel (kernel `scanpy_env`):

- **Fig 5A/B** — UMAP from the ATAC profile (IterativeLSI on a genome-wide TileMatrix),
  coloured by the RNA-defined cell state / pseudotime — cross-modal concordance
  (ARI ≈ 0.68 between ATAC clusters and RNA states; |ρ| ≈ 0.71 for pseudotime).
- **Fig 5C** — TF motif accessibility per cell state; **5D** — along pseudotime
  (root NE → differentiated tip); **5E** — TF motif accessibility vs RNA expression
  (examples: GRHL2 concordant; ZEB1/NFIC/ZBTB14 discordant).
- **Supplements** — ATAC QC, pseudotime segments (ATAC vs RNA UMAP), per-perturbation
  and state×perturbation accessibility, and TF-motif-accessibility-vs-mRNA. Legends
  in `Src/PaperFigures/Figure5_legends.md` (5A–E) and `Supplementary_figures.md` (S5A–H).

**Data policy (differs from the RNA side).** The small/medium ATAC result CSVs *are*
committed under `Src/ComboScreen/ATAC/CSV_Files/` (allowlisted in `.gitignore`) so the chromatin
panels reproduce without the ArchR data; the large regenerable matrices (`ATAC_LSI`,
`PertPeaksFC*`, the motif matrices, `PseudotimePeaksFC`) stay out. The UMAP /
pseudotime panels additionally read the external `ComboScreen_processed.h5ad`
(19 GB; set `ATAC_H5AD` to its path).
