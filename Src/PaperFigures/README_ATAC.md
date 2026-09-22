# Figure 5 — chromatin TF motif accessibility of the TF screen

Publication panels for **Figure 5** (single-cell ATAC / TF motif accessibility story) of the AbbasTFScreen
paper. One notebook per panel; each imports `paperfig_style.py`, loads its data, draws one
panel, writes `panels/<name>.pdf` (vector, editable text), and displays it inline in the
notebook. PDF only — no PNG files are written to disk.

## How to run

All notebooks use the **`scanpy_env`** Jupyter kernel (pandas/numpy/matplotlib/h5py/scanpy).
Open any notebook and run it, or regenerate + execute every panel from the command line:

```bash
# build the .ipynb files only
/home/eraslab1/miniconda3/bin/python _build_notebooks.py
# build AND execute all panels (writes PDFs/PNGs into panels/)
/home/eraslab1/miniconda3/bin/python _build_notebooks.py --run
```

`_build_notebooks.py` is the single source of truth for panel code — edit a panel there and
rerun, or edit the generated notebook directly (both stay in sync only if you regenerate).

## Panels

### Main — Figure 5
| notebook | panel |
|---|---|
| `Fig5A_ATAC_UMAP_cellstate` | ATAC-profile UMAP (IterativeLSI on genome-wide TileMatrix) colored by RNA cell state; ARI(ATAC clusters, RNA states) reports cross-modal concordance |
| `Fig5B_ATAC_UMAP_pseudotime` | Same ATAC-only UMAP colored by RNA pseudotime — tests whether the differentiation trajectory is visible in chromatin |
| `Fig5C_StateTFActivity_heatmap` | TF motif accessibility per cell state (state vs rest), FDR-marked |
| `Fig5D_PseudotimeTF_lines` | TF motif accessibility trajectories across pseudotime — top rising/falling TFs (line view; full heatmap is Fig S5E) |
| `Fig5E_ATACvsRNA_examples` | TF motif accessibility vs RNA expression across pseudotime — ZEB1/GRHL2 (rising), NFIC/ZBTB14 (falling) concordance examples (full set Fig S5H) |

### Supplementary Figure S5

**S5A–C — ATAC QC** (from `FigS5_QC_ATAC.ipynb`; data via `_export_atac_qc.R`):

| panel | file | content |
|---|---|---|
| S5A | `FigS5A_QC_CellCalling` | TSS enrichment vs log10(unique fragments), 2D density over all 8.17M barcodes — cell calling |
| S5B | `FigS5B_QC_Distributions` | per-cell TSS enrichment, log10 fragments, FRIP (medians 11.8 / 6,734 / 0.29) |
| S5C | `FigS5C_QC_TSSprofile` | aggregate TSS enrichment profile — sharp peak + nucleosome shoulder |

*(Excluded on purpose: fragment-size / nucleosome QC — the width-1 insertion-site input makes
`nDiFrags = nMultiFrags = NucleosomeRatio = 0`, so those ArchR plots are degenerate.)*

**S5D–H — analysis** (individual notebooks in `_build_notebooks.py`):

| panel | notebook | content |
|---|---|---|
| S5D | `FigS5D_PseudotimeSegments_UMAP` | 16 pseudotime segments on the ATAC UMAP vs the RNA UMAP |
| S5E | `FigS5E_PseudotimeTFActivity_heatmap` | Full pseudotime TF motif accessibility heatmap (all sig TFs; main Fig 5D is the line view) |
| S5F | `FigS5F_PertTFActivity_heatmap` | TF motif accessibility per perturbation, pooled vs NTC, pseudotime peak set (few hits) |
| S5G | `FigS5G_StateXPert_TFActivity` | TF motif accessibility per state × perturbation (vs same-state NTC), grouped by state |
| S5H | `FigS5H_TFActivity_vs_mRNA` | TF motif accessibility vs RNA per segment — remaining individual TFs (examples in main Fig 5E); mRNA from `_export_mrna_by_segment.py` |

Regenerate QC:
```bash
Rscript _export_atac_qc.R                                           # arrow/project -> CSV_Files/QC/
/home/eraslab1/miniconda3/bin/python _build_qc_notebook.py --run   # CSV -> panels
```

## Data sources (read-only)
- `../../CSV_Files/*.csv` — TF motif accessibility coef/FDR matrices from pipeline stages 05/07/12/13/14.
- `../../CSV_Files/ATAC_UMAP.csv` — ATAC-only UMAP coords + RNA cell state + ATAC cluster,
  produced by `_compute_atac_umap.R` (IterativeLSI on a genome-wide TileMatrix → UMAP on the
  cell-state ArchR project). Regenerate: `Rscript _compute_atac_umap.R` (long ArchR job).
- `../../ComboScreen_processed.h5ad` — `obs/{multimodal_cell_names,dpt_pseudotime}` (RNA
  pseudotime is mapped onto ATAC cells via the shared multiome barcode).

## Shared style (`paperfig_style.py`)
Central palettes (states, diverging activity colormap), publication rcParams (7 pt sans,
vector-editable PDF text), the NE→Differentiated `STATE_ORDER`, data loaders, and
`savepanel()`. Change the look in one place and rerun to restyle every panel consistently.
