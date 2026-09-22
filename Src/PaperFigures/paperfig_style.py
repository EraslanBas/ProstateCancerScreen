"""
paperfig_style.py  —  shared style + data helpers for the AbbasTFScreen Figure 5 panels.

Every panel notebook starts with `from paperfig_style import *`. This module owns:
  * publication matplotlib rcParams (vector-editable PDF text, consistent fonts/sizes)
  * the project palettes (cell states, TF classes, diverging activity colormap)
  * the canonical NE->Differentiated trajectory order of the cell states
  * path constants + small loaders for the h5ad and the CSV_Files result matrices
  * savepanel(): writes each panel to panels/<name>.pdf AND panels/<name>.png

Design decisions (keep consistent across all panels):
  - Cells:        one modality (ATAC-informed) TF-activity story; states colored identically everywhere.
  - Diverging:    TF activity coef uses a symmetric blue-white-red map centered at 0.
  - FDR marker:   '*' = FDR <= FDR_THRESHOLD (0.1), the threshold used across the pipeline.
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

# ----------------------------------------------------------------------------- paths
# This file lives in <repo>/Src/PaperFigures (shared RNA+ATAC figures). The ATAC
# pipeline outputs (result CSVs, motif inputs) live under <repo>/Src/ComboScreen/ATAC.
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "ComboScreen", "ATAC"))
CSV_DIR      = os.path.join(PROJECT_ROOT, "CSV_Files")
DATA_DIR     = os.path.join(PROJECT_ROOT, "Data")
# The RNA/multiome AnnData is a large EXTERNAL input (not in the repo); override with
# ATAC_H5AD if it lives elsewhere. Only the UMAP/pseudotime panels (5B, S5D) read it.
H5AD         = os.environ.get("ATAC_H5AD",
                              "/home/eraslab1/Projects/AbbasTFScreen/ComboScreen_processed.h5ad")
PANEL_DIR    = os.path.join(os.path.dirname(__file__), "panels")
os.makedirs(PANEL_DIR, exist_ok=True)

FDR_THRESHOLD = 0.1          # a TF activity call is "significant" if FDR <= this (pipeline-wide)

# ----------------------------------------------------------------------------- rcParams
def set_style():
    mpl.rcParams.update({
        # vector text stays editable in Illustrator/Inkscape (not outlined)
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Nimbus Sans", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "figure.dpi": 120, "savefig.dpi": 400,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "axes.titlepad": 4,
    })

set_style()

# ----------------------------------------------------------------------------- palettes
# Cell states, ordered NE (root) -> Differentiated (tip). Colors: cool (NE) -> warm (diff).
STATE_ORDER = ["Neuroendocrine", "Intermediate-3", "Intermediate-1",
               "Intermediate-2", "Differentiated-2", "Differentiated_1"]
STATE_COLORS = {
    "Neuroendocrine":  "#3B4CC0",   # deep blue  (NE root)
    "Intermediate-3":  "#6F9BD8",
    "Intermediate-1":  "#9DC3B8",
    "Intermediate-2":  "#E7D08B",
    "Differentiated-2":"#E38A46",
    "Differentiated_1":"#B40426",   # deep red   (differentiated tip)
}
# Data keys above are the h5ad `final_label` values; DISPLAY names use the NEPC nomenclature.
STATE_MAP = {
    "Neuroendocrine":   "NEPC-N",
    "Intermediate-1":   "NEPC-A&N-1",
    "Intermediate-2":   "NEPC-A&N-2",
    "Intermediate-3":   "NEPC-A&N-3",
    "Differentiated_1": "NEPC-A-1",
    "Differentiated-2": "NEPC-A-2",
}
def state_label(s):
    """final_label -> NEPC display name (identity if unmapped)."""
    return STATE_MAP.get(s, s)

# The 13 screened TFs (perturbation targets).
SCREENED_TFS = ["ASCL1", "KLF14", "NEUROD1", "NEUROG1", "NR3C1", "SIM1", "TET2",
                "TWIST1", "VSX1", "ZNF385A", "ZNF547", "ZNF660", "ZNF776"]

# Diverging colormap for TF-activity coefficients (blue = less active, red = more active).
# 3-stop blue-white-red, matching 08_PlotTFActivity_TFScreen.ipynb: colorRampPalette(c("#2166AC","white","#B2182B")).
ACTIVITY_CMAP = LinearSegmentedColormap.from_list(
    "activity", ["#2166AC", "#FFFFFF", "#B2182B"])
PSEUDOTIME_CMAP = "viridis"

import re
# HOCOMOCO motif-model name -> official gene symbol, for consistent labels across panels.
# Irregular cases listed explicitly; regular HOCOMOCO conventions handled by the patterns.
_TF_ALIAS = {
    "GCR":"NR3C1","ANDR":"AR","PRGR":"PGR","NF2L2":"NFE2L2","THA":"THRA","THB":"THRB",
    "NDF1":"NEUROD1","NGN1":"NEUROG1","TWST1":"TWIST1","HTF4":"TCF3","TFE2":"TCF3","ITF2":"TCF4",
    "KAISO":"ZBTB33","COT1":"NR2F1","COT2":"NR2F2","SRBP1":"SREBF1","SRBP2":"SREBF2",
    "HEN1":"NHLH1","BMAL1":"ARNTL","P63":"TP63","P53":"TP53","P73":"TP73","PIT1":"POU1F1",
    "COE1":"EBF1","EVI1":"MECOM","HXB13":"HOXB13","HXA10":"HOXA10","HXA13":"HOXA13","TYY1":"YY1",
}
_TF_PAT = [
    (re.compile(r"^ZN(\d+[A-Z]?)$"),  lambda m: "ZNF" + m.group(1)),   # ZN317  -> ZNF317
    (re.compile(r"^ZBT(\d+[A-Z]?)$"), lambda m: "ZBTB" + m.group(1)),  # ZBT14  -> ZBTB14
    (re.compile(r"^PO(\d)F(\d)$"),    lambda m: f"POU{m.group(1)}F{m.group(2)}"),  # PO3F1 -> POU3F1
    (re.compile(r"^AP2([A-Z])$"),     lambda m: "TFAP2" + m.group(1)),  # AP2A  -> TFAP2A
]
def tf_label(name):
    """HOCOMOCO motif name -> gene symbol for display. TF groups and already-symbols pass through.
    Apply this to EVERY displayed TF name so panels use identical symbols (e.g. ZBT14 -> ZBTB14)."""
    if not isinstance(name, str) or name.startswith("TFGroup"):
        return name
    if name in _TF_ALIAS:
        return _TF_ALIAS[name]
    for pat, fn in _TF_PAT:
        m = pat.match(name)
        if m:
            return fn(m)
    return name

def tf_labels(names):
    return [tf_label(n) for n in names]

# Pseudotime-segment palette — shared with FigS5D so segment colors match across panels.
def segment_cmap(nseg=16):
    return plt.get_cmap("turbo", nseg)

def sym_norm(mat):
    """Symmetric TwoSlopeNorm centered at 0 for a coefficient matrix (robust to all-zero)."""
    v = np.nanmax(np.abs(np.asarray(mat, dtype=float)))
    v = v if (np.isfinite(v) and v > 0) else 1.0
    return TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v)

# ----------------------------------------------------------------------------- loaders
def load_matrix(name):
    """Load a CSV_Files matrix (rows=TF, cols=condition) with the first column as index."""
    return pd.read_csv(os.path.join(CSV_DIR, name), index_col=0)

def load_obs_categorical(f, col):
    """Decode an anndata categorical obs column from an open h5py file -> np.array[str]."""
    g = f["obs"][col]
    cats = [c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]]
    codes = g["codes"][:]
    return np.array([cats[c] if c >= 0 else "NA" for c in codes]), cats

def rna_to_atac_cell(x):
    """RNA/multimodal cell name -> ArchR ATAC cell name (shared 2-level-ligation barcode).
    250804_2lig_10G.ACGG...  ->  TFScreen#250804_2lig_ATAC_10G_ACGG...-1  (matches the R pipeline)."""
    left, seq = x.rsplit(".", 1)
    samp, well = left.rsplit("_", 1)
    return f"TFScreen#{samp}_ATAC_{well}_{seq}-1"

def load_obs_by_atac(col):
    """Map an h5ad obs column onto ATAC cell names via the shared-barcode transform.
    Returns dict {atac_cellName: value} (value col may be plain or anndata-categorical)."""
    import h5py
    with h5py.File(H5AD, "r") as f:
        mm = [x.decode() if isinstance(x, bytes) else x
              for x in f["obs"]["multimodal_cell_names"][:]]
        if isinstance(f["obs"][col], type(f["obs"])) and "categories" in f["obs"][col]:
            val, _ = load_obs_categorical(f, col)
        else:
            val = f["obs"][col][:]
    return {rna_to_atac_cell(m): v for m, v in zip(mm, val)}

def load_umap_obs(cols=("final_label",)):
    """Return (umap[N,2], {col: np.array}) reading only what is needed from the h5ad."""
    import h5py
    out = {}
    with h5py.File(H5AD, "r") as f:
        umap = f["obsm"]["X_umap"][:]
        for c in cols:
            if isinstance(f["obs"][c], type(f["obs"])) and "categories" in f["obs"][c]:
                out[c], _ = load_obs_categorical(f, c)
            else:
                out[c] = f["obs"][c][:]
    return umap, out

# ----------------------------------------------------------------------------- saving
def savepanel(fig, name):
    """Write a panel to panels/<name>.pdf (vector only) and display it inline in the notebook."""
    pdf = os.path.join(PANEL_DIR, name + ".pdf")
    fig.savefig(pdf)
    print(f"saved  {os.path.relpath(pdf, PROJECT_ROOT)}")
    plt.show()          # render inline (embeds a preview in the notebook; no .png file written)
    return pdf

def significant_submatrix(coef, fdr, thr=FDR_THRESHOLD):
    """Coef restricted to rows(TF)/cols(condition) that have >=1 FDR<=thr entry; aligned fdr too."""
    coef, fdr = coef.copy(), fdr.reindex_like(coef)
    sig = fdr <= thr
    r = sig.any(axis=1); c = sig.any(axis=0)
    return coef.loc[r, c], fdr.loc[r, c], sig.loc[r, c]
