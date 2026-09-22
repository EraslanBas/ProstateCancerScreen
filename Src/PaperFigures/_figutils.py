"""
Shared utilities for the manuscript figure notebooks (ComboScreen / NEPC screen).

Every figure notebook does ``import _figutils as fu`` and then calls the helpers
below so that paths, palettes, plot styling and label conventions live in exactly
one place.

The canonical object is ``Data/ComboScreen_processed.h5ad``, (re)built by
``00_Prepare_data.ipynb`` straight from the raw ``ComboScreen.h5ad``. One object is
used for every figure; analyses pull whichever representation they need:
  - ``.layers['counts']``  raw UMI counts            (QC, re-normalisation)
  - ``.X`` and ``.raw``    log1p-normalised expression (gene plots, scoring;
                           ``use_raw=True`` keeps working)
  - ``.obsm['X_pca']``, ``.obsm['X_umap']``, ``.obsm['X_diffmap']``  embeddings
                           (scaling is done only transiently to compute PCA/UMAP)
  - ``.obs['final_label']``       original sub-state label
  - ``.obs['state']``             3-state manuscript label (NEPC-N / NEPC-A&N / NEPC-A)
  - ``.obs['perturbation']``      e.g. 'ASCL1', 'ASCL1+TWIST1', 'NTC'
  - ``.obs['perturbation_clean']`` perturbation with '+NTC' collapsed away
  - ``.obs['perturbation_time']`` perturbation x {day04, day10}
  - ``.obs['dpt_pseudotime']``    diffusion pseudotime (Doxo1-rooted)
  - trajectory signature scores + Doxo1program_score
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns  # noqa: F401  (imported for notebook convenience)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR = Path("/home/eraslab1/Projects/AbbasScreen")
DATA_DIR = PROJECT_DIR / "Data"

RAW_H5AD = DATA_DIR / "ComboScreen.h5ad"                  # 449k cells, raw counts + metadata
PROCESSED_H5AD = DATA_DIR / "ComboScreen_processed.h5ad"  # canonical object (built by 00_Prepare)

FIG_DIR = PROJECT_DIR / "Src" / "PaperFigures" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Doxo-1 differentiated signature (DE gene list); 'names' column holds the genes.
DOXO1_SIGNATURE_CSV = PROJECT_DIR / "Src" / "ComboScreen" / "Doxo_1_differentiated.DEGs.csv"

# ---------------------------------------------------------------------------
# Biological constants
# ---------------------------------------------------------------------------
# One-hot guide columns present in adata.obs (NTC is the control).
GUIDE_COLS = ["ASCL1", "KLF14", "NEUROD1", "NEUROG1", "NR3C1", "NTC", "SIM1",
              "TET2", "TWIST1", "VSX1", "ZNF385A", "ZNF547", "ZNF660", "ZNF776"]
TARGET_GENES = [g for g in GUIDE_COLS if g != "NTC"]

# Original sub-state labels as they appear in obs['final_label'] (trajectory order).
# NB: spelling is inconsistent in the source data ('Differentiated_1' uses an
# underscore, 'Differentiated-2' a hyphen) -- kept verbatim so it matches obs.
RAW_STATE_ORDER = ["Neuroendocrine", "Intermediate-1", "Intermediate-2",
                   "Intermediate-3", "Differentiated_1", "Differentiated-2"]

# Manuscript naming: collapse the sub-states into three biological states along
# the neuroendocrine (NEPC-N) -> mixed (NEPC-A&N) -> adeno-like (NEPC-A) axis.
STATE_MAP = {
    "Neuroendocrine":   "NEPC-N",
    "Intermediate-1":   "NEPC-A&N-1",
    "Intermediate-2":   "NEPC-A&N-2",
    "Intermediate-3":   "NEPC-A&N-3",
    "Differentiated_1": "NEPC-A-1",
    "Differentiated-2": "NEPC-A-2",
}
# Trajectory order: NE (root) -> differentiated.
STATE_ORDER = ["NEPC-N", "NEPC-A&N-1", "NEPC-A&N-2", "NEPC-A&N-3", "NEPC-A-1", "NEPC-A-2"]
TRAJECTORY_ORDER = STATE_ORDER          # backwards-compatible alias

# Convenience groupings (the two collapsed families).
INTERMEDIATE_STATES = ["NEPC-A&N-1", "NEPC-A&N-2", "NEPC-A&N-3"]
DIFFERENTIATED_STATES = ["NEPC-A-1", "NEPC-A-2"]
TIMEPOINT_ORDER = ["day04", "day10"]

CELL_STATE_COL = "state"                # canonical manuscript label (added in 00_Prepare)
RAW_STATE_COL = "final_label"           # original sub-state label
NTC_LABEL = "NTC"

# Trajectory signature scores already present in obs (upstream pipeline).
STATE_SCORE_COLS = ["Neuroendocrine_score", "Intermediate-1_score",
                    "Intermediate-2_score", "Intermediate-3_score",
                    "Differentiated-1_score", "Differentiated-2_score",
                    "Doxo-1-differentiated_score"]

# ---------------------------------------------------------------------------
# Palettes  (colour-blind friendly, ordered along the trajectory)
# ---------------------------------------------------------------------------
STATE_PALETTE = {
    "NEPC-N":     "#f781bf",   # pink  (neuroendocrine)
    # intermediate (brown), shades light -> dark in PSEUDOTIME order: A&N-3 < A&N-1 < A&N-2
    "NEPC-A&N-3": "#d9a26c",   # light brown  (earliest pseudotime)
    "NEPC-A&N-1": "#a65628",   # brown
    "NEPC-A&N-2": "#6b3410",   # dark brown   (latest)
    # differentiated (purple), light -> dark
    "NEPC-A-1":   "#bf8fd0",   # light purple
    "NEPC-A-2":   "#984ea3",   # purple
}
TIME_PALETTE = {"day04": "#4C72B0", "day10": "#DD8452"}
PHASE_PALETTE = {"G1": "#BBBBBB", "S": "#66C2A5", "G2M": "#FC8D62"}

# ---------------------------------------------------------------------------
# Plot theme
# ---------------------------------------------------------------------------
def set_theme():
    """Publication-style defaults: editable text in PDF/SVG, Arial, no frames."""
    sc.set_figure_params(dpi=120, dpi_save=300, fontsize=11,
                         frameon=False, vector_friendly=True)
    mpl.rcParams.update({
        "pdf.fonttype": 42,      # TrueType -> editable in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load(which="processed", backed=None):
    """Read a project AnnData object. which in {processed, raw}."""
    paths = {"processed": PROCESSED_H5AD, "raw": RAW_H5AD}
    if which not in paths:
        raise ValueError(f"which must be one of {list(paths)}")
    return sc.read_h5ad(paths[which], backed=backed)


# ---------------------------------------------------------------------------
# Perturbation labelling
# ---------------------------------------------------------------------------
def combine_onehot(adata, cols=GUIDE_COLS, sep="+"):
    """Collapse the one-hot guide columns into a single 'A+B' perturbation label."""
    def _row(row):
        active = [c for c in cols if row[c] == 1]
        return sep.join(active) if active else "None"
    return adata.obs[cols].apply(_row, axis=1)


def clean_perturbation_labels(pert):
    """Drop the NTC token from combinatorial labels: 'ASCL1+NTC' -> 'ASCL1'.

    A pure-control cell ('NTC') stays 'NTC'.
    """
    def _rm(label):
        parts = [p for p in str(label).split("+") if p != "NTC"]
        return "+".join(parts) if parts else "NTC"
    return pd.Series(np.asarray(pert), dtype=object).map(_rm).values


def apply_state_naming(adata, raw_col=RAW_STATE_COL, new_col=CELL_STATE_COL):
    """Add the 3-state manuscript label (`state`) by collapsing `final_label`.

    Neuroendocrine -> NEPC-N, Intermediate-* -> NEPC-A&N,
    Differentiated_* -> NEPC-A. Returns the present states in trajectory order.
    """
    mapped = adata.obs[raw_col].astype(str).map(STATE_MAP)
    if mapped.isna().any():
        missing = sorted(adata.obs[raw_col].astype(str)[mapped.isna()].unique())
        raise ValueError(f"Unmapped {raw_col} values (add to STATE_MAP): {missing}")
    present = [s for s in STATE_ORDER if s in set(mapped)]
    adata.obs[new_col] = pd.Categorical(mapped, categories=present, ordered=True)
    return present


def ordered_states(adata, col=CELL_STATE_COL):
    """Manuscript states that are actually present, in trajectory order."""
    present = set(adata.obs[col].astype(str).unique())
    return [s for s in STATE_ORDER if s in present]


def set_state_categories(adata, col=CELL_STATE_COL):
    """(Re)derive the manuscript `state` column from `final_label` and order it.

    Always recomputes from `final_label` so the current STATE_MAP/naming takes
    effect on load, without rebuilding the processed object.
    """
    apply_state_naming(adata, new_col=col)
    return ordered_states(adata, col)


# ---------------------------------------------------------------------------
# Pseudotime state regions
# ---------------------------------------------------------------------------
def state_pseudotime_bounds(adata, col=CELL_STATE_COL, pt_key="dpt_pseudotime"):
    """Pseudotime boundaries between consecutive states (midpoints of medians).

    Returns (order, edges, medians): `order` is the present states in trajectory
    order, `edges` has len(order)+1 values [min, b01, b12, ..., max] defining the
    pseudotime span assigned to each state, `medians` the per-state median pt.
    """
    order = ordered_states(adata, col)
    pt = adata.obs[pt_key].values
    med = adata.obs.groupby(col, observed=True)[pt_key].median().reindex(order)
    edges = [float(np.min(pt))]
    for i in range(len(order) - 1):
        edges.append(float((med.iloc[i] + med.iloc[i + 1]) / 2))
    edges.append(float(np.max(pt)))
    return order, edges, med


def shade_state_regions(ax, adata, alpha=0.13, label=False, legend=True, label_y=0.98,
                        fontsize=7, col=CELL_STATE_COL, pt_key="dpt_pseudotime",
                        legend_loc="upper left", legend_bbox=(1.02, 1.0),
                        legend_title="Cell state"):
    """Shade the pseudotime axis of `ax` into coloured cell-state regions.

    With ``legend=True`` (default) the state names are shown as a colour-patch
    legend on the right of the axis (added as a separate artist so it coexists
    with the plot's own data legend). Set ``label=True`` for inline band labels.
    """
    import matplotlib.patches as mpatches
    order, edges, _ = state_pseudotime_bounds(adata, col=col, pt_key=pt_key)
    for i, s in enumerate(order):
        ax.axvspan(edges[i], edges[i + 1], color=STATE_PALETTE[s],
                   alpha=alpha, lw=0, zorder=0)
        if label:
            ax.text((edges[i] + edges[i + 1]) / 2, label_y, s,
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=fontsize, color=STATE_PALETTE[s])
    if legend:
        handles = [mpatches.Patch(facecolor=STATE_PALETTE[s], label=s) for s in order]
        leg = ax.legend(handles=handles, title=legend_title, loc=legend_loc,
                        bbox_to_anchor=legend_bbox, frameon=False, fontsize=fontsize,
                        title_fontsize=fontsize + 1)   # title scales with fontsize
        ax.add_artist(leg)          # keep it when the plot calls ax.legend() for data
    return order, edges


def rolling_trend(x, y, frac=0.10, n_grid=200, min_window=50, min_periods_frac=0.2):
    """Rolling-window local mean ± SD of y along sorted x, sampled on a grid.

    A LOWESS-style smoother: sort by x, take a centred window of
    ``max(min_window, frac*n)`` consecutive cells, and report the running mean
    and SD, interpolated onto ``n_grid`` evenly spaced points spanning the full
    data range. Returns (grid, mean, std).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    order = np.argsort(x, kind="mergesort")
    xs, ys = x[order], y[order]
    w = int(min(max(int(frac * n), min_window), n))
    mp = max(10, int(min_window * min_periods_frac))
    s = pd.Series(ys)
    rmean = s.rolling(w, center=True, min_periods=mp).mean().to_numpy()
    rstd = s.rolling(w, center=True, min_periods=mp).std().to_numpy()
    # Make xs strictly increasing for np.interp (break exact pseudotime ties).
    xs = xs + np.arange(n) * 1e-12
    grid = np.linspace(xs.min(), xs.max(), n_grid)
    good = ~np.isnan(rmean)
    gm = np.interp(grid, xs[good], rmean[good])
    gs = np.interp(grid, xs[good], np.nan_to_num(rstd[good]))
    return grid, gm, gs


def pseudotime_trend(x, y, frac=0.10, n_grid=200, min_window=50, min_cells=5, z=1.96):
    """Local mean + CI band of y along x using a **fixed cell-count** sliding window.

    Same smoother as :func:`rolling_trend` (sort by x; centred window of
    ``max(min_window, frac*n)`` consecutive cells), but it also returns the 95% CI
    of the running mean: ``mean ± z * SE`` with ``SE = sd / sqrt(count)``. The
    window must hold at least ``min_cells`` cells (min_periods); estimates from
    fewer are dropped (only the very edges). Values are interpolated onto
    ``n_grid`` points so the line stays continuous. Returns ``(grid, mean, lo, hi, n)``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nan = np.full(n_grid, np.nan)
    n = x.size
    if n == 0:
        return np.linspace(0, 1, n_grid), nan.copy(), nan.copy(), nan.copy(), np.zeros(n_grid)
    order = np.argsort(x, kind="mergesort")
    xs, ys = x[order], y[order]
    w = int(min(max(int(frac * n), min_window), n))
    s = pd.Series(ys)
    rmean = s.rolling(w, center=True, min_periods=min_cells).mean().to_numpy()
    rstd = s.rolling(w, center=True, min_periods=min_cells).std().to_numpy()
    rcnt = s.rolling(w, center=True, min_periods=min_cells).count().to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        rse = rstd / np.sqrt(rcnt)
    xs = xs + np.arange(n) * 1e-12                      # strictly increasing for interp
    grid = np.linspace(xs.min(), xs.max(), n_grid)
    good = ~np.isnan(rmean) & (rcnt >= min_cells)
    if not good.any():
        return grid, nan.copy(), nan.copy(), nan.copy(), np.zeros(n_grid)
    gm = np.interp(grid, xs[good], rmean[good])
    gse = np.interp(grid, xs[good], np.nan_to_num(rse[good]))
    gcnt = np.interp(grid, xs[good], rcnt[good])
    return grid, gm, gm - z * gse, gm + z * gse, gcnt


def pseudotime_trend_fixedwindow(x, y, n_grid=200, bandwidth=None, frac=0.05,
                                 min_cells=3, xlim=None):
    """Local mean of y using a **fixed-width window along x** (not a cell count).

    At each of ``n_grid`` points spanning ``xlim`` (default the data range), the
    mean is taken over the cells whose x is within ``bandwidth`` (default
    ``frac * range``) of it. The line **disconnects** (NaN) over any interval that
    holds fewer than ``min_cells`` cells -- so a group that has no cells in a
    pseudotime stretch simply has a gap there. Returns ``(grid, mean, n)``.

    Pass ``xlim=(pt.min(), pt.max())`` so every group is drawn on the same axis
    and the gaps line up.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    nan = np.full(n_grid, np.nan)
    lo, hi = xlim if xlim is not None else (
        (float(np.nanmin(x)), float(np.nanmax(x))) if x.size else (0.0, 1.0))
    grid = np.linspace(lo, hi, n_grid)
    if x.size == 0:
        return grid, nan.copy(), np.zeros(n_grid)
    h = bandwidth if bandwidth is not None else frac * (hi - lo)
    order = np.argsort(x, kind="mergesort")
    xs, ys = x[order], y[order]
    csum = np.concatenate([[0.0], np.cumsum(ys)])
    left = np.searchsorted(xs, grid - h, side="left")
    right = np.searchsorted(xs, grid + h, side="right")
    cnt = (right - left).astype(float)
    s = csum[right] - csum[left]
    mean = nan.copy()
    ok = cnt >= min_cells
    with np.errstate(invalid="ignore", divide="ignore"):
        mean[ok] = (s / cnt)[ok]
    return grid, mean, cnt


def plot_pseudotime_trend(ax, x, y, color, label=None, lw=1.8, band=True,
                          band_alpha=0.18, zorder=3, **kw):
    """Draw a density-aware pseudotime trend (mean line + CI band) on ``ax``.

    Thin wrapper over :func:`pseudotime_trend`; gaps in unsupported regions are
    left blank (NaN), so the line breaks where the data can't support it.
    """
    g, m, lo, hi, _ = pseudotime_trend(x, y, **kw)
    line, = ax.plot(g, m, color=color, lw=lw, label=label, zorder=zorder)
    if band:
        ax.fill_between(g, lo, hi, color=color, alpha=band_alpha, lw=0, zorder=zorder - 1)
    return g, m, lo, hi


# ---------------------------------------------------------------------------
# Pseudotime backbone (the main path the trajectory traces through the embedding)
# ---------------------------------------------------------------------------
def pseudotime_backbone(adata, basis="X_umap", pt_key="dpt_pseudotime",
                        n_bins=30, smooth=1, robust=True):
    """The main trajectory path through an embedding, ordered by pseudotime.

    Slice the cells into ``n_bins`` equal-occupancy pseudotime windows
    (quantile edges, so every slice holds ~the same number of cells), take a
    robust centre of each slice in the embedding (median by default), and join
    those centres into one curve running root -> end. This makes the pseudotime
    axis used by the other panels *visible* as a path across the UMAP.

    ``smooth`` applies a centred rolling mean (in bins) to de-jitter the curve.
    Returns (path, pt) -- ``path`` is (n_kept, 2) ordered low->high pseudotime,
    ``pt`` the mean pseudotime of each point (for colouring the line).
    """
    E = np.asarray(adata.obsm[basis][:, :2], dtype=float)
    t = np.asarray(adata.obs[pt_key].values, dtype=float)
    edges = np.quantile(t, np.linspace(0, 1, n_bins + 1))
    edges[-1] = np.nextafter(edges[-1], np.inf)     # include the max in the last bin
    binid = np.clip(np.digitize(t, edges) - 1, 0, n_bins - 1)

    pts, ptc = [], []
    for b in range(n_bins):
        m = binid == b
        if not m.any():
            continue
        c = np.median(E[m], axis=0) if robust else E[m].mean(axis=0)
        pts.append(c); ptc.append(float(t[m].mean()))
    pts = np.asarray(pts); ptc = np.asarray(ptc)

    if smooth and smooth > 0 and len(pts) > 2 * smooth + 1:
        w = 2 * int(smooth) + 1
        sx = pd.Series(pts[:, 0]).rolling(w, center=True, min_periods=1).mean().to_numpy()
        sy = pd.Series(pts[:, 1]).rolling(w, center=True, min_periods=1).mean().to_numpy()
        pts = np.column_stack([sx, sy])
    return pts, ptc


def plot_pseudotime_backbone(adata, ax=None, basis="X_umap", pt_key="dpt_pseudotime",
                             color_by="pseudotime", cmap="viridis", point_size=3,
                             scatter_alpha=0.6, lw=5, line_color="crimson",
                             arrow_color=None, n_arrows=6, n_bins=30, smooth=1,
                             robust=True, mark_ends=True, state_col=None):
    """Embedding scatter with the pseudotime *backbone* path drawn on top.

    ``color_by='pseudotime'`` paints the cloud viridis; ``'state'`` paints it by
    manuscript cell state. The backbone always carries arrowheads pointing root
    -> end. ``line_color`` controls the path: a colour name draws a solid line
    (use this over a pseudotime cloud, where a pseudotime-coloured line would be
    invisible); ``None`` colours the line itself by pseudotime (nice over a state
    cloud). Returns ``(ax, info)`` with ``info['cbar']`` = the mappable to hand to
    ``fig.colorbar`` (the pseudotime cloud, or the line if it carries the scale).
    """
    from matplotlib.collections import LineCollection
    state_col = state_col or CELL_STATE_COL
    arrow_color = arrow_color or (line_color if line_color else "black")
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    E = np.asarray(adata.obsm[basis][:, :2], dtype=float)
    t = np.asarray(adata.obs[pt_key].values, dtype=float)

    # --- background cloud
    cloud = None
    if color_by == "state":
        for s in ordered_states(adata, state_col):
            m = adata.obs[state_col].astype(str).values == s
            ax.scatter(E[m, 0], E[m, 1], s=point_size, color=STATE_PALETTE[s],
                       alpha=scatter_alpha, lw=0, rasterized=True, label=s)
    else:
        cloud = ax.scatter(E[:, 0], E[:, 1], c=t, cmap=cmap, s=point_size, lw=0,
                           alpha=scatter_alpha, rasterized=True)

    # --- backbone path, with a white halo underneath for contrast
    path, ptc = pseudotime_backbone(adata, basis=basis, pt_key=pt_key,
                                    n_bins=n_bins, smooth=smooth, robust=robust)
    ax.plot(path[:, 0], path[:, 1], color="white", lw=lw + 3, zorder=4,
            solid_capstyle="round", solid_joinstyle="round")
    line = None
    if line_color is None:                      # path coloured by pseudotime
        seg = np.stack([path[:-1], path[1:]], axis=1)
        line = LineCollection(seg, cmap=cmap, zorder=5, capstyle="round")
        line.set_array(0.5 * (ptc[:-1] + ptc[1:])); line.set_linewidth(lw)
        ax.add_collection(line)
    else:                                       # solid contrasting path
        ax.plot(path[:, 0], path[:, 1], color=line_color, lw=lw, zorder=5,
                solid_capstyle="round", solid_joinstyle="round")

    # --- arrowheads along the path (direction = increasing pseudotime)
    if n_arrows:
        idx = np.linspace(1, len(path) - 1, n_arrows + 2)[1:-1].astype(int)
        for k in idx:
            ax.annotate("", xy=path[k], xytext=path[k - 1],
                        arrowprops=dict(arrowstyle="-|>", color=arrow_color, lw=0,
                                        mutation_scale=18), zorder=6)
    if mark_ends:
        ax.scatter(*path[0], s=130, marker="o", color="white", edgecolor="black",
                   lw=1.5, zorder=7)
        ax.scatter(*path[-1], s=200, marker="*", color="white", edgecolor="black",
                   lw=1.5, zorder=7)

    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"{basis.replace('X_', '').upper()}1")
    ax.set_ylabel(f"{basis.replace('X_', '').upper()}2")
    info = {"path": path, "pt": ptc, "line": line, "cloud": cloud,
            "cbar": cloud if cloud is not None else line}
    return ax, info


# ---------------------------------------------------------------------------
# Doxo1 double-KO interaction model (shared by the Fig2 score plot and Fig3)
# ---------------------------------------------------------------------------
def interaction_doxo1_pairs(adata, score="Doxo1program_score", genes=None,
                            select_substates=("Intermediate-1", "Intermediate-2",
                                              "Intermediate-3", "Differentiated_1",
                                              "Differentiated-2"),
                            fit_states=("NEPC-A&N-1", "NEPC-A&N-2", "NEPC-A&N-3",
                                        "NEPC-A-1", "NEPC-A-2"), step1_fdr=0.3,
                            one_sided=True):
    """Two-step Doxo1 double-KO selection (per timepoint).

    Step 1 (SELECT) -- for each perturbation, a one-sided Mann-Whitney U test
    (perturbation > NTC) on the Doxo1 score WITHIN each non-NE sub-state in
    ``select_substates``; a perturbation is kept if significant in >=1 sub-state
    (BH FDR < ``step1_fdr``). This is the conservative per-perturbation test
    (NOT the over-powered pooled model).

    Step 2 (DECOMPOSE) -- fit ``score ~ sum(gene_i) + sum(gene_i:gene_j)`` within
    ``fit_states`` and report, for every pair, each effect's beta, p-value and
    BH-FDR: ``beta_g1/pval_g1/fdr_g1``, ``beta_g2/pval_g2/fdr_g2``,
    ``beta_interaction/pval_interaction/fdr_interaction`` and
    ``total_effect/pval_total/fdr_total``. Base-effect and interaction p-values are
    two-sided Wald; the total-effect p is one-sided (increase) when ``one_sided``.
    FDR is BH-corrected across pairs separately for each effect.

    Returns (pairs_by_day, tables_by_day): pairs = step-1-significant double KOs
    (ordered by total effect); tables = the step-2 decomposition for all pairs.
    """
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests
    from scipy.stats import mannwhitneyu
    from itertools import combinations

    genes = list(genes or TARGET_GENES)
    obs = adata.obs
    tp = obs["time_point"].astype(str).values
    sub = obs[RAW_STATE_COL].astype(str).values
    pert = obs["perturbation_clean"].astype(str).values
    Xall = obs[genes].astype(float)
    yall = obs[score].astype(float).values
    fit_subs = [r for r in RAW_STATE_ORDER if STATE_MAP.get(r) in fit_states]

    pairs_by_day, tables_by_day = {}, {}
    for day in TIMEPOINT_ORDER:
        dmask = tp == day

        # ---- Step 1: per-perturbation one-sided MWU within each non-NE sub-state
        ntc = dmask & (pert == NTC_LABEL)
        rows = []
        for p in sorted(set(pert[dmask]) - {NTC_LABEL}):
            pm = dmask & (pert == p)
            for s in select_substates:
                a = yall[pm & (sub == s)]
                b = yall[ntc & (sub == s)]
                if len(a) < 3 or len(b) < 3:
                    continue
                _, pv = mannwhitneyu(a, b, alternative="greater")
                rows.append({"perturbation": p, "sub_state": s, "pval": pv})
        s1 = pd.DataFrame(rows)
        if len(s1):
            s1["fdr"] = multipletests(s1["pval"], method="fdr_bh")[1]
            sig_perts = s1[s1.fdr < step1_fdr].groupby("perturbation")["fdr"].min().index.tolist()
        else:
            sig_perts = []
        sig_doubles = [p for p in sig_perts if "+" in p]

        # ---- Step 2: interaction-model decomposition (fit on fit_states)
        fmask = dmask & np.isin(sub, fit_subs)
        X = Xall[fmask].reset_index(drop=True)
        y = yall[fmask]
        inter = []
        for g1, g2 in combinations(genes, 2):
            v = X[g1] * X[g2]
            if v.sum() > 0:
                X[f"{g1}:{g2}"] = v
                inter.append((g1, g2))
        model = sm.OLS(y, sm.add_constant(X)).fit()
        names = list(model.params.index)
        trows = []
        for g1, g2 in inter:
            col = f"{g1}:{g2}"
            cz = np.zeros(len(names))
            for t in (g1, g2, col):
                cz[names.index(t)] = 1.0
            tt = model.t_test(cz)
            eff = float(tt.effect); p2 = float(tt.pvalue)
            pv = (p2 / 2 if eff > 0 else 1 - p2 / 2) if one_sided else p2
            # base-effect / interaction p-values are the model's two-sided Wald p;
            # the total-effect p (pval_total) is one-sided (increase) when one_sided.
            trows.append({"pair": f"{g1}+{g2}", "gene_1": g1, "gene_2": g2,
                          "beta_g1": float(model.params[g1]), "pval_g1": float(model.pvalues[g1]),
                          "beta_g2": float(model.params[g2]), "pval_g2": float(model.pvalues[g2]),
                          "beta_interaction": float(model.params[col]),
                          "pval_interaction": float(model.pvalues[col]),
                          "total_effect": eff, "pval_total": pv})
        tab = pd.DataFrame(trows)
        for _c in ("g1", "g2", "interaction", "total"):       # BH-FDR per effect across pairs
            tab[f"fdr_{_c}"] = (multipletests(tab[f"pval_{_c}"], method="fdr_bh")[1]
                                if len(tab) else [])
        tab = tab.sort_values("total_effect", ascending=False).reset_index(drop=True)

        rank = tab.set_index("pair")["total_effect"] if len(tab) else pd.Series(dtype=float)
        sig_doubles = [p for p in sig_doubles if p in set(tab.get("pair", []))]
        sig_doubles = sorted(sig_doubles, key=lambda x: rank.get(x, -np.inf), reverse=True)
        pairs_by_day[day] = sig_doubles
        tables_by_day[day] = tab
    return pairs_by_day, tables_by_day


def spline_pseudotime_de(adata, score="Doxo1program_score", pert_col="perturbation_clean",
                         ref=NTC_LABEL, df_spline=5, min_cells=30, fdr=0.05,
                         require_increase=True):
    """Per-perturbation spline-interaction test vs NTC along pseudotime, per day.

    For each perturbation, fit (on NTC + that perturbation's cells, full pseudotime
    range)::

        full:    score ~ bs(pseudotime, df) * C(cond)
        reduced: score ~ bs(pseudotime, df)

    and compare with a nested F-test (statsmodels ``anova_lm``). A significant
    test means the score-vs-pseudotime *curve* differs from NTC (the tradeSeq
    ``conditionTest`` analogue for a single score). BH-correct across perturbations
    within each timepoint.

    Returns (sig_by_day, tables_by_day): sig = perturbations with FDR<fdr (and,
    if ``require_increase``, mean score > NTC); tables = full per-perturbation stats.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    from statsmodels.stats.multitest import multipletests

    obs = adata.obs
    tp = obs["time_point"].astype(str).values
    sig_by_day, tables_by_day = {}, {}
    for day in TIMEPOINT_ORDER:
        sub = obs.loc[tp == day, [pert_col, score, "dpt_pseudotime"]].copy()
        sub.columns = ["pert", "y", "pt"]
        sub["pert"] = sub["pert"].astype(str)
        ntc = sub[sub.pert == ref]
        rows = []
        for p in sorted(set(sub.pert) - {ref}):
            pp = sub[sub.pert == p]
            if len(pp) < min_cells:
                continue
            dd = pd.concat([ntc.assign(cond=0), pp.assign(cond=1)], ignore_index=True)
            try:
                full = smf.ols(f"y ~ bs(pt, df={df_spline}) * C(cond)", data=dd).fit()
                red = smf.ols(f"y ~ bs(pt, df={df_spline})", data=dd).fit()
                pval = float(anova_lm(red, full)["Pr(>F)"].iloc[1])
            except Exception:
                pval = np.nan
            rows.append({"perturbation": p, "n_pert": int(len(pp)),
                         "delta_mean": float(pp.y.mean() - ntc.y.mean()), "pval": pval})
        r = pd.DataFrame(rows).dropna(subset=["pval"])
        r["fdr"] = multipletests(r["pval"], method="fdr_bh")[1] if len(r) else []
        r = r.sort_values("fdr").reset_index(drop=True)
        keep = r[r.fdr < fdr]
        if require_increase:
            keep = keep[keep.delta_mean > 0]
        sig_by_day[day] = keep.sort_values("fdr")["perturbation"].tolist()
        tables_by_day[day] = r
    return sig_by_day, tables_by_day


def wilcoxon_doxo1_state_time(adata, score="Doxo1program_score", pert_col="perturbation_clean",
                              state_col=None, ref=NTC_LABEL, fdr=0.05, min_cells=10,
                              require_increase=True):
    """Single-cell one-sided Wilcoxon (perturbation > NTC) on the Doxo1 score,
    run SEPARATELY within each (timepoint x cell state); BH-corrected across
    perturbations within each (timepoint, state) group.

    Returns (sig_union, table): sig_union = perturbations significant (FDR<fdr,
    and increasing if require_increase) in >=1 (timepoint, state) test; table =
    the full per-(timepoint, state, perturbation) stats.
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    state_col = state_col or CELL_STATE_COL
    obs = adata.obs
    tp = obs["time_point"].astype(str).values
    sv = obs[state_col].astype(str).values
    pv = obs[pert_col].astype(str).values
    y = obs[score].astype(float).values
    states = [s for s in STATE_ORDER if s in set(sv)]
    out = []
    for day in TIMEPOINT_ORDER:
        for s in states:
            base = (tp == day) & (sv == s)
            ntc = y[base & (pv == ref)]
            if len(ntc) < min_cells:
                continue
            rows = []
            for p in sorted(set(pv[base]) - {ref}):
                x = y[base & (pv == p)]
                if len(x) < min_cells:
                    continue
                _, pval = mannwhitneyu(x, ntc, alternative="greater")
                rows.append({"timepoint": day, "state": s, "perturbation": p,
                             "n_pert": int(len(x)),
                             "delta_mean": float(x.mean() - ntc.mean()), "pval": float(pval)})
            d = pd.DataFrame(rows)
            if len(d):
                d["fdr"] = multipletests(d["pval"], method="fdr_bh")[1]  # BH within (day, state)
                out.append(d)
    cols = ["timepoint", "state", "perturbation", "n_pert", "delta_mean", "pval", "fdr"]
    res = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=cols)
    m = res["fdr"] < fdr
    if require_increase:
        m &= res["delta_mean"] > 0
    sig_union = sorted(res.loc[m, "perturbation"].unique())
    return sig_union, res


# ---------------------------------------------------------------------------
# Hierarchical (mixed) model: perturbation effect on Doxo1, state = random effect
# ---------------------------------------------------------------------------
def mixedlm_doxo1_state(adata, score="Doxo1program_score", pert_col="perturbation_clean",
                        state_col=None, ref=NTC_LABEL, fdr=0.05, min_cells=30,
                        require_increase=True, reml=True):
    """Per-timepoint linear mixed model ``score ~ C(perturbation) + (1 | state)``.

    Tests which perturbations raise the Doxo1 score vs NTC while **correcting for
    cell state as a random intercept**: every state gets its own baseline and each
    perturbation's fixed-effect coefficient is the within-state shift vs NTC,
    pooled across states (so a perturbation that scores high only by relocating
    cells to later states no longer counts -- the state random effect absorbs that
    baseline). One-sided Wald test (coef > 0) on each perturbation coefficient,
    BH-corrected across perturbations within each timepoint.

    NTC is the fixed-effect reference; only NTC + perturbations with >= ``min_cells``
    cells in that timepoint enter the model. NB the random effect has only ~6 levels
    (the states), so it is a modest baseline correction, not a well-powered variance
    component; a fixed ``C(state)`` term gives a near-identical adjustment.

    Returns (sig_by_day, tables_by_day): sig = perturbations with FDR<``fdr`` (and
    coef>0 if ``require_increase``); tables = full per-perturbation stats.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests

    state_col = state_col or CELL_STATE_COL
    obs = adata.obs
    tp = obs["time_point"].astype(str).values
    sig_by_day, tables_by_day = {}, {}
    cols = ["timepoint", "perturbation", "n_pert", "coef", "se", "z", "pval", "fdr"]
    for day in TIMEPOINT_ORDER:
        d = obs.loc[tp == day, [pert_col, score, state_col]].copy()
        d.columns = ["pert", "y", "state"]
        d["pert"] = d["pert"].astype(str)
        d["state"] = d["state"].astype(str)
        d["y"] = d["y"].astype(float)

        counts = d["pert"].value_counts()
        keep = [p for p in counts.index if p == ref or counts[p] >= min_cells]
        if ref not in keep or len(keep) < 2:
            tables_by_day[day] = pd.DataFrame(columns=cols); sig_by_day[day] = []
            continue
        d = d[d["pert"].isin(keep)].copy()
        d["pert"] = pd.Categorical(d["pert"],
                                   categories=[ref] + sorted(set(d["pert"]) - {ref}))
        try:
            mdf = smf.mixedlm("y ~ C(pert)", d, groups=d["state"]).fit(
                reml=reml, method="lbfgs")
        except Exception:
            tables_by_day[day] = pd.DataFrame(columns=cols); sig_by_day[day] = []
            continue

        rows = []
        for name in mdf.params.index:
            if not name.startswith("C(pert)[T."):
                continue
            p = name[len("C(pert)[T."):-1]
            coef, se = float(mdf.params[name]), float(mdf.bse[name])
            p2 = float(mdf.pvalues[name])                       # two-sided Wald
            pv = (p2 / 2 if coef > 0 else 1 - p2 / 2)           # one-sided: increase
            rows.append({"timepoint": day, "perturbation": p,
                         "n_pert": int((d["pert"] == p).sum()),
                         "coef": coef, "se": se, "z": coef / se if se > 0 else 0.0,
                         "pval": pv})
        r = pd.DataFrame(rows)
        r["fdr"] = multipletests(r["pval"], method="fdr_bh")[1] if len(r) else []
        r = r.sort_values("coef", ascending=False).reset_index(drop=True)
        m = r["fdr"] < fdr
        if require_increase:
            m &= r["coef"] > 0
        sig_by_day[day] = r.loc[m, "perturbation"].tolist()
        tables_by_day[day] = r
    return sig_by_day, tables_by_day


# ---------------------------------------------------------------------------
# Load a precomputed per-perturbation DE (pdex) table into beta / FDR matrices
# ---------------------------------------------------------------------------
def _bh_down_columns(P):
    """BH-FDR down each column independently (across rows), NaN-aware. P is 2-D array."""
    A = np.asarray(P, dtype=float)
    valid = ~np.isnan(A)
    nper = valid.sum(axis=0).astype(float)                 # tests per column
    order = np.argsort(np.where(valid, A, np.inf), axis=0, kind="mergesort")
    Asort = np.take_along_axis(A, order, axis=0)
    vsort = np.take_along_axis(valid, order, axis=0)
    pos = np.cumsum(vsort, axis=0).astype(float)           # rank among valid (1..nper)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = Asort * nper[None, :] / pos
    q = np.where(vsort, q, np.inf)
    q = np.minimum.accumulate(q[::-1], axis=0)[::-1]        # monotone (step-up)
    q = np.clip(q, 0, 1)
    q = np.where(vsort, q, np.nan)
    out = np.empty_like(A)
    np.put_along_axis(out, order, q, axis=0)
    return out


def bh_with_m(p, m):
    """BH step-up q-values using ``m`` as the number of tests rather than ``len(p)``.

    This is the correction behind the ``fdr`` column of the ``*_fdr8000.csv`` tables:
    BH over a fixed universe of ``m`` genes. NB it is applied to all ~23k tested
    genes, so with m < len(p) it is anti-conservative relative to a plain BH over
    what was actually tested; see the note in ``Figure5_legends.md``.

    ``_add_fdr8000.py`` keeps its own copy of this so it can run without importing
    scanpy; keep the two in step.
    """
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p, kind="mergesort")
    qs = p[order] * m / np.arange(1, n + 1)
    qs = np.clip(np.minimum.accumulate(qs[::-1])[::-1], 0, 1)   # monotone (step-up)
    q = np.empty(n)
    q[order] = qs
    return q


def load_perturbation_de(path, day, perts=None, fdr_across="table", n_genes=8000,
                         chunksize=1_000_000):
    """Read a precomputed per-perturbation DE table (pdex CSV) into beta/p/FDR matrices.

    The table has columns ``target`` (e.g. ``'NEUROG1+SIM1_day10'``), ``feature``
    (gene), ``fold_change`` (target mean / NTC mean), ``reference_mean`` (NTC mean),
    ``p_value`` and ``fdr``. **beta = log2(fold_change)** (vs NTC). Returns
    ``(beta, pval, fdr)`` DataFrames indexed by perturbation x gene.

    ``perts`` selects which perturbations to return (``None`` = all that day).
    ``fdr_across`` chooses the multiple-testing model:
      - ``'table'``  : the file's ``fdr`` (BH across ALL ~23k genes per perturbation);
      - ``'genes'``  : **BH per perturbation across the ~``n_genes`` most-expressed
        genes** (rank by NTC ``reference_mean``) -- the systematic "midway" that
        corrects over a realistic tested-gene universe (~8000) instead of 23k or 104;
        genes outside the universe get FDR = NaN;
      - ``'perturbations'`` : BH per gene across all perturbations tested (loads all).
    """
    cols = ["target", "feature", "fold_change", "reference_mean", "p_value", "fdr"]
    load_all = perts is None or fdr_across == "perturbations"
    parts = []
    if load_all:
        suffix = f"_{day}"
        for ch in pd.read_csv(path, usecols=cols, chunksize=chunksize):
            hit = ch[ch["target"].astype(str).str.endswith(suffix)]
            if len(hit):
                parts.append(hit)
        df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)
        df["pert"] = df["target"].str.replace(suffix, "", regex=False)
    else:
        want = {f"{p}_{day}": p for p in perts}
        for ch in pd.read_csv(path, usecols=cols, chunksize=chunksize):
            hit = ch[ch["target"].isin(want)]
            if len(hit):
                parts.append(hit)
        df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)
        df["pert"] = df["target"].map(want)

    df["log2fc"] = np.log2(df["fold_change"].clip(lower=1e-9))
    beta = df.pivot_table(index="pert", columns="feature", values="log2fc")
    pval = df.pivot_table(index="pert", columns="feature", values="p_value")
    if fdr_across == "perturbations":                       # BH down columns = across perturbations
        fdr = pd.DataFrame(_bh_down_columns(pval.values), index=pval.index, columns=pval.columns)
    elif fdr_across == "genes":                             # BH per perturbation over a gene universe
        ref = df.groupby("feature")["reference_mean"].mean()
        univ = (ref.sort_values(ascending=False).head(n_genes).index
                if n_genes else pval.columns)
        PU = pval.reindex(columns=univ)
        fdrU = pd.DataFrame(_bh_down_columns(PU.values.T).T, index=PU.index, columns=PU.columns)
        fdr = fdrU.reindex(columns=pval.columns)            # genes outside the universe -> NaN
    else:                                                   # "table"
        fdr = df.pivot_table(index="pert", columns="feature", values="fdr")

    if perts is not None:                                   # subset rows to requested perturbations
        keep = [p for p in perts if p in beta.index]
        beta, pval, fdr = beta.loc[keep], pval.loc[keep], fdr.loc[keep]
    return beta, pval, fdr


# ---------------------------------------------------------------------------
# Joint per-gene OLS of expression on a one-hot perturbation design (NTC = base)
# ---------------------------------------------------------------------------
def perturbation_gene_effects(adata, day, pert_col="perturbation_clean", ref=NTC_LABEL,
                              min_cells=20, min_frac=0.05, use_raw=False,
                              time_key="time_point"):
    """Joint OLS of every gene's expression on a one-hot **perturbation** design.

    On day ``day`` cells, build a single design ``X = [intercept | one-hot(perturbation)]``
    with **NTC as the base level** (its mean is the intercept), keeping perturbations
    with >= ``min_cells`` cells (cells of rarer perturbations are dropped, not pooled
    into NTC). For every gene detected in >= ``min_frac`` of the kept cells, fit
    ``expr ~ X`` in one vectorised solve (limma-style), so each perturbation's
    coefficient is its mean-expression shift vs NTC with the within-group pooled
    variance shared across all perturbations (a one-way ANOVA per gene). Two-sided
    Wald p; BH-FDR **across genes within each perturbation** (the DE convention).

    Returns ``(beta, pval, fdr)`` DataFrames indexed by perturbation (non-NTC) x gene;
    per-group cell counts are in ``beta.attrs['n']``.
    """
    import scipy.sparse as sp
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    obs = adata.obs
    daymask = obs[time_key].astype(str).values == day
    pv_all = obs[pert_col].astype(str)
    counts = pv_all[daymask].value_counts()
    kept = [p for p in counts.index if p != ref and counts[p] >= min_cells]
    cell_mask = daymask & pv_all.isin([ref] + kept).values
    sub = adata[cell_mask]
    pv = sub.obs[pert_col].astype(str)

    D = pd.get_dummies(pv).reindex(columns=kept, fill_value=0).astype(float)
    X = np.column_stack([np.ones(len(D)), D.values])      # const + one-hot perturbations
    n, p = X.shape

    Y = sub.raw.X if use_raw else sub.X
    var_names = sub.raw.var_names if use_raw else sub.var_names
    Y = sp.csc_matrix(Y)
    det = np.asarray((Y > 0).sum(axis=0)).ravel() / Y.shape[0]
    gmask = det >= min_frac
    Y = Y[:, gmask]
    genes = np.asarray(var_names)[gmask]

    XtX = X.T @ X
    XtX_inv = np.linalg.inv(XtX)
    XtY = (Y.T @ X).T                                     # (p x m)
    beta = XtX_inv @ XtY                                  # (p x m)
    yty = np.asarray(Y.multiply(Y).sum(axis=0)).ravel()
    rss = np.maximum(yty - np.einsum("pm,pm->m", beta, XtY), 0.0)
    dof = n - p
    sigma2 = rss / dof
    se = np.sqrt(np.maximum(np.outer(np.diag(XtX_inv), sigma2), 0.0))   # (p x m)
    with np.errstate(divide="ignore", invalid="ignore"):
        tval = np.where(se > 0, beta / se, 0.0)
    pval = 2 * stats.t.sf(np.abs(tval), dof)

    B = pd.DataFrame(beta[1:], index=kept, columns=genes)             # drop intercept row
    P = pd.DataFrame(pval[1:], index=kept, columns=genes)
    F = P.copy()
    for pt in kept:                                                   # BH across genes per pert
        F.loc[pt, :] = multipletests(P.loc[pt].to_numpy(), method="fdr_bh")[1]
    B.attrs["n"] = {ref: int((pv == ref).sum()), **{q: int((pv == q).sum()) for q in kept}}
    return B, P, F


# ---------------------------------------------------------------------------
# Per-combo OLS Doxo1 interaction + total effect (plain statsmodels, no state term)
# ---------------------------------------------------------------------------
def interaction_doxo1_ols(adata, combos, day, score="Doxo1program_score", states=None,
                          pert_col="perturbation_clean", time_key="time_point",
                          min_cells=10):
    """Per-combo OLS interaction + total effect on the Doxo1 score (no state correction).

    For each ``g1+g2`` combo, on day ``day`` cells within ``states`` (default = the
    non-NE compartment, intermediate + differentiated), fit two plain statsmodels
    OLS models and read the coefficients / p-values straight off them:

      decomposition:  ``Doxo1 ~ g1 + g2 + g1:g2``   (NTC + cells with only g1/g2 guides)
                      -> beta_g1, beta_g2, beta_interaction (each with its two-sided p)
      total:          ``Doxo1 ~ combo``              (NTC + the double-KO cells)
                      -> total_effect = the combo coefficient (with its two-sided p)

    The decomposition design is the saturated 2x2 {NTC, g1-only, g2-only, double}, so
    its implied total (β_g1+β_g2+β_g1:g2) equals the separate total model's coefficient;
    the total model is kept separate so its p-value is read directly from statsmodels.
    BH-FDR is applied across combos per effect.

    Returns one row per combo with beta_/pval_/fdr_ for g1, g2, interaction, total.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests

    states = list(states or (INTERMEDIATE_STATES + DIFFERENTIATED_STATES))
    obs = adata.obs
    base = ((obs[time_key].astype(str).values == day) &
            obs[CELL_STATE_COL].astype(str).isin(states).values)
    pert = obs[pert_col].astype(str).values
    y_all = obs[score].astype(float).values
    rows = []
    for combo in combos:
        if "+" not in combo:
            continue
        g1, g2 = combo.split("+")[:2]
        others = [g for g in TARGET_GENES if g not in (g1, g2)]
        rec = {"pair": combo, "gene_1": g1, "gene_2": g2, "n_double": 0, "n_total": 0,
               "beta_g1": np.nan, "pval_g1": np.nan, "beta_g2": np.nan, "pval_g2": np.nan,
               "beta_interaction": np.nan, "pval_interaction": np.nan,
               "total_effect": np.nan, "pval_total": np.nan}

        # (1) decomposition model: NTC + cells whose only target guides are g1/g2
        dmask = base & (obs[others].sum(axis=1).values == 0)
        d = pd.DataFrame({"y": y_all[dmask],
                          "g1": (obs[g1].values[dmask] == 1).astype(float),
                          "g2": (obs[g2].values[dmask] == 1).astype(float)})
        d["g1g2"] = d["g1"] * d["g2"]
        rec["n_double"] = int(d["g1g2"].sum())
        if rec["n_double"] >= min_cells:
            m = smf.ols("y ~ g1 + g2 + g1g2", data=d).fit()
            rec.update(beta_g1=float(m.params["g1"]), pval_g1=float(m.pvalues["g1"]),
                       beta_g2=float(m.params["g2"]), pval_g2=float(m.pvalues["g2"]),
                       beta_interaction=float(m.params["g1g2"]),
                       pval_interaction=float(m.pvalues["g1g2"]))
            # (2) separate total-effect model: NTC vs the double-KO group
            tmask = base & ((pert == NTC_LABEL) | (pert == combo))
            t = pd.DataFrame({"y": y_all[tmask],
                              "combo": (pert[tmask] == combo).astype(float)})
            rec["n_total"] = int(t["combo"].sum())
            if rec["n_total"] >= min_cells:
                mt = smf.ols("y ~ combo", data=t).fit()
                rec.update(total_effect=float(mt.params["combo"]),
                           pval_total=float(mt.pvalues["combo"]))
        rows.append(rec)

    df = pd.DataFrame(rows)
    for c in ("g1", "g2", "interaction", "total"):
        df[f"fdr_{c}"] = np.nan
        ok = df[f"pval_{c}"].notna()
        if ok.any():
            df.loc[ok, f"fdr_{c}"] = multipletests(df.loc[ok, f"pval_{c}"], method="fdr_bh")[1]
    return df


# ---------------------------------------------------------------------------
# Mixed-effects Doxo1 interaction model (cell state as random intercept)
# ---------------------------------------------------------------------------
def interaction_doxo1_mixed(adata, combos, day, score="Doxo1program_score",
                            pert_onehot=True, state_col=None, time_key="time_point",
                            min_cells=10, reml=True):
    """Per-combo 2-gene interaction model on the Doxo1 score, **state = random intercept**.

    For each ``g1+g2`` combo, fit (on day ``day`` cells whose only target guides are
    g1 and/or g2; NTC = reference)::

        Doxo1 ~ g1 + g2 + g1:g2 + (1 | state)

    so the single-gene **base effects** (beta_g1, beta_g2) and the **interaction**
    (beta_g1:g2) are estimated while a per-state random intercept absorbs the
    differentiation baseline (composition correction). Same shape as the Fig2f OLS
    but mixed-effects. BH-corrects the interaction p-values across combos.

    Returns a DataFrame (one row per combo) with beta_/p_ for g1, g2, interaction,
    plus ``fdr_interaction`` and ``n_double``.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests

    state_col = state_col or CELL_STATE_COL
    obs = adata.obs
    daymask = obs[time_key].astype(str).values == day
    rows = []
    for combo in combos:
        if "+" not in combo:
            continue
        g1, g2 = combo.split("+")[:2]
        others = [g for g in TARGET_GENES if g not in (g1, g2)]
        keep = daymask & (obs[others].sum(axis=1).values == 0)
        d = pd.DataFrame({
            "y": obs[score].astype(float).values[keep],
            "g1": (obs[g1].values[keep] == 1).astype(float),
            "g2": (obs[g2].values[keep] == 1).astype(float),
            "state": obs[state_col].astype(str).values[keep]})
        d["g1g2"] = d["g1"] * d["g2"]
        n_double = int(d["g1g2"].sum())
        base = {"pair": combo, "gene_1": g1, "gene_2": g2, "n_double": n_double,
                "beta_g1": np.nan, "beta_g2": np.nan, "beta_interaction": np.nan,
                "p_g1": np.nan, "p_g2": np.nan, "p_interaction": np.nan}
        if n_double >= min_cells:
            try:
                m = smf.mixedlm("y ~ g1 + g2 + g1g2", d, groups=d["state"]).fit(
                    reml=reml, method="lbfgs")
                base.update(beta_g1=float(m.params["g1"]), beta_g2=float(m.params["g2"]),
                            beta_interaction=float(m.params["g1g2"]),
                            p_g1=float(m.pvalues["g1"]), p_g2=float(m.pvalues["g2"]),
                            p_interaction=float(m.pvalues["g1g2"]))
            except Exception:
                pass
        rows.append(base)
    df = pd.DataFrame(rows)
    df["fdr_interaction"] = np.nan
    ok = df["p_interaction"].notna()
    if ok.any():
        df.loc[ok, "fdr_interaction"] = multipletests(df.loc[ok, "p_interaction"],
                                                       method="fdr_bh")[1]
    return df


def combo_total_effect_mixed(adata, combos, day, score="Doxo1program_score",
                             pert_col="perturbation_clean", state_col=None,
                             ref=NTC_LABEL, time_key="time_point", min_cells=10, reml=True):
    """Total double-KO effect per combo from a **group** mixed model, state = random intercept.

    One model on day ``day``, NTC + the combo groups::

        Doxo1 ~ C(perturbation) + (1 | state)        (NTC = reference)

    Each combo's coefficient is its **total effect** (combo vs NTC) adjusted for the
    state random intercept; reported with a two-sided Wald p-value and BH-FDR across
    combos. This is the "model each combo as its own group" total, separate from the
    gene decomposition in :func:`interaction_doxo1_mixed`.

    Returns a DataFrame: perturbation, n, total_effect, se, pval, fdr.
    """
    import statsmodels.formula.api as smf
    from statsmodels.stats.multitest import multipletests

    state_col = state_col or CELL_STATE_COL
    obs = adata.obs
    daymask = obs[time_key].astype(str).values == day
    pert = obs[pert_col].astype(str).values
    d = pd.DataFrame({"y": obs[score].astype(float).values,
                      "pert": pert,
                      "state": obs[state_col].astype(str).values})[daymask].copy()
    counts = d["pert"].value_counts()
    keep = [ref] + [c for c in combos if counts.get(c, 0) >= min_cells]
    d = d[d["pert"].isin(keep)].copy()
    if d["pert"].nunique() < 2 or (d["pert"] == ref).sum() < min_cells:
        return pd.DataFrame(columns=["perturbation", "n", "total_effect", "se", "pval", "fdr"])
    d["pert"] = pd.Categorical(d["pert"], categories=[ref] + sorted(set(keep) - {ref}))
    m = smf.mixedlm("y ~ C(pert)", d, groups=d["state"]).fit(reml=reml, method="lbfgs")
    rows = []
    for name in m.params.index:
        if not name.startswith("C(pert)[T."):
            continue
        p = name[len("C(pert)[T."):-1]
        rows.append({"perturbation": p, "n": int((d["pert"] == p).sum()),
                     "total_effect": float(m.params[name]), "se": float(m.bse[name]),
                     "pval": float(m.pvalues[name])})       # two-sided Wald
    r = pd.DataFrame(rows)
    r["fdr"] = multipletests(r["pval"], method="fdr_bh")[1] if len(r) else []
    return r


# ---------------------------------------------------------------------------
# Gene-level 2x2 interaction / epistasis model (expr ~ g1 + g2 + g1:g2)
# ---------------------------------------------------------------------------
def interaction_de_2gene(adata, g1="NEUROG1", g2="SIM1", states=None, use_raw=False,
                         min_frac=0.05, min_cells=10, adjust_time=True,
                         time_key="time_point"):
    """Per-gene 2x2 interaction model ``expr ~ g1 + g2 + g1:g2`` (NTC = reference).

    Cells kept: NTC + cells whose ONLY target guides are ``g1`` and/or ``g2``
    (every other target one-hot == 0), within ``states`` (default = the non-NE
    compartment, intermediate + differentiated, where the Doxo1 epistasis lives).
    For every gene detected in >= ``min_frac`` of those cells, OLS
    ``expr ~ const + g1 + g2 + g1:g2 (+ day)`` is fit by a SINGLE vectorised
    matrix solve over all genes (limma-style), so 10k+ genes cost one lstsq, not
    10k statsmodels fits. Expression is the log-normalised ``.X`` (or ``.raw``).

    Because the 2x2 design is saturated, the coefficients equal empirical
    group-mean differences: ``beta_g1`` = (g1-only - NTC), ``beta_g2`` =
    (g2-only - NTC), ``beta_total`` = ``beta_g1+beta_g2+beta_int`` = (double - NTC),
    and ``beta_int`` = the epistasis (double - additive expectation). Each term is
    tested two-sided and BH-FDR-corrected across genes.

    Returns a DataFrame indexed by gene with ``beta_/p_/fdr_`` for ``g1, g2, int,
    total`` plus per-group mean expression; group sizes are in ``.attrs['n']``.
    """
    import scipy.sparse as sp
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    states = list(states or (INTERMEDIATE_STATES + DIFFERENTIATED_STATES))
    others = [g for g in TARGET_GENES if g not in (g1, g2)]
    obs = adata.obs
    cell_mask = ((obs[others].sum(axis=1).values == 0) &
                 obs[CELL_STATE_COL].astype(str).isin(states).values)
    sub = adata[cell_mask]
    a1 = (sub.obs[g1].values == 1).astype(float)
    a2 = (sub.obs[g2].values == 1).astype(float)

    groups = {"NTC": (a1 == 0) & (a2 == 0), g1: (a1 == 1) & (a2 == 0),
              g2: (a2 == 1) & (a1 == 0), f"{g1}+{g2}": (a1 == 1) & (a2 == 1)}
    n = {k: int(v.sum()) for k, v in groups.items()}
    small = {k: c for k, c in n.items() if c < min_cells}
    if small:
        raise ValueError(f"groups with < {min_cells} cells: {small} (n={n}); "
                         "relax `states`/`min_cells`.")

    # Design matrix (saturated 2x2, optional additive day term).
    cols, names = [np.ones(len(a1)), a1, a2, a1 * a2], ["const", g1, g2, f"{g1}:{g2}"]
    if adjust_time:
        day = (sub.obs[time_key].astype(str).values == TIMEPOINT_ORDER[-1]).astype(float)
        if 0 < day.sum() < len(day):
            cols.append(day); names.append(TIMEPOINT_ORDER[-1])
    X = np.column_stack(cols)
    p = X.shape[1]

    # Expression matrix (cells x genes), filtered to detected genes.
    Y = sub.raw.X if use_raw else sub.X
    var_names = sub.raw.var_names if use_raw else sub.var_names
    Y = sp.csc_matrix(Y)
    det = np.asarray((Y > 0).sum(axis=0)).ravel() / Y.shape[0]
    gmask = det >= min_frac
    Y = Y[:, gmask]
    genes = np.asarray(var_names)[gmask]

    # Vectorised OLS: beta = (X'X)^-1 X'Y over all genes at once.
    XtX = X.T @ X
    XtX_inv = np.linalg.inv(XtX)
    XtY = (Y.T @ X).T                                  # (p x m)
    beta = XtX_inv @ XtY                               # (p x m)
    yty = np.asarray(Y.multiply(Y).sum(axis=0)).ravel()
    rss = np.maximum(yty - np.einsum("pm,pm->m", beta, XtY), 0.0)
    df = Y.shape[0] - p
    sigma2 = rss / df

    idx = {nm: i for i, nm in enumerate(names)}
    def term(contrast):
        c = np.asarray(contrast, float)
        est = c @ beta
        se = np.sqrt(np.maximum(sigma2 * (c @ XtX_inv @ c), 0.0))
        with np.errstate(divide="ignore", invalid="ignore"):
            tval = np.where(se > 0, est / se, 0.0)
        return est, 2 * stats.t.sf(np.abs(tval), df)

    res = pd.DataFrame(index=pd.Index(genes, name="gene"))
    contrasts = {"g1": {g1: 1}, "g2": {g2: 1}, "int": {f"{g1}:{g2}": 1},
                 "total": {g1: 1, g2: 1, f"{g1}:{g2}": 1}}
    for label, spec in contrasts.items():
        c = np.zeros(p)
        for nm, val in spec.items():
            c[idx[nm]] = val
        est, pv = term(c)
        res[f"beta_{label}"] = est
        res[f"p_{label}"] = pv
        res[f"fdr_{label}"] = multipletests(pv, method="fdr_bh")[1]
    for gname, m in groups.items():
        res[f"mean_{gname}"] = np.asarray(Y[m].mean(axis=0)).ravel()
    res.attrs["n"] = n
    res.attrs["genes2"] = (g1, g2)
    return res


def volcano_plot(ax, df, effect, fdr_thr=0.05, eff_thr=0.1, top_n=12,
                 up_color="#d62728", down_color="#1f77b4", title=None,
                 label=True, fontsize=6, highlight=None,
                 highlight_color="#111111", xlabel="Effect size (Δ log-expr vs NTC)"):
    """Volcano of one effect column of ``df`` (``beta_``/``p_``/``fdr_`` + ``effect``).

    x = effect size, y = -log10(p); genes with FDR<``fdr_thr`` and
    |effect|>``eff_thr`` are coloured by sign and the ``top_n`` most significant are
    labelled. Works both for an interaction-model term (``g1``/``g2``/``int``/
    ``total``) and for a per-perturbation DE table, where the effect is the log2
    fold-change vs NTC; pass ``xlabel`` to say which.

    ``highlight`` is an optional set of genes (e.g. the senescence panel of 4L) to
    ring and always label, so the same genes can be followed across panels.
    Returns the significant genes.
    """
    x = df[f"beta_{effect}"].values
    pv = df[f"p_{effect}"].values
    q = df[f"fdr_{effect}"].values
    y = -np.log10(np.clip(pv, 1e-300, None))
    sig = (q < fdr_thr) & (np.abs(x) >= eff_thr)
    genes = df.index.values

    ax.scatter(x[~sig], y[~sig], s=4, color="0.8", lw=0, rasterized=True)
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=9, color=up_color, lw=0)
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=9, color=down_color, lw=0)
    ax.axvline(0, color="k", lw=0.5)
    for s in (-eff_thr, eff_thr):
        ax.axvline(s, color="0.6", ls=":", lw=0.6)

    hit = np.zeros(len(genes), dtype=bool)
    if highlight:
        hit = np.isin(genes, list(highlight)) & sig
        if hit.any():
            ax.scatter(x[hit], y[hit], s=34, facecolor="none",
                       edgecolor=highlight_color, lw=1.1, zorder=4)

    if label and sig.any():
        idx = np.where(sig)[0]
        top = list(idx[np.argsort(pv[idx])[:top_n]]) + list(np.where(hit)[0])
        texts = []
        for i in dict.fromkeys(top):                 # de-duplicate, keep order
            texts.append(ax.text(x[i], y[i], genes[i], fontsize=fontsize, zorder=5,
                                 color=highlight_color if hit[i] else "black",
                                 fontweight="bold" if hit[i] else "normal",
                                 ha="left" if x[i] >= 0 else "right"))
        # The significant cloud is dense, so static offsets overlap badly. Repel the
        # labels off each other and off the points where adjustText is installed;
        # fall back to a small fixed nudge when it is not.
        try:
            from adjustText import adjust_text
            adjust_text(texts, x=x[sig], y=y[sig], ax=ax,
                        expand_points=(1.3, 1.5), expand_text=(1.2, 1.4),
                        force_text=(0.4, 0.8), force_points=(0.2, 0.4),
                        arrowprops=dict(arrowstyle="-", color="0.5", lw=0.4))
        except ImportError:
            for t in texts:
                t.set_position((t.get_position()[0], t.get_position()[1] + 0.03))
    ax.set_xlabel(xlabel)
    ax.set_ylabel("-log10 p")
    if title:
        ax.set_title(title, fontsize=10)
    return df.index[sig].tolist()


# ---------------------------------------------------------------------------
# Manuscript figure numbering
# ---------------------------------------------------------------------------
# Single source of truth mapping each notebook's internal panel name to the name
# it is published under. `savefig` consults this, so notebooks keep their own
# (stable) names and only this table changes when the numbering changes.
# Anything absent from the table is written under its own name unchanged.
#
#   Figure 4      -> main figure panels
#   FigS01-S03    -> supplementary group 1: data quality control
#   FigS04-S17    -> supplementary group 2: complementary to the presented analysis
#                    (e.g. the day04 counterpart of a day10 main panel)
PAPER_PANELS = {
    # ---- Figure 4 (main) -------------------------------------------------
    "pseudotime_a2_umap_doxo1":                          "Fig4B_umap_doxo1",
    "pseudotime_a1_backbone":                 "Fig4C_pseudotime_backbone",
    "overview_b_umap_timepoint":                       "Fig4D_umap_timepoint",
    "overview_a_umap_cellstate":                       "Fig4E_umap_cellstate",
    "pseudotime_a3_doxo1_vs_pseudotime":                 "Fig4F_doxo1_vs_pseudotime",
    "pseudotime_c_state_enrichment_heatmap_day04":       "Fig4G_state_enrichment_heatmap_day04",
    "pseudotime_c_state_enrichment_heatmap_day10":       "Fig4H_state_enrichment_heatmap_day10",
    "pseudotime_e1_doxo1_vs_pseudotime_rollcells_day10": "Fig4J_doxo1_vs_pseudotime_rollcells_day10",
    "interaction_heatmap_day10":             "Fig4K_interaction_heatmap_day10",
    "geneeffects_perturbation_gene_heatmap_day10":       "Fig4L_perturbation_gene_heatmap_day10",

    # ---- Figure 5 (main) -------------------------------------------------
    "geneeffects_volcano_NEUROG1+SIM1_day10": "Fig5F_volcano_NEUROG1_SIM1_day10",
    "geneeffects_senescence_interaction_NEUROG1+SIM1_day10": "Fig5G_senescence_interaction_NEUROG1_SIM1_day10",

    # ---- Supplementary group 1: data quality control ---------------------
    "overview_QC_depth":                   "FigS01_QC_sequencing_depth",
    "overview_QC_panels":                  "FigS02_QC_library_composition_panels",
    "overview_QC_cells_per_perturbation":  "FigS03_QC_cells_per_perturbation",

    # ---- Supplementary group 2: complementary to the presented analysis --
    "overview_c_umap_phase":                             "FigS04_Comp_umap_cellcycle_phase",
    # not published: plain pseudotime UMAP, superseded by Figure 4C.
    # None = savefig skips the file entirely.
    "pseudotime_a_umap_pseudotime":                        None,
    "overview_e_umap_state_scores":                      "FigS05_Comp_umap_state_signature_scores",
    "overview_f_marker_dotplot":                         "FigS06_Comp_state_marker_dotplot",
    "pseudotime_b_by_state":                    "FigS07_Comp_pseudotime_by_state",
    # not published (dropped from the figure set):
    "pseudotime_a4_umap_cdkn1a":                           None,
    "pseudotime_a5_cdkn1a_vs_pseudotime":                  "FigS08_Comp_CDKN1A_vs_pseudotime",
    "overview_d_composition_by_time":                    "FigS09_Comp_composition_control_vs_perturbed",
    "pseudotime_c_timepoint_state_barplots":               "FigS10_Comp_timepoint_state_composition",
    # not published (dropped from the figure set):
    "pseudotime_d_density_enriched_byday":      None,
    "pseudotime_e1_doxo1_vs_pseudotime_rollcells_day04":   "FigS11_Comp_doxo1_vs_pseudotime_day04",
    # not published (dropped from the figure set):
    "pseudotime_e1_doxo1_vs_pseudotime_fixedpt_day04":     None,
    # not published (dropped from the figure set):
    "pseudotime_e1_doxo1_vs_pseudotime_fixedpt_day10":     None,
    "pseudotime_e_doxo1_vs_pseudotime":                    "FigS16_Comp_doxo1_vs_pseudotime_wilcoxon_hits",
    # not published (dropped from the figure set):
    "pseudotime_e0_mixedlm_doxo1_state":                   None,
    "interaction_heatmap_day04":               "FigS12_Comp_interaction_heatmap_day04",
    "interaction_state_vs_doxo1_day04":                    "FigS13_Comp_state_vs_doxo1_day04",
    "interaction_state_vs_doxo1_day10":                    "FigS14_Comp_state_vs_doxo1_day10",
    "pseudotime_f_NEUROG1_SIM1_interaction":               "FigS17_Comp_NEUROG1_SIM1_interaction",
    "geneeffects_perturbation_gene_heatmap_day04":         "FigS15_Comp_perturbation_gene_heatmap_day04",
}


# ---------------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------------
def savefig(name, fig=None, formats=("pdf",), tight=True, **kw):
    """Write a figure to FIG_DIR as vector PDF (pass `formats` for other types)."""
    fig = fig or plt.gcf()
    if tight:
        fig.tight_layout()
    out = PAPER_PANELS.get(name, name)          # published name, if this panel has one
    if out is None:                             # explicitly not part of the figure set
        print(f"skipped {name} (not published)")
        return
    for ext in formats:
        fig.savefig(FIG_DIR / f"{out}.{ext}", dpi=300, bbox_inches="tight", **kw)
    label = f"{out}  [{name}]" if out != name else name
    print(f"saved {label} -> {', '.join(formats)}")
