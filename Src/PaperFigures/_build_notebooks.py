#!/usr/bin/env python
"""
_build_notebooks.py — (re)generate the Figure 5 panel notebooks for AbbasTFScreen.

Each panel is one notebook that imports paperfig_style, loads its data, draws ONE
panel, and writes panels/<name>.pdf + .png. This builder emits clean, fully
editable .ipynb files (kernel: scanpy_env) so you can open and tweak any panel
directly; rerun this script only if you want to regenerate them from scratch.

    /home/eraslab1/miniconda3/envs/scanpy_env/bin/python _build_notebooks.py          # build
    /home/eraslab1/miniconda3/envs/scanpy_env/bin/python _build_notebooks.py --run    # build + execute
"""
import sys, os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell

HERE = os.path.dirname(os.path.abspath(__file__))
KERNEL = "scanpy_env"

HDR = "from paperfig_style import *\nimport numpy as np, pandas as pd, matplotlib.pyplot as plt\n"

# ============================================================================ PANELS
PANELS = {}
def panel(name, *cells): PANELS[name] = cells

# ---------------------------------------------------------------- Fig 5A  ATAC UMAP x cell state
panel("Fig5A_ATAC_UMAP_cellstate",
("md", "# Fig 5A — ATAC-profile UMAP colored by (RNA-defined) cell state\n"
       "UMAP built from the **ATAC chromatin profile alone** (IterativeLSI on a genome-wide "
       "500 bp TileMatrix — features independent of the RNA labels), colored by the RNA-derived "
       "`cellState`. If cells of the same state group together here, the two modalities are "
       "**concordant**. Concordance is quantified by the Adjusted Rand Index (ARI) between the "
       "unsupervised ATAC clusters and the RNA states."),
("code", HDR),
("code",
"from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score\n"
"df = pd.read_csv(os.path.join(CSV_DIR, 'ATAC_UMAP.csv'))\n"
"# ATAC clusters: use the precomputed column shipped in ATAC_UMAP.csv; if absent (i.e. running\n"
"# from the full LSI on the source VM), Leiden-cluster ATAC_LSI.csv (162 MB, not in the repo).\n"
"if 'ATAC_cluster' not in df.columns:\n"
"    import scanpy as sc, anndata as ad\n"
"    lsi = pd.read_csv(os.path.join(CSV_DIR, 'ATAC_LSI.csv'), index_col=0).loc[df['cellName']]\n"
"    A = ad.AnnData(np.zeros((len(df), 1), dtype='float32')); A.obsm['X_lsi'] = lsi.values.astype('float32')\n"
"    sc.pp.neighbors(A, use_rep='X_lsi', n_neighbors=15)\n"
"    sc.tl.leiden(A, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)\n"
"    df['ATAC_cluster'] = A.obs['leiden'].values\n"
"df['ATAC_cluster'] = df['ATAC_cluster'].astype(str)\n"
"mask = df['cellState'].notna() & (df['cellState'] != 'NA')\n"
"ari = adjusted_rand_score(df['cellState'][mask], df['ATAC_cluster'][mask])\n"
"nmi = normalized_mutual_info_score(df['cellState'][mask], df['ATAC_cluster'][mask])\n"
"df = df[mask]\n"
"print('ATAC cells:', len(df), '| Leiden clusters:', df['ATAC_cluster'].nunique(),\n"
"      '| ARI(state,ATACcluster)=%.2f  NMI=%.2f' % (ari, nmi))\n"
"fig, ax = plt.subplots(figsize=(3.4, 3.2))\n"
"for s in STATE_ORDER:\n"
"    m = df['cellState'] == s\n"
"    ax.scatter(df['UMAP1'][m], df['UMAP2'][m], s=1.4, c=STATE_COLORS[s],\n"
"               label=state_label(s), linewidths=0, rasterized=True)\n"
"ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('ATAC-UMAP1'); ax.set_ylabel('ATAC-UMAP2')\n"
"for sp in ['left','bottom']: ax.spines[sp].set_visible(False)\n"
"ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), markerscale=4,\n"
"          handletextpad=0.3, labelspacing=0.35, borderpad=0)\n"
"savepanel(fig, 'Fig5A_ATAC_UMAP_cellstate')\n"),
)

# ---------------------------------------------------------------- Fig 5B  ATAC UMAP x pseudotime
panel("Fig5B_ATAC_UMAP_pseudotime",
("md", "# Fig 5B — ATAC-profile UMAP colored by (RNA-defined) pseudotime\n"
       "Same ATAC-only UMAP, colored by the RNA `dpt_pseudotime` mapped onto ATAC cells via the "
       "shared multiome barcode. The black→white line is the **RNA main trajectory drawn in "
       "chromatin space** — the running-median ATAC-UMAP position across 80 pseudotime bins, "
       "root (NE, dot) → tip (Differentiated, arrow). A smooth gradient with the trajectory "
       "tracking it means the differentiation path is **also visible in chromatin**. Spearman "
       "correlation of pseudotime with ATAC-UMAP position summarizes the gradient."),
("code", HDR),
("code",
"from scipy.stats import spearmanr\n"
"df = pd.read_csv(os.path.join(CSV_DIR, 'ATAC_UMAP.csv'))\n"
"pt_map = load_obs_by_atac('dpt_pseudotime')          # {atac_cellName: pseudotime}\n"
"df['pt'] = df['cellName'].map(pt_map)\n"
"df = df[df['pt'].notna()]\n"
"print('cells with pseudotime:', len(df))\n"
"# orient by whichever UMAP axis best tracks pseudotime (just for the reported number)\n"
"rho = max(abs(spearmanr(df['pt'], df['UMAP1'])[0]), abs(spearmanr(df['pt'], df['UMAP2'])[0]))\n"
"fig, ax = plt.subplots(figsize=(3.6, 3.2))\n"
"o = np.argsort(df['pt'].values)                      # draw high pt on top\n"
"sc = ax.scatter(df['UMAP1'].values[o], df['UMAP2'].values[o], c=df['pt'].values[o],\n"
"                cmap=PSEUDOTIME_CMAP, s=1.4, linewidths=0, rasterized=True)\n"
"# main RNA trajectory drawn in chromatin space: running-median ATAC-UMAP position\n"
"# across 80 RNA-pseudotime bins, connected root(NE) -> tip(Differentiated)\n"
"u = df[['UMAP1','UMAP2']].values; ptv = df['pt'].values\n"
"NB = 80\n"
"oo = np.argsort(ptv, kind='mergesort'); tb = np.empty(len(ptv), int)\n"
"tb[oo] = (np.arange(len(ptv)) * NB) // len(ptv)\n"
"cx = np.array([np.median(u[tb==b, 0]) for b in range(NB)])\n"
"cy = np.array([np.median(u[tb==b, 1]) for b in range(NB)])\n"
"rm = lambda a, w=5: np.array([np.median(a[max(0,i-w//2):i+w//2+1]) for i in range(len(a))])\n"
"sx, sy = rm(cx), rm(cy)\n"
"ax.plot(sx, sy, color='black', lw=2.6, solid_capstyle='round', zorder=5)\n"
"ax.plot(sx, sy, color='white', lw=0.9, solid_capstyle='round', zorder=6)\n"
"ax.scatter(sx[0], sy[0], s=34, color='black', edgecolor='white', lw=0.8, zorder=7)  # root\n"
"ax.annotate('', xy=(sx[-1], sy[-1]), xytext=(sx[-4], sy[-4]),\n"
"            arrowprops=dict(arrowstyle='-|>', color='black', lw=2), zorder=7)         # direction\n"
"cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02); cb.set_label('Pseudotime based on RNA profile')\n"
"cb.outline.set_linewidth(0.5)\n"
"ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('ATAC-UMAP1'); ax.set_ylabel('ATAC-UMAP2')\n"
"for sp in ['left','bottom']: ax.spines[sp].set_visible(False)\n"
"savepanel(fig, 'Fig5B_ATAC_UMAP_pseudotime')\n"),
)

# ---------------------------------------------------------------- Fig 5C  state TF motif accessibility
panel("Fig5C_StateTFActivity_heatmap",
("md", "# Fig 5C — TF motif accessibility per cell state\n"
       "Faithful reproduction of **08_PlotTFActivity Part A**. OLS coefficient of per-peak Log2FC "
       "(each state vs the rest) on the grouped HOCOMOCO motif matrix. Rows = ALL TFs significant "
       "(FDR≤0.1) in ≥1 state, ordered by the state where they peak (then strongest-first within a "
       "state); columns = states in NE→Differentiated order. `*` = FDR≤0.1; symmetric 3-stop "
       "blue-white-red scale with breaks ±max|coef|."),
("code", HDR),
("code",
"# --- 08_PlotTFActivity Part A (cells 3-4), replicated ------------------------\n"
"coef = load_matrix('StateTFActivity_coef.csv')\n"
"fdr  = load_matrix('StateTFActivity_FDR.csv')\n"
"cols = [s for s in STATE_ORDER if s in coef.columns]\n"
"coef, fdr = coef[cols], fdr[cols]\n"
"sigTF = coef.index[(fdr <= FDR_THRESHOLD).sum(axis=1) >= 1]        # sig in >=1 state (no cap)\n"
"mat = coef.loc[sigTF]; star = (fdr.loc[sigTF] <= FDR_THRESHOLD).values\n"
"# order rows by WHERE they peak (max.col), then strongest first within that state\n"
"peak    = mat.values.argmax(axis=1)\n"
"peakval = mat.values[np.arange(len(mat)), peak]\n"
"order   = np.lexsort((-peakval, peak))\n"
"mat = mat.iloc[order]; star = star[order]\n"
"print(len(sigTF), 'sig TFs |  states:', ' -> '.join(cols))\n"
"lim = np.abs(mat.values).max()\n"
"fig, ax = plt.subplots(figsize=(3.2, max(4.0, 0.16*len(mat))), layout='constrained')\n"
"im = ax.imshow(mat.values, aspect='auto', cmap=ACTIVITY_CMAP, vmin=-lim, vmax=lim)\n"
"ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels([state_label(c) for c in mat.columns], rotation=90)\n"
"ax.set_yticks(range(len(mat))); ax.set_yticklabels(tf_labels(mat.index))\n"
"for i in range(mat.shape[0]):\n"
"    for j in range(mat.shape[1]):\n"
"        if star[i, j]: ax.text(j, i, '*', ha='center', va='center', fontsize=8, color='black')\n"
"ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]\n"
"# horizontal colorbar at the bottom (fits the tall heatmap better than a side bar)\n"
"cb = fig.colorbar(im, ax=ax, location='bottom', fraction=0.03, pad=0.05, aspect=35)\n"
"cb.set_label('TF motif accessibility (coef)'); cb.outline.set_linewidth(0.5)\n"
"ax.set_title('TF motif accessibility per cell state\\n(NE \\u2192 Differentiated;  * FDR\\u22640.1)', fontsize=8)\n"
"savepanel(fig, 'Fig5C_StateTFActivity_heatmap')\n"),
)

# ---------------------------------------------------------------- Supp C  pseudotime TF motif accessibility (full heatmap)
panel("FigS5E_PseudotimeTFActivity_heatmap",
("md", "# Fig S5E — TF motif accessibility across pseudotime (full heatmap)\n"
       "Full heatmap companion to main-figure **Fig 5D**; based on **08_PlotTFActivity Part C**. "
       "Per pseudotime segment (vs the NE root `PT_seg01`), OLS coef of peak Log2FC on the grouped "
       "motif matrix. Rows = ALL TFs significant (FDR≤0.1) in **≥1 segment**, ordered by Spearman "
       "trend across segments — **rising at top, falling at bottom**; columns run root NE → "
       "Differentiated tip. `*` = FDR≤0.1."),
("code", HDR),
("code",
"from scipy.stats import spearmanr\n"
"# --- 08_PlotTFActivity Part C (cells 10-11), replicated ----------------------\n"
"pcoef = load_matrix('PseudotimeTFActivity_coef.csv')\n"
"pfdr  = load_matrix('PseudotimeTFActivity_FDR.csv')\n"
"seg = sorted(pcoef.columns, key=lambda c: int(c.replace('PT_seg', '')))\n"
"pcoef, pfdr = pcoef[seg], pfdr[seg]\n"
"MIN_SIG = 1   # supplementary full heatmap: keep TFs significant in >=1 segment\n"
"sigTF = pcoef.index[(pfdr <= FDR_THRESHOLD).sum(axis=1) >= MIN_SIG]\n"
"pmat = pcoef.loc[sigTF]; pstar = (pfdr.loc[sigTF] <= FDR_THRESHOLD).values\n"
"idx = np.arange(pmat.shape[1])\n"
"trend = np.array([spearmanr(idx, pmat.values[i])[0] for i in range(len(pmat))])\n"
"order = np.argsort(-trend)                                         # rising -> falling\n"
"pmat = pmat.iloc[order]; pstar = pstar[order]\n"
"print(len(sigTF), 'sig TFs across', pmat.shape[1], 'pseudotime segments (root NE -> Diff tip)')\n"
"lim = np.abs(pmat.values).max()\n"
"fig, ax = plt.subplots(figsize=(3.8, max(4.0, 0.14*len(pmat))), layout='constrained')\n"
"im = ax.imshow(pmat.values, aspect='auto', cmap=ACTIVITY_CMAP, vmin=-lim, vmax=lim)\n"
"ax.set_xticks(range(0, pmat.shape[1], 2)); ax.set_xticklabels(range(1, pmat.shape[1]+1, 2))\n"
"ax.set_xlabel('pseudotime segment (root NE \\u2192 Differentiated tip)')\n"
"ax.set_yticks(range(len(pmat))); ax.set_yticklabels(tf_labels(pmat.index))\n"
"for i in range(pmat.shape[0]):\n"
"    for j in range(pmat.shape[1]):\n"
"        if pstar[i, j]: ax.text(j, i, '*', ha='center', va='center', fontsize=7, color='black')\n"
"ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]\n"
"# horizontal colorbar at the bottom (fits the tall heatmap better than a side bar)\n"
"cb = fig.colorbar(im, ax=ax, location='bottom', fraction=0.03, pad=0.06, aspect=35)\n"
"cb.set_label('TF motif accessibility (coef)'); cb.outline.set_linewidth(0.5)\n"
"ax.set_title('TF motif accessibility across pseudotime\\n(root NE \\u2192 Differentiated tip;  * FDR\\u22640.1)', fontsize=8)\n"
"savepanel(fig, 'FigS5E_PseudotimeTFActivity_heatmap')\n"),
)

# ---------------------------------------------------------------- Supp F  perturbation TF motif accessibility
panel("FigS5F_PertTFActivity_heatmap",
("md", "# Fig S5F — TF motif accessibility per perturbation (pooled vs NTC)\n"
       "Faithful reproduction of **08_PlotTFActivity Part D** (heatmap). Each perturbation vs "
       "pooled NTC, no state stratification, on the **pseudotime peak set** (114k peaks) → "
       "getMarkerFeatures Log2FC → OLS on the grouped motifs. Shown: the nonzero FDR-significant "
       "coefficient sub-matrix (TFs × perturbations), rows/cols hierarchically clustered, each "
       "cell labelled with its coefficient."),
("code", HDR),
("code",
"from scipy.cluster.hierarchy import linkage, leaves_list\n"
"# --- 08_PlotTFActivity Part D (cell 15), replicated (pseudotime-peak perturbation result) ---\n"
"pcf = load_matrix('PertTFActivity_ptpeaks_coefFDRsig.csv')\n"
"pcf = pcf.loc[(pcf != 0).any(axis=1), (pcf != 0).any(axis=0)]\n"
"print('nonzero sub-matrix:', pcf.shape, '(TF x perturbation)')\n"
"def order(df, axis):                                    # pheatmap default: complete linkage\n"
"    X = df.values if axis == 0 else df.values.T\n"
"    if X.shape[0] < 3: return np.arange(X.shape[0])\n"
"    return leaves_list(linkage(X, method='complete'))\n"
"pcf = pcf.iloc[order(pcf, 0), order(pcf, 1)]\n"
"lim = np.abs(pcf.values).max()\n"
"fig, ax = plt.subplots(figsize=(max(3.2, 0.5*pcf.shape[1]+1.5), max(2.6, 0.3*pcf.shape[0]+1)))\n"
"im = ax.imshow(pcf.values, aspect='auto', cmap=ACTIVITY_CMAP, vmin=-lim, vmax=lim)\n"
"ax.set_xticks(range(pcf.shape[1])); ax.set_xticklabels(pcf.columns, rotation=90)\n"
"ax.set_yticks(range(pcf.shape[0])); ax.set_yticklabels(tf_labels(pcf.index))\n"
"for i in range(pcf.shape[0]):\n"
"    for j in range(pcf.shape[1]):\n"
"        v = pcf.values[i, j]\n"
"        if v != 0: ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=6, color='black')\n"
"ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]\n"
"cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03); cb.set_label('TF motif accessibility (coef)')\n"
"cb.outline.set_linewidth(0.5)\n"
"ax.set_title('Perturbation vs NTC (pseudotime peaks)\\nFDR-significant TF motif accessibility', fontsize=8)\n"
"savepanel(fig, 'FigS5F_PertTFActivity_heatmap')\n"),
)

# ---------------------------------------------------------------- Supp G  state x pert calls
panel("FigS5G_StateXPert_TFActivity",
("md", "# Fig S5G — TF motif accessibility per (state × perturbation)\n"
       "The stratified analysis (08_PlotTFActivity Part B): each perturbation vs same-state NTC. "
       "Nonzero FDR-significant coefficient sub-matrix (TF × condition). Columns are **grouped by "
       "cell state** (trajectory order; perturbation within), with a colored strip and state "
       "labels above and separators between blocks; the x-axis shows only the perturbation. TF "
       "rows hierarchically clustered. Symmetric 3-stop scale."),
("code", HDR),
("code",
"import itertools\n"
"from scipy.cluster.hierarchy import linkage, leaves_list\n"
"from matplotlib.patches import Rectangle\n"
"cf = load_matrix('TFActivity_coefFDRsig_stateXpert.csv')\n"
"cf = cf.loc[(cf != 0).any(axis=1), (cf != 0).any(axis=0)]        # all conditions with >=1 sig TF\n"
"cs = [c.split('___')[0] for c in cf.columns]\n"
"present = [s for s in STATE_ORDER if s in set(cs)]\n"
"# columns grouped by state (trajectory order); within each state, cluster the conditions\n"
"col_ord = []\n"
"for s in present:\n"
"    idx = [j for j in range(cf.shape[1]) if cs[j] == s]\n"
"    if len(idx) > 2:\n"
"        sub = leaves_list(linkage(cf.iloc[:, idx].values.T, method='complete'))\n"
"        idx = [idx[k] for k in sub]\n"
"    col_ord.extend(idx)\n"
"cf = cf.iloc[:, col_ord]\n"
"cs   = [c.split('___')[0] for c in cf.columns]\n"
"pert = [c.split('___')[1] for c in cf.columns]\n"
"ri = leaves_list(linkage(cf.values, method='complete')) if cf.shape[0] >= 3 else np.arange(cf.shape[0])\n"
"cf = cf.iloc[ri]\n"
"n = cf.shape[1]; lim = np.abs(cf.values).max()\n"
"print('nonzero sub-matrix:', cf.shape, '| states:', present)\n"
"fig = plt.figure(figsize=(max(5.5, 0.17*n + 1.8), max(3.6, 0.14*cf.shape[0] + 1.4)))\n"
"gs = fig.add_gridspec(2, 1, height_ratios=[1, 24], hspace=0.04)\n"
"axs = fig.add_subplot(gs[0]); axh = fig.add_subplot(gs[1])\n"
"im = axh.imshow(cf.values, aspect='auto', cmap=ACTIVITY_CMAP, vmin=-lim, vmax=lim)\n"
"axh.set_xticks(range(n)); axh.set_xticklabels(pert, rotation=90, fontsize=5)\n"
"axh.set_yticks(range(cf.shape[0])); axh.set_yticklabels(tf_labels(cf.index), fontsize=5)\n"
"axh.set_xlim(-0.5, n-0.5); axh.tick_params(length=0)\n"
"[s.set_visible(False) for s in axh.spines.values()]\n"
"# ---- state annotation strip above (NEPC display names) ----\n"
"axs.set_xlim(-0.5, n-0.5); axs.set_ylim(0, 1); axs.axis('off')\n"
"start = 0\n"
"for s, grp in itertools.groupby(cs):\n"
"    cnt = sum(1 for _ in grp); center = start + cnt/2 - 0.5\n"
"    axs.add_patch(Rectangle((start-0.5, 0), cnt, 1, color=STATE_COLORS.get(s, 'grey'), ec='white', lw=0.8))\n"
"    axs.text(center, 1.7, state_label(s), ha='center', va='bottom', fontsize=5, rotation=0)\n"
"    if start > 0: axh.axvline(start-0.5, color='white', lw=1.6)\n"
"    start += cnt\n"
"# legend: color -> NEPC state name (reliable reference for narrow blocks)\n"
"from matplotlib.patches import Patch\n"
"handles = [Patch(color=STATE_COLORS[s], label=state_label(s)) for s in present]\n"
"axh.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.14, 1.0), fontsize=5.5,\n"
"           handlelength=1, borderpad=0.3, title='cell state', title_fontsize=6)\n"
"cb = fig.colorbar(im, ax=[axs, axh], fraction=0.03, pad=0.02); cb.set_label('coef (FDR≤0.1)')\n"
"cb.outline.set_linewidth(0.5)\n"
"savepanel(fig, 'FigS5G_StateXPert_TFActivity')\n"),
)

# ---------------------------------------------------------------- Supp D  pt segments UMAP
panel("FigS5D_PseudotimeSegments_UMAP",
("md", "# Fig S5D — Pseudotime segments: ATAC vs RNA embedding\n"
       "The 16 adaptive even-trajectory pseudotime segments (equal-width DPT steps of 0.02 merged "
       "forward to ≥2,000 cells — the segments used to re-call peaks in 09_PseudotimePeaks) shown "
       "on the **ATAC-profile UMAP (left)** and the **RNA UMAP (right)**. Segments are defined "
       "once from pseudotime and colored identically in both panels, so a concordant layout shows "
       "the trajectory bins occupy matching territories in chromatin and transcriptome."),
("code", HDR),
("code",
"BIN_WIDTH, BIN_FLOOR = 0.02, 2000\n"
"# --- ATAC cells: ATAC-UMAP coords + RNA pseudotime mapped by shared barcode ---\n"
"atac = pd.read_csv(os.path.join(CSV_DIR, 'ATAC_UMAP.csv'))\n"
"atac['pt'] = atac['cellName'].map(load_obs_by_atac('dpt_pseudotime'))\n"
"atac = atac[atac['pt'].notna()].copy()\n"
"# --- define segments once from ATAC pseudotime (same scheme as 09_PseudotimePeaks) ---\n"
"pt = atac['pt'].values\n"
"maxb = int(np.floor(pt.max()/BIN_WIDTH))\n"
"b = np.minimum((pt//BIN_WIDTH).astype(int), maxb)\n"
"cnt = np.bincount(b, minlength=maxb+1)\n"
"newid = np.zeros(maxb+1, int); cur=1; acc=0\n"
"for i in range(maxb+1):\n"
"    newid[i]=cur; acc+=cnt[i]\n"
"    if acc>=BIN_FLOOR: cur+=1; acc=0\n"
"if acc>0 and cur>1: newid[newid==cur]=cur-1\n"
"atac['seg'] = newid[b]\n"
"nseg = int(atac['seg'].max())\n"
"# --- RNA cells: multiome UMAP + pseudotime, SAME pt->segment mapping ---\n"
"rna_umap, obs = load_umap_obs(('dpt_pseudotime',))\n"
"rna_pt = np.asarray(obs['dpt_pseudotime'], float)\n"
"valid = np.isfinite(rna_pt)\n"
"rb = np.minimum((rna_pt[valid]//BIN_WIDTH).astype(int), maxb)\n"
"rna_seg = newid[rb]; rna_u = rna_umap[valid]\n"
"print('segments:', nseg, '| ATAC cells:', len(atac), '| RNA cells:', int(valid.sum()))\n"
"cmap = plt.get_cmap('turbo', nseg)\n"
"fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.4))\n"
"o = np.argsort(atac['seg'].values)\n"
"axs[0].scatter(atac['UMAP1'].values[o], atac['UMAP2'].values[o], c=atac['seg'].values[o],\n"
"               cmap=cmap, vmin=0.5, vmax=nseg+0.5, s=1.2, linewidths=0, rasterized=True)\n"
"axs[0].set_title(f'ATAC UMAP \\u2014 {nseg} pseudotime segments', fontsize=8)\n"
"axs[0].set_xlabel('ATAC-UMAP1'); axs[0].set_ylabel('ATAC-UMAP2')\n"
"o2 = np.argsort(rna_seg)\n"
"sc = axs[1].scatter(rna_u[o2,0], rna_u[o2,1], c=rna_seg[o2],\n"
"                    cmap=cmap, vmin=0.5, vmax=nseg+0.5, s=1.2, linewidths=0, rasterized=True)\n"
"axs[1].set_title(f'RNA UMAP \\u2014 {nseg} pseudotime segments', fontsize=8)\n"
"axs[1].set_xlabel('RNA-UMAP1'); axs[1].set_ylabel('RNA-UMAP2')\n"
"cb = fig.colorbar(sc, ax=axs, fraction=0.03, pad=0.02); cb.set_label('pseudotime segment (1 = NE root)')\n"
"cb.outline.set_linewidth(0.5)\n"
"for ax in axs:\n"
"    ax.set_xticks([]); ax.set_yticks([])\n"
"    for sp in ['left','bottom']: ax.spines[sp].set_visible(False)\n"
"savepanel(fig, 'FigS5D_PseudotimeSegments_UMAP')\n"),
)

# ---------------------------------------------------------------- Fig 5D  pseudotime TF trend lines (main)
panel("Fig5D_PseudotimeTF_lines",
("md", "# Fig 5D — TF motif accessibility across pseudotime (rising vs falling)\n"
       "Based on **08_PlotTFActivity Part C** (line view). Coefficient trajectory across "
       "pseudotime segments (root NE=1 → Differentiated tip) of the top 10 **rising** (solid) and "
       "top 10 **falling** (dashed) TFs by Spearman trend. Coefficient = TF motif accessibility of each "
       "segment vs. the NE root; the full heatmap of all significant TFs is Fig S5E."),
("code", HDR),
("code",
"# --- 08_PlotTFActivity Part C line view (cell 12), replicated ----------------\n"
"pcoef = load_matrix('PseudotimeTFActivity_coef.csv')\n"
"tr = pd.read_csv(os.path.join(CSV_DIR, 'PseudotimeTFActivity_trend.csv'))\n"
"tr = tr[tr['nSigSeg'] >= 1]\n"
"seg = sorted(pcoef.columns, key=lambda c: int(c.replace('PT_seg', '')))\n"
"pcoef = pcoef[seg]; x = np.arange(1, len(seg)+1)\n"
"riseTF = [g for g in tr.sort_values('trend', ascending=False)['TF'] if g in pcoef.index][:10]\n"
"fallTF = [g for g in tr.sort_values('trend')['TF'] if g in pcoef.index][:10]\n"
"fig, ax = plt.subplots(figsize=(6.0, 3.4))\n"
"ax.axhline(0, color='grey', lw=0.6)\n"
"cmap = plt.get_cmap('tab20')\n"
"for i, g in enumerate(riseTF):\n"
"    ax.plot(x, pcoef.loc[g].values, '-',  lw=1.0, color=cmap(i % 20), label=f'{tf_label(g)} (rise)')\n"
"    ax.plot(x, pcoef.loc[g].values, '.',  ms=3, color=cmap(i % 20))\n"
"for i, g in enumerate(fallTF):\n"
"    ax.plot(x, pcoef.loc[g].values, '--', lw=1.0, color=cmap((i+10) % 20), label=f'{tf_label(g)} (fall)')\n"
"    ax.plot(x, pcoef.loc[g].values, '.',  ms=3, color=cmap((i+10) % 20))\n"
"ax.set_xticks(x)\n"
"ax.set_xlabel('pseudotime segment (root NE = 1 \\u2192 Differentiated tip)')\n"
"ax.set_ylabel('Relative motif accessibility compared to the root')\n"
"ax.set_title('Top TFs changing across pseudotime', fontsize=8)\n"
"ax.legend(fontsize=4.5, ncol=2, loc='center left', bbox_to_anchor=(1.005, 0.5),\n"
"          handlelength=1.4, labelspacing=0.25)\n"
"savepanel(fig, 'Fig5D_PseudotimeTF_lines')\n"),
)

# ---------------------------------------------------------------- Fig 5E  ATAC-vs-RNA examples (main)
panel("Fig5E_ATACvsRNA_examples",
("md", "# Fig 5E — TF motif accessibility vs RNA expression (examples)\n"
       "Example TFs from the pseudotime trajectory: two with **rising** TF motif accessibility (ZEB1, GRHL2) "
       "and two **falling** (NFIC, ZBTB14). Each point is a pseudotime segment (colored NE root → "
       "differentiated tip, same palette as Fig S5D); x = TF motif accessibility coefficient (vs NE root), "
       "y = mean mRNA. GRHL2 is concordant (activity and expression rise together); ZEB1/NFIC/ZBTB14 "
       "are discordant. Full set in Fig S5H."),
("code", HDR),
("code",
"from scipy.stats import spearmanr\n"
"coef = load_matrix('PseudotimeTFActivity_coef.csv')\n"
"seg = sorted([c for c in coef.columns if c.startswith('PT_seg')], key=lambda c: int(c.replace('PT_seg','')))\n"
"coef = coef[seg]\n"
"mrna = load_matrix('mRNA_by_pseudotime_segment.csv'); mrna = mrna[[c for c in seg if c in mrna.columns]]\n"
"mrna = mrna[~mrna.index.duplicated(keep='first')]\n"
"MAIN = [('ZEB1','rise'), ('GRHL2','rise'), ('NFIC','fall'), ('ZBT14','fall')]  # coef rownames\n"
"segidx = np.arange(1, len(seg)+1); cmap = segment_cmap(len(seg))\n"
"fig, axs = plt.subplots(2, 2, figsize=(3.7, 3.5), layout='constrained')\n"
"axs = axs.ravel()\n"
"for ax, (t, d) in zip(axs, MAIN):\n"
"    g = tf_label(t); x = coef.loc[t].values.astype(float); y = mrna.loc[g].values.astype(float)\n"
"    rho = spearmanr(x, y)[0]\n"
"    sc = ax.scatter(x, y, c=segidx, cmap=cmap, vmin=0.5, vmax=len(seg)+0.5, s=13, edgecolor='k', linewidths=0.25)\n"
"    ax.set_title(f'{g}  ({d}, \\u03c1={rho:.2f})', fontsize=6.5)\n"
"    ax.axvline(0, color='grey', lw=0.4, ls='--')\n"
"    ax.set_xlabel('TF motif accessibility (vs root)', fontsize=5.5); ax.tick_params(labelsize=5)\n"
"for ax in (axs[0], axs[2]): ax.set_ylabel('mean mRNA', fontsize=5.5)\n"
"cb = fig.colorbar(sc, ax=axs.tolist(), fraction=0.03, pad=0.02)\n"
"cb.set_label('pseudotime segment (1 = NE root)', fontsize=5.5); cb.ax.tick_params(labelsize=5)\n"
"savepanel(fig, 'Fig5E_ATACvsRNA_examples')\n"),
)

# ---------------------------------------------------------------- Supp H  TF motif accessibility vs RNA (rest)
panel("FigS5H_TFActivity_vs_mRNA",
("md", "# Fig S5H — TF motif accessibility vs RNA expression across pseudotime\n"
       "For each Fig 5D TF with a single gene (TFGroup_N motif groups excluded), a scatter of its "
       "TF motif accessibility coefficient (vs NE root) against its mean mRNA level, one point per "
       "pseudotime segment (colored NE root → differentiated tip). Spearman ρ per TF quantifies "
       "chromatin–transcriptome concordance. Concordant TFs (e.g. GRHL2, ρ≈+0.9) track together; "
       "shared-motif TFs (e.g. ZEB1) diverge because motif activity reflects a whole E-box family, "
       "not the single gene's expression."),
("code", HDR),
("code",
"from scipy.stats import spearmanr\n"
"coef = load_matrix('PseudotimeTFActivity_coef.csv')\n"
"seg = sorted([c for c in coef.columns if c.startswith('PT_seg')], key=lambda c: int(c.replace('PT_seg','')))\n"
"coef = coef[seg]\n"
"tr = pd.read_csv(os.path.join(CSV_DIR, 'PseudotimeTFActivity_trend.csv')); tr = tr[tr['nSigSeg'] >= 1]\n"
"rise = [g for g in tr.sort_values('trend', ascending=False)['TF'] if g in coef.index][:10]\n"
"fall = [g for g in tr.sort_values('trend')['TF'] if g in coef.index][:10]\n"
"mrna = load_matrix('mRNA_by_pseudotime_segment.csv'); mrna = mrna[[c for c in seg if c in mrna.columns]]\n"
"mrna = mrna[~mrna.index.duplicated(keep='first')]\n"
"# individual TFs only (motif groups discarded); alias the few HOCOMOCO names -> gene symbols\n"
"ALIAS = {'ZN317':'ZNF317','ZBT14':'ZBTB14','NDF1':'NEUROD1','NGN1':'NEUROG1','TWST1':'TWIST1',\n"
"         'GCR':'NR3C1','ANDR':'AR','PRGR':'PGR','NF2L2':'NFE2L2','THA':'THRA','THB':'THRB'}\n"
"def to_gene(t):\n"
"    if t.startswith('TFGroup_'): return None       # discard motif groups\n"
"    if t in mrna.index: return t\n"
"    return ALIAS.get(t) if ALIAS.get(t) in mrna.index else None\n"
"MAIN = {'ZEB1','GRHL2','NFIC','ZBT14'}   # shown in main Fig 5E; the rest go here\n"
"pairs = [(t, to_gene(t)) for t in rise + fall if t not in MAIN]\n"
"tfs = [(t, g) for t, g in pairs if g is not None]\n"
"dropped = [t for t, g in pairs if g is None and not t.startswith('TFGroup_')]\n"
"print('individual TFs plotted:', len(tfs), '|', [g for _, g in tfs])\n"
"print('dropped (absent from RNA reference):', dropped)\n"
"segidx = np.arange(1, len(seg)+1); cmap = segment_cmap(len(seg))\n"
"ncol = 3; nrow = int(np.ceil(len(tfs)/ncol))\n"
"fig, axs = plt.subplots(nrow, ncol, figsize=(2.2*ncol, 2.1*nrow), layout='constrained')\n"
"axs = np.atleast_1d(axs).ravel()\n"
"for k, (t, g) in enumerate(tfs):\n"
"    ax = axs[k]; x = coef.loc[t].values.astype(float); y = mrna.loc[g].values.astype(float)\n"
"    rho = spearmanr(x, y)[0]\n"
"    sc = ax.scatter(x, y, c=segidx, cmap=cmap, vmin=0.5, vmax=len(seg)+0.5, s=20, edgecolor='k', linewidths=0.3)\n"
"    d = 'rise' if t in rise else 'fall'\n"
"    ax.set_title(f'{g}  ({d}, \\u03c1={rho:.2f})', fontsize=7)\n"
"    ax.axvline(0, color='grey', lw=0.4, ls='--')\n"
"    ax.set_xlabel('TF motif accessibility (coef vs root)', fontsize=6)\n"
"    ax.set_ylabel('mean mRNA', fontsize=6); ax.tick_params(labelsize=5)\n"
"for k in range(len(tfs), len(axs)): axs[k].axis('off')\n"
"cb = fig.colorbar(sc, ax=axs.tolist(), fraction=0.015, pad=0.02); cb.set_label('pseudotime segment', fontsize=6)\n"
"fig.suptitle('TF motif accessibility vs RNA expression across pseudotime', fontsize=8)\n"
"savepanel(fig, 'FigS5H_TFActivity_vs_mRNA')\n"),
)

def build():
    made = []
    for name, cells in PANELS.items():
        nb = new_notebook()
        nb.cells = [new_markdown_cell(s) if k == "md" else new_code_cell(s) for k, s in cells]
        nb.metadata["kernelspec"] = {"name": KERNEL, "display_name": KERNEL, "language": "python"}
        nb.metadata["language_info"] = {"name": "python"}
        path = os.path.join(HERE, name + ".ipynb")
        nbf.write(nb, path)
        made.append(path)
        print("wrote", os.path.basename(path))
    return made

def run(paths):
    from nbclient import NotebookClient
    for p in paths:
        print("executing", os.path.basename(p), "...", flush=True)
        nb = nbf.read(p, as_version=4)
        NotebookClient(nb, timeout=1200, kernel_name=KERNEL, resources={"metadata": {"path": HERE}}).execute()
        nbf.write(nb, p)
    print("all notebooks executed")

if __name__ == "__main__":
    paths = build()
    if "--run" in sys.argv:
        run(paths)
