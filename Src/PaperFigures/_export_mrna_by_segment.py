#!/usr/bin/env python
"""
_export_mrna_by_segment.py — mean RNA expression per pseudotime segment, per gene.

For the ATAC-vs-RNA concordance panels: computes the mean (normalized) expression of
every gene within each of the 16 pseudotime segments, using the SAME adaptive
even-trajectory binning as the ATAC analysis (09_PseudotimePeaks). Segments are defined
from the ATAC cells' pseudotime (so they align 1:1 with the ATAC TF-activity coefficient
columns PT_seg01..16), then applied to all RNA cells.

Output: CSV_Files/mRNA_by_pseudotime_segment.csv  (genes x PT_seg01..N, mean expression)

Run:  /home/eraslab1/miniconda3/envs/scanpy_env/bin/python _export_mrna_by_segment.py
"""
import os, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
from paperfig_style import H5AD, CSV_DIR, load_obs_by_atac

BIN_WIDTH, BIN_FLOOR = 0.02, 2000

# --- define segments from the ATAC cells' pseudotime (same as 09 / FigS5D) ---
atac = pd.read_csv(os.path.join(CSV_DIR, "ATAC_UMAP.csv"))
atac_pt = atac["cellName"].map(load_obs_by_atac("dpt_pseudotime")).dropna().values
maxb = int(np.floor(atac_pt.max() / BIN_WIDTH))
b = np.minimum((atac_pt // BIN_WIDTH).astype(int), maxb)
cnt = np.bincount(b, minlength=maxb + 1)
newid = np.zeros(maxb + 1, int); cur = 1; acc = 0
for i in range(maxb + 1):
    newid[i] = cur; acc += cnt[i]
    if acc >= BIN_FLOOR: cur += 1; acc = 0
if acc > 0 and cur > 1: newid[newid == cur] = cur - 1
nseg = int(newid.max())
print(f"segments: {nseg}")

# --- load the RNA matrix + assign each RNA cell to a segment via its pseudotime ---
with h5py.File(H5AD, "r") as f:
    shape = tuple(f["X"].attrs["shape"])
    X = csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=shape)
    genes = [g.decode() if isinstance(g, bytes) else g for g in f["var/gene_symbols"][:]]
    pt = np.asarray(f["obs/dpt_pseudotime"][:], float)
print(f"RNA matrix: {shape[0]} cells x {shape[1]} genes")

valid = np.isfinite(pt)
seg = np.zeros(len(pt), int)
seg[valid] = newid[np.minimum((pt[valid] // BIN_WIDTH).astype(int), maxb)]  # 1..nseg; 0 = drop

# --- mean expression per segment = (onehot^T @ X) / counts, via sparse matmul ---
keep = seg > 0
rows = np.where(keep)[0]
onehot = csr_matrix((np.ones(len(rows)), (rows, seg[rows] - 1)), shape=(shape[0], nseg))
seg_sum = np.asarray((onehot.T @ X).todense())          # nseg x genes
seg_n = np.asarray(onehot.sum(axis=0)).ravel()          # cells per segment
mean = seg_sum / seg_n[:, None]

out = pd.DataFrame(mean.T, index=genes, columns=[f"PT_seg{n:02d}" for n in range(1, nseg + 1)])
out.index.name = "gene"
path = os.path.join(CSV_DIR, "mRNA_by_pseudotime_segment.csv")
out.to_csv(path)
print(f"cells/segment: {seg_n.astype(int)}")
print(f"wrote {path}  ({out.shape[0]} genes x {out.shape[1]} segments)")
