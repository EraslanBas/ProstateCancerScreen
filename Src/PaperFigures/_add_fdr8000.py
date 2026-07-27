"""Add an 8000-test-corrected FDR to the pdex DE tables.

For each perturbation (target), recompute BH-FDR on its gene p-values but with the
number of tests fixed at M=8000 (a realistic tested-gene universe), instead of the
~23k genes used in the original `fdr`. Writes new files `Day{04,10}DEGs_fdr8000.csv`
with an added column `fdr_8000` and the `fdr` column overwritten with these values;
the originals are left untouched.
"""
import numpy as np
import pandas as pd

M = 8000
SRC = "/home/eraslab1/Projects/AbbasScreen/Src/ComboScreen"


def bh_with_m(p, m):
    """BH step-up q-values using `m` as the number of tests (not len(p))."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p, kind="mergesort")
    qs = p[order] * m / np.arange(1, n + 1)
    qs = np.minimum.accumulate(qs[::-1])[::-1]      # enforce monotone (step-up)
    qs = np.clip(qs, 0, 1)
    q = np.empty(n)
    q[order] = qs
    return q


for day in ["04", "10"]:
    src = f"{SRC}/Day{day}DEGs.csv"
    df = pd.read_csv(src)
    print(f"Day{day}: read {df.shape}, {df['target'].nunique()} perturbations")
    df["fdr_8000"] = df.groupby("target")["p_value"].transform(lambda s: bh_with_m(s.values, M))
    df["fdr"] = df["fdr_8000"]
    out = f"{SRC}/Day{day}DEGs_fdr8000.csv"
    df.to_csv(out, index=False)
    # quick CDKN1A sanity
    cd = df[df.feature == "CDKN1A"]
    print(f"  wrote {out} | CDKN1A: fdr<0.1 in {(cd.fdr < 0.1).sum()} / {len(cd)} perturbations")
