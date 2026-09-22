"""Build ``Data/ComboScreen.h5ad``, the object every downstream step starts from.

This is the first step this repository can run. Its two inputs are produced by an
external pipeline that is NOT part of this repository (see README, "External
inputs"):

  * ``NEPC_sec_screen_RNA_scanpy_raw_counts.h5ad``
        raw UMI counts, all cells, all genes.
  * ``NEPC_merged_sec_screen_RNA_scanpy_final_singlets.HVG_downstream.h5ad``
        the same cells after the upstream QC / clustering / scoring, carrying the
        35 ``obs`` columns the figures rely on: the 14 guide one-hots, ``time_point``,
        ``final_label``, the six state signature scores plus
        ``Doxo-1-differentiated_score``, ``phase``, ``UMI_counts``, ``gene_number``.

What it does (unchanged from the original, which lived commented out in cell 1 of
``01_Data_outlook.ipynb``): read the raw counts, drop genes seen in fewer than
``MIN_CELLS`` cells, copy every ``obs`` column across from the downstream object,
and write the result.

Expected output: 449,267 cells x 23,423 genes, 35 obs columns.

Usage::

    python 00_Build_ComboScreen.py
    python 00_Build_ComboScreen.py --data-dir /path/to/Data
"""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

# Genes must be detected in at least this many cells to be kept.
MIN_CELLS = 1000

RAW_COUNTS = "NEPC_sec_screen_RNA_scanpy_raw_counts.h5ad"
HVG_DOWNSTREAM = "NEPC_merged_sec_screen_RNA_scanpy_final_singlets.HVG_downstream.h5ad"
OUTPUT = "ComboScreen.h5ad"

DEFAULT_DATA_DIR = Path(__file__).resolve().parents[2] / "Data"


def build(data_dir: Path, overwrite: bool = False) -> Path:
    raw_fn, meta_fn = data_dir / RAW_COUNTS, data_dir / HVG_DOWNSTREAM
    out_fn = data_dir / OUTPUT

    for fn in (raw_fn, meta_fn):
        if not fn.exists():
            raise FileNotFoundError(
                f"{fn} not found. This file comes from the external upstream "
                "pipeline; see the README section 'External inputs'."
            )
    if out_fn.exists() and not overwrite:
        raise FileExistsError(f"{out_fn} already exists; pass --overwrite to rebuild.")

    print(f"reading {raw_fn.name}")
    adata = sc.read_h5ad(raw_fn)
    print(f"  raw: {adata.shape}")

    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    print(f"  after filter_genes(min_cells={MIN_CELLS}): {adata.shape}")

    print(f"reading {meta_fn.name}")
    adata_1 = sc.read_h5ad(meta_fn)

    # The copy below aligns on the obs index, so cells missing from adata_1 would be
    # filled with NaN rather than raising. Report the overlap so that stays visible.
    shared = adata.obs_names.intersection(adata_1.obs_names)
    print(f"  obs overlap: {len(shared)} / {adata.n_obs} cells")
    if len(shared) < adata.n_obs:
        print(f"  WARNING: {adata.n_obs - len(shared)} cells absent from "
              f"{meta_fn.name}; their metadata will be NaN.")

    for elem in adata_1.obs.columns:
        adata.obs[elem] = adata_1.obs[elem]
    print(f"  copied {len(adata_1.obs.columns)} obs columns")

    adata.write(out_fn)
    print(f"wrote {out_fn}  ({adata.shape[0]} cells x {adata.shape[1]} genes, "
          f"{adata.obs.shape[1]} obs columns)")
    return out_fn


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR,
                    help=f"directory holding the inputs (default: {DEFAULT_DATA_DIR})")
    ap.add_argument("--overwrite", action="store_true",
                    help="rebuild even if ComboScreen.h5ad already exists")
    args = ap.parse_args()
    build(args.data_dir, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
