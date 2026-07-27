from __future__ import annotations

from libraries import *
from parameters import *

import os
from pathlib import Path
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from cnmf import cNMF



def _set_thread_env(blas_threads: int) -> None:
    # Limit BLAS/OpenMP threads per worker (critical for performance)
    os.environ["OMP_NUM_THREADS"] = str(blas_threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(blas_threads)
    os.environ["MKL_NUM_THREADS"] = str(blas_threads)
    os.environ["VECLIB_MAXIMUM_THREADS"] = str(blas_threads)
    os.environ["NUMEXPR_NUM_THREADS"] = str(blas_threads)


def _factorize_one(args: tuple[str, str, int, int, int]) -> int:
    out_dir, name, worker_i, total_workers, blas_threads = args
    _set_thread_env(blas_threads)
    cnmf_obj = cNMF(output_dir=out_dir, name=name)
    cnmf_obj.factorize(worker_i=worker_i, total_workers=total_workers)
    return worker_i


def run_cnmf_parallel(
    adata: ad.AnnData,
    out_dir: str | Path,
    name: str,
    k_list=(8, 10, 12, 14),
    k_final: int = 10,
    n_iter: int = 100,
    seed: int = 14,
    density_threshold: float = 0.01,
    n_workers: int = 32,
    blas_threads_per_worker: int = 1,
    programs_csv_prefix: str | Path | None = None,
) -> dict:
    """
    Faster cNMF on a single large node by parallelizing factorize() across workers.

    Recommended settings for 128 cores:
      - n_workers: 32 or 64
      - blas_threads_per_worker: 1
    """

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if "counts" not in adata.layers:
        raise KeyError(f"adata.layers['counts'] not found. Available layers: {list(adata.layers.keys())}")

    X = adata.layers["counts"]

    # cNMF requires non-negative counts
    if sp.issparse(X):
        if (X.data < 0).any():
            raise ValueError("adata.layers['counts'] contains negative values; cNMF requires non-negative counts.")
    else:
        X_arr = np.asarray(X)
        if (X_arr < 0).any():
            raise ValueError("adata.layers['counts'] contains negative values; cNMF requires non-negative counts.")

    # Minimal AnnData for cNMF: counts in .X, keep obs_names/var_names
    # Avoid heavy copies; this is usually enough for cNMF I/O
    adata_counts = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=adata.obs_names),
        var=pd.DataFrame(index=adata.var_names),
    )

    counts_fn = out_dir / f"{name}.counts.h5ad"
    adata_counts.write_h5ad(counts_fn)

    # Prepare (single process)
    _set_thread_env(blas_threads_per_worker)  # also limit threads in main proc
    cnmf_obj = cNMF(output_dir=str(out_dir), name=name)
    cnmf_obj.prepare(
        counts_fn=str(counts_fn),
        components=np.array(list(k_list), dtype=int),
        n_iter=int(n_iter),
        seed=int(seed),
    )

    # Factorize in parallel (multiple processes)
    n_workers = int(n_workers)
    tasks = [(str(out_dir), name, i, n_workers, int(blas_threads_per_worker)) for i in range(n_workers)]

    # Use spawn for safety on some systems; fork is faster on Linux but can inherit bad state
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=n_workers) as pool:
        for wid in pool.imap_unordered(_factorize_one, tasks, chunksize=1):
            pass  # could print progress if you want

    # Combine + consensus (single process)
    cnmf_obj = cNMF(output_dir=str(out_dir), name=name)
    cnmf_obj.combine()
    cnmf_obj.k_selection_plot()
    cnmf_obj.consensus(k=int(k_final), density_threshold=float(density_threshold))

    usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
        K=int(k_final),
        density_threshold=float(density_threshold),
    )

    # Label outputs
    cells = adata.obs_names.to_list()
    genes = adata.var_names.to_list()
    programs = [f"{name}_k{k_final}_p{i+1}" for i in range(int(k_final))]

    usage_df = pd.DataFrame(usage, index=cells, columns=programs)
    spectra_tpm_df = pd.DataFrame(spectra_tpm, index=programs, columns=genes)
    spectra_scores_df = pd.DataFrame(spectra_scores, index=programs, columns=genes)

    # Save programs
    if programs_csv_prefix is None:
        programs_csv_prefix = out_dir / f"{name}_k{k_final}_programs"
    else:
        programs_csv_prefix = Path(programs_csv_prefix)

    programs_tpm_csv = programs_csv_prefix.with_suffix("")  # strip suffix if any
    programs_tpm_csv = Path(str(programs_tpm_csv) + "_tpm.csv")
    programs_scores_csv = Path(str(programs_csv_prefix.with_suffix("")) + "_scores.csv")

    spectra_tpm_df.to_csv(programs_tpm_csv)
    spectra_scores_df.to_csv(programs_scores_csv)

    return {
        "usage_df": usage_df,
        "spectra_tpm_df": spectra_tpm_df,
        "spectra_scores_df": spectra_scores_df,
        "top_genes": top_genes,
        "paths": {
            "counts_fn": str(counts_fn),
            "out_dir": str(out_dir),
            "programs_tpm_csv": str(programs_tpm_csv),
            "programs_scores_csv": str(programs_scores_csv),
        },
    }


def main():
    os.getcwd()
    os.chdir(projectDir)

    adata = sc.read("./Data/adata_20K.h5ad")


    results = run_cnmf_parallel(
        adata,
        out_dir="./cnmf_out",
        name="my_run",
        k_list=(8, 10, 12),
        k_final=10,
        n_iter=100,
        n_workers=64,
        blas_threads_per_worker=1,
    )
    print(results["paths"])

if __name__ == "__main__":
    main()