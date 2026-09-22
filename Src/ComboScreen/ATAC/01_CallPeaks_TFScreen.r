#!/usr/bin/env Rscript
# 01_CallPeaks_TFScreen.r
# -----------------------------------------------------------------------------
# Adapted from 01_CallPeaks.r (Osaka Omics project) for the AbbasTFScreen project.
#
# What changed vs the original:
#   * Builds the ArchRProject from THIS project's arrow(s): the single merged
#     master arrow assembled from the 24 per-chromosome arrows by
#     Src/parallel_arrows/03_merge_arrows.r  (the per-chromosome arrows cannot be
#     used directly — every cell would live in 24 arrows; see that script's header).
#   * Attaches per-cell PERTURBATION labels read directly from
#     ComboScreen_processed.h5ad (obs/perturbation_clean), mapped onto ArchR's
#     ATAC cell names by a pure string transform of the RNA/multimodal cell name
#     (the 2-level-ligation barcode is shared across modalities; verified 100% of
#     names map into the valid ATAC barcode set):
#         250804_2lig_10G.ACGG...  ->  TFScreen#250804_2lig_ATAC_10G_ACGG...-1
#   * Calls reproducible peaks grouped by perturbation (the screen's TF/guide),
#     the analogue of the original's "compoundDosage".
#
# Prerequisite: the merged master arrow must exist and contain all 24 chromosomes.
#   Rscript Src/parallel_arrows/03_merge_arrows.r \
#       /data/AbbasTFScreen/parallel_build/build \
#       /data/AbbasTFScreen/TFScreen.arrow \
#       chr1,chr2,...,chr22,chrX,chrY 24
#
# Usage:  Rscript Src/01_CallPeaks_TFScreen.r
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ArchR)
  library(rhdf5)
})

# ---- paths / params ---------------------------------------------------------
# NOTE: Conf.R sets atacFilesDir="/data/abbastfscreen", which does not exist on
# this machine; the arrows live under /data/AbbasTFScreen. Paths are set
# explicitly here so the script is self-contained.
DATA_DIR     <- "/data/AbbasTFScreen"
MASTER_ARROW <- file.path(DATA_DIR, "TFScreen.arrow")            # from 03_merge_arrows.r
PROJ_DIR     <- file.path(DATA_DIR, "TFScreenATAC")              # ArchRProject output
H5AD         <- "/home/eraslab1/Projects/AbbasTFScreen/ComboScreen_processed.h5ad"

# Grouping column in the h5ad obs. "perturbation_clean" (92 groups, NTC-paired
# singletons collapsed) is the analysis-ready grouping; switch to "perturbation"
# (105 raw groups) if you want the raw labels.
GROUP_COL           <- "perturbation_clean"
MIN_CELLS_PER_GROUP <- 100        # drop perturbations with fewer cells than this
THREADS             <- 90

setwd(DATA_DIR)
addArchRThreads(threads = THREADS)
addArchRGenome("hg38")

# ---- helpers: read an anndata obs column (categorical or plain) from .h5ad ----
read_obs <- function(h5ad, col) {
  base <- paste0("obs/", col)
  cats <- tryCatch(h5read(h5ad, paste0(base, "/categories")), error = function(e) NULL)
  if (!is.null(cats)) {                     # anndata categorical: categories + 0-based codes
    codes <- as.integer(h5read(h5ad, paste0(base, "/codes")))
    out <- rep(NA_character_, length(codes))
    ok  <- codes >= 0                       # -1 == NaN
    out[ok] <- as.character(cats)[codes[ok] + 1L]
    return(out)
  }
  as.vector(h5read(h5ad, base))             # plain dataset
}

# ---- 0. sanity-check the merged master arrow (must hold all 24 chromosomes) --
if (!file.exists(MASTER_ARROW)) {
  stop("Master arrow not found: ", MASTER_ARROW,
       "\n  -> run Src/parallel_arrows/03_merge_arrows.r first (see header).")
}
stopifnot(h5read(MASTER_ARROW, "Class") == "Arrow")
ls_arrow <- h5ls(MASTER_ARROW); h5closeAll()
fchr <- ls_arrow$name[ls_arrow$group == "/Fragments"]
message("Master arrow: ", length(fchr), " chromosomes -> ", paste(sort(fchr), collapse = ", "))
if (length(fchr) < 24) {
  stop("Master arrow has only ", length(fchr), " chromosomes; re-run ",
       "03_merge_arrows.r after all 24 per-chromosome arrows exist.")
}

# ---- 1. ArchRProject from the single merged arrow ---------------------------
proj <- ArchRProject(
  ArrowFiles      = MASTER_ARROW,
  outputDirectory = PROJ_DIR,
  copyArrows      = FALSE,
  showLogo        = FALSE
)
message("Project created: ", nCells(proj), " cells")

# ---- 2. attach per-cell perturbation (+ covariates) from the h5ad -----------
# obs index (multimodal_cell_names) -> ArchR ATAC cell name via shared-barcode transform.
mm <- as.character(h5read(H5AD, "obs/multimodal_cell_names"))     # 250804_2lig_10G.ACGG...
rna_to_atac_cell <- function(x) {
  left <- sub("\\.[^.]*$", "", x)     # 250804_2lig_10G
  seq  <- sub("^[^.]*\\.", "", x)     # ACGGAGGCGGCTCTGGAGAC
  samp <- sub("_[^_]*$", "", left)    # 250804_2lig
  well <- sub("^.*_", "", left)       # 10G
  paste0("TFScreen#", samp, "_ATAC_", well, "_", seq, "-1")
}
meta <- data.frame(
  cellName     = rna_to_atac_cell(mm),
  perturbation = read_obs(H5AD, GROUP_COL),
  stringsAsFactors = FALSE
)
for (col in c("final_label", "phase", "time_point", "conditions")) {
  meta[[col]] <- tryCatch(read_obs(H5AD, col), error = function(e) NA_character_)
}
h5closeAll()

m <- match(proj$cellNames, meta$cellName)                        # metadata row per project cell
message("Cells matched to h5ad metadata: ", sum(!is.na(m)), " / ", nCells(proj))

# ArchRProject is S4: only `$<-` works (not `[[<-`), so assign each column explicitly.
proj$perturbation <- meta$perturbation[m]
proj$final_label  <- meta$final_label[m]
proj$phase        <- meta$phase[m]
proj$time_point   <- meta$time_point[m]
proj$conditions   <- meta$conditions[m]

# keep only cells with a usable perturbation label
proj <- proj[!is.na(proj$perturbation) & nzchar(as.character(proj$perturbation)), ]
message("Cells with a perturbation label: ", nCells(proj))

# ---- 3. drop perturbations with too few cells to form pseudobulk replicates --
tab  <- table(proj$perturbation)
keep <- names(tab)[tab >= MIN_CELLS_PER_GROUP]
message("Perturbations kept (>= ", MIN_CELLS_PER_GROUP, " cells): ",
        length(keep), " / ", length(tab),
        "  [cells kept: ", sum(tab[keep]), "]")
proj <- proj[as.character(proj$perturbation) %in% keep, ]

# ---- 4. group coverages -> reproducible peak set -> peak matrix --------------
# Pseudobulk replicates per perturbation. Params are permissive to accommodate
# smaller perturbation groups (the original used minReplicates=5/reproducibility=5
# for a few large dose groups). Tune for your group sizes.
proj <- addGroupCoverages(
  ArchRProj     = proj,
  groupBy       = "perturbation",
  minCells      = 40,
  maxCells      = 500,
  minReplicates = 3,
  maxReplicates = 5,
  threads       = THREADS
)

pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
  ArchRProj       = proj,
  groupBy         = "perturbation",
  reproducibility = "2",          # a peak must recur in >=2 pseudobulk replicates
  maxPeaks        = 300000,
  minCells        = 40,
  cutOff          = 0.01,
  method          = "q",
  threads         = THREADS,
  pathToMacs2     = pathToMacs2
)

proj <- addPeakMatrix(proj, threads = THREADS)

# ---- 5. save + export the peak set ------------------------------------------
saveArchRProject(ArchRProj = proj, outputDirectory = PROJ_DIR, load = FALSE)

peaks <- data.frame(getPeakSet(proj))
write.csv(peaks, file.path(PROJ_DIR, "PeakSet_perturbation.csv"), row.names = FALSE)

message("DONE. Peaks: ", nrow(peaks),
        " | groups: ", length(keep),
        " | project saved to ", PROJ_DIR)
