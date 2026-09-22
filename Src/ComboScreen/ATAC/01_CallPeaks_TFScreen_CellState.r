#!/usr/bin/env Rscript
# 01_CallPeaks_TFScreen_CellState.r
# -----------------------------------------------------------------------------
# Sibling of 01_CallPeaks_TFScreen.r, but calls reproducible peaks grouped by
# CELL STATE instead of perturbation.
#
# Cell state is read from ComboScreen_processed.h5ad obs. Two equivalent 6-way
# partitions are available (same cells, different naming):
#     final_label : Neuroendocrine / Intermediate-1..3 / Differentiated_1 / Differentiated-2
#     state       : NEPC-N / NEPC-A&N-1..3 / NEPC-A-1 / NEPC-A-2
# Default is final_label; switch GROUP_COL to "state", or to "current_leiden_1.5"
# (21 Leiden clusters) for a finer grouping.
#
# Cells are keyed onto ArchR's ATAC cell names by the shared-barcode string
# transform of the multimodal/RNA cell name (verified 100% coverage):
#     250804_2lig_10G.ACGG...  ->  TFScreen#250804_2lig_ATAC_10G_ACGG...-1
#
# Prerequisite: the merged master arrow (all 24 chromosomes) from 03_merge_arrows.r.
# Output goes to a SEPARATE project dir so it does not clobber the perturbation run.
#
# Usage:  Rscript Src/01_CallPeaks_TFScreen_CellState.r
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ArchR)
  library(rhdf5)
})

# ---- paths / params ---------------------------------------------------------
DATA_DIR     <- "/data/AbbasTFScreen"
MASTER_ARROW <- file.path(DATA_DIR, "TFScreen.arrow")               # from 03_merge_arrows.r
PROJ_DIR     <- file.path(DATA_DIR, "TFScreenATAC_cellstate")       # SEPARATE from perturbation run
H5AD         <- "/home/eraslab1/Projects/AbbasTFScreen/ComboScreen_processed.h5ad"

# Cell-state column in the h5ad obs: "final_label" (6 states, default),
# "state" (same 6, NEPC naming), or "current_leiden_1.5" (21 clusters).
GROUP_COL           <- "final_label"
MIN_CELLS_PER_GROUP <- 100        # drop states with fewer cells than this
THREADS             <- 90

setwd(DATA_DIR)
addArchRThreads(threads = THREADS)
addArchRGenome("hg38")

# ---- helper: read an anndata obs column (categorical or plain) from .h5ad ----
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
       "\n  -> run Src/parallel_arrows/03_merge_arrows.r first.")
}
stopifnot(h5read(MASTER_ARROW, "Class") == "Arrow")
ls_arrow <- h5ls(MASTER_ARROW); h5closeAll()
fchr <- ls_arrow$name[ls_arrow$group == "/Fragments"]
message("Master arrow: ", length(fchr), " chromosomes -> ", paste(sort(fchr), collapse = ", "))
if (length(fchr) < 24) {
  stop("Master arrow has only ", length(fchr), " chromosomes; re-run 03_merge_arrows.r.")
}

# ---- 1. ArchRProject from the single merged arrow ---------------------------
proj <- ArchRProject(
  ArrowFiles      = MASTER_ARROW,
  outputDirectory = PROJ_DIR,
  copyArrows      = FALSE,
  showLogo        = FALSE
)
message("Project created: ", nCells(proj), " cells")

# ---- 2. attach per-cell CELL STATE (+ covariates) from the h5ad -------------
mm <- as.character(h5read(H5AD, "obs/multimodal_cell_names"))
rna_to_atac_cell <- function(x) {
  left <- sub("\\.[^.]*$", "", x)
  seq  <- sub("^[^.]*\\.", "", x)
  samp <- sub("_[^_]*$", "", left)
  well <- sub("^.*_", "", left)
  paste0("TFScreen#", samp, "_ATAC_", well, "_", seq, "-1")
}
meta <- data.frame(
  cellName  = rna_to_atac_cell(mm),
  cellState = read_obs(H5AD, GROUP_COL),
  stringsAsFactors = FALSE
)
for (col in c("perturbation_clean", "phase", "time_point", "conditions")) {
  meta[[col]] <- tryCatch(read_obs(H5AD, col), error = function(e) NA_character_)
}
h5closeAll()

m <- match(proj$cellNames, meta$cellName)
message("Cells matched to h5ad metadata: ", sum(!is.na(m)), " / ", nCells(proj))

# ArchRProject is S4: only `$<-` works (not `[[<-`), so assign each column explicitly.
proj$cellState          <- meta$cellState[m]
proj$perturbation_clean <- meta$perturbation_clean[m]
proj$phase              <- meta$phase[m]
proj$time_point         <- meta$time_point[m]
proj$conditions         <- meta$conditions[m]

# keep only cells with a usable cell-state label
proj <- proj[!is.na(proj$cellState) & nzchar(as.character(proj$cellState)), ]
message("Cells with a cell-state label: ", nCells(proj))

# ---- 3. drop states with too few cells (all 6 pass at 100; matters for leiden) --
tab  <- table(proj$cellState)
keep <- names(tab)[tab >= MIN_CELLS_PER_GROUP]
message("States kept (>= ", MIN_CELLS_PER_GROUP, " cells): ", length(keep), " / ", length(tab),
        "  [cells kept: ", sum(tab[keep]), "]")
proj <- proj[as.character(proj$cellState) %in% keep, ]

# ---- 4. group coverages -> reproducible peak set -> peak matrix --------------
# States are large (thousands of cells each), so replicates form easily; raise
# reproducibility if you want stricter peaks.
proj <- addGroupCoverages(
  ArchRProj     = proj,
  groupBy       = "cellState",
  minCells      = 40,
  maxCells      = 1000,
  minReplicates = 3,
  maxReplicates = 6,
  threads       = THREADS
)

pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
  ArchRProj       = proj,
  groupBy         = "cellState",
  reproducibility = "2",
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
write.csv(peaks, file.path(PROJ_DIR, "PeakSet_cellstate.csv"), row.names = FALSE)

message("DONE. Peaks: ", nrow(peaks),
        " | states: ", length(keep),
        " | project saved to ", PROJ_DIR)
