#!/usr/bin/env Rscript
# _compute_atac_umap.R — ATAC-only embedding for the cross-modal concordance panels.
# -----------------------------------------------------------------------------
# Builds a UMAP from the ATAC chromatin profile alone, using a GENOME-WIDE
# TileMatrix (500 bp bins) so the features are INDEPENDENT of the RNA-derived
# cell-state labels — making "do ATAC and RNA agree?" a fair test (not circular,
# as it would be if we embedded on the state-called PeakMatrix).
#
# Steps: addTileMatrix -> addIterativeLSI(TileMatrix) -> addClusters -> addUMAP.
# Exports per-cell ATAC UMAP coords + RNA cell state + ATAC cluster to CSV; the
# RNA pseudotime is joined later in Python from the h5ad.
#
# Usage:  Rscript Src/PaperFigures/_compute_atac_umap.R
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ArchR) })
addArchRThreads(threads = 32)
addArchRGenome("hg38")
set.seed(1)

PROJ <- "/data/AbbasTFScreen/TFScreenATAC_cellstate"
OUT  <- "/home/eraslab1/Projects/AbbasTFScreen/CSV_Files/ATAC_UMAP.csv"

proj <- loadArchRProject(PROJ, showLogo = FALSE)
message(">> loaded project: ", nCells(proj), " cells; matrices: ",
        paste(getAvailableMatrices(proj), collapse = ", "))

# ---- 1. genome-wide 500 bp TileMatrix (unbiased features) --------------------
if (!("TileMatrix" %in% getAvailableMatrices(proj))) {
  message(">> addTileMatrix (500 bp) ...")
  proj <- addTileMatrix(proj, tileSize = 500, binarize = TRUE, force = TRUE)
}
message(">> matrices now: ", paste(getAvailableMatrices(proj), collapse = ", "))

# ---- 2. IterativeLSI on the TileMatrix --------------------------------------
# iterations = 1: single-iteration LSI on the top accessible tiles. ArchR only
# runs its internal clustering step BETWEEN iterations (it calls Seurat::FindClusters,
# which is NOT installed here), so iterations=1 skips clustering entirely and avoids
# that dependency. Clustering for the concordance metric is done in Python (Leiden).
message(">> addIterativeLSI (TileMatrix, iterations=1; clustering done in Python) ...")
proj <- addIterativeLSI(
  ArchRProj   = proj,
  useMatrix   = "TileMatrix",
  name        = "IterativeLSI",
  iterations  = 1,
  varFeatures = 25000,
  dimsToUse   = 1:30,
  force = TRUE
)

# ---- 3. UMAP from the ATAC LSI (no clustering needed) ------------------------
message(">> addUMAP ...")
proj <- addUMAP(proj, reducedDims = "IterativeLSI", name = "UMAP",
                nNeighbors = 30, minDist = 0.3, metric = "cosine", force = TRUE)

saveArchRProject(ArchRProj = proj, outputDirectory = PROJ, load = FALSE)

# ---- 4. export coords + labels + the LSI matrix -----------------------------
um <- getEmbedding(proj, embedding = "UMAP", returnDF = TRUE)
colnames(um) <- c("UMAP1", "UMAP2")
cc <- getCellColData(proj)
df <- data.frame(
  cellName  = rownames(um),
  UMAP1     = um$UMAP1,
  UMAP2     = um$UMAP2,
  cellState = as.character(cc[rownames(um), "cellState"]),
  stringsAsFactors = FALSE
)
write.csv(df, OUT, row.names = FALSE)

# LSI reduced dims -> Python does Leiden clustering on this for the ARI concordance metric
lsi <- getReducedDims(proj, reducedDims = "IterativeLSI")           # cells x dims
lsi <- lsi[df$cellName, , drop = FALSE]
write.csv(data.frame(cellName = rownames(lsi), lsi, check.names = FALSE),
          sub("ATAC_UMAP.csv$", "ATAC_LSI.csv", OUT), row.names = FALSE)
message(">> DONE. wrote ", OUT, " (", nrow(df), " cells) + ATAC_LSI.csv (",
        ncol(lsi), " LSI dims)")
