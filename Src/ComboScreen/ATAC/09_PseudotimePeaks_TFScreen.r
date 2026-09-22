#!/usr/bin/env Rscript
# 09_PseudotimePeaks_TFScreen.r
# -----------------------------------------------------------------------------
# Second exercise: call peaks on PSEUDOTIME SEGMENTS (not cell states).
# Bin cells into N_SEG quantile pseudotime segments and call reproducible peaks
# grouped by segment -> a trajectory-defined peak set + PeakMatrix, written to a
# NEW project dir so the cell-state project is untouched.
#
# Uses the SAME DPT trajectory stored in ComboScreen_processed.h5ad
# (obs/dpt_pseudotime), mapped onto ATAC cells by the verified barcode transform.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ArchR); library(rhdf5) })

DATA_DIR <- "/data/AbbasTFScreen"
MASTER   <- file.path(DATA_DIR, "TFScreen.arrow")
PROJ_DIR <- file.path(DATA_DIR, "TFScreenATAC_pseudotime")
H5AD     <- "/home/eraslab1/Projects/AbbasTFScreen/ComboScreen_processed.h5ad"
# Adaptive EVEN-TRAJECTORY binning: equal-width pseudotime steps of BIN_WIDTH,
# adjacent bins merged forward until each has >= BIN_FLOOR cells. This tiles the
# trajectory evenly (fine bins across the populated NE->transition range, sparse
# differentiated tip merged) instead of the quantile scheme that piled ~half the
# bins into the NE root. Yields ~16 bins, all >= ~2,300 cells.
BIN_WIDTH <- 0.02
BIN_FLOOR <- 2000
THREADS   <- 40

setwd(DATA_DIR)
addArchRThreads(THREADS); addArchRGenome("hg38")

# ---- ArchRProject from the master arrow -------------------------------------
proj <- ArchRProject(ArrowFiles = MASTER, outputDirectory = PROJ_DIR,
                     copyArrows = FALSE, showLogo = FALSE)
message("Project cells: ", nCells(proj))

# ---- attach DPT pseudotime (h5ad obs -> ATAC cells) -------------------------
mm <- as.character(h5read(H5AD, "obs/multimodal_cell_names"))
pt <- as.numeric(h5read(H5AD, "obs/dpt_pseudotime"))
h5closeAll()
r2a <- function(x) {
  left <- sub("\\.[^.]*$", "", x); seq <- sub("^[^.]*\\.", "", x)
  samp <- sub("_[^_]*$", "", left); well <- sub("^.*_", "", left)
  paste0("TFScreen#", samp, "_ATAC_", well, "_", seq, "-1")
}
proj$pseudotime <- pt[match(proj$cellNames, r2a(mm))]
message("cells with pseudotime: ", sum(!is.na(proj$pseudotime)), " / ", nCells(proj))
proj <- proj[!is.na(proj$pseudotime), ]

# ---- adaptive even-trajectory segments; PT_seg01 = root (lowest pt) ----------
pt   <- proj$pseudotime
maxb <- floor(max(pt) / BIN_WIDTH)
b    <- pmin(floor(pt / BIN_WIDTH), maxb)          # 0..maxb equal-width fine bins
cnt  <- tabulate(b + 1L, nbins = maxb + 1L)        # cells per fine bin (pt order)
newid <- integer(maxb + 1L); cur <- 1L; acc <- 0L  # merge forward until >= BIN_FLOOR
for (i in seq_len(maxb + 1L)) {
  newid[i] <- cur; acc <- acc + cnt[i]
  if (acc >= BIN_FLOOR) { cur <- cur + 1L; acc <- 0L }
}
if (acc > 0 && cur > 1L) newid[newid == cur] <- cur - 1L   # fold trailing small bin
seg <- newid[b + 1L]
lvl <- sort(unique(seg))
proj$ptSeg <- factor(sprintf("PT_seg%02d", seg), levels = sprintf("PT_seg%02d", lvl))
tb <- table(proj$ptSeg)
message("even-trajectory segments: ", length(lvl),
        "  cells/segment range: ", min(tb), "-", max(tb))

# ---- reproducible peaks grouped by pseudotime segment -----------------------
proj <- addGroupCoverages(proj, groupBy = "ptSeg", minCells = 40, maxCells = 1000,
                          minReplicates = 3, maxReplicates = 5, threads = THREADS)
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(proj, groupBy = "ptSeg", reproducibility = "2",
                               maxPeaks = 300000, minCells = 40, cutOff = 0.01,
                               method = "q", threads = THREADS, pathToMacs2 = pathToMacs2)
proj <- addPeakMatrix(proj, threads = THREADS)
saveArchRProject(proj, outputDirectory = PROJ_DIR, load = FALSE)

peaks <- data.frame(getPeakSet(proj))
write.csv(peaks, file.path(PROJ_DIR, "PeakSet_pseudotime.csv"), row.names = FALSE)
message("DONE. Pseudotime peak set: ", nrow(peaks), " peaks -> ", PROJ_DIR)
