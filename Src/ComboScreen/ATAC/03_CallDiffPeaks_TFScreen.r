#!/usr/bin/env Rscript
# 03_CallDiffPeaks_TFScreen.r
# -----------------------------------------------------------------------------
# Adapted from 03_CallDiffPeaks.r (Osaka Omics) for the AbbasTFScreen project.
#
# Design (as requested): peaks are defined by CELL STATE (the dominant source of
# variation) — run 01_CallPeaks_TFScreen_CellState.r first so the project has a
# cell-state ReproduciblePeakSet + PeakMatrix. Here we then, WITHIN each cell
# state, test the differential peak occupancy of every perturbation against the
# NTC (non-targeting control) cells of the SAME state — the direct analogue of
# the original's "each compound vs NOSTIM".
#
# Comparison groups use a combined key  <cellState>___<perturbation_clean>, so a
# test like  Neuroendocrine___ASCL1  vs  Neuroendocrine___NTC  is automatically
# restricted to Neuroendocrine cells for both foreground and background.
#
# Output: one RDS per (state, perturbation) under
#   <PROJ>/DiffPeaks/<state>/<perturbation>_vs_NTC.rds
# plus a combined AllDiffPeaks_stateXpert.rds.
#
# Usage:  Rscript Src/03_CallDiffPeaks_TFScreen.r
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ArchR)
  library(parallel)
})

# ---- paths / params ---------------------------------------------------------
DATA_DIR  <- "/data/AbbasTFScreen"
PROJ_DIR  <- file.path(DATA_DIR, "TFScreenATAC_cellstate")  # from 01_CallPeaks_TFScreen_CellState.r
OUT_DIR   <- file.path(PROJ_DIR, "DiffPeaks")

STATE_COL <- "cellState"           # cell-state column attached during peak calling
PERT_COL  <- "perturbation_clean"  # perturbation column attached as a covariate
CONTROL   <- "NTC"                 # background / non-targeting control group
MIN_TEST_CELLS <- 40               # min cells in BOTH fg and bg to run a test
NCORES    <- 40                    # comparisons run in parallel (each test single-threaded)

addArchRThreads(threads = 8)       # for load/setup; per-test threads set to 1 below
addArchRGenome("hg38")
setwd(DATA_DIR)

# ---- load the cell-state project (must contain PeakMatrix + perturbation_clean)
myProj <- loadArchRProject(path = PROJ_DIR, force = FALSE, showLogo = FALSE)

cc <- getCellColData(myProj)
if (!(STATE_COL %in% colnames(cc)) || !(PERT_COL %in% colnames(cc))) {
  stop("Project is missing '", STATE_COL, "' and/or '", PERT_COL,
       "'. Run 01_CallPeaks_TFScreen_CellState.r first (it attaches both).")
}
if (!("PeakMatrix" %in% getAvailableMatrices(myProj))) {
  stop("PeakMatrix not found in project. Run 01_CallPeaks_TFScreen_CellState.r ",
       "(addReproduciblePeakSet + addPeakMatrix) first.")
}

# ---- combined grouping key  <state>___<perturbation> ------------------------
myProj$stateXpert <- paste0(as.character(cc[[STATE_COL]]), "___",
                            as.character(cc[[PERT_COL]]))
cc      <- getCellColData(myProj)
counts  <- table(as.character(cc$stateXpert))
states  <- sort(unique(as.character(cc[[STATE_COL]])))
perts   <- sort(unique(as.character(cc[[PERT_COL]])))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- one differential test: perturbation vs same-state NTC ------------------
generateDiffPeaks <- function(state, pert) {
  fg <- paste0(state, "___", pert)
  bg <- paste0(state, "___", CONTROL)
  if (is.na(counts[fg]) || counts[fg] < MIN_TEST_CELLS) return(invisible(FALSE))
  if (is.na(counts[bg]) || counts[bg] < MIN_TEST_CELLS) return(invisible(FALSE))

  tryCatch({
    markerTest <- getMarkerFeatures(
      ArchRProj  = myProj,
      useMatrix  = "PeakMatrix",
      groupBy    = "stateXpert",
      testMethod = "wilcoxon",
      bias       = c("TSSEnrichment", "log10(nFrags)"),
      useGroups  = fg,
      bgdGroups  = bg,
      threads    = 1                       # parallelism is across comparisons (mclapply)
    )
    res <- rowData(markerTest)
    res$Log2FC <- unlist(assays(markerTest)$Log2FC)
    res$Pval   <- unlist(assays(markerTest)$Pval)
    res$FDR    <- unlist(assays(markerTest)$FDR)
    res <- res[!is.na(res$Pval) & !is.na(res$FDR), ]
    # peakIdent "chr_start_end" — the universal join key across the downstream
    # matrix / motif / regression stages (matches PERCISTRA/EpigeneticInhibitors).
    res$peakIdent    <- paste0(res$seqnames, "_", res$start, "_", res$end)
    res$State        <- state
    res$Perturbation <- pert
    res$Condition    <- fg          # <state>___<perturbation>, the matrix column label
    res$Background   <- bg
    res$nFg <- as.integer(counts[fg])
    res$nBg <- as.integer(counts[bg])

    sdir <- file.path(OUT_DIR, gsub("[^A-Za-z0-9._+-]", "_", state))
    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, file.path(sdir, paste0(gsub("[^A-Za-z0-9._+-]", "_", pert),
                                        "_vs_", CONTROL, ".rds")))
    TRUE
  }, error = function(e) {
    message("FAIL ", fg, " vs ", bg, ": ", conditionMessage(e)); FALSE
  })
}

# ---- build the list of viable (state, perturbation) comparisons -------------
combos <- list()
for (s in states) for (p in perts) {
  if (p == CONTROL) next
  fg <- paste0(s, "___", p); bg <- paste0(s, "___", CONTROL)
  if (!is.na(counts[fg]) && counts[fg] >= MIN_TEST_CELLS &&
      !is.na(counts[bg]) && counts[bg] >= MIN_TEST_CELLS) {
    combos[[length(combos) + 1]] <- c(s, p)
  }
}
message("Viable comparisons (perturbation vs same-state NTC, >= ", MIN_TEST_CELLS,
        " cells each): ", length(combos))

# ---- run them in parallel ---------------------------------------------------
ok <- mclapply(combos, function(sp) generateDiffPeaks(sp[[1]], sp[[2]]),
               mc.cores = NCORES)
message("Completed: ", sum(unlist(ok) == TRUE), " / ", length(combos))

# ---- combine into one table -------------------------------------------------
allf <- list.files(OUT_DIR, pattern = "_vs_NTC\\.rds$", recursive = TRUE, full.names = TRUE)
if (length(allf)) {
  allDiff <- do.call(rbind, lapply(allf, readRDS))
  saveRDS(allDiff, file.path(OUT_DIR, "AllDiffPeaks_stateXpert.rds"))
  message("Combined: ", nrow(allDiff), " rows across ", length(allf),
          " comparisons -> ", file.path(OUT_DIR, "AllDiffPeaks_stateXpert.rds"))
}
message("DONE. Per-comparison RDS under ", OUT_DIR)
