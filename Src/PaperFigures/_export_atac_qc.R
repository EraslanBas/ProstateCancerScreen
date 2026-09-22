#!/usr/bin/env Rscript
# _export_atac_qc.R — extract ATAC QC data (ArchR) -> CSV for the Python QC panel.
# -----------------------------------------------------------------------------
# Follows the repo pattern (R computes, dumps CSV; a notebook plots). Produces
# only the QC quantities that are MEANINGFUL for this dataset. The input is a
# single-base insertion-site fragments file (width-1), so ArchR's fragment-size /
# nucleosome metrics are degenerate (nDiFrags = nMultiFrags = NucleosomeRatio = 0,
# nMonoFrags = nFrags) and are deliberately NOT exported.
#
# Outputs (CSV_Files/QC/):
#   prefilter_tss_frags.csv  all barcodes (subsampled): nFrags, TSSEnrichment, Keep
#                            -> the canonical TSS-vs-unique-fragments cell-calling plot
#   cellstate_qc.csv         analysis cells: cellState + per-cell TSS/nFrags/ratios
#   tss_profile.csv          aggregate TSS enrichment profile (ArchR plotTSSEnrichment)
#   counts.csv               filtering waterfall (barcodes -> passQC -> labelled -> analysis)
#
# Usage:  Rscript Src/PaperFigures/_export_atac_qc.R
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(ArchR); library(rhdf5) })
addArchRThreads(8); addArchRGenome("hg38")
set.seed(1)

DATA_DIR  <- "/data/AbbasTFScreen"
ARROW     <- file.path(DATA_DIR, "TFScreen.arrow")
PREFILT   <- file.path(DATA_DIR, "QualityControl/TFScreen/TFScreen-Pre-Filter-Metadata.rds")
STATEPROJ <- file.path(DATA_DIR, "TFScreenATAC_cellstate")
OUT       <- "/home/eraslab1/Projects/AbbasTFScreen/CSV_Files/QC"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

counts <- list()

# ---- 1. pre-filter metadata: all barcodes -> TSS-vs-frags cell-calling plot ---
pf <- as.data.frame(readRDS(PREFILT))
counts$barcodes_total <- nrow(pf)
counts$barcodes_keep  <- sum(pf$Keep == 1, na.rm = TRUE)
message("pre-filter barcodes: ", nrow(pf), " | Keep: ", counts$barcodes_keep)
# subsample for a plottable density scatter, but keep ALL retained (Keep==1) cells
keep_idx <- which(pf$Keep == 1)
rest_idx <- which(!(seq_len(nrow(pf)) %in% keep_idx))
n_rest   <- min(length(rest_idx), 150000)
samp     <- c(keep_idx, sample(rest_idx, n_rest))
pfs <- pf[samp, c("nFrags", "TSSEnrichment", "Keep")]
pfs <- pfs[pfs$nFrags >= 1, ]                       # log10 axis
write.csv(pfs, file.path(OUT, "prefilter_tss_frags.csv"), row.names = FALSE)
message("wrote prefilter_tss_frags.csv (", nrow(pfs), " points; all Keep + ", n_rest, " sampled background)")

# ---- 2. analysis cells: per-cell QC + cell state (from the cell-state project) --
proj <- loadArchRProject(STATEPROJ, showLogo = FALSE)
cc <- as.data.frame(getCellColData(proj))
num <- intersect(c("cellState","TSSEnrichment","nFrags","ReadsInTSS","ReadsInPromoter",
                   "ReadsInBlacklist","PromoterRatio","BlacklistRatio","FRIP"), colnames(cc))
qc <- cc[, num, drop = FALSE]
qc$log10nFrags <- log10(qc$nFrags)
write.csv(qc, file.path(OUT, "cellstate_qc.csv"), row.names = TRUE)
counts$cells_analysis <- nrow(qc)
counts$cells_passQC   <- length(h5read(ARROW, "Metadata/PassQC"))
message("wrote cellstate_qc.csv (", nrow(qc), " analysis cells, ",
        length(unique(qc$cellState)), " states)")

# ---- 3. aggregate TSS enrichment profile (hallmark ATAC QC) -------------------
tss <- tryCatch(
  plotTSSEnrichment(ArchRProj = proj, groupBy = "Sample", returnDF = TRUE),
  error = function(e) { message("plotTSSEnrichment returnDF failed: ", conditionMessage(e)); NULL })
if (!is.null(tss)) {
  write.csv(as.data.frame(tss), file.path(OUT, "tss_profile.csv"), row.names = FALSE)
  message("wrote tss_profile.csv (", nrow(tss), " rows)")
}

# ---- 4. filtering waterfall counts ------------------------------------------
labelled <- sum(!is.na(cc$cellState))
cdf <- data.frame(
  stage = c("all barcodes", "pass QC (arrow)", "with cell-state label", "analysis set"),
  cells = c(counts$barcodes_total, counts$cells_passQC, labelled, counts$cells_analysis))
write.csv(cdf, file.path(OUT, "counts.csv"), row.names = FALSE)
print(cdf)
message("DONE. QC CSVs -> ", OUT)
