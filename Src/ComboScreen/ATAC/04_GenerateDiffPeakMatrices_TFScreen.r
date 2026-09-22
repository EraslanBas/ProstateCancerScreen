#!/usr/bin/env Rscript
# 04_GenerateDiffPeakMatrices_TFScreen.r
# -----------------------------------------------------------------------------
# Assemble the per-comparison differential-peak RDS files (from
# 03_CallDiffPeaks_TFScreen.r) into peaks x condition matrices, where a
# "condition" is <cellState>___<perturbation> (each vs same-state NTC).
# Mirrors PERCISTRA build_diff_peak_matrix / EpigeneticInhibitors AllPeaksFCs.
#
#   AllPeaksFCs_stateXpert.csv   : peaks x conditions, cells = Log2FC (0 if the
#                                  peak was not returned for that condition)
#   AllPeaksPVals_stateXpert.csv : peaks x conditions, cells = Wilcoxon Pval (NA if absent)
#
# Rows keyed by peakIdent "chr_start_end"; column order = condition labels.
# The FC matrix is the regression target for stage 05 (TF-activity).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

PROJ <- "/data/AbbasTFScreen/TFScreenATAC_cellstate"
DIFF <- file.path(PROJ, "DiffPeaks")
OUT  <- file.path(DIFF, "CSVFiles")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- peak universe (the 69,042 cell-state peaks), in peak-set order ---------
ps <- fread(file.path(PROJ, "PeakSet_cellstate.csv"))
peakIdent <- paste0(ps$seqnames, "_", ps$start, "_", ps$end)
message("Peak universe: ", length(peakIdent))

# ---- read every per-condition RDS, fill Log2FC / Pval columns ---------------
rds <- list.files(DIFF, pattern = "_vs_NTC\\.rds$", recursive = TRUE, full.names = TRUE)
message("Differential RDS files: ", length(rds))
if (length(rds) == 0) stop("No *_vs_NTC.rds under ", DIFF, " — run 03_CallDiffPeaks first.")

conds <- character(length(rds))
fcL <- vector("list", length(rds))
pvL <- vector("list", length(rds))
for (i in seq_along(rds)) {
  r <- as.data.frame(readRDS(rds[i]))
  conds[i] <- if ("Condition" %in% names(r) && length(r$Condition)) r$Condition[1]
              else sub("_vs_NTC\\.rds$", "", basename(rds[i]))
  idx <- match(r$peakIdent, peakIdent)
  ok  <- !is.na(idx)
  fc <- numeric(length(peakIdent))              # default 0
  pv <- rep(NA_real_, length(peakIdent))
  fc[idx[ok]] <- r$Log2FC[ok]
  pv[idx[ok]] <- r$Pval[ok]
  fcL[[i]] <- fc; pvL[[i]] <- pv
}

fcMat <- do.call(cbind, fcL); rownames(fcMat) <- peakIdent; colnames(fcMat) <- conds
pvMat <- do.call(cbind, pvL); rownames(pvMat) <- peakIdent; colnames(pvMat) <- conds

write.csv(fcMat, file.path(OUT, "AllPeaksFCs_stateXpert.csv"))
write.csv(pvMat, file.path(OUT, "AllPeaksPVals_stateXpert.csv"))

# quick per-state summary of how many conditions each state contributes
st <- sub("___.*$", "", conds)
message("Conditions per state:")
print(table(st))
message("DONE. FC matrix ", nrow(fcMat), " peaks x ", ncol(fcMat), " conditions -> ",
        file.path(OUT, "AllPeaksFCs_stateXpert.csv"))
