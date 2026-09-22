#!/usr/bin/env Rscript
# 05_TFActivityRegression_TFScreen.r
# -----------------------------------------------------------------------------
# TF-ACTIVITY per (cell state x perturbation) vs NTC, via motif regression —
# the same idea as EpigeneticInhibitors 08_PredictTFGroupsOfConditions /
# PERCISTRA fit_ols_union: for each condition, regress the per-peak Log2FC
# vector (perturbation vs same-state NTC) on the peaks x TF binary motif matrix.
# The fitted OLS coefficient for TF t = how much peaks carrying TF t's motif move
# under the perturbation, controlling for all other motifs:
#     coef > 0  => TF t's target peaks OPEN   => TF t MORE active than in NTC
#     coef < 0  => TF t's target peaks CLOSE  => TF t LESS active than in NTC
#
# Inputs (stage 02 + 04):
#   MotifMatrix_hocomoco.csv       peaks x TF (0/1)
#   AllPeaksFCs_stateXpert.csv     peaks x condition (Log2FC)
# Outputs (peaks-> conditions collapsed to TF x condition):
#   TFActivity_coef_stateXpert.csv       TF x condition OLS coefficient
#   TFActivity_pval_stateXpert.csv       "         "     t-test p-value
#   TFActivity_FDR_stateXpert.csv        "         "     BH-FDR (per condition)
#   TFActivity_coefFDRsig_stateXpert.csv coef, zeroed where FDR > threshold
#   TFActivity_significant_long.csv      tidy: State, Perturbation, TF, coef, pval, FDR, direction
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(reshape2)
})

PROJ <- "/data/AbbasTFScreen/TFScreenATAC_cellstate"
CSV  <- file.path(PROJ, "DiffPeaks", "CSVFiles")                    # big matrices (in /data)
PROJECT_CSV <- "/home/eraslab1/Projects/AbbasTFScreen/CSV_Files"    # grouped motif matrix + small outputs (in repo)
dir.create(PROJECT_CSV, showWarnings = FALSE, recursive = TRUE)
FDR_THRESHOLD <- 0.1   # a TF is "active" in a condition if BH-FDR <= this

# HOCOMOCO uses UniProt-style names for some TFs; map the notable ones back to
# gene symbols for readability (extend as needed).
HOCOMOCO_ALIAS <- c(GCR = "NR3C1", ANDR = "AR", AP2A = "TFAP2A", AP2B = "TFAP2B",
                    AP2C = "TFAP2C", PRGR = "PGR", ESR1 = "ESR1", NF2L2 = "NFE2L2",
                    THA = "THRA", THB = "THRB", RARA = "RARA")

# ---- load motif matrix (X) and FC matrix (Y) --------------------------------
# Grouped motif matrix from 06_ClusterTFs_TFScreen.ipynb (TF groups collapse
# collinear motifs; ungrouped TFs kept as-is). Switch back to
# "MotifMatrix_hocomoco.csv" to regress on individual (ungrouped) TFs.
motif <- as.matrix(fread(file.path(PROJECT_CSV, "MotifMatrixTFGroups_hocomoco.csv")), rownames = 1)
fc    <- as.matrix(fread(file.path(CSV, "AllPeaksFCs_stateXpert.csv")),  rownames = 1)
pv    <- as.matrix(fread(file.path(CSV, "AllPeaksPVals_stateXpert.csv")), rownames = 1)

# Keep only SIGNIFICANT peak changes as the regression signal (as in
# EpigeneticInhibitors 02_GenerateDifferentialPeakMatrices): Log2FC where the
# per-peak Wilcoxon Pval < SIG_PVAL, 0 otherwise. Regressing on raw Log2FC for
# all peaks drowns the motif signal in non-significant noise (p-value rate ~ null).
SIG_PVAL <- 0.05
stopifnot(all(dim(fc) == dim(pv)))
fc[is.na(pv) | pv >= SIG_PVAL] <- 0
cat("kept significant peak-changes (Pval<", SIG_PVAL, "): ",
    round(100 * mean(fc != 0), 2), "% of matrix nonzero\n", sep = "")

common <- intersect(rownames(motif), rownames(fc))
message("Peaks: motif=", nrow(motif), " fc=", nrow(fc), " common=", length(common))
motif <- motif[common, , drop = FALSE]
fc    <- fc[common, , drop = FALSE]
motif <- motif[, colSums(motif) > 0, drop = FALSE]   # drop TFs with no motif hits
n <- nrow(motif); p <- ncol(motif); nc <- ncol(fc)
message("Regression: ", n, " peaks, ", p, " TFs, ", nc, " conditions")

# ---- per-condition OLS in closed form (shared design matrix) ----------------
X   <- cbind(Intercept = 1, motif)                    # n x (p+1)
XtX <- crossprod(X)
XtXinv <- tryCatch(solve(XtX), error = function(e) {
  message("XtX not invertible (", e$message, ") — adding small ridge 1e-6.")
  solve(XtX + diag(1e-6, ncol(XtX)))
})
B   <- XtXinv %*% crossprod(X, fc)                    # (p+1) x nc coefficients
res <- fc - X %*% B                                   # residuals n x nc
dfree  <- n - ncol(X)
sigma2 <- colSums(res^2) / dfree                      # per condition
SE   <- outer(sqrt(diag(XtXinv)), sqrt(sigma2))       # (p+1) x nc
Tval <- B / SE
Pval <- 2 * pt(-abs(Tval), df = dfree)

# drop intercept row
keep <- rownames(B) != "Intercept"
coef <- B[keep, , drop = FALSE]
pval <- Pval[keep, , drop = FALSE]
fdr  <- apply(pval, 2, p.adjust, method = "fdr")      # BH per condition (column)
dimnames(fdr) <- dimnames(coef)

# rename HOCOMOCO aliases -> gene symbols
rn <- rownames(coef); hit <- rn %in% names(HOCOMOCO_ALIAS)
rn[hit] <- HOCOMOCO_ALIAS[rn[hit]]
rownames(coef) <- rownames(pval) <- rownames(fdr) <- rn

coefSig <- coef; coefSig[fdr > FDR_THRESHOLD] <- 0

# ---- write matrices ---------------------------------------------------------
write.csv(coef,    file.path(PROJECT_CSV, "TFActivity_coef_stateXpert.csv"))
write.csv(pval,    file.path(PROJECT_CSV, "TFActivity_pval_stateXpert.csv"))
write.csv(fdr,     file.path(PROJECT_CSV, "TFActivity_FDR_stateXpert.csv"))
write.csv(coefSig, file.path(PROJECT_CSV, "TFActivity_coefFDRsig_stateXpert.csv"))

# ---- tidy long table of significant TF activities ---------------------------
mC <- melt(coef); colnames(mC) <- c("TF", "Condition", "coef")
mP <- melt(pval); mF <- melt(fdr)
long <- cbind(mC, pval = mP$value, FDR = mF$value)
long$State        <- sub("___.*$", "", long$Condition)
long$Perturbation <- sub("^.*___", "", long$Condition)
long$direction    <- ifelse(long$coef > 0, "more_active", "less_active")
sig <- long[long$FDR <= FDR_THRESHOLD, c("State","Perturbation","TF","coef","pval","FDR","direction")]
sig <- sig[order(sig$State, sig$Perturbation, -abs(sig$coef)), ]
write.csv(sig, file.path(PROJECT_CSV, "TFActivity_significant_long.csv"), row.names = FALSE)

message("Significant (FDR<=", FDR_THRESHOLD, ") TF-activity calls: ", nrow(sig),
        " across ", length(unique(paste(sig$State, sig$Perturbation))), " conditions")
message("  more_active: ", sum(sig$direction == "more_active"),
        " | less_active: ", sum(sig$direction == "less_active"))
message("DONE. Outputs in ", PROJECT_CSV)
