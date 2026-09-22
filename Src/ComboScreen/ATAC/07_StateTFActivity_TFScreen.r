#!/usr/bin/env Rscript
# 07_StateTFActivity_TFScreen.r
# -----------------------------------------------------------------------------
# Which TF (groups) are more/less active in EACH cell state.
#   1. differential peak occupancy per cell state (one-vs-rest getMarkerFeatures)
#   2. regress each state's peak Log2FC on the grouped HOCOMOCO motif matrix
#      -> OLS coefficient per TF = TF activity in that state vs the rest
#         (coef > 0 : TF's motif-peaks open in the state -> more active; < 0 : less)
#
# Well-powered (6 states, thousands of cells each) — a sanity check before the
# pseudotime-trajectory version. Reuses the existing cell-state peak set/PeakMatrix.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({ library(ArchR); library(data.table) })

PROJ        <- "/data/AbbasTFScreen/TFScreenATAC_cellstate"
PROJECT_CSV <- "/home/eraslab1/Projects/AbbasTFScreen/CSV_Files"
CSV         <- file.path(PROJ, "DiffPeaks", "CSVFiles")
FDR_THRESHOLD <- 0.1
THREADS     <- 40
addArchRThreads(THREADS); addArchRGenome("hg38")

HOCOMOCO_ALIAS <- c(GCR="NR3C1", ANDR="AR", AP2A="TFAP2A", AP2B="TFAP2B", AP2C="TFAP2C",
                    PRGR="PGR", NF2L2="NFE2L2", THA="THRA", THB="THRB",
                    NDF1="NEUROD1", NGN1="NEUROG1", TWST1="TWIST1", ANDR="AR")

# ---- 1. per-state differential peaks (each state vs the rest) ---------------
proj <- loadArchRProject(PROJ, showLogo = FALSE)
stopifnot("cellState" %in% colnames(getCellColData(proj)))
message("States: ", paste(names(table(getCellColData(proj)$cellState)), collapse=", "))

markers <- getMarkerFeatures(
  ArchRProj  = proj,
  useMatrix  = "PeakMatrix",
  groupBy    = "cellState",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)")
)
rd <- rowData(markers)
peakIdent <- paste0(rd$seqnames, "_", rd$start, "_", rd$end)
fc <- as.matrix(assays(markers)[["Log2FC"]]); rownames(fc) <- peakIdent
pv <- as.matrix(assays(markers)[["Pval"]]);   rownames(pv) <- peakIdent
colnames(fc) <- colnames(pv) <- colnames(markers)   # state names
write.csv(fc, file.path(CSV, "StatePeaksFC.csv"))
write.csv(pv, file.path(CSV, "StatePeaksPval.csv"))
message("state diff-peak matrix: ", nrow(fc), " peaks x ", ncol(fc), " states")

# ---- 2. regression: Log2FC ~ grouped motif matrix, one OLS per state --------
motif <- as.matrix(fread(file.path(PROJECT_CSV, "MotifMatrixTFGroups_hocomoco.csv")), rownames = 1)
common <- intersect(rownames(motif), rownames(fc))
motif <- motif[common, , drop = FALSE]
fc    <- fc[common, , drop = FALSE]
motif <- motif[, colSums(motif) > 0, drop = FALSE]
n <- nrow(motif); p <- ncol(motif)
message("regression: ", n, " peaks, ", p, " TFs, ", ncol(fc), " states")

X   <- cbind(Intercept = 1, motif)
XtXinv <- solve(crossprod(X))
B   <- XtXinv %*% crossprod(X, fc)
res <- fc - X %*% B
dfree  <- n - ncol(X)
sigma2 <- colSums(res^2) / dfree
SE   <- outer(sqrt(diag(XtXinv)), sqrt(sigma2))
Tval <- B / SE
Pval <- 2 * pt(-abs(Tval), df = dfree)

keep <- rownames(B) != "Intercept"
coef <- B[keep, , drop = FALSE]; pval <- Pval[keep, , drop = FALSE]
fdr  <- apply(pval, 2, p.adjust, method = "fdr"); dimnames(fdr) <- dimnames(coef)

rn <- rownames(coef); hit <- rn %in% names(HOCOMOCO_ALIAS); rn[hit] <- HOCOMOCO_ALIAS[rn[hit]]
rownames(coef) <- rownames(pval) <- rownames(fdr) <- rn
coefSig <- coef; coefSig[fdr > FDR_THRESHOLD] <- 0

write.csv(coef,    file.path(PROJECT_CSV, "StateTFActivity_coef.csv"))
write.csv(fdr,     file.path(PROJECT_CSV, "StateTFActivity_FDR.csv"))
write.csv(coefSig, file.path(PROJECT_CSV, "StateTFActivity_coefFDRsig.csv"))

# ---- 3. report: top more/less active TFs per state --------------------------
cat("\n==== TF activity per cell state (raw-Log2FC regression) ====\n")
for (s in colnames(coef)) {
  o <- order(-coef[, s])
  up <- head(o, 8); dn <- head(rev(o), 8)
  sig <- which(fdr[, s] <= FDR_THRESHOLD)
  cat("\n--- ", s, "  (", length(sig), " FDR<=", FDR_THRESHOLD, ") ---\n", sep="")
  cat("  MORE active: ", paste0(rownames(coef)[up], "(", round(coef[up, s], 2), ")", collapse=", "), "\n")
  cat("  LESS active: ", paste0(rownames(coef)[dn], "(", round(coef[dn, s], 2), ")", collapse=", "), "\n")
}
message("\nDONE. Outputs in ", PROJECT_CSV)
