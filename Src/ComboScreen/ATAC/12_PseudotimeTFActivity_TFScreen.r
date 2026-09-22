#!/usr/bin/env Rscript
# 12_PseudotimeTFActivity_TFScreen.r
# TF activity ALONG the trajectory: differential peak occupancy of each pseudotime
# segment vs the ROOT segment (PT_seg01), then regress each segment's peak Log2FC
# on the grouped HOCOMOCO motif matrix -> TF activity per pseudotime bin.
# coef > 0 : TF's motif-peaks open relative to the NE root -> TF gains activity along trajectory.

suppressPackageStartupMessages({ library(ArchR); library(data.table); library(parallel) })
PROJ <- "/data/AbbasTFScreen/TFScreenATAC_pseudotime"
CSVF <- "/home/eraslab1/Projects/AbbasTFScreen/CSV_Files"
FDR_THRESHOLD <- 0.1; THREADS <- 40
addArchRThreads(THREADS); addArchRGenome("hg38")

HOCOMOCO_ALIAS <- c(GCR="NR3C1", ANDR="AR", AP2A="TFAP2A", AP2B="TFAP2B", AP2C="TFAP2C",
                    PRGR="PGR", NF2L2="NFE2L2", THA="THRA", THB="THRB",
                    NDF1="NEUROD1", NGN1="NEUROG1", TWST1="TWIST1")

proj <- loadArchRProject(PROJ, showLogo = FALSE)
lvl  <- levels(getCellColData(proj)$ptSeg)
root <- lvl[1]
segs <- setdiff(lvl, root)
message("segments: ", length(lvl), "  root=", root)

# ArchR getMarkerFeatures' bias-matcher calls table(groups[idY]); if the grouping
# column is a FACTOR, table() returns a bin for every level (all 16 segments),
# extending obsbgd to length 16 while estbgd stays length 1 -> cor() throws
# "incompatible dimensions". 03_CallDiffPeaks_TFScreen.r works only because its
# stateXpert grouping is a plain character. So expose ptSeg as CHARACTER here.
proj$ptSegChr <- as.character(getCellColData(proj)$ptSeg)

# ---- differential: each non-root segment vs the root segment ----------------
# SAME as 03_CallDiffPeaks_TFScreen.r (perturbation vs same-state NTC): one
# getMarkerFeatures Wilcoxon test per segment, useGroups=<segment> vs
# bgdGroups=<root>, ArchR doing its bias-matching + scaleTo depth normalization.
# Per-peak Log2FC(segment / root) is the regression signal. Isolated per segment
# (threads=1), run a few in parallel (mc.cores kept low to avoid arrow-connection
# read errors seen at high parallelism).
diffSeg <- function(seg) {
  tryCatch({
    mt <- getMarkerFeatures(
      ArchRProj  = proj, useMatrix = "PeakMatrix", groupBy = "ptSegChr",
      testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups  = seg, bgdGroups = root, threads = 1)
    rd <- rowData(mt)
    l2 <- unlist(assays(mt)$Log2FC)
    names(l2) <- paste0(rd$seqnames, "_", rd$start, "_", rd$end)
    l2
  }, error = function(e) { message("FAIL ", seg, " vs ", root, ": ", conditionMessage(e)); NULL })
}
res <- mclapply(segs, diffSeg, mc.cores = 6)
names(res) <- segs
res <- res[!vapply(res, is.null, logical(1))]             # drop any that errored
stopifnot(length(res) > 0)
pid <- names(res[[1]])
fc  <- vapply(res, function(v) v[pid], numeric(length(pid)))
rownames(fc) <- pid                                       # peaks x (surviving) non-root segments
lfc <- cbind(0, fc); colnames(lfc)[1] <- root             # prepend root = 0 baseline
lfc <- lfc[, intersect(lvl, colnames(lfc)), drop = FALSE] # order along trajectory
write.csv(lfc, file.path(CSVF, "PseudotimePeaksFC.csv"))  # full trajectory incl root (=0)
message("segment getMarkerFeatures Log2FC vs root: ", nrow(fc), " peaks x ", ncol(fc), " segments",
        if (length(res) < length(segs)) paste0("  (", length(segs) - length(res), " segment(s) errored)") else "")

# ---- regression: Log2FC ~ grouped motif matrix, one OLS per segment ---------
motif <- as.matrix(fread(file.path(CSVF, "MotifMatrixTFGroups_hocomoco_pseudotime.csv")), rownames = 1)
common <- intersect(rownames(motif), rownames(fc))
motif <- motif[common, , drop = FALSE]; fc <- fc[common, , drop = FALSE]
motif <- motif[, colSums(motif) > 0, drop = FALSE]
n <- nrow(motif)
X <- cbind(Intercept = 1, motif); XtXinv <- solve(crossprod(X))
B <- XtXinv %*% crossprod(X, fc); res <- fc - X %*% B
dfree <- n - ncol(X); sigma2 <- colSums(res^2) / dfree
SE <- outer(sqrt(diag(XtXinv)), sqrt(sigma2)); Tval <- B / SE
Pval <- 2 * pt(-abs(Tval), df = dfree)
keep <- rownames(B) != "Intercept"
coef <- B[keep, , drop = FALSE]; pval <- Pval[keep, , drop = FALSE]
fdr  <- apply(pval, 2, p.adjust, method = "fdr"); dimnames(fdr) <- dimnames(coef)
rn <- rownames(coef); h <- rn %in% names(HOCOMOCO_ALIAS); rn[h] <- HOCOMOCO_ALIAS[rn[h]]
rownames(coef) <- rownames(pval) <- rownames(fdr) <- rn

# prepend the root segment as the 0 baseline -> full trajectory of length(lvl)
z <- matrix(0, nrow(coef), 1, dimnames = list(rownames(coef), root))
coefF <- cbind(z, coef); fdrF <- cbind(matrix(1, nrow(coef), 1, dimnames=list(rownames(coef), root)), fdr)

write.csv(coefF, file.path(CSVF, "PseudotimeTFActivity_coef.csv"))
write.csv(fdrF,  file.path(CSVF, "PseudotimeTFActivity_FDR.csv"))
coefSig <- coefF; coefSig[fdrF > FDR_THRESHOLD] <- 0
write.csv(coefSig, file.path(CSVF, "PseudotimeTFActivity_coefFDRsig.csv"))

# ---- report: TFs rising / falling along the trajectory ----------------------
# trend = slope of coef across ordered segments (Spearman with segment index)
idx <- seq_len(ncol(coefF))
trend <- apply(coefF, 1, function(y) suppressWarnings(cor(idx, y, method = "spearman")))
nSig <- rowSums(fdrF <= FDR_THRESHOLD)
tab <- data.frame(TF = rownames(coefF), trend = round(trend, 2),
                  maxAbsCoef = round(apply(abs(coefF), 1, max), 3), nSigSeg = nSig)
tab <- tab[order(-tab$trend), ]
write.csv(tab, file.path(CSVF, "PseudotimeTFActivity_trend.csv"), row.names = FALSE)
cat("\n==== TFs RISING along trajectory (NE -> Differentiated) ====\n")
print(head(tab[tab$nSigSeg > 0, ], 15), row.names = FALSE)
cat("\n==== TFs FALLING along trajectory ====\n")
print(head(tab[tab$nSigSeg > 0, ][order(tab$trend[tab$nSigSeg>0]), ], 15), row.names = FALSE)
message("\nDONE. Outputs in ", CSVF)
