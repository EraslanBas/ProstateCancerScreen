#!/usr/bin/env Rscript
# 13_PertTFActivity_TFScreen.r
# -----------------------------------------------------------------------------
# THIRD pass: TF activity per PERTURBATION, POOLED across all cell states
# (i.e. NOT stratified by state, unlike 03/05). For each perturbation we test
# its differential peak occupancy vs the NTC cells over the whole population,
# then regress the per-peak Log2FC on the grouped HOCOMOCO motif matrix -> which
# TFs are more/less active under each perturbation, marginalizing over state.
#
# Same machinery as 03_CallDiffPeaks_TFScreen.r + 07 regression, but grouping by
# perturbation_clean (character!) with bgdGroups = NTC. Reuses the cell-state
# peak set / PeakMatrix and MotifMatrixTFGroups_hocomoco.csv.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ArchR); library(data.table); library(parallel); library(reshape2)
})

PROJ <- "/data/AbbasTFScreen/TFScreenATAC_cellstate"
CSVF <- "/home/eraslab1/Projects/AbbasTFScreen/CSV_Files"
CONTROL <- "NTC"; MIN_TEST_CELLS <- 40
FDR_THRESHOLD <- 0.1; THREADS <- 40
addArchRThreads(THREADS); addArchRGenome("hg38")

HOCOMOCO_ALIAS <- c(GCR="NR3C1", ANDR="AR", AP2A="TFAP2A", AP2B="TFAP2B", AP2C="TFAP2C",
                    PRGR="PGR", NF2L2="NFE2L2", THA="THRA", THB="THRB",
                    NDF1="NEUROD1", NGN1="NEUROG1", TWST1="TWIST1")

# ---- load project; expose perturbation as CHARACTER grouping ----------------
proj <- loadArchRProject(PROJ, showLogo = FALSE)
cc <- getCellColData(proj)
stopifnot("perturbation_clean" %in% colnames(cc))
proj$pertChr <- as.character(cc$perturbation_clean)   # character -> avoids the
                                                      # factor table() bias-match bug
cnt <- table(proj$pertChr)
stopifnot(CONTROL %in% names(cnt) && cnt[CONTROL] >= MIN_TEST_CELLS)
perts <- setdiff(names(cnt)[cnt >= MIN_TEST_CELLS], CONTROL)
message("NTC cells: ", cnt[CONTROL], " | perturbations vs NTC (>=", MIN_TEST_CELLS,
        " cells): ", length(perts))

# ---- differential: each perturbation vs pooled NTC --------------------------
diffPert <- function(p) {
  tryCatch({
    mt <- getMarkerFeatures(
      ArchRProj  = proj, useMatrix = "PeakMatrix", groupBy = "pertChr",
      testMethod = "wilcoxon", bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups  = p, bgdGroups = CONTROL, threads = 1)
    rd <- rowData(mt)
    l2 <- unlist(assays(mt)$Log2FC)
    names(l2) <- paste0(rd$seqnames, "_", rd$start, "_", rd$end)
    l2
  }, error = function(e) { message("FAIL ", p, " vs ", CONTROL, ": ", conditionMessage(e)); NULL })
}
res <- mclapply(perts, diffPert, mc.cores = 6)
names(res) <- perts
res <- res[!vapply(res, is.null, logical(1))]
stopifnot(length(res) > 0)
pid <- names(res[[1]])
fc  <- vapply(res, function(v) v[pid], numeric(length(pid)))
rownames(fc) <- pid                                    # peaks x perturbations (Log2FC vs NTC)
write.csv(fc, file.path(CSVF, "PertPeaksFC.csv"))
message("perturbation Log2FC vs NTC: ", nrow(fc), " peaks x ", ncol(fc), " perturbations",
        if (length(res) < length(perts)) paste0("  (", length(perts) - length(res), " errored)") else "")

# ---- regression: Log2FC ~ grouped motif matrix, one OLS per perturbation ----
motif <- as.matrix(fread(file.path(CSVF, "MotifMatrixTFGroups_hocomoco.csv")), rownames = 1)
common <- intersect(rownames(motif), rownames(fc))
motif <- motif[common, , drop = FALSE]; fc <- fc[common, , drop = FALSE]
motif <- motif[, colSums(motif) > 0, drop = FALSE]
n <- nrow(motif)
message("regression: ", n, " peaks, ", ncol(motif), " TFs, ", ncol(fc), " perturbations")

X <- cbind(Intercept = 1, motif); XtXinv <- solve(crossprod(X))
B <- XtXinv %*% crossprod(X, fc); resid <- fc - X %*% B
dfree <- n - ncol(X); sigma2 <- colSums(resid^2) / dfree
SE <- outer(sqrt(diag(XtXinv)), sqrt(sigma2)); Tval <- B / SE
Pval <- 2 * pt(-abs(Tval), df = dfree)
keep <- rownames(B) != "Intercept"
coef <- B[keep, , drop = FALSE]; pval <- Pval[keep, , drop = FALSE]
fdr  <- apply(pval, 2, p.adjust, method = "fdr"); dimnames(fdr) <- dimnames(coef)
rn <- rownames(coef); h <- rn %in% names(HOCOMOCO_ALIAS); rn[h] <- HOCOMOCO_ALIAS[rn[h]]
rownames(coef) <- rownames(pval) <- rownames(fdr) <- rn
coefSig <- coef; coefSig[fdr > FDR_THRESHOLD] <- 0

write.csv(coef,    file.path(CSVF, "PertTFActivity_coef.csv"))
write.csv(fdr,     file.path(CSVF, "PertTFActivity_FDR.csv"))
write.csv(coefSig, file.path(CSVF, "PertTFActivity_coefFDRsig.csv"))

# ---- tidy long table of significant calls -----------------------------------
mC <- melt(coef); colnames(mC) <- c("TF", "Perturbation", "coef")
mF <- melt(fdr)
long <- cbind(mC, FDR = mF$value)
long$direction <- ifelse(long$coef > 0, "more_active", "less_active")
sig <- long[long$FDR <= FDR_THRESHOLD, c("Perturbation", "TF", "coef", "FDR", "direction")]
sig <- sig[order(sig$Perturbation, -abs(sig$coef)), ]
write.csv(sig, file.path(CSVF, "PertTFActivity_significant_long.csv"), row.names = FALSE)
message("significant (FDR<=", FDR_THRESHOLD, "): ", nrow(sig), " across ",
        length(unique(sig$Perturbation)), " perturbations  | more=",
        sum(sig$direction == "more_active"), " less=", sum(sig$direction == "less_active"))

# ---- report -----------------------------------------------------------------
cat("\n==== TF activity per perturbation (pooled vs NTC) ====\n")
for (p in colnames(coef)) {
  sg <- which(fdr[, p] <= FDR_THRESHOLD); if (length(sg) == 0) next
  o <- order(-coef[, p])
  cat("\n--- ", p, "  (", length(sg), " sig) ---\n", sep = "")
  cat("  MORE: ", paste0(rownames(coef)[head(o, 6)], "(", round(coef[head(o, 6), p], 2), ")", collapse = ", "), "\n")
  cat("  LESS: ", paste0(rownames(coef)[head(rev(o), 6)], "(", round(coef[head(rev(o), 6), p], 2), ")", collapse = ", "), "\n")
}
message("\nDONE. Outputs in ", CSVF)
