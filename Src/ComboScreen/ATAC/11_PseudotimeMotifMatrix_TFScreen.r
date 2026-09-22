#!/usr/bin/env Rscript
# 11_PseudotimeMotifMatrix_TFScreen.r
# HOCOMOCO motif matrix for the PSEUDOTIME peak set, collapsed by the SAME TF
# grouping used for the cell-state analysis (TFGroups_hocomoco.csv) so TF columns
# stay comparable across the two exercises.
# Output: CSV_Files/MotifMatrixTFGroups_hocomoco_pseudotime.csv (rownames chr_start_end)

suppressPackageStartupMessages({
  library(ArchR); library(TFBSTools); library(motifmatchr)
  library(SummarizedExperiment); library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38); library(data.table)
})
PROJ <- "/data/AbbasTFScreen/TFScreenATAC_pseudotime"
CSVF <- "/home/eraslab1/Projects/AbbasTFScreen/CSV_Files"
HOCO <- "/home/eraslab1/Projects/AbbasTFScreen/Data/HOCOMOCOv11_core_pwms_HUMAN_mono.txt"
addArchRThreads(8); addArchRGenome("hg38")

read_hocomoco_pwms <- function(fp) {
  L <- readLines(fp); hi <- grep("^>", L); be <- c(hi, length(L) + 1); out <- list()
  for (i in seq_along(hi)) {
    nm <- sub("^>", "", L[hi[i]])
    M <- do.call(rbind, lapply(L[(hi[i]+1):(be[i+1]-1)], function(x) as.numeric(strsplit(trimws(x), "\t")[[1]])))
    colnames(M) <- c("A", "C", "G", "T")   # REQUIRED: PWMatrix needs A/C/G/T rownames after t()
    o <- tryCatch(TFBSTools::PWMatrix(ID=nm, name=nm, profileMatrix=t(as.matrix(M)),
                  bg=c(A=.25,C=.25,G=.25,T=.25)), error=function(e) NULL)
    if (inherits(o, "PWMatrix")) out[[nm]] <- o
  }
  do.call(TFBSTools::PWMatrixList, out)
}
pwms <- read_hocomoco_pwms(HOCO)

proj <- loadArchRProject(PROJ, showLogo = FALSE)
gr <- getPeakSet(proj)
pid <- paste0(seqnames(gr), "_", start(gr), "_", end(gr)); names(gr) <- pid
message("pseudotime peaks: ", length(gr))

mm  <- matchMotifs(pwms, gr, genome = BSgenome.Hsapiens.UCSC.hg38, out = "matches")
hit <- as.matrix(assay(mm)); storage.mode(hit) <- "integer"
tf  <- vapply(colnames(hit), function(x) strsplit(x, "_")[[1]][1], character(1))
tfs <- unique(tf)
mat <- vapply(tfs, function(t){ c <- which(tf==t)
  if (length(c)==1) hit[,c] else as.integer(rowSums(hit[,c,drop=FALSE])>0) }, integer(nrow(hit)))
rownames(mat) <- pid; colnames(mat) <- tfs

# collapse by the cell-state TF grouping (grouped -> TFGroup_N max; ungrouped kept)
grp <- fread(file.path(CSVF, "TFGroups_hocomoco.csv"))
out <- as.data.frame(mat, check.names = FALSE)
for (g in setdiff(unique(grp$TF_group), "NA")) {
  members <- intersect(grp[grp$TF_group == g, TF_name], colnames(out))
  if (length(members) == 0) next
  out[[g]] <- if (length(members) == 1) out[[members]] else apply(out[, members, drop=FALSE], 1, max)
  out[, setdiff(members, g)] <- NULL
}
o <- file.path(CSVF, "MotifMatrixTFGroups_hocomoco_pseudotime.csv")
write.csv(out, o)
message("DONE. grouped motif matrix ", nrow(out), " peaks x ", ncol(out), " -> ", o)
