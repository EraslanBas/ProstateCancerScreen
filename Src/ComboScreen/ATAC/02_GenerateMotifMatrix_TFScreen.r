#!/usr/bin/env Rscript
# 02_GenerateMotifMatrix_TFScreen.r
# -----------------------------------------------------------------------------
# Build the peaks x TFs binary motif-hit matrix for the cell-state peak set
# (the 69,042 reproducible peaks), using HOCOMOCO v11 core mono PWMs — the same
# motif set and method as PERCISTRA's notebooks/ATAC_preprocessing/02_GenerateMotifMatrix.ipynb
# (R/motifs.R: read_hocomoco_pwms + match_motifs_binary), replicated inline here
# so this script is self-contained.
#
# This is the design matrix X for the TF-activity regression (stage 05): for each
# condition (<cellState>___<perturbation>) we regress the peak Log2FC vector on
# these motif columns; the fitted coefficient per TF = TF activity vs NTC.
#
# Runs matchMotifs on the peak GRanges directly (does NOT modify the project).
# Output: <PROJ>/DiffPeaks/CSVFiles/MotifMatrix_hocomoco.csv  (rownames = chr_start_end)
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ArchR)
  library(TFBSTools)
  library(motifmatchr)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

PROJ_DIR <- "/data/AbbasTFScreen/TFScreenATAC_cellstate"
OUT_DIR  <- file.path(PROJ_DIR, "DiffPeaks", "CSVFiles")
HOCOMOCO <- "/home/eraslab1/Projects/AbbasTFScreen/Data/HOCOMOCOv11_core_pwms_HUMAN_mono.txt"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
addArchRThreads(threads = 8)
addArchRGenome("hg38")

# ---- read HOCOMOCO PWMs into a PWMatrixList (PERCISTRA read_hocomoco_pwms) ---
read_hocomoco_pwms <- function(file_path) {
  lines <- readLines(file_path)
  header_idx <- grep("^>", lines)
  if (length(header_idx) == 0) stop("No motif headers (^>) in ", file_path)
  block_ends <- c(header_idx, length(lines) + 1)
  pwm_list <- list()
  for (i in seq_along(header_idx)) {
    start <- header_idx[i]; end <- block_ends[i + 1] - 1
    motif_name <- sub("^>", "", lines[start])
    pwm_lines  <- lines[(start + 1):end]
    pwm_matrix <- do.call(rbind, lapply(pwm_lines, function(line) {
      as.numeric(strsplit(trimws(line), "\t")[[1]])
    }))
    colnames(pwm_matrix) <- c("A", "C", "G", "T")
    pwm_obj <- tryCatch(TFBSTools::PWMatrix(
      ID = motif_name, name = motif_name,
      profileMatrix = t(as.matrix(pwm_matrix)),
      bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)),
      error = function(e) { warning(motif_name, ": ", e$message); NULL })
    if (inherits(pwm_obj, "PWMatrix")) pwm_list[[motif_name]] <- pwm_obj
  }
  do.call(TFBSTools::PWMatrixList, pwm_list)
}

pwms <- read_hocomoco_pwms(HOCOMOCO)
message("HOCOMOCO PWMs: ", length(pwms))

# ---- peak set as GRanges keyed by peakIdent ---------------------------------
proj <- loadArchRProject(PROJ_DIR, showLogo = FALSE)
gr <- getPeakSet(proj)
peakIdent <- paste0(seqnames(gr), "_", start(gr), "_", end(gr))
names(gr) <- peakIdent
message("Peaks: ", length(gr))

# ---- match motifs -> peaks x motifs binary, strip TF names ------------------
mm  <- matchMotifs(pwms, subject = gr, genome = BSgenome.Hsapiens.UCSC.hg38, out = "matches")
hit <- as.matrix(assay(mm))                     # logical peaks x motifs
storage.mode(hit) <- "integer"
# strip HOCOMOCO suffix: "ASCL1_HUMAN.H11MO.0.A" -> "ASCL1"
tf_per_motif <- vapply(colnames(hit), function(x) strsplit(x, "_")[[1]][1], character(1))

# collapse motifs sharing a stripped TF name (max) -> unique TF columns
tfs <- unique(tf_per_motif)
mat <- vapply(tfs, function(t) {
  cols <- which(tf_per_motif == t)
  if (length(cols) == 1L) hit[, cols] else as.integer(rowSums(hit[, cols, drop = FALSE]) > 0L)
}, integer(nrow(hit)))
rownames(mat) <- peakIdent
colnames(mat) <- tfs

message("Motif matrix: ", nrow(mat), " peaks x ", ncol(mat), " TFs; ",
        "mean hits/peak = ", round(mean(rowSums(mat)), 1))
# HOCOMOCO uses UniProt-style names for some TFs (NR3C1->GCR, AR->ANDR, TFAP2A->AP2A);
# report which screened TFs appear (by common symbol AND HOCOMOCO alias).
screened <- c("ASCL1","KLF14","NEUROD1","NEUROG1","NR3C1","SIM1","TET2","TWIST1",
              "VSX1","ZNF385A","ZNF547","ZNF660","ZNF776")
message("Screened TFs present (direct symbol): ",
        paste(intersect(screened, colnames(mat)), collapse = ", "))

out <- file.path(OUT_DIR, "MotifMatrix_hocomoco.csv")
write.csv(mat, out)
message("DONE. Wrote ", out)
