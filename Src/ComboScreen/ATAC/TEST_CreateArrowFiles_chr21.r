#!/usr/bin/env Rscript
# Quick verification run: chr21 only, offset disabled, to confirm fragments survive.
source("Main.R")
source("Utilities.R")
source("Conf.R")
library(ArchR)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

# --- FIX for ArchR 1.0.2 bug ---
# On the tabix/fragments read path, .createArrow calls .tmpToArrow WITHOUT
# forwarding minFragSize, so it uses the hardcoded default of 10 and discards
# every fragment shorter than 10bp. Our insertion-site fragments are width 1,
# so all are dropped. Lower the default to 0 so they are retained.
.fixedTmpToArrow <- ArchR:::.tmpToArrow
formals(.fixedTmpToArrow)$minFragSize <- 0
assignInNamespace(".tmpToArrow", .fixedTmpToArrow, ns = "ArchR")

setwd("/data/AbbasTFScreen/test_chr21")

ArrowFiles <- createArrowFiles(
  inputFiles      = c(TFScreenChr21 = "/data/AbbasTFScreen/test_chr21/TFScreen_chr21.tsv.bgz"),
  sampleNames     = "TFScreenChr21",
  minTSS          = 0,
  minFrags        = 0,
  addTileMat      = FALSE,
  addGeneScoreMat = FALSE,
  force           = TRUE,
  threads         = 8,
  minFragSize     = 0,
  maxFragSize     = 20000000,
  offsetPlus      = 0,
  offsetMinus     = 0,
  nChunk          = 1
)

cat("\n==== TEST RESULT ====\n")
meta <- readRDS("QualityControl/TFScreenChr21/TFScreenChr21-Pre-Filter-Metadata.rds")
cat("nCells:", nrow(meta), "\n")
cat("Total nFrags:", sum(meta$nFrags), "\n")
cat("Median nFrags:", median(meta$nFrags), "\n")
cat("Max nFrags:", max(meta$nFrags), "\n")
cat("Cells with >0 frags:", sum(meta$nFrags > 0), "\n")
