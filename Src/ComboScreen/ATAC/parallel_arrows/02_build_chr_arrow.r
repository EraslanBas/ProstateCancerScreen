#!/usr/bin/env Rscript
# 02_build_chr_arrow.r
#
# Build ONE valid ArchR arrow from a single-chromosome fragment file. This is
# the exact path proven to work in TEST_CreateArrowFiles_chr21.r (89M fragments
# survived for chr21), generalized so it can be launched once per chromosome in
# parallel processes.
#
# Each chromosome builds into its own work dir (so QC/tmp files never collide),
# and all use the SAME sampleName "TFScreen" so barcodes are keyed identically
# (TFScreen#<barcode>) for the later merge.
#
# Usage:
#   Rscript 02_build_chr_arrow.r <chr> <in_bgz> <work_dir> [threads] [valid_barcodes_file]
#
# valid_barcodes_file (optional): one ATAC barcode per line (matching the 4th
#   column of the fragment file, e.g. "250804_2lig_ATAC_10G_ACGG...-1"). When
#   given, only these barcodes are retained — i.e. the ~291K real multimodal
#   cells from SelectedCells.csv, not the millions of empty droplets. Generate
#   it from SelectedCells.csv (see parallel_arrows/README or run_all.sh).
#
# Produces <work_dir>/TFScreen.arrow containing that chromosome's fragments
# (plus empty placeholder groups for the other chromosomes, as ArchR always does).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: 02_build_chr_arrow.r <chr> <in_bgz> <work_dir> [threads] [valid_barcodes_file]")
chr      <- args[[1]]
in_bgz   <- args[[2]]
work_dir <- args[[3]]
threads  <- if (length(args) >= 4) as.integer(args[[4]]) else 4L
vb_file  <- if (length(args) >= 5 && nzchar(args[[5]])) args[[5]] else NULL

validBarcodes <- NULL
if (!is.null(vb_file)) {
  vb <- readLines(vb_file)
  vb <- vb[nzchar(vb)]
  validBarcodes <- setNames(list(vb), "TFScreen")
  message("[", chr, "] restricting to ", length(vb), " valid barcodes from ", vb_file)
}

suppressPackageStartupMessages(library(ArchR))
addArchRThreads(threads = threads)
addArchRGenome("hg38")

# --- FIX for ArchR 1.0.2 bug ---
# On the tabix read path, .createArrow calls .tmpToArrow WITHOUT forwarding
# minFragSize, so it uses the hardcoded default of 10 and discards every
# fragment shorter than 10bp. Our insertion-site fragments are width 1, so all
# would be dropped. Lower the default to 0 so they are retained.
.fixedTmpToArrow <- ArchR:::.tmpToArrow
formals(.fixedTmpToArrow)$minFragSize <- 0
assignInNamespace(".tmpToArrow", .fixedTmpToArrow, ns = "ArchR")

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

message("[", chr, "] building arrow in ", work_dir, " from ", in_bgz)

ArrowFiles <- createArrowFiles(
  inputFiles      = setNames(in_bgz, "TFScreen"),
  sampleNames     = "TFScreen",
  validBarcodes   = validBarcodes,
  minTSS          = 0,
  minFrags        = 0,
  maxFrags        = 1e9,
  addTileMat      = FALSE,
  addGeneScoreMat = FALSE,
  force           = TRUE,
  threads         = threads,
  minFragSize     = 0,
  maxFragSize     = 20000000,
  offsetPlus      = 0,
  offsetMinus     = 0,
  # nChunk splits each chromosome's tabix read into this many genomic sub-chunks.
  # nChunk=1 (one giant chunk) hangs on the largest chromosome (chr1, ~232M
  # fragments) — its single chunk crosses a size where ArchR's per-chunk
  # processing degrades pathologically (chr2 at 205M in one chunk was just under
  # it). Splitting into 16 keeps every sub-chunk well below that threshold and
  # lets threads process them in parallel. Harmless for smaller chromosomes.
  nChunk          = 16,
  # Each input file holds exactly ONE chromosome that we explicitly want, so
  # exclude nothing here (the default would drop chrY).
  excludeChr      = character(0)
)

# --- per-chromosome validation: this chr must be non-empty ---
suppressPackageStartupMessages(library(rhdf5))
af  <- file.path(work_dir, "TFScreen.arrow")
stopifnot(file.exists(af))
if (h5read(af, "Class") != "Arrow") stop("[", chr, "] not a valid Arrow (Class != 'Arrow')")
nfrag <- sum(h5read(af, paste0("Fragments/", chr, "/RGLengths")))
h5closeAll()
cat(sprintf("[%s] DONE  arrow=%s  fragments=%s\n", chr, af, format(nfrag, big.mark = ",")))
if (nfrag == 0) stop("[", chr, "] built arrow has 0 fragments — FAILED")
