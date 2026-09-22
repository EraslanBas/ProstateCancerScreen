#!/usr/bin/env Rscript
# 03_merge_arrows.r
#
# Assemble per-chromosome arrows (built in parallel by 02_build_chr_arrow.r)
# into ONE valid ArchR arrow for the single logical sample "TFScreen".
#
# ArchR keys every cell as <sample>#<barcode> and assumes one sample == one
# arrow, so per-chromosome arrows cannot be used as a project directly (the same
# cell would live in 24 arrows). We therefore copy each chromosome's HDF5
# fragment group into one master arrow and recompute the per-cell /Metadata
# across all chromosomes.
#
#  - Fragments: copied verbatim (Ranges / RGLengths / RGValues) per chromosome.
#  - nFrags / nMonoFrags / nDiFrags / nMultiFrags: recomputed exactly from the
#    run-length data + fragment widths (no lossy approximation).
#  - TSSEnrichment / ReadsInTSS: recomputed genome-wide with ArchR's own
#    ArchR:::.fastTSSEnrichment on the assembled arrow (identical code path to
#    createArrowFiles), so the covariate matches a native build.
#
# Usage:
#   Rscript 03_merge_arrows.r <perchr_root> <out_arrow> [chrs_csv] [threads]
#     <perchr_root> : dir containing <chr>/TFScreen.arrow subdirs
#     <out_arrow>   : path of the master arrow to write
#     <chrs_csv>    : optional comma list (default chr1..22,chrX)
#     <threads>     : threads for TSS counting (default 16)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: 03_merge_arrows.r <perchr_root> <out_arrow> [chrs_csv] [threads]")
perchr_root <- args[[1]]
out_arrow   <- args[[2]]
chrs <- if (length(args) >= 3 && nzchar(args[[3]])) {
  strsplit(args[[3]], ",")[[1]]
} else {
  c(paste0("chr", 1:22), "chrX", "chrY")
}
threads <- if (length(args) >= 4) as.integer(args[[4]]) else 16L

suppressPackageStartupMessages({
  library(ArchR); library(rhdf5)
})
addArchRThreads(threads = threads)
addArchRGenome("hg38")

SAMPLE <- "TFScreen"
arrow_path <- function(chr) file.path(perchr_root, chr, "TFScreen.arrow")

# Only keep chromosomes whose per-chr arrow exists and is populated.
chrs <- Filter(function(chr) {
  af <- arrow_path(chr)
  if (!file.exists(af)) { message("skip ", chr, " (no arrow)"); return(FALSE) }
  n <- sum(h5read(af, paste0("Fragments/", chr, "/RGLengths"))); h5closeAll()
  if (n == 0) { message("skip ", chr, " (0 fragments)"); return(FALSE) }
  TRUE
}, chrs)
if (length(chrs) == 0) stop("No populated per-chromosome arrows found under ", perchr_root)
message("Merging chromosomes: ", paste(chrs, collapse = ", "))

# ---- 1. Fresh master arrow skeleton -----------------------------------------
if (file.exists(out_arrow)) file.remove(out_arrow)
h5closeAll()
h5createFile(out_arrow)
h5write("Arrow", out_arrow, "Class")
h5write(as.character(packageVersion("ArchR")), out_arrow, "ArchRVersion")  # native arrows carry this
h5createGroup(out_arrow, "Metadata")
h5createGroup(out_arrow, "Fragments")
h5write(SAMPLE, out_arrow, "Metadata/Sample")

# ---- 2. Copy each chromosome's fragments + accumulate per-cell counts --------
# nFrags / nucleosome bins are accumulated into a global per-barcode vector.
bc_levels <- character(0)      # global barcode universe (no sample prefix)
nFrags_g  <- numeric(0)
nMono_g   <- numeric(0)
nDi_g     <- numeric(0)
nMulti_g  <- numeric(0)
NUC <- 147

ensure_levels <- function(bc) {
  new <- bc[!(bc %in% bc_levels)]
  if (length(new)) {
    bc_levels <<- c(bc_levels, new)
    z <- numeric(length(new))
    nFrags_g <<- c(nFrags_g, z); nMono_g <<- c(nMono_g, z)
    nDi_g <<- c(nDi_g, z); nMulti_g <<- c(nMulti_g, z)
  }
}

for (chr in chrs) {
  af <- arrow_path(chr)
  g  <- paste0("Fragments/", chr)
  ranges    <- h5read(af, paste0(g, "/Ranges"))      # N x 2 : start, width
  rgLengths <- as.numeric(h5read(af, paste0(g, "/RGLengths")))
  rgValues  <- as.character(h5read(af, paste0(g, "/RGValues")))
  h5closeAll()

  # Copy fragment group verbatim into master.
  h5createGroup(out_arrow, g)
  # Match ArchR's NATIVE on-disk layout exactly, or .fastTSSEnrichment ->
  # .fastTSSCounts fails coercing run-lengths ("no method for coercing array to
  # LLint"): RGLengths must be INTEGER N x 1 and RGValues STRING N x 1 (NOT a
  # 1-D float / 1-D string). Ranges is already INTEGER N x 2.
  h5write(ranges,                                      out_arrow, paste0(g, "/Ranges"))
  h5write(matrix(as.integer(rgLengths),   ncol = 1L),  out_arrow, paste0(g, "/RGLengths"))
  h5write(matrix(as.character(rgValues),  ncol = 1L),  out_arrow, paste0(g, "/RGValues"))

  # Per-barcode fragment count on this chr (one run per barcode, ArchR-sorted).
  ensure_levels(rgValues)
  idx <- match(rgValues, bc_levels)
  nFrags_g[idx] <- nFrags_g[idx] + rgLengths

  # Nucleosome bin per fragment from width: 1=mono(<147), 2=di, 3=multi.
  widths <- ranges[, 2]
  wbin <- trunc(widths / NUC) + 1; wbin[wbin > 3] <- 3
  # Expand run-length barcode -> per-fragment barcode index, then tabulate bins.
  frag_bc_idx <- rep(idx, times = rgLengths)
  for (b in 1:3) {
    sel <- which(wbin == b)
    if (length(sel)) {
      tb <- tabulate(frag_bc_idx[sel], nbins = length(bc_levels))
      if (b == 1) nMono_g  <- nMono_g  + tb
      if (b == 2) nDi_g    <- nDi_g    + tb
      if (b == 3) nMulti_g <- nMulti_g + tb
    }
  }
  message(sprintf("  copied %s : %s frags, %s cells (cum cells=%s)",
                  chr, format(sum(rgLengths), big.mark=","),
                  format(length(rgValues), big.mark=","),
                  format(length(bc_levels), big.mark=",")))
}

cellNames <- paste0(SAMPLE, "#", bc_levels)
nCells <- length(cellNames)
message(sprintf("Total cells=%s  total frags=%s",
                format(nCells, big.mark=","), format(sum(nFrags_g), big.mark=",")))

# ---- 3. Genome-wide TSS enrichment via ArchR's own routine ------------------
message("Computing TSS enrichment (ArchR:::.fastTSSEnrichment) ...")
TSS <- getArchRGenome(geneAnnotation = TRUE)$TSS
tss <- ArchR:::.fastTSSEnrichment(
  TSS = TSS, ArrowFile = out_arrow, cellNames = cellNames,
  threads = threads, prefix = "(merge)", logFile = NULL
)
tssScores <- as.numeric(tss$tssScores)
tssReads  <- as.numeric(tss$tssReads)
tssScores[is.na(tssScores)] <- 0
tssReads[is.na(tssReads)]   <- 0

# ---- 4. Write /Metadata (parallel per-cell vectors) -------------------------
nan <- rep(NaN, nCells); zero <- rep(0, nCells)
meta <- list(
  # BARE barcodes (no "Sample#"): ArchR prepends the sample prefix at load time,
  # so storing them prefixed here yields double-prefixed cellNames (TFScreen#TFScreen#...).
  CellNames        = bc_levels,
  nFrags           = nFrags_g,
  nMonoFrags       = nMono_g,
  nDiFrags         = nDi_g,
  nMultiFrags      = nMulti_g,
  TSSEnrichment    = tssScores,
  ReadsInTSS       = tssReads,
  ReadsInPromoter  = zero,
  ReadsInBlacklist = zero,
  PromoterRatio    = nan,
  BlacklistRatio   = nan,
  NucleosomeRatio  = (nDi_g + nMulti_g) / pmax(nMono_g, 1),
  PassQC           = rep(1, nCells),     # minTSS=0, minFrags=0 -> all pass
  Date             = as.character(Sys.Date()),
  Completed        = "TRUE"
)
for (nm in names(meta)) h5write(meta[[nm]], out_arrow, paste0("Metadata/", nm))
h5closeAll()

# ---- 5. Validate -------------------------------------------------------------
stopifnot(h5read(out_arrow, "Class") == "Arrow")
present <- h5ls(out_arrow); h5closeAll()
fchr <- present$name[present$group == "/Fragments"]
cat("\n==== MERGE COMPLETE ====\n")
cat("Arrow:", out_arrow, "\n")
cat("Chromosomes:", paste(sort(fchr), collapse=", "), "\n")
cat("Cells:", format(nCells, big.mark=","), " Total fragments:",
    format(sum(nFrags_g), big.mark=","), "\n")
cat("Median nFrags:", median(nFrags_g),
    " Median TSS:", median(tssScores), "\n")
