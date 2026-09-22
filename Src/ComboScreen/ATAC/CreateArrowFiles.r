source("Main.R")
source("Utilities.R")
source("Conf.R")
library(ArchR)

addArchRThreads(threads = 30)
addArchRGenome("hg38")

# --- FIX for ArchR 1.0.2 bug (tabix/fragments read path) ---
# .createArrow calls .tmpToArrow WITHOUT forwarding minFragSize on the tabix
# path, so it uses the hardcoded default of 10 and silently discards every
# fragment shorter than 10bp. This input is a single-base insertion-site file
# (all fragments are width 1), so without this patch ALL fragments are dropped
# (verified: produces an arrow with 0 fragments). Lower the default to 0 so the
# width-1 insertions are retained.
.fixedTmpToArrow <- ArchR:::.tmpToArrow
formals(.fixedTmpToArrow)$minFragSize <- 0
assignInNamespace(".tmpToArrow", .fixedTmpToArrow, ns = "ArchR")

setwd("/data/AbbasTFScreen")


out_dir <- "/data/AbbasTFScreen/ArchR_ArrowFiles"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ArrowFiles <- createArrowFiles(
  inputFiles      = c(TFScreen = "/data/AbbasTFScreen/NEPC_sec_screen_ATAC_fragments_for_archR.sorted.5col.archr.tsv.bgz"),
  sampleNames     = "TFScreen",
  minTSS          = 0,
  minFrags        = 0,
  maxFrags        = 1e9,
  addTileMat      = FALSE,
  addGeneScoreMat = FALSE,
  force           = TRUE,
  threads         = 30,
  minFragSize=0,
  maxFragSize=20000000,
  offsetPlus      = 0,
  offsetMinus     = 0,
  nChunk=2
)

# --- POST-BUILD VALIDATION ---------------------------------------------------
# The Jun 2026 build "completed successfully" but silently produced an arrow in
# which 10 of 24 chromosomes (incl. chr1 & chr2) had 0 fragments: the filtered-
# arrow step found temp chunk files for only 23 of 48 chunks and wrote empty
# chromosomes for the rest, with no error raised. createArrowFiles will NOT
# catch this, so we verify per-chromosome fragment counts ourselves and abort
# loudly if any expected chromosome is empty.
suppressPackageStartupMessages(library(rhdf5))

expectedChr <- paste0("chr", c(1:22, "X", "Y"))

for (af in ArrowFiles) {
  stopifnot(file.exists(af))
  if (h5read(af, "Class") != "Arrow") {
    stop("Validation FAILED: ", af, " is not a valid Arrow file (Class != 'Arrow').")
  }

  ls    <- h5ls(af)
  chrs  <- unique(ls$name[ls$group == "/Fragments"])
  chrs  <- chrs[!grepl("Info", chrs)]

  fragByChr <- vapply(expectedChr, function(ch) {
    if (!(ch %in% chrs)) return(0)
    sum(h5read(af, paste0("Fragments/", ch, "/RGLengths")))
  }, numeric(1))
  h5closeAll()

  emptyChr <- names(fragByChr)[fragByChr == 0]

  cat("\n==== Arrow validation:", af, "====\n")
  print(fragByChr)
  cat("Total fragments:", format(sum(fragByChr), big.mark = ","), "\n")

  if (length(emptyChr) > 0) {
    stop("Validation FAILED for ", af, ": ", length(emptyChr),
         " chromosome(s) have 0 fragments -> ", paste(emptyChr, collapse = ", "),
         "\nThe arrow is INCOMPLETE and must be rebuilt (check disk/temp-chunk writes).")
  }
}

cat("\n==== Arrow validation PASSED: all", length(expectedChr),
    "chromosomes contain fragments ====\n")