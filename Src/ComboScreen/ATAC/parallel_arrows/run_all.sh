#!/usr/bin/env bash
# run_all.sh — full parallel arrow build for ALL barcodes, ALL chromosomes.
#
# Pipeline (all detached-safe, logged under $ROOT/logs):
#   1. split the 73 GB fragments file into one bgz per chromosome  (parallel tabix)
#   2. build one ArchR arrow per chromosome                        (parallel, bounded)
#   3. merge per-chromosome arrows into one master TFScreen.arrow
#
# No validBarcodes filtering — every droplet barcode is kept; valid-cell
# filtering happens downstream. Launch detached:
#   cd Src/parallel_arrows && nohup bash run_all.sh > /data/AbbasTFScreen/run_all.log 2>&1 &
set -uo pipefail

SRC=/data/AbbasTFScreen/NEPC_sec_screen_ATAC_fragments_for_archR.sorted.5col.archr.tsv.bgz
ROOT=/data/AbbasTFScreen/parallel_build
FRAGS=$ROOT/frags
BUILD=$ROOT/build
LOGS=$ROOT/logs
MASTER=/data/AbbasTFScreen/TFScreen.arrow
HERE="$(cd "$(dirname "$0")" && pwd)"

SPLIT_JOBS=12          # concurrent tabix extractions
BUILD_JOBS=12          # concurrent arrow builds — filtered to ~291K cells each
                       # build is light (~1-2GB), so memory is no longer the limit.
BUILD_THREADS=4        # threads per build
MERGE_THREADS=24
MIN_FREE_GB=30         # memory guard (mostly moot once filtered)
# Restrict every build to the ~291K real cells (from ComboScreen_processed.h5ad),
# mapped to ATAC barcodes. Empty droplets are dropped at build time -> each build
# goes from ~38GB/hours to ~1-2GB/minutes. Leave empty ("") to keep all barcodes.
VB_FILE=/data/AbbasTFScreen/valid_barcodes_combo.txt
                       # RAM (free+reclaimable, the `available` column) is free.
                       # Each big-chromosome build peaks ~40GB.

# Block until at least MIN_FREE_GB of RAM is available, so we never start a
# build into a low-memory state and trigger the OOM killer (which would
# silently corrupt an arrow — the exact failure mode we are avoiding).
wait_for_mem() {
  local waited=0 avail
  while :; do
    avail=$(free -g | awk '/Mem:/{print $7}')   # 'available' column
    if (( avail >= MIN_FREE_GB )); then break; fi
    if (( waited % 60 == 0 )); then
      echo "  [mem-guard] available ${avail}GB < ${MIN_FREE_GB}GB — holding next build ($(date +%H:%M:%S))"
    fi
    sleep 15; waited=$((waited+15))
  done
}

# Interleaved big<->small so each concurrent group of 6 mixes large and small
# chromosomes (steadier memory + a steady drip of completions).
# Remaining-to-build first, interleaved big<->small so each pair mixes a large
# and a smaller chromosome. Already-built ones (trailing) are skipped via .done.
CHRS=( chr1 chr19 chr2 chr17 chr3 chr16 chr5 chr15 chr6 chr14 chr7 chr13 \
       chrX chr12 chr8 chr11 chr9 \
       chr4 chr18 chr20 chr21 chr22 chrY chr10 )

mkdir -p "$FRAGS" "$BUILD" "$LOGS"
echo "=== $(date) START full parallel arrow build ==="

# ---- Step 1: split by chromosome (parallel) ---------------------------------
echo "--- Step 1: split by chromosome ($SPLIT_JOBS parallel) ---"
if ! bash "$HERE/01_split_by_chr.sh" "$SRC" "$FRAGS" "$SPLIT_JOBS" > "$LOGS/01_split.log" 2>&1; then
  echo "SPLIT FAILED — see $LOGS/01_split.log"; tail -5 "$LOGS/01_split.log"; exit 1
fi
echo "split done: $(ls "$FRAGS"/*.tsv.bgz 2>/dev/null | wc -l) files"

# ---- Step 2: build one arrow per chromosome (bounded parallel) --------------
echo "--- Step 2: build arrows ($BUILD_JOBS parallel x $BUILD_THREADS threads) ---"
running=0
for chr in "${CHRS[@]}"; do
  infile="$FRAGS/TFScreen_${chr}.tsv.bgz"
  if [[ ! -s "$infile" ]]; then echo "  WARN missing $infile — skipping $chr"; continue; fi
  wdir="$BUILD/$chr"; mkdir -p "$wdir"
  # Resumability: skip chromosomes already built successfully.
  if [[ -f "$wdir/.done" && -s "$wdir/TFScreen.arrow" ]]; then
    echo "  [skip] $chr already built"; continue
  fi
  wait_for_mem                      # <-- memory guard before each launch
  ( Rscript "$HERE/02_build_chr_arrow.r" "$chr" "$infile" "$wdir" "$BUILD_THREADS" "$VB_FILE" \
      > "$LOGS/build_${chr}.log" 2>&1 \
      && touch "$wdir/.done" && echo "  [ok] $chr" \
      || echo "  [FAIL] $chr (see $LOGS/build_${chr}.log)" ) &
  running=$((running+1))
  if (( running >= BUILD_JOBS )); then wait -n; running=$((running-1)); fi
done
wait
echo "all builds finished"

# Report per-chromosome fragment counts.
echo "--- per-chromosome arrow status ---"
ok=0
for chr in "${CHRS[@]}"; do
  af="$BUILD/$chr/TFScreen.arrow"
  if [[ -s "$af" ]]; then echo "  $chr : arrow present ($(du -h "$af" | cut -f1))"; ok=$((ok+1));
  else echo "  $chr : MISSING"; fi
done
echo "$ok / ${#CHRS[@]} chromosomes built"

# ---- Step 3: merge ----------------------------------------------------------
echo "--- Step 3: merge -> $MASTER ---"
chrs_csv=$(IFS=,; echo "${CHRS[*]}")
Rscript "$HERE/03_merge_arrows.r" "$BUILD" "$MASTER" "$chrs_csv" "$MERGE_THREADS" \
   > "$LOGS/03_merge.log" 2>&1 \
   && echo "MERGE OK" || { echo "MERGE FAILED (see $LOGS/03_merge.log)"; exit 1; }

tail -12 "$LOGS/03_merge.log"
echo "=== $(date) DONE. Master arrow: $MASTER ==="
