#!/usr/bin/env bash
# rebuild_missing_serial.sh — rebuild the arrows that the parallel run lost to
# the OOM-killer, STRICTLY ONE CHROMOSOME AT A TIME.
#
# Background: run_all.sh launched 12 concurrent builds assuming the validBarcodes
# filter kept each at ~1-2GB. In reality ArchR's tabix reader loads the ENTIRE
# chromosome (all droplet barcodes) into RAM before filtering, so each build
# peaked at 40-123GB. Twelve at once blew past the 850GB / no-swap machine and
# the kernel OOM-killed the biggest workers (always the largest chromosomes),
# which surfaced as the misleading "No fragments found! Possible error with
# validBarcodes!" error.
#
# Fix: serial builds. A single big-chromosome worker peaks ~120GB, far under
# 850GB, so memory can never be the limiting factor. Resumable via .done, so a
# re-run only builds what is still missing.
#
# Launch detached:
#   cd Src/parallel_arrows && nohup bash rebuild_missing_serial.sh \
#       > /data/AbbasTFScreen/rebuild_serial.log 2>&1 &
set -uo pipefail

ROOT=/data/AbbasTFScreen/parallel_build
FRAGS=$ROOT/frags
BUILD=$ROOT/build
LOGS=$ROOT/logs
HERE="$(cd "$(dirname "$0")" && pwd)"

BUILD_THREADS=4        # proven-good thread count; only ONE build runs at a time
VB_FILE=/data/AbbasTFScreen/valid_barcodes_combo.txt

# The 11 chromosomes the parallel run lost to OOM, largest split-file first so
# the longest jobs run earliest. Any already-built ones are skipped below.
CHRS=( chr1 chr2 chr3 chr6 chr12 chr8 chr17 chr14 chr19 chr13 chr20 )

mkdir -p "$LOGS"
echo "=== $(date) START serial rebuild of ${#CHRS[@]} missing chromosomes ==="
free -g | awk '/Mem:/{print "  available RAM at start: "$7"GB"}'

built=0; skipped=0; failed=0
for chr in "${CHRS[@]}"; do
  infile="$FRAGS/TFScreen_${chr}.tsv.bgz"
  wdir="$BUILD/$chr"

  if [[ ! -s "$infile" ]]; then
    echo "  [MISS-SPLIT] $chr — no split file $infile, skipping"; failed=$((failed+1)); continue
  fi
  if [[ -f "$wdir/.done" && -s "$wdir/TFScreen.arrow" ]]; then
    echo "  [skip] $chr already built"; skipped=$((skipped+1)); continue
  fi

  mkdir -p "$wdir"
  echo "--- $(date +%H:%M:%S) building $chr (serial, ${BUILD_THREADS} threads) ---"
  avail=$(free -g | awk '/Mem:/{print $7}')
  echo "    available RAM before build: ${avail}GB"

  # Foreground: the loop does not advance until THIS build finishes. One at a time.
  if Rscript "$HERE/02_build_chr_arrow.r" "$chr" "$infile" "$wdir" "$BUILD_THREADS" "$VB_FILE" \
        > "$LOGS/build_${chr}.log" 2>&1; then
    touch "$wdir/.done"
    sz=$(du -h "$wdir/TFScreen.arrow" 2>/dev/null | cut -f1)
    echo "  [ok] $chr -> $sz"; built=$((built+1))
  else
    echo "  [FAIL] $chr (see $LOGS/build_${chr}.log)"; tail -3 "$LOGS/build_${chr}.log"
    failed=$((failed+1))
  fi
done

echo "=== $(date) DONE. built=$built skipped=$skipped failed=$failed ==="
echo "--- full per-chromosome arrow status ---"
ok=0; ALL=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
            chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY )
for chr in "${ALL[@]}"; do
  af="$BUILD/$chr/TFScreen.arrow"
  if [[ -s "$af" ]]; then echo "  $chr : present ($(du -h "$af" | cut -f1))"; ok=$((ok+1));
  else echo "  $chr : MISSING"; fi
done
echo "$ok / ${#ALL[@]} chromosomes now built"
