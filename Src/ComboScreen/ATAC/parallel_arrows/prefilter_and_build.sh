#!/usr/bin/env bash
# prefilter_and_build.sh — rebuild the 11 missing arrows the fast, robust way.
#
# Why: ArchR's tabix reader loads EVERY droplet's fragments into RAM before it
# applies validBarcodes, so the big chromosomes crawl for hours and peak at tens
# of GB. Here we pre-shrink each per-chromosome fragment file to just the 291,301
# valid cells (a single streaming awk pass on column 4), which cuts each file
# ~10x. ArchR then reads a small file -> each arrow builds in minutes with tiny
# memory, so neither speed nor OOM is a concern.
#
# Strictly ONE chromosome at a time. Resumable: skips chromosomes already built
# (.done + arrow) and reuses any filtered file already present.
#
# Launch DETACHED so it survives this shell / session ending:
#   cd Src/parallel_arrows
#   setsid bash prefilter_and_build.sh </dev/null \
#       > /data/AbbasTFScreen/prefilter_build.log 2>&1 &
#   disown
set -uo pipefail

ROOT=/data/AbbasTFScreen/parallel_build
FRAGS=$ROOT/frags               # full per-chr fragment files (all droplets)
FILT=$ROOT/frags_filtered       # valid-cell-only per-chr fragment files
BUILD=$ROOT/build
LOGS=$ROOT/logs
HERE="$(cd "$(dirname "$0")" && pwd)"

VB_FILE=/data/AbbasTFScreen/valid_barcodes_combo.txt   # 291,301 ATAC barcodes
BUILD_THREADS=4                 # per-build ArchR threads (only one build at a time)
FILT_THREADS=8                  # bgzip (de)compression threads for the filter step

# The 11 chromosomes lost to the earlier OOM run, largest-first.
CHRS=( chr1 chr2 chr3 chr6 chr12 chr8 chr17 chr14 chr19 chr13 chr20 )

mkdir -p "$FILT" "$LOGS"
echo "=== $(date) START pre-filter + serial build of ${#CHRS[@]} chromosomes ==="
free -g | awk '/Mem:/{print "  available RAM at start: "$7"GB"}'

built=0; skipped=0; failed=0
for chr in "${CHRS[@]}"; do
  infile="$FRAGS/TFScreen_${chr}.tsv.bgz"
  filt="$FILT/TFScreen_${chr}.filt.tsv.bgz"
  wdir="$BUILD/$chr"

  if [[ ! -s "$infile" ]]; then
    echo "  [MISS-SPLIT] $chr — no split file $infile"; failed=$((failed+1)); continue
  fi
  if [[ -f "$wdir/.done" && -s "$wdir/TFScreen.arrow" ]]; then
    echo "  [skip] $chr already built"; skipped=$((skipped+1)); continue
  fi

  # ---- Step A: pre-filter to valid barcodes (resumable, atomic via .part) ----
  if [[ -s "$filt" && -s "$filt.tbi" ]]; then
    echo "  [filt-skip] $chr already filtered ($(du -h "$filt"|cut -f1))"
  else
    echo "--- $(date +%H:%M:%S) filtering $chr (col4 in $VB_FILE) ---"
    # awk: load the barcode whitelist first (NR==FNR), then stream fragments and
    # keep only rows whose 4th column (the ATAC barcode) is in the whitelist.
    if bgzip -cd -@ "$FILT_THREADS" "$infile" \
         | awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($4 in keep)' "$VB_FILE" - \
         | bgzip -c -@ "$FILT_THREADS" > "$filt.part" \
       && mv "$filt.part" "$filt" \
       && tabix -p bed -f "$filt"; then
      echo "  [filt-ok] $chr -> $(du -h "$filt"|cut -f1) (from $(du -h "$infile"|cut -f1))"
    else
      echo "  [FAIL-FILT] $chr"; rm -f "$filt.part"; failed=$((failed+1)); continue
    fi
  fi

  # ---- Step B: build arrow from the small filtered file (serial) -------------
  mkdir -p "$wdir"
  echo "--- $(date +%H:%M:%S) building $chr from filtered file (${BUILD_THREADS} threads) ---"
  avail=$(free -g | awk '/Mem:/{print $7}')
  echo "    available RAM before build: ${avail}GB"
  # validBarcodes still passed as a harmless double-check; the file is already
  # filtered, so the read is small and fast regardless.
  if Rscript "$HERE/02_build_chr_arrow.r" "$chr" "$filt" "$wdir" "$BUILD_THREADS" "$VB_FILE" \
        > "$LOGS/build_${chr}.log" 2>&1; then
    touch "$wdir/.done"
    echo "  [ok] $chr -> $(du -h "$wdir/TFScreen.arrow"|cut -f1)"; built=$((built+1))
  else
    echo "  [FAIL] $chr (see $LOGS/build_${chr}.log)"; tail -4 "$LOGS/build_${chr}.log"
    failed=$((failed+1))
  fi
done

echo "=== $(date) DONE. built=$built skipped=$skipped failed=$failed ==="
echo "--- full per-chromosome arrow status ---"
ok=0; ALL=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
            chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY )
for chr in "${ALL[@]}"; do
  af="$BUILD/$chr/TFScreen.arrow"
  if [[ -s "$af" ]]; then echo "  $chr : present ($(du -h "$af"|cut -f1))"; ok=$((ok+1));
  else echo "  $chr : MISSING"; fi
done
echo "$ok / ${#ALL[@]} chromosomes now built"
