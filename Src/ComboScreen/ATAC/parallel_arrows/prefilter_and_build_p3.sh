#!/usr/bin/env bash
# prefilter_and_build_p3.sh — rebuild the missing arrows, 3 CHROMOSOMES AT A TIME.
#
# Same pipeline as prefilter_and_build.sh (pre-filter each per-chr fragment file
# to the 291,301 valid cells, then build the arrow from the small file), but runs
# up to BUILD_JOBS chromosomes concurrently. Safe on memory: each filtered build
# peaks ~40GB, so 3 at once ~120GB, far under the 850GB machine.
#
# Resumable: skips chromosomes already built (.done + arrow) and reuses any
# filtered file already present (e.g. chr1's filter, done earlier).
#
# Launch DETACHED so it survives this shell / session ending:
#   cd Src/parallel_arrows
#   setsid bash prefilter_and_build_p3.sh </dev/null \
#       > /data/AbbasTFScreen/prefilter_build_p3.log 2>&1 &
#   disown
set -uo pipefail

ROOT=/data/AbbasTFScreen/parallel_build
FRAGS=$ROOT/frags
FILT=$ROOT/frags_filtered
BUILD=$ROOT/build
LOGS=$ROOT/logs
HERE="$(cd "$(dirname "$0")" && pwd)"

VB_FILE=/data/AbbasTFScreen/valid_barcodes_combo.txt
BUILD_JOBS=3          # chromosomes in flight at once (user-requested cap)
BUILD_THREADS=8       # ArchR threads per build. Only chr1 remains (nothing else
                      # contends), so give it 8 threads to parallelize its 16
                      # read sub-chunks (see nChunk in 02_build_chr_arrow.r).
FILT_THREADS=8        # bgzip (de)compression threads per filter

# Largest-first, so the 3 heaviest (chr1/chr2/chr3) run together in the first wave.
CHRS=( chr1 chr2 chr3 chr6 chr12 chr8 chr17 chr14 chr19 chr13 chr20 )

mkdir -p "$FILT" "$LOGS"
echo "=== $(date) START 3-way parallel pre-filter+build of ${#CHRS[@]} chromosomes ==="
free -g | awk '/Mem:/{print "  available RAM at start: "$7"GB"}'

running=0
for chr in "${CHRS[@]}"; do
  wdir="$BUILD/$chr"
  if [[ -f "$wdir/.done" && -s "$wdir/TFScreen.arrow" ]]; then
    echo "  [skip] $chr already built"; continue
  fi

  # One chromosome = one detached subshell: filter (if needed) then build.
  (
    infile="$FRAGS/TFScreen_${chr}.tsv.bgz"
    filt="$FILT/TFScreen_${chr}.filt.tsv.bgz"

    if [[ ! -s "$infile" ]]; then echo "  [MISS-SPLIT] $chr"; exit 1; fi

    # ---- filter to valid barcodes (resumable, atomic) ----
    if [[ -s "$filt" && -s "$filt.tbi" ]]; then
      echo "  [filt-skip] $chr ($(du -h "$filt"|cut -f1))"
    else
      echo "  [$(date +%H:%M:%S)] filtering $chr"
      if bgzip -cd -@ "$FILT_THREADS" "$infile" \
           | awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($4 in keep)' "$VB_FILE" - \
           | bgzip -c -@ "$FILT_THREADS" > "$filt.part" \
         && mv "$filt.part" "$filt" && tabix -p bed -f "$filt"; then
        echo "  [filt-ok] $chr -> $(du -h "$filt"|cut -f1) (from $(du -h "$infile"|cut -f1))"
      else
        echo "  [FAIL-FILT] $chr"; rm -f "$filt.part"; exit 1
      fi
    fi

    # ---- build arrow from the small filtered file ----
    mkdir -p "$wdir"
    echo "  [$(date +%H:%M:%S)] building $chr"
    if Rscript "$HERE/02_build_chr_arrow.r" "$chr" "$filt" "$wdir" "$BUILD_THREADS" "$VB_FILE" \
          > "$LOGS/build_${chr}.log" 2>&1; then
      touch "$wdir/.done"
      echo "  [ok] $chr -> $(du -h "$wdir/TFScreen.arrow"|cut -f1)"
    else
      echo "  [FAIL] $chr (see $LOGS/build_${chr}.log)"; tail -4 "$LOGS/build_${chr}.log" | sed 's/^/      /'
      exit 1
    fi
  ) &

  running=$((running+1))
  if (( running >= BUILD_JOBS )); then wait -n; running=$((running-1)); fi
done
wait

echo "=== $(date) DONE (3-way parallel run) ==="
echo "--- full per-chromosome arrow status ---"
ok=0; ALL=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
            chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY )
for chr in "${ALL[@]}"; do
  af="$BUILD/$chr/TFScreen.arrow"
  if [[ -s "$af" ]]; then echo "  $chr : present ($(du -h "$af"|cut -f1))"; ok=$((ok+1));
  else echo "  $chr : MISSING"; fi
done
echo "$ok / ${#ALL[@]} chromosomes now built"
