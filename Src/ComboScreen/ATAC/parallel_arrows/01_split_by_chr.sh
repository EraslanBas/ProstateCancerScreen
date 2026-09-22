#!/usr/bin/env bash
# 01_split_by_chr.sh
#
# Parallelize the slow part of arrow creation: the single-threaded tabix
# decompression of the 73 GB fragments file. We extract each main chromosome
# into its own bgzipped + tabix-indexed fragment file, all in parallel. Each
# tabix range query uses the .tbi index to seek directly to that chromosome's
# blocks, so the 24 jobs together read the file ~once but with 24x the
# throughput of ArchR's serial reader.
#
# Output: <OUTDIR>/TFScreen_<chr>.tsv.bgz (+ .tbi) for chr1..22, X.
# chrY/chrM are excluded (ArchR's default excludeChr for hg38 drops chrM/chrY).
#
# Usage: bash 01_split_by_chr.sh [SRC_BGZ] [OUTDIR] [NJOBS]

# NOTE: deliberately NOT using `set -e`. With bounded-parallel `wait -n`, a
# background job returning non-zero (e.g. a harmless SIGPIPE) would otherwise
# abort the whole loop and orphan running extractions.
set -uo pipefail

SRC="${1:-/data/AbbasTFScreen/NEPC_sec_screen_ATAC_fragments_for_archR.sorted.5col.archr.tsv.bgz}"
OUTDIR="${2:-/data/AbbasTFScreen/perchr_fragments}"
NJOBS="${3:-24}"

CHRS=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
       chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY )

mkdir -p "$OUTDIR"
echo "Splitting $SRC -> $OUTDIR  (${#CHRS[@]} chromosomes, $NJOBS parallel)"

extract_one() {
  local chr="$1" src="$2" outdir="$3"
  local out="$outdir/TFScreen_${chr}.tsv.bgz"
  # Skip if already done (valid bgz + index present) — lets a re-run reuse work.
  if [[ -s "$out" && -s "$out.tbi" ]]; then
    echo "  [skip] $chr already present ($(du -h "$out" | cut -f1))"
    return 0
  fi
  # tabix range-query the chromosome, re-bgzip, index. Atomic via .part.
  if tabix "$src" "$chr" | bgzip -c > "${out}.part" && mv "${out}.part" "$out" \
       && tabix -p bed -f "$out"; then
    echo "  [done] $chr -> $(du -h "$out" | cut -f1)"
  else
    echo "  [FAIL] $chr (extraction error)"; rm -f "${out}.part"
    return 1
  fi
}
export -f extract_one

# Run with a bounded number of parallel jobs (plain bash, no GNU parallel dep).
pids=()
running=0
for chr in "${CHRS[@]}"; do
  extract_one "$chr" "$SRC" "$OUTDIR" &
  pids+=($!)
  running=$((running+1))
  if (( running >= NJOBS )); then
    wait -n
    running=$((running-1))
  fi
done
wait

# Verify every requested chromosome produced a non-empty indexed file.
missing=()
for chr in "${CHRS[@]}"; do
  out="$OUTDIR/TFScreen_${chr}.tsv.bgz"
  if [[ ! -s "$out" || ! -s "$out.tbi" ]]; then missing+=("$chr"); fi
done
echo "Split summary: $(( ${#CHRS[@]} - ${#missing[@]} )) / ${#CHRS[@]} chromosomes present."
ls -la "$OUTDIR"/*.tsv.bgz 2>/dev/null
if (( ${#missing[@]} > 0 )); then
  echo "ERROR: missing splits -> ${missing[*]}"
  exit 1
fi
echo "All chromosome splits complete."
