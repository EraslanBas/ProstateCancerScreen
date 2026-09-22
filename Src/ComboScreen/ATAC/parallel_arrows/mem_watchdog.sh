#!/usr/bin/env bash
# mem_watchdog.sh — OOM safety net for an already-running build phase.
#
# The bounded build loop in run_all.sh now has a pre-launch memory guard, but a
# job launched before that guard existed has none. This watchdog protects such a
# run WITHOUT killing anything (which would corrupt an arrow): when available RAM
# falls below HARD_FLOOR it SIGSTOPs the youngest build worker to halt further
# allocation; when RAM recovers past RESUME_FLOOR it SIGCONTs paused workers.
# A paused worker keeps its RSS but stops growing, so finishing builds can drain
# memory safely. SIGSTOP/SIGCONT is safe for CPU/HDF5 work — it just freezes.
#
# Usage: nohup bash mem_watchdog.sh > /data/AbbasTFScreen/mem_watchdog.log 2>&1 &
set -uo pipefail

HARD_FLOOR=${1:-12}     # GB available -> start pausing
RESUME_FLOOR=${2:-35}   # GB available -> resume paused
TAG='02_build_chr_arrow.r'

declare -A PAUSED=()

worker_pids() {  # heavy build worker R procs (the actual ~40GB ones), newest first
  for d in /proc/[0-9]*; do
    pid=${d#/proc/}
    [[ -r "$d/cmdline" ]] || continue
    if tr '\0' ' ' < "$d/cmdline" 2>/dev/null | grep -q "$TAG"; then
      rss=$(awk '/^VmRSS:/{print $2}' "$d/status" 2>/dev/null)
      [[ -n "${rss:-}" ]] && (( rss > 5000000 )) && echo "$pid"   # >~5GB = real worker
    fi
  done | sort -rn   # highest pid first ~ youngest
}

echo "=== mem_watchdog start $(date)  HARD_FLOOR=${HARD_FLOOR}GB RESUME_FLOOR=${RESUME_FLOOR}GB ==="
while :; do
  mapfile -t workers < <(worker_pids)
  # Exit once the build phase is over (no workers and merge running or run_all gone).
  if (( ${#workers[@]} == 0 )); then
    if ! pgrep -f 'bash run_all.sh' >/dev/null 2>&1 || pgrep -f '03_merge_arrows.r' >/dev/null 2>&1; then
      echo "$(date) no build workers left — watchdog exiting"; break
    fi
    sleep 10; continue
  fi

  avail=$(free -g | awk '/Mem:/{print $7}')

  # Resume paused workers if memory recovered.
  if (( avail >= RESUME_FLOOR )) && (( ${#PAUSED[@]} > 0 )); then
    for pid in "${!PAUSED[@]}"; do
      kill -CONT "$pid" 2>/dev/null && echo "$(date) [resume] pid=$pid (avail=${avail}GB)"
      unset 'PAUSED[$pid]'
    done
  fi

  # Pause youngest running worker if memory critical.
  if (( avail < HARD_FLOOR )); then
    for pid in "${workers[@]}"; do
      if [[ -z "${PAUSED[$pid]:-}" ]]; then
        kill -STOP "$pid" 2>/dev/null && PAUSED[$pid]=1 \
          && echo "$(date) [PAUSE] pid=$pid (avail=${avail}GB < ${HARD_FLOOR}GB) — halting allocation"
        break   # pause one at a time, re-check next loop
      fi
    done
  fi

  sleep 10
done
echo "=== mem_watchdog end $(date) ==="
