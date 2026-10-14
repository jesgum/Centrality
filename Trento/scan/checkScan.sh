#!/bin/bash
# Progress of a (running) runScan.sh: number of finished, running and failed scan points,
# plus the number of trento events done for the running points.
# Settings (OUTDIR, NEVENTS, DO_* steps) are read from runScan.sh, so run this with the
# same runScan.sh configuration as the scan.
#   ./checkScan.sh        summary
#   ./checkScan.sh -v     also list the running and failed points
# Use e.g. `watch -n 60 ./checkScan.sh` to keep it updated.
set -euo pipefail

SCAN_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RUNSCAN="$SCAN_DIR/runScan.sh"

VERBOSE=0
case ${1:-} in
  "") ;;
  -v | --verbose) VERBOSE=1 ;;
  *)
    echo "Usage: $0 [-v|--verbose]" >&2
    exit 1
    ;;
esac

# Read a "VAR=value" setting from the configuration block of runScan.sh
getSetting() {
  sed -n -E "s/^$1=\"?([^\" #]*)\"?.*/\1/p" "$RUNSCAN" | head -1
}
OUTDIR=$(getSetting OUTDIR)
NEVENTS=$(getSetting NEVENTS)
[[ "$OUTDIR" = /* ]] || OUTDIR="$SCAN_DIR/$OUTDIR"

# The last enabled step decides which file marks a finished point
FINAL_FILE=""
[[ $(getSetting DO_TRENTO) -eq 1 ]] && FINAL_FILE="PbPb.dat"
[[ $(getSetting DO_SAVE) -eq 1 ]] && FINAL_FILE="outfile.root"
[[ $(getSetting DO_FIT) -eq 1 ]] && FINAL_FILE="trentofit.root"
[[ $(getSetting DO_CALIB) -eq 1 ]] && FINAL_FILE="trentofit_calib.root"
if [[ -z "$FINAL_FILE" ]]; then
  echo "ERROR: no step enabled in $RUNSCAN" >&2
  exit 1
fi

if [[ ! -d "$OUTDIR" ]]; then
  echo "No output directory $OUTDIR yet, scan not started"
  exit 0
fi

# Total number of points from the scan ranges ("min max nsteps"), same logic as runScan.sh
# (`runScan.sh --list` is too slow for large grids)
SCAN_MODE=$(getSetting SCAN_MODE)
nTotal=0
for par in norm p k w x d; do
  nsteps=$(sed -n -E "s/^SCAN_$par=\"[^\"]*[[:space:]]([0-9]+)[[:space:]]*\".*/\1/p" "$RUNSCAN" | head -1)
  [[ -n "$nsteps" ]] || continue
  if [[ $SCAN_MODE == "grid" ]]; then
    nTotal=$(( (nTotal == 0 ? 1 : nTotal) * nsteps ))
  else
    nTotal=$((nTotal + nsteps))
  fi
done
((nTotal > 0)) || nTotal=1

# Points with a running trento/root.exe process: the working directory of that process
declare -A running=()
for pid in $(pgrep -u "$USER" -x 'trento|root.exe' || true); do
  cwd=$(readlink "/proc/$pid/cwd" 2>/dev/null) || continue
  [[ "$cwd" == "$OUTDIR"/* ]] && running[$cwd]=1
done

nStarted=0
nDone=0
failed=()
for dir in "$OUTDIR"/*/; do
  [[ -d "$dir" ]] || continue
  dir=${dir%/}
  nStarted=$((nStarted + 1))
  if [[ -s "$dir/$FINAL_FILE" ]]; then
    nDone=$((nDone + 1))
  elif [[ -z "${running[$dir]:-}" ]]; then
    failed+=("$dir")
  fi
done
nRunning=${#running[@]}
nFailed=${#failed[@]}

printf 'Scan progress (%s)\n' "$OUTDIR"
printf '  done:     %d / %d (%.1f%%)\n' $nDone "$nTotal" "$(awk -v a=$nDone -v b="$nTotal" 'BEGIN { print (b > 0) ? 100 * a / b : 0 }')"
printf '  running:  %d\n' $nRunning
printf '  waiting:  %d\n' $((nTotal - nStarted))
printf '  failed:   %d%s\n' $nFailed "$( ((nFailed > 0)) && echo '  (or between two steps, rerun to confirm)')"

# Seconds -> "1d 02h 03m"
formatDuration() {
  local s=$1
  if ((s >= 86400)); then
    printf '%dd %02dh %02dm' $((s / 86400)) $((s % 86400 / 3600)) $((s % 3600 / 60))
  else
    printf '%dh %02dm' $((s / 3600)) $((s % 3600 / 60))
  fi
}

# Running time from the runScan.sh process, ETA from the rate of points finished since it started
# (points finished by an earlier run in the same OUTDIR are not counted in the rate)
scanPid=$(pgrep -o -u "$USER" -f "^[^ ]*bash [^ ]*runScan\.sh( |$)" || true)
if [[ -n "$scanPid" ]]; then
  elapsed=$(ps -o etimes= -p "$scanPid" | tr -d ' ')
  now=$(date +%s)
  tStart=$((now - elapsed))
  nDoneNow=$(find "$OUTDIR" -mindepth 2 -maxdepth 2 -name "$FINAL_FILE" -newermt "@$tStart" | wc -l)
  nLeft=$((nTotal - nDone))
  printf '  started:  %s (running for %s)\n' "$(date -d "@$tStart" '+%a %d %b %H:%M')" "$(formatDuration "$elapsed")"
  if ((nDoneNow > 0)); then
    eta=$((elapsed * nLeft / nDoneNow))
    printf '  ETA:      %s (in %s, %.1f points/h)\n' "$(date -d "@$((now + eta))" '+%a %d %b %H:%M')" \
      "$(formatDuration $eta)" "$(awk -v n=$nDoneNow -v t="$elapsed" 'BEGIN { print 3600 * n / t }')"
  else
    echo "  ETA:      available once the first point has finished"
  fi
else
  echo "  runScan.sh is not running"
fi

if [[ $VERBOSE -eq 1 ]]; then
  if ((nRunning > 0)); then
    echo "Running points:"
    for dir in "${!running[@]}"; do
      if [[ ! -s "$dir/outfile.root" && -f "$dir/PbPb.dat" ]]; then
        # trento writes one line per event
        printf '  %s  trento %d / %d events\n' "${dir##*/}" "$(wc -l <"$dir/PbPb.dat")" "$NEVENTS"
      else
        printf '  %s  %s\n' "${dir##*/}" "$(ls -t "$dir"/*.log 2>/dev/null | head -1 | xargs -r basename)"
      fi
    done | sort
  fi
  if ((nFailed > 0)); then
    echo "Failed points (check the logs in these directories):"
    printf '  %s\n' "${failed[@]##*/}"
  fi
fi
