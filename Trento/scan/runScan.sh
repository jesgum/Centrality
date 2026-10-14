#!/bin/bash
# Trento scan: for each scan point
#   1) run trento                 -> PbPb.dat
#   2) saveTrento.cc              -> outfile.root          (hEntropy, hNpNc, hNsNp, hNsNc)
#   3) runTrentoFit.cc            -> trentofit.root/.pdf   (Trento + NBD fit to data)
#   4) runTrentoCalib.cc          -> trentofit_calib.root  (centrality calibration, Npart/Ncoll)
# Every scan point gets its own directory, e.g. $OUTDIR/p_0.5 (one-at-a-time) or
# $OUTDIR/norm_18.5_p_0.5_... (grid). `./runScan.sh --list` only prints the scan points,
# `./runScan.sh -y` skips the confirmation prompt.
#
# Needs trento and ROOT with the O2Physics multiplicity tools (multGlauberNBDFitter, multCalibrator)
# in the environment, e.g. run inside `alienv enter O2Physics/latest`.
set -euo pipefail

#### Configuration ##################################################
OUTDIR="output"
NEVENTS=1000000
NJOBS=160 # number of scan points processed in parallel

# Trento parameters (default values, the scanned parameters are overwritten per point)
PROJ="Pb"
norm=18.5 # normalization
p=0.007   # reduced-thickness
k=0.918   # fluctuation
w=0.4     # nucleon-width [fm]
x=6.82    # cross-section [fm^2]
d=1.27    # nucleon-min-dist [fm]
GRID_MAX=10
GRID_STEP=0.2

# Scan range per parameter as "min max nsteps" (nsteps points, min and max included),
# leave empty to keep the default above. All empty -> single run with the defaults.
SCAN_norm="100 200 20"
SCAN_p="-0.1 0.1 10"
SCAN_k="1.8 2.2 40"
SCAN_w=""
SCAN_x=""
SCAN_d="0 0.8 40"
# oneAtATime: vary one parameter at a time, the others at their default
# grid:       every combination of the values above (number of points multiplies quickly!)
SCAN_MODE="grid"

# runTrentoFit: data to fit
DATA_FILE="data/AR_544122.root"
DATA_HIST="centrality-study/hFT0C_BCs"

# Steps to run (set to 0 to e.g. redo only the fit/calibration on existing output)
DO_TRENTO=1
DO_SAVE=1
DO_FIT=1
DO_CALIB=1
KEEP_DAT=1 # keep the (~100 MB per 1M events) trento text output
#####################################################################

SCAN_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[[ "$DATA_FILE" = /* ]] || DATA_FILE="$SCAN_DIR/$DATA_FILE"
[[ "$OUTDIR" = /* ]] || OUTDIR="$SCAN_DIR/$OUTDIR"

SCAN_PARAMS="norm p k w x d"

# Expand "min max nsteps" into the list of values
scanValues() {
  local range="SCAN_$1"
  [[ -n "${!range}" ]] || return 0
  read -r min max nsteps extra <<<"${!range}"
  if [[ -z "$nsteps" || -n "$extra" ]] || ! [[ "$nsteps" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: SCAN_$1=\"${!range}\" should be \"min max nsteps\" with nsteps >= 1" >&2
    return 1
  fi
  awk -v lo="$min" -v hi="$max" -v n="$nsteps" 'BEGIN {
    for (i = 0; i < n; i++) printf "%g ", (n == 1) ? lo : lo + i * (hi - lo) / (n - 1)
  }'
}
for par in $SCAN_PARAMS; do
  values=$(scanValues "$par") || exit 1
  printf -v "SCAN_VALUES_$par" '%s' "$values"
done

# Build the list of scan points, each one a string like "p=0.5 d=0.4"
points=()
if [[ $SCAN_MODE == "oneAtATime" ]]; then
  for par in $SCAN_PARAMS; do
    vals="SCAN_VALUES_$par"
    for val in ${!vals}; do
      points+=("$par=$val")
    done
  done
elif [[ $SCAN_MODE == "grid" ]]; then
  points=("")
  for par in $SCAN_PARAMS; do
    vals="SCAN_VALUES_$par"
    [[ -n "${!vals}" ]] || continue
    newPoints=()
    for pt in "${points[@]}"; do
      for val in ${!vals}; do
        newPoints+=("${pt:+$pt }$par=$val")
      done
    done
    points=("${newPoints[@]}")
  done
  [[ -n "${points[0]}" ]] || points=()
else
  echo "ERROR: unknown SCAN_MODE '$SCAN_MODE' (use oneAtATime or grid)" >&2
  exit 1
fi
((${#points[@]} > 0)) || points=("default")

pointName() {
  local name=${1// /_}
  echo "${name//=/_}"
}

LIST_ONLY=0
ASSUME_YES=0
for arg in "$@"; do
  case $arg in
    --list) LIST_ONLY=1 ;;
    -y | --yes) ASSUME_YES=1 ;;
    *)
      echo "Usage: $0 [--list] [-y|--yes]" >&2
      exit 1
      ;;
  esac
done

if [[ $LIST_ONLY -eq 1 ]]; then
  for pt in "${points[@]}"; do
    echo "$(pointName "$pt")"
  done
  echo "${#points[@]} scan point(s)"
  exit 0
fi

if [[ $ASSUME_YES -eq 0 ]]; then
  while true; do
    read -r -p "This scan will spawn ${#points[@]} total jobs with $NJOBS running at a time. Are you sure you want to proceed? [yes/no] " answer || answer="no"
    case ${answer,,} in
      y | yes) break ;;
      n | no)
        echo "Aborted"
        exit 0
        ;;
      *) echo "Please answer yes/y or no/n" ;;
    esac
  done
fi

if [[ $DO_TRENTO -eq 1 ]] && ! command -v trento >/dev/null; then
  echo "ERROR: trento not found in PATH" >&2
  exit 1
fi

writeConfig() {
  cat >"$1" <<EOF
number-events = $NEVENTS
normalization = $norm
reduced-thickness = $p
fluctuation = $k
nucleon-width = $w
cross-section = $x
nucleon-min-dist = $d
ncoll = 1
grid-max = $GRID_MAX
grid-step = $GRID_STEP
EOF
}

# root.exe does not always return an error code, so check that the expected output exists
requireFile() {
  if [[ ! -s "$1" ]]; then
    echo "ERROR: $PWD/$1 missing after $2, check $PWD/$2.log" >&2
    return 1
  fi
}

runPoint() {
  local name=$1
  local dir="$OUTDIR/$name"
  mkdir -p "$dir"
  cd "$dir"
  trap 'echo "[$name] FAILED, see logs in $dir" >&2' ERR

  if [[ $DO_TRENTO -eq 1 ]]; then
    echo "[$name] running trento ($NEVENTS events)"
    writeConfig trento.conf
    trento "$PROJ" "$PROJ" -c trento.conf >PbPb.dat 2>trento.log
    requireFile PbPb.dat trento
  fi
  if [[ $DO_SAVE -eq 1 ]]; then
    echo "[$name] saveTrento"
    root.exe -q -b -l "$SCAN_DIR/saveTrento.cc(\"PbPb.dat\", \"outfile.root\")" &>saveTrento.log
    requireFile outfile.root saveTrento
    [[ $KEEP_DAT -eq 1 ]] || rm -f PbPb.dat
  fi
  if [[ $DO_FIT -eq 1 ]]; then
    echo "[$name] runTrentoFit"
    root.exe -q -b -l "$SCAN_DIR/runTrentoFit.cc(\"$DATA_FILE\", \"$DATA_HIST\", \"outfile.root\", \"trentofit.root\")" &>runTrentoFit.log
    requireFile trentofit.root runTrentoFit
  fi
  if [[ $DO_CALIB -eq 1 ]]; then
    echo "[$name] runTrentoCalib"
    root.exe -q -b -l "$SCAN_DIR/runTrentoCalib.cc(\"trentofit.root\", \"outfile.root\")" &>runTrentoCalib.log
    requireFile trentofit_calib.root runTrentoCalib
  fi
  echo "[$name] done -> $dir"
}

echo "Running ${#points[@]} scan point(s) with $NJOBS parallel job(s)"
nFailed=0
for pt in "${points[@]}"; do
  while (($(jobs -rp | wc -l) >= NJOBS)); do
    wait -n || nFailed=$((nFailed + 1))
  done
  # subshell so the scanned parameters do not leak into the other points
  (
    if [[ $pt != "default" ]]; then
      for assignment in $pt; do
        printf -v "${assignment%%=*}" '%s' "${assignment#*=}"
      done
    fi
    runPoint "$(pointName "$pt")"
  ) &
done

while true; do
  wait -n && continue
  status=$?
  ((status == 127)) && break # no jobs left
  nFailed=$((nFailed + 1))
done
echo "Scan finished, output in $OUTDIR"
if ((nFailed > 0)); then
  echo "WARNING: $nFailed scan point(s) failed" >&2
  exit 1
fi
