#!/bin/bash
# Run all plotting macros in Macros/ (each also writes a .root next to its figures).
# Not included on purpose: uploadToCCDB.cc (uploads to CCDB), performSurgery.cc
# (modifies centrality-run-*.root in place), printEventSelectionParams.cc (prints only).
#
# Usage: ./doAllMacros.sh            run everything
#        ./doAllMacros.sh drawNpart  run only macros whose name contains "drawNpart"
# Logs go to logs/<macro>.log

cd "$(dirname "$0")" || exit 1
mkdir -p logs

macros=(
  # systematics
  'doSystematicsOO.cc'
  'doSystematicsNeNe.cc'
  'doSystematicsPbPb.cc'
  # light ion
  'drawLightIonCalibration.cc'
  'drawLightIonAllCalibration.cc'
  'drawLightIonCollisionSelections.cc'
  'drawLightIonSameSelections.cc'
  'drawLightIonMonteCarlo.cc'
  'drawLightIonGlauberFit.cc'
  'drawLightIonGlauberAnchorPoints.cc(false)'
  'drawLightIonGlauberAnchorPoints.cc(true)'
  'drawLightIonNpartNcollVsAnchorPoint.cc'
  'drawLightIonNpartNcollVsRun.cc'
  'drawNpartNcollLightIon.cc'
  'drawNpartNcollTH2.cc'
  # Pb-Pb
  'drawPbPbNpartNcollVsAnchorPoint.cc'
  'drawCollisionVsBC.cc'
  'drawNpartNcoll.cc'
  'drawNpartNcollVsRun.cc'
  'drawNpartNcollTrento.cc'
  'drawGlauberVsTrento.cc'
  'drawMultipleGlaubers.cc'
  'drawMultipleTrento.cc'
  'drawTGlauberMC.cc'
  'drawNchNancVsCentrality.cc'
  # other
  'calculateNGloToFT0C.cc'
  'drawCalibrationMC.cc'
  'drawRejectNoFoundBC.cc'
  'drawSelectionVariation.cc'
  'drawSummaryPlots.cc'
)

failed=()
for m in "${macros[@]}"; do
  [[ -n "$1" && "$m" != *"$1"* ]] && continue
  log="logs/$(echo "$m" | tr -c 'A-Za-z0-9_.\n' '_').log"
  echo ">>> $m"
  if root -l -b -q "$m" > "$log" 2>&1 && ! grep -q "^Error\|error:" "$log"; then
    echo "    ok"
  else
    echo "    FAILED (see $log)"
    failed+=("$m")
  fi
done

echo
if ((${#failed[@]})); then
  echo "${#failed[@]} macro(s) failed:"
  printf '  %s\n' "${failed[@]}"
  exit 1
fi
echo "All macros ran successfully."
