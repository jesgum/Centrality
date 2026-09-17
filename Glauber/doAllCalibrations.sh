#!/bin/bash

echo "Will now run all relevant calibrations"

ANCHOR=90
GLOBAL_NORM=kFALSE
RECO_PASSES="LHC23_pass5 LHC24_pass3 LHC25_pass1"
JOBS=40

process_file() {
  FILE="$1"
  ANCHOR="$2"
  GLOBAL_NORM="$3"

  echo "Processing file ${FILE}"
  BASE=$(basename "$FILE" .root)
  RUN=$(echo "$BASE" | sed 's/AR_\([0-9]*\)_.*/\1/')

  if echo "$BASE" | grep -q '_BCs$'; then
    HISTTYPE=$(echo "$BASE" | sed 's/.*_fixedMu_h\(.*\)_BCs/\1/')
    EST="${HISTTYPE} Amplitude"
  else
    HISTTYPE=$(echo "$BASE" | sed 's/.*_fixedMu_hN\(.*\)/\1/')
    EST="$HISTTYPE"
  fi

  # root.exe -q -b "runCalibrationStandard.cc(\"${FILE}\",${ANCHOR},${GLOBAL_NORM})"
  root.exe -q -b "drawSummaryPlotsPbPb.cc(\"${FILE}\",\"${RUN}\",\"${EST}\",${GLOBAL_NORM})"
}

export -f process_file

for RECO_PASS in $RECO_PASSES; do
  FILES=$(ls ../AnalysisResults/${RECO_PASS}/*glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root)
  printf '%s\n' $FILES | parallel -j "$JOBS" process_file {} "$ANCHOR" "$GLOBAL_NORM"
done