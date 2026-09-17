#!/bin/bash

INPUT_FILE="$1"
HISTOGRAM="$2"
ANCESTOR_MODE="$3"
RECO_PASS="$4"

LOCK_FILE="jobLock_Analysis_${RECO_PASS}_$(basename "${INPUT_FILE}")_${HISTOGRAM}_${ANCESTOR_MODE}"

touch "$LOCK_FILE"

root.exe -q -b runGlauFit.cc\(\"${INPUT_FILE}\",500,\"${HISTOGRAM}\",${ANCESTOR_MODE}\)

rm -f "$LOCK_FILE"