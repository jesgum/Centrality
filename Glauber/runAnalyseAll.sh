#!/bin/bash

COUNTER=0
JOBS=32

RECO_PASSES="LHC23_pass5 LHC24_pass3 LHC25_pass1"

mkdir -p ../logs
mkdir -p ../results

for RECO_PASS in $RECO_PASSES; do
  FILES="$(ls ../AnalysisResults/${RECO_PASS}/*.root)"
  for FILE in $FILES; do
    while [ "$(ls jobLock_Analysis_* 2>/dev/null | wc -l)" -ge "$JOBS" ]; do
      echo "Sleeping. At counter: ${COUNTER}..."
      sleep 5
    done

    BASENAME=$(basename "$FILE")

    echo "[Starting processing: ${RECO_PASS}/${BASENAME}]"
    screen -d -m -L -Logfile "../logs/log_analysis_${RECO_PASS}_${BASENAME}_FT0C_2.txt" ./analyse.sh "$FILE" hFT0C_BCs 2 "$RECO_PASS"
    sleep 0.1
    let COUNTER=COUNTER+1

  done
done