#!/bin/bash

# AR_567905.root

root.exe runTrentoCalib.cc -b -q
root.exe drawTrentoPlots.cc -b -q

wait
echo "done"