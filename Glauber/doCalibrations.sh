#!/bin/bash

AR_PATH="../AnalysisResults/"
PBPB_2025="${AR_PATH}LHC25_pass1/"
PBPB_2024="${AR_PATH}LHC24_pass3/"
PBPB_2023="${AR_PATH}LHC23_pass5/"
NeNe_2025="${AR_PATH}LHC25af_pass2/"
OO_2025="${AR_PATH}LHC25ae_pass2_small/"

# TString lInputFileName = "results/AR_544122_glauberNBD_ancestorMode2_hFT0C_BCs.root"
# double anchorPointPercentage = 90.0
# double matchRange = 200.0
# bool doNpartNcoll = false

# AR_567905.root

root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_567905_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)