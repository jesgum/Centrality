#!/bin/bash

AR_PATH="../AnalysisResults/"
PBPB_2025="${AR_PATH}LHC25_pass1/"
PBPB_2024="${AR_PATH}LHC24_pass3/"
PBPB_2023="${AR_PATH}LHC23_pass5/"
NeNe_2025="${AR_PATH}LHC25af_pass2/"
OO_2025="${AR_PATH}LHC25ae_pass2_small/"

# TString lInputFileName = "AnalysisResultsLHC24ar.root"
# Double_t lFitRange = /*350*/500.
# TString histogramName = "hFT0C_BCs"
# int ancestorMode = 2
# Bool_t lFreek = kFALSE
# Bool_t use_dMu_dNanc = kFALSE
# Bool_t lFreef = kFALSE
# Float_t lfvalue = 0.800

# AR_567905.root

root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_567905.root\"\,500,\"hFT0C_BCs\"\,2\)
