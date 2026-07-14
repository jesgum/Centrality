#!/bin/bash

AR_PATH="../AnalysisResults/"
PBPB_2025="${AR_PATH}LHC25_pass1/"
PBPB_2024="${AR_PATH}LHC24_pass3/"
PBPB_2023="${AR_PATH}LHC23_pass5/"
NeNe_2025="${AR_PATH}LHC25af_pass2/"
OO_2025="${AR_PATH}LHC25ae_pass2/"

# TString lInputFileName = "AnalysisResultsLHC24ar.root"
# Double_t lFitRange = /*350*/500.
# TString histogramName = "hFT0C_BCs"
# int ancestorMode = 2
# Bool_t lFreek = kFALSE
# Bool_t use_dMu_dNanc = kFALSE
# Bool_t lFreef = kFALSE
# Float_t lfvalue = 0.800

###  High statistics runs
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_567905.root\"\,500,\"hFT0C_BCs\"\,2\)
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568242.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568146.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568067.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568409.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568117.root\"\,500,\"hFT0C_BCs\"\,2\)&

###  Low statistics runs
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_567911.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_567912.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568367.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568522.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_567905.root\"\,500,\"hFT0C_BCs\"\,2\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568524.root\"\,500,\"hFT0C_BCs\"\,2\)&

###  Basehistos test
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568117.root\"\,500,\"hFT0C_BCs\"\,\"basehistos_Pb\"\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568117.root\"\,500,\"hFT0C_BCs\"\,\"basehistos_PbHN\"\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568117.root\"\,500,\"hFT0C_BCs\"\,\"basehistos_Pbpn\"\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568117.root\"\,500,\"hFT0C_BCs\"\,\"basehistos_Pbstar\"\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568117.root\"\,500,\"hFT0C_BCs\"\,\"basehistos_Pbpnrw\"\)&
# root.exe -q -b runGlauFit.cc\(\"${PBPB_2025}AR_568117.root\"\,500,\"hFT0C_BCs\"\,\"basehistos_Pbrw\"\)&

###  OO: Fitting BC vs collisions
root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564374.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564356.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564359.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564373.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\)&

# root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564374.root\"\,\"hFT0M_BCs\"\,\"hFT0M_BCs\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564356.root\"\,\"hFT0M_BCs\"\,\"hFT0M_BCs\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564359.root\"\,\"hFT0M_BCs\"\,\"hFT0M_BCs\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"${OO_2025}AR_564373.root\"\,\"hFT0M_BCs\"\,\"hFT0M_BCs\"\)&

wait
echo "Done!"