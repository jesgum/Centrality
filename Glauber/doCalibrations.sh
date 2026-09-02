#!/bin/bash

AR_PATH="../AnalysisResults/"
PBPB_2025="${AR_PATH}LHC25_pass1/"
PBPB_2024="${AR_PATH}LHC24_pass3/"
PBPB_2023="${AR_PATH}LHC23_pass5/"
NeNe_2025="${AR_PATH}LHC25af_pass2/"
OO_2025="${AR_PATH}LHC25ae_pass2/"

# TString lInputFileName = "results/AR_544122_glauberNBD_ancestorMode2_hFT0C_BCs.root"
# double anchorPointPercentage = 90.0
# double matchRange = 200.0
# bool doNpartNcoll = false

# AR_567905.root

# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_567905_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568242_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568146_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568067_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568409_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568117_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&

# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_567911_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_567912_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568367_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568522_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_567905_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568524_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs.root\"\,90.0,200.0,kTRUE\)&

# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568117_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb.root\"\,90.0,200.0,kTRUE\,\"basehistos_Pb\"\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568117_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_PbHN.root\"\,90.0,200.0,kTRUE\,\"basehistos_PbHN\"\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568117_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbpn.root\"\,90.0,200.0,kTRUE\,\"basehistos_Pbpn\"\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568117_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbstar.root\"\,90.0,200.0,kTRUE\,\"basehistos_Pbstar\"\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568117_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbpnrw.root\"\,90.0,200.0,kTRUE\,\"basehistos_Pbpnrw\"\)&
# root.exe -q -b runCalibration.cc\(\"${PBPB_2025}AR_568117_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pbrw.root\"\,90.0,200.0,kTRUE\,\"basehistos_Pbrw\"\)&

DO_NPART_NCOLL=kTRUE

# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&


# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_BCs.root\",${DO_NPAR/T_NCOLL}\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_Collisions.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_BCs.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&

# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_BCs.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_Collisions.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_BCs.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&

# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_BCs.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_Collisions.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_BCs.root\"\,${ANCHOR_POINT_PERCENTAGE},200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&

# root.exe -b -q ../Macros/drawSummaryPlots.cc
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_Collisions.root\"\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_Collisions.root\"\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_Collisions.root\"\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_BCs.root\"\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_BCs.root\"\)&
# root.exe -q -b runCalibrationLightIon.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_BCs.root\"\)&

#!/bin/bash

# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,95,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,90,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,80,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,60,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,50,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,95,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,90,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,80,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,60,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,50,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&


# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564387_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564400_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564414_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564430_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564445_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564356_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564359_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564373_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564374_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564387_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564400_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564414_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564430_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"${OO_2025}AR_564445_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&


root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2_extra2/AR_564374_glauberNBD_hFT0M_Collisions_vtxZ.root\",${DO_NPART_NCOLL}\)&
root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2_extra2/AR_564374_glauberNBD_hFT0M_Collisions_nopileup.root\",${DO_NPART_NCOLL}\)&
root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2_extra2/AR_564374_glauberNBD_hFT0M_Collisions_nosel.root\",${DO_NPART_NCOLL}\)&


wait
echo "done"