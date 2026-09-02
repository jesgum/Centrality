#!/bin/bash

DO_NPART_NCOLL=kTRUE
AR_PATH=""

OO_2025=""
root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\,0\)&
root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\,1\)&
wait
root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\",kTRUE\)&
root.exe -q -b 'drawSummaryPlots.cc("../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root","564374")'
root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions_0.root\",kTRUE\)&
root.exe -q -b 'drawSummaryPlots.cc("../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions_0.root","564374")'
root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions_1.root\",kTRUE\)&
root.exe -q -b 'drawSummaryPlots.cc("../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions_1.root","564374")'

# TString lInputFileName
# TString lRunNumber = "000000",
# TString lXaxisTitle = "FT0C Amplitude",
# Bool_t lMode = kFALSE

# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\,1\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564356.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564359.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564373.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564387.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564400.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564414.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564430.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564445.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_bc_vtxZ\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2_systematics/AR_564374.root\"\,\"hFT0M_BCs\"\,\"hFT0M_BCs\"\,\"\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2_extra2/AR_564374.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"\"\)&
# root.exe -q -b runGlauFitOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374.root\"\,\"hFT0M_BCs\"\,\"hFT0M_BCs\"\,\"\"\)&

# root.exe -q -b runGlauFitNeNe.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_evsel\"\)&
# root.exe -q -b runGlauFitNeNe.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564472.root\"\,\"hFT0M_Collisions\"\,\"hFT0M_Collisions\"\,\"_evsel\"\)&
# root.exe -q -b runGlauFitNeNe.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468.root\"\,\"hFT0M_BCs\"\,\"hFT0M_BCs\"\,\"\"\)&

# This is what we'll compare with
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions_0.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions_1.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564356_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564359_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564373_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564387_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&


# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564400_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564414_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564430_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564445_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationOO.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&

# root.exe -q -b runCalibrationNeNe.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibrationNeNe.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564472_glauberNBD_hFT0M_Collisions.root\",${DO_NPART_NCOLL}\)&

# Anchor point comparison
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\"\,50,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\"\,60,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\"\,80,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2/AR_564374_glauberNBD_hFT0M_Collisions.root\"\,90,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&

# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,50,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,60,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,80,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,90,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,91,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,92,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&


# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,93,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,94,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,95,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,96,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,97,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,98,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_Collisions.root\"\,99,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&

# Comparison to Col with BC selections

# Comparison to BC selections
# root.exe -q -b runCalibrationLightIon.cc\(\"../AnalysisResults/LHC25ae_pass2_systematics/AR_564374_glauberNBD_hFT0M_BCs.root\",${DO_NPART_NCOLL}\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2_systematics/AR_564374_glauberNBD_hFT0M_BCs.root\"\,50,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2_systematics/AR_564374_glauberNBD_hFT0M_BCs.root\"\,60,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2_systematics/AR_564374_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2_systematics/AR_564374_glauberNBD_hFT0M_BCs.root\"\,80,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25ae_pass2_systematics/AR_564374_glauberNBD_hFT0M_BCs.root\"\,90,200.0,${DO_NPART_NCOLL}\,\"basehistos_Oho2_Oho2_omega0.3\"\)&

# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_BCs.root\"\,50,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_BCs.root\"\,60,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_BCs.root\"\,70,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_BCs.root\"\,80,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&
# root.exe -q -b runCalibration.cc\(\"../AnalysisResults/LHC25af_pass2_systematics/AR_564468_glauberNBD_hFT0M_BCs.root\"\,90,200.0,${DO_NPART_NCOLL}\,\"basehistos_Ne_Ne_omega0.3\"\)&

echo "done"
# pingme "done"