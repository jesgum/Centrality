#!/bin/bash

# AR_567905.root


# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_A_d0.0.root\"\,\"comsolpc/A_d0.0/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_A_d0.4.root\"\,\"comsolpc/A_d0.4/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_A_d0.8.root\"\,\"comsolpc/A_d0.8/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_B_d0.0.root\"\,\"comsolpc/B_d0.0/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_B_d0.4.root\"\,\"comsolpc/B_d0.4/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_B_d0.8.root\"\,\"comsolpc/B_d0.8/outfile.root\"\)&

# root.exe -q -b runTrentoCalib.cc\(\"trentofit.root\"\,\"outfileNature.root\"\)
root.exe -q -b drawTrentoPlots.cc\(\"trentofit.root\"\,\"trentofit_calib.root\"\)

wait
echo "done"