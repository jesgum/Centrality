#!/bin/bash

# AR_567905.root


# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_A_d0.0.root\"\,\"comsolpc/A_d0.0/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_A_d0.4.root\"\,\"comsolpc/A_d0.4/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_A_d0.8.root\"\,\"comsolpc/A_d0.8/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_B_d0.0.root\"\,\"comsolpc/B_d0.0/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_B_d0.4.root\"\,\"comsolpc/B_d0.4/outfile.root\"\)&
# root.exe -q -b runTrentoCalib.cc\(\"comsolpc/trentofit_B_d0.8.root\"\,\"comsolpc/B_d0.8/outfile.root\"\)&

# root.exe -q -b runTrentoCalib.cc\(\"trentofit_Jesper2.root\"\,\"outfileJesper.root\"\)
# root.exe -q -b drawTrentoPlots.cc\(\"trentofit_Jesper2.root\"\,\"trentofit_Jesper2_calib.root\"\)
# root.exe -q -b runTrentoCalib.cc\(\"trentofit_Nature.root\"\,\"outfileNature.root\"\)
# root.exe -q -b drawTrentoPlots.cc\(\"trentofit_Nature.root\"\,\"trentofit_Nature_calib.root\"\)
root.exe -q -b runTrentoCalib.cc\(\"trentofit_New.root\"\,\"outfileNature.root\"\)
root.exe -q -b drawTrentoPlots.cc\(\"trentofit_New.root\"\,\"trentofit_New_calib.root\"\)

wait
echo "done"