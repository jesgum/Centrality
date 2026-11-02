// Same plot as drawMultipleGlaubers.cc for Trento + NBD fits: data and fits with a peripheral inset on top,
// data / fit with an inset at the bottom. The inputs are Trento calibration files (output of runTrentoCalib.cc).

#include "drawMultipleGlaubers.cc"

void drawMultipleTrento()
{
  // ---- fits to compare: {calibration file, legend label, color} ----
  const std::vector<GlauberInput> inputs = {
    { "../Trento/trentofit_Jesper2_calib.root", "T_{R}ENTo Run 1 params", kRed + 1 },
    { "../Trento/trentofit_New_calib.root", "T_{R}ENTo Run 2 params", kBlue + 1 },
    { "../Trento/madhatter/norm_166.667_p_-0.1_k_2.03158_d_0.294737/trentofit_calib.root", "Paramter scan", kGreen + 2 },
  };
  drawMultipleFits(inputs, "MultipleTrento.pdf");
}
