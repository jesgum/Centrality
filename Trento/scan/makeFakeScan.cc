#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TString.h"
#include "TSystem.h"

#include <cmath>
#include <fstream>
#include <iostream>

// Fake output of runScan.sh to test drawScan.cc: a norm x p x w grid of scan points, each with outfile.root,
// trentofit.root, trentofit_calib.root and trento.conf. <Npart>/<Ncoll> are the reference scaled by a
// parameter dependent factor, so the combined deviation spans ~0-40%.
//   root -l -b -q makeFakeScan.cc
//   root -l 'drawScan.cc("output/fake", "", "", "Reference", 5, 20.f)'

namespace
{
// y -> y * scale(cent), keeping the x binning and the per-bin y distribution shape
TH2* scaleY(const TH2* ref, const char* name, double dev0, double devSlope)
{
  TH2* hist = static_cast<TH2*>(ref->Clone(name));
  hist->Reset();
  for (int ix = 1; ix <= ref->GetNbinsX(); ++ix) {
    const double cent = ref->GetXaxis()->GetBinCenter(ix);
    const double scale = 1 + dev0 + devSlope * cent / 100;
    for (int iy = 1; iy <= ref->GetNbinsY(); ++iy) {
      const double content = ref->GetBinContent(ix, iy);
      if (content > 0) {
        hist->Fill(cent, ref->GetYaxis()->GetBinCenter(iy) * scale, content);
      }
    }
  }
  return hist;
}
} // namespace

void makeFakeScan(TString outDir = "output/fake",
                  TString reference = "data/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root")
{
  TH1::AddDirectory(kFALSE);
  TFile* fileRef = TFile::Open(reference, "READ");
  if (!fileRef || fileRef->IsZombie()) {
    std::cerr << "Error: cannot open reference file " << reference << std::endl;
    return;
  }
  TH2* refNpart = dynamic_cast<TH2*>(fileRef->Get("h2dNpart"));
  TH2* refNcoll = dynamic_cast<TH2*>(fileRef->Get("h2dNcoll"));
  if (!refNpart || !refNcoll) {
    std::cerr << "Error: h2dNpart/h2dNcoll not found in " << reference << std::endl;
    return;
  }

  TRandom3 rnd(42);
  // multiplicity distribution: falling exponential with a smeared knee
  TF1 shape("shape", "exp(-x/[0]) * 0.5 * TMath::Erfc((x - [1]) / [2])", 0, 30000);
  const double norms[] = { 16, 18, 20 };
  const double ps[] = { -0.5, 0, 0.5 };
  const double ws[] = { 0.4, 0.6, 0.8, 1.0 };
  int nPoints = 0;
  for (double norm : norms) {
    for (double p : ps) {
      for (double w : ws) {
        const TString name = Form("norm_%g_p_%g_w_%g", norm, p, w);
        const TString dir = outDir + "/" + name;
        gSystem->mkdir(dir, kTRUE);

        // deviation from the reference grows with |p|, |w - 0.6| and |norm - 18|; Ncoll more sensitive than Npart
        const double dev0 = 0.25 * p + 0.3 * (w - 0.6) + 0.01 * (norm - 18) + rnd.Gaus(0, 0.005);
        const double devSlope = 0.1 * p - 0.05 * (w - 0.6);

        // outfile.root: entropy distribution
        TH1D hEntropy("hEntropy", "", 1000, 0, 500 * norm / 18);
        for (int i = 0; i < 100000; ++i) {
          hEntropy.Fill(hEntropy.GetXaxis()->GetXmax() * std::pow(rnd.Rndm(), 1.5 + p));
        }
        TFile fileOut(dir + "/outfile.root", "RECREATE");
        hEntropy.Write();
        fileOut.Close();

        // trentofit.root: data, NBD fit, fit parameters and range
        TH1D hData("hV0MUltraFine", "", 3000, 0, 30000);
        TH1D hFit("hGlauber", "", 3000, 0, 30000);
        shape.SetParameters(4000, 20000, 1500);
        for (int ibin = 1; ibin <= hData.GetNbinsX(); ++ibin) {
          const double x = hData.GetBinCenter(ibin);
          const double expected = 1e3 * shape.Eval(x);
          hData.SetBinContent(ibin, rnd.Poisson(expected));
          hData.SetBinError(ibin, std::sqrt(expected));
        }
        shape.SetParameters(4000 * (1 + 0.5 * dev0), 20000 * (1 + 0.3 * dev0), 1500);
        for (int ibin = 1; ibin <= hFit.GetNbinsX(); ++ibin) {
          hFit.SetBinContent(ibin, 1e3 * shape.Eval(hFit.GetBinCenter(ibin)));
        }
        TH1D hPars("hGlauberParameters", "", 7, 0.5, 7.5);
        hPars.SetBinContent(1, 30 * norm / 18 * (1 + dev0)); // mu
        hPars.SetBinContent(2, 0.8 + 0.1 * w);               // k
        hPars.SetBinContent(7, 120);                         // ndf
        hPars.SetBinContent(6, 120 * (1 + 20 * dev0 * dev0)); // chi2
        TH1D hFitRange("hGlauberFitRange", "", 2, 0.5, 2.5);
        hFitRange.SetBinContent(1, 2000);
        hFitRange.SetBinContent(2, 25000);
        TFile fileFit(dir + "/trentofit.root", "RECREATE");
        hData.Write();
        hFit.Write();
        hPars.Write();
        hFitRange.Write();
        fileFit.Close();

        // trentofit_calib.root: anchor point and Npart/Ncoll vs. centrality
        TH1D hAnchor("hAnchorPoint", "", 1, 0, 1);
        hAnchor.SetBinContent(1, 1500 + 100 * p);
        TH2* h2dNpart = scaleY(refNpart, "h2dNpart", dev0, devSlope);
        TH2* h2dNcoll = scaleY(refNcoll, "h2dNcoll", 2 * dev0, 2 * devSlope);
        TFile fileCalib(dir + "/trentofit_calib.root", "RECREATE");
        hAnchor.Write();
        h2dNpart->Write();
        h2dNcoll->Write();
        fileCalib.Close();
        delete h2dNpart;
        delete h2dNcoll;

        std::ofstream config((dir + "/trento.conf").Data());
        config << "projectile = Pb\nprojectile = Pb\nnumber-events = 1000000\n"
               << "normalization = " << norm << "\nreduced-thickness = " << p << "\nnucleon-width = " << w << "\n";
        ++nPoints;
      }
    }
  }
  fileRef->Close();
  std::cout << "Wrote " << nPoints << " fake scan points to " << outDir << std::endl;
}
