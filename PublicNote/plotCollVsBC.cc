#include <cstdio>

#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TString.h"

#include "plotHelpers.h"

void plotCollVsBC()
{
  plot::init();

  for (int year : { 2024, 2025 }) {
    TFile* f = TFile::Open(Form("Files/hCollVsBC_Zoom_Pb_%d.root", year));
    TCanvas* c = (TCanvas*)f->Get("canvCollVsBC_Zoom");
    c->SetName(Form("canvCollVsBC_Zoom_%d", year));
    c->Draw();
    plot::removeStats(c);

    printf("\n=== %d ===\n", year);
    plot::printCanvasContents(c);

    // Edit the plot here, e.g.:
    TH1* h = plot::get<TH1>(c, "hCoT0_Zoom");  // frame of the upper pad; same as plot::getFrame(c, "padUp2")
    h->GetYaxis()->SetTitle("Normalised count");
    plot::removeLatex(c, "Normalization range: [500, 50000]");


    // Move the legend (it sits in the upper pad; coordinates are NDC of that pad)
    TPad* padUp = plot::get<TPad>(c, "padUp2");
    TLegend* leg = plot::get<TLegend>(padUp, "TPave");
    printf("legend at x = [%.2f, %.2f], y = [%.2f, %.2f]\n", leg->GetX1NDC(), leg->GetX2NDC(), leg->GetY1NDC(), leg->GetY2NDC());
    leg->SetX1NDC(0.64);
    leg->SetX2NDC(0.97);
    leg->SetY1NDC(0.4);
    leg->SetY2NDC(0.6);
    leg->SetFillColorAlpha(0, 0);
    padUp->Modified();

    c->cd();
    TLatex* lat = new TLatex();
    lat->SetTextSize(0.05);
    lat->DrawLatexNDC(0.16, 0.9, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    lat->SetTextFont(42);
    lat->SetTextSize(0.04);
    
    if (year == 2024) {
      plot::removeLatex(c, "Pb#minusPb 2024");
      lat->DrawLatexNDC(0.16, 0.84, "Data-taking year: 2024");
    }
    if (year == 2025) {
      plot::removeLatex(c, "Pb#minusPb 2025");
      lat->DrawLatexNDC(0.16, 0.84, "Data-taking year: 2025");
    }

    c->SaveAs(Form("Figures/hCollVsBC_Zoom_PbPb_%d.pdf", year));
  }
}
