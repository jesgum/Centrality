#include <cstdio>
#include <map>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TString.h"

#include "plotHelpers.h"

void plotAnchorPointStability()
{
  plot::init();
  for (int year : { 2023, 2024, 2025, 2026 }) {
    TFile* f = TFile::Open(Form("Files/anchor_point_vs_run_%d.root", year));
    TCanvas* c = (TCanvas*)f->Get("c_anchor");
    c->SetName(Form("c_anchor_%d", year));
    c->Draw();
    plot::removeStats(c);

    printf("\n=== %d ===\n", year);
    plot::printCanvasContents(c);

    // Edit the plot here, e.g.:
    TH1F* h = (TH1F*)c->GetPrimitive("hFrame");
    // h->SetMinimum(210);
    h->GetYaxis()->SetRangeUser(140, 360);
    plot::setBoxHeights(c, 140, 360);
    // if (year == 2023) { ... }  // year-specific tweaks

    c->cd();
    TLatex* lat = new TLatex(); 
    lat->SetTextSize(0.06);
    lat->DrawLatexNDC(0.1, 0.72, "ALICE Pb#minusPb #sqrt{#it{s}_{NN}} = 5.36 TeV");
    lat->SetTextSize(0.05);
    lat->SetTextFont(42);
    lat->DrawLatexNDC(0.1, 0.65, Form("Data-taking year: %d", year));
    c->SaveAs(Form("Figures/anchor_point_vs_run_%d.pdf", year));
  }
}
