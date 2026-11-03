#include <cstdio>

#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"

#include "plotHelpers.h"

// Pb-Pb <Npart> and <Ncoll> vs centrality from TGlauberMC for different Pb nuclear
// profiles (Pb, Pb*, Pbpn, PbHN, Pbrw, Pbpnrw), redrawn as is.
// Input: Files/h<var>TGlauberMC.root (from Macros/drawTGlauberMC.cc), canvas
// canv_h<var>TGlauberMC (1000x900) with sub-pads
//   padTop_h<var>TGlauberMC     hFrameTop_* (frame), one hAvg<var>_local per profile, legend, label
//   padBottom_h<var>TGlauberMC  hFrameBottom_* (frame), line at 1, hRatio_* per profile (ratio to Pb)

void plotTGlauberMC()
{
  plot::init();

  for (const auto& var : { "Npart", "Ncoll" }) {
    const TString fileName = Form("Files/h%sTGlauberMC.root", var);
    TFile* f = TFile::Open(fileName);
    if (!f || f->IsZombie()) {
      printf("cannot open %s, skipping\n", fileName.Data());
      continue;
    }
    TCanvas* c = (TCanvas*)f->Get(Form("canv_h%sTGlauberMC", var));
    c->Draw();

    printf("\n=== %s TGlauberMC ===\n", var);
    plot::printCanvasContents(c);

    // ALICE label in bold (62); TLatex objects are unnamed, so find them by type
    TPad* padTop = plot::get<TPad>(c, Form("padTop_h%sTGlauberMC", var));
    for (const auto& obj : *padTop->GetListOfPrimitives()) {
      if (obj->InheritsFrom(TLatex::Class())) {
        ((TLatex*)obj)->SetTextFont(62);
      }
    }
    padTop->Modified();

    // x-axis title (on the frame of the bottom pad)
    TPad* padBot = plot::get<TPad>(c, Form("padBottom_h%sTGlauberMC", var));
    plot::getFrame(padBot)->GetXaxis()->SetTitle("FT0C Centrality (%)");
    padBot->Modified();

    c->SaveAs(Form("Figures/h%sTGlauberMC_PbPb.pdf", var));
  }
}
