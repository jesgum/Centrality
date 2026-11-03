#include <cstdio>
#include <utility>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"

#include "plotHelpers.h"

void plothGlauber()
{
  plot::init();

  // { file tag, canvas name }
  for (auto [tag, canvName] : { std::pair{ "Col", "canvCol" }, std::pair{ "Bcs", "canvBcs" } }) {
    TFile* f = TFile::Open(Form("Files/hGlauber%s_OO.root", tag));
    TCanvas* c = (TCanvas*)f->Get(canvName);
    c->Draw();
    plot::removeStats(c);

    printf("\n=== %s ===\n", tag);
    plot::printCanvasContents(c);

    // Edit the plot here, e.g.:
    // The canvas is split into pads; each pad has its own frame (first histogram drawn)
    TPad* padMain = plot::find<TPad>(c, Form("%s_1", canvName));  // top: data + Glauber fits
    TPad* padRatio = plot::find<TPad>(c, Form("%s_2", canvName)); // bottom: data / fit
    TPad* padMainInset = plot::find<TPad>(c, "padTop");           // inset in the top pad
    TPad* padRatioInset = plot::find<TPad>(c, "padBot");          // inset in the bottom pad

    TH1* frameMain = plot::getFrame(padMain);
    TH1* frameRatio = plot::getFrame(padRatio);
    TH1* frameMainInset = plot::getFrame(padMainInset);
    TH1* frameRatioInset = plot::getFrame(padRatioInset);

    frameMain->SetTitle("");
    frameRatio->SetTitle("");
    frameMainInset->SetTitle("");
    frameRatioInset->SetTitle("");
    // frameMain->GetYaxis()->SetTitle("Counts");
    // frameRatio->GetYaxis()->SetRangeUser(0.6, 1.5);
    // if (TString(tag) == "Col") { ... }  // collisions-only tweaks

    c->cd();
    TLatex* lat = new TLatex();
    lat->SetTextSize(0.04);
    lat->DrawLatexNDC(0.1, 0.9, "ALICE OO #sqrt{#it{s}_{NN}} = 5.36 TeV");
    lat->SetTextFont(42);
    lat->SetTextSize(0.03);
    if (TString(tag) == "Col") {
      lat->DrawLatexNDC(0.1, 0.86, "Collision based");
      plot::removeLatex(c, "Collisions");
    } else {
      lat->DrawLatexNDC(0.1, 0.86, "Bunch Crossing based");
      plot::removeLatex(c, "Bunch crossings");
      frameMain->SetMaximum(9e+6);

    }

    c->SaveAs(Form("Figures/hGlauber%s_OO.pdf", tag));
  }
}
