#include <cstdio>

#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"

#include "plotHelpers.h"

// Impact-parameter distributions per 10% centrality class, OO (left) and Pb-Pb (right),
// redrawn as is.
// Input: Files/hImpactParameter.root (from
// service-work/centrality-study-4/DrawBExampleJesper01.C), canvas c1 (1600x600):
//   padL  OO:    hCentrality_left_<0..90> (frame: 0-10%), hCentFull_left (0-100%), label
//   padR  Pb-Pb: hCentrality_right_<0..90> (frame: 0-10%), hCentFull_right, label, legend

void plotImpactParameter()
{
  plot::init();

  TFile* f = TFile::Open("Files/hImpactParameter.root");
  if (!f || f->IsZombie()) {
    printf("cannot open Files/hImpactParameter.root\n");
    return;
  }
  TCanvas* c = (TCanvas*)f->Get("c1");
  c->SetName("canvImpactParameter");
  c->Draw();

  printf("\n=== impact parameter ===\n");
  plot::printCanvasContents(c);

  // ALICE labels: fix the text, keep position, size and font from the file.
  // TLatex objects are unnamed; their text is stored as the title.
  for (const auto& padName : { "padL", "padR" }) {
    TPad* pad = plot::get<TPad>(c, padName);
    for (const auto& obj : *pad->GetListOfPrimitives()) {
      if (obj->InheritsFrom(TLatex::Class())) {
        TLatex* lat = (TLatex*)obj;
        TString text = lat->GetTitle();
        text.ReplaceAll("#it{s}_{#it{NN}}", "#it{s}_{NN}"); // roman NN as in the other figures
        lat->SetTitle(text);
        // lat->SetTextSize(0.05);
        // lat->SetX(0.2); lat->SetY(0.9); // NDC
      }
    }
    pad->Modified();
  }

  c->SaveAs("Figures/hImpactParameter.pdf");
}
