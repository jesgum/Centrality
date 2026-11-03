#include <cstdio>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPad.h"
#include "TString.h"

#include "plotHelpers.h"

// O-O <Npart> and <Ncoll> vs FT0C centrality for different oxygen nuclear-structure
// models (noTrento) and Glauber Oho2 vs TRENTo (trento).
// Input: Files/<var>_<tag>.root, one canvas (c0 or c1) with two sub-pads
//   <canvas>_1  <Npart>/<Ncoll> per model + legend (frame is the first model)
//   <canvas>_2  ratio to Oho2, omega = 0.3

// The single canvas stored in the file, whatever its name
TCanvas* getCanvas(TFile* f)
{
  for (const auto& key : *f->GetListOfKeys()) {
    TObject* obj = ((TKey*)key)->ReadObj();
    if (obj->InheritsFrom(TCanvas::Class())) {
      return (TCanvas*)obj;
    }
  }
  return nullptr;
}

// Move a legend loaded from a file, in NDC of its pad
void setLegendPosition(TLegend* leg, double x1, double y1, double x2, double y2)
{
  leg->SetX1NDC(x1);
  leg->SetY1NDC(y1);
  leg->SetX2NDC(x2);
  leg->SetY2NDC(y2);
}

void plotNpartNcollNuclearStructure()
{
  plot::init();

  // Files/<var>_<tag>.root
  for (const auto& tag : { "noTrento", "trento" }) {
    for (const auto& var : { "npart", "ncoll" }) {
      const TString fileName = Form("Files/%s_%s.root", var, tag);
      TFile* f = TFile::Open(fileName);
      if (!f || f->IsZombie()) {
        printf("cannot open %s, skipping\n", fileName.Data());
        continue;
      }
      TCanvas* c = getCanvas(f);
      if (!c) {
        printf("no canvas in %s, skipping\n", fileName.Data());
        continue;
      }
      c->SetName(Form("canv_%s_%s", var, tag));
      c->Draw();
      plot::removeStats(c);

      // Canvas size in pixels (original 796x722 or 796x772)
      c->SetCanvasSize(800, 800);
      
      printf("\n=== %s %s ===\n", var, tag);
      plot::printCanvasContents(c);

      const bool isNpart = TString(var) == "npart";
      const bool isTrento = TString(tag) == "trento";

      // The two sub-pads are the first two primitives of the canvas
      TPad* padTop = (TPad*)c->GetListOfPrimitives()->At(0);
      TPad* padBot = (TPad*)c->GetListOfPrimitives()->At(1);

      // Margins belong to the sub-pads, not the canvas (nothing is drawn in the canvas
      // itself). Keep left/right equal in both pads so the x-axes stay aligned.
      // Originals: top pad 0.14, 0.02, 0.001, 0.02; bottom pad 0.14, 0.02, 0.195, 0.001
      const double leftMargin = 0.14;
      const double rightMargin = 0.04;
      plot::setMargins(padTop, leftMargin, rightMargin, 0.001, 0.04); // left, right, bottom, top
      plot::setMargins(padBot, leftMargin, rightMargin, 0.195, 0.001);

      // Upper pad: one curve per model, the first one is the frame
      TH1* hFrame = plot::getFrame(padTop);
      hFrame->GetYaxis()->SetTitle(isNpart ? "#LT#it{N}_{part}#GT" : "#LT#it{N}_{coll}#GT");
      // hFrame->GetYaxis()->SetRangeUser(0, 30);
      padTop->Modified();

      TLegend* leg = plot::get<TLegend>(padTop, "TPave");
      for (const auto& obj : *leg->GetListOfPrimitives()) {
        TLegendEntry* entry = (TLegendEntry*)obj;
        if (TString(entry->GetLabel()) == "Trento") {
          entry->SetLabel("T_{R}ENTo");
        }
      }
      // Legend position in NDC of the upper pad (x1, y1, x2, y2); the values are the
      // originals from the file, change them to move or resize the legend
      padTop->cd();
      TLatex* lat = new TLatex();
      lat->SetTextSize(0.09);
      if (isTrento) {
        lat->SetTextAlign(31);
        const float legOffsetY = -0.12;
        setLegendPosition(leg, 0.62, 0.7 + legOffsetY, 0.95, 0.923 + legOffsetY);
        const float latPosX = 0.94;
        const float latPosY = 0.855;
        lat->DrawLatexNDC(latPosX, latPosY, "ALICE OO #sqrt{#it{s}_{NN}} = 5.36 TeV");
      } else {
        hFrame->SetMaximum(hFrame->GetMaximum() * 1.14);
        setLegendPosition(leg, 0.72, 0.2345, 0.95, 0.943);
        const float latPosX = 0.165;
        const float latPosY = 0.855;
        lat->DrawLatexNDC(latPosX, latPosY, "ALICE OO #sqrt{#it{s}_{NN}} = 5.36 TeV");
      }
      // leg->SetTextSize(0.05);
      padTop->Modified();

      // Thicker lines for the nuclear-structure models (curves, ratios, bands and the
      // legend entries, which keep their own copy of the line attributes)
      if (!isTrento) {
        const int lineWidth = 2;
        for (TPad* pad : { padTop, padBot }) {
          for (TH1* h : plot::getHistograms(pad)) {
            h->SetLineWidth(lineWidth);
          }
          pad->Modified();
        }
        for (const auto& obj : *leg->GetListOfPrimitives()) {
          ((TLegendEntry*)obj)->SetLineWidth(lineWidth);
        }
      }

      // Lower pad: ratio to Oho2, omega = 0.3
      TH1* hRatio = plot::getFrame(padBot);
      hRatio->GetXaxis()->SetTitle("FT0C Centrality (%)");
      hRatio->GetYaxis()->SetTitle("Ratio to Oho2, #omega = 0.3");
      // hRatio->GetYaxis()->SetRangeUser(0.9, 1.1);
      if (isTrento && isNpart) {
        hRatio->SetMaximum(1.039); // the 1.04 label would be cut at the pad edge
      }
      plot::drawHorizontalLine(padBot, 1.); // line at unity, behind the ratios
      padBot->Modified();


      c->Modified();
      c->Update();
      c->SaveAs(Form("Figures/%s_%s_OO.pdf", var, tag));
    }
  }
}
