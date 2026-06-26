#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TGraph.h"

#include <iostream>
#include <vector>
#include <string>  

#include "colorManager.h"

void styleTH1(TH1F* hist, int color)
{
  hist->SetLineColor(color);
  hist->GetXaxis()->SetRangeUser(0, 14000);
  hist->GetXaxis()->SetTitle("FT0M Amplitude");
  hist->GetYaxis()->SetTitle("Centrality precentile");
}

void styleTGraph(TGraph* gr, int color)
{
  gr->SetMarkerStyle(kFullCircle);
  gr->SetMarkerColor(color);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Centrality percentile");
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitle("FT0M Amplitude");
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetRangeUser(0, 100);
}

void drawLightIonCalibration()
{
    gStyle->SetOptStat(0);
    TFile* fCalibColl = new TFile("../AnalysisResults/LHC25ae_pass2_small/AR_564374_calibration_hFT0M_Collisions_LightIonDef.root", "read");
    TFile* fCalibBCs = new TFile("../AnalysisResults/LHC25ae_pass2_small/AR_564374_calibration_hFT0M_BCs_LightIonDef.root", "read");
    TFile* fCalibAncColl = new TFile("../AnalysisResults/LHC25ae_pass2_small/AR_564374_calibration_hFT0M_Collisions_Anchor90.root", "read");
    TFile* fCalibAncBCs = new TFile("../AnalysisResults/LHC25ae_pass2_small/AR_564374_calibration_hFT0M_BCs_Anchor90.root", "read");

    TH1F* hCalibColl    = (TH1F*)fCalibColl->Get("hCalib");
    TH1F* hCalibBCs     = (TH1F*)fCalibBCs->Get("hCalib");
    TH1F* hCalibAncColl = (TH1F*)fCalibAncColl->Get("hCalib");
    TH1F* hCalibAncBCs  = (TH1F*)fCalibAncBCs->Get("hCalib");

    // --- Helper lambda to swap axes of a TH1F into a TH2F ---
    // The original TH1F has: X = FT0M amplitude, Y = calibration value (bin content)
    // The swapped TH2F has: X = calibration value range, Y = FT0M amplitude
    auto swapAxes = [](TH1F* h, const char* newName) -> TH2F* {
        int    nBins = h->GetNbinsX();
        double xMin  = h->GetXaxis()->GetXmin();
        double xMax  = h->GetXaxis()->GetXmax();
        double yMin  = h->GetMinimum();
        double yMax  = h->GetMaximum();

        TH2F* hSwap = new TH2F(newName, h->GetTitle(),
                               100, yMin, yMax,   // new X = old content (calib value)
                               nBins, xMin, xMax); // new Y = old X (FT0M amplitude)

        hSwap->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle());
        hSwap->GetYaxis()->SetTitle(h->GetXaxis()->GetTitle());

        for (int ix = 1; ix <= nBins; ix++) {
            double content = h->GetBinContent(ix);
            double error   = h->GetBinError(ix);
            if (content == 0) continue;
            double oldX = h->GetBinCenter(ix);
            int newBinX = hSwap->GetXaxis()->FindBin(content);
            int newBinY = hSwap->GetYaxis()->FindBin(oldX);
            hSwap->SetBinContent(newBinX, newBinY, 1); // mark the point
            (void)error; // preserve if needed
        }
        return hSwap;
    };

    // --- Alternatively: swap into a TGraph (cleaner for a calibration curve) ---
    auto swapToGraph = [](TH1F* h, const char* newName) -> TGraph* {
        int nPoints = 0;
        for (int ix = 1; ix <= h->GetNbinsX(); ix++)
            if (h->GetBinContent(ix) != 0) nPoints++;

        TGraph* g = new TGraph(nPoints);
        g->SetName(newName);
        g->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle()); // old content -> new X
        g->GetYaxis()->SetTitle(h->GetXaxis()->GetTitle()); // old X -> new Y

        int iPoint = 0;
        for (int ix = 1; ix <= h->GetNbinsX(); ix++) {
            double content = h->GetBinContent(ix);
            if (content == 0) continue;
            double oldX = h->GetBinCenter(ix);
            g->SetPoint(iPoint++, content, oldX); // swap: content -> X, oldX -> Y
        }
        return g;
    };

    TGraph* gSwapColl    = swapToGraph(hCalibColl,    "gSwapColl");
    TGraph* gSwapBCs     = swapToGraph(hCalibBCs,     "gSwapBCs");
    TGraph* gSwapAncColl = swapToGraph(hCalibAncColl, "gSwapAncColl");
    TGraph* gSwapAncBCs  = swapToGraph(hCalibAncBCs,  "gSwapAncBCs");

    ColorManager cm(4);
    styleTH1(hCalibColl,    cm.getColor(0));
    styleTH1(hCalibBCs,     cm.getColor(1));
    styleTH1(hCalibAncColl, cm.getColor(2));
    styleTH1(hCalibAncBCs,  cm.getColor(3));

    gSwapColl->SetLineColor(cm.getColor(0));    gSwapColl->SetLineWidth(2);
    gSwapBCs->SetLineColor(cm.getColor(1));     gSwapBCs->SetLineWidth(2);
    gSwapAncColl->SetLineColor(cm.getColor(2)); gSwapAncColl->SetLineWidth(2);
    gSwapAncBCs->SetLineColor(cm.getColor(3));  gSwapAncBCs->SetLineWidth(2);

    // --- Original plot ---
    TCanvas* canvCalib = new TCanvas("canvCalib", "", 1200, 800);
    canvCalib->SetTicks(1, 1);
    canvCalib->SetTopMargin(0.03);
    canvCalib->SetLeftMargin(0.13);
    canvCalib->SetRightMargin(0.03);
    canvCalib->SetBottomMargin(0.13);
    hCalibColl->Draw();
    hCalibBCs->Draw("same");
    hCalibAncColl->Draw("same");
    hCalibAncBCs->Draw("same");

    TLegend* leg = new TLegend(0.6, 0.67, 0.95, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0, 0);
    leg->AddEntry(hCalibColl,    "Collision data slicing",  "l");
    leg->AddEntry(hCalibAncColl, "Collision data 90% anchor", "l");
    leg->AddEntry(hCalibBCs,     "BC data slicing",         "l");
    leg->AddEntry(hCalibAncBCs,  "BC data 90% anchor",      "l");
    leg->Draw();
    canvCalib->SaveAs("hLightIonCalibration.pdf");

TCanvas* canvSwap = new TCanvas("canvSwap", "", 1200, 1000);

// --- Pads ---
TPad* padMain  = new TPad("padMain",  "", 0, 0.35, 1, 1.0);
TPad* padRatio = new TPad("padRatio", "", 0, 0.00, 1, 0.35);

padMain->SetTopMargin(0.04);
padMain->SetBottomMargin(0);
padMain->SetLeftMargin(0.13);
padMain->SetRightMargin(0.03);
padMain->SetTicks(1, 1);
padMain->SetLogy();
padMain->Draw();

padRatio->SetTopMargin(0);
padRatio->SetBottomMargin(0.25);
padRatio->SetLeftMargin(0.13);
padRatio->SetRightMargin(0.03);
padRatio->SetTicks(1, 1);
padRatio->Draw();

// --- Main pad ---
padMain->cd();

styleTGraph(gSwapColl,    cm.getColor(0));
styleTGraph(gSwapBCs,     cm.getColor(1));
styleTGraph(gSwapAncColl, cm.getColor(2));
styleTGraph(gSwapAncBCs,  cm.getColor(3));

gSwapColl->Draw("AP");
gSwapColl->GetXaxis()->SetLabelSize(0);  // hide X labels on main pad
gSwapColl->GetYaxis()->SetTitle("FT0M amplitude");
gSwapColl->GetYaxis()->SetTitleSize(0.05);
gSwapColl->GetYaxis()->SetTitleOffset(1.2);
gSwapBCs->Draw("P same");
gSwapAncColl->Draw("P same");
gSwapAncBCs->Draw("P same");

TLegend* legSwap = new TLegend(0.6, 0.65, 0.95, 0.9);
legSwap->SetBorderSize(0);
legSwap->SetFillColorAlpha(0, 0);
legSwap->AddEntry(gSwapColl,    "Collision data slicing",    "p");
legSwap->AddEntry(gSwapAncColl, "Collision data 90% anchor", "p");
legSwap->AddEntry(gSwapBCs,     "BC data slicing",           "p");
legSwap->AddEntry(gSwapAncBCs,  "BC data 90% anchor",        "p");
legSwap->Draw();

// --- Build ratio graphs: gSwapBCs/gSwapColl, gSwapAncColl/gSwapColl, gSwapAncBCs/gSwapColl ---
auto makeRatio = [](TGraph* gNum, TGraph* gDen, const char* name) -> TGraph* {
    int n = gNum->GetN();
    TGraph* gRatio = new TGraph(n);
    gRatio->SetName(name);
    for (int i = 0; i < n; i++) {
        double x, yNum, yDen;
        gNum->GetPoint(i, x, yNum);
        // Evaluate denominator at the same x by finding closest point
        yDen = gDen->Eval(x);
        double ratio = (yDen != 0) ? yNum / yDen : 0;
        gRatio->SetPoint(i, x, ratio);
    }
    return gRatio;
};

TGraph* gRatioBCs     = makeRatio(gSwapBCs,     gSwapColl, "gRatioBCs");
TGraph* gRatioAncColl = makeRatio(gSwapAncColl, gSwapColl, "gRatioAncColl");
TGraph* gRatioAncBCs  = makeRatio(gSwapAncBCs,  gSwapColl, "gRatioAncBCs");

styleTGraph(gRatioBCs,     cm.getColor(1));
styleTGraph(gRatioAncColl, cm.getColor(2));
styleTGraph(gRatioAncBCs,  cm.getColor(3));

// --- Ratio pad ---
padRatio->cd();

// Reference line at 1
double xMin = gSwapColl->GetXaxis()->GetXmin();
double xMax = gSwapColl->GetXaxis()->GetXmax();
TLine* lineUnity = new TLine(xMin, 1.0, 100, 1.0);
lineUnity->SetLineColor(cm.getColor(0));
lineUnity->SetLineStyle(2);
lineUnity->SetLineWidth(2);

gRatioBCs->Draw("AP");
gRatioBCs->GetXaxis()->SetTitle("Centrality percentile");
gRatioBCs->GetYaxis()->SetTitle("Ratio to coll. slicing");
gRatioBCs->GetYaxis()->SetNdivisions(505);
gRatioBCs->GetYaxis()->SetTitleSize(0.09);
gRatioBCs->GetYaxis()->SetTitleOffset(0.6);
gRatioBCs->GetYaxis()->SetLabelSize(0.08);
gRatioBCs->GetYaxis()->SetRangeUser(0.45, 1.05);
gRatioBCs->GetXaxis()->SetRangeUser(0, 100);
gRatioBCs->GetXaxis()->SetTitleSize(0.09);
gRatioBCs->GetXaxis()->SetTitleOffset(1.1);
gRatioBCs->GetXaxis()->SetLabelSize(0.08);
gRatioAncColl->Draw("P same");
gRatioAncBCs->Draw("P same");
lineUnity->Draw();

canvSwap->cd();
canvSwap->SaveAs("hLightIonCalibration_swapped.pdf");
}