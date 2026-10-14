#include <iostream>
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include <vector>
#include "Includes/colorManager.h"

using namespace std;

// Compares Glauber+NBD and TRENTo fits to the same data.
// Inputs are calibration files (output of Glauber/runCalibration.cc and
// Trento/runTrentoCalib.cc), which contain hData and the fitted model curve.
// The first enabled entry in the list provides the data histogram that is drawn.
enum class ModelType { Glauber,
                       Trento };

struct ModelInput {
  Bool_t enabled;
  TString fileName;
  TString label;
  ModelType type;
  Int_t color = kBlack; // assigned from the palette of its model type
};

void drawGlauberVsTrento()
{
  // ---- switch groups of histograms on/off ----
  Bool_t lDrawGlauber = kTRUE; // all entries with ModelType::Glauber
  Bool_t lDrawTrento = kTRUE;  // all entries with ModelType::Trento
  Bool_t lDrawPythia = kTRUE;  // scaled PYTHIA 8 Angantyr N_ch
  Bool_t lDrawNGlo = kFALSE;    // scaled N global tracks
  Bool_t lRebinToMult = kFALSE;
  Int_t lLineWidth = 1;

  // ---- fits to compare: {enabled, calibration file, legend label, model type} ----
  // set the first field to kFALSE to drop a single entry
  // colors are assigned automatically: one palette per model type (see below)
  std::vector<ModelInput> lModels = {
    { kTRUE, "../AnalysisResults/LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root", "Glauber mode 2", ModelType::Glauber },
    { kFALSE, "../AnalysisResults/LHC23_pass5/AR_544122_calibration_ancestorMode1_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root", "Glauber mode 1", ModelType::Glauber },
    { kFALSE, "../AnalysisResults/LHC23_pass5/AR_544122_calibration_ancestorMode0_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root", "Glauber mode 0", ModelType::Glauber },
    { kTRUE, "../Trento/trentofit_Jesper2_calib.root", "T_{R}ENTo v1 params", ModelType::Trento },
    { kTRUE, "../Trento/trentofit_New_calib.root", "T_{R}ENTo v2 params", ModelType::Trento },
  };
  TString lRunNumber = "544122";
  TString lOutputName = "GlauberVsTrento.pdf";
  TString lXaxisTitle = "FT0C Amplitude";

  // ---- normalisation: curves are scaled to the plotted data in [lNormMin, lNormMax] ----
  // PYTHIA/NGlo are always normalised, the model fits only if lNormaliseFits
  // (otherwise they keep the normalisation from the fit)
  Double_t lNormMin = 15000;
  Double_t lNormMax = 35000;
  Bool_t lNormaliseFits = kTRUE;

  // ---- multiplicity distributions mapped onto the FT0C axis ----
  // x: multiplied by lScaleTracksToFT0C (from calculateNGloToFT0C.cc, same factor for both)
  Double_t lScaleTracksToFT0C = 21.3344;

  // ---- rebin data and fits to roughly the PYTHIA/NGlo bin width ----
  // one N_ch unit spans lScaleTracksToFT0C in FT0C amplitude; contents are scaled by
  // 1/rebin so everything stays in counts per original data bin
  // lRebinToMult only affects the top pad; Data/PYTHIA and Data/NGlo ratios always use
  // data rebinned by lRebinToMultFactor (Data/Fit ratios use rebinRatio)
  Int_t lRebinToMultFactor = -1; // <0: use the nearest integer to lScaleTracksToFT0C

  TString lPythiaFileName = "../Trento/pythia/events_pythia8_angantyr_FT0.root";
  TString lPythiaHistName = "hChargedPartMid";
  TString lPythiaLabel = "PYTHIA N_{ch} (scaled)";

  TString lNGloFileName = "../AnalysisResults/LHC23_pass5/AR_544122.root";
  TString lNGloHistName = "centrality-study/hNGlobalTracks";
  TString lNGloLabel = "N_{global tracks} (scaled)";

  // ---- palettes (from Includes/colorManager.h): one per group, dark -> light within a group ----
  ColorManager::Palette lPaletteGlauber = ColorManager::Palette::Blue;
  ColorManager::Palette lPaletteTrento = ColorManager::Palette::Orange;
  ColorManager::Palette lPaletteMult = ColorManager::Palette::Green; // PYTHIA and NGlo

  // keep only the enabled entries
  std::vector<ModelInput> lAllModels = lModels;
  lModels.clear();
  for (auto& model : lAllModels) {
    Bool_t lGroupEnabled = (model.type == ModelType::Glauber) ? lDrawGlauber : lDrawTrento;
    if (model.enabled && lGroupEnabled)
      lModels.push_back(model);
  }
  if (lModels.empty()) {
    cout << "No fits enabled: at least one is needed to provide the data histogram" << endl;
    return;
  }

  const Int_t nModels = lModels.size();

  Int_t nGlauber = 0, nTrento = 0;
  for (auto& model : lModels)
    (model.type == ModelType::Glauber ? nGlauber : nTrento)++;
  ColorManager lColorsGlauber(TMath::Max(1, nGlauber), lPaletteGlauber);
  ColorManager lColorsTrento(TMath::Max(1, nTrento), lPaletteTrento);
  ColorManager lColorsMult(2, lPaletteMult);
  Int_t iGlauber = 0, iTrento = 0;
  for (auto& model : lModels) {
    model.color = (model.type == ModelType::Glauber) ? lColorsGlauber.getColor(iGlauber++) : lColorsTrento.getColor(iTrento++);
  }
  Int_t lPythiaColor = lColorsMult.getColor(0);
  Int_t lNGloColor = lColorsMult.getColor(1);

  int rebinRatio = 20;
  Double_t lFitRangeMax = 54000;

  std::vector<TH1F*> hData(nModels), hFit(nModels), hRatio(nModels), hRatioFine(nModels);

  // sum of the bin contents in [lNormMin, lNormMax]
  auto IntegralInNormRange = [&](TH1* h) -> Double_t {
    return h->Integral(h->FindBin(lNormMin), h->FindBin(lNormMax));
  };

  for (Int_t im = 0; im < nModels; im++) {
    TFile* file = new TFile(lModels[im].fileName, "READ");
    if (!file || file->IsZombie()) {
      cout << "Could not open " << lModels[im].fileName << endl;
      return;
    }
    file->ls();

    hData[im] = (TH1F*)file->Get("hData");
    if (!hData[im])
      hData[im] = (TH1F*)file->Get("hV0MUltraFine");

    hFit[im] = (TH1F*)file->Get("hGlauber");
    if (!hFit[im])
      hFit[im] = (TH1F*)file->Get("hGlauberFine");

    if (!hData[im] || !hFit[im]) {
      cout << "Data or fit histogram not found in " << lModels[im].fileName << endl;
      return;
    }
    hData[im]->SetName(Form("hData_%i", im));
    hFit[im]->SetName(Form("hFit_%i", im));

    // before the ratios are made, so Data/Model shows the normalised fit
    // (both are in counts per their own bin, so the integrals compare directly)
    if (lNormaliseFits) {
      Double_t lNormFactor = IntegralInNormRange(hData[0]) / IntegralInNormRange(hFit[im]);
      cout << lModels[im].label << ": fit scaled by " << lNormFactor << endl;
      hFit[im]->Scale(lNormFactor);
    }

    // ratio is computed against the data stored in the same file, so binning always matches
    hRatio[im] = (TH1F*)hData[im]->Clone(Form("hRatio_%i", im));
    TH1F* hFitRebinned = (TH1F*)hFit[im]->Clone(Form("hFitRebinned_%i", im));
    // rebin both models to the same final bin width as the reference (first entry) data
    Int_t lRebin = TMath::Max(1, TMath::Nint(rebinRatio * hData[0]->GetBinWidth(1) / hData[im]->GetBinWidth(1)));
    Int_t lRebinFit = TMath::Max(1, TMath::Nint(lRebin * hData[im]->GetBinWidth(1) / hFit[im]->GetBinWidth(1)));
    hRatio[im]->Rebin(lRebin);
    hFitRebinned->Rebin(lRebinFit);
    hRatio[im]->Divide(hFitRebinned);

    // unrebinned ratio for the zoomed inset (fit only rebinned if it is finer than the data)
    hRatioFine[im] = (TH1F*)hData[im]->Clone(Form("hRatioFine_%i", im));
    TH1F* hFitFine = (TH1F*)hFit[im]->Clone(Form("hFitFine_%i", im));
    hFitFine->Rebin(TMath::Max(1, TMath::Nint(hData[im]->GetBinWidth(1) / hFit[im]->GetBinWidth(1))));
    hRatioFine[im]->Divide(hFitFine);

    // express fit in counts per bin of the plotted (first entry) data histogram
    hFit[im]->Scale(hData[0]->GetBinWidth(1) / hFit[im]->GetBinWidth(1));
    hFit[im]->SetLineColor(lModels[im].color);
    hFit[im]->SetMarkerColor(lModels[im].color);
    hFit[im]->SetLineWidth(lLineWidth);
  }

  // stretch x-axis by lScaleTracksToFT0C, then normalise to the plotted data in [lNormMin, lNormMax]
  auto GetScaledToFT0C = [&](TString lFileName, TString lHistName, TString lNewName, Int_t lColor) -> TH1F* {
    TFile* fileIn = new TFile(lFileName, "READ");
    TH1* hRaw = fileIn->IsZombie() ? 0x0 : (TH1*)fileIn->Get(lHistName);
    if (!hRaw) {
      cout << lHistName << " not found in " << lFileName << endl;
      return 0x0;
    }
    const Int_t nBins = hRaw->GetNbinsX();
    TH1F* hScaled = new TH1F(lNewName, "", nBins,
                             hRaw->GetXaxis()->GetXmin() * lScaleTracksToFT0C,
                             hRaw->GetXaxis()->GetXmax() * lScaleTracksToFT0C);
    for (Int_t ii = 1; ii <= nBins; ii++)
      hScaled->SetBinContent(ii, hRaw->GetBinContent(ii));

    // contents in counts per bin of the plotted data
    hScaled->Scale(IntegralInNormRange(hData[0]) / IntegralInNormRange(hScaled) * hData[0]->GetBinWidth(1) / hScaled->GetBinWidth(1));

    hScaled->SetLineColor(lColor);
    hScaled->SetLineWidth(lLineWidth);
    // hScaled->SetLineStyle(2);
    return hScaled;
  };

  TH1F* hPythia = lDrawPythia ? GetScaledToFT0C(lPythiaFileName, lPythiaHistName, "hPythia", lPythiaColor) : 0x0;
  TH1F* hNGlo = lDrawNGlo ? GetScaledToFT0C(lNGloFileName, lNGloHistName, "hNGlo", lNGloColor) : 0x0;

  TH1F* hV0M = hData[0];

  if (lRebinToMultFactor < 0)
    lRebinToMultFactor = TMath::Max(1, TMath::Nint(lScaleTracksToFT0C));

  // Data/PYTHIA and Data/NGlo: always w.r.t. data rebinned to the multiplicity bin width.
  // Both are in counts per original data bin; the multiplicity curve is evaluated at the
  // centre of each rebinned data bin since its stretched bins do not align with the data bins
  TH1F* hDataMult = (TH1F*)hData[0]->Clone("hDataMult");
  hDataMult->Rebin(lRebinToMultFactor);
  hDataMult->Scale(1. / lRebinToMultFactor);
  auto GetRatioToMult = [&](TH1F* hMult, TString lNewName) -> TH1F* {
    if (!hMult)
      return 0x0;
    TH1F* hRatioMult = (TH1F*)hDataMult->Clone(lNewName);
    for (Int_t ii = 1; ii <= hRatioMult->GetNbinsX(); ii++) {
      Double_t lMult = hMult->Interpolate(hRatioMult->GetBinCenter(ii));
      Double_t lRatio = lMult > 0 ? hDataMult->GetBinContent(ii) / lMult : 0;
      Double_t lError = lMult > 0 ? hDataMult->GetBinError(ii) / lMult : 0;
      hRatioMult->SetBinContent(ii, lRatio);
      hRatioMult->SetBinError(ii, lError);
    }
    hRatioMult->SetLineColor(hMult->GetLineColor());
    hRatioMult->SetMarkerColor(hMult->GetLineColor());
    hRatioMult->SetLineStyle(hMult->GetLineStyle());
    hRatioMult->SetLineWidth(lLineWidth);
    return hRatioMult;
  };
  TH1F* hRatioPythia = GetRatioToMult(hPythia, "hRatioPythia");
  TH1F* hRatioNGlo = GetRatioToMult(hNGlo, "hRatioNGlo");

  if (lRebinToMult) {
    cout << "Rebinning data and fits by " << lRebinToMultFactor << " (scaled by 1/" << lRebinToMultFactor << ")" << endl;
    std::vector<TH1F*> lToRebin = hFit;
    lToRebin.push_back(hV0M); // keep last: its original bin width is the reference above
    for (auto h : lToRebin) {
      // same final bin width for all; a fit with different binning gets its own factor
      Int_t lRebinHist = TMath::Max(1, TMath::Nint(lRebinToMultFactor * hV0M->GetBinWidth(1) / h->GetBinWidth(1)));
      h->Rebin(lRebinHist);
      h->Scale(1. / lRebinHist);
    }
  }

  // transparent background outside the axes; insets keep an opaque frame so the
  // curves of the pad underneath do not show through the zoomed plot
  auto MakeTransparent = [](TVirtualPad* pad, Bool_t lOpaqueFrame = kFALSE) {
    pad->SetFillStyle(0);
    pad->SetFillColorAlpha(0, 0);
    pad->SetFrameFillStyle(lOpaqueFrame ? 1001 : 0);
    pad->SetFrameFillColor(kWhite);
    pad->SetBorderMode(0);
    pad->SetFrameBorderMode(0);
  };

  TCanvas* c1 = new TCanvas("c1", "", 900, 900);
  MakeTransparent(c1);
  c1->Divide(1, 2);
  MakeTransparent(c1->cd(1));
  MakeTransparent(c1->cd(2));

  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1, 1);
  c1->cd(1)->SetPad(0, 0.5, 1, 1);
  c1->cd(2)->SetPad(0, 0.0, 1, .5);

  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.05);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.08);

  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.05);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.08);
  c1->cd(2)->SetTicks(1, 1);
  c1->cd(1);

  Double_t lYminTop = 0.25;
  hV0M->GetXaxis()->SetRangeUser(0, lFitRangeMax);
  hV0M->GetYaxis()->SetRangeUser(lYminTop, hV0M->GetMaximum() * 20);
  hV0M->SetLineColor(kBlack);
  hV0M->SetMarkerStyle(20);
  hV0M->SetMarkerColor(kBlack);
  hV0M->SetMarkerSize(0.5);
  hV0M->GetYaxis()->SetTitleSize(0.07);
  hV0M->GetYaxis()->SetLabelSize(0.05);
  hV0M->GetYaxis()->SetTitle("Count");
  hV0M->GetYaxis()->SetTitleOffset(0.9);
  hV0M->GetXaxis()->SetLabelSize(0.05);
  hV0M->GetXaxis()->SetTitleSize(0.06);
  hV0M->GetXaxis()->SetTitle(lXaxisTitle.Data());
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->SetStats(0);
  hV0M->SetLineWidth(lLineWidth);
  hV0M->SetTitle("");
  hV0M->Draw("hist");
  for (Int_t im = 0; im < nModels; im++)
    hFit[im]->Draw("same hist");
  if (hPythia)
    hPythia->Draw("same hist");
  if (hNGlo)
    hNGlo->Draw("same hist");

  TLatex* latPer = new TLatex();
  latPer->SetNDC();
  latPer->SetTextFont(42);
  latPer->SetTextSize(0.06);
  latPer->DrawLatex(0.12, 0.88, "ALICE Pb-Pb 5.36 TeV");
  // legend height grows with the number of entries
  Int_t nLegEntries = 1 + nModels + (hPythia ? 1 : 0) + (hNGlo ? 1 : 0);
  Double_t lLegBottom = 0.85 - 0.04 * nLegEntries;
  latPer->SetTextSize(0.05);
  // latPer -> DrawLatex(0.12, lLegBottom - 0.07, Form("Run %s", lRunNumber.Data()));

  TLegend* legO = new TLegend(0.12, lLegBottom, 0.42, 0.85);
  legO->SetTextSize(0.035);
  legO->SetBorderSize(0);
  legO->SetFillColorAlpha(0, 0);
  legO->AddEntry(hV0M, "Data", "l");
  for (Int_t im = 0; im < nModels; im++)
    legO->AddEntry(hFit[im], lModels[im].label.Data(), "l");
  if (hPythia)
    legO->AddEntry(hPythia, lPythiaLabel.Data(), "l");
  if (hNGlo)
    legO->AddEntry(hNGlo, lNGloLabel.Data(), "l");
  legO->Draw();

  // inset: zoom on peripheral region
  TPad* inset1 = new TPad("inset1", "", 0.575724, 0.441648, 0.924, 0.961098);
  inset1->SetPad(0.575724, 0.441648, 0.924, 0.961098);
  MakeTransparent(inset1, kTRUE);
  inset1->cd();

  // Double_t lXmaxInset1 = 590*lFitRangeMax/59000.;
  Double_t lXmaxInset1 = 300;
  Double_t lYminInset1 = 2e+2;
  Double_t lYmaxInset1 = 2e+7;
  TH1F* hV0Mzoomed = (TH1F*)hV0M->Clone("hV0Mzoomed");
  hV0Mzoomed->GetXaxis()->SetRangeUser(0, lXmaxInset1);
  hV0Mzoomed->GetXaxis()->SetNdivisions(404);
  hV0Mzoomed->GetXaxis()->SetTitleSize(0.1);
  hV0Mzoomed->GetYaxis()->UnZoom();
  hV0Mzoomed->SetMinimum(lYminInset1);
  hV0Mzoomed->SetMaximum(lYmaxInset1);
  hV0Mzoomed->GetYaxis()->SetTitleSize(0.1);
  hV0Mzoomed->GetYaxis()->SetTitleOffset(1.1);
  hV0Mzoomed->GetXaxis()->SetLabelSize(0.06);
  hV0Mzoomed->GetYaxis()->SetLabelSize(0.06);
  hV0Mzoomed->Draw("hist");
  for (Int_t im = 0; im < nModels; im++)
    hFit[im]->Draw("same hist");
  if (hPythia)
    hPythia->Draw("same hist");
  if (hNGlo)
    hNGlo->Draw("same hist");

  c1->cd(1);

  inset1->SetLogy();
  inset1->SetTicks(1, 1);
  inset1->SetBottomMargin(0.194);
  inset1->SetRightMargin(0.03);
  inset1->SetTopMargin(0.085);
  inset1->SetLeftMargin(0.11);

  inset1->Draw("same");

  c1->cd(2);
  Float_t lLoRangeRatio = 0.3;
  Float_t lHiRangeRatio = 2.495;

  for (Int_t im = 0; im < nModels; im++) {
    hRatio[im]->GetYaxis()->SetTitle("Data/Model");
    hRatio[im]->GetXaxis()->SetTitle(lXaxisTitle.Data());
    hRatio[im]->GetYaxis()->SetTitleSize(0.055);
    hRatio[im]->GetYaxis()->SetTitleOffset(0.7);
    hRatio[im]->GetXaxis()->SetTitleSize(0.055);
    hRatio[im]->GetYaxis()->SetLabelSize(0.045);
    hRatio[im]->GetXaxis()->SetLabelSize(0.045);
    hRatio[im]->GetYaxis()->SetRangeUser(lLoRangeRatio, lHiRangeRatio);
    hRatio[im]->GetXaxis()->SetRangeUser(0, lFitRangeMax);
    hRatio[im]->SetMarkerStyle(20);
    hRatio[im]->SetMarkerColor(lModels[im].color);
    hRatio[im]->SetLineColor(lModels[im].color);
    hRatio[im]->SetMarkerSize(.5);
    hRatio[im]->SetLineWidth(lLineWidth);
    hRatio[im]->SetStats(0);
    hRatio[im]->SetTitle("");
  }

  hRatio[0]->Draw("hist");

  TLine* line = new TLine(0, 1, lFitRangeMax, 1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray + 1);
  line->Draw();

  for (Int_t im = 0; im < nModels; im++)
    hRatio[im]->Draw("same hist");
  if (hRatioPythia)
    hRatioPythia->Draw("same hist");
  if (hRatioNGlo)
    hRatioNGlo->Draw("same hist");

  // inset: zoom on peripheral region of the ratio
  Double_t lXmaxInset2 = 1150 * lFitRangeMax / 59000.;
  TPad* inset2 = new TPad("inset2", "", 0.104677, 0.47032, 0.558775, 0.965753);
  inset2->SetPad(0.104677, 0.47032, 0.558775, 0.965753);
  MakeTransparent(inset2, kTRUE);
  inset2->SetTicks(1, 1);
  inset2->SetBottomMargin(0.194);
  inset2->SetRightMargin(0.03);
  inset2->SetTopMargin(0.045);
  inset2->SetLeftMargin(0.11);
  inset2->cd();

  std::vector<TH1F*> hRatioZoomed(nModels);
  for (Int_t im = 0; im < nModels; im++) {
    hRatioZoomed[im] = (TH1F*)hRatioFine[im]->Clone(Form("hRatioZoomed_%i", im));
    hRatioZoomed[im]->SetMarkerStyle(20);
    hRatioZoomed[im]->SetMarkerColor(lModels[im].color);
    hRatioZoomed[im]->SetLineColor(lModels[im].color);
    hRatioZoomed[im]->SetMarkerSize(.5);
    hRatioZoomed[im]->SetLineWidth(lLineWidth);
    hRatioZoomed[im]->SetStats(0);
    hRatioZoomed[im]->SetTitle("");
    hRatioZoomed[im]->GetXaxis()->SetRangeUser(0, lXmaxInset2);
    hRatioZoomed[im]->GetYaxis()->SetRangeUser(0.4, 1.55);
    hRatioZoomed[im]->GetYaxis()->SetTitleSize(0.1);
    hRatioZoomed[im]->GetXaxis()->SetTitleSize(0.1);
    hRatioZoomed[im]->GetYaxis()->SetTitleOffset(1.1);
    hRatioZoomed[im]->GetXaxis()->SetLabelSize(0.06);
    hRatioZoomed[im]->GetYaxis()->SetLabelSize(0.06);
  }
  hRatioZoomed[0]->Draw("hist");

  TLine* line2 = new TLine(0, 1, lXmaxInset2, 1);
  line2->SetLineStyle(7);
  line2->SetLineColor(kGray + 1);
  line2->Draw();

  for (Int_t im = 0; im < nModels; im++)
    hRatioZoomed[im]->Draw("hist same");
  // multiplicity ratios stay at the rebinned-data granularity, also in the zoom
  if (hRatioPythia)
    hRatioPythia->Draw("hist same");
  if (hRatioNGlo)
    hRatioNGlo->Draw("hist same");

  c1->cd(2);
  inset2->Draw();

  c1->SaveAs(lOutputName.Data());

  // standalone version of the top-pad inset (peripheral region)
  TCanvas* c2 = new TCanvas("c2", "", 900, 700);
  MakeTransparent(c2);
  c2->SetLogy();
  c2->SetTicks(1, 1);
  c2->SetBottomMargin(0.12);
  c2->SetRightMargin(0.05);
  c2->SetTopMargin(0.04);
  c2->SetLeftMargin(0.12);

  TH1F* hV0MInset = (TH1F*)hV0Mzoomed->Clone("hV0MInset");
  hV0MInset->GetXaxis()->SetNdivisions(510);
  hV0MInset->GetXaxis()->SetTitleSize(0.045);
  hV0MInset->GetYaxis()->SetTitleSize(0.045);
  hV0MInset->GetXaxis()->SetLabelSize(0.04);
  hV0MInset->GetYaxis()->SetLabelSize(0.04);
  hV0MInset->GetXaxis()->SetTitleOffset(1.1);
  hV0MInset->GetYaxis()->SetTitleOffset(1.3);
  hV0MInset->GetYaxis()->SetTickLength(0.03);
  hV0MInset->Draw("hist");
  for (Int_t im = 0; im < nModels; im++)
    hFit[im]->Draw("same hist");
  if (hPythia)
    hPythia->Draw("same hist");
  if (hNGlo)
    hNGlo->Draw("same hist");

  TLatex* latInset = new TLatex();
  latInset->SetNDC();
  latInset->SetTextFont(42);
  latInset->SetTextSize(0.045);
  latInset->DrawLatex(0.16, 0.89, "ALICE Pb-Pb 5.36 TeV");
  latInset->SetTextSize(0.04);
  // latInset->DrawLatex(0.16, 0.84, Form("Run %s", lRunNumber.Data()));

  TLegend* legInset = new TLegend(0.50, 0.93 - 0.045 * nLegEntries, 0.92, 0.93);
  legInset->SetBorderSize(0);
  legInset->SetFillColorAlpha(0, 0);
  legInset->SetTextSize(0.04);
  legInset->AddEntry(hV0MInset, "Data", "l");
  for (Int_t im = 0; im < nModels; im++)
    legInset->AddEntry(hFit[im], lModels[im].label.Data(), "l");
  if (hPythia)
    legInset->AddEntry(hPythia, lPythiaLabel.Data(), "l");
  if (hNGlo)
    legInset->AddEntry(hNGlo, lNGloLabel.Data(), "l");
  legInset->Draw();

  TString lOutputNameInset = lOutputName;
  lOutputNameInset.ReplaceAll(".pdf", "_inset.pdf");
  c2->SaveAs(lOutputNameInset.Data());
}
