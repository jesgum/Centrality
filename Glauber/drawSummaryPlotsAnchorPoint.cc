TH1D *GetCentralityHisto(TH1F *hMult, TH1 *hCalib, Double_t lCent1, Double_t lCent2, TString lName = "h"){
  TH1D *h = (TH1D*) hMult->Clone(lName.Data());
  
  //cout<<"lLo = "<<lLo<<", lHi = "<<lHi<<endl;
  for(Int_t ii=1; ii<hMult->GetNbinsX()+1; ii++) {
    float centrality = hCalib->GetBinContent( hCalib->FindBin(hMult->GetBinCenter(ii)));
    if( centrality < lCent1 ) h->SetBinContent(ii, 0);
    if( centrality > lCent2 ) h->SetBinContent(ii, 0);
  }
  return h;
}

Double_t GetMidRange(TH1* histo, Double_t loweringEdge = 56000){
  // this calculates the average of rising edge vs lowering edge.
  Double_t risingEdge = -1;
  for(int ii=1; ii<histo->GetNbinsX(); ii++){
    if(histo->GetBinContent(ii)>0){
      risingEdge = histo->GetBinCenter(ii);
      break;
    }
  }
  for(int ii=histo->FindBin(loweringEdge); ii>0; ii--){
    if(histo->GetBinContent(ii)>0){
      loweringEdge = histo->GetBinCenter(ii);
      break;
    }
  }
  return 0.5*(risingEdge+loweringEdge);
}

Double_t GetBoundaryForPercentile( TH1 *histo, Double_t lPercentileRequested ) {
  //This function returns the boundary for a specific percentile.
  Double_t lReturnValue = 0.0;
  Double_t lPercentile = 100.0 - lPercentileRequested;
  
  const Long_t lNBins = histo->GetNbinsX();
  Double_t lCountDesired = lPercentile * histo->GetEntries()/100;
  Long_t lCount = 0;
  for(Long_t ibin=1;ibin<lNBins;ibin++){
    lCount += histo->GetBinContent(ibin);
    if(lCount >= lCountDesired){
      //Found bin I am looking for!
      Double_t lWidth = histo->GetBinWidth(ibin);
      Double_t lLeftPercentile = 100.*(lCount - histo->GetBinContent(ibin))/histo->GetEntries();
      Double_t lRightPercentile = 100.*lCount / histo->GetEntries();
      
      Double_t lProportion = (lPercentile - lLeftPercentile)/(lRightPercentile-lLeftPercentile);
      
      lReturnValue = histo->GetBinLowEdge(ibin) + lProportion*lWidth;
      break;
    }
  }
  return lReturnValue;
}

Double_t GetAnchorPoint(TH1F* hGlauber_, Double_t anchorPointPercentage) {
  // Same approach as Trento/runTrentoFit.cc: build the normalised cumulative
  // distribution of the Glauber + NBD model alone and find where it crosses
  // the fraction of events below the anchor point (e.g. 10% for a 90% anchor)
  const Double_t lFractionBelow = 1. - anchorPointPercentage / 100.;
  const Int_t lNBins = hGlauber_->GetNbinsX();
  const Double_t lTotal = hGlauber_->Integral(1, lNBins);
  if (lTotal <= 0) {
    cout << "Glauber histogram " << hGlauber_->GetName() << " is empty, cannot determine anchor point!" << endl;
    return -1;
  }

  // Cumulative evaluated at the upper bin edges, starting from 0 at the lower edge of bin 1
  Double_t lCumuPrev = 0;
  for (Int_t ii = 1; ii < lNBins + 1; ii++) {
    Double_t lCumu = lCumuPrev + hGlauber_->GetBinContent(ii) / lTotal;
    if (lCumu > lFractionBelow) {
      // Linear interpolation between the bin edges
      Double_t lFrac = (lFractionBelow - lCumuPrev) / (lCumu - lCumuPrev);
      Double_t lAnchorPoint = hGlauber_->GetBinLowEdge(ii) + lFrac * hGlauber_->GetBinWidth(ii);
      cout << "Anchor point for " << hGlauber_->GetName() << " determined to be: " << lAnchorPoint << endl;
      return lAnchorPoint;
    }
    lCumuPrev = lCumu;
  }
  return -1;
}

void DrawAnchorLine(Double_t lAnchorPoint, Double_t lYLo, Double_t lYHi, Double_t lXMax, Int_t lColor){
  // Vertical dashed line at the anchor point, skipped if outside the drawn x range
  if(lAnchorPoint < 0 || lAnchorPoint > lXMax) return;
  TLine *l = new TLine(lAnchorPoint, lYLo, lAnchorPoint, lYHi);
  l->SetLineColor(lColor);
  l->SetLineStyle(2);
  l->SetLineWidth(2);
  l->Draw();
}

void drawSummaryPlotsAnchorPoint( TString lCalibFileName = "../AnalysisResults/LHC23_pass5/AR_544122_calibration_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb_Anchor90.0.root",
                                  TString lInputFileName = "../AnalysisResults/LHC23_pass5/AR_544122_glauberNBD_ancestorMode2_fixedK_fixedMu_hFT0C_BCs_basehistos_Pb.root")
{
  TString lXaxisTitle = "FT0C Amplitude";
  TString lRunNumber = "544122";
  TFile *file1 = new TFile(lInputFileName.Data(), "READ");
  TFile *file2 = new TFile(lCalibFileName.Data(), "READ");

  TH1F *hV0M = (TH1F*) file2->Get("hData");
  if(!hV0M) hV0M = (TH1F*) file2->Get("hV0MUltraFine");

  TH1F* hFittedGlauber = (TH1F*) file1->Get("hGlauber");
  TH1F* hCalibGlauber = (TH1F*) file2->Get("hGlauber");
  if(!hFittedGlauber || !hCalibGlauber){
    cout<<"Glauber histogram(s) not found! Stopping!"<<endl;
    return;
  }

  // hCalibGlauber is only hFittedGlauber times a global scale factor, so its normalised
  // cumulative (and thus anchor) is identical to the fit. What the calibration actually
  // uses is the stitched histogram (calib Glauber below the anchor, data above)
  TH1F* hStitched = (TH1F*) file2->Get("hStitched");
  TH1F* hStoredAnchor = (TH1F*) file2->Get("hAnchorPoint");
  const float anchorPointFit = GetAnchorPoint(hFittedGlauber, 90);
  const float anchorPointCalib = hStitched ? GetAnchorPoint(hStitched, 90) : -1;
  cout << "Fit anchor (Glauber + NBD fit only): "<<anchorPointFit<<endl;
  cout << "Calib anchor (stitched Glauber + data): "<<anchorPointCalib<<endl;
  if(hStoredAnchor) cout << "Calib anchor stored by runCalibration: "<<hStoredAnchor->GetBinContent(1)<<endl;


  hFittedGlauber->SetName("hFittedGlauber");
  hCalibGlauber->SetName("hCalibGlauber");

  // The fit is done on the rebinned hV0M, so the stored fit function histogram
  // is normalised per coarse bin: bring it to the fine binning of the data
  TH1F *hV0MCoarse = (TH1F*) file1->Get("hV0M");
  if(hV0MCoarse){
    Double_t lBinRatio = hV0MCoarse->GetBinWidth(1) / hV0M->GetBinWidth(1);
    cout<<"Scaling fitted Glauber by 1/"<<lBinRatio<<" (fit binning vs data binning)"<<endl;
    hFittedGlauber->Scale(1./lBinRatio);
  }

  int rebinRatio = 20;
  TH1F *hData = (TH1F*) hV0M->Clone("hDataRebinned");
  hData->Rebin(rebinRatio);
  TH1F *hRatioCalib = (TH1F*) hData->Clone("hRatioCalib");
  TH1F *hRatioFit = (TH1F*) hData->Clone("hRatioFit");
  TH1F *hTmp = (TH1F*) hCalibGlauber->Clone("hTmpCalib");
  hTmp->Rebin(rebinRatio);
  hRatioCalib->Divide(hTmp);
  hTmp = (TH1F*) hFittedGlauber->Clone("hTmpFit");
  hTmp->Rebin(rebinRatio);
  hRatioFit->Divide(hTmp);

  Double_t lFitRangeMax = GetBoundaryForPercentile(hV0M, 0.002);
  cout<<"Fit range max estimated from histogram: "<<lFitRangeMax<<endl;

  TCanvas *c1 = new TCanvas("c1", "", 900,900);
  c1->SetFrameFillStyle(0);
  c1->SetFillStyle(0);
  c1->Divide(1,2);

  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1,1);
  c1->cd(1)->SetPad(0,0.5,1,1);
  c1->cd(2)->SetPad(0,0.0,1,.5);

  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.05);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.08);

  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.05);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.08);
  c1->cd(2)->SetTicks(1,1);
  c1->cd(1);

  hV0M->GetXaxis()->SetRangeUser(0,lFitRangeMax);
  const Double_t lYLoMain = 0.25;
  const Double_t lYHiMain = hV0M->GetMaximum()*3;
  hV0M->GetYaxis()->SetRangeUser(lYLoMain,lYHiMain);
  hV0M->SetLineColor(kBlack);
  hV0M->SetMarkerStyle(20);
  hV0M->SetMarkerColor(kBlack);
  hV0M->SetMarkerSize(0.5);
  hV0M->GetYaxis()->SetTitleSize(0.07);
  hV0M->GetYaxis()->SetLabelSize(0.05);
  hV0M->GetYaxis()->SetTitle("Counts");
  hV0M->GetYaxis()->SetTitleOffset(0.9);
  hV0M->GetXaxis()->SetLabelSize(0.05);
  hV0M->GetXaxis()->SetTitleSize(0.06);
  hV0M->GetXaxis()->SetTitle(lXaxisTitle.Data());
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->SetStats(0);
  hV0M->SetTitle("");
  hV0M->Draw("hist");

  hFittedGlauber->SetLineColor(kBlue);
  hFittedGlauber->SetMarkerColor(kBlue);
  hCalibGlauber->SetLineColor(kRed);
  hCalibGlauber->SetMarkerColor(kRed);
  hFittedGlauber->Draw("same hist");
  hCalibGlauber->Draw("same hist");
  DrawAnchorLine(anchorPointFit, lYLoMain, lYHiMain, lFitRangeMax, kBlue);
  DrawAnchorLine(anchorPointCalib, lYLoMain, lYHiMain, lFitRangeMax, kRed);

  TLatex *latPer = new TLatex();
  latPer->SetNDC();
  latPer->SetTextFont(42);
  latPer->SetTextSize(0.06);
  latPer-> DrawLatex(0.12, 0.88, "ALICE Pb-Pb 5.36 TeV");
  latPer->SetTextSize(0.062);
  // latPer -> DrawLatex(0.12,0.62, Form("Run %s", lRunNumber.Data()));

  TLegend *legO = new TLegend(0.12, 0.53, 0.44, 0.85);
  legO->SetBorderSize(0);
  legO->SetFillColorAlpha(0, 0);
  legO->SetTextSize(0.048);
  legO->AddEntry(hV0M, "Data", "l");
  legO->AddEntry(hFittedGlauber, "Glauber + NBD fit", "l");
  legO->AddEntry(hCalibGlauber, "Glauber + NBD renormalised", "l");
  TLine *lLegAnchorFit = new TLine();
  lLegAnchorFit->SetLineColor(kBlue);
  lLegAnchorFit->SetLineStyle(2);
  lLegAnchorFit->SetLineWidth(2);
  TLine *lLegAnchorCalib = (TLine*) lLegAnchorFit->Clone();
  lLegAnchorCalib->SetLineColor(kRed);
  legO->AddEntry(lLegAnchorFit, Form("90%% anchor (fit): %.2f", anchorPointFit), "l");
  legO->AddEntry(lLegAnchorCalib, Form("90%% anchor (renormalised): %.2f", anchorPointCalib), "l");
  legO->Draw();

  // create extra inset pad
  TPad *inset1 = new TPad("inset1", "",0.575724, 0.441648, 0.924, 0.961098);
  inset1->SetPad(0.575724, 0.441648, 0.924, 0.961098);
  inset1->cd();

  TH1F *hV0Mzoomed = (TH1F*) hV0M->Clone("hV0Mzoomed");
  hV0Mzoomed -> GetXaxis()->SetRangeUser(0,590*lFitRangeMax/59000.);
  hV0Mzoomed -> GetXaxis()->SetNdivisions(404);
  hV0Mzoomed -> GetXaxis()->SetTitleSize(0.1);
  hV0Mzoomed -> GetYaxis()->UnZoom();
  const Double_t lYLoInset1 = hCalibGlauber->GetMaximum() * 0.05;
  const Double_t lYHiInset1 = hV0Mzoomed->GetMaximum() * 5.5;
  hV0Mzoomed -> GetYaxis()->SetRangeUser(lYLoInset1, lYHiInset1);
  hV0Mzoomed -> GetYaxis()->SetTitleSize(0.1);
  hV0Mzoomed -> GetYaxis()->SetTitleOffset(1.1);
  hV0Mzoomed -> GetXaxis()->SetLabelSize(0.06);
  hV0Mzoomed -> GetYaxis()->SetLabelSize(0.06);
  hV0Mzoomed->Draw();
  hFittedGlauber->Draw("same hist");
  hCalibGlauber->Draw("same hist");
  DrawAnchorLine(anchorPointFit, lYLoInset1, lYHiInset1, 590*lFitRangeMax/59000., kBlue);
  DrawAnchorLine(anchorPointCalib, lYLoInset1, lYHiInset1, 590*lFitRangeMax/59000., kRed);

  c1->cd(1);
  inset1->SetLogy();
  inset1->SetTicks(1,1);
  inset1->SetBottomMargin(0.194);
  inset1->SetRightMargin(0.03);
  inset1->SetTopMargin(0.085);
  inset1->SetLeftMargin(0.11);
  inset1->Draw("same");

  c1->cd(2);
  Float_t lLoRangeRatio = 0.62;
  Float_t lHiRangeRatio = 2.495;

  for(TH1F *h : {hRatioFit, hRatioCalib}){
    h->GetYaxis()->SetTitle("Data/Fit");
    h->GetXaxis()->SetTitle(lXaxisTitle.Data());
    h->GetYaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetTitleOffset(0.7);
    h->GetXaxis()->SetTitleSize(0.055);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetRangeUser(lLoRangeRatio,lHiRangeRatio);
    h->GetXaxis()->SetRangeUser(0,lFitRangeMax);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(.5);
    h->SetLineWidth(1);
    h->SetStats(0);
    h->SetTitle("");
  }
  hRatioFit->SetMarkerColor(kBlue);
  hRatioFit->SetLineColor(kBlue);
  hRatioCalib->SetMarkerColor(kRed);
  hRatioCalib->SetLineColor(kRed);

  hRatioFit->Draw("hist");
  TLine *line = new TLine(0,1,lFitRangeMax,1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray+1);
  line->Draw();
  hRatioFit->Draw("same hist");
  hRatioCalib->Draw("same hist");
  DrawAnchorLine(anchorPointFit, lLoRangeRatio, lHiRangeRatio, lFitRangeMax, kBlue);
  DrawAnchorLine(anchorPointCalib, lLoRangeRatio, lHiRangeRatio, lFitRangeMax, kRed);

  TH1F *hRatioFitZoomed = (TH1F*) hRatioFit->Clone("hRatioFitZoomed");
  TH1F *hRatioCalibZoomed = (TH1F*) hRatioCalib->Clone("hRatioCalibZoomed");
  hRatioFitZoomed->GetXaxis()->SetRangeUser(0,1150*lFitRangeMax/59000.);
  hRatioCalibZoomed->GetXaxis()->SetRangeUser(0,1150*lFitRangeMax/59000.);

  TPad *inset2 = new TPad("inset2", "",0.104677, 0.47032, 0.558775, 0.965753);
  inset2->SetPad(0.104677, 0.47032, 0.558775, 0.965753);
  inset2->SetTicks(1,1);
  inset2->SetBottomMargin(0.194);
  inset2->SetRightMargin(0.03);
  inset2->SetTopMargin(0.045);
  inset2->SetLeftMargin(0.11);
  inset2->cd();
  hRatioFitZoomed -> GetYaxis()->SetRangeUser(lLoRangeRatio, 1.55);
  hRatioFitZoomed -> GetYaxis()->SetTitleSize(0.1);
  hRatioFitZoomed -> GetXaxis()->SetTitleSize(0.1);
  hRatioFitZoomed -> GetYaxis()->SetTitleOffset(1.1);
  hRatioFitZoomed -> GetXaxis()->SetLabelSize(0.06);
  hRatioFitZoomed -> GetYaxis()->SetLabelSize(0.06);
  hRatioFitZoomed->Draw("hist");

  TLine *line2 = new TLine(0,1,1150*lFitRangeMax/59000.,1);
  line2->SetLineStyle(7);
  line2->SetLineColor(kGray+1);
  line2->Draw();
  hRatioFitZoomed->Draw("hist same");
  hRatioCalibZoomed->Draw("hist same");
  DrawAnchorLine(anchorPointFit, lLoRangeRatio, 1.55, 1150*lFitRangeMax/59000., kBlue);
  DrawAnchorLine(anchorPointCalib, lLoRangeRatio, 1.55, 1150*lFitRangeMax/59000., kRed);

  c1->cd(2);
  inset2->Draw();

  c1->SaveAs("hAnchorPointShift.pdf");

  // Second canvas: standalone version of the zoomed-in top pad
  const Double_t lZoomXMax = 590*lFitRangeMax/59000.;
  TCanvas *c2 = new TCanvas("c2", "", 900, 600);
  c2->SetFrameFillStyle(0);
  c2->SetFillStyle(0);
  c2->SetLogy();
  c2->SetTicks(1,1);
  c2->SetBottomMargin(0.12);
  c2->SetRightMargin(0.03);
  c2->SetTopMargin(0.03);
  c2->SetLeftMargin(0.10);

  TH1F *hV0MzoomedFull = (TH1F*) hV0M->Clone("hV0MzoomedFull");
  hV0MzoomedFull->GetXaxis()->SetRangeUser(0,lZoomXMax);
  hV0MzoomedFull->GetYaxis()->SetRangeUser(lYLoInset1, lYHiInset1);
  hV0MzoomedFull->GetXaxis()->SetTitleSize(0.05);
  hV0MzoomedFull->GetXaxis()->SetLabelSize(0.04);
  hV0MzoomedFull->GetYaxis()->SetTitleSize(0.05);
  hV0MzoomedFull->GetYaxis()->SetLabelSize(0.04);
  hV0MzoomedFull->GetYaxis()->SetTitleOffset(0.95);
  hV0MzoomedFull->Draw("hist");
  hFittedGlauber->Draw("same hist");
  hCalibGlauber->Draw("same hist");
  DrawAnchorLine(anchorPointFit, lYLoInset1, lYHiInset1 / 2e+2, lZoomXMax, kBlue);
  DrawAnchorLine(anchorPointCalib, lYLoInset1, lYHiInset1 / 2e+2, lZoomXMax, kRed);

  latPer->SetTextSize(0.045);
  latPer->DrawLatex(0.14, 0.89, "ALICE Pb#minusPb 5.36 TeV");

  TLegend *legZoom = (TLegend*) legO->Clone("legZoom");
  legZoom->SetX1NDC(0.5);
  legZoom->SetX2NDC(0.92);
  legZoom->SetY1NDC(0.62);
  legZoom->SetY2NDC(0.94);
  legZoom->SetTextSize(0.038);
  legZoom->Draw();

  c2->SaveAs("hAnchorPointShift_zoom.pdf");
}
