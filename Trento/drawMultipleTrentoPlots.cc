#include <iostream>
#include "TNamed.h"
// For Run Ranges functionality
#include <map>
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TStopwatch.h"
#include "TVirtualFitter.h"


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

Double_t GetCentralityBoundary(TH1 *hCalib, Double_t lCentrality){
  // hCalib maps multiplicity -> centrality (decreasing with multiplicity);
  // return the multiplicity at which the centrality crosses lCentrality
  for(Int_t ii=2; ii<hCalib->GetNbinsX()+1; ii++){
    if( hCalib->GetBinContent(ii-1) > lCentrality && hCalib->GetBinContent(ii) < lCentrality )
      return hCalib->GetBinLowEdge(ii);
  }
  return -1;
}

void drawMultipleTrentoPlots( TString lCalibFileNameNature  = "trentofit_New_calib.root",
                              TString lCalibFileNameJesper2 = "trentofit_Jesper2_calib.root") {
  TString lRunNumber = "572183";
  TString lXaxisTitle = "FT0C Amplitude";
  
  const Int_t nModels = 2;
  TString lCalibFileName[nModels] = {lCalibFileNameNature, lCalibFileNameJesper2};
  TString lModelName[nModels] = {"Nature", "Jesper2"};
  Int_t lColor[nModels] = {kBlue, kRed};
  
  int rebinRatio = 20;
  Double_t lFitRangeMax = 54000;
  
  TH1F *hV0M = 0x0;
  TH1F *hGlauberFunctionHyperFine[nModels];
  TH1F *hRatio[nModels];
  TH1F *hCalib[nModels];
  Double_t lBoundary[nModels][21];  // multiplicity at centrality = 0, 5, ..., 100%
  Double_t lMidRange[nModels][20];  // label position of the 5%-wide intervals
  
  for(Int_t im=0; im<nModels; im++){
    TFile *file = new TFile(lCalibFileName[im], "READ");
    file->ls();
    
    // data is identical in both files: take it from the first one
    if(im==0){
      hV0M = (TH1F*) file->Get("hData");
      if(!hV0M) {
        cout<<"hData not found"<<endl;
        hV0M = (TH1F*) file->Get("hV0MUltraFine");
      }
    }
    
    hGlauberFunctionHyperFine[im] = (TH1F*) file->Get("hGlauber");
    if (!hGlauberFunctionHyperFine[im]) {
      hGlauberFunctionHyperFine[im] = (TH1F*) file->Get("hGlauberFine");
    }
    hGlauberFunctionHyperFine[im]->SetName(Form("hGlauber_bc_%s", lModelName[im].Data()));
    hGlauberFunctionHyperFine[im]->SetLineColor(lColor[im]);
    hGlauberFunctionHyperFine[im]->SetMarkerColor(lColor[im]);
    hGlauberFunctionHyperFine[im]->SetLineWidth(2);
    
    hRatio[im] = (TH1F*) hV0M->Clone(Form("hRatio_%s", lModelName[im].Data()));
    TH1F *hGlauberFunction = (TH1F*) hGlauberFunctionHyperFine[im]->Clone(Form("hGlauberFunction_%s", lModelName[im].Data()));
    hRatio[im]->Rebin(rebinRatio);
    hGlauberFunction->Rebin(rebinRatio);
    hRatio[im]->Divide(hGlauberFunction);
    
    hCalib[im] = (TH1F*) file->Get("hCalib");
    
    for(Int_t ii=0; ii<21; ii++){
      lBoundary[im][ii] = GetCentralityBoundary(hCalib[im], ii*5);
    }
    for(Int_t ii=0; ii<20; ii++){
      TH1D *hSlice = GetCentralityHisto( hV0M, hCalib[im], ii*5, ii*5+5, Form("hSlice_%s_%i", lModelName[im].Data(), ii) );
      lMidRange[im][ii] = GetMidRange(hSlice, lFitRangeMax);
    }
  }
  
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
  
  Double_t lYminTop = 0.25;
  hV0M->GetXaxis()->SetRangeUser(0,lFitRangeMax);
  hV0M->GetYaxis()->SetRangeUser(lYminTop,hV0M->GetMaximum()*20);
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
  hV0M->GetXaxis()->SetTitle("FT0C Amplitude");
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->SetStats(0);
  hV0M->SetTitle("");
  hV0M->Draw("hist");
  
  TLatex *latPer = new TLatex();
  latPer->SetNDC();
  latPer->SetTextFont(42);
  latPer->SetTextSize(0.04);
  
  TLatex *latPerAngle = new TLatex();
  latPerAngle->SetTextFont(42);
  latPerAngle->SetTextAlign(22);
  latPerAngle->SetTextSize(0.030);
  latPerAngle->SetTextAngle(90);
  
  // 5% boundaries: lines from the bottom of the frame up to the data
  TLine *lBoundaryLine = new TLine();
  lBoundaryLine->SetLineStyle(2);
  auto DrawBoundaries = [&](Double_t lYmin, Bool_t lUpToData, Double_t lYmax, Double_t lXmax){
    for(Int_t im=0; im<nModels; im++){
      lBoundaryLine->SetLineColor(lColor[im]);
      for(Int_t ii=1; ii<20; ii++){
        Double_t x = lBoundary[im][ii];
        if( x<0 || x>lXmax ) continue;
        Double_t y = lUpToData ? hV0M->GetBinContent(hV0M->FindBin(x)) : lYmax;
        lBoundaryLine->DrawLine(x, lYmin, x, y);
      }
    }
  };
  
  DrawBoundaries(lYminTop, kTRUE, 0, lFitRangeMax);
  
  Double_t lLabelYTop[nModels] = {1, 12};
  for(Int_t im=0; im<nModels; im++){
    latPerAngle->SetTextColor(lColor[im]);
    for(Int_t ii=0; ii<20; ii=ii+1){
      if( ii*5<60 )
        latPerAngle->DrawLatex( lMidRange[im][ii], lLabelYTop[im], Form("%i-%i%%", ii*5, ii*5+5));
    }
  }
  hV0M->Draw("same hist");
  for(Int_t im=0; im<nModels; im++) hGlauberFunctionHyperFine[im]->Draw("same hist");
  
  latPer->SetTextSize(0.06);
  latPer-> DrawLatex(0.12, 0.88, "ALICE Pb-Pb 5.36 TeV");
  latPer->SetTextSize(0.033);
  latPer -> DrawLatex(0.12,0.835, "5%-wide intervals shown");
  latPer->SetTextSize(0.062);
  latPer -> DrawLatex(0.12,0.62, Form("Run %s", lRunNumber.Data()));
  
  TLegend *legO = new TLegend(0.12, 0.68, 0.42, 0.821);
  legO->SetBorderSize(0);
  legO->SetFillColorAlpha(0, 0);
  legO->AddEntry(hGlauberFunctionHyperFine[1], "T_{R}ENTo params from Run 1 fit", "l");
  legO->AddEntry(hGlauberFunctionHyperFine[0], "T_{R}ENTo params from Run 2 fit", "l");
  legO->Draw();
  
  // create extra inset pad
  
  TPad *inset1 = new TPad("inset1", "",0.575724, 0.441648, 0.924, 0.961098);
  inset1->SetPad(0.575724, 0.441648, 0.924, 0.961098);
  inset1->cd();
  
  Double_t lXmaxInset1 = 590*lFitRangeMax/59000.;
  Double_t lYminInset1 = hGlauberFunctionHyperFine[0]->GetMaximum() * 0.01;
  TH1F *hV0Mzoomed = (TH1F*) hV0M->Clone("hV0Mzoomed");
  hV0Mzoomed -> GetXaxis()->SetRangeUser(0,lXmaxInset1);
  hV0Mzoomed -> GetXaxis()->SetNdivisions(404);
  hV0Mzoomed -> GetXaxis()->SetTitleSize(0.1);
  hV0Mzoomed -> GetYaxis()->UnZoom();
  hV0Mzoomed -> GetYaxis()->SetRangeUser(lYminInset1, hV0Mzoomed->GetMaximum() * 5.5);
  hV0Mzoomed -> GetYaxis()->SetTitleSize(0.1);
  hV0Mzoomed -> GetYaxis()->SetTitleOffset(1.1); 
  hV0Mzoomed -> GetXaxis()->SetLabelSize(0.06);
  hV0Mzoomed -> GetYaxis()->SetLabelSize(0.06);
  hV0Mzoomed->Draw("hist");
  
  DrawBoundaries(lYminInset1, kTRUE, 0, lXmaxInset1);
  
  latPerAngle->SetTextSize(0.037);
  latPerAngle->SetTextAngle(0);
  Double_t lLabelYInset1[nModels] = {0.018, 0.04};
  for(Int_t im=0; im<nModels; im++){
    latPerAngle->SetTextColor(lColor[im]);
    for(Int_t ii=0; ii<20; ii=ii+1){
      if( ii*5>77 )
        latPerAngle->DrawLatex( lMidRange[im][ii], hGlauberFunctionHyperFine[0]->GetMaximum() * lLabelYInset1[im], Form("%i-%i%%", ii*5, ii*5+5));
    }
  }
  latPerAngle->SetTextSize(0.033);
  latPerAngle->SetTextAngle(90);
  
  hV0Mzoomed->Draw("same hist");
  for(Int_t im=0; im<nModels; im++) hGlauberFunctionHyperFine[im]->Draw("same hist");
  
  c1->cd(1);
  
  inset1->SetLogy();
  inset1->SetTicks(1,1);
  inset1->SetBottomMargin(0.194);
  inset1->SetRightMargin(0.03);
  inset1->SetTopMargin(0.085);
  inset1->SetLeftMargin(0.11);
  
  inset1->Draw("same");
  
  
  c1->cd(2);
  Float_t lLoRangeRatio = 0.3;
  Float_t lHiRangeRatio = 2.495;
  
  for(Int_t im=0; im<nModels; im++){
    hRatio[im]->GetYaxis()->SetTitle("Data/Fit");
    hRatio[im]->GetXaxis()->SetTitle(lXaxisTitle.Data());
    hRatio[im]->GetYaxis()->SetTitleSize(0.055);
    hRatio[im]->GetYaxis()->SetTitleOffset(0.7);
    hRatio[im]->GetXaxis()->SetTitleSize(0.055);
    hRatio[im]->GetYaxis()->SetLabelSize(0.045);
    hRatio[im]->GetXaxis()->SetLabelSize(0.045);
    hRatio[im]->GetYaxis()->SetRangeUser(lLoRangeRatio,lHiRangeRatio);
    hRatio[im]->GetXaxis()->SetRangeUser(0,lFitRangeMax);
    hRatio[im]->SetMarkerStyle(20);
    hRatio[im]->SetMarkerColor(lColor[im]);
    hRatio[im]->SetLineColor(lColor[im]);
    hRatio[im]->SetMarkerSize(.5);
    hRatio[im]->SetLineWidth(1);
    hRatio[im]->SetStats(0);
    hRatio[im]->SetTitle("");
  }
  
  hRatio[0]->Draw("hist");
  
  TLine *line = new TLine(0,1,lFitRangeMax,1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray+1);
  line->Draw();
  
  DrawBoundaries(lLoRangeRatio, kFALSE, 0.78, lFitRangeMax);
  
  Double_t lLabelYRatio[nModels] = {0.43, 0.64};
  latPerAngle->SetTextSize(0.026);
  for(Int_t im=0; im<nModels; im++){
    latPerAngle->SetTextColor(lColor[im]);
    for(Int_t ii=0; ii<20; ii=ii+1){
      if( ii*5<60 )
        latPerAngle->DrawLatex( lMidRange[im][ii], lLabelYRatio[im], Form("%i-%i%%", ii*5, ii*5+5));
    }
  }
  latPerAngle->SetTextSize(0.033);
  
  for(Int_t im=0; im<nModels; im++) hRatio[im]->Draw("same hist");
  
  Double_t lXmaxInset2 = 1150*lFitRangeMax/59000.;
  TPad *inset2 = new TPad("inset2", "",0.104677, 0.47032, 0.558775, 0.965753);
  inset2->SetPad(0.104677, 0.47032, 0.558775, 0.965753);
  inset2->SetTicks(1,1);
  inset2->SetBottomMargin(0.194);
  inset2->SetRightMargin(0.03);
  inset2->SetTopMargin(0.045);
  inset2->SetLeftMargin(0.11);
  inset2->cd();
  
  TH1F *hRatioZoomed[nModels];
  for(Int_t im=0; im<nModels; im++){
    hRatioZoomed[im] = (TH1F*) hRatio[im]->Clone(Form("hRatioZoomed_%s", lModelName[im].Data()));
    hRatioZoomed[im] -> GetXaxis()->SetRangeUser(0,lXmaxInset2);
    hRatioZoomed[im] -> GetYaxis()->SetRangeUser(0.4, 1.55);
    hRatioZoomed[im] -> GetYaxis()->SetTitleSize(0.1);
    hRatioZoomed[im] -> GetXaxis()->SetTitleSize(0.1);
    hRatioZoomed[im] -> GetYaxis()->SetTitleOffset(1.1);
    hRatioZoomed[im] -> GetXaxis()->SetLabelSize(0.06);
    hRatioZoomed[im] -> GetYaxis()->SetLabelSize(0.06);
  }
  hRatioZoomed[0]->Draw("hist");
  
  TLine *line2 = new TLine(0,1,lXmaxInset2,1);
  line2->SetLineStyle(7);
  line2->SetLineColor(kGray+1);
  line2->Draw();
  
  DrawBoundaries(0.4, kFALSE, 0.85, lXmaxInset2);
  
  latPerAngle->SetTextSize(0.04);
  latPerAngle->SetTextAngle(45);
  Double_t lLabelYInset2[nModels] = {0.52, 0.72};
  for(Int_t im=0; im<nModels; im++){
    latPerAngle->SetTextColor(lColor[im]);
    for(Int_t ii=0; ii<20; ii=ii+1){
      if( ii*5>70 )
        latPerAngle->DrawLatex( lMidRange[im][ii]+10, lLabelYInset2[im], Form("%i-%i%%", ii*5, ii*5+5));
    }
  }
  latPerAngle->SetTextSize(0.033);
  latPerAngle->SetTextAngle(90);
  
  for(Int_t im=0; im<nModels; im++) hRatioZoomed[im]->Draw("hist same");
  
  c1->cd(2);
  inset2->Draw();
  
  c1->SaveAs("trentofit_Nature_vs_Jesper2.pdf");
}
