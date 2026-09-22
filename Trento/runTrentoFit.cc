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

#include "multCalibrator.h"
#include "multGlauberNBDFitter.h"

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

//________________________________________________________________
Double_t FastIntegrate(TF1 *f1, Double_t a, Double_t b, Int_t n = 5){
  //Do fast integration with N sampling points
  const Int_t nc = n;
  Double_t x[nc], y[nc];
  Double_t lWidth = (b-a)/((double)(n-1));
  for(Int_t ii=0; ii<n; ii++){
    x[ii] = a + ((double)(ii))*lWidth;
    y[ii] = f1->Eval( x[ii] );
  }
  //Now go via trapezoids, please (this probably has a name)
  Double_t lIntegral = 0;
  for(Int_t ii=0; ii<n-1; ii++){
    lIntegral += 0.5*lWidth*(y[ii]+y[ii+1]);
  }
  return lIntegral/(b-a);
}

void runTrentoFit(Double_t lFitRange = 132.5, Double_t lPlotXrange = 54000)
{
  // This macro does an NBD + Trento fit to ALICE data
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1);
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // Acquire data to start
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  TFile* file = new TFile("../AnalysisResults/LHC23_pass5/AR_544122.root", "READ");
  TH1D* hV0M = (TH1D*)file->Get("centrality-study/hFT0C_BCs");
  TH1D* hV0MUltraFine = (TH1D*)hV0M->Clone("hV0MUltraFine");
  hV0M->SetName("hData");
  hV0M->SetMinimum(2e-1);
  hV0M->Rebin(20);

  TCanvas* c1 = new TCanvas("c1", "", 1600, 1200);
  c1->Divide(1, 2);
  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1, 1);
  c1->cd(1)->SetPad(0, 0.5, 1, 1);
  c1->cd(2)->SetPad(0, 0.0, 1, .5);

  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.02);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.07);

  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.02);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.07);
  c1->cd(2)->SetTicks(1, 1);
  c1->cd(1);

  hV0M->GetXaxis()->SetRangeUser(0, lPlotXrange);
  hV0M->SetTitle("");
  hV0M->SetLineColor(kBlack);
  hV0M->SetMarkerSize(0.5);
  hV0M->GetYaxis()->SetTitleSize(0.07);
  hV0M->GetYaxis()->SetLabelSize(0.05);
  hV0M->GetYaxis()->SetTitle("Count");
  hV0M->GetYaxis()->SetTitleOffset(0.5);
  hV0M->GetXaxis()->SetLabelSize(0.05);
  hV0M->GetXaxis()->SetTitleSize(0.06);
  hV0M->GetXaxis()->SetTitle("V0M Amplitude");
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->Draw("hist E");

  // Get Entropy, please

  TFile* fileTrento = new TFile("outfileJesper.root", "READ");
  TH1D* hEntropy = (TH1D*)fileTrento->Get("hEntropy");

  gStyle->SetOptStat(0);
  TCanvas* c1ent = new TCanvas("c1ent", "", 600, 600);
  c1ent->SetTicks(1, 1);
  c1ent->SetTopMargin(0.02);
  c1ent->SetBottomMargin(0.15);
  c1ent->SetLeftMargin(0.16);
  c1ent->SetRightMargin(0.05);
  c1ent->SetLogy();

  hEntropy->GetXaxis()->SetTitle("Entropy");
  hEntropy->GetXaxis()->SetTitleOffset(1.35);
  hEntropy->GetYaxis()->SetTitle("Frequency");
  hEntropy->GetYaxis()->SetTitleSize(0.050);
  hEntropy->GetXaxis()->SetTitleSize(0.050);
  hEntropy->GetYaxis()->SetLabelSize(0.041);
  hEntropy->GetXaxis()->SetLabelSize(0.041);
  hEntropy->Draw("");

  c1->cd(1);

  multGlauberNBDFitter* g = new multGlauberNBDFitter("fitter");
  g->InitTrentoNBD();
  g->SetNSources(hEntropy);
  g->SetInputV0M(hV0M);
  
  Double_t lFitRangeMax = GetBoundaryForPercentile(hV0M, 0.008);
  lFitRange = 0.012*GetBoundaryForPercentile(hV0M, 0.01);
  
  g->SetFitRange(lFitRange, lFitRangeMax);
  cout<<"Fit range: ["<<lFitRange <<", "<<lFitRangeMax<<"]"<<endl;
  // Step 3: go for it ...

  TString lString = "REM0";
  g->SetFitOptions(lString.Data());
  g->SetFitNpx(100000);
  TF1* fitfunc = g->GetTrentoNBD();

  Double_t guessedMu = lFitRangeMax/53968.4 * 0.175*3.53971e+02;
  Double_t guessedK = 1.5f;
  Double_t guessedF = 0.800f;
  // Double_t guessedF = 0;
  Double_t guessedDmu = 0;

  fitfunc->SetParameter(0, guessedMu); // mu value
  fitfunc->SetParameter(1, guessedK);  // k value
  fitfunc->FixParameter(2, guessedF);  // k value
  fitfunc->FixParameter(4, guessedDmu);  // k value

  g->DoFit();
  hV0M->GetXaxis()->SetRangeUser(0, lPlotXrange);

  for (Int_t ii = 0; ii < 50; ii++) {
    cout << "Sample function at " << ii * 1000 << ": " << fitfunc->Eval(ii * 1000) << endl;
  }

  // Do a ratio plot
  TH1D* hGlauber = (TH1D*)hV0MUltraFine->Clone("hGlauber");
  TH1D* hRatio = (TH1D*)hV0MUltraFine->Clone("hRatio");
  hGlauber->Reset();

  cout<<"Calculating glauber function histogram with the same binning as data input... please wait..."<<endl;
  for(Int_t ii=1; ii<hGlauber->GetNbinsX()+1; ii++){
    Double_t lFuncVal = FastIntegrate( fitfunc, hGlauber->GetBinLowEdge(ii), hGlauber->GetBinLowEdge(ii+1), 4);
    hGlauber->SetBinContent(ii, lFuncVal);
    if(ii%500==0){
      cout<<"At integration #"<<ii<<"/"<<hGlauber->GetNbinsX()+1<<"..."<<endl;
    }
  }
  cout<<"Glauber function evaluated. Should go quickly now."<<endl;
  hGlauber->SetLineColor(kRed);
  hGlauber->Draw("same");
  hRatio->Divide(hGlauber);

  c1->cd(2);
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetXaxis()->SetTitle("FT0C Amplitude");
  hRatio->GetYaxis()->SetRangeUser(0.45, 1.55);
  hRatio->SetLineColor(kRed);
  hRatio->SetMarkerSize(0.5);

  hRatio->Draw("hist");
  TLine* line = new TLine(0, 1, lPlotXrange, 1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray + 1);
  line->Draw();
  hRatio->Draw("same hist");

  c1->cd(1);
  TLatex* lat = new TLatex();
  lat->SetNDC();
  lat->DrawLatex(0.1, 0.53 - 0.11, "Pb-Pb 5.36 TeV Trento + NBD fit (2026)");
  lat->SetTextFont(42);
  lat->DrawLatex(0.1, 0.47 - 0.11, Form("Fit range: [%.1f, %.1f]", lFitRange, lFitRangeMax));
  lat->DrawLatex(0.1, 0.41 - 0.11, Form("#Chi^{2}/ndf: %.1f / %i = %.3f", fitfunc->GetChisquare(), fitfunc->GetNDF(), fitfunc->GetChisquare() / ((Double_t)(fitfunc->GetNDF()))));
  lat->DrawLatex(0.1, 0.35 - 0.11, Form("Fit options: %s", lString.Data()));
  lat->DrawLatex(0.1, 0.29 - 0.11, Form("NBD #mu: %.3f", fitfunc->GetParameter(0)));
  lat->DrawLatex(0.1, 0.23 - 0.11, Form("NBD k: %.3f", fitfunc->GetParameter(1)));

  // return;

  // Now extract hyper-fine cumulative function from fit function
  const Long_t lSamplePoints = 1e+5; // because, because. Just because.
  Double_t lDelta = 50000. / ((Double_t)(lSamplePoints));
  Double_t lX[lSamplePoints], lY[lSamplePoints];

  lX[0] = 0;
  lY[0] = 0;
  for (Long_t ii = 1; ii < lSamplePoints; ii++) {
    if (ii % 50000 == 0) {
      cout << "At sample #" << ii << endl;
    }
    lX[ii] = ((Double_t)ii) * lDelta;
    lY[ii] = lY[ii - 1] + fitfunc->Eval(lX[ii]);
  }
  Int_t lFirstPointAbove = -1;
  for (Long_t ii = 1; ii < lSamplePoints; ii++) {
    lX[ii] = ((Double_t)ii) * lDelta;
    lY[ii] = lY[ii] / lY[lSamplePoints - 1]; // Normalize
    if (lFirstPointAbove < 0 && lY[ii] > 0.1) {
      lFirstPointAbove = ii;
    }
  }

  TGraph* gr = new TGraph(lSamplePoints, lX, lY);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.3);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);

  TCanvas* cCumu = new TCanvas("cCumu", "", 800, 600);
  cCumu->SetTicks(1, 1);
  gr->Draw("ALP");

  // Determine anchor point based on hyperfinely binned TGraph, please
  Double_t lAnchorPoint;
  Double_t lFrac = (0.1 - gr->GetY()[lFirstPointAbove - 1]) / (gr->GetY()[lFirstPointAbove] - gr->GetY()[lFirstPointAbove - 1]);
  lAnchorPoint = gr->GetX()[lFirstPointAbove - 1] + lFrac * (gr->GetX()[lFirstPointAbove] - gr->GetX()[lFirstPointAbove - 1]);
  cout << "Anchor point determined to be: " << lAnchorPoint << endl;
  Double_t lFracAnchoredOut = gr->Eval(lFitRange);

  c1->cd(1);
  lat->DrawLatex(0.1, 0.11 - 0.05, Form("90%% anchor point: %.3f (percentage above fit lower limit: %.1f%%)", lAnchorPoint, 100. * (1 - lFracAnchoredOut)));

  c1->SaveAs("TrentoNBD.pdf");

  TFile* outfile = new TFile("trentofit.root", "recreate");
  hV0M->Write();
  hV0MUltraFine->Write();
  hGlauber->Write();
  hRatio->Write();

  // save the glauber parameters explicitly
  TH1D *hGlauberParameters = new TH1D("hGlauberParameters", "", 10,0,10);
  TH1D *hGlauberFitRange = new TH1D("hGlauberFitRange", "", 10,0,10);
  
  //fitfunc
  hGlauberParameters -> SetBinContent( 1, fitfunc -> GetParameter(0));
  hGlauberParameters -> SetBinContent( 2, fitfunc -> GetParameter(1));
  hGlauberParameters -> SetBinContent( 3, fitfunc -> GetParameter(2));
  hGlauberParameters -> SetBinContent( 4, fitfunc -> GetParameter(3));
  hGlauberParameters -> SetBinContent( 5, fitfunc -> GetParameter(4));
  hGlauberParameters -> Write();
  
  Double_t lLoRangeGlauber, lHiRangeGlauber;
  fitfunc->GetRange(lLoRangeGlauber, lHiRangeGlauber);
  hGlauberFitRange->SetBinContent(1, lLoRangeGlauber); 
  hGlauberFitRange->SetBinContent(2, lHiRangeGlauber);
  hGlauberFitRange->Write();
  outfile->Write();
  outfile->Close();
}