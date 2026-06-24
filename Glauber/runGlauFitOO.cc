#include "multCalibrator.h"
#include "multGlauberNBDFitter.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TTree.h"
#include "TStopwatch.h"

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

Double_t GetBoundaryForPercentile( TH1 *histo, Double_t lPercentileRequested ) {
  //This function returns the boundary for a specific percentile.
  Double_t lReturnValue = 0.0;
  Double_t lPercentile = 100.0 - lPercentileRequested;
  
  const Long_t lNBins = histo->GetNbinsX();
  Double_t lCountDesired = lPercentile * histo->Integral(1, histo->GetNbinsX())/100.;
  Long_t lCount = 0;
  for(Long_t ibin=1;ibin<lNBins;ibin++){
    lCount += histo->GetBinContent(ibin);
    if(lCount >= lCountDesired){
      //Found bin I am looking for!
      Double_t lWidth = histo->GetBinWidth(ibin);
      Double_t lLeftPercentile = 100.*(lCount - histo->GetBinContent(ibin))/histo->Integral(1, histo->GetNbinsX());
      Double_t lRightPercentile = 100.*lCount / histo->GetEntries();
      
      Double_t lProportion = (lPercentile - lLeftPercentile)/(lRightPercentile-lLeftPercentile);
      
      lReturnValue = histo->GetBinLowEdge(ibin) + lProportion*lWidth;
      break;
    }
  }
  return lReturnValue;
}

void runGlauFitOO(TString lInputFileName = "AnalysisResults-25ae.root", TString type = "hNPVContributors", TString histogramName = "hFT0C_BCs", TString nuclearProfile = "Oho2", float omega = 0.3){
  cout<<"Starting!"<<endl;
  TFile *file = new TFile(lInputFileName.Data(), "READ");
  file->ls();
  if(!file) cout<<"Problem with file!"<<endl;
  TH1F *hV0Mfine = 0x0;
  
  // regulate type
  // type options: hFT0C_Collisions hFT0M_Collisions hFV0A_Collisions hNPVContributors hNMFTTracks hNGlobalTracks
  TString XaxisLabel = "FT0C Amplitude";
  
  Double_t guessedMu = 0, guessedNorm = 0;
  int rebinFactor = 5;
  if(type == "hNPVContributors"){
    guessedMu = 1;
    guessedNorm = 1.60575e+08;
    XaxisLabel = "N PV contributors";
  }
  if(type == "hNGlobalTracks"){
    guessedMu = 5;
    guessedNorm = 1.60575e+08;
    XaxisLabel = "N global tracks";
  }
  if(type == "hNMFTTracks"){
    guessedMu = 5;
    guessedNorm = 1.60575e+08;
    XaxisLabel = "N MFT tracks";
  }
  // FIT detector series
  if(type == "hFT0C_Collisions"){
    rebinFactor = 20;
    guessedMu = 75;
    guessedNorm = 2.95339e+09;
  }
  if(type == "hFT0M_Collisions"){
    rebinFactor = 20;
    guessedMu = 150;
    guessedNorm = 2.95339e+09;
    XaxisLabel = "FT0M amplitude";
  }
  if(type == "hFT0M_BCs"){
    rebinFactor = 20;
    guessedMu = 150;
    guessedNorm = 2.95339e+09;
    XaxisLabel = "FT0M amplitude";
  }
  if(type == "hFT0MOuterA_Collisions"){
    rebinFactor = 20;
    guessedMu = 150;
    guessedNorm = 2.95339e+09;
    XaxisLabel = "FT0M amplitude";
  }
  if(type == "hFV0A_Collisions"){
    rebinFactor = 20;
    guessedMu = 500;
    guessedNorm = 2.95339e+09;
    XaxisLabel = "FV0A amplitude";
  }
  
  
  //hV0Mfine = (TH1F *) file -> Get(Form("centrality-study/Run_564359/%s", "hFT0C_Collisions"));
  // hV0Mfine = (TH1F *) file -> Get(Form("centrality-study/Run_%i/%s", runNumber, type.Data())); // hNPVContributors
  hV0Mfine = (TH1F *) file -> Get(Form("centrality-study/%s", histogramName.Data()));

  hV0Mfine->SetTitle("");
  
  // disregard bin zero
  cout<<"Received bin zero content: "<< hV0Mfine ->GetBinContent(0)<<", will set to zero..."<<endl;
  hV0Mfine ->SetBinContent(0, 0);
  
  cout<<"Input histogram has been received successfully! Information: "<<endl;
  
  cout<<"Counts: "<<hV0Mfine->GetEntries()<<endl;
  cout<<"NbinsX: "<<hV0Mfine->GetNbinsX()<<endl;
  cout<<"MaxX: "<<hV0Mfine->GetBinLowEdge(hV0Mfine->GetNbinsX()+1)<<endl;
  

  cout<<"Creating rebinned histogram with rebin factor: "<<rebinFactor<<endl;
  TH1F *hV0M = (TH1F*) hV0Mfine->Clone("hV0M");
  TH1F *hV0MUltraFine = (TH1F*) hV0Mfine->Clone("hV0MUltraFine");
  hV0M->SetName("hV0M");
  hV0M->SetTitle("");
  hV0M->Rebin(rebinFactor);

  TCanvas *c1 = new TCanvas("c1", "", 1300,900);
  c1->SetFrameFillStyle(0);
  c1->SetFillStyle(0);
  c1->Divide(1,2);
  c1->cd(1)->SetFrameFillStyle(0);
  c1->cd(1)->SetFillStyle(0);
  c1->cd(2)->SetFrameFillStyle(0);
  c1->cd(2)->SetFillStyle(0);
  
  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1,1);
  c1->cd(1)->SetPad(0,0.5,1,1);
  c1->cd(2)->SetPad(0,0.0,1,.5);
  
  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.25);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.07);
  
  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.25);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.07);
  c1->cd(2)->SetTicks(1,1);
  c1->cd(1);
  
  Double_t lFitRangeMax = GetBoundaryForPercentile(hV0Mfine, 0.008);
  cout<<"Fit range max estimated from histogram: "<<lFitRangeMax<<endl;
  
  float lFitRange = 0.012*GetBoundaryForPercentile(hV0Mfine, 0.01);
  if(lFitRangeMax<10000) lFitRange = 0.02*GetBoundaryForPercentile(hV0Mfine, 0.01);
  cout<<"Fit range low: "<<lFitRange<<endl;
  
  hV0M->GetXaxis()->SetRangeUser(0,lFitRangeMax);
  hV0M->GetYaxis()->SetRangeUser(0.25,hV0M->GetMaximum()*3);
  hV0M->SetLineColor(kBlack);
  hV0M->SetMarkerStyle(20);
  hV0M->SetMarkerColor(kBlack);
  hV0M->SetMarkerSize(0.5);
  hV0M->GetYaxis()->SetTitleSize(0.07);
  hV0M->GetYaxis()->SetLabelSize(0.05);
  hV0M->GetYaxis()->SetTitle("Count");
  hV0M->GetYaxis()->SetTitleOffset(0.5);
  hV0M->GetXaxis()->SetLabelSize(0.05);
  hV0M->GetXaxis()->SetTitleSize(0.06);
  hV0M->GetXaxis()->SetTitle("FT0A+C Amplitude");
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->SetStats(0);
  hV0M->Draw("E");
  
  
  //=======================================================
  
  TString baseHistoIdentifier = Form("basehistos_%s_%s_omega%.1f.root", nuclearProfile.Data(), nuclearProfile.Data(), omega);
  
  // TFile *fOutput = new TFile(Form("glauberfit_%i_%s_%s_%.1f.root", runNumber, type.Data(), nuclearProfile.Data(), omega), "RECREATE");
  
  TString lProcessedFileName = lInputFileName.Data();
  lProcessedFileName.ReplaceAll("ARs/", "results/");
  lProcessedFileName.ReplaceAll(".root", Form("_glauberNBD_%s.root", histogramName.Data()));
  TFile *fOutput = new TFile(lProcessedFileName.Data(), "RECREATE");

  
  //Stand back! Imma gonna do GLAUBER FITTIN'
  multGlauberNBDFitter *g = new multGlauberNBDFitter("lglau");
  g->SetAncestorMode(2);
  
  //Step 1: open the (Npart, Ncoll) pair information, provide
  TFile *fbasefile = new TFile(baseHistoIdentifier.Data(),"READ");
  TH2D *hNpNc = (TH2D*) fbasefile->Get("hNpNc");
  if(!hNpNc){
    cout<<"Problem with Npart Ncoll correlation!"<<endl;
    return;
  }
  
  // return to proper scope
  fOutput->cd();
  g->SetNpartNcollCorrelation(hNpNc);
  g->SetInputV0M(hV0M);
  g->SetFitRange(lFitRange,lFitRangeMax);
  
  TString lString = "REM0";
  g->SetFitOptions(lString.Data());
  g->SetFitNpx(100000);
  TF1 *fitfunc = g->GetGlauberNBD();
  


  cout<<"Guessed GlauberNBD mu value: "<<guessedMu<<endl;
  fitfunc->SetParameter(0,guessedMu); //mu value
  fitfunc->SetParLimits(0,0.25*guessedMu,2000);
  
  fitfunc->SetParameter(1,1.5); //k value
  fitfunc->FixParameter(2,0.8); // f value
  
  fitfunc->SetParameter(3, guessedNorm);
  
  //dMu/dNanc
  fitfunc->FixParameter(4,0);
  
  //dk/dNanc
  fitfunc->FixParameter(5,0);
  //fitfunc->SetParameter(5,1.63590e-03);
  
  //d2Mu/dNanc2
  fitfunc->FixParameter(6,0.0);
  //fitfunc->SetParameter(6,4.02271e-05);

  fitfunc->FixParameter(7,0.0);
  //fitfunc->SetParameter(7,-1.24349e-06);
  
  g->InitializeNpNc();
  g->InitAncestor();
  
  cout<<"WILL NOW ATTEMPT GLAUBER FIT"<<endl;
  cout<<"This will take a while. Please wait..."<<endl;
  Int_t lFitStatus = 0;
  lFitStatus = g->DoFit();
  Int_t lAttempts = 1;
  while(lAttempts < 10 && lFitStatus==0){
    cout<<"Attempting fit again ("<<lAttempts<<" attempt)..." <<endl;
    lFitStatus = g->DoFit();
  }
  cout<<"Final fit status: "<<lFitStatus<<endl;
  
  //Do a ratio plot
  TH1D *hGlauber = (TH1D*) hV0M->Clone("hGlauber");
  TH1D *hGlauberFine = (TH1D*) hV0MUltraFine->Clone("hGlauberFine");
  TH1D *hRatio = (TH1D*) hV0M->Clone("hRatio");
  hGlauber->Reset();
  
  c1->cd(1);

  
  cout<<"Calculating glauber function histogram with the same binning as data input... please wait..."<<endl;
  for(Int_t ii=1; ii<hGlauber->GetNbinsX()+1; ii++){
    Double_t lFuncVal = FastIntegrate( fitfunc, hGlauber->GetBinLowEdge(ii), hGlauber->GetBinLowEdge(ii+1), 4);
    hGlauber->SetBinContent(ii, lFuncVal);
    if(ii%500==0){
      cout<<"At integration #"<<ii<<"/"<<hGlauber->GetNbinsX()+1<<"..."<<endl;
    }
  }
  for(Int_t ii=1; ii<hGlauberFine->GetNbinsX()+1; ii++){
    Double_t lFuncVal = FastIntegrate( fitfunc, hGlauberFine->GetBinLowEdge(ii), hGlauberFine->GetBinLowEdge(ii+1), 4);
    hGlauberFine->SetBinContent(ii, lFuncVal);
    if(ii%500==0){
      cout<<"Fine calculation ... At integration #"<<ii<<"/"<<hGlauberFine->GetNbinsX()+1<<"..."<<endl;
    }
  }
  cout<<"Glauber function evaluated. Should go quickly now."<<endl;
  
  hGlauber->SetLineColor(kRed);
  hGlauber->SetLineWidth(1);
  hGlauber->Draw("same");
  
  c1->cd(2);
  Float_t lLoRangeRatio = 0.35;
  Float_t lHiRangeRatio = 1.65;
  hRatio->Divide(hGlauber);
  //hRatio->Scale(static_cast<float>(rebinFactor));
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetXaxis()->SetTitle(XaxisLabel.Data());
  hRatio->GetYaxis()->SetTitleSize(0.055);
  hRatio->GetYaxis()->SetTitleOffset(0.7);
  hRatio->GetXaxis()->SetTitleSize(0.055);
  hRatio->GetYaxis()->SetLabelSize(0.045);
  hRatio->GetXaxis()->SetLabelSize(0.045);
  hRatio->GetYaxis()->SetRangeUser(lLoRangeRatio,lHiRangeRatio);
  hRatio->GetXaxis()->SetRangeUser(0,lFitRangeMax);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerColor(kGray+2);
  hRatio->SetLineColor(kGray+2);
  //hRatio->SetMarkerSize(1.0);
  hRatio->SetMarkerSize(.7);
  hRatio->SetStats(0);
  //hRatioWide->SetStats(0);
  
  hRatio->Draw("hist");
  
  // draw 90% and 80% lines
  
  TLine *line = new TLine(0,1,lFitRangeMax,1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray+1);
  line->Draw();
  
  double boundary90 = GetBoundaryForPercentile(hGlauber, 90);
  double boundary80 = GetBoundaryForPercentile(hGlauber, 80);
  double boundary70 = GetBoundaryForPercentile(hGlauber, 70);
  double boundary60 = GetBoundaryForPercentile(hGlauber, 60);
  
  cout<<"boundary 90 = "<<boundary90<<endl;
  cout<<"boundary 80 = "<<boundary80<<endl;
  
  TLine *line90 = new TLine(boundary90,0.85,boundary90,1.15);
  TLine *line80 = new TLine(boundary80,0.85,boundary80,1.15);
  TLine *line70 = new TLine(boundary70,0.85,boundary70,1.15);
  TLine *line60 = new TLine(boundary60,0.85,boundary60,1.15);
  line90->SetLineStyle(7);
  line80->SetLineStyle(7);
  line70->SetLineStyle(7);
  line60->SetLineStyle(7);
  line90->SetLineColor(kGray+1);
  line80->SetLineColor(kGray+1);
  line70->SetLineColor(kGray+1);
  line60->SetLineColor(kGray+1);
  line90->Draw();
  line80->Draw();
  line70->Draw();
  line60->Draw();
  
  c1->cd(1);
  TLatex *lat = new TLatex();
  lat->SetNDC();
  Float_t lPosText = 0.76;
  Float_t lYShift = 0.25;
  lat->SetTextSize(0.042);

  // save the glauber parameters explicitly
  TH1D *hGlauberParameters = new TH1D("hGlauberParameters", "", 10,0,10);
  TH1D *hGlauberFitRange = new TH1D("hGlauberFitRange", "", 10,0,10);
  
  // Large-scale comparison of integrals
  cout<<"Ratio of integrals: "<<hV0MUltraFine->Integral(1,hV0MUltraFine->GetNbinsX())/hGlauber->Integral(1,hGlauber->GetNbinsX())<<endl;
  
  //fitfunc
  hGlauberParameters -> SetBinContent( 1, fitfunc -> GetParameter(0));
  hGlauberParameters -> SetBinContent( 2, fitfunc -> GetParameter(1));
  hGlauberParameters -> SetBinContent( 3, fitfunc -> GetParameter(2));
  hGlauberParameters -> SetBinContent( 4, fitfunc -> GetParameter(3));
  hGlauberParameters -> SetBinContent( 5, fitfunc -> GetParameter(4));
  hGlauberParameters -> Write();
  
  TH1D *hRatioOfIntegrals = new TH1D("hRatioOfIntegrals", "", 1, -0.5f, 0.5f);
  hRatioOfIntegrals->SetBinContent(1, hV0MUltraFine->Integral(1,hV0MUltraFine->GetNbinsX())/hGlauber->Integral(1,hGlauber->GetNbinsX()));
  hRatioOfIntegrals->Write();
  
  Double_t lLoRangeGlauber, lHiRangeGlauber;
  fitfunc->GetRange(lLoRangeGlauber, lHiRangeGlauber);
  hGlauberFitRange->SetBinContent(1, lLoRangeGlauber);
  hGlauberFitRange->SetBinContent(2, lHiRangeGlauber);
  hGlauberFitRange->Write();

  hV0MUltraFine->Write();
  hV0Mfine->Write();
  hV0M->Write();
  
  hGlauberFine->Write();
  hGlauber->Write();
  hRatio->Write();
  fOutput->Write();
  
  c1->SaveAs(Form("%s-fit.pdf", type.Data()));
  
}




