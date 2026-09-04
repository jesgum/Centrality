#include "TNamed.h"
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

#include <iostream>

using namespace std;
class NBDFitter : public TNamed {
public:
  //basic functionality
  NBDFitter(const char * name, const char * title = "Glauber+NBD fitter");
  ~NBDFitter();
  
  //Master fitter function
  Double_t ProbDistrib(Double_t *x, Double_t *par);
  
  //Do Fit: where everything happens
  Bool_t DoFit();
  
  //Set input characteristics: the 2D plot with Npart, Nanc
  Bool_t SetNsources(TH1 *hNsources);
  
  //Set main input to be fitted (the V0M distribution)
  Bool_t SetInputV0M(TH1 *hV0M);
  
  //Interface to get funtions if asked to
  TF1 *GetNBD();
  TF1 *GetGlauberNBD();
  
  //Helper
  Bool_t InitializeNpNc();
  
  //Interface to set vals
  void SetMu ( Double_t lVal ) {fMu = lVal;}
  void Setk ( Double_t lVal ) {fk = lVal;}
  void SetNorm ( Double_t lVal ) {fnorm = lVal;}
  
  //Interface to get vals
  Double_t GetMu ()   {return fMu;}
  Double_t Getk ()    {return fk;}
  Double_t GetNorm () {return fnorm;}
  
  void SetFitRange  (Double_t lMin, Double_t lMax);
  void SetFitOptions(TString lOpt);
  
  //void    Print(Option_t *option="") const;
  
private:
  //This function serves as the (analytical) NBD
  TF1 *fNBD;
  
  //This function is the key fitting function
  TF1 *fMasterFit;
  
  //Reference histo
  TH1 *fhNsources; //basic ancestor distribution
  TH1 *fhV0M; //basic ancestor distribution
  
  //The actual output: mu, k, norm
  Double_t fMu;
  Double_t fk;
  Double_t fnorm;
  
  TString fFitOptions;
  
  ClassDef(NBDFitter, 1);
};
//________________________________________________________________
NBDFitter::NBDFitter(const char * name, const char * title): TNamed(name,title),
fNBD(0x0),
fMasterFit(0x0),
fhNsources(0x0),
fMu(45),
fk(1.5),
fnorm(100),
fFitOptions("R0")
{
  //NBD
  fNBD = new TF1("fNBD","ROOT::Math::negative_binomial_pdf(x,[0],[1])",0,800);
  
  //master function
  fMasterFit = new TF1("fMasterFit", this, &NBDFitter::ProbDistrib,
                       0, 50000, 3 , "NBDFitter", "ProbDistrib");
  fMasterFit->SetParameter(0,fMu);
  fMasterFit->SetParameter(1,fk);
  fMasterFit->SetParameter(2,fnorm);
}
//________________________________________________________________
NBDFitter::~NBDFitter() {
  // Destructor
  if (fNBD) {
    delete fNBD;
    fNBD = 0x0;
  }
  if (fMasterFit) {
    delete fMasterFit;
    fMasterFit = 0x0;
  }
  if (fhNsources) {
    delete fhNsources;
    fhNsources = 0x0;
  }
}
//______________________________________________________
Double_t NBDFitter::ProbDistrib(Double_t *x, Double_t *par)
//Master fitter function
{
  //Double_t lMultValue = TMath::Floor(x[0]+0.5);
  Double_t lMultValue = x[0];
  Double_t lProbability = 0.0;
  //______________________________________________________
  //Actually ealuate function
  for(Long_t iNanc = 1; iNanc<fhNsources->GetNbinsX()+1; iNanc++){
    Double_t lNsources = fhNsources->GetBinCenter(iNanc);
    Double_t lThisMu = ((Double_t)lNsources)*par[0];
    Double_t lThisk = ((Double_t)lNsources)*par[1];
    Double_t lpval = TMath::Power(1+lThisMu/lThisk,-1);
    fNBD->SetParameter(1,lThisk);
    fNBD->SetParameter(0,lpval);
    Double_t lMult = fNBD->Eval(lMultValue);
    lProbability += fhNsources->GetBinContent(fhNsources->FindBin(iNanc))*lMult;
  }
  //______________________________________________________
  return par[2]*lProbability;
}
//________________________________________________________________
Bool_t NBDFitter::SetNsources(TH1 *hNsources){
  Bool_t lReturnValue = kTRUE;
  if( hNsources ){
    fhNsources = (TH1*) hNsources;
  }else{
    lReturnValue = kFALSE;
  }
  return lReturnValue;
}
//________________________________________________________________
Bool_t NBDFitter::SetInputV0M(TH1 *hV0M){
  Bool_t lReturnValue = kTRUE;
  if( hV0M ){
    fhV0M = (TH1*) hV0M;
  }else{
    lReturnValue = kFALSE;
  }
  return lReturnValue;
}
//________________________________________________________________
TF1 *NBDFitter::GetNBD(){
  return fNBD;
}
//________________________________________________________________
TF1 *NBDFitter::GetGlauberNBD(){
  return fMasterFit;
}
//________________________________________________________________
void NBDFitter::SetFitRange( Double_t lMin, Double_t lMax){
  fMasterFit -> SetRange(lMin, lMax);
}
//________________________________________________________________
void NBDFitter::SetFitOptions(TString lOpt){
  fFitOptions = lOpt;
}
//________________________________________________________________
Bool_t NBDFitter::DoFit(){
  //Try very hard, please
  TVirtualFitter::SetMaxIterations(5000000);
  
  TStopwatch* timer = new TStopwatch();
  timer->Start ( kTRUE );
  cout<<"---> Now fitting, please wait..."<<endl;
  
  fMasterFit->SetNpx(1000);
  fhV0M->Fit("fMasterFit",fFitOptions.Data());
  
  timer->Stop();
  Double_t lTotalTime = timer->RealTime();
  cout<<"---> Fitting took "<<lTotalTime<<" seconds"<<endl;
  
  fMu   = fMasterFit -> GetParameter(0);
  fk    = fMasterFit -> GetParameter(1);
  fnorm = fMasterFit -> GetParameter(2);
  
  return kTRUE;
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
//________________________________________________________________

void doTrentoFit(Double_t lFitRange = 132.5, Double_t lPlotXrange = 45000){
  //This macro does an NBD + Trento fit to ALICE data
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1);
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  // Acquire data to start
  //+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
  //TFile *file = new TFile("~/Downloads/calibrator/LEGO/LHC15o/AnalysisResults.root", "READ");
  //TFile *file = new TFile("gagliardi/AnalysisResultsAll.root", "READ");
  TFile *file = new TFile("AnalysisResults_244918.root", "READ");
  
  TTree *fTree = (TTree*) file->Get("MultSelection/fTreeEvent");
  fTree->Print();
  
  TProfile *fProVtxZ = new TProfile("fProVtxZ","", 10,-10,10);
  fTree->Draw("(fAmplitude_V0A+fAmplitude_V0C):fEvSel_VtxZ>>fProVtxZ", "fnContributors>1&&fEvSel_Triggered&&TMath::Abs(fEvSel_VtxZ)<10", "goff");
  
  TCanvas *c0 = new TCanvas("c0", "", 800,600);
  fProVtxZ->Draw();
  
  TF1 *f1 = new TF1("f1", "[0]*(1+[1]*x+[2]*x*x+[3]*x*x*x)", -10,10);
  fProVtxZ->Fit("f1", "REM0");
  f1->Draw("same");
  
  TH1D *hV0M = new TH1D("hV0M", "", 1200, 0, 60000);
  TString lExpression = Form("(fAmplitude_V0A+fAmplitude_V0C)/(1+(%.10f)*fEvSel_VtxZ+(%.10f)*fEvSel_VtxZ*fEvSel_VtxZ+(%.10f)*fEvSel_VtxZ*fEvSel_VtxZ*fEvSel_VtxZ)",
                             f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3));
  //lExpression = "(fAmplitude_V0A+fAmplitude_V0C)";
  fTree->Draw(Form("%s>>hV0M",lExpression.Data()), "fnContributors>1&&fEvSel_Triggered&&TMath::Abs(fEvSel_VtxZ)<10", "goff");
  
  TCanvas *c1 = new TCanvas("c1", "", 1600,900);
  c1->Divide(1,2);
  c1->cd(1);
  c1->cd(1)->SetLogy();
  c1->cd(1)->SetTicks(1,1);
  c1->cd(1)->SetPad(0,0.5,1,1);
  c1->cd(2)->SetPad(0,0.0,1,.5);
  
  c1->cd(1)->SetBottomMargin(0.001);
  c1->cd(1)->SetRightMargin(0.02);
  c1->cd(1)->SetTopMargin(0.02);
  c1->cd(1)->SetLeftMargin(0.07);
  
  c1->cd(2)->SetBottomMargin(0.14);
  c1->cd(2)->SetRightMargin(0.02);
  c1->cd(2)->SetTopMargin(0.001);
  c1->cd(2)->SetLeftMargin(0.07);
  c1->cd(2)->SetTicks(1,1);
  c1->cd(1);
  
  hV0M->GetXaxis()->SetRangeUser(0,39000);
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
  hV0M->GetXaxis()->SetTitle("V0M Amplitude");
  hV0M->GetYaxis()->SetTickLength(0.015);
  hV0M->Draw("E");
  
  //Get Entropy, please
  
  TFile *fileTrento = new TFile("Trento_PbPb_5p02_converted.root", "READ");
  TTree *fTreeTrento = (TTree*) fileTrento->Get("fTree");
  
  TH1D *hEntropy = new TH1D("hEntropy", "", 700,0,7e+2);
  fTreeTrento->Draw("fEntropy/100.>>hEntropy", "", "goff");
  
  gStyle->SetOptStat(0);
  TCanvas *c1ent = new TCanvas("c1ent", "", 600,600);
  c1ent->SetTicks(1,1);
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
  hEntropy->Draw("") ;
  
  c1->cd(1);
  NBDFitter *g = new NBDFitter("fBazooka");
  g->SetNsources( hEntropy );
  g->SetInputV0M( hV0M );
  
  Double_t lFitRangeMax = 43000;
  
  g->SetFitRange(lFitRange,lFitRangeMax);
  //Step 3: go for it ...
  
  
  TString lString = "REM0";
  g->SetFitOptions(lString.Data());
  
  TF1 *fitfunc = g->GetGlauberNBD();
  
  fitfunc->SetParameter(0,50);
  fitfunc->SetParameter(1,30);
  //fitfunc->SetParameter(2,1);
  //
  //fitfunc->FixParameter(2,0.8);
  
  g->DoFit();
  
  hV0M->GetXaxis()->SetRangeUser(0,lPlotXrange);
  fitfunc->Draw("same");
  
  
  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  fitfunc->SetLineColor(kRed);
  fitfunc->Draw("");
  
  for (Int_t ii=0; ii<50; ii++){
    cout<<"Sample function at "<<ii*1000<<": "<<fitfunc->Eval(ii*1000)<<endl;
  }
  
  
  //Do a ratio plot
  TH1D *hRatio = (TH1D*) hV0M->Clone("hRatio");
  TH1D *hRatioWide = (TH1D*) hV0M->Clone("hRatio");
  hRatioWide->Rebin(20);
  hRatioWide->Scale(1./20.);
  //hRatio->Scale(1.,"width") ;
  
  for(Int_t ii=1; ii<hRatio->GetNbinsX()+1; ii++){
    Double_t lRatio = hRatio->GetBinContent(ii);
    Double_t lFuncVal = fitfunc->Eval( hRatio->GetBinCenter(ii) );
    
    if ( lRatio < 2 ){
      hRatio->SetBinContent(ii,-1);
      hRatio->SetBinError(ii,1e-9);
      continue;
    }
    lRatio /= lFuncVal;
    hRatio->SetBinContent(ii, lRatio);
    hRatio->SetBinError(ii, hRatio->GetBinError(ii)/lFuncVal);
  }
  
  for(Int_t ii=1; ii<hRatioWide->GetNbinsX()+1; ii++){
    Double_t lRatio = hRatioWide->GetBinContent(ii);
    Double_t lFuncValCenter = fitfunc->Eval( hRatioWide->GetBinCenter(ii) );
    Double_t lFuncVal = FastIntegrate( fitfunc, hRatioWide->GetBinLowEdge(ii), hRatioWide->GetBinLowEdge(ii+1), 100 );
    cout<<"center: "<<lFuncValCenter<<" fastintegrate: "<<lFuncVal<<endl;
    if ( lRatio < 2 ){
      hRatioWide->SetBinContent(ii,-1);
      hRatioWide->SetBinError(ii,1e-9);
      continue;
    }
    lRatio /= lFuncVal;
    hRatioWide->SetBinContent(ii, lRatio);
    hRatioWide->SetBinError(ii, hRatioWide->GetBinError(ii)/lFuncVal);
  }
  
  c1->cd(2);
  hRatio->GetYaxis()->SetTitle("Data/Fit");
  hRatio->GetYaxis()->SetRangeUser(0.45,1.55);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerSize(0.5);
  
  hRatio->Draw("hist");
  hRatioWide->SetLineColor(kBlue);
  hRatioWide->SetMarkerColor(kBlue);
  hRatioWide->SetMarkerStyle(20);
  hRatioWide->SetMarkerSize(0.7);
  hRatioWide->SetLineWidth(2);
  
  TLine *line = new TLine(0,1,lPlotXrange,1);
  line->SetLineStyle(7);
  line->SetLineColor(kGray+1);
  line->Draw();
  
  TLine *lFitRangeLine = new TLine( lFitRange, 0.5, lFitRange, 1.5) ;
  lFitRangeLine->SetLineColor(kRed);
  lFitRangeLine->Draw();
  
  TH1D *hRatioGrayed = (TH1D*) hRatio->Clone("hRatioGrayed");
  hRatioGrayed->SetMarkerColor(kGray+1);
  hRatioGrayed->SetLineColor(kGray+1);
  hRatioGrayed->Draw("same");
  
  hRatio->SetLineWidth(1);
  hRatio->Draw("same hist");
  hRatioWide->Draw("same");
  
  c1->cd(1);
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->DrawLatex(0.1,0.53-0.11, "Pb-Pb 5.02 TeV Trento + NBD fit (2015)");
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.47-0.11, Form("Fit range: %.1f-%.1f", lFitRange, lFitRangeMax) );
  lat->DrawLatex(0.1,0.41-0.11, Form("#Chi^{2}/ndf: %.1f / %i = %.3f", fitfunc->GetChisquare(), fitfunc->GetNDF(), fitfunc->GetChisquare() / ((Double_t)(fitfunc->GetNDF() ) ) ) );
  lat->DrawLatex(0.1,0.35-0.11, Form("Fit options: %s", lString.Data() ) );
  lat->DrawLatex(0.1,0.29-0.11, Form("NBD #mu: %.3f", fitfunc->GetParameter(0) ) );
  lat->DrawLatex(0.1,0.23-0.11, Form("NBD k: %.3f", fitfunc->GetParameter(1) ) );
  
  
  //return;
  
  //Now extract hyper-fine cumulative function from fit function
  const Long_t lSamplePoints = 1e+5; //because, because. Just because.
  
  Double_t lDelta = 50000./((Double_t)(lSamplePoints));
  
  Double_t lX[lSamplePoints], lY[lSamplePoints];
  
  lX[0] = 0; lY[0] = 0;
  for(Long_t ii=1; ii<lSamplePoints; ii++){
      if(ii%50000==0) cout<<"At sample #"<<ii<<endl;
      lX[ii] = ((Double_t) ii) * lDelta ;
      lY[ii] = lY[ii-1] + fitfunc->Eval(lX[ii]);
  }
  Int_t lFirstPointAbove = -1;
  for(Long_t ii=1; ii<lSamplePoints; ii++){
      lX[ii] = ((Double_t) ii) * lDelta ;
      lY[ii] = lY[ii] / lY[lSamplePoints-1] ; //Normalize
      if( lFirstPointAbove<0 && lY[ii]>0.1 ) lFirstPointAbove=ii;
  }
  
  TGraph *gr = new TGraph(lSamplePoints, lX, lY);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.3);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);
  
  TCanvas *cCumu = new TCanvas("cCumu", "", 800,600);
  cCumu->SetTicks(1,1);
  gr->Draw("ALP") ;
  
  //Determine anchor point based on hyperfinely binned TGraph, please
  Double_t lAnchorPoint;
  Double_t lFrac = (0.1 - gr->GetY()[lFirstPointAbove-1])/(gr->GetY()[lFirstPointAbove] - gr->GetY()[lFirstPointAbove-1]);
  lAnchorPoint = gr->GetX()[lFirstPointAbove-1] + lFrac*(gr->GetX()[lFirstPointAbove] - gr->GetX()[lFirstPointAbove-1]);
  cout<<"Anchor point determined to be: "<<lAnchorPoint<<endl;
  
  /*
   TH1D *hCumulative = new TH1D("hCumulative", "", lSamplePoints, 0, 50000);
  cout<<"Calculating cumulative..."<<endl;
  hCumulative->SetBinContent(1,0);
  for(Long_t ii=1; ii<lSamplePoints; ii++){
      if(ii%50000==0) cout<<"At sample #"<<ii<<endl;
      hCumulative->SetBinContent(ii+1, hCumulative->GetBinContent(ii) +
      lDelta*0.5*(fitfunc->Eval(hCumulative->GetBinLowEdge(ii)) ) );
  }
  
  TCanvas *cCumu = new TCanvas("cCumu", "", 800,600);
  cCumu->SetTicks(1,1);
  hCumulative->Scale(1./hCumulative->GetBinContent( lSamplePoints-1 ));
  hCumulative->Draw();

  Double_t lAnchorPoint = 0;
  //Find precise 90% anchor point
  for(Long_t ii=1; ii<lSamplePoints; ii++){
      if( hCumulative->GetBinContent(ii) > 0.1){
          cout<<"Located anchor point for 90% anchoring!"<<endl;
          cout<<"It's at approximately: "<<hCumulative->GetBinCenter(ii)<<endl;
          cout<<"Will now smoothen it slightly for you, hang on"<<endl;
          Double_t lVal2 = hCumulative->GetBinLowEdge(ii);
          Double_t lVal1 = hCumulative->GetBinLowEdge(ii-1);
          
          //Proportion:
          Double_t lFrac = (10.-hCumulative->GetBinContent(ii-1))/
          (hCumulative->GetBinContent(11)-hCumulative->GetBinContent(ii-1));
          lAnchorPoint = lVal1 + lFrac*(lVal2-lVal1);
          cout<<"Better value: "<<lAnchorPoint<<endl;
          
          break;
      }
  }
   */
  
  Double_t lFracAnchoredOut = gr->Eval(lFitRange);
  /* DEPRECATED
  for(Long_t ii=1; ii<lSamplePoints; ii++){
      if( hCumulative->GetBinCenter(ii) > lFitRange){
          cout<<"Located equivalent percentage of hadronic cross section"<<endl;
          //Proportion:
          Double_t lVal2 = hCumulative->GetBinContent(ii);
          Double_t lVal1 = hCumulative->GetBinContent(ii-1);
          Double_t lFrac = (lFitRange-hCumulative->GetBinCenter(ii-1))/
          (hCumulative->GetBinCenter(11)-hCumulative->GetBinCenter(ii-1));
          lFracAnchoredOut = lVal1 + lFrac*(lVal2-lVal1);
          cout<<"Fraction Anchored out: "<<lFracAnchoredOut<<endl;
          cout<<"Anchor percentile for fit: "<<1-lFracAnchoredOut<<endl;
          break;
      }
  }
   */
  
  c1->cd(1);
  lat->DrawLatex(0.1,0.11-0.05, Form("90%% anchor point: %.3f (percentage above fit lower limit: %.1f%%)", lAnchorPoint, 100.*(1-lFracAnchoredOut) ) );
  
  
  c1->SaveAs("TrentoNBD.pdf");
} 
