// void runCalibration(TString lInputFileName = "results/AR_544122_glauberNBD_ancestorMode2_hFT0C_BCs.root")
void runCalibrationOO(TString lInputFileName = "results/AR_544122_glauberNBD_ancestorMode2_hFT0C_BCs.root", Bool_t doNpartNcoll = kFALSE) {
  TString nuclearProfile = "Oho2";
  // TString nuclearProfile = "Ne";
  float omega = 0.3;
  TFile *file = new TFile(lInputFileName.Data(), "READ");
  file->ls();
  
  TH1F *hData = (TH1F*) file->Get("hV0MUltraFine");
  
  TH1F *hGlauberParameters = (TH1F*) file->Get("hGlauberParameters");
  if(!hGlauberParameters) cout<<"No hGlauberParameters"<<endl;
  TH1F *hGlauberFitRange = (TH1F*) file->Get("hGlauberFitRange");
  if(!hGlauberFitRange) cout<<"No hGlauberFitRange"<<endl;
  hData->SetName("hData");
  TH1F *hStitched = (TH1F*) hData->Clone("hStitched");
  TH1F *hFit = (TH1F*) file->Get("hGlauberFine");
  if(!hFit) cout<<"No hFit"<<endl;
  
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->SetLeftMargin(0.17);
  c1->SetBottomMargin(0.17);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.05);
  c1->SetTicks(1,1);
  c1->SetLogz();
  c1->SetFrameFillStyle(0);
  c1->SetFillStyle(0);
  
  /*
    double matchRange = 50; // first guess
    cout<<"Data bin width: "<<hData->GetBinWidth(1)<<endl;
    cout<<"Fit bin width: "<<hFit->GetBinWidth(1)<<endl;
    cout<<"Match range to use: "<<matchRange<<endl;
    
    double anchorPoint = -1;
    
    // doing partial integration up to certain point for anchor point
    for(int ii=1; ii<hData->GetNbinsX()+1; ii++){
      // renormalize data curve
      int bin1 = ii+1;
      int bin2 = hData->FindBin( hData->GetBinLowEdge(ii+1) + matchRange + 1e-3 );
      double matchRangeData = hData -> Integral( bin1, bin2);
      double matchRangeFit = hFit -> Integral( bin1, bin2);
      
      hFit->Scale(matchRangeData/matchRangeFit);
      
      double integralFit = hFit->Integral(1,ii);
      double integralData = hData->Integral(ii+1,hData->GetNbinsX()+1);
      double integralAll = integralFit+integralData;
      
      cout<<"at bin #"<<ii<<", integrated up to "<<hData->GetBinLowEdge(ii+1)<<" fraction above this value is: "<<integralData/integralAll<<endl;
      anchorPoint = hData->GetBinLowEdge(ii+1);

      if(integralData/integralAll<.9) break;
    }
    
    for(int ii=1; ii<hData->GetNbinsX()+1; ii++){
      // renormalize data curve
      if( hData->GetBinCenter(ii) < anchorPoint) hStitched->SetBinContent(ii, hFit->GetBinContent(ii));
    }
    
    cout<<"Anchor point determined to be: "<<anchorPoint<<endl;
    cout<<"Preparing stitched histogram ... "<<endl;
  */

  hFit->SetLineColor(kRed);
  hStitched->SetLineColor(kBlue);
  
  hData->GetYaxis()->SetTitleSize(0.055);
  hData->GetXaxis()->SetTitleSize(0.055);
  hData->GetYaxis()->SetLabelSize(0.04);
  hData->GetXaxis()->SetLabelSize(0.04);
  hData->SetTitle("");
  hData->Draw("hist");
  hFit->Draw("hist same");
  hStitched->Draw("hist same");
  
  //All fine, let's try the calibrator
  multCalibrator *lCalib = new multCalibrator("lCalib");
  lCalib->SetAnchorPointPercentage(100.0f);
  lCalib->SetAnchorPointRaw(-1e-6);
  
  //Set standard Pb-Pb boundaries
  lCalib->SetStandardOnePercentBoundaries();
  
  TString calibFileName = lInputFileName.Data();
  calibFileName.ReplaceAll("glauberNBD", "calibration");
  calibFileName.ReplaceAll(".root", "_LightIonDef.root");
  TFile *fileCalib = new TFile(calibFileName.Data(), "RECREATE");
  
  TH1F *hCalib = lCalib->GetCalibrationHistogram(hStitched, "hCalib");
  
  TCanvas *c2 = new TCanvas("c2", "", 800, 600);
  c2->SetLeftMargin(0.17);
  c2->SetBottomMargin(0.17);
  c2->SetRightMargin(0.15);
  c2->SetTopMargin(0.05);
  c2->SetTicks(1,1);
//  c2->SetLogz();
  c2->SetFrameFillStyle(0);
  c2->SetFillStyle(0);
  
  hCalib->GetYaxis()->SetTitleSize(0.055);
  hCalib->GetXaxis()->SetTitleSize(0.055);
  hCalib->GetYaxis()->SetLabelSize(0.04);
  hCalib->GetXaxis()->SetLabelSize(0.04);
  hCalib->SetTitle("");
  hCalib->Draw();

  if (doNpartNcoll) {
    cout<<"Will now attempt to calculate % -> Np, Nc map..."<<endl;
    
    TProfile *hProfileNpart = new TProfile("hProfileNpart", "", 100, 0, 100);
    TProfile *hProfileNcoll = new TProfile("hProfileNcoll", "", 100, 0, 100);
    TH2F *h2dNpart = new TH2F("h2dNpart", "", 100, 0, 100, 500, -0.5f, 499.5f);
    TH2F *h2dNcoll = new TH2F("h2dNcoll", "", 100, 0, 100, 3000, -0.5f, 2999.5);
    TH2F *h2dEcc = new TH2F("h2dEcc", "", 100, 0, 100, 100, 0.0f, 1.0f);
    TH2F *h2dB = new TH2F("h2dB", "", 100, 0, 100, 100, 0.0f, 20.0f);
    TH2F *h2dNancestors = new TH2F("h2dNancestors", "", 100, 0, 100, 900, -0.5f, 899.5f);
    
    
    // Replay
    multGlauberNBDFitter *g = new multGlauberNBDFitter("lglau");
    TF1 *fitfunc = g->GetGlauberNBD();
    
    //Step 1: open the (Npart, Ncoll) pair information, provide
    TFile *fbasefile = new TFile(Form("basehistos_%s_%s_omega%.1f.root", nuclearProfile.Data(), nuclearProfile.Data(), omega),"READ");
    TH2D *hNpNc = (TH2D*) fbasefile->Get("hNpNc");
    TH3D *hNpNcEcc = (TH3D*) fbasefile->Get("hNpNcEcc");
    TH3D *hNpNcB = (TH3D*) fbasefile->Get("hNpNcB");
    g->SetNpartNcollCorrelation(hNpNc);
    g->InitializeNpNc();
    
    fitfunc->SetParameter(0, hGlauberParameters->GetBinContent(1));
    fitfunc->SetParameter(1, hGlauberParameters->GetBinContent(2));
    fitfunc->SetParameter(2, hGlauberParameters->GetBinContent(3));
    fitfunc->SetParameter(3, hGlauberParameters->GetBinContent(4));
    fitfunc->SetParameter(4, hGlauberParameters->GetBinContent(5));
    
    //  hGlauberParameters -> SetBinContent( 1, fitfunc -> GetParameter(0));
    //  hGlauberParameters -> SetBinContent( 2, fitfunc -> GetParameter(1));
    //  hGlauberParameters -> SetBinContent( 3, fitfunc -> GetParameter(2));
    //  hGlauberParameters -> SetBinContent( 4, fitfunc -> GetParameter(3));
    //  hGlauberParameters -> SetBinContent( 5, fitfunc -> GetParameter(4));
    
    Double_t lMax = hData->GetBinLowEdge( hData->GetNbinsX() + 1);
    g->CalculateAvNpNc( hProfileNpart, hProfileNcoll, h2dNpart, h2dNcoll, hCalib, 0, lMax , hNpNcEcc, h2dEcc, hNpNcB, h2dB, h2dNancestors );
    
    fileCalib->cd();
    
    hProfileNpart->Write();
    hProfileNcoll->Write();
    h2dNpart->Write();
    h2dNcoll->Write();
  }
  
  hData->Write();
  hCalib->Write();
  hStitched->Write();
  hFit->Write();
  fileCalib->Write();
  
  
}
