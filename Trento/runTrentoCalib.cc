// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
/// \file saveCorrelation.C
/// \brief
/// \author ALICE
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
/// @brief function to calibrate centrality
/// @param lInputFileName name of input file.
/// @param anchorPointPercentage anchor point percentage to use
/// @param matchRange width of region in which data/glauber matching is to be done in rolling anchoring test
/// @param doNpartNcoll wether or not to attempt calculating Npart, Ncoll in centrality bins
void runTrentoCalib(TString lInputFileName = "trentofit.root", TString lBasefile = "outfile.root")
{

  double anchorPointPercentage = 90.f;
  Bool_t globalNormalisation = kFALSE;
  double matchRange = 200.0;
  bool doNpartNcoll = true;

  TFile* file = new TFile(lInputFileName.Data(), "READ");
  file->ls();

  TH1F* hData = (TH1F*)file->Get("hV0MUltraFine");
  if (!hData) {
    hData = (TH1F*)file->Get("hData");
  }

  TH1F* hGlauberParameters = (TH1F*)file->Get("hGlauberParameters");
  TH1F* hGlauberFitRange = (TH1F*)file->Get("hGlauberFitRange");
  hData->SetName("hData");
  TH1F* hStitched = (TH1F*)hData->Clone("hStitched");
  TH1F* hFit = (TH1F*)file->Get("hGlauber");

  TCanvas* c1 = new TCanvas("c1", "", 800, 600);
  c1->SetLeftMargin(0.17);
  c1->SetBottomMargin(0.17);
  c1->SetRightMargin(0.15);
  c1->SetTopMargin(0.05);
  c1->SetTicks(1, 1);
  c1->SetLogz();
  c1->SetFrameFillStyle(0);
  c1->SetFillStyle(0);

  cout << "Data bin width: " << hData->GetBinWidth(1) << endl;
  cout << "Fit bin width: " << hFit->GetBinWidth(1) << endl;
  cout << "Match range to use: " << matchRange << endl;

  //____________________________________________
  double anchorPointFraction = anchorPointPercentage / 100.f;
  double anchorPoint = -1; // the anchor point value in raw

  //____________________________________________
  // doing partial integration up to certain point for finding anchor point bin
  for (int ii = 1; ii < hData->GetNbinsX() + 1; ii++) {
    // renormalize data curve
    int bin1 = ii + 1;
    int bin2 = (globalNormalisation) ? hData->GetNbinsX() : hData->FindBin(hData->GetBinLowEdge(ii + 1) + matchRange + 1e-3);
    double matchRangeData = hData->Integral(bin1, bin2);
    double matchRangeFit = hFit->Integral(bin1, bin2);

    // rescale fit to match in the vicinity of the region we're at
    hFit->Scale(matchRangeData / matchRangeFit);

    double integralFit = hFit->Integral(1, ii);
    double integralData = hData->Integral(ii + 1, hData->GetNbinsX() + 1);
    double integralAll = integralFit + integralData;

    cout << "at bin #" << ii << ", integrated up to " << hData->GetBinLowEdge(ii + 1) << " fraction above this value is: " << integralData / integralAll << endl;
    anchorPoint = hData->GetBinLowEdge(ii + 1);

    if (integralData / integralAll < anchorPointFraction) {
      break;
    }
  }

  //____________________________________________
  for (int ii = 1; ii < hData->GetNbinsX() + 1; ii++) {
    // renormalize data curve
    if (hData->GetBinCenter(ii) < anchorPoint) {
      hStitched->SetBinContent(ii, hFit->GetBinContent(ii));
    }
  }

  cout << "Anchor point determined to be: " << anchorPoint << endl;
  cout << "Preparing stitched histogram ... " << endl;

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

  // All fine, let's try the calibrator
  multCalibrator* lCalib = new multCalibrator("lCalib");
  lCalib->SetAnchorPointPercentage(100.0f);
  lCalib->SetAnchorPointRaw(-1e-6);

  // Set standard Pb-Pb boundaries
  lCalib->SetStandardOnePercentBoundaries();

  TString calibFileName = lInputFileName.Data();
  calibFileName.ReplaceAll(".root", "_calib.root");
  TFile* fileCalib = new TFile(calibFileName.Data(), "RECREATE");
  // TFile* fileCalib = new TFile(outfile, "RECREATE");

  TH1F* hCalib = lCalib->GetCalibrationHistogram(hStitched, "hCalib");

  TCanvas* c2 = new TCanvas("c2", "", 800, 600);
  c2->SetLeftMargin(0.17);
  c2->SetBottomMargin(0.17);
  c2->SetRightMargin(0.15);
  c2->SetTopMargin(0.05);
  c2->SetTicks(1, 1);
  //  c2->SetLogz();
  c2->SetFrameFillStyle(0);
  c2->SetFillStyle(0);

  hCalib->GetYaxis()->SetTitleSize(0.055);
  hCalib->GetXaxis()->SetTitleSize(0.055);
  hCalib->GetYaxis()->SetLabelSize(0.04);
  hCalib->GetXaxis()->SetLabelSize(0.04);
  hCalib->SetTitle("");
  hCalib->Draw();

  fileCalib->cd();

  hData->Write();
  hCalib->Write();
  hStitched->Write();
  hFit->Write();

  if (doNpartNcoll) {
    cout << "Will now attempt to calculate % -> Np, Nc map..." << endl;

    TProfile* hProfileNpart = new TProfile("hProfileNpart", "", 100, 0, 100);
    TProfile* hProfileNcoll = new TProfile("hProfileNcoll", "", 100, 0, 100);
    TH2D* h2dNpart = new TH2D("h2dNpart", "", 100, 0, 100, 500, -0.5f, 499.5f);
    TH2D* h2dNcoll = new TH2D("h2dNcoll", "", 100, 0, 100, 3000, -0.5f, 2999.5);

    // Replay
    multGlauberNBDFitter* g = new multGlauberNBDFitter("lglau");
    g->InitTrentoNBD();
    TF1* fitfunc = g->GetTrentoNBD();

    // Step 1: open the (Npart, Ncoll) pair information, provide
    TFile* fbasefile = new TFile(lBasefile.Data(), "READ");
    // TFile* fbasefile = new TFile("outfileJesper.root", "READ");
    TH1F* hNSources = (TH1F*)fbasefile->Get("hEntropy");
    TH2D* hNsNp = (TH2D*)fbasefile->Get("hNsNp");
    TH2D* hNsNc = (TH2D*)fbasefile->Get("hNsNc");

    Double_t lMax = hData->GetBinLowEdge(hData->GetNbinsX() + 1);
    fitfunc->SetParameter(0, hGlauberParameters->GetBinContent(1));
    fitfunc->SetParameter(1, hGlauberParameters->GetBinContent(2));
    fitfunc->SetParameter(2, hGlauberParameters->GetBinContent(3));
    fitfunc->SetParameter(3, hGlauberParameters->GetBinContent(4));
    fitfunc->SetParameter(4, hGlauberParameters->GetBinContent(5));

    TF1* fNBD = g->GetNBD();
    const Double_t muPerSrc = hGlauberParameters->GetBinContent(1);
    const Double_t kPerSrc = hGlauberParameters->GetBinContent(2);

    const TAxis* xCent = hProfileNpart->GetXaxis(); // assumes h2dNpart/h2dNcoll/hProfileNcoll share this x binning
    const int nCB = xCent->GetNbins() + 2;          // incl. under/overflow
    const int nNp = hNsNp->GetNbinsY();
    const int nNc = hNsNc->GetNbinsY();

    // 1) Centrality bin of each multiplicity: source-independent, compute once
    std::vector<int> centBin(lMax, 0);
    for (int mult = 1; mult < lMax; ++mult) {
      Double_t cent = hCalib->GetBinContent(hCalib->FindBin(mult));
      centBin[mult] = xCent->FindFixBin(cent);
    }

    // Accumulators: weight[centBin][yBin]
    std::vector<double> wNp(nCB * (nNp + 1), 0.0), wNc(nCB * (nNc + 1), 0.0);
    std::vector<double> pCent(nCB), cNp(nNp + 1), cNc(nNc + 1);

    for (int isrc = 1; isrc <= hNSources->GetNbinsX(); ++isrc) {
      const Double_t pSource = hNSources->GetBinContent(isrc);
      if (pSource <= 0) {
        continue;
      }

      const Double_t nSources = hNSources->GetBinCenter(isrc);
      const Double_t thisMu = nSources * muPerSrc;
      const Double_t thisK = nSources * kPerSrc;
      const Double_t pVal = 1.0 / (1.0 + thisMu / thisK);
      fNBD->SetParameter(0, pVal);  // NBDPar::p
      fNBD->SetParameter(1, thisK); // NBDPar::k
      // 2) Sum P(mult | src) per centrality bin
      std::fill(pCent.begin(), pCent.end(), 0.0);
      for (int mult = 1; mult < lMax; ++mult)
        pCent[centBin[mult]] += fNBD->Eval(mult);

      // Cache the Np/Nc columns for this source
      Double_t sumNp = 0.0, sumNc = 0.0;
      for (int i = 1; i <= nNp; ++i)
        sumNp += (cNp[i] = hNsNp->GetBinContent(isrc, i));
      for (int i = 1; i <= nNc; ++i)
        sumNc += (cNc[i] = hNsNc->GetBinContent(isrc, i));

      // 3) Outer product into plain arrays
      for (int cb = 0; cb < nCB; ++cb) {
        if (pCent[cb] == 0)
          continue;
        if (sumNp > 0) {
          const double a = pSource * pCent[cb] / sumNp;
          double* row = &wNp[cb * (nNp + 1)];
          for (int i = 1; i <= nNp; ++i)
            row[i] += a * cNp[i];
        }
        if (sumNc > 0) {
          const double a = pSource * pCent[cb] / sumNc;
          double* row = &wNc[cb * (nNc + 1)];
          for (int i = 1; i <= nNc; ++i)
            row[i] += a * cNc[i];
        }
      }
    }

    // 4) One Fill per non-empty (cent, y) cell
    for (int cb = 0; cb < nCB; ++cb) {
      const Double_t cent = xCent->GetBinCenter(cb);
      for (int i = 1; i <= nNp; ++i) {
        const double w = wNp[cb * (nNp + 1) + i];
        if (w == 0)
          continue;
        const Double_t nPart = hNsNp->GetYaxis()->GetBinCenter(i);
        hProfileNpart->Fill(cent, nPart, w);
        h2dNpart->Fill(cent, nPart, w);
      }
      for (int i = 1; i <= nNc; ++i) {
        const double w = wNc[cb * (nNc + 1) + i];
        if (w == 0)
          continue;
        const Double_t nColl = hNsNc->GetYaxis()->GetBinCenter(i);
        hProfileNcoll->Fill(cent, nColl, w);
        h2dNcoll->Fill(cent, nColl, w);
      }
    }

    fileCalib->cd();
    hProfileNpart->Write();
    hProfileNcoll->Write();
    h2dNpart->Write();
    h2dNcoll->Write();
  }

  TH1F* hMatchRange = new TH1F("hMatchRange", "", 1, 0, 1);
  TH1F* hAnchorPoint = new TH1F("hAnchorPoint", "", 1, 0, 1);
  hMatchRange->SetBinContent(1, matchRange);
  hMatchRange->Write();
  hAnchorPoint->SetBinContent(1, anchorPoint);
  hAnchorPoint->Write();
  fileCalib->Write();
}