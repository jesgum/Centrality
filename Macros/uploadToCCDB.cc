#include <gsl/span>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <TSystem.h>
#include <TFile.h>
#include <TH1F.h>
#include <TFormula.h>
#include <TKey.h>
#include "Framework/Logger.h"
#include "CCDB/BasicCCDBManager.h"

Double_t GetBoundaryForPercentile(TH1* histo, Double_t lPercentileRequested)
{
  // This function returns the boundary for a specific percentile.
  Double_t lReturnValue = 0.0;
  Double_t lPercentile = 100.0 - lPercentileRequested;

  const Long_t lNBins = histo->GetNbinsX();
  Double_t lCountDesired = lPercentile * histo->GetEntries() / 100;
  Long_t lCount = 0;
  for (Long_t ibin = 1; ibin < lNBins; ibin++) {
    lCount += histo->GetBinContent(ibin);
    if (lCount >= lCountDesired) {
      // Found bin I am looking for!
      Double_t lWidth = histo->GetBinWidth(ibin);
      Double_t lLeftPercentile = 100. * (lCount - histo->GetBinContent(ibin)) / histo->GetEntries();
      Double_t lRightPercentile = 100. * lCount / histo->GetEntries();

      Double_t lProportion = (lPercentile - lLeftPercentile) / (lRightPercentile - lLeftPercentile);

      lReturnValue = histo->GetBinLowEdge(ibin) + lProportion * lWidth;
      break;
    }
  }
  return lReturnValue;
}

void uploadToCCDB()
{
  // Read run numbers from file
  std::vector<int> runs = {564356, 564359, 564373, 564374, 564387, 564400, 564414, 564430, 564445};

  const char* ccdbUser = "Users/j/jekarlss";
  const char* ccdbPath = "Test2/LightIon";

  int nRuns = 9;
  cout << "Processing " << nRuns << " runs..." << endl;

  o2::ccdb::CcdbApi ccdb_api;
  ccdb_api.init("https://alice-ccdb.cern.ch");
  std::map<string, string> metadataRCT, headers;

  for (int ii = 0; ii < runs.size(); ii++) {
    // for(int ii=0; ii<5; ii++){
    //  Get the desired timestamps from the CCDB, no need to suffer with independent code
    headers = ccdb_api.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runs[ii]), metadataRCT, -1);
    int64_t tsSOR = atol(headers["SOR"].c_str());
    int64_t tsEOR = atol(headers["EOR"].c_str());

    // safety margins
    tsSOR = tsSOR - 300;
    tsEOR = tsEOR + 300;

    cout << "Run " << runs[ii] << " SOR " << tsSOR << " EOR " << tsEOR << endl;

    // Interface with local objects to be uploaded
    cout << "Now performing dedicated object upload... please wait..." << endl;

    // Constructing entirety of list to be uploaded for this run
    TFile* ccdbfile = new TFile(Form("centrality-run-%d.root", runs[ii]), "read");
    TList* listHistograms = new TList();

    TIter next(ccdbfile->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
      TObject* obj = key->ReadObj(); // actually reads/creates the object
      if (obj) listHistograms->Add(obj);
    }

    //__________________________________________________________________
    // Upload calibration, please

    cout << "Defining metadata for this run..." << endl;
    map<string, string> metadata; // can be empty
    metadata.insert(std::pair{ "Description", Form("pass2-based calib for run %s: update with extra estimators", Form("%i", runs[ii])) });
    metadata.insert(std::pair{ "Author", "David Dobrigkeit Chinellato & Jesper Karlsson Gumprecht" });

    // Send off to CCDB
    cout << "Attempting CCDB upload..." << endl;
    try {
      ccdb_api.storeAsTFileAny(listHistograms, Form("%s/%s", ccdbUser, ccdbPath), metadata, tsSOR, tsEOR);
      //      ccdb_api.storeAsTFileAny(listHistograms, "Centrality/Estimators", metadata, tsSOR, tsEOR);
    } catch (std::exception const& e) {
      LOG(fatal) << "Failed at CCDB submission!";
    }
    cout << "Finished with upload of run " << runs[ii] << " update! " << endl;
  }
  cout << "Done!" << endl;
}