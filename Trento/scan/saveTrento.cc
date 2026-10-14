#include "TrentoEvent.h"

#include "TH2D.h"
#include "TFile.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

using Array3D = std::array<double, 3>;

TH2D* createTH2D(const char* name, const char* title, Array3D xaxis, Array3D yaxis)
{
  return new TH2D(name, title, xaxis[0], xaxis[1], xaxis[2], yaxis[0], yaxis[1], yaxis[2]);
}

TH1D* createTH1D(const char* name, const char* title, Array3D axis)
{
  return new TH1D(name, title, axis[0], axis[1], axis[2]);
}

/// @param datFile trento text output to read
/// @param outfileName output file with histograms used by runTrentoFit/runTrentoCalib
void saveTrento(const char* datFile = "PbPb.dat", const char* outfileName = "outfile.root")
{
  std::vector<trento::Event> events = trento::parse::text(datFile);
  // entropy axis range from the events, it scales with the trento normalization
  double multMax = 0;
  for (const auto& ev : events) {
    multMax = std::max(multMax, ev.mult());
  }
  multMax = std::ceil(1.05 * multMax / 100.) * 100.;
  std::cout << "Entropy axis range taken from events: [0, " << multMax << "]" << std::endl;
  Array3D multAxis = {5000, 0, multMax};
  Array3D npartAxis = {500, -0.5, 499.5};
  Array3D ncollAxis = {2500, -0.5, 2499.5};

  TH1D* hMult = createTH1D("hMult", "hMult;Multiplicity;Count", multAxis);

  TH2D* hMultVsNpart = createTH2D("hMultVsNpart", "hMultVsNpart;Multiplicity;Npart", multAxis, npartAxis);
  TH2D* hMultVsNcoll = createTH2D("hMultVsNcoll", "hMultVsNcoll;Multiplicity;Ncoll", multAxis, ncollAxis);
  TH2D* hNpartVsNColl = createTH2D("hNpartVsNColl", "hNpartVsNColl;Npart;Ncoll", npartAxis, ncollAxis);
  
  for (size_t iev = 0; iev < events.size(); ++iev) {
    trento::Event& ev = events[iev];
    hMult->Fill(ev.mult());
    hMultVsNpart->Fill(ev.mult(), ev.npart());
    hMultVsNcoll->Fill(ev.mult(), ev.ncoll());
    hNpartVsNColl->Fill(ev.npart(), ev.ncoll());
  }

  TFile* outfile = new TFile(outfileName, "recreate");
  hMult->GetXaxis()->SetTitle("Entropy");
  hMult->GetYaxis()->SetTitle("Frequency");
  hMult->SetName("hEntropy");
  hMult->Write();
  hNpartVsNColl->SetName("hNpNc");
  hNpartVsNColl->Write();
  hMultVsNpart->SetName("hNsNp");
  hMultVsNpart->Write();
  hMultVsNcoll->SetName("hNsNc");
  hMultVsNcoll->Write();
  outfile->Write();
  outfile->Close();
}