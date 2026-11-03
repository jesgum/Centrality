#ifndef SAVEHELPER_H_
#define SAVEHELPER_H_

#include <cstdio>
#include <filesystem>

#include "TDirectory.h"
#include "TROOT.h"
#include "TFile.h"
#include "TFrame.h"
#include "TPad.h"
#include "TString.h"

namespace save {

// Write all objects drawn in a pad (and its sub-pads) to the file
inline void writePrimitives(TPad* pad, TFile* f)
{
  for (const auto& obj : *pad->GetListOfPrimitives()) {
    if (obj->InheritsFrom(TFrame::Class())) {
      continue;
    }
    if (obj->InheritsFrom(TPad::Class())) {
      writePrimitives((TPad*)obj, f);
      continue;
    }
    f->WriteTObject(obj, obj->GetName());
  }
}

// Save the canvas and its objects to a .root file next to the figure,
// e.g. "Figures/hNpart.pdf" -> "Figures/hNpart.root", so the plot can be
// re-opened and restyled in another macro.
inline void canvasToRoot(TPad* canv, const char* figureName)
{
  TString rootName = figureName;
  Ssiz_t dot = rootName.Last('.');
  if (dot != kNPOS) {
    rootName.Remove(dot);
  }
  rootName += ".root";

  // Never clobber an input: refuse if a file with this name is currently open
  const auto target = std::filesystem::weakly_canonical(rootName.Data());
  for (const auto& obj : *gROOT->GetListOfFiles()) {
    if (std::filesystem::weakly_canonical(obj->GetName()) == target) {
      printf("save::canvasToRoot: %s is an open input file, not overwriting it\n", rootName.Data());
      return;
    }
  }

  TDirectory::TContext ctx; // restore gDirectory afterwards
  TFile* f = TFile::Open(rootName, "RECREATE");
  if (!f || f->IsZombie()) {
    printf("save::canvasToRoot: cannot create %s\n", rootName.Data());
    return;
  }
  f->WriteTObject(canv, canv->GetName());
  writePrimitives(canv, f);
  f->Close();
  delete f;
}

} // namespace save

#endif // SAVEHELPER_H_
