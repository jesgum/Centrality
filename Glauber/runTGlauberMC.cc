#include <iostream>
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "runglauber_v3.2.C"

void runTGlauberMC()
{
  gSystem->Load("libMathMore");
  gSystem->CompileMacro("runglauber_v3.2.C","-kfo");

  Int_t nevents = 1000000;
  Double_t signn = 68.21;
  Double_t mind = 0.4;
  
  runAndSaveNtuple(nevents, "Pb",    "Pb",    signn, mind);
  runAndSaveNtuple(nevents, "Pbrw",  "Pbrw",  signn, mind);
  runAndSaveNtuple(nevents, "Pb*",   "Pb*",   signn, mind);
  runAndSaveNtuple(nevents, "PbHN",  "Pb",    signn, mind);
  runAndSaveNtuple(nevents, "Pbpn",  "Pbpn",  signn, mind);
  runAndSaveNtuple(nevents, "Pbnrw", "Pbnrw", signn, mind);
}