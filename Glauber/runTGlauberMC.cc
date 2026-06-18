#include <iostream>
#include "TSystem.h"
#include "TString.h"
#include "TROOT.h"
#include "runglauber_v3.2.C"

void runTGlauberMC(const char* sysA, const char* sysB)
{
  gSystem->Load("libMathMore");
  gSystem->CompileMacro("runglauber_v3.2.C","-kfo");

  Int_t nevents = 1000000;
  Double_t signn = 68.21;
  Double_t mind = 0.4;
  
  runAndSaveNtuple(nevents, sysA, sysB, signn, mind);
}