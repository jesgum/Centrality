#!/bin/bash

root.exe -b -q runTGlauberMC.cc\(\"Pb\",\"\Pb\"\)&
root.exe -b -q runTGlauberMC.cc\(\"Pbrw\",\"\Pbrw\"\)&
root.exe -b -q runTGlauberMC.cc\(\"Pb*\",\"\Pb*\"\)&
root.exe -b -q runTGlauberMC.cc\(\"PbHN\",\"\PbHN\"\)&
root.exe -b -q runTGlauberMC.cc\(\"Pbpn\",\"\Pbpn\"\)&
root.exe -b -q runTGlauberMC.cc\(\"Pbpnrw\",\"\Pbnrw\"\)&

wait
pingme "done"