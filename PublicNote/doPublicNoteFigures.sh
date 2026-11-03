#!/bin/bash
# Produce all public-note figures from the .root files in Files/ into Figures/

cd "$(dirname "$0")" || exit 1
mkdir -p Figures

root -l -b -q plotAnchorPointStability.cc
root -l -b -q plotNpartNcoll.cc
root -l -b -q plothGlauber.cc
root -l -b -q plotCollVsBC.cc
root -l -b -q plotNpartNcollVsRun.cc
root -l -b -q plotNpartNcollNuclearStructure.cc
root -l -b -q plotAnchorPointShift.cc
root -l -b -q plotNpartNcollOOVsRun.cc
root -l -b -q plotImpactParameter.cc
root -l -b -q plotTGlauberMC.cc
