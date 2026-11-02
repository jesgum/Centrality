// Same plots as drawScan.cc, but only for the scan points listed in a scan<tag>_best.txt written by drawScan.cc,
// so the full scan directory does not have to be read and ranked again:
//   <prefix>_Fits.pdf, <prefix>_Entropy.pdf, <prefix>_Npart.pdf, <prefix>_Ncoll.pdf, <prefix>_Parameters.pdf
// with <prefix> the name of the list file without .txt, e.g. scan_best_Fits.pdf. The plots are written next to
// the list file, e.g. into Trento/madhatter for ../madhatter/scan_best.txt.
//
// The scan point directories in the list are absolute paths on the machine that ran drawScan.cc. A directory
// that does not exist here is looked up next to the list file instead (where Trento/madhatter/fetchBest.sh puts it),
// and the reference calibration file also in data/.
//
// Usage: root -l -b -q 'drawBest.cc+("../madhatter/scan_best.txt")'

#include "drawScan.cc"

/// @param listFile scan<tag>_best.txt from drawScan.cc
/// @param reference overrides the reference stored in the list: a calibration file (*.root) or a scan point directory
/// @param referenceLabel legend entry of an overriding reference
void drawBest(TString listFile = "../madhatter/scan_best.txt", TString reference = "", TString referenceLabel = "Reference")
{
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);

  std::ifstream list(listFile.Data());
  if (!list) {
    std::cerr << "Error: cannot open " << listFile << std::endl;
    return;
  }

  // the path as listed if it exists, otherwise the file or directory of that name next to the list file,
  // or in data/ (where the reference calibration file is kept when running from Trento/scan)
  const TString listDir = gSystem->GetDirName(listFile);
  auto localPath = [&](const TString& path) -> TString {
    const TString name = gSystem->BaseName(path);
    for (const TString& candidate : { path, listDir + "/" + name, "data/" + name }) {
      if (!gSystem->AccessPathName(candidate)) {
        return candidate;
      }
    }
    return path;
  };

  // "# reference <tab> path <tab> label", then "directory <tab> deviation (%) <tab> label" per scan point
  std::vector<ScanPoint> points;
  std::vector<std::pair<std::string, double>> ranking; // as in the list: most accurate first
  TString listReference, listReferenceLabel;
  std::string line;
  while (std::getline(list, line)) {
    TObjArray* tokens = TString(line).Tokenize("\t");
    const int n = tokens->GetEntries();
    const TString first = n > 0 ? tokens->At(0)->GetName() : "";
    if (first == "# reference" && n >= 2) {
      listReference = tokens->At(1)->GetName();
      listReferenceLabel = n >= 3 ? tokens->At(2)->GetName() : "Reference";
    } else if (n >= 1 && !first.BeginsWith("#")) {
      ScanPoint point;
      point.name = gSystem->BaseName(first);
      point.label = n >= 3 ? TString(tokens->At(2)->GetName()) : makeLabel(point.name);
      if (loadPoint(localPath(first), point)) {
        points.push_back(point);
        ranking.emplace_back(point.name, n >= 2 ? TString(tokens->At(1)->GetName()).Atof() / 100 : 0);
      } else {
        std::cout << "Skipping " << first << ": could not read its output" << std::endl;
      }
    }
    delete tokens;
  }
  if (points.empty()) {
    std::cerr << "Error: no readable scan points in " << listFile << std::endl;
    return;
  }

  if (reference.IsNull()) {
    reference = listReference;
    referenceLabel = listReferenceLabel;
  }
  if (reference.IsNull()) {
    std::cerr << "Warning: no reference in " << listFile << ", using " << points[0].name << std::endl;
  }

  // The reference first, then the scan points in the order of the list
  std::vector<ScanPoint> centPoints;
  const bool referenceFromFile = reference.EndsWith(".root");
  if (!reference.IsNull()) {
    reference = localPath(reference);
  }
  if (!reference.IsNull()) {
    ScanPoint ref;
    ref.name = referenceFromFile ? "reference" : gSystem->BaseName(reference);
    ref.label = referenceFromFile ? referenceLabel : makeLabel(ref.name);
    if (!(referenceFromFile ? loadReference(reference, ref) : loadPoint(reference, ref))) {
      std::cerr << "Error: cannot read the reference " << reference << std::endl;
      return;
    }
    centPoints.push_back(ref);
  }
  centPoints.insert(centPoints.end(), points.begin(), points.end());

  printSummary(centPoints, referenceFromFile);
  printf("Deviation from %s (as stored in %s):\n", centPoints[0].label.Data(), listFile.Data());
  for (const auto& entry : ranking) {
    printf("  %-40s %9.2f%%\n", entry.first.c_str(), 100 * entry.second);
  }

  TString prefix = listDir + "/" + gSystem->BaseName(listFile);
  prefix.ReplaceAll(".txt", "");
  if (!drawPlots(centPoints, referenceFromFile, prefix)) {
    return;
  }
  drawParametersVsRank(ranking, -1, prefix + "_Parameters.pdf");
}
