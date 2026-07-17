#include "Includes/EventSelectionParams.h"
#include "ALICE3/Core/FlatTrackSmearer.h"
#include "ALICE3/Core/FlatLutEntry.h"

#include <map>

void printEventSelectionParams()
{
  EventSelectionParams* par = nullptr;
  constexpr long Timestamp = 175172992900; // Run 564374


  // CCDB boilerplate declarations
  const char* ccdbUrl = "http://alice-ccdb.cern.ch";
  o2::ccdb::CcdbApi ccdb;
  ccdb.init("https://alice-ccdb.cern.ch");
  std::map<std::string, std::string> metadata; // can be empty

  par = ccdb.retrieveFromTFileAny<EventSelectionParams>("EventSelection/EventSelectionParams", metadata, Timestamp);


  std::cout << "par->fZNABBlower: " << par->fZNABBlower << std::endl;;
  std::cout << "par->fZNABBupper: " << par->fZNABBupper << std::endl;;
  
}