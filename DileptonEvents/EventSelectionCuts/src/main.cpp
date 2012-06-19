#include "EventSelectionCuts.h"
#include "SelectBestCandidates.h"

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Not enough arguments. Usage:" << std::endl
              << "EventSelectionCuts input_ntuple analysis_name" 
              << std::endl;
    return 1;
  }
  EventSelectionCuts basf_cuts(argv[1], argv[2]);
  basf_cuts.processNtuple();
  basf_cuts.saveNewNtuple();
  
  SelectBestCandidates best_candidate_selection(
      TString(argv[1]) + TString(".postbasf.root"), argv[2]);
  best_candidate_selection.processNtuple();
  best_candidate_selection.isolateBestCandidates();
  best_candidate_selection.saveNewNtuple();
  
  return 0;
}