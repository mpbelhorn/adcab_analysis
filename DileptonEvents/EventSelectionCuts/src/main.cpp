#include "EventSelectionCuts.h"

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
  return 0;
}
