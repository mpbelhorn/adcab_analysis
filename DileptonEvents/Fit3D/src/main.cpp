#include "Fit3D.h"

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Not enough arguments. Needs a file and analysis name!" 
              << std::endl;
    return 1;
  }
  
  // Process ntuple: generate plots and dataset.
  Fit3D psum_loa_dz_fit(argv[1], argv[2]);
  psum_loa_dz_fit.processNtuple();
  
  // Generate models.
  
 
  // Fit data.
  
  return 0;
  
}
