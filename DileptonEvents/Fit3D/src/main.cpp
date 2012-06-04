#include "Fit3D.h"

int main(int argc, char *argv[])
{
  if (argc != 3) {
    std::cout << "Not enough arguments. Usage:" << std::endl
              << "Fit3D input_ntuple analysis_name" 
              << std::endl;
    return 1;
  }
  
  // Process ntuple: generate plots and dataset.
  Fit3D psum_loa_dz_fit(argv[1], argv[2],
      "psum", 2.2, 5, "loa", -1, 1, "deltaz", -.5, .5);
  std::cout << "Passed init." << std::endl;
  psum_loa_dz_fit.processNtuple();
  psum_loa_dz_fit.data_set_->Print();
  // Generate models.
  
  
  // Fit data.
  
  return 0;
  
}
