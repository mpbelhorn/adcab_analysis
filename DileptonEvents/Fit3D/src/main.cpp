#include "Fit3D.h"

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Not enough arguments. Usage:" << std::endl
              << "Fit3D input_ntuple analysis_name run_mode" 
              << std::endl;
    return 1;
  }
  
  // Process ntuple: generate plots and dataset.
  Fit3D psum_loa_dz_fit(argv[1], argv[2],
      "|p_{0}| + |p_{1}| (GeV/c)", 2.0, 5.25,
      "Cos(#Theta_{ll})", -0.80, 0.98,
      "#Delta z", 0, .2);

  // Process raw data.
  psum_loa_dz_fit.setCreateDataSet(true);
  psum_loa_dz_fit.processNtuple();

  // Save processed data.
  psum_loa_dz_fit.drawHistograms();
  psum_loa_dz_fit.saveHistograms();
  psum_loa_dz_fit.saveDataSet();
  
  // Generate models.
  if (*argv[3] == 'g') {
    psum_loa_dz_fit.generateModels();
  }
  else if (*argv[3] == 'f') {
    psum_loa_dz_fit.fitData();
  }
  else {
    std::cout << "Fit/generate flags not set! Third argument should be "
              << "(g)enerate or (f)it." << std::endl;
  }
  
  return 0;
}
