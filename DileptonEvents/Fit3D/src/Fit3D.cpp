#include "Fit3D.h"
#include "EventSelectors.h"
#include <iostream>

Fit3D::Fit3D(TTree *tree)
  : DileptonEvents(tree)
{
  std::cout << "Initiliazing Fit3D class." << std::endl;
}

Fit3D::~Fit3D()
{
  // Blank for now.
}

void Fit3D::ntupleLoopCore()
{
  x_var = l0_pcm + l1_pcm;
  y_var = cos_thta;
  z_var = l0_ip_dz - l1_ip_dz;
  
  AnalysisSelectors cuts(*this);  
  bool select_parent[] = { 
      cuts.is_signal_Bs() || cuts.is_signal_Bdu(),
      cuts.is_correct_wrong(),
      cuts.is_wrong_wrong(),
      cuts.is_continuum()};
  
  for (int i = 0; i < 4; i++) {
    std::cout << select_parent[i] << ", ";
  }
  std::cout << std::endl;
  std::cout << x_var << ", " << y_var << ", " << z_var << std::endl;
}
