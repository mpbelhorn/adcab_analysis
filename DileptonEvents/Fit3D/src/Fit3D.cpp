#include "Fit3D.h"
#include "EventSelectors.h"
#include <iostream>

Fit3D::Fit3D(TString input_ntuple_file, TString analysis_name)
  : DileptonEvents(input_ntuple_file, analysis_name)
{
  std::cout << "Initiliazing Fit3D class." << std::endl;
}

Fit3D::~Fit3D()
{
  // Blank for now.
}

void Fit3D::ntupleLoopCore()
{
  x_value_ = l0_pcm + l1_pcm;
  y_value_ = cos_thta;
  z_value_ = l0_ip_dz - l1_ip_dz;
  
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
  std::cout << x_value_ << ", " << y_value_ << ", " << z_value_ << std::endl;
}

void Fit3D::fillDataSet(const double &x_var,
                                const double &y_var,
                                const int &component)
{
  if (flag_output_a_dataset_ && data_set_) {
    x_var_->setVal(x_var);
    y_var_->setVal(y_var);
    component_->setIndex(component);
    event_sign_->setIndex(evt_sign);
    event_species_->setIndex(typ_tru);
    data_set_->add(RooArgSet(
        *x_var_, *y_var_, *component_, *event_sign_, *event_species_));
  }
}