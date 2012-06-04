#include "Fit3D.h"
#include "EventSelectors.h"
#include <iostream>

Fit3D::Fit3D(
    TString input_ntuple_file,
    TString analysis_name,
    TString x_axis_label,
    double min_x_bin_edge,
    double max_x_bin_edge,
    TString y_axis_label,
    double min_y_bin_edge,
    double max_y_bin_edge,
    TString z_axis_label,
    double min_z_bin_edge,
    double max_z_bin_edge)
  : DileptonEvents(input_ntuple_file, analysis_name)
{
  std::cout << "Initiliazing Fit3D class." << std::endl;
  // Define the data set elements and add them to the dataset.
  x_variable_ = new RooRealVar(
      "x_variable", x_axis_label, min_x_bin_edge, max_x_bin_edge);
  y_variable_ = new RooRealVar(
      "y_variable", y_axis_label, min_y_bin_edge, max_y_bin_edge);
  z_variable_ = new RooRealVar(
      "z_variable", z_axis_label, min_z_bin_edge, max_z_bin_edge);
  std::cout << "Adding new columns to dataset." << std::endl;
  RooArgSet analysis_variables(*x_variable_, *y_variable_, *z_variable_);
  data_set_->printArgs(cout);
  std::cout << std::endl;
  recreateDataSet(analysis_variables);
  data_set_->printArgs(cout);
  std::cout << std::endl;
}

Fit3D::~Fit3D()
{
  std::cout << "Starting destructor." << std::endl;
  // TODO member pointers need to be deleted, but these lines give segfault.
  // delete x_variable_;
  // delete y_variable_;
  // delete y_variable_;
}

void Fit3D::ntupleLoopCore()
{
  x_value_ = l0_pcm + l1_pcm;
  y_value_ = cos_thta;
  if (evt_sign != 0) {
    z_value_ = abs(l0_ip_dz - l1_ip_dz);
  } else {
    z_value_ = l0_ip_dz - l1_ip_dz;
  }
  
  AnalysisSelectors cuts(*this);  
  bool select_parent[] = { 
      cuts.is_signal_Bs() || cuts.is_signal_Bdu(),
      cuts.is_correct_wrong(),
      cuts.is_wrong_wrong(),
      cuts.is_continuum()};
      
  fillDataSet(cuts.signal_type());
}

void Fit3D::fillDataSet(const int &component)
{ 
  // WARNING - Missing columns will silently not be added to the data set!
  if (flag_output_a_dataset_ && data_set_) {
    x_variable_->setVal(x_value_);
    y_variable_->setVal(y_value_);
    z_variable_->setVal(z_value_);
    component_->setIndex(component);
    event_sign_->setIndex(evt_sign);
    event_species_->setIndex(typ_tru);
    data_set_->add(
        RooArgSet(
            *x_variable_,
            *y_variable_,
            *z_variable_,
            *component_,
            *event_sign_,
            *event_species_));
  }
}

void DileptonEvents::drawHistograms()
{
  TH1* hist = data_set_->createHistogram("name", x_variable_);
  
  
}