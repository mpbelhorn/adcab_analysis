#ifndef Fit3D_h
#define Fit3D_h

#include "DileptonEvents.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"

class Fit3D : public DileptonEvents {
 public:
  Fit3D(
      TString input_ntuple_file="default.root",
      TString analysis_name="UnspecifiedAnalysis");
  ~Fit3D();

  void ntupleLoopCore();
  void initializeDataSet(
      const TString &x_variable_name,
      const TString &y_variable_name,
      const TString &z_variable_name,
      const TString &x_axis_label,
      const double &min_bin_x, const double &max_bin_x, const int &num_bins_x,
      const TString &y_axis_label,
      const double &min_bin_y, const double &max_bin_y, const int &num_bins_y,
      const TString &z_axis_label,
      const double &min_bin_z, const double &max_bin_z, const int &num_bins_z);
  void fillDataSet(
      const double &x_var,
      const double &y_var,
      const double &z_var,
      const int &component);
  void saveDataSet();
  
  double x_value_;
  double y_value_;
  double z_value_;
  
  RooRealVar* x_variable_;
  RooRealVar* y_variable_;
  RooRealVar* z_variable_;

};

#endif
