#ifndef Fit3D_h
#define Fit3D_h

#include "DileptonEvents.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include <vector>

typedef vector<TH1F> HistogramListOverVariables;
typedef vector<HistogramListOverVariables> HistogramListOverComponents;
typedef vector<HistogramListOverComponents> HistogramListOverEventSigns;
typedef vector<HistogramListOverEventSigns> HistogramListOverEventSpecies;

class Fit3D : public DileptonEvents {
 public:
  Fit3D(
      TString input_ntuple_file="default.root",
      TString analysis_name="UnspecifiedAnalysis",
      TString x_axis_label = "x variable",
      double min_x_bin_edge = 0,
      double max_x_bin_edge = 1,
      TString y_axis_label = "y variable",
      double min_y_bin_edge = 0,
      double max_y_bin_edge = 1,
      TString z_axis_label = "z variable",
      double min_z_bin_edge = 0,
      double max_z_bin_edge = 1);
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
  void fillDataSet(const int &component);
  void fillHistograms(const int& component);
  void drawHistograms();
  void saveHistograms(const TString& filename = "histograms.root");
  void generateModels();
  
  float x_value_;
  float y_value_;
  float z_value_;
  
  // Species, sign, component, variable
  HistogramListOverEventSpecies histograms_;
  
  RooRealVar* x_variable_;
  RooRealVar* y_variable_;
  RooRealVar* z_variable_;
};

#endif
