#ifndef Fit3D_h
#define Fit3D_h

#include "DileptonEvents.h"
#include "EventSelectors.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TH3.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

typedef vector<TH3D> HistogramListOverComponents;
typedef vector<HistogramListOverComponents> HistogramListOverEventSigns;
typedef vector<HistogramListOverEventSigns> HistogramListOverEventSpecies;

class Fit3D : public DileptonEvents {
 public:
  Fit3D(
      TString input_ntuple_file = "default.root",
      TString analysis_name = "UnspecifiedAnalysis",
      TString x_axis_label = "x variable",
      TString x_axis_unit = "",
      double min_x_bin_edge = 0,
      double max_x_bin_edge = 1,
      TString y_axis_label = "y variable",
      TString y_axis_unit = "",
      double min_y_bin_edge = 0,
      double max_y_bin_edge = 1,
      TString z_axis_label = "z variable",
      TString z_axis_unit = "",
      double min_z_bin_edge = 0,
      double max_z_bin_edge = 1);
  ~Fit3D();

  void ntupleLoopCore(const int& entry_id = 0);
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
  void generateModels(const int& interpolation_order = 0);
  void fitData(
      const TString& filename = "data.root",
      const TString& data_set = "data_set");
  void plotAsymmetry();
  void plotFitAccuracy(
      const RooDataSet& mc_data,
      const RooFitResult& fit);
  void plotFitProjection(
      const RooRealVar &independant_variable,
      const RooDataSet &data,
      const RooFitResult& fit,
      const RooAbsPdf &model,
      const RooAbsPdf &bs_pdf,
      const RooAbsPdf &bd_pdf,
      const RooAbsPdf &cw_pdf,
      const RooAbsPdf &ww_pdf,
      const RooAbsPdf &cn_pdf,
      const TString &filename);
  
  float x_value_;
  float y_value_;
  float z_value_;
  
  // Species, sign, component, variable
  HistogramListOverEventSpecies histograms_;
  
  RooRealVar* x_variable_;
  RooRealVar* y_variable_;
  RooRealVar* z_variable_;
  
  RooFitResult* pp_fit_results_;
  RooFitResult* nn_fit_results_;
  
  TString output_path_;
  
};

#endif
