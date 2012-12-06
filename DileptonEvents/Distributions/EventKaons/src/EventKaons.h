#ifndef EventKaons_h
#define EventKaons_h

#include "DileptonEvents.h"
#include "EventSelectors.h"
#include "TLatex.h"
#include "TH1.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

typedef vector<TH1D> ComponentHistograms;
typedef vector<ComponentHistograms> KaonSignHistograms;
typedef vector<KaonSignHistograms> EventKaonHistograms;
typedef vector<EventKaonHistograms> EventSignHistograms;
typedef vector<EventSignHistograms> SpeciesHistograms;

class EventKaons : public DileptonEvents {
 public:
  EventKaons(
      TString input_ntuple_file = "default.root",
      TString analysis_name = "UnspecifiedAnalysis",
      TString x_axis_label = "x variable",
      TString x_axis_unit = "",
      double min_x_bin_edge = 0,
      double max_x_bin_edge = 1);
  ~EventKaons();

  void ntupleLoopCore(const int& entry_id = 0);
  void fillHistograms(
      const float& event_kaons,
      const int& number_of_k_minus,
      const int& component);
  void drawHistograms();
  void saveHistograms(const TString& filename = "histograms.root");
  
  float x_value_;
  
  // Species, sign, component, variable
  SpeciesHistograms histograms_;
  int max_event_kaons_;
  TString output_path_;
  
};

#endif
