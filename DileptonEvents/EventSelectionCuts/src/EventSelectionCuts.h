#ifndef EventSelectionCuts_h
#define EventSelectionCuts_h

#include "DileptonEvents.h"

#include <vector>
#include <iostream>
struct BasfCutValues {
  static const double dilepton_opening_angle_min = -0.80; //-0.85;
  static const double dilepton_opening_angle_max =  0.98; // 0.95;
  static const double fox_wolfram_r2_max = 0.36; //0.4;
  static const double hadronb_min = 10;
  static const double jpsi_veto_dielectron_delta_m_min = -0.15; // GeV/C^2
  static const double jpsi_veto_dielectron_delta_m_max =  0.05; // GeV/C^2
  static const double jpsi_veto_dimuon_delta_m_min = -0.05; // GeV/C^2
  static const double jpsi_veto_dimuon_delta_m_max =  0.05; // GeV/C^2
  static const double lepton_cm_momentum_min = 1.01; //1.1; // GeV/c
  static const double lepton_cm_momentum_max =  2.56; //2.5; // GeV/c
  static const double muon_id_chi2_per_klm_layers_hit_max = 3.5;
  static const double lab_frame_polar_angle_cosine_min = -0.707;
  static const double lab_frame_polar_angle_cosine_max =  0.866;

  static const double dz_samesign_min = 0.00; // cm.
  static const double dz_samesign_max = 0.18; // cm.
};

class EventSelectionCuts : public DileptonEvents {
 public:
  EventSelectionCuts(
      TString input_ntuple_file,
      TString analysis_name);
  ~EventSelectionCuts();

  void ntupleLoopCore(const int& entry_id);
  void saveNewNtuple();
  
  bool passes_hadronb();
  bool passes_fox_wolfram_r2();
  bool passes_jpsi_veto();
  bool passes_cm_momentum();
  bool passes_lepton_opening_angle();
  
  TFile* output_file_;
  TTree* output_ntuple_;
  BasfCutValues cuts_;
};

#endif