#include "EventSelectionCuts.h"

EventSelectionCuts::EventSelectionCuts(
    TString input_ntuple_file,
    TString analysis_name)
  : DileptonEvents(input_ntuple_file, analysis_name)
{
  std::cout << "Initiliazing EventSelectionCuts class." << std::endl;
  // For some reason the output file needs to be made before the cloned
  //   TTree. I don't know why or how this is true.
  TString output_file_name = input_ntuple_file_ + TString(".postbasf.root");
  output_file_ = new TFile(output_file_name, "recreate");
  output_ntuple_ = fChain->CloneTree(0);
}


EventSelectionCuts::~EventSelectionCuts()
{
  // delete output_ntuple_;
}

void EventSelectionCuts::ntupleLoopCore(const int& entry_id)
{
  if (passes_hadronb() &&
      passes_fox_wolfram_r2() &&
      passes_jpsi_veto() &&
      passes_cm_momentum() &&
      passes_lepton_opening_angle()) {
    output_ntuple_->Fill();
  }
}

void EventSelectionCuts::saveNewNtuple()
{
  output_ntuple_->Write();
  output_file_->Close();
}

bool EventSelectionCuts::passes_hadronb()
{
  return hadronb >= cuts_.hadronb_min;
}

bool EventSelectionCuts::passes_fox_wolfram_r2()
{
  return fw_r2 < cuts_.fox_wolfram_r2_max;
}

bool EventSelectionCuts::passes_jpsi_veto()
{
  double delta_mass = inv_mass - 3.096;
  
  if (abs(typ_asn) == 1) {
    return (delta_mass < cuts_.jpsi_veto_dielectron_delta_m_min ||
        cuts_.jpsi_veto_dielectron_delta_m_max < delta_mass);
  } else if (abs(typ_asn) == 2) {
    return (delta_mass < cuts_.jpsi_veto_dimuon_delta_m_min ||
        cuts_.jpsi_veto_dimuon_delta_m_max < delta_mass);
  } else {
    return true;
  }
}

bool EventSelectionCuts::passes_cm_momentum()
{
  bool l0_pcm_above_low_cut = (l0_pcm > cuts_.lepton_cm_momentum_min);
  bool l1_pcm_above_low_cut = (l1_pcm > cuts_.lepton_cm_momentum_min);
  bool l0_pcm_below_high_cut = (l0_pcm < cuts_.lepton_cm_momentum_max);
  bool l1_pcm_below_high_cut = (l1_pcm < cuts_.lepton_cm_momentum_max);
  bool l0_accepted = (l0_pcm_above_low_cut && l0_pcm_below_high_cut);
  bool l1_accepted = (l1_pcm_above_low_cut && l1_pcm_below_high_cut);
  bool event_accepted = (l0_accepted && l1_accepted);
  return event_accepted;
}

bool EventSelectionCuts::passes_lepton_opening_angle()
{
  return (cuts_.dilepton_opening_angle_min < cos_thta &&
      cos_thta < cuts_.dilepton_opening_angle_max);
}