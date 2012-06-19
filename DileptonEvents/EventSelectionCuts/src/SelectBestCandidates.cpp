#include "SelectBestCandidates.h"

SelectBestCandidates::SelectBestCandidates(
    TString input_ntuple_file,
    TString analysis_name)
  : DileptonEvents(input_ntuple_file, analysis_name)
{
  std::cout << "Initiliazing SelectBestCandidates class." << std::endl;
  // For some reason the output file needs to be made before the cloned
  //   TTree. I don't know why or how this is true.
  TString output_file_name = input_ntuple_file_ + TString(".best.root");
  output_file_ = new TFile(output_file_name, "recreate");
  output_ntuple_ = fChain->CloneTree(0);
}


SelectBestCandidates::~SelectBestCandidates()
{
  // delete output_ntuple_;
}

void SelectBestCandidates::ntupleLoopCore(const int& entry_id)
{
  TString event_candidate_id =
      floatToString(stm_no) + 
      floatToString(exp_no) + 
      floatToString(run_no) + 
      floatToString(evt_no);
  
  EventCandidate event_candidate;
  event_candidate.entry = entry_id;
  event_candidate.momentum = l0_pcm + l1_pcm;
  
  unique_events_[event_candidate_id].push_back(event_candidate);
}

void SelectBestCandidates::isolateBestCandidates()
{
  vector<int> best_candidate_entries;
  UniqueEvents::iterator unique_event;
  for (unique_event = unique_events_.begin();
      unique_event != unique_events_.end(); ++unique_event) {
    EventCandidates &event_candidates = unique_event->second;
    int number_of_candidates = event_candidates.size();

    EventCandidate best_candidate;
    best_candidate.entry = event_candidates.front().entry;
    best_candidate.momentum = event_candidates.front().momentum;
    if (number_of_candidates > 1) {
      EventCandidates::iterator ith_candidate;
      for (ith_candidate = event_candidates.begin();
          ith_candidate != event_candidates.end(); ++ith_candidate) {
        if (ith_candidate->momentum > best_candidate.momentum) {
          best_candidate.entry = ith_candidate->entry;
          best_candidate.momentum = ith_candidate->momentum;
        }
      }
    }
    best_candidate_entries.push_back(best_candidate.entry);
  }
  
  // TTree filling is very slow if entries are not sequentially added.
  std::sort(best_candidate_entries.begin(), best_candidate_entries.end());
  for (vector<int>::iterator i = best_candidate_entries.begin();
      i != best_candidate_entries.end(); ++i) {
    fChain->GetEntry(*i);
    output_ntuple_->Fill();
  }
}

void SelectBestCandidates::saveNewNtuple()
{
  output_ntuple_->Write();
}

TString SelectBestCandidates::floatToString(const float& value)
{
  stringstream ss (stringstream::in | stringstream::out);
  ss << value;
  return TString(ss.str());
}