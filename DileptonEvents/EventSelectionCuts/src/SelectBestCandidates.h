#ifndef SelectBestCandidates_h
#define SelectBestCandidates_h

#include "DileptonEvents.h"

#include <sstream>
#include <vector>
#include <map>
#include <algorithm>

struct EventCandidate{
  int entry;
  int momentum;
};

typedef vector<EventCandidate > EventCandidates;
typedef map<TString, EventCandidates> UniqueEvents;

class SelectBestCandidates : public DileptonEvents {
 public:
  SelectBestCandidates(
      TString input_ntuple_file,
      TString analysis_name);
  ~SelectBestCandidates();

  void ntupleLoopCore(const int& entry_id);
  void isolateBestCandidates();
  void saveNewNtuple();
  TString floatToString(const float& value);
  
  TFile* output_file_;
  TTree* output_ntuple_;
  
  UniqueEvents unique_events_;
  
};

#endif