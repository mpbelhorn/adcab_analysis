#ifndef EventSelectionCuts_h
#define EventSelectionCuts_h

#include "DileptonEvents.h"

#include <vector>
#include <iostream>

class EventSelectionCuts : public DileptonEvents {
 public:
  EventSelectionCuts();
  ~EventSelectionCuts();

  void ntupleLoopCore();
};

#endif
