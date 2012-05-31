#ifndef Fit3D_h
#define Fit3D_h

#include "DileptonEvents.h"

class Fit3D : public DileptonEvents {
 public:
  Fit3D(TTree *tree = 0);
  ~Fit3D();

  void ntupleLoopCore();
  
  double x_var;
  double y_var;
  double z_var;

};

#endif
