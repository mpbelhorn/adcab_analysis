#ifndef Fit3D_h
#define Fit3D_h

#include "DileptonEvents.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"

class Fit3D : public DileptonEvents {
 public:
  Fit3D(TTree *tree = 0);
  ~Fit3D();

  void ntupleLoopCore();
  
  double x_var;
  double y_var;
  double z_var;
  
  RooRealVar* x_var_;
  RooRealVar* y_var_;
  RooCategory* component_;
  RooCategory* event_sign_;
  RooCategory* event_species_;
  RooDataSet* data_set_;

};

#endif
