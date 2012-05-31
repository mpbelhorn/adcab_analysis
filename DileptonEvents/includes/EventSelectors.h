#ifndef ANALYSISSELECTORS_H
#define ANALYSISSELECTORS_H

#include <cmath>
#include "DileptonEvents.h"
#include <ParticleIdCodes.h>

//______________________________________________________________________________
// Veiled attributes are the particle attributes such as species that can be
//   either assigned (as in veiled or unknown to be true) and the truthtable
//   values if from MC. This class exists as a way to switch between tru and asn
//   values based on a script switch (information_switch). As of 2011/07/27 it
//   has not been implemented.
class VeiledAttributes {
 public:
  VeiledAttributes() {}
  VeiledAttributes( int l0_or_l1, int asn_or_tru, DileptonEvents& entry );
  ~VeiledAttributes() {}
  
  bool electron();
  bool muon();

 private:
  DileptonEvents * en_;
  int l0_or_l1_;
  int asn_or_tru_;
  int situationId_;
};



//______________________________________________________________________________
class ParticleAttributes {
 public:
  ParticleAttributes() {}
  ParticleAttributes(int l0_or_l1, DileptonEvents& entry);
  ~ParticleAttributes() {}
  
  VeiledAttributes asn;
  VeiledAttributes tru;
  
  bool lepton();
  bool mom( int momId );
  bool source_Bs();
  bool source_Bdu();
  bool source_K();
  bool source_D();
  bool source_pi();
  bool source_gamma();
  bool source_rho();
  bool source_Jpsi();
  
  bool is_signal();
  
 private:
  DileptonEvents * en_;
  int l0_or_l1_;
  PythiaPidCodes pid_;
};

//______________________________________________________________________________
class EventAttributes {
 public:
  EventAttributes() {}
  EventAttributes( int asn_or_tru, DileptonEvents& entry );
  ~EventAttributes() {}
  
  bool dielectron();
  bool dimuon();
  bool elmu();
  char type();

 private:
  DileptonEvents * en_;
  int asn_or_tru_;
};

//______________________________________________________________________________
class AnalysisSelectors {
 public:
  AnalysisSelectors(DileptonEvents& entry);
  ~AnalysisSelectors() {}
  
  ParticleAttributes l0;
  ParticleAttributes l1;
  
  EventAttributes asn;
  EventAttributes tru;
  
  bool sameSign();
  bool oppositeSign();
  bool dilepton();
  bool isSignal();
  bool is_signal_Bs();
  bool is_signal_Bdu();
  bool is_correct_wrong();
  bool is_wrong_wrong();
  bool is_continuum();
  char species_type();
  int signal_type();
  char event_sign();

 private:
  DileptonEvents * en_;
};

#endif
