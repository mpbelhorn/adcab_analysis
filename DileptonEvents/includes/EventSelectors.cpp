#include "EventSelectors.h"

VeiledAttributes::VeiledAttributes( int l0_or_l1, int asn_or_tru,
    DileptonEvents& entry )
{
  en_ = &entry;
  l0_or_l1_ = l0_or_l1;
  asn_or_tru_ = asn_or_tru;
  
  // l0_or_l1 and asn_or_tru can be 0 or 1. To make the switch cases unique,
  //   asn_or_tru is multiplied by 2 so the cases are (lepton flag first):
  //   00 = 0, 10 = 1, 01 = 2, and 11 = 3.
  situationId_ = l0_or_l1_ + ( 2 * asn_or_tru_ );
}

bool VeiledAttributes::electron()
{
  switch ( situationId_ ) {
    case 0: // l0 and asn
      return ( std::abs( (*en_).l0_idasn ) == 11 );
    case 1: // l1 and asn
      return ( std::abs( (*en_).l1_idasn ) == 11 );
    case 2: // l0 and tru
      return ( std::abs( (*en_).l0_idtru ) == 11 );
    case 3: // l1 and tru
      return ( std::abs( (*en_).l1_idtru ) == 11 );
    default:
      return false;
  }
}

bool VeiledAttributes::muon()
{
  switch ( situationId_ ) {
    case 0: // l0 and asn
      return ( std::abs( (*en_).l0_idasn ) == 13 );
    case 1: // l1 and asn
      return ( std::abs( (*en_).l1_idasn ) == 13 );
    case 2: // l0 and tru
      return ( std::abs( (*en_).l0_idtru ) == 13 );
    case 3: // l1 and tru
      return ( std::abs( (*en_).l1_idtru ) == 13 );
    default:
      return false;
  }
}

//______________________________________________________________________________

ParticleAttributes::ParticleAttributes( int l0_or_l1, DileptonEvents& entry )
{
  en_ = &entry;
  l0_or_l1_ = l0_or_l1;
  asn = VeiledAttributes( l0_or_l1, 0, entry );
  tru = VeiledAttributes( l0_or_l1, 1, entry );
}

bool ParticleAttributes::lepton()
{
  return (tru.electron() || tru.muon());
}

bool ParticleAttributes::mom( int momId )
{
  switch ( l0_or_l1_ ) {
    case 0: // l0
      return (std::abs( (*en_).l0_idmom ) == momId );
    case 1: // l1
      return (std::abs( (*en_).l1_idmom ) == momId );
    default:
      return false;
  }
}

bool ParticleAttributes::source_Bs()
{
  return ( mom(pid_.Bs0) || mom(pid_.Bs0_star) );
}

bool ParticleAttributes::source_Bdu()
{
  return ( mom(pid_.B0) ||
           mom(pid_.B0_star) ||
           mom(pid_.B_plus) ||
           mom(pid_.B_plus_star) );
}

bool ParticleAttributes::source_K()
{
  return ( mom(pid_.K0) ||
           mom(pid_.K0_star) ||
           mom(pid_.K_plus) ||
           mom(pid_.K_plus_star) );
}

bool ParticleAttributes::source_D()
{
  return ( mom(pid_.D_plus) ||
           mom(pid_.D_plus_star) ||
           mom(pid_.D0) ||
           mom(pid_.D0_star) ||
           mom(pid_.Ds_plus) ||
           mom(pid_.Ds_plus_star) );
}

bool ParticleAttributes::source_pi()
{
  return ( mom(pid_.pi0) || mom(pid_.pi_plus) );
}

bool ParticleAttributes::source_gamma()
{
  return ( mom(pid_.gamma) || mom(pid_.v_gamma) );
}

bool ParticleAttributes::source_rho()
{
  return ( mom(pid_.rho0) || mom(pid_.rho_plus) );
}

bool ParticleAttributes::source_Jpsi()
{
  return ( mom(pid_.J_psi) );
}

bool ParticleAttributes::is_signal()
{
  return ((source_Bs() || source_Bdu()) && lepton());
}

//______________________________________________________________________________

EventAttributes::EventAttributes( int asn_or_tru, DileptonEvents& entry )
{
  en_ = &entry;
  asn_or_tru_ = asn_or_tru;
}

bool EventAttributes::dielectron()
{
  switch ( asn_or_tru_ ) {
    case 0: // asn
      return ( ( std::abs( (*en_).l0_idasn ) == 11 ) &&
               ( std::abs( (*en_).l1_idasn ) == 11 ) );
    case 1: // tru
      return ( ( std::abs( (*en_).l0_idtru ) == 11 ) &&
               ( std::abs( (*en_).l1_idtru ) == 11 ) );
    default:
      return false;
  }
}

bool EventAttributes::dimuon()
{
  switch ( asn_or_tru_ ) {
    case 0: // asn
      return ( ( std::abs( (*en_).l0_idasn ) == 13 ) &&
               ( std::abs( (*en_).l1_idasn ) == 13 ) );
    case 1: // tru
      return ( ( std::abs( (*en_).l0_idtru ) == 13 ) &&
               ( std::abs( (*en_).l1_idtru ) == 13 ) );
    default:
      return false;
  }
}

bool EventAttributes::elmu()
{
  switch ( asn_or_tru_ ) {
    case 0: // asn
      return ( ( ( std::abs( (*en_).l0_idasn ) == 11 ) &&
                 ( std::abs( (*en_).l1_idasn ) == 13 ) ) || 
               ( ( std::abs( (*en_).l0_idasn ) == 13 ) &&
                 ( std::abs( (*en_).l1_idasn ) == 11 ) ) );
    case 1: // tru
      return ( ( ( std::abs( (*en_).l0_idtru ) == 11 ) &&
                 ( std::abs( (*en_).l1_idtru ) == 13 ) ) || 
               ( ( std::abs( (*en_).l0_idtru ) == 13 ) &&
                 ( std::abs( (*en_).l1_idtru ) == 11 ) ) );
    default:
      return false;
  }
}

char EventAttributes::type()
{
  char event_type = 0;
  if (dielectron()) {
    event_type = 1;
  } else if (dimuon()) {
    event_type = 2;
  } else if (elmu()) {
    event_type = 3;
  }
  return event_type;
}

//______________________________________________________________________________

AnalysisSelectors::AnalysisSelectors(DileptonEvents& entry)
{ 
  en_ = &entry;
  l0 = ParticleAttributes( 0, entry );
  l1 = ParticleAttributes( 1, entry );
  asn = EventAttributes( 0, entry );
  tru = EventAttributes( 1, entry );
}

bool AnalysisSelectors::sameSign()
{
  return ((*en_).l0_chrge == (*en_).l1_chrge);
}

bool AnalysisSelectors::oppositeSign()
{
  return !(sameSign());
}

bool AnalysisSelectors::dilepton()
{
  return ( l0.lepton() && l1.lepton() );
}

bool AnalysisSelectors::isSignal()
{
  return ( l0.is_signal() && l1.is_signal() );
}

bool AnalysisSelectors::is_signal_Bs() 
{
  return ( l0.is_signal() && l1.is_signal() &&
      l1.source_Bs() && l0.source_Bs() );
}

bool AnalysisSelectors::is_signal_Bdu() 
{
  return ( l0.is_signal() && l1.is_signal() &&
      l1.source_Bdu() && l0.source_Bdu() );
}

bool AnalysisSelectors::is_correct_wrong() 
{
  return ( ( !(l0.is_signal()) && l1.is_signal()  ) ||
      (l0.is_signal() && !(l1.is_signal())));
}

bool AnalysisSelectors::is_wrong_wrong() 
{
  return (!(l0.is_signal()) && !(l1.is_signal()) && !(is_continuum()));
}

bool AnalysisSelectors::is_continuum()
{
  return (*en_).is_cntnm;
}

char AnalysisSelectors::species_type()
{ 
  return ((asn.type() * 10) + (tru.type()));
}

// Give a unique integer to each type of signal level an event can contain.
int AnalysisSelectors::signal_type()
{
  // Raw data will be categorized as WW so an ID of 0 should never occur outside
  //   an error.
  int signal_type_id = 0;
  
  if (is_signal_Bs()) {
    signal_type_id = 1;
  }
  if (is_signal_Bdu()) {
    signal_type_id = 2;
  }
  if (is_correct_wrong()) {
    signal_type_id = 3;
  }
  if (is_wrong_wrong()) {
    signal_type_id = 4;
  }
  if (is_continuum()) {
    signal_type_id = 5;
  }
  
  return signal_type_id;
}

char AnalysisSelectors::event_sign()
{
  char charge_sum = (*en_).l0_chrge + (*en_).l1_chrge;
  return charge_sum / 2;
}
