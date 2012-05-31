//
//******************************************************************************
// Filename: ParticleIdCodes.h
// Version: 2011.03.25.A
// Author: M.P. Belhorn
// Description: A global Pythia PID codes map.
//******************************************************************************

#ifndef PARTICLEIDCODES_H
#define PARTICLEIDCODES_H
#include <map>

// A structure of Pythia Particle ID Codes.
// Antiparticles (if they exist) are opposite (negative) sign of the related
//   particle code.
struct PythiaPidCodes {
  static const int e            = 11;
  static const int mu           = 13;
  static const int gamma        = 22;
  static const int v_gamma      = 10022; // Virtual photon.
  
  static const int pi0          = 111;
  static const int pi_plus      = 211;
  
  static const int eta          = 221;
  static const int eta_prime    = 331;
  
  static const int rho0         = 113;
  static const int rho_plus     = 213;
  static const int omega        = 223;
  static const int phi          = 333;
  
  static const int K0           = 311;
  static const int K0_star      = 313;
  static const int K_plus       = 321;
  static const int K_plus_star  = 323;
  
  static const int D_plus       = 411;
  static const int D_plus_star  = 413;
  static const int D0           = 421;
  static const int D0_star      = 423;
  static const int Ds_plus      = 431;
  static const int Ds_plus_star = 433;
  static const int J_psi        = 443;
  
  static const int B0           = 511;
  static const int B0_star      = 513;
  static const int B_plus       = 521;
  static const int B_plus_star  = 523;
  static const int Bs0          = 531;
  static const int Bs0_star     = 533;
};


#endif
