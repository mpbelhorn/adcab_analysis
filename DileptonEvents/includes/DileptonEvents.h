#ifndef DileptonEvents_h
#define DileptonEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

class DileptonEvents {
 public :
  TTree *fChain;  // Pointer to the analysis TTree or TChain.
  Int_t fCurrent; // Current Tree number in a TChain.

  // Declaration of N-tuple column variables.
  Float_t stm_no;
  Float_t exp_no;
  Float_t run_no;
  Float_t evt_no;
  Float_t is_mc;
  Float_t is_cntnm;
  Float_t fw_r2;
  Float_t hadronb;
  Float_t cm_enrgy;
  Float_t n_kaons;
  Float_t typ_asn;
  Float_t typ_tru;
  Float_t evt_sign;
  Float_t cos_thta;
  Float_t inv_mass;
  Float_t l0_chrge;
  Float_t l1_chrge;
  Float_t l0_idasn;
  Float_t l1_idasn;
  Float_t l0_idtru;
  Float_t l1_idtru;
  Float_t l0_idmom;
  Float_t l1_idmom;
  Float_t l0_plab;
  Float_t l1_plab;
  Float_t l0_pcm;
  Float_t l1_pcm;
  Float_t l0_e_cm;
  Float_t l1_e_cm;
  Float_t l0_pcm_x;
  Float_t l1_pcm_x;
  Float_t l0_pcm_y;
  Float_t l1_pcm_y;
  Float_t l0_pcm_z;
  Float_t l1_pcm_z;
  Float_t l0_ip_dr;
  Float_t l1_ip_dr;
  Float_t l0_ip_dz;
  Float_t l1_ip_dz;
  Float_t l0_svdr;
  Float_t l1_svdr;
  Float_t l0_svdz;
  Float_t l1_svdz;

  // List of branches.
  TBranch *b_stm_no;
  TBranch *b_exp_no;
  TBranch *b_run_no;
  TBranch *b_evt_no;
  TBranch *b_is_mc;
  TBranch *b_is_cntnm;
  TBranch *b_fw_r2;
  TBranch *b_hadronb;
  TBranch *b_cm_enrgy;
  TBranch *b_n_kaons;
  TBranch *b_typ_asn;
  TBranch *b_typ_tru;
  TBranch *b_evt_sign;
  TBranch *b_cos_thta;
  TBranch *b_inv_mass;
  TBranch *b_l0_chrge;
  TBranch *b_l1_chrge;
  TBranch *b_l0_idasn;
  TBranch *b_l1_idasn;
  TBranch *b_l0_idtru;
  TBranch *b_l1_idtru;
  TBranch *b_l0_idmom;
  TBranch *b_l1_idmom;
  TBranch *b_l0_plab;
  TBranch *b_l1_plab;
  TBranch *b_l0_pcm;
  TBranch *b_l1_pcm;
  TBranch *b_l0_e_cm;
  TBranch *b_l1_e_cm;
  TBranch *b_l0_pcm_x;
  TBranch *b_l1_pcm_x;
  TBranch *b_l0_pcm_y;
  TBranch *b_l1_pcm_y;
  TBranch *b_l0_pcm_z;
  TBranch *b_l1_pcm_z;
  TBranch *b_l0_ip_dr;
  TBranch *b_l1_ip_dr;
  TBranch *b_l0_ip_dz;
  TBranch *b_l1_ip_dz;
  TBranch *b_l0_svdr;
  TBranch *b_l1_svdr;
  TBranch *b_l0_svdz;
  TBranch *b_l1_svdz;

   DileptonEvents(TTree *tree=0);
   virtual ~DileptonEvents();
   virtual Int_t    cut(Long64_t entry);
   virtual Int_t    getEntry(Long64_t entry);
   virtual Long64_t loadTree(Long64_t entry);
   virtual void     init(TTree *tree);
   virtual void     processNtuple();
   virtual void     ntupleLoopCore();
   virtual Bool_t   notify();
   virtual void     show(Long64_t entry = -1);
};

#endif
