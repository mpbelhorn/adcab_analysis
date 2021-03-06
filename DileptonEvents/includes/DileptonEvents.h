#ifndef DileptonEvents_h
#define DileptonEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TCut.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include <iostream>
#include <vector>

class EventTags {
 public:
  EventTags()
  {
    species_.clear();
    signs_.clear();
    components_.clear();
    colors_.clear();
    
    species_.push_back("er");
    species_.push_back("elel");
    species_.push_back("elmu");
    species_.push_back("mumu");
    
    signs_.push_back("nn");
    signs_.push_back("pn");
    signs_.push_back("pp");
    
    components_.push_back("er");
    components_.push_back("bs");
    components_.push_back("bd");
    components_.push_back("cw");
    components_.push_back("ww");
    components_.push_back("cn");
    
     colors_.push_back(1);
     colors_.push_back(91);
     colors_.push_back(2);
     colors_.push_back(8);
     colors_.push_back(4);
     colors_.push_back(7);
    
  }
  ~EventTags() {}
  
  vector<TString> species_;
  vector<TString> signs_;
  vector<TString> components_;
  vector<int> colors_;
};

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
  Float_t n_k_min;
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
  TBranch *b_n_k_min;
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

  DileptonEvents(
      TString input_ntuple_file="default.root",
      TString analysis_name="UnspecifiedAnalysis");
  virtual ~DileptonEvents();
  virtual int  cut(Long64_t entry);
  virtual Int_t  getEntry(Long64_t entry);
  virtual Long64_t loadTree(Long64_t entry);
  virtual void  init(TTree *tree);
  virtual void  processNtuple();
  // TODO: replace argument below with something like GetIndex for 
  //   derived classes that need it.
  virtual void  ntupleLoopCore(const int& entry_id = 0);
  virtual bool  notify();
  virtual void  show(Long64_t entry = -1);
  
  TCut pp_events_cut_;
  TCut pn_events_cut_;
  TCut nn_events_cut_;
  TCut bs_events_cut_;
  TCut bd_events_cut_;
  TCut cw_events_cut_;
  TCut ww_events_cut_;
  TCut cn_events_cut_;
  
  void setCreateDataSet(const bool& yes_or_no);
  void recreateDataSet(RooArgSet& variables);
  void saveDataSet(const TString& filename = "data.root");
  
  TString input_ntuple_file_;
  TString analysis_name_;
  bool flag_output_a_dataset_;
  
  RooCategory* component_;
  RooCategory* event_sign_;
  RooCategory* event_species_;
  RooDataSet* data_set_;
  
  const EventTags tags_;
  
};

#endif
