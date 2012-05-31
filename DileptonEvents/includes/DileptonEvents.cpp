#include "DileptonEvents.h"

DileptonEvents::DileptonEvents(TTree *tree)
{
  // If parameter tree is not specified (or zero), connect default file or
  //     make a blank one.
  if (tree == 0) {
    std::cout
        << "Attempting to load default file: Adcab.GNMC.basf.root"
        << std::endl;
    TFile *f = new TFile(
        "/home/matt/research/belle/adcab/root/ntuples/Adcab.GNMC.basf.root",
        "READ");
    if (!f || !f->IsOpen()) {
      std::cout
          << "Error accessing root file. Creating dummy ntuple."
          << std::endl;
      f = new TFile("Adcab.GNMC.basf.root");
    }
    f->GetObject("h3", tree);
  }
  init(tree);
}

DileptonEvents::~DileptonEvents()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t DileptonEvents::getEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t DileptonEvents::loadTree(Long64_t entry)
{
  // Set the environment to read one entry.
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    notify();
  }
  return centry;
}

void DileptonEvents::init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  //   a new tree or chain. Typically here the branch addresses and branch
  //   pointers of the tree will be set.

  // Set branch addresses and branch pointers.
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("stm_no", &stm_no, &b_stm_no);
  fChain->SetBranchAddress("exp_no", &exp_no, &b_exp_no);
  fChain->SetBranchAddress("run_no", &run_no, &b_run_no);
  fChain->SetBranchAddress("evt_no", &evt_no, &b_evt_no);
  fChain->SetBranchAddress("is_mc", &is_mc, &b_is_mc);
  fChain->SetBranchAddress("is_cntnm", &is_cntnm, &b_is_cntnm);
  fChain->SetBranchAddress("fw_r2", &fw_r2, &b_fw_r2);
  fChain->SetBranchAddress("hadronb", &hadronb, &b_hadronb);
  fChain->SetBranchAddress("cm_enrgy", &cm_enrgy, &b_cm_enrgy);
  fChain->SetBranchAddress("n_kaons", &n_kaons, &b_n_kaons);
  fChain->SetBranchAddress("typ_asn", &typ_asn, &b_typ_asn);
  fChain->SetBranchAddress("typ_tru", &typ_tru, &b_typ_tru);
  fChain->SetBranchAddress("evt_sign", &evt_sign, &b_evt_sign);
  fChain->SetBranchAddress("cos_thta", &cos_thta, &b_cos_thta);
  fChain->SetBranchAddress("inv_mass", &inv_mass, &b_inv_mass);
  fChain->SetBranchAddress("l0_chrge", &l0_chrge, &b_l0_chrge);
  fChain->SetBranchAddress("l1_chrge", &l1_chrge, &b_l1_chrge);
  fChain->SetBranchAddress("l0_idasn", &l0_idasn, &b_l0_idasn);
  fChain->SetBranchAddress("l1_idasn", &l1_idasn, &b_l1_idasn);
  fChain->SetBranchAddress("l0_idtru", &l0_idtru, &b_l0_idtru);
  fChain->SetBranchAddress("l1_idtru", &l1_idtru, &b_l1_idtru);
  fChain->SetBranchAddress("l0_idmom", &l0_idmom, &b_l0_idmom);
  fChain->SetBranchAddress("l1_idmom", &l1_idmom, &b_l1_idmom);
  fChain->SetBranchAddress("l0_plab", &l0_plab, &b_l0_plab);
  fChain->SetBranchAddress("l1_plab", &l1_plab, &b_l1_plab);
  fChain->SetBranchAddress("l0_pcm", &l0_pcm, &b_l0_pcm);
  fChain->SetBranchAddress("l1_pcm", &l1_pcm, &b_l1_pcm);
  fChain->SetBranchAddress("l0_e_cm", &l0_e_cm, &b_l0_e_cm);
  fChain->SetBranchAddress("l1_e_cm", &l1_e_cm, &b_l1_e_cm);
  fChain->SetBranchAddress("l0_pcm_x", &l0_pcm_x, &b_l0_pcm_x);
  fChain->SetBranchAddress("l1_pcm_x", &l1_pcm_x, &b_l1_pcm_x);
  fChain->SetBranchAddress("l0_pcm_y", &l0_pcm_y, &b_l0_pcm_y);
  fChain->SetBranchAddress("l1_pcm_y", &l1_pcm_y, &b_l1_pcm_y);
  fChain->SetBranchAddress("l0_pcm_z", &l0_pcm_z, &b_l0_pcm_z);
  fChain->SetBranchAddress("l1_pcm_z", &l1_pcm_z, &b_l1_pcm_z);
  fChain->SetBranchAddress("l0_ip_dr", &l0_ip_dr, &b_l0_ip_dr);
  fChain->SetBranchAddress("l1_ip_dr", &l1_ip_dr, &b_l1_ip_dr);
  fChain->SetBranchAddress("l0_ip_dz", &l0_ip_dz, &b_l0_ip_dz);
  fChain->SetBranchAddress("l1_ip_dz", &l1_ip_dz, &b_l1_ip_dz);
  fChain->SetBranchAddress("l0_svdr", &l0_svdr, &b_l0_svdr);
  fChain->SetBranchAddress("l1_svdr", &l1_svdr, &b_l1_svdr);
  fChain->SetBranchAddress("l0_svdz", &l0_svdz, &b_l0_svdz);
  fChain->SetBranchAddress("l1_svdz", &l1_svdz, &b_l1_svdz);
  notify();
}

Bool_t DileptonEvents::notify()
{
  // The Notify() function is called when a new file is opened. This
  //   can be either for a new TTree in a TChain or when when a new TTree
  //   is started when using PROOF. It is normally not necessary to make
  //   changes to the generated code, but the routine can be extended by
  //   the user if needed. The return value is currently not used.

   return kTRUE;
}

void DileptonEvents::show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t DileptonEvents::cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

void DileptonEvents::processNtuple()
{
  // This is the loop skeleton where:
  //     jth_entry is the global entry number in the chain
  //     ith_entry is the entry number in the current Tree
  // Note that the argument to GetEntry must be:
  //   jth_entry for TChain::GetEntry
  //   ith_entry for TTree::GetEntry and TBranch::GetEntry
  //
  // To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2:
  //    replace line:
  //    fChain->GetEntry(jth_entry);       //read all branches
  //    by:
  //    b_branchname->GetEntry(ith_entry); //read only this branch

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jth_entry = 0; jth_entry < nentries; jth_entry++) {
      Long64_t ith_entry = loadTree(jth_entry);
      if (ith_entry < 0) break;
      nb = fChain->GetEntry(jth_entry);   nbytes += nb;
      // if (Cut(ith_entry) < 0) continue;
      ntupleLoopCore();
   }
}

void DileptonEvents::ntupleLoopCore()
{
  // Intentionally blank.
  // This is the set of commands run on each entry in the Ntuple. It is
  //   intended that this is overridden by derived classes.
}