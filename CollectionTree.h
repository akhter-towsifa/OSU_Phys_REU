//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun  4 13:09:13 2019 by ROOT version 6.06/06
// from TTree CollectionTree/
// found on file: bbww_x1000_s170.root
//////////////////////////////////////////////////////////

#ifndef CollectionTree_h
#define CollectionTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class CollectionTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           n_jet;
   Float_t         jet_pt[14];   //[n_jet]
   Float_t         jet_eta[14];   //[n_jet]
   Float_t         jet_phi[14];   //[n_jet]
   Float_t         jet_e[14];   //[n_jet]
   Float_t         jet_Jvt[14];   //[n_jet]
   Int_t           jet_btagged[14];   //[n_jet]
   Int_t           n_electrons;
   Float_t         electrons_pt[1];   //[n_electrons]
   Float_t         electrons_eta[1];   //[n_electrons]
   Float_t         electrons_phi[1];   //[n_electrons]
   Float_t         electrons_e[1];   //[n_electrons]
   Float_t         electrons_etcone20[1];   //[n_electrons]
   Float_t         electrons_ptcone30[1];   //[n_electrons]
   Int_t           n_muons;
   Float_t         muons_pt[1];   //[n_muons]
   Float_t         muons_eta[1];   //[n_muons]
   Float_t         muons_phi[1];   //[n_muons]
   Float_t         muons_e[1];   //[n_muons]
   Float_t         muons_etcone20[1];   //[n_muons]
   Float_t         muons_ptcone30[1];   //[n_muons]
   Float_t         met;
   Float_t         met_x;
   Float_t         met_y;
   Float_t         met_phi;
   Float_t         met_sumet;
   Float_t         mc_event_weight;
   Float_t         weight;
   Float_t         av_int_per_xing;
   Float_t         num_pv;
   Int_t           run_number;
   Int_t           lumi_block;
   Int_t           mc_channel_number;
   Int_t           event_number;

   // List of branches
   TBranch        *b_n_jet;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_e;   //!
   TBranch        *b_jet_Jvt;   //!
   TBranch        *b_jet_btagged;   //!
   TBranch        *b_n_electrons;   //!
   TBranch        *b_electrons_pt;   //!
   TBranch        *b_electrons_eta;   //!
   TBranch        *b_electrons_phi;   //!
   TBranch        *b_electrons_e;   //!
   TBranch        *b_electrons_etcone20;   //!
   TBranch        *b_electrons_ptcone30;   //!
   TBranch        *b_n_muons;   //!
   TBranch        *b_muons_pt;   //!
   TBranch        *b_muons_eta;   //!
   TBranch        *b_muons_phi;   //!
   TBranch        *b_muons_e;   //!
   TBranch        *b_muons_etcone20;   //!
   TBranch        *b_muons_ptcone30;   //!
   TBranch        *b_met;   //!
   TBranch        *b_met_x;   //!
   TBranch        *b_met_y;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_sumet;   //!
   TBranch        *b_mc_event_weight;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_av_int_per_xing;   //!
   TBranch        *b_num_pv;   //!
   TBranch        *b_run_number;   //!
   TBranch        *b_lumi_block;   //!
   TBranch        *b_mc_channel_number;   //!
   TBranch        *b_event_number;   //!

   CollectionTree(TTree *tree=0);
   virtual ~CollectionTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CollectionTree_cxx
CollectionTree::CollectionTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("bbww_x1000_s170.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("bbww_x1000_s170.root");
      }
      f->GetObject("CollectionTree",tree);

   }
   Init(tree);
}

CollectionTree::~CollectionTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CollectionTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CollectionTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CollectionTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("n_jet", &n_jet, &b_n_jet);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_e", jet_e, &b_jet_e);
   fChain->SetBranchAddress("jet_Jvt", jet_Jvt, &b_jet_Jvt);
   fChain->SetBranchAddress("jet_btagged", jet_btagged, &b_jet_btagged);
   fChain->SetBranchAddress("n_electrons", &n_electrons, &b_n_electrons);
   fChain->SetBranchAddress("electrons_pt", &electrons_pt, &b_electrons_pt);
   fChain->SetBranchAddress("electrons_eta", &electrons_eta, &b_electrons_eta); 
   fChain->SetBranchAddress("electrons_phi", &electrons_phi, &b_electrons_phi);
   fChain->SetBranchAddress("electrons_e", &electrons_e, &b_electrons_e);
   fChain->SetBranchAddress("electrons_etcone20", &electrons_etcone20, &b_electrons_etcone20);
   fChain->SetBranchAddress("electrons_ptcone30", &electrons_ptcone30, &b_electrons_ptcone30);
   fChain->SetBranchAddress("n_muons", &n_muons, &b_n_muons);
   fChain->SetBranchAddress("muons_pt", &muons_pt, &b_muons_pt);
   fChain->SetBranchAddress("muons_eta", &muons_eta, &b_muons_eta);
   fChain->SetBranchAddress("muons_phi", &muons_phi, &b_muons_phi);
   fChain->SetBranchAddress("muons_e", &muons_e, &b_muons_e);
   fChain->SetBranchAddress("muons_etcone20", &muons_etcone20, &b_muons_etcone20);
   fChain->SetBranchAddress("muons_ptcone30", &muons_ptcone30, &b_muons_ptcone30);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("met_x", &met_x, &b_met_x);
   fChain->SetBranchAddress("met_y", &met_y, &b_met_y);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_sumet", &met_sumet, &b_met_sumet);
   fChain->SetBranchAddress("mc_event_weight", &mc_event_weight, &b_mc_event_weight);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("av_int_per_xing", &av_int_per_xing, &b_av_int_per_xing);
   fChain->SetBranchAddress("num_pv", &num_pv, &b_num_pv);
   fChain->SetBranchAddress("run_number", &run_number, &b_run_number);
   fChain->SetBranchAddress("lumi_block", &lumi_block, &b_lumi_block);
   fChain->SetBranchAddress("mc_channel_number", &mc_channel_number, &b_mc_channel_number);
   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   Notify();
}

Bool_t CollectionTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CollectionTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CollectionTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CollectionTree_cxx
