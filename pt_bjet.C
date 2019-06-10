#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TH1F.h"
#include "TH1I.h"


/*This code stores all the jets that have a b particle in it. Then it selects
two of the most energetic b-jets [pt_bjet].*/

void pt_bjet(){
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 1000;
  double w = 1;
  
  auto h_jet_pt_btag = new TH1I("jet_btagged","Jet p_{T} b-tagged;X",nbins,0.0, 2.4e6);
  auto h_jet_pt = new TH1F("jet_pt","Jet_P_{T}; X; Counts",nbins,0.0,2.4e6);
  auto file = TFile::Open("bbww_x1000_s170.root");
  
  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         jet_pt[14];   //[n_jet]
  Int_t           jet_btagged[14];   //[n_jet]
  
  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt", jet_pt);
  ctree->SetBranchAddress("jet_btagged",jet_btagged);
  ctree->SetBranchStatus("*",0);
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  ctree->SetBranchStatus("jet_btagged",1);

  int nentries = ctree->GetEntries();
  TFile f("pt_bjet.root","recreate");
  
// This part creates the jet_pt histogram of all jets.    
  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    for(uint j = 0; j< n_jet; j++){
      h_jet_pt->Fill( jet_pt[j],w ); 
    }
  }

  int64_t bj[14] = {};		//stores the index of events with b-jet(s).
  int64_t pt_bjet_all = {};		//stores the pt of the events with greater than 2 b-jets.
  int64_t pt_bjet_2jets = {};		//stores the pt of the events with 2 b-jets with the greatest energy.
// This part creates the jet_pt histogram of all the b-jets.  
  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    int nbtags = 0; 
    for (int64_t j=0; j<n_jet; j++) {
      if(jet_btagged[j]>0){         
	bj[nbtags]=j;
        nbtags++;
	h_jet_pt_btag->Fill(jet_pt[j], w);
	cout << "\t pt_bjet \t" << jet_pt[j] << "\tj\t" << j << endl;
	
	if (j >=2) {
	  pt_bjet = jet_pt[j];
	  if (j > 2) {
	    for (int k=0; k<j+1; j++){
	      if (pt_bjet[k]
	    }
	  }
	  else {
	    pt_bjet_2jets = pt_bjet
	}
      }
    }
  }
  
  
  
  f.Write();
  f.Close();
}