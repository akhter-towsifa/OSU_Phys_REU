#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TH1F.h"
#include "TH1I.h"
#include "TROOT.h"
#include "TStyle.h"
#include "iostream"


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
  int pt_bjet_all[14] = {};	//stores the pt of the events with greater than 2 b-jets.
  int64_t pt_bjet_2jets[14] = {};//stores the pt of the events with 2 b-jets with the greatest energy.

// This part creates the jet_pt histogram of all the b-jets.  
  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    int nbtags = 0; 
    
    for (int64_t j=0; j<n_jet; j++) {
      if(jet_btagged[j]>0){         
	bj[nbtags]=j;
        nbtags++;
	h_jet_pt_btag->Fill(jet_pt[j], w);	
      }
      pt_bjet_all[i] = jet_pt[j];

	
      int n = sizeof(pt_bjet_all)/sizeof(pt_bjet_all[0]);
      //cout << "pt_bjet_all\t\t" << pt_bjet_all[i] << endl;
      sort(pt_bjet_all, pt_bjet_all+n, greater<int>());
    }
    //cout << "sorted pt_bjet_all\t\t" << i << "\t\t" << pt_bjet_all[i] << endl;
    
    int size_pt_bjet_all = sizeof(pt_bjet_all)/sizeof(pt_bjet_all[0]);
    //cout << "\t\tsize of pt_bjet\t\t" << size_pt_bjet_all << endl;
    
    
    
    if (size_pt_bjet_all == 2){
      pt_bjet_2jets[i] = pt_bjet_all[i];
      cout << "pt bjet 2 jets\t\t" << pt_bjet_2jets[i] << endl;
    }
    else if (size_pt_bjet_all > 2){
      for (int k=0; k< 2 ; k++){
	pt_bjet_2jets[i] = pt_bjet_all[k];
      }
      //cout << "pt bjet 2 jets\t\t\t" << pt_bjet_2jets[i] << endl;
    }
    
  }

  
  
  
  f.Write();
  f.Close();
}