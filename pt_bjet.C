#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TH1F.h"
#include "TH1I.h"
#include "TROOT.h"
#include "TStyle.h"
#include "iostream"
#include "array"
#include <bits/stdc++.h>
//#include "TLorentzVector"
#include "TMath.h"
#include <math.h>
#include <cmath>

/*This code stores all the jets that have a b particle in it. Then it selects
two of the most energetic b-jets [pt_bjet].*/

void pt_bjet(){
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 1000;
  double w = 1;
  
  auto h_jet_pt_btag 	= new TH1I("jet_btagged","Jet p_{T} b-tagged;X",	nbins,0.0, 2.4e6);
  auto h_jet_pt 	= new TH1F("jet_pt","Jet_P_{T}; X; Counts",		nbins,0.0, 2.4e6);
  auto h_M_inv 		= new TH1F("M_inv", "Invariant Mass; Mass[GeV]; Counts",nbins,0.0, 1.0e7);
  auto file = TFile::Open("bbww_x1000_s170.root");

  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         jet_pt[14];   	//[n_jet]
  Int_t           jet_btagged[14];   	//[n_jet]
  Float_t         jet_eta[14];   	//[n_jet]
  Float_t         jet_phi[14];  	//[n_jet]
  
  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt", jet_pt);
  ctree->SetBranchAddress("jet_btagged",jet_btagged);
  ctree->SetBranchAddress("jet_eta", jet_eta);
  ctree->SetBranchAddress("jet_phi", jet_phi);

  ctree->SetBranchStatus("*",0);
  
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  ctree->SetBranchStatus("jet_btagged",1);
  ctree->SetBranchStatus("jet_eta",1);
  ctree->SetBranchStatus("jet_phi",1);


  int nentries = ctree->GetEntries();
  TFile f("mass_invariant.root","recreate");
  

  int64_t bj[14] = {};			//stores the index of events with b-jet(s).
  int pt_bjet_all[20000][2];		//stores the pt of the events with greater than 2 b-jets.
  int64_t pt_bjet_2jets[20000][1] = {};	//stores the pt of the events with 2 b-jets with the greatest energy.
  
  float phi_bjet_all[20000][1] = {};	//stores the phi angles of all of the jets
  float eta_bjet_all[20000][1] = {};  	//stores the eta angles of all of the jets
  
  float phi_bjet_2jets[20000][1] = {};	//stores the phi angles of the 2 jets
  float eta_bjet_2jets[20000][1] = {};	//stores the eta angles of the 2 jets 
  
  float x[20000] = {};			//stores the (P_x,0 	+ P_x,1) 	values
  float y[20000] = {};			//stores the (P_y,0 	+ P_y,1) 	values
  float z[20000] = {};			//stores the (P_z,0 	+ P_z,1) 	values
  float t[20000] = {};			//stores the (P_total,1 + P_total_2) 	values
  
  float M_inv[20000] = {};		//stores the invariant mass values
  
//The following loop creates a list pt_bjet_2jets that contains 2 b-jets  

  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    //cout << "pt_bjet_all sorted\t\t" << i << endl;
    int nbtags = 0;
    //cout << "size of pt_bjet_all\t\t" << sizeof(pt_bjet_all) << endl;
    
    for (int j=0; j<n_jet; j++) {
      
      if(jet_btagged[j]>0){         
	bj[nbtags]=j;
        nbtags++;
	h_jet_pt_btag->Fill(jet_pt[j], w);	
      }
      pt_bjet_all[i][j] = jet_pt[j];
      phi_bjet_all[i][j] = jet_phi[j];
      eta_bjet_all[i][j] = jet_eta[j];
      
//      int n = sizeof(pt_bjet_all)/sizeof(pt_bjet_all[0]);
//      sort(pt_bjet_all, pt_bjet_all+n, greater<int>());
//      cout << "\t\t\t\t\t" << pt_bjet_all[i][j] << endl;
	    
    }


    if (sizeof(pt_bjet_all) >= 2) {
      pt_bjet_2jets[i][0] = pt_bjet_all[i][0];
      pt_bjet_2jets[i][1] = pt_bjet_all[i][1];

      phi_bjet_2jets[i][0] = phi_bjet_all[i][0];
      phi_bjet_2jets[i][1] = phi_bjet_all[i][1];
      
      eta_bjet_2jets[i][0] = phi_bjet_all[i][0];
      eta_bjet_2jets[i][1] = phi_bjet_all[i][1];
      
      //cout << "pt_bjet_2jets\t" << i << "\t" << pt_bjet_2jets[i][0] << "\t" << phi_bjet_2jets[i][0] << "\t" << eta_bjet_2jets[i][0] << endl;
      //cout << "pt_bjet_2jets\t" << " " << "\t" << pt_bjet_2jets[i][1] << "\t" << phi_bjet_2jets[i][1] << "\t" << eta_bjet_2jets[i][1] << endl;
      
    }
    
    x[i] = pt_bjet_2jets[i][0] * cos(phi_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * cos(phi_bjet_2jets[i][1]);
    y[i] = pt_bjet_2jets[i][0] * sin(phi_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * sin(phi_bjet_2jets[i][1]);
    z[i] = pt_bjet_2jets[i][0] * sinh(eta_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * sinh(eta_bjet_2jets[i][1]);
    t[i] = sqrt( pow(pt_bjet_2jets[i][0] * cos(phi_bjet_2jets[i][0]) , 2)
		+pow(pt_bjet_2jets[i][0] * sin(phi_bjet_2jets[i][0]) , 2)
		+pow(pt_bjet_2jets[i][0] * sinh(eta_bjet_2jets[i][0]) , 2)
	   )
	  +sqrt( pow(pt_bjet_2jets[i][1] * cos(phi_bjet_2jets[i][1]) , 2)
		+pow(pt_bjet_2jets[i][1] * sin(phi_bjet_2jets[i][1]) , 2)
		+pow(pt_bjet_2jets[i][1] * sinh(eta_bjet_2jets[i][1]) , 2)
	   );
	  
    M_inv[i] = sqrt ( pow(t[i] , 2) - pow(x[i] , 2) - pow(y[i] , 2) - pow(z[i] , 2) );
    h_M_inv->Fill(M_inv[i], w);
    
	  
  }
  
  h_M_inv->Write();
  
  f.Write();
  f.Close();
}



