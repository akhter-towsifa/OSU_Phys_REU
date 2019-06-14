#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "TH1F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "iostream"
#include "array"
#include <bits/stdc++.h>
//#include "TLorentzVector"
#include "TMath.h"
#include <math.h>
#include <cmath>
#include "TCanvas.h"

/*This code stores all the jets that have a b particle in it. Then it selects
two of the most energetic b-jets [pt_bjet].*/

void invariant_mass(){
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 100;
  int nbins_1 = 1000;
  int w_2 = 1;

  auto h_jet_pt 	= new TH1F("jet_pt","Jet_P_{T}; X; Counts",		nbins_1,0.0, 2.4e6);
  auto h_jet_pt_btag 	= new TH1I("jet_btagged","Jet p_{T} b-tagged;X",	nbins_1,0.0, 2.4e6);
  auto h_M_inv 		= new TH1F("M_inv", "Invariant Mass; Mass[GeV]; Counts",nbins,0.0, 2.5e2);
  auto file = TFile::Open("bbww_x1000_s170.root");

  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         jet_pt[14];   	//[n_jet]
  Int_t           jet_btagged[14];   	//[n_jet]
  Float_t         jet_eta[14];   	//[n_jet]
  Float_t         jet_phi[14];  	//[n_jet]
  Float_t		w;
  
  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt", jet_pt);
  ctree->SetBranchAddress("jet_btagged",jet_btagged);
  ctree->SetBranchAddress("jet_eta", jet_eta);
  ctree->SetBranchAddress("jet_phi", jet_phi);
  ctree->SetBranchAddress("weight", &w);
  
  ctree->SetBranchStatus("*",0);
  
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  ctree->SetBranchStatus("jet_btagged",1);
  ctree->SetBranchStatus("jet_eta",1);
  ctree->SetBranchStatus("jet_phi",1);
  ctree->SetBranchStatus("weight", 1);

  int nentries = ctree->GetEntries();
  TFile f("mass_invariant.root","recreate");
  

  int64_t bj[14] = {0};			//stores the index of events with b-jet(s).
  float pt_bjet_all[20000][2];		//stores the pt of the events with greater than 2 b-jets in a 2-dimensional array.
  float pt_bjet_2jets[20000][2] = {0};//stores the pt of the events with 2 b-jets with the greatest energy in a 2-dimensional array.
  
  float phi_bjet_all[20000][2] = {0};	//stores the phi angles of all of the jets in a 2-dimensional array
  float eta_bjet_all[20000][2] = {0};  //stores the eta angles of all of the jets in a 2-dimensional array
  
  float phi_bjet_2jets[20000][2] = {0};	//stores the phi angles of the 2 jets in a 2-dimensional array
  float eta_bjet_2jets[20000][2] = {0};	//stores the eta angles of the 2 jets in a 2-dimensional array
 
  float x[20000] = {0};			//stores the (P_x,0 	+ P_x,1) 	values
  float y[20000] = {0};			//stores the (P_y,0 	+ P_y,1) 	values
  float z[20000] = {0};			//stores the (P_z,0 	+ P_z,1) 	values
  float t[20000] = {0};			//stores the (P_total,1 + P_total_2) 	values
  
  float M_inv[20000] = {0};		//stores the invariant mass values
  
//The following loop creates a list pt_bjet_2jets that contains 2 b-jets  

  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    //cout << "pt_bjet_all sorted\t\t" << i << endl;
    int nbtags = 0;
    //cout << "size of pt_bjet_all\t\t" << sizeof(pt_bjet_all) << endl;
    
    for (int j=0; j<n_jet; j++) {
      
      if(jet_btagged[j]>0){         	
       
	h_jet_pt_btag->Fill(jet_pt[j], w_2);	
	
	//cout << "\t\t" << jet_pt[j] << endl;
      
	pt_bjet_all[i][nbtags] = jet_pt[j];
	phi_bjet_all[i][nbtags] = jet_phi[j];
	eta_bjet_all[i][nbtags] = jet_eta[j];
	 nbtags++;
      } // btagging
      
	    
    } // jets in an event
 

    if (nbtags>=2) {
      
      pt_bjet_2jets[i][0] = pt_bjet_all[i][0];
      pt_bjet_2jets[i][1] = pt_bjet_all[i][1];
      
      phi_bjet_2jets[i][0] = phi_bjet_all[i][0];
      phi_bjet_2jets[i][1] = phi_bjet_all[i][1];
      
      eta_bjet_2jets[i][0] = phi_bjet_all[i][0];
      eta_bjet_2jets[i][1] = phi_bjet_all[i][1];
      
      //cout << "Entry: " << i << ";\tn_b-tagged: " << nbtags << ";\tpt1 = " << pt_bjet_2jets[i][0] << ";\tpt2 = " << pt_bjet_2jets[i][1] << endl << endl;
            
    
    
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
	  
    M_inv[i] = (sqrt ( pow(t[i] , 2) - pow(x[i] , 2) - pow(y[i] , 2) - pow(z[i] , 2) )) / 1000;
    h_M_inv->Fill(M_inv[i], w);

    }
    
//this loop creates all the jet_pt histogram
    for(int k = 0; k< n_jet; k++){
      h_jet_pt->Fill( jet_pt[k], w_2);
    }
    
    
  }
  
 
  h_M_inv      ->Write();
  h_jet_pt     ->Write();
  h_jet_pt_btag->Write();
  
  c2 = new TCanvas("c2", "Particle Jets from PP collision", 1200, 400);
  c2 -> Divide (3,1);
  
  c2->cd(1);
  h_jet_pt->Draw();
  
  c2->cd(2);
  h_jet_pt_btag->Draw();
  
  c2->cd(3);
  h_M_inv->Draw();
  
  f.Write();
  f.Close();
}



