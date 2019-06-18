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

//this code graphs the invariant mass of the bb-bar from the tt-bar.

void ttbar_1() {
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 100;
  int w = 1;
  
  auto h_M_inv 		= new TH1F("M_inv", "Invariant Mass; Mass[GeV]; Counts",nbins,0.0, 6.0e2);
  auto file = TFile::Open("ttbar.root");
  
  TTree* nominalTree = (TTree*)file->Get("nominal");
  
  vector<float>   *jet_pt = 0;
  vector<float>   *jet_eta = 0;
  vector<float>   *jet_phi = 0;
  vector<char>    *jet_btagged;
  /*vector<float>   *ljet_pt;
  vector<float>   *ljet_eta;
  vector<float>   *ljet_phi;
  vector<float>   *ljet_m;
  vector<float>   *tjet_pt;
  vector<float>   *tjet_eta;
  vector<float>   *tjet_phi;
  vector<char>    *tjet_btagged;*/
  
  nominalTree->SetBranchAddress("jet_pt", &jet_pt);
  nominalTree->SetBranchAddress("jet_eta", &jet_eta);
  nominalTree->SetBranchAddress("jet_phi", &jet_phi);
  nominalTree->SetBranchAddress("jet_isbtagged_MV2c10_85", &jet_btagged);
  /*nominalTree->SetBranchAddress("ljet_pt", &ljet_pt);
  nominalTree->SetBranchAddress("ljet_eta", &ljet_eta);
  nominalTree->SetBranchAddress("ljet_phi", &ljet_phi);
  nominalTree->SetBranchAddress("ljet_m", &ljet_m);
  nominalTree->SetBranchAddress("tjet_pt", &tjet_pt);
  nominalTree->SetBranchAddress("tjet_eta", &tjet_eta);
  nominalTree->SetBranchAddress("tjet_phi", &tjet_phi);
  nominalTree->SetBranchAddress("tjet_isbtagged_MV2c10_85", &tjet_btagged);*/
  
  nominalTree->SetBranchStatus("*",0);
  
  nominalTree->SetBranchStatus("jet_pt",1);
  nominalTree->SetBranchStatus("jet_eta",1);
  nominalTree->SetBranchStatus("jet_phi",1);
  nominalTree->SetBranchStatus("jet_isbtagged_MV2c10_85",1);
  /*nominalTree->SetBranchStatus("ljet_pt",1);
  nominalTree->SetBranchStatus("ljet_eta",1);
  nominalTree->SetBranchStatus("ljet_phi",1);
  nominalTree->SetBranchStatus("ljet_m",1);
  nominalTree->SetBranchStatus("tjet_pt",1);
  nominalTree->SetBranchStatus("tjet_eta",1);
  nominalTree->SetBranchStatus("tjet_phi",1);
  nominalTree->SetBranchStatus("tjet_isbtagged_MV2c10_85",1);*/
  
  int nentries = nominalTree->GetEntries();
  TFile f("ttbar_invariant_mass_1.root", "recreate");
  
  cout << "beginning of the array definition" << endl;
  
  float pt_bjet_all[2000000][2];		//stores the pt of the events with greater than 2 b-jets in a 2-dimensional array.
  float pt_bjet_2jets[2000000][2] = {0};//stores the pt of the events with 2 b-jets with the greatest energy in a 2-dimensional array.
  
  float phi_bjet_all[2000000][2] = {0};	//stores the phi angles of all of the jets in a 2-dimensional array
  float eta_bjet_all[2000000][2] = {0};  //stores the eta angles of all of the jets in a 2-dimensional array
  
  float phi_bjet_2jets[2000000][2] = {0};	//stores the phi angles of the 2 jets in a 2-dimensional array
  float eta_bjet_2jets[2000000][2] = {0};	//stores the eta angles of the 2 jets in a 2-dimensional array
  cout << "end of array definitions\t" << endl;
  
  float x[2000000] = {0};
  float y[2000000] = {0};
  float z[2000000] = {0};
  float m[2000000] = {0};
  float M_inv[2000000] = {0};
  
  cout << "Number of entries: " << nentries << endl;

  for (int i=0; i<nentries; i++) {
    nominalTree->GetEntry(i);
    cout << "we're in the loop\t" << i << endl;
    int n_bjets=0;
    int n_jets = jet_pt->size();
    float pt_bjet_1 = 0;
    float pt_bjet_2 = 0;

    
    for (int j=0; j<n_jets; j++) {
            
      if ((int)(*jet_btagged)[j]==1){
	cout << "\tB-tagging: " << (int)(*jet_btagged)[j] << endl;
	
	pt_bjet_all[i][n_bjets] = (float)(*jet_pt)[j];
	phi_bjet_all[i][n_bjets] = (float)(*jet_phi)[j];
	eta_bjet_all[i][n_bjets] = (float)(*jet_eta)[j];
	
	n_bjets++;
      }
      cout << "Event(" << i << "), n_bjets = " << n_bjets << endl;
    }
    
    if (n_bjets>=2) {
      
      pt_bjet_2jets[i][0] = pt_bjet_all[i][0];
      pt_bjet_2jets[i][1] = pt_bjet_all[i][1];
      
      phi_bjet_2jets[i][0] = phi_bjet_all[i][0];
      phi_bjet_2jets[i][1] = phi_bjet_all[i][1];
      
      eta_bjet_2jets[i][0] = phi_bjet_all[i][0];
      eta_bjet_2jets[i][1] = phi_bjet_all[i][1];
      


    
    cout << "Event(" << i << "), n_bjets = " << n_bjets << endl;
    
    x[i] = pt_bjet_2jets[i][0] * cos(phi_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * cos(phi_bjet_2jets[i][1]);
    y[i] = pt_bjet_2jets[i][0] * sin(phi_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * sin(phi_bjet_2jets[i][1]);
    z[i] = pt_bjet_2jets[i][0] * sinh(eta_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * sinh(eta_bjet_2jets[i][1]);
    m[i] = sqrt( pow(pt_bjet_2jets[i][0] * cos(phi_bjet_2jets[i][0]) , 2)
		+pow(pt_bjet_2jets[i][0] * sin(phi_bjet_2jets[i][0]) , 2)
		+pow(pt_bjet_2jets[i][0] * sinh(eta_bjet_2jets[i][0]) , 2)
	   )
	  +sqrt( pow(pt_bjet_2jets[i][1] * cos(phi_bjet_2jets[i][1]) , 2)
		+pow(pt_bjet_2jets[i][1] * sin(phi_bjet_2jets[i][1]) , 2)
		+pow(pt_bjet_2jets[i][1] * sinh(eta_bjet_2jets[i][1]) , 2)
	   );
      
    M_inv[i] = (sqrt ( pow(m[i] , 2) - pow(x[i] , 2) - pow(y[i] , 2) - pow(z[i] , 2) ));
    h_M_inv->Fill(M_inv[i], w);
    //cout << "Event(" << i << "), n_bjets = " << n_bjets << endl;
    }    
  }
  
   h_M_inv ->Write();

  
  f.Write();
  f.Close(); 
  
}