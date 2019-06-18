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

void ttbar_2() {
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 100;
  int w = 1;
  
//   h1 = new TH1F("h1", "Invariant Mass", nbins, 0, 600);
  
  auto h_M_inv = new TH1F("M_inv", "Invariant Mass; Mass[GeV]; Counts",nbins,0.0, 6.0e5);
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
  
  int nentries = nominalTree->GetEntry();
  
  cout << "Number of entries: " << nentries << endl;
  return;
  
  float x[] = {0};
  float y[] = {0};
  float z[] = {0};
  float m[] = {};
  float M_inv[] = {};
  
//   float pt_bjet_all[981137][2] = {0};		//stores all the pt bjets with more than 1 btags
//   float pt_bjet_2jets[981137][2] = {0};		//stores the 2 highest energy bjets 
//   
//   float phi_bjet_all[981137][2] = {0};		//stores phi all the pt bjets with more than 1 btags
//   float phi_bjet_2jets[981137][2] = {0};		//stores phi for the 2 highest energy bjets
//   
//   float eta_bjet_all[981137][2] = {0};		//stores eta for all the pt bjets with more than 1 btags
//   float eta_bjet_2jets[981137][2] = {0};		//stores eta for the 2 highest energy bjets
//   
  
  for (int i=0; i<nentries; i++) {
    nominalTree->GetEntry(i);
  
    int n_bjets=0;
    int n_jets = jet_pt->size();
    
    cout << "\nNumber of jets: " << n_jets << endl;
    
    for (int j=0; j<n_jets; j++) {
            
      if ((int)(*jet_btagged)[j]==1){

// 	pt_bjet_all[i][n_bjets] = (float)(*jet_pt)[j];
// 	phi_bjet_all[i][n_bjets] = (float)(*jet_phi)[j];
// 	eta_bjet_all[i][n_bjets] = (float)(*jet_eta)[j];
	
	n_bjets++;
      }

    }
    cout << "Event(" << i << "), n_bjets = " << n_bjets << endl;

//     if (n_bjets>=2)
//     {
//       pt_bjet_2jets[i][0] = pt_bjet_all[i][0];
//       pt_bjet_2jets[i][1] = pt_bjet_all[i][1];
//       
//       phi_bjet_2jets[i][0] = phi_bjet_all[i][0];
//       phi_bjet_2jets[i][1] = phi_bjet_all[i][1];
//       
//       eta_bjet_2jets[i][0] = phi_bjet_all[i][0];
//       eta_bjet_2jets[i][1] = phi_bjet_all[i][1];
//       
//       x[i] = pt_bjet_2jets[i][0] * cos(phi_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * cos(phi_bjet_2jets[i][1]);
//       y[i] = pt_bjet_2jets[i][0] * sin(phi_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * sin(phi_bjet_2jets[i][1]);
//       z[i] = pt_bjet_2jets[i][0] * sinh(eta_bjet_2jets[i][0]) + pt_bjet_2jets[i][1] * sinh(eta_bjet_2jets[i][1]);
//       m[i] = sqrt( pow(pt_bjet_2jets[i][0] * cos(phi_bjet_2jets[i][0]) , 2)
// 		+pow(pt_bjet_2jets[i][0] * sin(phi_bjet_2jets[i][0]) , 2)
// 		+pow(pt_bjet_2jets[i][0] * sinh(eta_bjet_2jets[i][0]) , 2)
// 	   )
// 	  +sqrt( pow(pt_bjet_2jets[i][1] * cos(phi_bjet_2jets[i][1]) , 2)
// 		+pow(pt_bjet_2jets[i][1] * sin(phi_bjet_2jets[i][1]) , 2)
// 		+pow(pt_bjet_2jets[i][1] * sinh(eta_bjet_2jets[i][1]) , 2)
// 	   );
//       
//       M_inv[i] = (sqrt ( pow(m[i] , 2) - pow(x[i] , 2) - pow(y[i] , 2) - pow(z[i] , 2) ));
//       h_M_inv->Fill(M_inv[i], w);
// //       h1->Fill( M_inv[i]);
//       
//     }
//     
  }
  
//   c1 = new TCanvas("h1", "Invariant Mass", 1200, 400);
//   c1->Divide(3,1);
//   c1->cd(1);
//   h1->Draw();
  
  file->Close();
  
  TFile f("ttbar_invariant_mass_2.root", "RECREATE");
  h_M_inv ->Write();

  
  //f.Write();
  //f.Close(); 
  
}