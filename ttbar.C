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

void ttbar() {
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 100;
  int nbins_1 = 1000;
  int w = 1;
  
  auto h_jet_pt 	= new TH1F("jet_pt","Jet_P_{T}; X; Counts",		nbins_1,0.0, 2.4e6);
  auto h_jet_pt_btag 	= new TH1I("jet_btagged","Jet p_{T} b-tagged;X",	nbins_1,0.0, 2.4e6);
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
  TFile f("ttbar_invariant_mass.root", "recreate");
  
  float x[] = {};
  float y[] = {};
  float z[] = {};
  float m[] = {};
  
  cout << "Number of entries: " << nentries << endl;
  
  for (int i=0; i<nentries; i++) {
    nominalTree->GetEntry(i);

    int n_bjets=0;
    int n_jets = jet_pt->size();
    
    //cout << "\nNumber of jets: " << n_jets << endl;
    
    for (int j=0; j<n_jets; j++) {
            
      if ((int)(*jet_btagged)[j]==1){
	//cout << "\tB-tagging: " << (int)(*jet_btagged)[0] << endl;
	//float pt_bjet_1 = jet_pt[i];
	n_bjets++;
      }
      //cout << "Event(" << i << "), n_bjets = " << n_bjets << endl;
    }
//     
//     if (jet_pt[i] >= 420 && abs(jet_eta[i]) <= 2){
//       float pt_bjet_2 = jet_pt[i];
//       n_bjets++;
//     }
//     
//     if (n_bjets == 2){
//       x[i] = pt_bjet_1 * cos(jet_phi[i]) + pt_bjet_2 * cos(jet_phi[i]) ;
//       y[i] = pt_bjet_1 * sin(jet_phi[i]) + pt_bjet_2 * sin(jet_phi[i]) ;
//       z[i] = pt_bjet_1 * sinh(jet_eta[i]) + pt_bjet_2 * sinh(jet_eta[i]) ;
//       m[i] = sqrt(pow(pt_bjet_1 * cos(jet_phi[i]) ,2)
// 		  pow(pt_bjet_1 * sin(jet_phi[i]) ,2)
// 		  pow(pt_bjet_1 * sinh(jet_eta[i]) ,2)
// 	     )
// 	    +sqrt(pow(pt_bjet_2 * cos(jet_phi[i]) ,2)
// 		  pow(pt_bjet_2 * sin(jet_phi[i]) ,2)
// 		  pow(pt_bjet_2 * sinh(jet_eta[i]) ,2)
// 	     );
//       M_inv[i] = (sqrt ( pow(m[i] , 2) - pow(x[i] , 2) - pow(y[i] , 2) - pow(z[i] , 2) ));
//       h_M_inv->Fill(M_inv[i], w);
//     }
  
    
  }
  
//  h_M_inv ->Write();

  
  f.Write();
  f.Close(); 
  
}