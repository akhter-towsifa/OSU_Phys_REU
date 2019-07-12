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
#include "TMath.h"
#include <math.h>
#include <cmath>
#include "TCanvas.h"
#include "TLegend.h" 
#include "TLorentzVector.h"

//This code finds the large jet and (jet) b-tagged Higgs Boson mass using the bbww_x2000_s1500 file

void a_1(){
  gStyle -> SetOptStat("nemr");
  int nbins = 1000;
  int w = 1;				//weight
  
  auto h_Minv = new TH1F("M_{inv}", "Invariant Mass of the Higgs Boson and W from large jet; Mass [MeV]", nbins, 0.0, 2.0e5);
    
  auto file = TFile::Open("bbww_x2000_s1500.root");
  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         jet_pt[14];		//[n_jet]
  Int_t           jet_btagged[14];	//[n_jet]
  Float_t         jet_eta[14];		//[n_jet]
  Float_t         jet_phi[14];		//[n_jet] 
  Float_t         jet_e[14];		//[n_jet]
  Int_t           n_ljet;
  Float_t         ljet_pt[9];   //[n_ljet]
  Float_t         ljet_eta[9];   //[n_ljet]
  Float_t         ljet_phi[9];   //[n_ljet]
  Float_t         ljet_e[9];   //[n_ljet]
  
  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt", jet_pt);
  ctree->SetBranchAddress("jet_btagged",jet_btagged);
  ctree->SetBranchAddress("jet_eta", jet_eta);
  ctree->SetBranchAddress("jet_phi", jet_phi);
  ctree->SetBranchAddress("jet_e", jet_e);
  ctree->SetBranchAddress("n_ljet", &n_ljet);
  ctree->SetBranchAddress("ljet_pt", ljet_pt);
  ctree->SetBranchAddress("ljet_eta", ljet_eta);
  ctree->SetBranchAddress("ljet_phi", ljet_phi);
  ctree->SetBranchAddress("ljet_e", ljet_e);
  
  ctree->SetBranchStatus("*",0);
  
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  ctree->SetBranchStatus("jet_btagged",1);
  ctree->SetBranchStatus("jet_eta",1);
  ctree->SetBranchStatus("jet_phi",1);
  ctree->SetBranchStatus("jet_e", 1);
  ctree->SetBranchStatus("n_ljet",1);
  ctree->SetBranchStatus("ljet_pt", 1);
  ctree->SetBranchStatus("ljet_eta", 1);
  ctree->SetBranchStatus("ljet_phi", 1);
  ctree->SetBranchStatus("ljet_e", 1);
  
  int nentries = ctree->GetEntries();
  cout << nentries << endl;
  
  for (int i=0; i<nentries; i++)			//loops over all the 20,000 entries.
  {
    ctree->GetEntry(i);
    
    for (int j=0; j<n_ljet; j++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
                  
      if (ljet_pt[j] >= 250e3 && abs(ljet_eta[j]) <= 2.0)
      {
	  h_Minv->Fill(ljet.M(), w);
      }
    }  
  } 
  
  file->Close();
  TFile f("a_1.root", "recreate");
  
  h_Minv -> Write();
  h_Minv -> Draw();
  
  f.Close();
} 
 
