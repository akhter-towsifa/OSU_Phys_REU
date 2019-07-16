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

//This code finds the large jet and Higgs Boson, S masses using the bbww_x2000_s1500 file (large jet: 1H, >1 W's)

void a_1_x2000_s1500(){
  gStyle -> SetOptStat("nemr");
  int nbins = 1000;
  int w = 1;				//weight
  
  auto h_Minv = new TH1F("M_{inv}", "Invariant Mass of the Higgs Boson and W from large jet; Mass [MeV]", nbins, 0.0, 2.0e5);
  auto h_Higgs = new TH1F("M_{inv} Higgs", "Invariant Mass of the Higgs Boson from large jet; Mass [MeV]", nbins, 0.0, 2.0e5);
  auto h_W1 = new TH1F("M_{inv} W1", "Invariant Mass of first W Boson from large jet; Mass [MeV]", nbins, 0.0, 2.0e5);
  auto h_W2 = new TH1F("M_{inv} W2", "Invariant Mass of second W Boson from large jet; Mass [MeV]", nbins, 0.0, 2.0e5);
  auto h_S = new TH1F("M_{inv} S", "Invariant Mass of the S particle from W's; Mass [MeV]", nbins, 0.0, 3.0e6);
  auto h_S_totalnumber = new TH1F("M_{inv} S number", "Invariant Mass of the S particle from W's; Mass [MeV]", nbins, 0.0, 3.0e6);
  auto h_X = new TH1F("M_{inv} X", "Invariant Mass of the X Particle from large jet Higgs and W; Mass [MeV]", nbins, 0.0, 3.0e6);
  auto h_X_totalnumber = new TH1F("M_{inv} X number", "Invariant Mass of the X Particle from large jet Higgs and W; Mass [MeV]", nbins, 0.0, 3.0e6);
  
  auto file = TFile::Open("bbww_x2000_s1500.root");
  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_ljet;
  Float_t         ljet_pt[9];   //[n_ljet]
  Float_t         ljet_eta[9];   //[n_ljet]
  Float_t         ljet_phi[9];   //[n_ljet]
  Float_t         ljet_e[9];   //[n_ljet]
  
  ctree->SetBranchAddress("n_ljet", &n_ljet);
  ctree->SetBranchAddress("ljet_pt", ljet_pt);
  ctree->SetBranchAddress("ljet_eta", ljet_eta);
  ctree->SetBranchAddress("ljet_phi", ljet_phi);
  ctree->SetBranchAddress("ljet_e", ljet_e);
  
  ctree->SetBranchStatus("*",0);
  
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
    TLorentzVector H;					//Higgs Lorentz Vector
    TLorentzVector W1;					//First W from ljet2
    TLorentzVector W2;					//Second W from ljet3
    TLorentzVector S;					//S Lorentz Vector
    TLorentzVector X;					//X Lorentz Vector

    
    int count_w = 0;					//counting the number of W's
    
    for (int j=0; j<n_ljet; j++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
                  
      if (ljet_pt[j] >= 250e3 && abs(ljet_eta[j]) <= 2.0)
      {	
	h_Minv->Fill(ljet.M(), w);
	
	if (ljet.M() > 105e3 && ljet.M() < 145e3)
	{
	  H.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	  h_Higgs->Fill(H.M(), w);
	}
	
	if (ljet.M() > 60e3 && ljet.M() < 100e3)
	{
	  count_w++;
	  if (count_w == 1)
	  {
	    W1.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	    h_W1->Fill(W1.M(), w);
	  }
	  else
	  {
	    W2.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	    h_W2->Fill(W2.M(), w);
	  }
	}
	S = W1 + W2;
	/*if(S.M()>0) */h_S->Fill(S.M(), w);
	if(S.M() > 1400e3 && S.M() < 1550e3) h_S_totalnumber->Fill(S.M(), w);
	
      }
    }  
    X = H + S;
    /*if (X.M() > 0) */h_X->Fill(X.M(), w);
    if (X.M() > 1850e3 && X.M() < 2100e3) h_X_totalnumber->Fill(X.M(), w);
    
  } 
  
  file->Close();
  TFile f("a_1.root", "recreate");
  
  h_Minv -> Write();
  h_Higgs -> Write();
  h_W1 -> Write();
  h_W2 -> Write();
  h_S -> Write();
  h_S_totalnumber ->  Write();
  h_X -> Write();
  h_X_totalnumber -> Write();
  
  c = new TCanvas("canvas", "M_Higgs_S_X", 2000, 1000);
  c -> Divide (3,2);
  
  c -> cd(1);
  h_Minv -> Draw();
  
  c -> cd(2);
  h_Higgs -> Draw();
  
  c -> cd(3);
  h_W1 -> Draw();
  
  c -> cd(4);
  h_W2 -> Draw();
  
  c -> cd(5);
  h_S -> Draw();
  
  c -> cd(6);
  h_X -> Draw();
  
  c1 = new TCanvas("canvas1", "M_S_X", 2000, 500);
  c1 -> Divide (2,1);
  
  c1 -> cd(1);
  h_S_totalnumber -> Draw();
  
  c1 -> cd(2);
  h_X_totalnumber -> Draw();
  
  f.Close();
} 
 
