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

//This code finds the large jet and (jet) b-tagged Higgs Boson mass

void a_0_W(){
  gStyle -> SetOptStat("nemr");
  int nbins = 1000;
  int w = 1;				//weight
  
  auto h_Higgs = new TH1F("M_{inv} Higgs", "Invariant Mass of the Higgs Boson from large jet with b-tagging; Mass [MeV]", nbins, 0.0, 2.0e5);
  auto h_W = new TH1F("M_{inv} W", "Invariant Mass of W Boson from large jet; Mass [MeV]", nbins, 0.0, 2.0e5);
  
  auto file = TFile::Open("bbww_x1000_s170.root");
  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         jet_pt[14];		//[n_jet]
  Int_t           jet_btagged[14];	//[n_jet]
  Float_t         jet_eta[14];		//[n_jet]
  Float_t         jet_phi[14];		//[n_jet] 
  Float_t         jet_e[14];		//[n_jet]
  Int_t           n_ljet;
  Float_t         ljet_pt[8];   	//[n_ljet]
  Float_t         ljet_eta[8];   	//[n_ljet]
  Float_t         ljet_phi[8];   	//[n_ljet]
  Float_t         ljet_e[8];   		//[n_ljet]
  
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
      
      if (ljet_pt[j] >= 250e3 && abs(ljet_eta[j]) <= 2.0 && ljet.M() > 105e3 && ljet.M() < 145e3)
      {
	for (int k=0; k<n_jet; k++)			//Loops over all regular jets
	{
	  TLorentzVector jet;				//Defines a Lorentz Vector "jet" that will store the regular jet information
	  jet.SetPtEtaPhiE(jet_pt[k], jet_eta[k], jet_phi[k], jet_e[k]);
	  
	  if (jet_btagged[k] >= 1)
	  {
	    float delta_r;				//This variable is the Delta R value between jet and ljet
	    delta_r = ljet.DeltaR(jet);
	    
	    if (delta_r < 0.2)
	    {
	      h_Higgs->Fill(ljet.M(), w);
	    }
	  } 
	}
      }
      
      else
      {
	if (ljet_pt[j] >= 250e3 && abs(ljet_eta[j]) <= 2.0 && ljet.M() > 160e3 /*60e3 && ljet.M() < 100e3*/)
	{
	  for (int l=0; l<n_jet; l++)			//Loops over all regular jets
	  {
	    TLorentzVector jet_0;				//Defines a Lorentz Vector "jet_0" that will store the regular jet information
	    jet_0.SetPtEtaPhiE(jet_pt[l], jet_eta[l], jet_phi[l], jet_e[l]);
	  
	    if (jet_btagged[l] == 0)
	    {
	      float delta_r;				//This variable is the Delta R value between jet and ljet
	      delta_r = ljet.DeltaR(jet_0);
	    
	      if (delta_r < 0.2)
	      {
		h_W->Fill(ljet.M(), w);
	      }
	    } 
	  }
	}
      }
    }    
  } 
  
  file->Close();
  TFile f("a_0_W.root", "recreate");
  
  h_Higgs -> Write();
  h_W -> Write();
  
  c = new TCanvas("canvas", "M_Higgs_W", 2000, 500);
  c -> Divide (2,1);
  
  c -> cd(1);
  h_Higgs -> Draw();
  
  c -> cd(2);
  h_W -> Draw();
  
//   Double_t norm = 1;
//   h_Higgs -> Scale(norm, "width");
//   c -> cd(2);
//   h_Higgs -> Draw();
  
  f.Close();
} 
