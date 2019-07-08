 
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

//This code deals with the boosted analysis. ljet 1-> H; ljet 2->W; ljet 3->W

void a_boosted_semiboosted_fullyboosted(){
  
  gStyle -> SetOptStat("nemr");
  int nbins = 100;
  int nbins_1 = 1000;
  int w = 1;				//weight
  
  auto h_X_boosted = new TH1F("M_{inv} X boosted", "Invariant Mass of the 'X' Particle; Mass[MeV]", nbins_1, 0.0, 2.0e6);
  auto h_X_semiboosted = new TH1F("M_{inv} X semiboosted", "Invariant Mass of the 'X' Particle; Mass[MeV]", nbins_1, 0.0, 2.0e6);
  auto h_X_fullyboosted = new TH1F("M_{inv} X fully-boosted", "Invariant Mass of the 'X' Particle; Mass[MeV]", nbins_1, 0.0, 2.0e6);
  auto file = TFile::Open("bbww_x1000_s170.root");
  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  
  Int_t           n_jet;
  Float_t         jet_pt[14];		//[n_jet]
  Int_t           jet_btagged[14];	//[n_jet]
  Float_t         jet_eta[14];		//[n_jet]
  Float_t         jet_phi[14];		//[n_jet] 
  Float_t         jet_e[14];		//[n_jet]
  Int_t           n_ljet;
  Float_t         ljet_pt[8];   //[n_ljet]
  Float_t         ljet_eta[8];   //[n_ljet]
  Float_t         ljet_phi[8];   //[n_ljet]
  Float_t         ljet_e[8];   //[n_ljet]
  
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
  for (int i=0; i<nentries; i++)
  {
    ctree->GetEntry(i);
    TLorentzVector A;
    TLorentzVector B;
    TLorentzVector C;
    TLorentzVector S;
    TLorentzVector Tot;
    
    float ljet_compare = 0;
    int count_w = 0;
    
    cout << "\nbeginning of " << i << " loop" << endl;
    for (int j=0; j<n_ljet; j++)		//Large jets
    {
      if (ljet_pt[j] >= 250000 && abs(ljet_eta[j]) <= 2.0)		//This loop collects the pt, phi, and eta of the highest jet within the large jets
      {
	float ljet_current = ljet_e[j];
	ljet_compare = max(ljet_compare, ljet_current);
	
	if (ljet_compare == ljet_current && count_w ==0 )
	{
	  cout << "Higgs: " << ljet_compare << endl;
	  A.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	}
	
	if (ljet_current < ljet_compare && count_w == 1)
	{
	  cout << "W1 : " << ljet_current << endl;
	  B.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	}
	
	if (ljet_current < ljet_compare &&  count_w == 2)
	{
	  cout << "W2 : " << ljet_current << endl;
	  C.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	}
	
	count_w++;
      } 
      S = B + C;
      Tot = A + S;
      /*if (S.M() > 0)*/ h_X_boosted->Fill(Tot.M(), w);
      if (S.M() > 160e3) h_X_fullyboosted->Fill(Tot.M(), w);
      
    }
    cout << "end of " << i << " loop\n" << endl;
  }

  for (int k=0; k<nentries; k++)
  {
    ctree->GetEntry(k);

    TLorentzVector HiggsBosonb;			//Higgs Boson from the two b jets (ljet 1 and ljet 2)
    TLorentzVector Jet1;			//1st ljet
    TLorentzVector Jet2;			//2nd ljet
    TLorentzVector X;				//X particle from the Higgs and large jet W
    TLorentzVector W_ljet;			//W particle from the l-jet 
    
    int nbtags = 0;
    int nwtags = 0;
    float ljet_compare = 0;
    int count_w = 0;
    
    for (int l=0; l<n_jet; l++)
    {
      if (jet_btagged[l]>0 && nbtags < 2)
      {
	
	if (nbtags == 0)
	{ 
	  Jet1.SetPtEtaPhiE(jet_pt[l], jet_eta[l], jet_phi[l], jet_e[l]);
	}
      
	if (nbtags == 1)
	{
	  Jet2.SetPtEtaPhiE(jet_pt[l], jet_eta[l], jet_phi[l], jet_e[l]);
	  HiggsBosonb = Jet1 + Jet2;
	}
	nbtags++;

      }
      if (ljet_pt[l] >= 250000 && abs(ljet_eta[l]) <= 2.0)		//This loop collects the pt, phi, and eta of the highest jet within the large jets
      {
	float ljet_current = ljet_e[l];
	ljet_compare = max(ljet_compare, ljet_current);
	if (ljet_current < ljet_compare && count_w == 1)
	{
	  W_ljet.SetPtEtaPhiE(ljet_pt[l], ljet_eta[l], ljet_phi[l], ljet_e[l]);
	}
	count_w++;
      }
      
      
      X = HiggsBosonb + W_ljet;
      /*if (W_ljet.M() >0 ) */h_X_semiboosted->Fill(X.M(), w);
    }
  }
  
  file->Close();
  TFile f("a_boosted_semiboosted_fullyboosted.root", "recreate");
  
  

  h_X_boosted->Write();
  h_X_semiboosted -> Write();
  h_X_fullyboosted -> Write();
  
  c = new TCanvas("c", "M_X", 1000, 500);
  c -> Divide (1,1);
  Double_t norm = 1;

  h_X_boosted->Scale(norm, "width");
  h_X_semiboosted -> Scale(norm, "width");
  h_X_fullyboosted -> Scale(norm, "width");
  
  c -> cd(1);
  h_X_boosted 		-> Draw();
  h_X_semiboosted 	-> Draw("same");
  h_X_fullyboosted 	-> Draw("same");
  
  h_X_boosted 		-> SetLineColor(kBlue+1);
  h_X_semiboosted 	-> SetLineColor(2);
  h_X_fullyboosted 	-> SetLineColor(3);
  
  auto legend = new TLegend(0.6,0.6,0.9,0.75);	// (x1, y1, x2, y2)
  legend->SetHeader("M_{inv} of 'X' from large jets");
  legend->AddEntry(h_X_boosted,"Boosted: 1 H + 2 W");
  legend->AddEntry(h_X_semiboosted,"Semi-Boosted: 1 H + 1 W");
  legend->AddEntry(h_X_fullyboosted, "Fully-Boosted: 1 H + WW jet (m_{ljet} > 160GeV)");
  legend->Draw();

  f.Close();
    
}
