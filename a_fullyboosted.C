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

//This code deals with the fully boosted analysis. H-> bb; ljet 1->W; ljet 2->W

void a_fullyboosted(){
  
  gStyle -> SetOptStat("nemr");
  int nbins = 100;
  int nbins_1 = 1000;
  int w = 1;				//weight
  
  auto h_M_inv_Higgs = new TH1F("M_{inv} Higgs", "Invariant Mass of Higgs Boson; Mass[MeV]; Counts", nbins,0.0, 4.0e5);
  auto h_M_inv_W1 = new TH1F("M_{inv} W1", "Invariant Mass of the first W Boson; Mass[MeV]; Counts", nbins, 0.0, 4.0e5);
  auto h_M_inv_W2 = new TH1F("M_{inv} W2", "Invariant Mass of the second W Boson; Mass[MeV]; Counts", nbins, 0.0, 4.0e5);
  auto h_M_inv_S = new TH1F("M_{inv} S", "Invariant Mass of the S particle; Mass[MeV]; Counts", nbins, 0.0, 4.0e5);
  auto h_M_inv_X = new TH1F("M_{inv} X", "Invariant Mass of the 'X' Particle; Mass[MeV]; Counts", nbins_1, 0.0, 2.0e6);
  
  auto file = TFile::Open("bbww_x1000_s170.root");
  TTree* t = (TTree*)file->Get("CollectionTree");
  
  
  Int_t           n_ljet;
  Float_t         ljet_pt[8];   //[n_ljet]
  Float_t         ljet_eta[8];   //[n_ljet]
  Float_t         ljet_phi[8];   //[n_ljet]
  Float_t         ljet_e[8];   //[n_ljet]
  
  t->SetBranchAddress("n_ljet", &n_ljet);
  t->SetBranchAddress("ljet_pt", ljet_pt);
  t->SetBranchAddress("ljet_eta", ljet_eta);
  t->SetBranchAddress("ljet_phi", ljet_phi);
  t->SetBranchAddress("ljet_e", ljet_e);
  
  t->SetBranchStatus("*",0);
  
  t->SetBranchStatus("n_ljet",1);
  t->SetBranchStatus("ljet_pt", 1);
  t->SetBranchStatus("ljet_eta", 1);
  t->SetBranchStatus("ljet_phi", 1);
  t->SetBranchStatus("ljet_e", 1);
  
  int nentries = t->GetEntries();
  cout << nentries << endl;
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
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
// 	cout << "\t\t\t\tcomp: " << ljet_compare << "\tcurr: " << ljet_current << endl;
	ljet_compare = max(ljet_compare, ljet_current);
	
	if (ljet_compare == ljet_current && count_w ==0 )
	{
	  cout << "Higgs: " << ljet_compare << endl;
	  A.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	  h_M_inv_Higgs->Fill(A.M(), w);
	}
	
	if (ljet_current < ljet_compare && count_w == 1)
	{
	  cout << "W1 : " << ljet_current << endl;
	  B.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	  h_M_inv_W1->Fill(B.M(), w);
	}
	
	if (ljet_current < ljet_compare &&  count_w == 2)
	{
	  cout << "W2 : " << ljet_current << endl;
	  C.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	  h_M_inv_W2->Fill(C.M(), w);
	}
	
	count_w++;
      } 
      S = B + C;
      if (S.M()> 0.0) h_M_inv_S->Fill(S.M(), w);
      Tot = A + B + C;
      h_M_inv_X->Fill(Tot.M(), w);
      
    }
    cout << "end of " << i << " loop\n" << endl;
  }

  
  file->Close();
  TFile f("a_fullyboosted.root", "recreate");
  

  h_M_inv_Higgs->Write();
  h_M_inv_W1->Write();
  h_M_inv_W2->Write();
  h_M_inv_S->Write();
  h_M_inv_X->Write();
  
  c = new TCanvas("c", "Fully Boosted Analysis", 2000, 1000);
  c -> Divide (3,2);
  
  c -> cd(1);
  h_M_inv_Higgs -> Draw();
  
  c -> cd(2);
  h_M_inv_W1 -> Draw();
  
  c -> cd(3);
  h_M_inv_W2 -> Draw();
  
  c -> cd(4);
  h_M_inv_S -> Draw();
  
  c -> cd(5);
  h_M_inv_X -> Draw();
  
  f.Close();
}
 
