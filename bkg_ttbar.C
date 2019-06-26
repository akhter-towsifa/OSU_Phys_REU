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
#include "TLorentzVector.h"
#include "TMath.h"
#include <math.h>
#include <cmath>
#include "TCanvas.h"
#include "TLatex.h"


//This code graphs the invariant mass of regular and large b-jets // i am using a ttbar monte carlo simulated data
void bkg_ttbar() {
  
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemr");
  int w = 1; 					//weight of the bins
  int nbin = 1000;				//bin size
  
  auto h_reg_M_inv = new TH1F("M_inv1", "Invariant Mass of Regular Jet; Mass [MeV]; Events", nbin, 0.0, 6.0e5);
  auto h_large_M_inv = new TH1F("M_inv2", "Invariant Mass of Large Jet; Mass [MeV]; Events", nbin, 0.0, 6.0e5);
  auto h_M_inv = new TH1F("M_inv", "Invariant Mass of t#bar{t}; Mass [MeV]; Events", nbin, 0.0, 2.0e6);
  auto file = TFile::Open("ttbar.root");	//Creates a pointer to the main data file
  TTree* t = (TTree*)file->Get("nominal");	//Creates a pointer to the tree
  
  vector<float> *jet_pt		=0;
  vector<char>	*jet_btagged	=0;
  vector<float>	*jet_eta	=0;
  vector<float>	*jet_phi	=0;
  vector<float> *jet_e		=0;
  vector<float> *ljet_pt	=0;
  vector<float> *ljet_eta	=0;
  vector<float> *ljet_phi	=0;
  vector<float> *ljet_m		=0;
  
  t->SetBranchAddress("jet_pt", 			&jet_pt);
  t->SetBranchAddress("jet_isbtagged_MV2c10_85", 	&jet_btagged);
  t->SetBranchAddress("jet_eta", 			&jet_eta);
  t->SetBranchAddress("jet_phi", 			&jet_phi);
  t->SetBranchAddress("jet_e", 				&jet_e);
  t->SetBranchAddress("ljet_pt", 			&ljet_pt);
  t->SetBranchAddress("ljet_eta", 			&ljet_eta);
  t->SetBranchAddress("ljet_phi",			&ljet_phi);
  t->SetBranchAddress("ljet_m", 			&ljet_m);
  
  int nentries = t->GetEntries();
  
  
  cout << "Number of entries: " << nentries << endl;
  
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int n_bjets = 0;
    int n_jet = jet_pt->size();
    int n_ljet = ljet_pt->size();
    TLorentzVector Minv1;
    TLorentzVector Minv2;
    TLorentzVector Minv_total;
    
    for (int j=0; j<n_jet; j++)
    {
      if ((int)(*jet_btagged)[j] == 1 && n_bjets < 1)
      {
	Minv1.SetPtEtaPhiE((*jet_pt)[j], (*jet_eta)[j], (*jet_phi)[j], (*jet_e)[j]);
	h_reg_M_inv->Fill(Minv1.M(), w);
	n_bjets++;
      } 
    }
    
    float ljet_m_compare = 0;
    int n_t_candidates=0;			//a step variable that counts the number of large jets or t quarks
    
    for (int k=0; k<n_ljet; k++)		//Large jets
    {
      if ((float)(*ljet_pt)[k] >= 420000 && abs((float)(*ljet_eta)[k]) <= 2.0 && n_bjets==1)		//This loop collects the pt, phi, and eta of the highest jet within the large jets
      {
	n_t_candidates++;
	float ljet_m_current = (float)(*ljet_m)[k];
	ljet_m_compare = max(ljet_m_compare, ljet_m_current);
	
	if (ljet_m_compare == ljet_m_current)
	{
	  Minv2.SetPtEtaPhiM((*ljet_pt)[k], (*ljet_eta)[k], (*ljet_phi)[k], (*ljet_m)[k]);
	  h_large_M_inv->Fill(Minv2.M(), w);
	}
      }
    }
    
    if (n_bjets==1 && n_t_candidates>0)
    {
      Minv_total = Minv1 + Minv2;
//       cout << i << "\tMass: " << Minv_total.M() << endl;
      h_M_inv->Fill(Minv_total.M(), w);
    }
    
    
  }
  TCanvas *c = new TCanvas("canvas", "b", 2000, 500);
  c -> Divide (3,1);
  c->cd(1);
  h_reg_M_inv->Draw();
  c->cd(2);
  h_large_M_inv->Draw();
  c->cd(3);
  h_M_inv->Draw();
    
  
  
  file->Close();
  TFile f("bkg_ttbar.root", "RECREATE");
  h_reg_M_inv->Write();
  h_large_M_inv->Write();
  h_M_inv->Write();
  f.Close();
  
}