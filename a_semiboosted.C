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

//This code deals with the semi-boosted analysis. H->bb (regular); ljet 1->W; regular jets->W

void a_semiboosted(){
  
  gStyle->SetOptStat("nemr");
  int nbins = 1000;
  int nbins_1 = 100;
  int w = 1;			//weight
  
  auto h_jet_pt 		= new TH1F("jet_pt","Jet_P_{T}; Mass [MeV]; Counts",		nbins,0.0, 1.0e6);
  auto h_jet_pt_btag 		= new TH1I("jet_btagged","Jet p_{T} b-tagged;Mass[MeV]",	nbins,0.0, 1.0e6);
  auto h_leading_bjet_1 	= new TH1F("h_leading_bjet_1", "First Leading b-Jet; Mass[MeV]; Counts", nbins, 0.0, 2.0e6);
  auto h_leading_bjet_2 	= new TH1F("h_leading_bjet_2", "Leading and Subleading b-Jets; Mass[MeV]; Counts", nbins, 0.0, 2.0e6);
  auto h_M_inv 			= new TH1F("M_inv", "Invariant Mass of Higgs Boson from the b-Jets; Mass[MeV]; Counts",	nbins,0.0, 4.0e5);
  auto h_w_pt_btag		= new TH1F("W_particle", "W P_{T}; Mass[MeV]; Counts", 		nbins, 0.0, 2.0e6);
  auto h_w_inv			= new TH1F("W1_inv", "Invariant mass of W; mass[MeV]; Counts", 	nbins, 0.0, 2.0e6);
  auto h_w2_inv			= new TH1F("W2_inv", "Invariant mass of 2nd W; mass[MeV]; Counts", 	nbins, 0.0, 2.0e6);
  auto h_M_inv_2 		= new TH1F("M_inv_2", "Invariant Mass of the S particle from the W bosons; Mass[MeV]; Counts",nbins,0.0, 2.0e6);
  auto h_X	 		= new TH1F("X_M_inv", "Invariant Mass of the X Particle; Mass[MeV]; Counts",nbins,0.0, 4.0e6);
  auto h_M_inv_W_ljet = new TH1F("M_{inv} W_{Ljet}", "Invariant Mass of the W Boson from Large Jet; Mass[MeV]; Counts", nbins_1, 0.0, 4.0e5);

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
//   cout << nentries << endl;
  
  for (int i=0; i<nentries; i++)
  {
    ctree->GetEntry(i);

    TLorentzVector HiggsBosonb;			//Higgs Boson from the two b jets
    TLorentzVector Jet1;			//1st b jet
    TLorentzVector Jet2;			//2nd b jet
    TLorentzVector W_1;				//1st W Boson
    TLorentzVector W_1_a;			//1st quark jet from the 1st W Boson
    TLorentzVector W_1_b;			//2nd quark jet from the 1st W Boson
    TLorentzVector W_2;				//2nd W Boson
    TLorentzVector W_2_a;			//1st quark jet from the 2nd W Boson 
    TLorentzVector W_2_b;			//2nd quark jet from the 2nd W Boson
    TLorentzVector HiggsBosonw;			//2nd Higgs Boson from the W bosons
    TLorentzVector X;				//X particle from the diHiggs
    TLorentzVector W_ljet;			//W particle from the l-jet 
    
    int nbtags = 0;
    int nwtags = 0;
    float ljet_compare = 0;
    int count_w = 0;
    
    for (int j=0; j<n_jet; j++)
    {
      if (jet_btagged[j]>0 && nbtags < 2)
      {
	
	if (nbtags == 0)
	{ 
	  h_leading_bjet_1->Fill(jet_pt[j]);			//Creates histogram for the leading b-jet
	  Jet1.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	}
      
	if (nbtags == 1)
	{
	  h_leading_bjet_2->Fill(jet_pt[j]);			//Creates histogram for the subleading b-jet
	  Jet2.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  HiggsBosonb = Jet1 + Jet2;
	  h_M_inv->Fill(HiggsBosonb.M(), w);
	}
	nbtags++;

      }
      
      else
      {
	h_w_pt_btag->Fill(jet_pt[j], w);
	if (nwtags == 0)
	{ 
// 	  cout << "\t\t1: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_1_a.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  
	}
      
	if (nwtags == 1)
	{
// 	  cout << "\t\t2: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_1_b.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  W_1 = W_1_a + W_1_b;
	  h_w_inv->Fill(W_1.M(), w);
	}
	
	if (nwtags == 2)
	{ 
// 	  cout << "\t\t3: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_2_a.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);	  
	}
      
	if (nwtags == 3)
	{
// 	  cout << "\t\t4: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_2_b.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  W_2 = W_2_a + W_2_b;
	  h_w2_inv->Fill(W_2.M(), w);
// 	  HiggsBosonw = W_1 + W_2;
// 	  h_M_inv_2->Fill(HiggsBosonw.M(), w);
	}
	
	nwtags++;
      }
      HiggsBosonw = W_1 + W_2;
      h_M_inv_2->Fill(HiggsBosonw.M(), w);
      
      if (ljet_pt[j] >= 250000 && abs(ljet_eta[j]) <= 2.0)		//This loop collects the pt, phi, and eta of the highest jet within the large jets
      {
	float ljet_current = ljet_e[j];
	ljet_compare = max(ljet_compare, ljet_current);
	if (ljet_current < ljet_compare && count_w == 1)
	{
// 	  cout << "W1 : " << ljet_current << endl;
	  W_ljet.SetPtEtaPhiE(ljet_pt[j], ljet_eta[j], ljet_phi[j], ljet_e[j]);
	  h_M_inv_W_ljet->Fill(W_ljet.M(), w);
	}
	count_w++;
      }
      
      X = HiggsBosonb + HiggsBosonw + W_ljet;
      h_X->Fill(X.M(), w);
    }
    
   
    for(int k = 0; k< n_jet; k++)
    {
      h_jet_pt->Fill( jet_pt[k], w);				//Creates all the jet_pt histogram
      if (jet_btagged[k]>0) h_jet_pt_btag->Fill(jet_pt[k], w);	//Creates the jet_pt btagged histogram
    }
  }
  
  file->Close();
  TFile f("a_semiboosted.root", "recreate");
  

  h_M_inv	->Write();
  h_w_pt_btag	->Write();
  h_w_inv	->Write();
  h_w2_inv	->Write();
  h_M_inv_2	->Write();
  h_X		->Write();
  h_M_inv_W_ljet->Write();
  
  c = new TCanvas("c", "Semi-Boosted Analysis", 2000, 1000);
  c -> Divide (3,2);
  
  
  c -> cd(4);
  h_M_inv -> Draw();
  
//   c -> cd(5);
//   h_M_inv_2 -> Draw();
  
  c -> cd(6);
  h_X -> Draw();
  
  c -> cd(1);
  h_w_pt_btag -> Draw();
  
  c -> cd(2);
  h_w_inv -> Draw();
  
  c-> cd(3);
  h_w2_inv -> Draw();
  
  c -> cd(5);
  h_M_inv_W_ljet -> Draw();
  /////////////
  

  
  f.Close();
}
 
