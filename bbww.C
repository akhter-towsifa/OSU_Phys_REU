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

//This code is regarding all the bbww files

void bbww(){
  
  gStyle->SetOptStat("nemr");
  int nbins = 1000;
  int w = 1;				//weight
  
  auto h_jet_pt 		= new TH1F("jet_pt","Jet_P_{T}; Mass [MeV]; Counts",		nbins,0.0, 1.0e6);
  auto h_jet_pt_btag 		= new TH1I("jet_btagged","Jet p_{T} b-tagged;Mass[MeV]",	nbins,0.0, 1.0e6);
  auto h_leading_bjet_1 	= new TH1F("h_leading_bjet_1", "First Leading b-Jet; Mass[MeV]; Counts", nbins, 0.0, 2.0e6);
  auto h_leading_bjet_2 	= new TH1F("h_leading_bjet_2", "Leading and Subleading b-Jets; Mass[MeV]; Counts", nbins, 0.0, 2.0e6);
  auto h_M_inv 			= new TH1F("M_inv", "Invariant Mass of Higgs Boson from the b-Jets; Mass[MeV]; Counts",	nbins,0.0, 4.0e5);
  auto h_w_pt_btag		= new TH1F("W_particle", "W P_{T}; Mass[MeV]; Counts", 		nbins, 0.0, 2.0e6);
  auto h_w_inv			= new TH1F("W1_inv", "Invariant mass of W; mass[MeV]; Counts", 	nbins, 0.0, 2.0e6);
  auto h_w2_inv			= new TH1F("W2_inv", "Invariant mass of 2nd W; mass[MeV]; Counts", 	nbins, 0.0, 2.0e6);

  auto file = TFile::Open("bbww_x1000_s170.root");

  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         jet_pt[14];		//[n_jet]
  Int_t           jet_btagged[14];	//[n_jet]
  Float_t         jet_eta[14];		//[n_jet]
  Float_t         jet_phi[14];		//[n_jet]
  Float_t         jet_e[14];		//[n_jet]

  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt", jet_pt);
  ctree->SetBranchAddress("jet_btagged",jet_btagged);
  ctree->SetBranchAddress("jet_eta", jet_eta);
  ctree->SetBranchAddress("jet_phi", jet_phi);
  ctree->SetBranchAddress("jet_e", jet_e);
  
  ctree->SetBranchStatus("*",0);
  
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  ctree->SetBranchStatus("jet_btagged",1);
  ctree->SetBranchStatus("jet_eta",1);
  ctree->SetBranchStatus("jet_phi",1);
  ctree->SetBranchStatus("jet_e", 1);
  
  int nentries = ctree->GetEntries();
//   cout << nentries << endl;
  
  for (int i=0; i<nentries; i++)
  {
    ctree->GetEntry(i);

    TLorentzVector HiggsBosonb;
    TLorentzVector Jet1;
    TLorentzVector Jet2;
    TLorentzVector W_1;
    TLorentzVector W_1_a;
    TLorentzVector W_1_b;
    TLorentzVector W_2;
    TLorentzVector W_2_a;
    TLorentzVector W_2_b;

    cout << i << endl;
    
    int nbtags = 0;
    int nwtags = 0;
    
    for (int j=0; j<n_jet; j++)
    { cout << "\t" << j << endl;
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
	  cout << "\t\t1: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_1_a.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  
	}
      
	if (nwtags == 1)
	{
	  cout << "\t\t2: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_1_b.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  W_1 = W_1_a + W_1_b;
	  h_w_inv->Fill(W_1.M(), w);
	}
	
	if (nwtags == 2)
	{ 
	  cout << "\t\t3: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_2_a.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  
	}
      
	if (nwtags == 3)
	{
	  cout << "\t\t4: " << jet_pt[j] << "\t" << jet_eta[j] << "\t" << jet_phi[j] << "\t" << jet_e[j] << endl;
	  W_2_b.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
	  W_2 = W_2_a + W_2_b;
	  h_w2_inv->Fill(W_2.M(), w);
	}
	nwtags++;
      }
    }
    
//     if (nbtags >= 2)
//     {
//       HiggsBosonb = Jet1 + Jet2;
//       h_M_inv->Fill(HiggsBosonb.M(), w);
//     }
    
   
    for(int k = 0; k< n_jet; k++)
    {
      h_jet_pt->Fill( jet_pt[k], w);				//Creates all the jet_pt histogram
      if (jet_btagged[k]>0) h_jet_pt_btag->Fill(jet_pt[k], w);
    }
  }
  
  file->Close();
  TFile f("all_bbww.root", "recreate");
  
  h_jet_pt     ->Write();
  h_jet_pt_btag->Write();
  h_leading_bjet_1->Write();
  h_leading_bjet_2->Write();
  h_M_inv	->Write();
  h_w_pt_btag	->Write();
  h_w_inv	->Write();
  h_w2_inv	->Write();
  
  c = new TCanvas("c", "bbWW Particle Jets from PP Collision", 2000, 1000);
  c -> Divide (3,3);
  
  c -> cd(1);
  h_jet_pt -> Draw();
  
  c -> cd(4);
  h_jet_pt_btag -> Draw();
  
  c -> cd(7);
  h_leading_bjet_2 -> Draw();
  h_leading_bjet_1 -> Draw("same");
  h_leading_bjet_2 -> SetLineColor(2);
  auto legend = new TLegend(0.7,0.6,0.9,0.7);
  legend->SetHeader("b-Jets");
  legend->AddEntry(h_leading_bjet_1,"Leading b-Jet");
  legend->AddEntry("h_leading_bjet_2","Subleading b-Jet");
  legend->Draw();
  
  c -> cd(5);
  h_M_inv -> Draw();
  
  c -> cd(3);
  h_w_pt_btag -> Draw();
  
  c -> cd(6);
  h_w_inv -> Draw();
  
  c-> cd(9);
  h_w2_inv -> Draw();
  
  f.Close();
}
