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


//This code graphs the invariant mass of regular and large b-jets // i am using a ttbar monte carlo simulated data
void bkg_ttbar() {
  
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemr");
  int w = 1; 					//weight of the bins
  int nbin = 1000;				//bin size
  
  auto h_M_inv = new TH1F("M_inv", "Invariant Mass of Regular and Large Jets; Mass [GeV]; Events", nbin, 0.0, 2.0e6);
						//Creates a pointer to the histogram to be created.
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
  
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("jet_isbtagged_MV2c10_85", &jet_btagged);
  t->SetBranchAddress("jet_eta", &jet_eta);
  t->SetBranchAddress("jet_phi", &jet_phi);
  t->SetBranchAddress("jet_e", &jet_e);
  t->SetBranchAddress("ljet_pt", &ljet_pt);
  t->SetBranchAddress("ljet_eta", &ljet_eta);
  t->SetBranchAddress("ljet_phi", &ljet_phi);
  t->SetBranchAddress("ljet_m", &ljet_m);
  
  int nentries = t->GetEntries();
  
  
  cout << "Number of entries: " << nentries << endl;
  
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int n_jet = jet_pt->size();
    int n_ljet = ljet_pt->size();

    for (int j=0; j<n_jet; j++)
    {
      if ((int)(*jet_btagged)[j] == 1)
      {
	cout << "pt: " << (*jet_pt)[j] << "\teta: " << (*jet_eta)[j] << "\tphi: " << (*jet_phi)[j] 
	<< "\tE :" << (*jet_e)[j] << endl;
	TLorentzVector Minv1;
	Minv1.SetPtEtaPhiE((*jet_pt)[j], (*jet_eta)[j], (*jet_phi)[j], (*jet_e)[j]);
	cout << Minv1(0) << endl;
      }
//       cout << "passed" << j << endl;
      
    }
    
  }
  TCanvas *c = new TCanvas("canvas", "b", 750, 500);
  h_M_inv->Draw();
    
  
  
  file->Close();
  TFile f("bkg_ttbar.root", "RECREATE");
  h_M_inv->Write(); 
  f.Close();
  
}