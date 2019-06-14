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

//This code creates a histogram for the missing transverse energy.

void missing_transverse_energy() {
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 1000;
  int w = 1;
 
  auto h_MET = new TH1F("met", "Missing Transverse Energy from the Jets; x-axis; Counts", nbins,0.0, 4.5e5);
  auto file = TFile::Open("bbww_x1000_s170.root");
  
  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         met;
  Float_t         met_x;
  Float_t         met_y;
  Float_t         met_phi;
  Float_t         met_sumet;

  ctree->SetBranchAddress("n_jet", &n_jet);
//   ctree->SetBranchAddress("jet_pt", &jet_pt);
  ctree->SetBranchAddress("met", &met);
  ctree->SetBranchAddress("met_x", &met_x);
  ctree->SetBranchAddress("met_y", &met_y);
  ctree->SetBranchAddress("met_phi", &met_phi);
  ctree->SetBranchAddress("met_sumet", &met_sumet);
  
  ctree->SetBranchStatus("*",0);
  
  ctree->SetBranchStatus("n_jet", 1);
//   ctree->SetBranchStatus("jet_pt", 1);
  ctree->SetBranchStatus("met", 1);
  ctree->SetBranchStatus("met_x", 1);
  ctree->SetBranchStatus("met_y", 1);
  ctree->SetBranchStatus("met_phi", 1);
  ctree->SetBranchStatus("met_sumet", 1);
  
  int nentries = ctree->GetEntries();
  TFile f("missing_transverse_energy.root","recreate");
  
  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    
    for (int j=0; j<n_jet; j++){
      h_MET->Fill( met, w);
    }
    
  }
  
  h_MET->Write();
  
  f.Write();
  f.Close();  
}