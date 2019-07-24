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
#include "TGraph.h"
#include "TGraphPainter.h"
#include "TGraphErrors.h"

 
//This code is to find the selection efficiency for all the different mass point and all the boosted types (taking regular jet b-tagging into account).

void b_selectionEfficiency_3btags(){
  
  float smass[6];	//smass list collects all the 6 S mass; fraction1,2,3 lists collect the nsel/nentries ratio at each of these mass points 
  float fraction1a[6], fraction2a[6], fraction3a[6];
  float fraction1b[6], fraction2b[6], fraction3b[6];
  float fraction1c[6], fraction2c[6], fraction3c[6];
  float fraction1aerror[6], fraction2aerror[6], fraction3aerror[6];
  float fraction1berror[6], fraction2berror[6], fraction3berror[6];
  float fraction1cerror[6], fraction2cerror[6], fraction3cerror[6];
  
////////////////////////////////X 1000 S 170 //////////////////////////////////
  
  auto file0 = TFile::Open("bbww_x1000_s170.root");
  TTree* t0 = (TTree*)file0->Get("CollectionTree");
  
  Int_t           n_jet0;
  Float_t         jet_pt0[14];		//[n_jet]
  Int_t           jet_btagged0[14];	//[n_jet]
  Float_t         jet_eta0[14];		//[n_jet]
  Float_t         jet_phi0[14];		//[n_jet] 
  Float_t         jet_e0[14];		//[n_jet]
  Int_t           n_ljet0;
  Float_t         ljet_pt0[8];   	//[n_ljet]
  Float_t         ljet_eta0[8];   	//[n_ljet]
  Float_t         ljet_phi0[8];   	//[n_ljet]
  Float_t         ljet_e0[8];   	//[n_ljet]
  Int_t           n_tjet0;
  Float_t         tjet_pt0[28];   //[n_tjet]
  Float_t         tjet_eta0[28];   //[n_tjet]
  Float_t         tjet_phi0[28];   //[n_tjet]
  Float_t         tjet_e0[28];   //[n_tjet]
  Int_t           tjet_btagged0[28];   //[n_tjet]  
  
  t0->SetBranchAddress("n_jet", &n_jet0);
  t0->SetBranchAddress("jet_pt", jet_pt0);
  t0->SetBranchAddress("jet_btagged",jet_btagged0);
  t0->SetBranchAddress("jet_eta", jet_eta0);
  t0->SetBranchAddress("jet_phi", jet_phi0);
  t0->SetBranchAddress("jet_e", jet_e0);
  t0->SetBranchAddress("n_ljet", &n_ljet0);
  t0->SetBranchAddress("ljet_pt", ljet_pt0);
  t0->SetBranchAddress("ljet_eta", ljet_eta0);
  t0->SetBranchAddress("ljet_phi", ljet_phi0);
  t0->SetBranchAddress("ljet_e", ljet_e0);
  t0->SetBranchAddress("n_tjet", &n_tjet0);
  t0->SetBranchAddress("tjet_pt", tjet_pt0);
  t0->SetBranchAddress("tjet_eta", tjet_eta0);
  t0->SetBranchAddress("tjet_phi", tjet_phi0);
  t0->SetBranchAddress("tjet_e", tjet_e0);
  t0->SetBranchAddress("tjet_btagged", tjet_btagged0);  
  
  t0->SetBranchStatus("*",0);
  
  t0->SetBranchStatus("n_jet",1);
  t0->SetBranchStatus("jet_pt",1);
  t0->SetBranchStatus("jet_btagged",1);
  t0->SetBranchStatus("jet_eta",1);
  t0->SetBranchStatus("jet_phi",1);
  t0->SetBranchStatus("jet_e", 1);
  t0->SetBranchStatus("n_ljet",1);
  t0->SetBranchStatus("ljet_pt", 1);
  t0->SetBranchStatus("ljet_eta", 1);
  t0->SetBranchStatus("ljet_phi", 1);
  t0->SetBranchStatus("ljet_e", 1);
  t0->SetBranchStatus("n_tjet", 1);
  t0->SetBranchStatus("tjet_pt", 1);
  t0->SetBranchStatus("tjet_eta", 1);
  t0->SetBranchStatus("tjet_phi", 1);
  t0->SetBranchStatus("tjet_e", 1);
  t0->SetBranchStatus("tjet_btagged", 1);
  
  int nentries0 = t0->GetEntries();
  int nsel01a = 0;		//counting fully boosted events that are selected
  int nsel01b = 0;
  int nsel01c = 0;
  int nsel02a = 0;		//counting boosted events that are selected
  int nsel02b = 0;
  int nsel02c = 0;
  int nsel03a = 0;		//counting semi-boosted events that are selected
  int nsel03b = 0;
  int nsel03c = 0;
  
  cout << "\nnentries0: " << nentries0 << "\tX: 1000 GeV \tS: 170 Gev" << endl;
  for (int i=0; i<nentries0; i++)			//loops over all the 20,000 entries in file0
  {
    t0->GetEntry(i);
    int count01a = 0;		//counts the number of higgs within each Ljet
    int count01b = 0;
    int count01c = 0;
    int count02 = 0;		//counts the number of WW Ljets
    int count03 = 0;		//counts the number of W Ljets
    int count04 = 0;		//counts the number of regular W jets
    
    for (int i0=0; i0<n_ljet0; i0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt0[i0], ljet_eta0[i0], ljet_phi0[i0], ljet_e0[i0]);
      
      int btag = 0;
      for (int i1=0; i1<n_jet0; i1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt0[i1], jet_eta0[i1], jet_phi0[i1], jet_e0[i1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3 && jet_btagged0[i1]==0) count04++;
      }
      for (int i2=0; i2<n_tjet0; i2++)
      {
	TLorentzVector tjet;
	tjet.SetPtEtaPhiE(tjet_pt0[i2], tjet_eta0[i2], tjet_phi0[i2], tjet_e0[i2]);
	float delta_r;
	delta_r = tjet.DeltaR(ljet);
	if (delta_r<0.2 && tjet_btagged0[i2] > 0) btag++;
      }

      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) 
      {
	if (btag == 0) count01a++;
	else if (btag == 1) count01b++;
	else count01c++;
      }
      if (ljet.M()>145.0e3 && btag==0) count02++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3 && btag==0) count03++;
      
    }
    if (count01a > 0 && count02 > 0) nsel01a++;
    else if (count01a > 0 && count03 > 1) nsel02a++;
    else if (count01a > 0 && count03 > 0 && count04 > 0) nsel03a++;

    if (count01b > 0 && count02 > 0) nsel01b++;
    else if (count01b > 0 && count03 > 1) nsel02b++;
    else if (count01b > 0 && count03 > 0 && count04 > 0) nsel03b++;
    
    if (count01c > 0 && count02 > 0) nsel01c++;
    else if (count01c > 0 && count03 > 1) nsel02c++;
    else if (count01c > 0 && count03 > 0 && count04 > 0) nsel03c++;
  }
  cout << "nsel01: " << nsel01a << "\t" << nsel01b << "\t" << nsel01c << endl;
  cout << "nsel02: " << nsel02a << "\t" << nsel02b << "\t" << nsel02c << endl;
  cout << "nsel03: " << nsel03a << "\t" << nsel03b << "\t" << nsel03c << endl;
  file0->Close();

  
//////////////////////////////////X 2000 S 1500 //////////////////////////
  
  auto file1 = TFile::Open("bbww_x2000_s1500.root");
  TTree* t1 = (TTree*)file1->Get("CollectionTree");
  
  Int_t           n_jet1;
  Float_t         jet_pt1[14];		//[n_jet]
  Int_t           jet_btagged1[14];	//[n_jet]
  Float_t         jet_eta1[14];		//[n_jet]
  Float_t         jet_phi1[14];		//[n_jet] 
  Float_t         jet_e1[14];		//[n_jet]
  Int_t           n_ljet1;
  Float_t         ljet_pt1[9];   	//[n_ljet]
  Float_t         ljet_eta1[9];   	//[n_ljet]
  Float_t         ljet_phi1[9];   	//[n_ljet]
  Float_t         ljet_e1[9];   	//[n_ljet]
  Int_t           n_tjet1;
  Float_t         tjet_pt1[27];   //[n_tjet]
  Float_t         tjet_eta1[27];   //[n_tjet]
  Float_t         tjet_phi1[27];   //[n_tjet]
  Float_t         tjet_e1[27];   //[n_tjet]
  Int_t           tjet_btagged1[27];   //[n_tjet]  
  
  t1->SetBranchAddress("n_jet", &n_jet1);
  t1->SetBranchAddress("jet_pt", jet_pt1);
  t1->SetBranchAddress("jet_btagged",jet_btagged1);
  t1->SetBranchAddress("jet_eta", jet_eta1);
  t1->SetBranchAddress("jet_phi", jet_phi1);
  t1->SetBranchAddress("jet_e", jet_e1);
  t1->SetBranchAddress("n_ljet", &n_ljet1);
  t1->SetBranchAddress("ljet_pt", ljet_pt1);
  t1->SetBranchAddress("ljet_eta", ljet_eta1);
  t1->SetBranchAddress("ljet_phi", ljet_phi1);
  t1->SetBranchAddress("ljet_e", ljet_e1);
  t1->SetBranchAddress("n_tjet", &n_tjet1);
  t1->SetBranchAddress("tjet_pt", tjet_pt1);
  t1->SetBranchAddress("tjet_eta", tjet_eta1);
  t1->SetBranchAddress("tjet_phi", tjet_phi1);
  t1->SetBranchAddress("tjet_e", tjet_e1);
  t1->SetBranchAddress("tjet_btagged", tjet_btagged1);  
  
  t1->SetBranchStatus("*",0);
  
  t1->SetBranchStatus("n_jet",1);
  t1->SetBranchStatus("jet_pt",1);
  t1->SetBranchStatus("jet_btagged",1);
  t1->SetBranchStatus("jet_eta",1);
  t1->SetBranchStatus("jet_phi",1);
  t1->SetBranchStatus("jet_e", 1);
  t1->SetBranchStatus("n_ljet",1);
  t1->SetBranchStatus("ljet_pt", 1);
  t1->SetBranchStatus("ljet_eta", 1);
  t1->SetBranchStatus("ljet_phi", 1);
  t1->SetBranchStatus("ljet_e", 1);
  t1->SetBranchStatus("n_tjet", 1);
  t1->SetBranchStatus("tjet_pt", 1);
  t1->SetBranchStatus("tjet_eta", 1);
  t1->SetBranchStatus("tjet_phi", 1);
  t1->SetBranchStatus("tjet_e", 1);
  t1->SetBranchStatus("tjet_btagged", 1);
  
  int nentries1 = t1->GetEntries();
  int nsel11a = 0;		//counting fully boosted events that are selected
  int nsel11b = 0;
  int nsel11c = 0;
  int nsel12a = 0;		//counting boosted events that are selected
  int nsel12b = 0;
  int nsel12c = 0;
  int nsel13a = 0;		//counting semi-boosted events that are selected
  int nsel13b = 0;
  int nsel13c = 0;
  
  cout << "\nnentries1: " << nentries1 << "\tX: 2000 GeV \tS: 1500 Gev" << endl;
  for (int j=0; j<nentries1; j++)			//loops over all the 20,000 entries in file1
  {
    t1->GetEntry(j);
    int count11a = 0;
    int count11b = 0;
    int count11c = 0;
    int count12 = 0;
    int count13 = 0;
    int count14 = 0;
    
    for ( int j0=0; j0<n_ljet1; j0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt1[j0], ljet_eta1[j0], ljet_phi1[j0], ljet_e1[j0]);
      
      int btag = 0;
      for (int j1=0; j1<n_jet1; j1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt1[j1], jet_eta1[j1], jet_phi1[j1], jet_e1[j1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3 && jet_btagged1[j1]==0) count14++;
      }
      for (int j2=0; j2<n_tjet1; j2++)
      {
	TLorentzVector tjet;
	tjet.SetPtEtaPhiE(tjet_pt1[j2], tjet_eta1[j2], tjet_phi1[j2], tjet_e1[j2]);
	float delta_r;
	delta_r = tjet.DeltaR(ljet);
	if (delta_r<0.2 && tjet_btagged1[j2] > 0) btag++;
      }      
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3)
      {
	if (btag == 0) count11a++;
	else if (btag == 1) count11b++;
	else count11c++;
      }
      if (ljet.M()>145.0e3 && btag==0) count12++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3 && btag==0) count13++;
    }
    if (count11a > 0 && count12 > 0) nsel11a++;
    else if (count11a > 0 && count13 > 1) nsel12a++;
    else if (count11a > 0 && count13 > 0 && count14 > 0) nsel13a++;
    
    if (count11b > 0 && count12 > 0) nsel11b++;
    else if (count11b > 0 && count13 > 1) nsel12b++;
    else if (count11b > 0 && count13 > 0 && count14 > 0) nsel13b++;
    
    if (count11c > 0 && count12 > 0) nsel11c++;
    else if (count11c > 0 && count13 > 1) nsel12c++;
    else if (count11c > 0 && count13 > 0 && count14 > 0) nsel13c++;

  }
  float nsel_11a = nsel11a;
  float nsel_11b = nsel11b;
  float nsel_11c = nsel11c;
  float nsel_12a = nsel12a;
  float nsel_12b = nsel12b;
  float nsel_12c = nsel12c;
  float nsel_13a = nsel13a;
  float nsel_13b = nsel13b;
  float nsel_13c = nsel13c;
  
  cout << "nsel11: " << nsel_11a << "\t" << nsel_11b << "\t" << nsel_11c << endl;
  cout << "nsel12: " << nsel_12a << "\t" << nsel_12b << "\t" << nsel_12c << endl;
  cout << "nsel13: " << nsel_13a << "\t" << nsel_13b << "\t" << nsel_13c << endl;
  
  smass[5] = 1500.0;
  fraction1a[5] = nsel_11a/nentries1;
  fraction2a[5] = nsel_12a/nentries1;
  fraction3a[5] = nsel_13a/nentries1;
  
  fraction1b[5] = nsel_11b/nentries1;
  fraction2b[5] = nsel_12b/nentries1;
  fraction3b[5] = nsel_13b/nentries1;
  
  fraction1c[5] = nsel_11c/nentries1;
  fraction2c[5] = nsel_12c/nentries1;
  fraction3c[5] = nsel_13c/nentries1;
  
  fraction1aerror[5] = sqrt(fraction1a[5] * (1-fraction1a[5])) / sqrt(nentries1);
  fraction2aerror[5] = sqrt(fraction2a[5] * (1-fraction2a[5])) / sqrt(nentries1);
  fraction3aerror[5] = sqrt(fraction3a[5] * (1-fraction3a[5])) / sqrt(nentries1);
  
  fraction1berror[5] = sqrt(fraction1b[5] * (1-fraction1b[5])) / sqrt(nentries1);
  fraction2berror[5] = sqrt(fraction2b[5] * (1-fraction2b[5])) / sqrt(nentries1);
  fraction3berror[5] = sqrt(fraction3b[5] * (1-fraction3b[5])) / sqrt(nentries1);
  
  fraction1cerror[5] = sqrt(fraction1c[5] * (1-fraction1c[5])) / sqrt(nentries1);
  fraction2cerror[5] = sqrt(fraction2c[5] * (1-fraction2c[5])) / sqrt(nentries1);
  fraction3cerror[5] = sqrt(fraction3c[5] * (1-fraction3c[5])) / sqrt(nentries1);
  
  cout << fraction1a[5] << "\t" << fraction1b[5] << "\t" << fraction1c[5] << endl;
  cout << fraction2a[5] << "\t" << fraction2b[5] << "\t" << fraction1c[5] << endl;
  cout << fraction3a[5] << "\t" << fraction3b[5] << "\t" << fraction3c[5]<< endl;
  
  file1->Close();

////////////////////////////////////X 2000 S 170 ///////////////////////////
  
  auto file2 = TFile::Open("bbww_x2000_s170.root");
  TTree* t2 = (TTree*)file2->Get("CollectionTree");
  
  Int_t           n_jet2;
  Float_t         jet_pt2[14];		//[n_jet]
  Int_t           jet_btagged2[14];	//[n_jet]
  Float_t         jet_eta2[14];		//[n_jet]
  Float_t         jet_phi2[14];		//[n_jet] 
  Float_t         jet_e2[14];		//[n_jet]
  Int_t           n_ljet2;
  Float_t         ljet_pt2[9];   	//[n_ljet]
  Float_t         ljet_eta2[9];   	//[n_ljet]
  Float_t         ljet_phi2[9];   	//[n_ljet]
  Float_t         ljet_e2[9];   	//[n_ljet]
  Int_t           n_tjet2;
  Float_t         tjet_pt2[26];   //[n_tjet]
  Float_t         tjet_eta2[26];   //[n_tjet]
  Float_t         tjet_phi2[26];   //[n_tjet]
  Float_t         tjet_e2[26];   //[n_tjet]
  Int_t           tjet_btagged2[26];   //[n_tjet]
  
  t2->SetBranchAddress("n_jet", &n_jet2);
  t2->SetBranchAddress("jet_pt", jet_pt2);
  t2->SetBranchAddress("jet_btagged",jet_btagged2);
  t2->SetBranchAddress("jet_eta", jet_eta2);
  t2->SetBranchAddress("jet_phi", jet_phi2);
  t2->SetBranchAddress("jet_e", jet_e2);
  t2->SetBranchAddress("n_ljet", &n_ljet2);
  t2->SetBranchAddress("ljet_pt", ljet_pt2);
  t2->SetBranchAddress("ljet_eta", ljet_eta2);
  t2->SetBranchAddress("ljet_phi", ljet_phi2);
  t2->SetBranchAddress("ljet_e", ljet_e2);
  t2->SetBranchAddress("n_tjet", &n_tjet2);
  t2->SetBranchAddress("tjet_pt", tjet_pt2);
  t2->SetBranchAddress("tjet_eta", tjet_eta2);
  t2->SetBranchAddress("tjet_phi", tjet_phi2);
  t2->SetBranchAddress("tjet_e", tjet_e2);
  t2->SetBranchAddress("tjet_btagged", tjet_btagged2); 
  
  t2->SetBranchStatus("*",0);
  
  t2->SetBranchStatus("n_jet",1);
  t2->SetBranchStatus("jet_pt",1);
  t2->SetBranchStatus("jet_btagged",1);
  t2->SetBranchStatus("jet_eta",1);
  t2->SetBranchStatus("jet_phi",1);
  t2->SetBranchStatus("jet_e", 1);
  t2->SetBranchStatus("n_ljet",1);
  t2->SetBranchStatus("ljet_pt", 1);
  t2->SetBranchStatus("ljet_eta", 1);
  t2->SetBranchStatus("ljet_phi", 1);
  t2->SetBranchStatus("ljet_e", 1);
  t2->SetBranchStatus("n_tjet", 1);
  t2->SetBranchStatus("tjet_pt", 1);
  t2->SetBranchStatus("tjet_eta", 1);
  t2->SetBranchStatus("tjet_phi", 1);
  t2->SetBranchStatus("tjet_e", 1);
  t2->SetBranchStatus("tjet_btagged", 1);
  
  int nentries2 = t2->GetEntries();
  int nsel21a = 0;		//counting fully boosted events that are selected
  int nsel21b = 0;
  int nsel21c = 0;
  int nsel22a = 0;		//counting boosted events that are selected
  int nsel22b = 0;
  int nsel22c = 0;
  int nsel23a = 0;		//counting semi-boosted events that are selected
  int nsel23b = 0;
  int nsel23c = 0;
  
  cout << "\nnentries2: " << nentries2 << "\tX: 2000 GeV \tS: 170 Gev" << endl;
  for (int k=0; k<nentries2; k++)			//loops over all the 20,000 entries in file2
  {
    t2->GetEntry(k);
    int count21a = 0;
    int count21b = 0;
    int count21c = 0;
    int count22 = 0;
    int count23 = 0;
    int count24 = 0;
    for ( int k0=0; k0<n_ljet2; k0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt2[k0], ljet_eta2[k0], ljet_phi2[k0], ljet_e2[k0]);
      
      int btag = 0;
      for (int k1=0; k1<n_jet2; k1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt2[k1], jet_eta2[k1], jet_phi2[k1], jet_e2[k1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3 && jet_btagged2[k1]==0) count24++;
      }
      for (int k2=0; k2<n_tjet2; k2++)
      {
	TLorentzVector tjet;
	tjet.SetPtEtaPhiE(tjet_pt2[k2], tjet_eta2[k2], tjet_phi2[k2], tjet_e2[k2]);
	float delta_r;
	delta_r = tjet.DeltaR(ljet);
	if (delta_r<0.2 && tjet_btagged2[k2] > 0) btag++;
      }      
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3)
      {
	if (btag == 0) count21a++;
	else if (btag == 1) count21b++;
	else count21c++;
      }
      if (ljet.M()>145.0e3 && btag == 0) count22++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3 && btag == 0) count23++;
    }
    if (count21a > 0 && count22 > 0) nsel21a++;
    else if (count21a > 0 && count23 > 1) nsel22a++;
    else if (count21a > 0 && count23 > 0 && count24 > 0) nsel23a++;
    
    if (count21b > 0 && count22 > 0) nsel21b++;
    else if (count21b > 0 && count23 > 1) nsel22b++;
    else if (count21b > 0 && count23 > 0 && count24 > 0) nsel23b++;
    
    if (count21c > 0 && count22 > 0) nsel21c++;
    else if (count21c > 0 && count23 > 1) nsel22c++;
    else if (count21c > 0 && count23 > 0 && count24 > 0) nsel23c++;
  }
  float nsel_21a = nsel21a;
  float nsel_22a = nsel22a;
  float nsel_23a = nsel23a;
  
  float nsel_21b = nsel21b;
  float nsel_22b = nsel22b;
  float nsel_23b = nsel23b;
  
  float nsel_21c = nsel21c;
  float nsel_22c = nsel22c;
  float nsel_23c = nsel23c;
    
  cout << "nsel21: " << nsel_21a << "\t" << nsel_21b << "\t" << nsel_21c << endl;
  cout << "nsel22: " << nsel_22a << "\t" << nsel_22b << "\t" << nsel_22c << endl;
  cout << "nsel23: " << nsel_23a << "\t" << nsel_23b << "\t" << nsel_23c << endl;
  
  smass[0] = 170.0;
  fraction1a[0] = nsel_21a/nentries2;
  fraction2a[0] = nsel_22a/nentries2;
  fraction3a[0] = nsel_23a/nentries2;
  
  fraction1b[0] = nsel_21b/nentries2;
  fraction2b[0] = nsel_22b/nentries2;
  fraction3b[0] = nsel_23b/nentries2;
  
  fraction1c[0] = nsel_21c/nentries2;
  fraction2c[0] = nsel_22c/nentries2;
  fraction3c[0] = nsel_23c/nentries2;
  
  fraction1aerror[0] = sqrt(fraction1a[0] * (1-fraction1a[0])) / sqrt(nentries2);
  fraction2aerror[0] = sqrt(fraction2a[0] * (1-fraction2a[0])) / sqrt(nentries2);
  fraction3aerror[0] = sqrt(fraction3a[0] * (1-fraction3a[0])) / sqrt(nentries2);
  
  fraction1berror[0] = sqrt(fraction1b[0] * (1-fraction1b[0])) / sqrt(nentries2);
  fraction2berror[0] = sqrt(fraction2b[0] * (1-fraction2b[0])) / sqrt(nentries2);
  fraction3berror[0] = sqrt(fraction3b[0] * (1-fraction3b[0])) / sqrt(nentries2);
  
  fraction1cerror[0] = sqrt(fraction1c[0] * (1-fraction1c[0])) / sqrt(nentries2);
  fraction2cerror[0] = sqrt(fraction2c[0] * (1-fraction2c[0])) / sqrt(nentries2);
  fraction3cerror[0] = sqrt(fraction3c[0] * (1-fraction3c[0])) / sqrt(nentries2);
  
  cout << fraction1a[0] << "\t" << fraction1b[0] << "\t" << fraction1c[0] << endl;
  cout << fraction2a[0] << "\t" << fraction2b[0] << "\t" << fraction1c[0] << endl;
  cout << fraction3a[0] << "\t" << fraction3b[0] << "\t" << fraction3c[0]<< endl;
  
  file2->Close();
//////////////////////////////////X 2000 S 240///////////////
  
  auto file3 = TFile::Open("bbww_x2000_s240.root");
  TTree* t3 = (TTree*)file3->Get("CollectionTree");
  
  Int_t           n_jet3;
  Float_t         jet_pt3[15];		//[n_jet]
  Int_t           jet_btagged3[15];	//[n_jet]
  Float_t         jet_eta3[15];		//[n_jet]
  Float_t         jet_phi3[15];		//[n_jet] 
  Float_t         jet_e3[15];		//[n_jet]
  Int_t           n_ljet3;
  Float_t         ljet_pt3[10];   	//[n_ljet]
  Float_t         ljet_eta3[10];   	//[n_ljet]
  Float_t         ljet_phi3[10];   	//[n_ljet]
  Float_t         ljet_e3[10];   	//[n_ljet]
  Int_t           n_tjet3;
  Float_t         tjet_pt3[32];   //[n_tjet]
  Float_t         tjet_eta3[32];   //[n_tjet]
  Float_t         tjet_phi3[32];   //[n_tjet]
  Float_t         tjet_e3[32];   //[n_tjet]
  Int_t           tjet_btagged3[32];   //[n_tjet]
  
  t3->SetBranchAddress("n_jet", &n_jet3);
  t3->SetBranchAddress("jet_pt", jet_pt3);
  t3->SetBranchAddress("jet_btagged",jet_btagged3);
  t3->SetBranchAddress("jet_eta", jet_eta3);
  t3->SetBranchAddress("jet_phi", jet_phi3);
  t3->SetBranchAddress("jet_e", jet_e3);
  t3->SetBranchAddress("n_ljet", &n_ljet3);
  t3->SetBranchAddress("ljet_pt", ljet_pt3);
  t3->SetBranchAddress("ljet_eta", ljet_eta3);
  t3->SetBranchAddress("ljet_phi", ljet_phi3);
  t3->SetBranchAddress("ljet_e", ljet_e3);
  t3->SetBranchAddress("n_tjet", &n_tjet3);
  t3->SetBranchAddress("tjet_pt", tjet_pt3);
  t3->SetBranchAddress("tjet_eta", tjet_eta3);
  t3->SetBranchAddress("tjet_phi", tjet_phi3);
  t3->SetBranchAddress("tjet_e", tjet_e3);
  t3->SetBranchAddress("tjet_btagged", tjet_btagged3);
  
  t3->SetBranchStatus("*",0);
  
  t3->SetBranchStatus("n_jet",1);
  t3->SetBranchStatus("jet_pt",1);
  t3->SetBranchStatus("jet_btagged",1);
  t3->SetBranchStatus("jet_eta",1);
  t3->SetBranchStatus("jet_phi",1);
  t3->SetBranchStatus("jet_e", 1);
  t3->SetBranchStatus("n_ljet",1);
  t3->SetBranchStatus("ljet_pt", 1);
  t3->SetBranchStatus("ljet_eta", 1);
  t3->SetBranchStatus("ljet_phi", 1);
  t3->SetBranchStatus("ljet_e", 1);
  t3->SetBranchStatus("n_tjet", 1);
  t3->SetBranchStatus("tjet_pt", 1);
  t3->SetBranchStatus("tjet_eta", 1);
  t3->SetBranchStatus("tjet_phi", 1);
  t3->SetBranchStatus("tjet_e", 1);
  t3->SetBranchStatus("tjet_btagged", 1);
  
  int nentries3 = t3->GetEntries();
  int nsel31a = 0;		//counting fully boosted events that are selected
  int nsel31b = 0;
  int nsel31c = 0;
  int nsel32a = 0;		//counting boosted events that are selected
  int nsel32b = 0;
  int nsel32c = 0;
  int nsel33a = 0;		//counting semi-boosted events that are selected
  int nsel33b = 0;
  int nsel33c = 0;
  
  cout << "\nnentries3: " << nentries3 << "\tX: 2000 GeV \tS: 240 Gev" << endl;
  for (int l=0; l<nentries3; l++)			//loops over all the 20,000 entries in file3
  {
    t3->GetEntry(l);
    int count31a = 0;
    int count31b = 0;
    int count31c = 0;
    int count32 = 0;
    int count33 = 0;
    int count34 = 0;
    
    for ( int l0=0; l0<n_ljet3; l0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt3[l0], ljet_eta3[l0], ljet_phi3[l0], ljet_e3[l0]);
      
      int btag = 0;
      for (int l1=0; l1<n_jet3; l1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt3[l1], jet_eta3[l1], jet_phi3[l1], jet_e3[l1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count34++;
      }
      for (int l2=0; l2<n_tjet3; l2++)
      {
	TLorentzVector tjet;
	tjet.SetPtEtaPhiE(tjet_pt3[l2], tjet_eta3[l2], tjet_phi3[l2], tjet_e3[l2]);
	float delta_r;
	delta_r = tjet.DeltaR(ljet);
	if (delta_r<0.2 && tjet_btagged3[l2] > 0) btag++;
      }      
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) 
      {
	if (btag == 0) count31a++;
	else if (btag == 1) count31b++;
	else count31c++;
      }
      if (ljet.M()>145.0e3 && btag==0) count32++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3 && btag==0) count33++;
    }
    if (count31a > 0 && count32 > 0) nsel31a++;
    else if (count31a > 0 && count33 > 1) nsel32a++;
    else if (count31a > 0 && count33 > 0 && count34 > 0) nsel33a++;
    
    if (count31b > 0 && count32 > 0) nsel31b++;
    else if (count31b > 0 && count33 > 1) nsel32b++;
    else if (count31b > 0 && count33 > 0 && count34 > 0) nsel33b++;
    
    if (count31c > 0 && count32 > 0) nsel31c++;
    else if (count31c > 0 && count33 > 1) nsel32c++;
    else if (count31c > 0 && count33 > 0 && count34 > 0) nsel33c++;
  }
  float nsel_31a = nsel31a;
  float nsel_32a = nsel32a;
  float nsel_33a = nsel33a;
  
  float nsel_31b = nsel31b;
  float nsel_32b = nsel32b;
  float nsel_33b = nsel33b;
  
  float nsel_31c = nsel31c;
  float nsel_32c = nsel32c;
  float nsel_33c = nsel33c;
  
  cout << "nsel31: " << nsel_31a << "\t" << nsel_31b << "\t" << nsel_31c << endl;
  cout << "nsel32: " << nsel_32a << "\t" << nsel_32b << "\t" << nsel_32c << endl;
  cout << "nsel33: " << nsel_33a << "\t" << nsel_33b << "\t" << nsel_33c << endl;
  
  smass[1] = 240.0;
  fraction1a[1] = nsel_21a/nentries3;
  fraction2a[1] = nsel_22a/nentries3;
  fraction3a[1] = nsel_23a/nentries3;
  
  fraction1b[1] = nsel_21b/nentries3;
  fraction2b[1] = nsel_22b/nentries3;
  fraction3b[1] = nsel_23b/nentries3;
  
  fraction1c[1] = nsel_21c/nentries3;
  fraction2c[1] = nsel_22c/nentries3;
  fraction3c[1] = nsel_23c/nentries3;
  
  fraction1aerror[1] = sqrt(fraction1a[1] * (1-fraction1a[1])) / sqrt(nentries3);
  fraction2aerror[1] = sqrt(fraction2a[1] * (1-fraction2a[1])) / sqrt(nentries3);
  fraction3aerror[1] = sqrt(fraction3a[1] * (1-fraction3a[1])) / sqrt(nentries3);
  
  fraction1berror[1] = sqrt(fraction1b[1] * (1-fraction1b[1])) / sqrt(nentries3);
  fraction2berror[1] = sqrt(fraction2b[1] * (1-fraction2b[1])) / sqrt(nentries3);
  fraction3berror[1] = sqrt(fraction3b[1] * (1-fraction3b[1])) / sqrt(nentries3);
  
  fraction1cerror[1] = sqrt(fraction1c[1] * (1-fraction1c[1])) / sqrt(nentries3);
  fraction2cerror[1] = sqrt(fraction2c[1] * (1-fraction2c[1])) / sqrt(nentries3);
  fraction3cerror[1] = sqrt(fraction3c[1] * (1-fraction3c[1])) / sqrt(nentries3);
  
  cout << fraction1a[1] << "\t" << fraction1b[1] << "\t" << fraction1c[1] << endl;
  cout << fraction2a[1] << "\t" << fraction2b[1] << "\t" << fraction1c[1] << endl;
  cout << fraction3a[1] << "\t" << fraction3b[1] << "\t" << fraction3c[1]<< endl;
  
  file3->Close();
  
///////////////////////////////X 2000 S 400//////////////////////////////////////////////////
  
  auto file4 = TFile::Open("bbww_x2000_s400.root");
  TTree* t4 = (TTree*)file4->Get("CollectionTree");
  
  Int_t           n_jet4;
  Float_t         jet_pt4[15];		//[n_jet]
  Int_t           jet_btagged4[15];	//[n_jet]
  Float_t         jet_eta4[15];		//[n_jet]
  Float_t         jet_phi4[15];		//[n_jet] 
  Float_t         jet_e4[15];		//[n_jet]
  Int_t           n_ljet4;
  Float_t         ljet_pt4[10];   	//[n_ljet]
  Float_t         ljet_eta4[10];   	//[n_ljet]
  Float_t         ljet_phi4[10];   	//[n_ljet]
  Float_t         ljet_e4[10];   		//[n_ljet]
  Int_t           n_tjet4;
  Float_t         tjet_pt4[26];   //[n_tjet]
  Float_t         tjet_eta4[26];   //[n_tjet]
  Float_t         tjet_phi4[26];   //[n_tjet]
  Float_t         tjet_e4[26];   //[n_tjet]
  Int_t           tjet_btagged4[26];   //[n_tjet]
  
  t4->SetBranchAddress("n_jet", &n_jet4);
  t4->SetBranchAddress("jet_pt", jet_pt4);
  t4->SetBranchAddress("jet_btagged",jet_btagged4);
  t4->SetBranchAddress("jet_eta", jet_eta4);
  t4->SetBranchAddress("jet_phi", jet_phi4);
  t4->SetBranchAddress("jet_e", jet_e4);
  t4->SetBranchAddress("n_ljet", &n_ljet4);
  t4->SetBranchAddress("ljet_pt", ljet_pt4);
  t4->SetBranchAddress("ljet_eta", ljet_eta4);
  t4->SetBranchAddress("ljet_phi", ljet_phi4);
  t4->SetBranchAddress("ljet_e", ljet_e4);
  t4->SetBranchAddress("n_tjet", &n_tjet4);
  t4->SetBranchAddress("tjet_pt", tjet_pt4);
  t4->SetBranchAddress("tjet_eta", tjet_eta4);
  t4->SetBranchAddress("tjet_phi", tjet_phi4);
  t4->SetBranchAddress("tjet_e", tjet_e4);
  t4->SetBranchAddress("tjet_btagged", tjet_btagged4);
  
  t4->SetBranchStatus("*",0);
  
  t4->SetBranchStatus("n_jet",1);
  t4->SetBranchStatus("jet_pt",1);
  t4->SetBranchStatus("jet_btagged",1);
  t4->SetBranchStatus("jet_eta",1);
  t4->SetBranchStatus("jet_phi",1);
  t4->SetBranchStatus("jet_e", 1);
  t4->SetBranchStatus("n_ljet",1);
  t4->SetBranchStatus("ljet_pt", 1);
  t4->SetBranchStatus("ljet_eta", 1);
  t4->SetBranchStatus("ljet_phi", 1);
  t4->SetBranchStatus("ljet_e", 1);
  t4->SetBranchStatus("n_tjet", 1);
  t4->SetBranchStatus("tjet_pt", 1);
  t4->SetBranchStatus("tjet_eta", 1);
  t4->SetBranchStatus("tjet_phi", 1);
  t4->SetBranchStatus("tjet_e", 1);
  t4->SetBranchStatus("tjet_btagged", 1);
  
  int nentries4 = t4->GetEntries();
  int nsel41a = 0;		//counting fully boosted events that are selected
  int nsel41b = 0;
  int nsel41c = 0;
  int nsel42a = 0;		//counting boosted events that are selected
  int nsel42b = 0;
  int nsel42c = 0;
  int nsel43a = 0;		//counting semi-boosted events that are selected
  int nsel43b = 0;
  int nsel43c = 0;
  
  cout << "\nnentries4: " << nentries4 << "\tX: 2000 GeV \tS: 400 Gev" << endl;
  for (int m=0; m<nentries4; m++)			//loops over all the 20,000 entries in file4
  {
    t4->GetEntry(m);
    int count41a = 0;
    int count41b = 0;
    int count41c = 0;
    int count42 = 0;
    int count43 = 0;
    int count44 = 0;
    for ( int m0=0; m0<n_ljet4; m0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt4[m0], ljet_eta4[m0], ljet_phi4[m0], ljet_e4[m0]);

      int btag = 0;
      for (int m1=0; m1<n_jet4; m1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt4[m1], jet_eta4[m1], jet_phi4[m1], jet_e4[m1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count44++;
      }
      for (int m2=0; m2<n_tjet4; m2++)
      {
	TLorentzVector tjet;
	tjet.SetPtEtaPhiE(tjet_pt4[m2], tjet_eta4[m2], tjet_phi4[m2], tjet_e4[m2]);
	float delta_r;
	delta_r = tjet.DeltaR(ljet);
	if (delta_r<0.2 && tjet_btagged4[m2] > 0) btag++;
      }      
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3)
      {
	if (btag == 0) count41a++;
	else if (btag == 1) count41b++;
	else count41c++;
      }
      if (ljet.M()>145.0e3 && btag==0) count42++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3 && btag==0) count43++;      
    }
    
    if (count41a > 0 && count42 > 0) nsel41a++;
    else if (count41a > 0 && count43 > 1) nsel42a++;
    else if (count41a > 0 && count43 > 0 && count44 > 0) nsel43a++;
    
    if (count41b > 0 && count42 > 0) nsel41b++;
    else if (count41b > 0 && count43 > 1) nsel42b++;
    else if (count41b > 0 && count43 > 0 && count44 > 0) nsel43b++;
    
    if (count41c > 0 && count42 > 0) nsel41c++;
    else if (count41c > 0 && count43 > 1) nsel42c++;
    else if (count41c > 0 && count43 > 0 && count44 > 0) nsel43c++;
  }
  float nsel_41a = nsel41a;
  float nsel_42a = nsel42a;
  float nsel_43a = nsel43a;

  float nsel_41b = nsel41b;
  float nsel_42b = nsel42b;
  float nsel_43b = nsel43b;
  
  float nsel_41c = nsel41c;
  float nsel_42c = nsel42c;
  float nsel_43c = nsel43c;
  
  cout << "nsel41: " << nsel_41a << "\t" << nsel_41b << "\t" << nsel_41c << endl;
  cout << "nsel42: " << nsel_42a << "\t" << nsel_42b << "\t" << nsel_42c << endl;
  cout << "nsel43: " << nsel_43a << "\t" << nsel_43b << "\t" << nsel_43c << endl;
  
  smass[2] = 400.0;
  fraction1a[2] = nsel_41a/nentries4;
  fraction2a[2] = nsel_42a/nentries4;
  fraction3a[2] = nsel_43a/nentries4;
  
  fraction1b[2] = nsel_41b/nentries4;
  fraction2b[2] = nsel_42b/nentries4;
  fraction3b[2] = nsel_43b/nentries4;
  
  fraction1c[2] = nsel_41c/nentries4;
  fraction2c[2] = nsel_42c/nentries4;
  fraction3c[2] = nsel_43c/nentries4;
  
  fraction1aerror[2] = sqrt(fraction1a[2] * (1-fraction1a[2])) / sqrt(nentries4);
  fraction2aerror[2] = sqrt(fraction2a[2] * (1-fraction2a[2])) / sqrt(nentries4);
  fraction3aerror[2] = sqrt(fraction3a[2] * (1-fraction3a[2])) / sqrt(nentries4);

  fraction1berror[2] = sqrt(fraction1b[2] * (1-fraction1b[2])) / sqrt(nentries4);
  fraction2berror[2] = sqrt(fraction2b[2] * (1-fraction2b[2])) / sqrt(nentries4);
  fraction3berror[2] = sqrt(fraction3b[2] * (1-fraction3b[2])) / sqrt(nentries4);
  
  fraction1cerror[2] = sqrt(fraction1c[2] * (1-fraction1c[2])) / sqrt(nentries4);
  fraction2cerror[2] = sqrt(fraction2c[2] * (1-fraction2c[2])) / sqrt(nentries4);
  fraction3cerror[2] = sqrt(fraction3c[2] * (1-fraction3c[2])) / sqrt(nentries4);
  
  cout << fraction1a[2] << "\t" << fraction1b[2] << "\t" << fraction1c[2] << endl;
  cout << fraction2a[2] << "\t" << fraction2b[2] << "\t" << fraction1c[2] << endl;
  cout << fraction3a[2] << "\t" << fraction3b[2] << "\t" << fraction3c[2]<< endl;
  
  file4->Close();
  
//////////////////////////X 2000 S 750//////////////////////////
  
  auto file5 = TFile::Open("bbww_x2000_s750.root");
  TTree* t5 = (TTree*)file5->Get("CollectionTree");
  
  Int_t           n_jet5;
  Float_t         jet_pt5[14];		//[n_jet]
  Int_t           jet_btagged5[14];	//[n_jet]
  Float_t         jet_eta5[14];		//[n_jet]
  Float_t         jet_phi5[14];		//[n_jet] 
  Float_t         jet_e5[14];		//[n_jet]
  Int_t           n_ljet5;
  Float_t         ljet_pt5[10];   	//[n_ljet]
  Float_t         ljet_eta5[10];   	//[n_ljet]
  Float_t         ljet_phi5[10];   	//[n_ljet]
  Float_t         ljet_e5[10];   		//[n_ljet]
  Int_t           n_tjet5;
  Float_t         tjet_pt5[28];   //[n_tjet]
  Float_t         tjet_eta5[28];   //[n_tjet]
  Float_t         tjet_phi5[28];   //[n_tjet]
  Float_t         tjet_e5[28];   //[n_tjet]
  Int_t           tjet_btagged5[28];   //[n_tjet]
  
  t5->SetBranchAddress("n_jet", &n_jet5);
  t5->SetBranchAddress("jet_pt", jet_pt5);
  t5->SetBranchAddress("jet_btagged",jet_btagged5);
  t5->SetBranchAddress("jet_eta", jet_eta5);
  t5->SetBranchAddress("jet_phi", jet_phi5);
  t5->SetBranchAddress("jet_e", jet_e5);
  t5->SetBranchAddress("n_ljet", &n_ljet5);
  t5->SetBranchAddress("ljet_pt", ljet_pt5);
  t5->SetBranchAddress("ljet_eta", ljet_eta5);
  t5->SetBranchAddress("ljet_phi", ljet_phi5);
  t5->SetBranchAddress("ljet_e", ljet_e5);
  t5->SetBranchAddress("n_tjet", &n_tjet5);
  t5->SetBranchAddress("tjet_pt", tjet_pt5);
  t5->SetBranchAddress("tjet_eta", tjet_eta5);
  t5->SetBranchAddress("tjet_phi", tjet_phi5);
  t5->SetBranchAddress("tjet_e", tjet_e5);
  t5->SetBranchAddress("tjet_btagged", tjet_btagged5);
  
  t5->SetBranchStatus("*",0);
  
  t5->SetBranchStatus("n_jet",1);
  t5->SetBranchStatus("jet_pt",1);
  t5->SetBranchStatus("jet_btagged",1);
  t5->SetBranchStatus("jet_eta",1);
  t5->SetBranchStatus("jet_phi",1);
  t5->SetBranchStatus("jet_e", 1);
  t5->SetBranchStatus("n_ljet",1);
  t5->SetBranchStatus("ljet_pt", 1);
  t5->SetBranchStatus("ljet_eta", 1);
  t5->SetBranchStatus("ljet_phi", 1);
  t5->SetBranchStatus("ljet_e", 1);
  t5->SetBranchStatus("n_tjet", 1);
  t5->SetBranchStatus("tjet_pt", 1);
  t5->SetBranchStatus("tjet_eta", 1);
  t5->SetBranchStatus("tjet_phi", 1);
  t5->SetBranchStatus("tjet_e", 1);
  t5->SetBranchStatus("tjet_btagged", 1);
  
  int nentries5 = t5->GetEntries();
  int nsel51a = 0;		//counting fully boosted events that are selected
  int nsel51b = 0;
  int nsel51c = 0;
  int nsel52a = 0;		//counting boosted events that are selected
  int nsel52b = 0;
  int nsel52c = 0;
  int nsel53a = 0;		//counting semi-boosted events that are selected
  int nsel53b = 0;
  int nsel53c = 0;
  cout << "\nnentries5: " << nentries5 << "\tX: 2000 GeV \tS: 750 Gev" << endl;
  for (int n=0; n<nentries5; n++)			//loops over all the 20,000 entries in file5
  {
    t5->GetEntry(n);
    int count51a = 0;
    int count51b = 0;
    int count51c = 0;
    int count52 = 0;
    int count53 = 0;
    int count54 = 0;
    for ( int n0=0; n0<n_ljet5; n0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt5[n0], ljet_eta5[n0], ljet_phi5[n0], ljet_e5[n0]);

      int btag = 0;
      for (int n1=0; n1<n_jet5; n1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt5[n1], jet_eta5[n1], jet_phi5[n1], jet_e5[n1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count54++;
      }
      for (int n2=0; n2<n_tjet5; n2++)
      {
	TLorentzVector tjet;
	tjet.SetPtEtaPhiE(tjet_pt5[n2], tjet_eta5[n2], tjet_phi5[n2], tjet_e5[n2]);
	float delta_r;
	delta_r = tjet.DeltaR(ljet);
	if (delta_r<0.2 && tjet_btagged5[n2] > 0) btag++;
      }      
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3)
      {
	if (btag == 0) count51a++;
	else if (btag == 1) count51b++;
	else count51c++;
      }
      if (ljet.M()>145.0e3 && btag>0) count52++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3 && btag>0) count53++;
    }
    
    if (count51a > 0 && count52 > 0) nsel51a++;
    else if (count51a > 0 && count53 > 1) nsel52a++;
    else if (count51a > 0 && count53 > 0 && count54 > 0) nsel53a++;
    
    if (count51b > 0 && count52 > 0) nsel51b++;
    else if (count51b > 0 && count53 > 1) nsel52b++;
    else if (count51b > 0 && count53 > 0 && count54 > 0) nsel53b++;
    
    if (count51c > 0 && count52 > 0) nsel51c++;
    else if (count51c > 0 && count53 > 1) nsel52c++;
    else if (count51c > 0 && count53 > 0 && count54 > 0) nsel53c++;
    
  }
  float nsel_51a = nsel51a;
  float nsel_52a = nsel52a;
  float nsel_53a = nsel53a;
  
  float nsel_51b = nsel51b;
  float nsel_52b = nsel52b;
  float nsel_53b = nsel53b;
  
  float nsel_51c = nsel51c;
  float nsel_52c = nsel52c;
  float nsel_53c = nsel53c;
  
  cout << "nsel51: " << nsel_51a << "\t" << nsel_51b << "\t" << nsel_51c << endl;
  cout << "nsel52: " << nsel_52a << "\t" << nsel_52b << "\t" << nsel_52c << endl;
  cout << "nsel53: " << nsel_53a << "\t" << nsel_53b << "\t" << nsel_53c << endl;
  
  smass[3] = 750.0;
  fraction1a[3] = nsel_51a/nentries5;
  fraction2a[3] = nsel_52a/nentries5;
  fraction3a[3] = nsel_53a/nentries5;
  
  fraction1b[3] = nsel_51b/nentries5;
  fraction2b[3] = nsel_52b/nentries5;
  fraction3b[3] = nsel_53b/nentries5;
  
  fraction1c[3] = nsel_51c/nentries5;
  fraction2c[3] = nsel_52c/nentries5;
  fraction3c[3] = nsel_53c/nentries5;
  
  fraction1aerror[3] = sqrt(fraction1a[3] * (1-fraction1a[3])) / sqrt(nentries5);
  fraction2aerror[3] = sqrt(fraction2a[3] * (1-fraction2a[3])) / sqrt(nentries5);
  fraction3aerror[3] = sqrt(fraction3a[3] * (1-fraction3a[3])) / sqrt(nentries5);

  fraction1berror[3] = sqrt(fraction1b[3] * (1-fraction1b[3])) / sqrt(nentries5);
  fraction2berror[3] = sqrt(fraction2b[3] * (1-fraction2b[3])) / sqrt(nentries5);
  fraction3berror[3] = sqrt(fraction3b[3] * (1-fraction3b[3])) / sqrt(nentries5);
  
  fraction1cerror[3] = sqrt(fraction1c[3] * (1-fraction1c[3])) / sqrt(nentries5);
  fraction2cerror[3] = sqrt(fraction2c[3] * (1-fraction2c[3])) / sqrt(nentries5);
  fraction3cerror[3] = sqrt(fraction3c[3] * (1-fraction3c[3])) / sqrt(nentries5);
  
  cout << fraction1a[3] << "\t" << fraction1b[3] << "\t" << fraction1c[3] << endl;
  cout << fraction2a[3] << "\t" << fraction2b[3] << "\t" << fraction1c[3] << endl;
  cout << fraction3a[3] << "\t" << fraction3b[3] << "\t" << fraction3c[3]<< endl;
  
  file5->Close();
  
/////////////////////////////X 2000 S 1000//////////////////////
  
  auto file6 = TFile::Open("bbww_x2000_s1000.root");
  TTree* t6 = (TTree*)file6->Get("CollectionTree");
  
  Int_t           n_jet6;
  Float_t         jet_pt6[13];		//[n_jet]
  Int_t           jet_btagged6[13];	//[n_jet]
  Float_t         jet_eta6[13];		//[n_jet]
  Float_t         jet_phi6[13];		//[n_jet] 
  Float_t         jet_e6[13];		//[n_jet]
  Int_t           n_ljet6;
  Float_t         ljet_pt6[9];   	//[n_ljet]
  Float_t         ljet_eta6[9];   	//[n_ljet]
  Float_t         ljet_phi6[9];   	//[n_ljet]
  Float_t         ljet_e6[9];   	//[n_ljet]
  Int_t           n_tjet6;
  Float_t         tjet_pt6[27];   //[n_tjet]
  Float_t         tjet_eta6[27];   //[n_tjet]
  Float_t         tjet_phi6[27];   //[n_tjet]
  Float_t         tjet_e6[27];   //[n_tjet]
  Int_t           tjet_btagged6[27];   //[n_tjet]  
  
  t6->SetBranchAddress("n_jet", &n_jet6);
  t6->SetBranchAddress("jet_pt", jet_pt6);
  t6->SetBranchAddress("jet_btagged",jet_btagged6);
  t6->SetBranchAddress("jet_eta", jet_eta6);
  t6->SetBranchAddress("jet_phi", jet_phi6);
  t6->SetBranchAddress("jet_e", jet_e6);
  t6->SetBranchAddress("n_ljet", &n_ljet6);
  t6->SetBranchAddress("ljet_pt", ljet_pt6);
  t6->SetBranchAddress("ljet_eta", ljet_eta6);
  t6->SetBranchAddress("ljet_phi", ljet_phi6);
  t6->SetBranchAddress("ljet_e", ljet_e6);
  t6->SetBranchAddress("n_tjet", &n_tjet6);
  t6->SetBranchAddress("tjet_pt", tjet_pt6);
  t6->SetBranchAddress("tjet_eta", tjet_eta6);
  t6->SetBranchAddress("tjet_phi", tjet_phi6);
  t6->SetBranchAddress("tjet_e", tjet_e6);
  t6->SetBranchAddress("tjet_btagged", tjet_btagged6);
  
  t6->SetBranchStatus("*",0);
  
  t6->SetBranchStatus("n_jet",1);
  t6->SetBranchStatus("jet_pt",1);
  t6->SetBranchStatus("jet_btagged",1);
  t6->SetBranchStatus("jet_eta",1);
  t6->SetBranchStatus("jet_phi",1);
  t6->SetBranchStatus("jet_e", 1);
  t6->SetBranchStatus("n_ljet",1);
  t6->SetBranchStatus("ljet_pt", 1);
  t6->SetBranchStatus("ljet_eta", 1);
  t6->SetBranchStatus("ljet_phi", 1);
  t6->SetBranchStatus("ljet_e", 1);
  t6->SetBranchStatus("n_tjet", 1);
  t6->SetBranchStatus("tjet_pt", 1);
  t6->SetBranchStatus("tjet_eta", 1);
  t6->SetBranchStatus("tjet_phi", 1);
  t6->SetBranchStatus("tjet_e", 1);
  t6->SetBranchStatus("tjet_btagged", 1);
  
  int nentries6 = t6->GetEntries();
  int nsel61a = 0;		//counting fully boosted events that are selected
  int nsel61b = 0;
  int nsel61c = 0;
  int nsel62a = 0;		//counting boosted events that are selected
  int nsel62b = 0;
  int nsel62c = 0;
  int nsel63a = 0;		//counting semi-boosted events that are selected
  int nsel63b = 0;
  int nsel63c = 0;
  
  cout << "\nnentries6: " << nentries6 << "\tX: 2000 GeV \tS: 1000 Gev" << endl;
  for (int p=0; p<nentries6; p++)			//loops over all the 20,000 entries in file6
  {
    t6->GetEntry(p);
    int count61a = 0;
    int count61b = 0;
    int count61c = 0;
    int count62 = 0;
    int count63 = 0;
    int count64 = 0;
    
    for ( int p0=0; p0<n_ljet6; p0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt6[p0], ljet_eta6[p0], ljet_phi6[p0], ljet_e6[p0]);

      int btag = 0;
      for (int p1=0; p1<n_jet0; p1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt6[p1], jet_eta6[p1], jet_phi6[p1], jet_e6[p1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count64++;
      }
      for (int p2=0; p2<n_tjet6; p2++)
      {
	TLorentzVector tjet;
	tjet.SetPtEtaPhiE(tjet_pt6[p2], tjet_eta6[p2], tjet_phi6[p2], tjet_e6[p2]);
	float delta_r;
	delta_r = tjet.DeltaR(ljet);
	if (delta_r<0.2 && tjet_btagged0[p2] > 0) btag++;
      }      
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3)
      {
	if (btag == 0) count61a++;
	else if (btag == 1) count61b++;
	else count61c++;
      }
      if (ljet.M()>145.0e3 && btag==0) count62++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3 && btag==0) count63++;
    }
    
    if (count61a > 0 && count62 > 0) nsel61a++;
    else if (count61a > 0 && count63 > 1) nsel62a++;
    else if (count61a > 0 && count63 > 0 && count64 > 0) nsel63a++;
    
    if (count61b > 0 && count62 > 0) nsel61b++;
    else if (count61b > 0 && count63 > 1) nsel62b++;
    else if (count61b > 0 && count63 > 0 && count64 > 0) nsel63b++;
    
    if (count61c > 0 && count62 > 0) nsel61c++;
    else if (count61c > 0 && count63 > 1) nsel62c++;
    else if (count61c > 0 && count63 > 0 && count64 > 0) nsel63c++;
  }
  float nsel_61a = nsel61a;
  float nsel_62a = nsel62a;
  float nsel_63a = nsel63a;
  
  float nsel_61b = nsel61b;
  float nsel_62b = nsel62b;
  float nsel_63b = nsel63b;
  
  float nsel_61c = nsel61c;
  float nsel_62c = nsel62c;
  float nsel_63c = nsel63c;
  
  cout << "nsel61: " << nsel_61a << "\t" << nsel_61b << "\t" << nsel_61c << endl;
  cout << "nsel62: " << nsel_62a << "\t" << nsel_62b << "\t" << nsel_62c << endl;
  cout << "nsel63: " << nsel_63a << "\t" << nsel_63b << "\t" << nsel_63c << endl;
  
  smass[4] = 1000.0;
  fraction1a[4] = nsel_61a/nentries6;
  fraction2a[4] = nsel_62a/nentries6;
  fraction3a[4] = nsel_63a/nentries6;
  
  fraction1b[4] = nsel_61b/nentries6;
  fraction2b[4] = nsel_62b/nentries6;
  fraction3b[4] = nsel_63b/nentries6;
  
  fraction1c[4] = nsel_61c/nentries6;
  fraction2c[4] = nsel_62c/nentries6;
  fraction3c[4] = nsel_63c/nentries6;
  
  fraction1aerror[4] = sqrt(fraction1a[4] * (1-fraction1a[4])) / sqrt(nentries6);
  fraction2aerror[4] = sqrt(fraction2a[4] * (1-fraction2a[4])) / sqrt(nentries6);
  fraction3aerror[4] = sqrt(fraction3a[4] * (1-fraction3a[4])) / sqrt(nentries6);
  
  fraction1berror[4] = sqrt(fraction1b[4] * (1-fraction1b[4])) / sqrt(nentries6);
  fraction2berror[4] = sqrt(fraction2b[4] * (1-fraction2b[4])) / sqrt(nentries6);
  fraction3berror[4] = sqrt(fraction3b[4] * (1-fraction3b[4])) / sqrt(nentries6);
  
  fraction1cerror[4] = sqrt(fraction1c[4] * (1-fraction1c[4])) / sqrt(nentries6);
  fraction2cerror[4] = sqrt(fraction2c[4] * (1-fraction2c[4])) / sqrt(nentries6);
  fraction3cerror[4] = sqrt(fraction3c[4] * (1-fraction3c[4])) / sqrt(nentries6);
  
  cout << fraction1a[4] << "\t" << fraction1b[4] << "\t" << fraction1c[4] << endl;
  cout << fraction2a[4] << "\t" << fraction2b[4] << "\t" << fraction1c[4] << endl;
  cout << fraction3a[4] << "\t" << fraction3b[4] << "\t" << fraction3c[4]<< endl;
  
  file6->Close();
 
/////////////////////////////////////////////////////////////////////////////////
  
  TFile f("b_selectionEfficiency_3btags.root", "recreate");

  ca = new TCanvas("ca", "Event Selection Efficiency for different S mass points", 1200, 900);
  ca -> Divide(1,1);
  
  float ex[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  ca -> cd(1);
  gr1 = new TGraphErrors(6, smass, fraction1a, ex, fraction1aerror);
  gr2 = new TGraphErrors(6, smass, fraction2a, ex, fraction2aerror);
  gr3 = new TGraphErrors(6, smass, fraction3a, ex, fraction3aerror);
  
  gr1->SetFillColor(kWhite);
  gr1->SetLineColor(3);
  gr1->SetLineWidth(2);
  gr1->SetMarkerStyle(2);
  gr1->SetMarkerColor(3);
  gr1->SetTitle("Event Selection Efficiency at different S mass points for M_{X} = 2 TeV (with b-tagging)");
  gr1->GetXaxis()->SetTitle("Mass of the S Particle, M_{S} [GeV]");
  gr1->GetYaxis()->SetRangeUser(0.0, 0.3); 
  gr1->GetYaxis()->SetTitle("Fraction of Events");
  
  gr2->SetFillColor(kWhite);
  gr2->SetLineColor(2);
  gr2->SetLineWidth(2);
  gr2->SetMarkerStyle(2);
  gr2->SetMarkerColor(2);
  
  gr3->SetFillColor(kWhite);
  gr3->SetLineColor(4);
  gr3->SetLineWidth(2);
  gr3->SetMarkerStyle(2);
  gr3->SetMarkerColor(4);
  
  gr1->Draw();
  gr2->Draw("same");
  gr3->Draw("same");
  
  auto legend = new TLegend(0.6,0.6,0.9,0.8);	// (x1, y1, x2, y2)
  legend->AddEntry(gr1,"Fully Boosted: 1 H + WW jet");
  legend->AddEntry(gr2,"Boosted: 1 H + 2 W");
  legend->AddEntry(gr3, "Semi-Boosted: 1 H + 1 W + 1 W_{j j}");
  legend->Draw();
  
  f.Close();
  
}
