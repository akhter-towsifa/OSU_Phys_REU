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

 
//This code is to find the selection efficiency for all the different mass point and all the boosted types (without taking regular jet b-tagging into account).

void a_selectionEfficiency(){
  
  float smass[6], fraction1[6], fraction2[6], fraction3[6];	//smass list collects all the 6 S mass; fraction1,2,3 lists collect the nsel/nentries ratio at each of these mass points 
  float fraction1error[6], fraction2error[6], fraction3error[6];
  float totalEfficiency[6], totalEfficiencyError[6];
  
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
  
  int nentries0 = t0->GetEntries();
  int nsel01 = 0;		//counting fully boosted events that are selected
  int nsel02 = 0;		//counting boosted events that are selected
  int nsel03 = 0;		//counting semi-boosted events that are selected
  
  cout << "\nnentries0: " << nentries0 << "\tX: 1000 GeV \tS: 170 Gev" << endl;
  for (int i=0; i<nentries0; i++)			//loops over all the 20,000 entries in file0
  {
    t0->GetEntry(i);
    int count01 = 0;		//counts the number of higgs within each Ljet
    int count02 = 0;		//counts the number of WW Ljets
    int count03 = 0;		//counts the number of W Ljets
    int count04 = 0;		//counts the number of regular W jets
    
    for (int i0=0; i0<n_ljet0; i0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt0[i0], ljet_eta0[i0], ljet_phi0[i0], ljet_e0[i0]);
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) count01++;
      if (ljet.M()>145.0e3) count02++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3) count03++;
      
      int countb0 = 0;
      
      for (int i1=0; i1<n_jet0; i1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt0[i1], jet_eta0[i1], jet_phi0[i1], jet_e0[i1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count04++;
      }
    }
    if (count01 > 0 && count02 > 0) nsel01++;
    else if (count01 > 0 && count03 > 1) nsel02++;
    else if (count01 > 0 && count03 > 0 && count04 > 0) nsel03++;

  }
  cout << "nsel01: " << nsel01 << endl;
  cout << "nsel02: " << nsel02 << endl;
  cout << "nsel03: " << nsel03 << endl;
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
  
  int nentries1 = t1->GetEntries();
  int nsel11 = 0;		//counting fully boosted events that are selected
  int nsel12 = 0;		//counting boosted events that are selected
  int nsel13 = 0;		//counting semi-boosted events that are selected
  
  cout << "\nnentries1: " << nentries1 << "\tX: 2000 GeV \tS: 1500 Gev" << endl;
  for (int j=0; j<nentries1; j++)			//loops over all the 20,000 entries in file1
  {
    t1->GetEntry(j);
    int count11 = 0;
    int count12 = 0;
    int count13 = 0;
    int count14 = 0;
    
    for ( int j0=0; j0<n_ljet1; j0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt1[j0], ljet_eta1[j0], ljet_phi1[j0], ljet_e1[j0]);
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) count11++;
      if (ljet.M()>145.0e3) count12++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3) count13++;
      
      for (int j1=0; j1<n_jet1; j1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt1[j1], jet_eta1[j1], jet_phi1[j1], jet_e1[j1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count14++;
      }

    }
    if (count11 > 0 && count12 > 0) nsel11++;
    else if (count11 > 0 && count13 > 1) nsel12++;
    else if (count11 > 0 && count13 > 0 && count14 > 0) nsel13++;

  }
  float nsel_11 = nsel11;
  float nsel_12 = nsel12;
  float nsel_13 = nsel13;
  
  cout << "nsel11: " << nsel_11 << endl;
  cout << "nsel12: " << nsel_12 << endl;
  cout << "nsel13: " << nsel_13 << endl;
  
  smass[5] = 1500.0;
  fraction1[5] = nsel_11/nentries1;
  fraction2[5] = nsel_12/nentries1;
  fraction3[5] = nsel_13/nentries1;
  
  fraction1error[5] = sqrt(fraction1[5] * (1-fraction1[5])) / sqrt(nentries1);
  fraction2error[5] = sqrt(fraction2[5] * (1-fraction2[5])) / sqrt(nentries1);
  fraction3error[5] = sqrt(fraction3[5] * (1-fraction3[5])) / sqrt(nentries1);
  
  totalEfficiency[5] = fraction1[5] + fraction2[5] + fraction3[5];
  totalEfficiencyError[5] = sqrt(totalEfficiency[5] * (1-totalEfficiency[5])) / sqrt(nentries1);
  
  cout << fraction1[5] << endl;
  cout << fraction2[5] << endl;
  cout << fraction3[5] << endl;
  
  cout << "total efficiency: " << totalEfficiency[5] << " +/- " << totalEfficiencyError[5] << endl;
  
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
  
  int nentries2 = t2->GetEntries();
  int nsel21 = 0;		//counting fully boosted events that are selected
  int nsel22 = 0;		//counting boosted events that are selected
  int nsel23 = 0;		//counting semi-boosted events that are selected
  
  cout << "\nnentries2: " << nentries2 << "\tX: 2000 GeV \tS: 170 Gev" << endl;
  for (int k=0; k<nentries2; k++)			//loops over all the 20,000 entries in file2
  {
    t2->GetEntry(k);
    int count21 = 0;
    int count22 = 0;
    int count23 = 0;
    int count24 = 0;
    for ( int k0=0; k0<n_ljet2; k0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt2[k0], ljet_eta2[k0], ljet_phi2[k0], ljet_e2[k0]);
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) count21++;
      if (ljet.M()>145.0e3) count22++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3) count23++;
      
      for (int k1=0; k1<n_jet2; k1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt2[k1], jet_eta2[k1], jet_phi2[k1], jet_e2[k1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count24++;
      }
    }
    if (count21 > 0 && count22 > 0) nsel21++;
    else if (count21 > 0 && count23 > 1) nsel22++;
    else if (count21 > 0 && count23 > 0 && count24 > 0) nsel23++;
  }
  float nsel_21 = nsel21;
  float nsel_22 = nsel22;
  float nsel_23 = nsel23;
  
  cout << "nsel21: " << nsel_21 << endl;
  cout << "nsel22: " << nsel_22 << endl;
  cout << "nsel23: " << nsel_23 << endl;
  
  smass[0] = 170.0;
  fraction1[0] = nsel_21/nentries2;
  fraction2[0] = nsel_22/nentries2;
  fraction3[0] = nsel_23/nentries2;
  
  fraction1error[0] = sqrt(fraction1[0] * (1-fraction1[0])) / sqrt(nentries2);
  fraction2error[0] = sqrt(fraction2[0] * (1-fraction2[0])) / sqrt(nentries2);
  fraction3error[0] = sqrt(fraction3[0] * (1-fraction3[0])) / sqrt(nentries2);
  
  totalEfficiency[0] = fraction1[0] + fraction2[0] + fraction3[0];
  totalEfficiencyError[0] = sqrt(totalEfficiency[0] * (1-totalEfficiency[0])) / sqrt(nentries2);
  
  cout << fraction1[0] << endl;
  cout << fraction2[0] << endl;
  cout << fraction3[0] << endl;
  
  cout << "total efficiency: " << totalEfficiency[0] << " +/- " << totalEfficiencyError[0] << endl;
  
  file2->Close();
/////////////////////////////////////////////////
  
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
  
  int nentries3 = t3->GetEntries();
  int nsel31 = 0;		//counting fully boosted events that are selected
  int nsel32 = 0;		//counting boosted events that are selected
  int nsel33 = 0;		//counting semi-boosted events that are selected
  
  cout << "\nnentries3: " << nentries3 << "\tX: 2000 GeV \tS: 240 Gev" << endl;
  for (int l=0; l<nentries3; l++)			//loops over all the 20,000 entries in file3
  {
    t3->GetEntry(l);
    int count31 = 0;
    int count32 = 0;
    int count33 = 0;
    int count34 = 0;
    
    for ( int l0=0; l0<n_ljet3; l0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt3[l0], ljet_eta3[l0], ljet_phi3[l0], ljet_e3[l0]);
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) count31++;
      if (ljet.M()>145.0e3) count32++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3) count33++;
      
      for (int l1=0; l1<n_jet3; l1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt3[l1], jet_eta3[l1], jet_phi3[l1], jet_e3[l1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count34++;
      }
    }
    if (count31 > 0 && count32 > 0) nsel31++;
    else if (count31 > 0 && count33 > 1) nsel32++;
    else if (count31 > 0 && count33 > 0 && count34 > 0) nsel33++;
  }
  float nsel_31 = nsel31;
  float nsel_32 = nsel32;
  float nsel_33 = nsel33;
  
  cout << "nsel31: " << nsel_31 << endl;
  cout << "nsel32: " << nsel_32 << endl;
  cout << "nsel33: " << nsel_33 << endl;
  
  smass[1] = 240.0;
  fraction1[1] = nsel_31/nentries3;
  fraction2[1] = nsel_32/nentries3;
  fraction3[1] = nsel_33/nentries3;
  
  fraction1error[1] = sqrt(fraction1[1] * (1-fraction1[1])) / sqrt(nentries3);
  fraction2error[1] = sqrt(fraction2[1] * (1-fraction2[1])) / sqrt(nentries3);
  fraction3error[1] = sqrt(fraction3[1] * (1-fraction3[1])) / sqrt(nentries3);
  
  totalEfficiency[1] = fraction1[1] + fraction2[1] + fraction3[1];
  totalEfficiencyError[1] = sqrt(totalEfficiency[1] * (1-totalEfficiency[1])) / sqrt(nentries3);
  
  cout << fraction1[1] << endl;
  cout << fraction2[1] << endl;
  cout << fraction3[1] << endl;
  
  cout << "total efficiency: " << totalEfficiency[1] << " +/- " << totalEfficiencyError[1] << endl;
  
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
  
  int nentries4 = t4->GetEntries();
  int nsel41 = 0;		//counting fully boosted events that are selected
  int nsel42 = 0;		//counting boosted events that are selected
  int nsel43 = 0;		//counting semi-boosted events that are selected
  
  cout << "\nnentries4: " << nentries4 << "\tX: 2000 GeV \tS: 400 Gev" << endl;
  for (int m=0; m<nentries4; m++)			//loops over all the 20,000 entries in file4
  {
    t4->GetEntry(m);
    int count41 = 0;
    int count42 = 0;
    int count43 = 0;
    int count44 = 0;
    for ( int m0=0; m0<n_ljet4; m0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt4[m0], ljet_eta4[m0], ljet_phi4[m0], ljet_e4[m0]);
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) count41++;
      if (ljet.M()>145.0e3) count42++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3) count43++;
      
      for (int m1=0; m1<n_jet4; m1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt4[m1], jet_eta4[m1], jet_phi4[m1], jet_e4[m1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count44++;
      }
    }
    
    if (count41 > 0 && count42 > 0) nsel41++;
    else if (count41 > 0 && count43 > 1) nsel42++;
    else if (count41 > 0 && count43 > 0 && count44 > 0) nsel43++;
  }
  float nsel_41 = nsel41;
  float nsel_42 = nsel42;
  float nsel_43 = nsel43;
  
  cout << "nsel41: " << nsel_41 << endl;
  cout << "nsel42: " << nsel_42 << endl;
  cout << "nsel43: " << nsel_43 << endl;
  
  smass[2] = 400.0;
  fraction1[2] = nsel_41/nentries4;
  fraction2[2] = nsel_42/nentries4;
  fraction3[2] = nsel_43/nentries4;
  
  fraction1error[2] = sqrt(fraction1[2] * (1-fraction1[2])) / sqrt(nentries4);
  fraction2error[2] = sqrt(fraction2[2] * (1-fraction2[2])) / sqrt(nentries4);
  fraction3error[2] = sqrt(fraction3[2] * (1-fraction3[2])) / sqrt(nentries4);
  
  totalEfficiency[2] = fraction1[2] + fraction2[2] + fraction3[2];
  totalEfficiencyError[2] = sqrt(totalEfficiency[2] * (1-totalEfficiency[2])) / sqrt(nentries4);
  
  cout << fraction1[2] << endl;
  cout << fraction2[2] << endl;
  cout << fraction3[2] << endl;
  
  cout << "total efficiency: " << totalEfficiency[2] << " +/- " << totalEfficiencyError[2] << endl;
  
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
  
  int nentries5 = t5->GetEntries();
  int nsel51 = 0;		//counting fully boosted events that are selected
  int nsel52 = 0;		//counting boosted events that are selected
  int nsel53 = 0;		//counting semi-boosted events that are selected
  cout << "\nnentries5: " << nentries5 << "\tX: 2000 GeV \tS: 750 Gev" << endl;
  for (int n=0; n<nentries5; n++)			//loops over all the 20,000 entries in file5
  {
    t5->GetEntry(n);
    int count51 = 0;
    int count52 = 0;
    int count53 = 0;
    int count54 = 0;
    for ( int n0=0; n0<n_ljet5; n0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt5[n0], ljet_eta5[n0], ljet_phi5[n0], ljet_e5[n0]);
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) count51++;
      if (ljet.M()>145.0e3) count52++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3) count53++;
      
      for (int n1=0; n1<n_jet5; n1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt5[n1], jet_eta5[n1], jet_phi5[n1], jet_e5[n1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count54++;
      }
    }
    
    if (count51 > 0 && count52 > 0) nsel51++;
    else if (count51 > 0 && count53 > 1) nsel52++;
    else if (count51 > 0 && count53 > 0 && count54 > 0) nsel53++;
    
  }
  float nsel_51 = nsel51;
  float nsel_52 = nsel52;
  float nsel_53 = nsel53;
  
  cout << "nsel51: " << nsel_51 << endl;
  cout << "nsel52: " << nsel_52 << endl;
  cout << "nsel53: " << nsel_53 << endl;
  
  smass[3] = 750.0;
  fraction1[3] = nsel_51/nentries5;
  fraction2[3] = nsel_52/nentries5;
  fraction3[3] = nsel_53/nentries5;
  
  fraction1error[3] = sqrt(fraction1[3] * (1-fraction1[3])) / sqrt(nentries5);
  fraction2error[3] = sqrt(fraction2[3] * (1-fraction2[3])) / sqrt(nentries5);
  fraction3error[3] = sqrt(fraction3[3] * (1-fraction3[3])) / sqrt(nentries5);
  
  totalEfficiency[3] = fraction1[3] + fraction2[3] + fraction3[3];
  totalEfficiencyError[3] = sqrt(totalEfficiency[3] * (1-totalEfficiency[3])) / sqrt(nentries5);
  
  cout << fraction1[3] << endl;
  cout << fraction2[3] << endl;
  cout << fraction3[3] << endl;
  
  cout << "total efficiency: " << totalEfficiency[3] << " +/- " << totalEfficiencyError[3] << endl;
  
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
  
  int nentries6 = t6->GetEntries();
  int nsel61 = 0;		//counting fully boosted events that are selected
  int nsel62 = 0;		//counting boosted events that are selected
  int nsel63 = 0;		//counting semi-boosted events that are selected
  
  cout << "\nnentries6: " << nentries6 << "\tX: 2000 GeV \tS: 1000 Gev" << endl;
  for (int p=0; p<nentries6; p++)			//loops over all the 20,000 entries in file6
  {
    t6->GetEntry(p);
    int count61 = 0;
    int count62 = 0;
    int count63 = 0;
    int count64 = 0;
    
    for ( int p0=0; p0<n_ljet6; p0++)			//loops over all the large jets in each entry
    {
      TLorentzVector ljet;				//Defines a Lorentz vector "ljet" that stores the large jets in each iteration j over all large jets
      ljet.SetPtEtaPhiE(ljet_pt6[p0], ljet_eta6[p0], ljet_phi6[p0], ljet_e6[p0]);
      if (ljet.M()>105.0e3 && ljet.M()<145.0e3) count61++;
      if (ljet.M()>145.0e3) count62++;
      if (ljet.M()>60.0e3 && ljet.M()<100.0e3) count63++;
      
      for (int p1=0; p1<n_jet0; p1++)			//loops over all the regular jets within each large jet
      {
	TLorentzVector jet;
	jet.SetPtEtaPhiE(jet_pt6[p1], jet_eta6[p1], jet_phi6[p1], jet_e6[p1]);
	float delta_r;
	delta_r = jet.DeltaR(ljet);
	if (delta_r > 0.2 && jet.M()>60.0e3 && jet.M()<100.0e3) count64++;
      }
    }
    
    if (count61 > 0 && count62 > 0) nsel61++;
    else if (count61 > 0 && count63 > 1) nsel62++;
    else if (count61 > 0 && count63 > 0 && count64 > 0) nsel63++;
    
  }
  float nsel_61 = nsel61;
  float nsel_62 = nsel62;
  float nsel_63 = nsel63;
  
  cout << "nsel61: " << nsel_61 << endl;
  cout << "nsel62: " << nsel_62 << endl;
  cout << "nsel63: " << nsel_63 << endl;
  
  smass[4] = 1000.0;
  fraction1[4] = nsel_61/nentries6;
  fraction2[4] = nsel_62/nentries6;
  fraction3[4] = nsel_63/nentries6;
  
  fraction1error[4] = sqrt(fraction1[4] * (1-fraction1[4])) / sqrt(nentries6);
  fraction2error[4] = sqrt(fraction2[4] * (1-fraction2[4])) / sqrt(nentries6);
  fraction3error[4] = sqrt(fraction3[4] * (1-fraction3[4])) / sqrt(nentries6);
  
  totalEfficiency[4] = fraction1[4] + fraction2[4] + fraction3[4];
  totalEfficiencyError[4] = sqrt(totalEfficiency[4] * (1-totalEfficiency[4])) / sqrt(nentries6);
  
  cout << fraction1[4] << endl;
  cout << fraction2[4] << endl;
  cout << fraction3[4] << endl;
  
  cout << "total efficiency: " << totalEfficiency[4] << " +/- " << totalEfficiencyError[4] << endl;
  
  file6->Close();
 
/////////////////////////////////////////////////////////////////////////////////
  
  TFile f("a_selectionEfficiency.root", "recreate");

  ca = new TCanvas("ca", "Event Selection Efficiency for different S mass points", 1200, 900);
  ca -> SetFillColor(kWhite);
  float ex[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  gr1 = new TGraphErrors(6, smass, fraction1, ex, fraction1error);
  gr2 = new TGraphErrors(6, smass, fraction2, ex, fraction2error);
  gr3 = new TGraphErrors(6, smass, fraction3, ex, fraction3error);
  gr4 = new TGraphErrors(6, smass, totalEfficiency, ex, totalEfficiencyError);
  
  gr1->SetFillColor(kWhite);
  gr1->SetLineColor(3);
  gr1->SetLineWidth(3);
  gr1->SetMarkerStyle(2);
  gr1->SetMarkerColor(3);
  gr1->SetTitle("Event Selection Efficiency at different S mass points for M_{X} = 2 TeV (without b-tagging)");
  gr1->GetXaxis()->SetTitle("M_{S}    [GeV]");
  gr1->GetYaxis()->SetRangeUser(0.0, 0.6); 
  gr1->GetYaxis()->SetTitle("Fraction of Events");
  
  gr2->SetFillColor(kWhite);
  gr2->SetLineColor(2);
  gr2->SetLineWidth(3);
  gr2->SetMarkerStyle(2);
  gr2->SetMarkerColor(2);
  
  gr3->SetFillColor(kWhite);
  gr3->SetLineColor(4);
  gr3->SetLineWidth(3);
  gr3->SetMarkerStyle(2);
  gr3->SetMarkerColor(4);
  
  gr4->SetFillColor(kWhite);
  gr4->SetLineColor(1);
  gr4->SetLineWidth(3);
  gr4->SetLineStyle(10);
  gr4->SetMarkerStyle(1);
  gr4->SetMarkerColor(1);
  
  gr1->Draw();
  gr2->Draw("same");
  gr3->Draw("same");
  gr4->Draw("same");
  
  auto legend = new TLegend(0.6,0.6,0.9,0.8);	// (x1, y1, x2, y2)
  legend->SetFillColor(kWhite);
  legend->AddEntry(gr4, "Total Efficiency");
  legend->AddEntry(gr1,"Fully Boosted: 1 H + WW jet");
  legend->AddEntry(gr2,"Boosted: 1 H + 2 W");
  legend->AddEntry(gr3, "Semi-Boosted: 1 H + 1 W + 1 W_{j j}");
  legend->Draw();
  
  f.Close();
}