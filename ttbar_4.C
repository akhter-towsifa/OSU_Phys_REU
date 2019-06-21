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
//#include "TLorentzVector.h"
#include "TMath.h"
#include <math.h>
#include <cmath>
#include "TCanvas.h"


//This code graphs the invariant mass of regular and large b-jets // i am using a ttbar monte carlo simulated data
void ttbar_4() {
  
  
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
  vector<float> *ljet_pt	=0;
  vector<float> *ljet_eta	=0;
  vector<float> *ljet_phi	=0;
  vector<float> *ljet_m		=0;
  
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("jet_isbtagged_MV2c10_85", &jet_btagged);
  t->SetBranchAddress("jet_eta", &jet_eta);
  t->SetBranchAddress("jet_phi", &jet_phi);
  t->SetBranchAddress("ljet_pt", &ljet_pt);
  t->SetBranchAddress("ljet_eta", &ljet_eta);
  t->SetBranchAddress("ljet_phi", &ljet_phi);
  t->SetBranchAddress("ljet_m", &ljet_m);
  
  int nentries = t->GetEntries();
  
  float x;			// the temp pt_x components
  float y;			// the temp pt_y components
  float z;			// the temp pt_z components
  float m;			// the temp total pt components
  float M_inv;			// the temp calculated M_inv components

  cout << "Number of entries: " << nentries << endl;
  
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int n_bjets = 0;				//this int variable is like a step for counting number of b-jets within each entry.
    int n_jets = jet_pt->size();
    int n_ljets = ljet_pt->size();
    
    vector<float> pt_bjet_all {0, 0};		//stores temporary pt bjets with more than 1 btag in each loop
    vector<float> phi_bjet_all {0, 0};		//stores temporary phi angles with more than 1 btag in each loop
    vector<float> eta_bjet_all {0, 0};		//stores temporary eta angles with more than 1 btag in each loop
    
    for (int j=0; j<n_jets; j++)		//regular jets
    {
      if ((int)(*jet_btagged)[j]==1)
      {
	if ( n_bjets < 1)			//This loop collects the pt, phi, and eta of the highest jet within the regular jets
	{
	  pt_bjet_all[n_bjets] = (float)(*jet_pt)[j];
	  phi_bjet_all[n_bjets] = (float)(*jet_phi)[j];
	  eta_bjet_all[n_bjets] = (float)(*jet_eta)[j];
	
	  n_bjets++;
	}
      }
    }
//     cout << "n_bjets after the first for loop\t\t" << n_bjets << endl;

    float ljet_m_compare = 0;
    int n_t_candidates=0;			//a step variable that counts the number of large jets or t quarks
    for (int k=0; k<n_ljets; k++)		//Large jets
    {
      if ((float)(*ljet_pt)[k] >= 420000 && abs((float)(*ljet_eta)[k]) <= 2.0 && n_bjets==1)		//This loop collects the pt, phi, and eta of the highest jet within the large jets
      {
	n_t_candidates++;
	float ljet_m_current = (float)(*ljet_m)[k];
	ljet_m_compare = max(ljet_m_compare, ljet_m_current);
	
	if (ljet_m_compare == ljet_m_current)
	{
// 	cout << k << "\t\t" << (float)(*ljet_pt)[k] << "\t\t"  << abs((float)(*ljet_eta)[k]) << endl;
	  pt_bjet_all[n_bjets] = (float)(*ljet_pt)[k];
	  phi_bjet_all[n_bjets] = (float)(*ljet_phi)[k];
	  eta_bjet_all[n_bjets] = (float)(*ljet_eta)[k];
	}
      }
    }
    
  
//      cout << "n_bjets after the second for loop: " << n_bjets << "\t\tt-cndidates: " << n_t_candidates << endl;
    
    if (n_bjets == 1 && n_t_candidates > 0)			//This loop calculates the invariant mass where [0] refers to the regular jet and [1] refers to the large jet.
    {
      x = pt_bjet_all[0] * cos(phi_bjet_all[0]) + pt_bjet_all[1] * cos(phi_bjet_all[1]);
      y = pt_bjet_all[0] * sin (phi_bjet_all[0]) + pt_bjet_all[1] * sin (phi_bjet_all[1]);
      z = pt_bjet_all[0] * sinh(eta_bjet_all[0]) + pt_bjet_all[1] * sinh(eta_bjet_all[1]);
      m = sqrt( pow(pt_bjet_all[0] * cos(phi_bjet_all[0]) , 2)
		+pow(pt_bjet_all[0] * sin(phi_bjet_all[0]) , 2)
		+pow(pt_bjet_all[0] * sinh(eta_bjet_all[0]) , 2)
	   )
	  +sqrt( pow(pt_bjet_all[1] * cos(phi_bjet_all[1]) , 2)
		+pow(pt_bjet_all[1] * sin(phi_bjet_all[1]) , 2)
		+pow(pt_bjet_all[1] * sinh(eta_bjet_all[1]) , 2)
	   );
      M_inv = (sqrt ( pow(m , 2) - pow(x , 2) - pow(y , 2) - pow(z , 2) ));      

      if (M_inv <= 25000) 
      {
	cout << i << "\t\tMass calculated: " << M_inv << ";\tm: " << m << ";\tx: " << x << ";\ty: " << y << ";\tz: " << z <<"\n" << endl;
	h_M_inv->Fill(M_inv, w);
      }
    }
  }  
   
  TCanvas *c = new TCanvas("canvas", "canvas", 1000, 1000);
  h_M_inv->Draw();
    
  
  
  file->Close();
  TFile f("ttbar_invariant_mass_4.root", "RECREATE");
  h_M_inv->Write(); 
  f.Close();
   
}

//The lines below are to look at the values of each of the pt components............ not part of the code

// 	float px0 = pt_bjet_all[0] * cos(phi_bjet_all[0]);
// 	float px1 = pt_bjet_all[1] * cos(phi_bjet_all[1]);
// 	float py0 = pt_bjet_all[0] * sin(phi_bjet_all[0]);
// 	float py1 = pt_bjet_all[1] * sin(phi_bjet_all[1]);
// 	float pz0 = pt_bjet_all[0] * sinh(eta_bjet_all[0]);
// 	float pz1 = pt_bjet_all[1] * sinh(eta_bjet_all[1]);
// 	float temp_m = 2 * sqrt( pow(pt_bjet_all[0] * cos(phi_bjet_all[0]) , 2)
// 		+pow(pt_bjet_all[0] * sin(phi_bjet_all[0]) , 2)
// 		+pow(pt_bjet_all[0] * sinh(eta_bjet_all[0]) , 2)
// 	   ) * sqrt( pow(pt_bjet_all[1] * cos(phi_bjet_all[1]) , 2)
// 		+pow(pt_bjet_all[1] * sin(phi_bjet_all[1]) , 2)
// 		+pow(pt_bjet_all[1] * sinh(eta_bjet_all[1]) , 2)
// 	   );
// 	cout << "px0\t\t" << px0 << endl;
// 	cout << "px1\t\t" << px1 << endl;
// 	cout << "py0\t\t" << py0 << endl;
// 	cout << "py1\t\t" << py1 << endl;
// 	cout << "pz0\t\t" << pz0 << endl;
// 	cout << "pz1\t\t" << pz1 << endl;
// 	cout << "ptotal\t\t" << temp_m << endl;
/*

    #if defined(R__FAST_MATH)
      inline Bool_t TMath::IsNaN(Float_t M_inv);	//{return isnan(M_inv);}
      cout << i << endl;
    #endif*/