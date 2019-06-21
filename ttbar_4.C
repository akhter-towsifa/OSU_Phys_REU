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
  gStyle->SetOptStat("nemi");
  int w = 1; 					//weight of the bins
  int nbin = 1000;				//bin size
  
  auto h_M_inv = new TH1F("M_inv", "Invariant Mass of Regular and Large Jets; Mass [GeV]; Events", nbin, 0.0, 2.0e6);
						//Creates a pointer to the histogram to be created.
  auto file = TFile::Open("ttbar.root");	//Creates a pointer to the main data file
  TTree* t = (TTree*)file->Get("nominal");	//Creates a pointer to the tree
  
  TFile f("ttbar_invariant_mass_4.root", "RECREATE");
  
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
	if ( n_bjets < 1)
	{
	  pt_bjet_all[n_bjets] = (float)(*jet_pt)[j];
	  phi_bjet_all[n_bjets] = (float)(*jet_phi)[j];
	  eta_bjet_all[n_bjets] = (float)(*jet_eta)[j];
	}
	n_bjets++;
      }
    }
    
    n_bjets = 1;
    
    for (int k=0; k<n_ljets; k++)		//Large jets
    {
      if ((float)(*ljet_pt)[k] >= 420000 && abs((float)(*ljet_eta)[k]) <= 2.0)
      {
// 	cout << k << "\t\t" << (float)(*ljet_pt)[k] << "\t\t"  << abs((float)(*ljet_eta)[k]) << endl;
	pt_bjet_all[n_bjets] = (float)(*ljet_pt)[k];
	phi_bjet_all[n_bjets] = (float)(*ljet_phi)[k];
	eta_bjet_all[n_bjets] = (float)(*ljet_eta)[k];
      }
      n_bjets++;
    }

    if (n_bjets >= 2)
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
// 	cout << x << "\t\t" << y << "\t\t" << z << "\t\t" << m << endl;
      
      M_inv = (sqrt ( pow(m , 2) - pow(x , 2) - pow(y , 2) - pow(z , 2) ));      
      h_M_inv->Fill(M_inv, w);
//       cout << i << "\t\tMass calculated: " << M_inv << endl;
    }
    
  }
  
  file->Close();
  h_M_inv->Write(); 
   
}
