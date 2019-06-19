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
//#include "TLorentzVector"
#include "TMath.h"
#include <math.h>
#include <cmath>
#include "TCanvas.h"


//This code graphs the invariant mass of regular b-jets // i am using the same ttbar monte carlo simulated data as Merrick Lavinsky is.
void ttbar_4() {
  
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int w = 1; 					//weight of the bins
  int nbin = 1000;				//bin size
  
  auto file = TFile::Open("ttbar.root");		//Creates a pointer to the main data file
  TTree* t = (TTree*)file->Get("nominal");		//Creates a pointer to the tree
  auto h_M_inv = new TH1F("M_inv", "Invariant Mass; Mass[GeV]; Events", nbin, 0.0, 6.0e5);
						//Creates a pointer to the histogram to be created.
  TFile f("ttbar_invariant_mass_4.root", "RECREATE");
  
  vector<float> *jet_pt		=0;
  vector<char>	*jet_btagged	=0;
  vector<float>	*jet_eta	=0;
  vector<float>	*jet_phi	=0;
  
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("jet_isbtagged_MV2c10_85", &jet_btagged);
  t->SetBranchAddress("jet_eta", &jet_eta);
  t->SetBranchAddress("jet_phi", &jet_phi);
  
  int nentries = t->GetEntries();
  
  float x;
  float y;
  float z;
  float m;
  float M_inv;

  cout << "Number of entries: " << nentries << endl;
  
  for (int i=0; i<nentries; i++)
  {
    t->GetEntry(i);
    int n_bjets = 0;				//this int variable is like a step for counting number of b-jets within each entry.
    int n_jets = jet_pt->size();
    
    vector<float> pt_bjet_all {0, 0};		//stores temporary pt bjets with more than 1 btags in each loop
    vector<float> phi_bjet_all {0, 0};		//stores temporary phi angles with more than 1 btags in each loop
    vector<float> eta_bjet_all {0, 0};		//stores temporary eta angles with more than 1 btags in each loop
    
    for (int j=0; j<n_jets; j++)
    {
      if ((int)(*jet_btagged)[j]==1)
      {
	if ( n_bjets < 2)
	{
	  pt_bjet_all[n_bjets] = (float)(*jet_pt)[j];
	  phi_bjet_all[n_bjets] = (float)(*jet_phi)[j];
	  eta_bjet_all[n_bjets] = (float)(*jet_eta)[j];
	  
	}
	n_bjets++;
      }
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
	  
      
      M_inv = (sqrt ( pow(m , 2) - pow(x , 2) - pow(y , 2) - pow(z , 2) ));      
      h_M_inv->Fill(M_inv, w);

    }
    
  }
  
  file->Close();
  h_M_inv->Write(); 
  
  f.Write();
  f.Close();  
  
}
