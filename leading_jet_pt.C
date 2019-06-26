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


/*This code stores all the jets that have a b particle in it. Then it selects
two of the most energetic leading b-jets [pt_bjet].*/

void leading_jet_pt(){
  
  gStyle->SetOptStat("nemr");
  int nbins = 1000;
  int w_2 = 1;			//weight

  auto h_jet_pt 		= new TH1F("jet_pt","Jet_P_{T}; Mass [MeV]; Counts",		nbins,0.0, 1.0e6);
  auto h_jet_pt_btag 		= new TH1I("jet_btagged","Jet p_{T} b-tagged;Mass[MeV]",	nbins,0.0, 1.0e6);
  auto h_leading_bjet_1 	= new TH1F("h_leading_bjet_1", "First Leading b-Jet; Mass[MeV]; Counts", nbins, 0.0, 2.0e6);
  auto h_leading_bjet_2 	= new TH1F("h_leading_bjet_2", "Leading and Subleading b-Jets; Mass[MeV]; Counts", nbins, 0.0, 2.0e6);
  auto h_M_inv 			= new TH1F("M_inv", "Invariant Mass of Higgs Boson from the b-Jets; Mass[GeV]; Counts",	nbins,0.0, 200.0);
  auto h_w_pt_btag		= new TH1F("W_particle", "W P_{T}; Mass[MeV]; Counts", 		nbins, 0.0, 1.0e6);
  auto h_w_inv			= new TH1F("W1_inv", "Invariant mass of W; mass[GeV]; Counts", 	nbins, 0.0, 1000.0);
  auto h_w2_inv			= new TH1F("W2_inv", "Invariant mass of 2nd W; mass[GeV]; Counts", 	nbins, 0.0, 1800.0);
  
  auto file = TFile::Open("bbww_x1000_s170.root");

  TTree* ctree = (TTree*)file->Get("CollectionTree");
  
  Int_t           n_jet;
  Float_t         jet_pt[14];   	//[n_jet]
  Int_t           jet_btagged[14];   	//[n_jet]
  Float_t         jet_eta[14];   	//[n_jet]
  Float_t         jet_phi[14];  	//[n_jet]
  Float_t         jet_e[14];  		//[n_jet]
  Float_t		w;
  
  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt", jet_pt);
  ctree->SetBranchAddress("jet_btagged",jet_btagged);
  ctree->SetBranchAddress("jet_eta", jet_eta);
  ctree->SetBranchAddress("jet_phi", jet_phi);
  ctree->SetBranchAddress("jet_e", jet_e);
  ctree->SetBranchAddress("weight", &w);
  
  ctree->SetBranchStatus("*",0);
  
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  ctree->SetBranchStatus("jet_btagged",1);
  ctree->SetBranchStatus("jet_eta",1);
  ctree->SetBranchStatus("jet_phi",1);
  ctree->SetBranchStatus("jet_e", 1);
  ctree->SetBranchStatus("weight", 1);

  int nentries = ctree->GetEntries();
  
//The following loop creates a list pt_bjet_2jets that contains 2 b-jets  

  float pt_bjet_all[2] = {0};		//stores the pt of the events with greater than 2 b-jets in a 2-dimensional array.
  float phi_bjet_all[2] = {0};		//stores the phi angles of all of the jets in a 2-dimensional array
  float eta_bjet_all[2] = {0};  	//stores the eta angles of all of the jets in a 2-dimensional array
 
  float pt_w_all[20000][4] = {0};
  float phi_w_all[20000][4] = {0};
  float eta_w_all[20000][4] = {0};
 
  float x = 0;			//stores the (P_x,0 	+ P_x,1) 	values within each loop
  float y = 0;			//stores the (P_y,0 	+ P_y,1) 	values within each loop 
  float z = 0;			//stores the (P_z,0 	+ P_z,1) 	values within each loop
  float t = 0;			//stores the (P_total,1 + P_total_2) 	values within each loop
  float M_inv = 0;		//stores the invariant mass value within each loop
  
  float w_x = 0;
  float w_y = 0;
  float w_z = 0;
  float w_t = 0;
  float w_M_inv = 0;
  
  float w_2_x = 0;
  float w_2_y = 0;
  float w_2_z = 0;
  float w_2_t = 0;
  float w_2_M_inv = 0;
  
  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    int nbtags = 0;
    
    for (int j=0; j<n_jet; j++) {
      if(jet_btagged[j]>0) h_jet_pt_btag->Fill(jet_pt[j], w_2); 	//Creates histogram for b-jets
      
      if(jet_btagged[j]>0 && nbtags < 2){
	nbtags++;

	if (nbtags==1) h_leading_bjet_1->Fill(jet_pt[j]);		//Creates histogram for the leading b-jet
	if (nbtags==2) h_leading_bjet_2->Fill(jet_pt[j]);		//Creates histogram for the subleading b-jet

	pt_bjet_all[nbtags-1] = jet_pt[j];				
	phi_bjet_all[nbtags-1] = jet_phi[j];
	eta_bjet_all[nbtags-1] = jet_eta[j];
	
      } // btagging
      
      else
      {
	nbtags++;
	
	h_w_pt_btag->Fill(jet_pt[j], w_2);
	
	pt_w_all[i][nbtags-1] = jet_pt[j];
	phi_w_all[i][nbtags-1] = jet_phi[j];
	eta_w_all[i][nbtags-1] = jet_eta[j];
      }
	
    } // n_jets
    
      if (nbtags>=2) {

          
	x = pt_bjet_all[0] * cos(phi_bjet_all[0]) + pt_bjet_all[1] * cos(phi_bjet_all[1]);
	y = pt_bjet_all[0] * sin(phi_bjet_all[0]) + pt_bjet_all[1] * sin(phi_bjet_all[1]);
	z = pt_bjet_all[0] * sinh(eta_bjet_all[0]) + pt_bjet_all[1] * sinh(eta_bjet_all[1]);
	t = sqrt( pow(pt_bjet_all[0] * cos(phi_bjet_all[0]) , 2)
		    +pow(pt_bjet_all[0] * sin(phi_bjet_all[0]) , 2)
		    +pow(pt_bjet_all[0] * sinh(eta_bjet_all[0]) , 2)
	      )
	      +sqrt( pow(pt_bjet_all[1] * cos(phi_bjet_all[1]) , 2)
		    +pow(pt_bjet_all[1] * sin(phi_bjet_all[1]) , 2)
		    +pow(pt_bjet_all[1] * sinh(eta_bjet_all[1]) , 2)
	      );
	      
	M_inv = (sqrt ( pow(t , 2) - pow(x , 2) - pow(y , 2) - pow(z , 2) )) / 1000;
	h_M_inv->Fill(M_inv, w);
	
	////////Below is the inv mass calculation for W particles/////////
	
	w_x = pt_w_all[i][0] * cos(phi_w_all[i][0]) + pt_w_all[i][1] * cos(phi_w_all[i][1]);
	w_y = pt_w_all[i][0] * sin(phi_w_all[i][0]) + pt_w_all[i][1] * sin(phi_w_all[i][1]);
	w_z = pt_w_all[i][0] * sinh(eta_w_all[i][0]) + pt_w_all[i][1] * sinh(eta_w_all[i][1]);
	w_t = sqrt( pow(pt_w_all[i][0] * cos(phi_w_all[i][0]) , 2)
		    +pow(pt_w_all[i][0] * sin(phi_w_all[i][0]) , 2)
		    +pow(pt_w_all[i][0] * sinh(eta_w_all[i][0]) , 2)
	      )
	      +sqrt( pow(pt_w_all[i][1] * cos(phi_w_all[i][1]) , 2)
		    +pow(pt_w_all[i][1] * sin(phi_w_all[i][1]) , 2)
		    +pow(pt_w_all[i][1] * sinh(eta_w_all[i][1]) , 2)
	      );
	      
	w_M_inv = (sqrt ( pow(w_t , 2) - pow(w_x , 2) - pow(w_y, 2) - pow(w_z , 2) )) / 1000;
	h_w_inv->Fill(w_M_inv, w);
	
	
	
	w_2_x = pt_w_all[i][2] * cos(phi_w_all[i][2]) + pt_w_all[i][3] * cos(phi_w_all[i][3]);
	w_2_y = pt_w_all[i][2] * sin(phi_w_all[i][2]) + pt_w_all[i][3] * sin(phi_w_all[i][3]);
	w_2_z = pt_w_all[i][2] * sinh(eta_w_all[i][2]) + pt_w_all[i][3] * sinh(eta_w_all[i][3]);
	w_2_t = sqrt( pow(pt_w_all[i][2] * cos(phi_w_all[i][2]) , 2)
		    +pow(pt_w_all[i][2] * sin(phi_w_all[i][2]) , 2)
		    +pow(pt_w_all[i][2] * sinh(eta_w_all[i][2]) , 2)
	      )
	      +sqrt( pow(pt_w_all[i][3] * cos(phi_w_all[i][3]) , 2)
		    +pow(pt_w_all[i][3] * sin(phi_w_all[i][3]) , 2)
		    +pow(pt_w_all[i][3] * sinh(eta_w_all[i][3]) , 2)
	      );
	      
	w_2_M_inv = (sqrt ( pow(w_2_t , 2) - pow(w_2_x , 2) - pow(w_2_y, 2) - pow(w_2_z , 2) )) / 1000;
	h_w2_inv->Fill(w_2_M_inv, w);
	    
	} // jets in an event

    
//this loop creates all the jet_pt histogram
    for(int k = 0; k< n_jet; k++)
    {
      h_jet_pt->Fill( jet_pt[k], w_2);
    }
    
  }
  file->Close();
  TFile f("leading_jet_pt.root","recreate");

  h_jet_pt     ->Write();
  h_jet_pt_btag->Write();
  h_leading_bjet_1->Write();
  h_leading_bjet_2->Write();
  
  c2 = new TCanvas("c2", "Particle Jets from PP collision", 2000, 1000);
  c2 -> Divide (3,3);
  
  c2->cd(1);
  h_jet_pt->Draw();
  
  c2->cd(4);
  h_jet_pt_btag->Draw();
  
  c2->cd(7);
  h_leading_bjet_2->Draw();
  h_leading_bjet_1->Draw("same");
  h_leading_bjet_2->SetLineColor(2);
  auto legend = new TLegend(0.7,0.6,0.9,0.7);
  legend->SetHeader("b-Jets");
  legend->AddEntry(h_leading_bjet_1,"Leading b-Jet");
  legend->AddEntry("h_leading_bjet_2","Subleading b-Jet");
  legend->Draw();

  c2->cd(5);
  h_M_inv->Draw();
  
  c2->cd(3);
  h_w_pt_btag->Draw();
  
  c2->cd(6);
  h_w_inv->Draw();
  
  c2->cd(9);
  h_w2_inv->Draw();
 
  
  f.Close();
}