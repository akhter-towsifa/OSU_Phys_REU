#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TBranch.h"
#include "TMath.h"
#include "TImage.h"
#include "TCanvas.h"
#include "TString.h"
#include "TAttFill.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLatex.h"
#include <math.h>
#include <cmath>
#include <TChain.h>


void hist_btagged(){
  
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptStat("nemi");
  int nbins = 1000;
  double w = 1;
  
  auto h_jet_pt_btag = new TH1I("jet_btagged","Jet p_{T} b-tagged;X",nbins,0.0,5);
  auto h_jet_pt = new TH1F("jet_pt","Jet_P_{T}; X; Counts",nbins,0.0,2.4e6);
  auto file = TFile::Open("bbww_x1000_s170.root");

  TTree* ctree = (TTree*)file->Get("CollectionTree");

  Int_t           n_jet;
  Float_t         jet_pt[14];   //[n_jet]
  Int_t           jet_btagged[14];   //[n_jet]
  
  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt",jet_pt);
  ctree->SetBranchAddress("jet_btagged",jet_btagged);
  ctree->SetBranchStatus("*",0);
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  ctree->SetBranchStatus("jet_btagged",1);

  int nentries = ctree->GetEntries();
  
  TFile f("btagged_pt.root","recreate");

  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    for(uint j = 0; j< n_jet; j++){
      h_jet_pt->Fill( jet_pt[j],w ); 
    }
  }

  int64_t bj[] = {};	
  
  std::cout << "total number of entries: " << nentries << std::endl;
  
  for (int i=0; i<nentries;i++){
    ctree->GetEntry(i);
    int nbtags = 0; 
    for (int j=0; j<n_jet; j++) {
      if(jet_btagged[j]>0){         
	bj[nbtags]=j;
        nbtags++;

	if (nbtags == 1) {
	  //std::cout << "ptbj_" << nbtags << std::endl; 
	  int pt_bj1 = jet_pt[bj[0]];
	  //int pt_bj2 = jet_pt[bj[1]];
	  //std::cout <<"1"<< "\t" << pt_bj1 <<std::endl;
	}
	else if (nbtags == 2) {
	  int pt_bj2 = jet_pt[bj[1]];
	  //std::cout <<"2"<< "\t" << pt_bj2 <<std::endl;
	}
	else if (nbtags == 3) {
	  int pt_bj3 = jet_pt[bj[2]];
	  //std::cout << "3"<< "\t" <<pt_bj3 <<std::endl;
	}
	else if (nbtags == 4) {
	  int pt_bj4 = jet_pt[bj[3]];
	  //std::cout << "4"<< "\t" << pt_bj4 <<std::endl;
	}
	else if (nbtags == 5) {
	  int pt_bj5 = jet_pt[bj[4]];
	  //std::cout << "5"<< "\t" << pt_bj5 <<std::endl;
	}
	for (int m=1; m<6; m++) {
	  h_jet_pt_btag->Fill(jet_pt[m],w);
	}
      }
    }
  }

  h_jet_pt_btag->Write();
  
  f.Write(); 
  f.Close();
}

// << "\t" << pt_bj2 