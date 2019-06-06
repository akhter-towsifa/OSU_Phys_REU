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
  
  auto h_jet_pt_btag = new TH1I("jet_btagged","Jet p_{T} b-tagged;X-Axis",nbins,0.0,2.4e6);
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


  nbtags = 0; 
  Int_t bj[] = {};
  for (int j=0; j<nentries;j++){
    ctree->GetEntry(j);
    if(jet_btagged[j]>0){
      nbtags++;
      bj[nbtags]=j;
      std::cout<<"passed..."<<std::endl;
      h_jet_pt_btag->Fill(jet_pt[bj[nbtags]],w);
    }
  }

  h_jet_pt_btag->Write();
  
  f.Write(); 
  f.Close();
}
