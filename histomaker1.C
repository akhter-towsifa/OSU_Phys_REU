#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <map>
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


void histomaker1(){
  
  gStyle->SetOptStat(kFALSE);
  
  int nbins = 1000;
  double w = 1;
  
  auto h_jet_pt = new TH1F("jet_pt","Jet_P_{T}",nbins,0.0,2.4e6);
  auto file = TFile::Open("bbww_x1000_s170.root");
  std::cout<<"this part worked 1"<<std::endl;
  
  //TTreeReaderValue<std::vector<float>> jet_pt(myReader,"jet_pt");
  TTree* ctree = (TTree*)file->Get("CollectionTree");

  Int_t           n_jet;
  Float_t         jet_pt[14];   //[n_jet]
  ctree->SetBranchAddress("n_jet", &n_jet);
  ctree->SetBranchAddress("jet_pt",jet_pt);
  ctree->SetBranchStatus("*",0);
  ctree->SetBranchStatus("n_jet",1);
  ctree->SetBranchStatus("jet_pt",1);
  std::cout<<"this part worked 2"<<std::endl;
  float nentries = ctree->GetEntries();
  
  
  TFile f("test.root","recreate");
  
  std::cout<<"this part worked 3"<<std::endl;
  

  for (int i=0; i<nentries;i++){
    std::cout<<"passed..."<<std::endl;
    ctree->GetEntry(i);
    for(uint j = 0; j< n_jet; j++){
      h_jet_pt->Fill( jet_pt[j],w );
      
    }
  }

  
  h_jet_pt->Write();
  
  f.Write();
  f.Close(); 
  
  
}
