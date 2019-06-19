#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"

//This code graphs the invariant mass of regular b-jets // i am using the same ttbar monte carlo simulated data as Merrick Lavinsky is.
void ttbar_3() {
  
  int w = 1; 					//weight of the bins
  int nbin = 1000;				//bin size
  
  TFile *f = new TFile("ttbar.root");		//Creates a pointer to the main data file
  TTree *t = (TTree*)f->Get("nominal");		//Creates a pointer to the tree
  auto h = new TH1F("plot", "Invariant Mass of the regular b-jet; Mass [GeV]; Events", nbin, 0, 1000000);
						//Creates a pointer to the histogram to be created.
  
  vector<float> *jet_pt		=0;
  vector<char>	*jet_btagged	=0;
  vector<float>	*jet_eta	=0;
  vector<float>	*jet_phi	=0;
  
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("jet_isbtagged_MV2c10_85", &jet_btagged);
  t->SetBranchAddress("jet_eta", &jet_eta);
  t->SetBranchAddress("jet_phi", &jet_phi);
  
  int nentries = t->GetEntries();
  
  vector<vector<float> > pt_bjet_all(981139, vector<float>(2));		//stores all the pt bjets with more than 1 btags
  vector<vector<float> > phi_bjet_all(981139, vector<float>(2));	//stores all the phi angles with more than 1 btags
  vector<vector<float> > eta_bjet_all(981139, vector<float>(2));	//stores all the eta angles with more than 1 btags
  
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

//     cout << i << "\tNumber of jets: " << n_jets << endl;
    
    for (int j=0; j<n_jets; j++)
    {
      if ((int)(*jet_btagged)[j]==1)
      {
	pt_bjet_all[i][n_bjets] = (float)(*jet_pt)[j];
	phi_bjet_all[i][n_bjets] = (float)(*jet_phi)[j];
	eta_bjet_all[i][n_bjets] = (float)(*jet_eta)[j];
// 	cout << "passes this point" << endl;
	n_bjets++;
      }
    }
      
    if (n_bjets >= 2)
    {
      cout << i << endl;
//       cout << "n_jets\t\t" << n_jets << endl;
//       
//       cout << "\t\tThe first pt_bjet\t\t" << pt_bjet_all[i][0] << endl;
//       cout << "\t\tThe second pt_bjet\t\t" << pt_bjet_all[i][1] << endl;
//       
//       cout << "\t\tThe first phi_bjet\t\t" << phi_bjet_all[i][0] << endl;
//       cout << "\t\tThe second phi_bjet\t\t" << phi_bjet_all[i][1] << endl;
//       
//       cout << "\t\tThe first eta_bjet\t\t" << eta_bjet_all[i][0] << endl;
//       cout << "\t\tThe second eta_bjet\t\t" << eta_bjet_all[i][1] << endl;
      
//       cout<< "about to calculate m_inv" << endl;
     
      x = pt_bjet_all[i][0] * cos(phi_bjet_all[i][0]) + pt_bjet_all[i][1] * cos(phi_bjet_all[i][1]);
      y = pt_bjet_all[i][0] * sin (phi_bjet_all[i][0]) + pt_bjet_all[i][1] * sin (phi_bjet_all[i][1]);
      
//       cout << "\t\t1st parameter in z\t\t" <<  sinh(eta_bjet_all[i][0]) << endl;
//       cout << "\t\t2nd parameter in z\t\t" <<  sinh(eta_bjet_all[i][1]) << endl;
      
      z = pt_bjet_all[i][0] * sinh(eta_bjet_all[i][0]) + pt_bjet_all[i][1] * sinh(eta_bjet_all[i][1]);

//       cout << "another checkpoint" << endl;      
      m = sqrt( pow(pt_bjet_all[i][0] * cos(phi_bjet_all[i][0]) , 2)
		+pow(pt_bjet_all[i][0] * sin(phi_bjet_all[i][0]) , 2)
		+pow(pt_bjet_all[i][0] * sinh(eta_bjet_all[i][0]) , 2)
	   )
	  +sqrt( pow(pt_bjet_all[i][1] * cos(phi_bjet_all[i][1]) , 2)
		+pow(pt_bjet_all[i][1] * sin(phi_bjet_all[i][1]) , 2)
		+pow(pt_bjet_all[i][1] * sinh(eta_bjet_all[i][1]) , 2)
	   );
	  
      
      M_inv = (sqrt ( pow(m , 2) - pow(x , 2) - pow(y , 2) - pow(z , 2) ));
      
//       cout << "calculated M-inv\t\t"<< M_inv <<  endl;
      
      h->Fill(M_inv, w);
      
      cout << "end of this loop" << endl;
    }
    
  }
  
  h->Draw();
  TFile file("ttbar_invariant_mass_3.root", "RECREATE");
  
  file.Write();
  file.Close();
  
  
}
