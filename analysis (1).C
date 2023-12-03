#define analysis_cxx
#include "analysis.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L analysis.C
//      root> analysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

// At this point we declare the output root file with the histograms we need
   TFile *rfile  = new TFile( "OutputHistograms.root","RECREATE");
   TH1F *h_pt_mu = new TH1F("h_pt_mu","Muon pT ", 100,0.,100.);
   TH1F *h_eta_mu = new TH1F("h_eta_mu","Muon eta ", 54,-2.7,2.7);
   TH1F *h_phi_mu = new TH1F("h_phi_mu","Muon phi ", 50,-acos(-1.),acos(-1.));

   TH1F *h_pt_el = new TH1F("h_pt_el","Electron pT ", 100,0.,100.);
   TH1F *h_eta_el = new TH1F("h_eta_el","Electron eta ", 50,-2.5,2.5);
   TH1F *h_phi_el = new TH1F("h_phi_el","Electron phi ", 50,-acos(-1.),acos(-1.));


   TH1F *h_mZ_mu = new TH1F("h_mZ_mu","Invariant mass mu+mu-", 200,0,200.);
   TH1F *h_mZ_el = new TH1F("h_mZ_el","Invariant mass e+e- ", 200,0.,200.);

//------------------------------------------------------------------------------



   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
// At this point we loop over the collision events of the original tree
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

// Show analysis progress
    if(jentry%10000==0)std::cout << "Analysed a total of: " << jentry << " events" << std::endl;

// At this point we construct the variables we need and apply our criteria
// As an example we loop over all leptons in the event 
// and we fill histograms with their pT (electrons and muons separately) if their pT exceeds 25 GeV 
// All quantities in the tree are expressed in MeV
	  float MeV2GeV = 0.001;
//Scale factors to improve MC simulation with actual detector real data agreement
	  Float_t m_mcWeight = mcWeight;
	  //cout << m_mcWeight <<" "<< XSection<<" "<<SumWeights<<endl;
      float scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER*scaleFactor_PILEUP; 
      float weight =1.;  
//Total weight for MC event  = MC event weight  mcWeight * scale factor  
         
      if(mcWeight!=0.)weight = scaleFactor*(mcWeight/SumWeights)*(XSection*10000.);	 // weight used to fill all MC histograms 
	  
	  
	  int n_mu=0;
	  int n_el=0;
	  float mu_v[5][10];
	  float el_v[5][10];
	  //if(!trigE)continue;   // Trigger for electrons
	  if(!trigM)continue;   // Trigger for muons

	for (int i=0; i<lep_n;i++)
	    {
		   //electrons
		   if(lep_type->at(i)==11 && lep_pt->at(i)>25000.)
		   {
			   h_pt_el->Fill((lep_pt->at(i))*MeV2GeV, weight);
			   h_eta_el->Fill(lep_eta->at(i), weight);
			   h_phi_el->Fill(lep_phi->at(i), weight);
			   el_v[n_el][0]=(lep_pt->at(i))*MeV2GeV; // lepton pT
			   el_v[n_el][1]=el_v[n_el][0]*cos(lep_phi->at(i)); //lepton px
			   el_v[n_el][2]=el_v[n_el][0]*sin(lep_phi->at(i)); //lepton py
				// calculate polar angle theta from pseudorapidity eta=-ln(tan(theta/2))
			   float theta=2*atan(exp(-lep_eta->at(i)));
			   el_v[n_el][3]=el_v[n_el][4]/tan(theta);  //lepton pz
			   el_v[n_el][4]=(lep_E->at(i))*MeV2GeV;    //lepton E 
			   el_v[n_el][5]=lep_charge->at(i); //lepton charge
			   el_v[n_el][6]=(lep_ptcone30->at(i))*MeV2GeV;   //track isolation
			   el_v[n_el][7]=(lep_etcone20->at(i))*MeV2GeV;   //calorimeter isolation
			   el_v[n_el][8]=lep_tracksigd0pvunbiased->at(i);  // Impact parameter significance in the R-phi plane
			   el_v[n_el][9]=lep_z0->at(i);  // lepton extrpolated position on the Z axis
			   n_el++; 			   
		   }
		   
		   
		   //muons
		   if(lep_type->at(i)==13 && lep_pt->at(i)>25000.)
		   {
			   h_pt_mu->Fill(lep_pt->at(i)/1000., weight);
			   h_eta_mu->Fill(lep_eta->at(i), weight);
			   h_phi_mu->Fill(lep_phi->at(i), weight);
			   mu_v[n_mu][0]=(lep_pt->at(i))*MeV2GeV; // lepton pT
			   mu_v[n_mu][1]=el_v[n_el][0]*cos(lep_phi->at(i)); //lepton px
			   mu_v[n_mu][2]=el_v[n_el][0]*sin(lep_phi->at(i)); //lepton py
				// calculate polar angle theta from pseudorapidity eta=-ln(tan(theta/2))
			   float theta=2*atan(exp(-lep_eta->at(i)));
			   mu_v[n_mu][3]=el_v[n_el][4]/tan(theta);  //lepton pz
			   mu_v[n_mu][4]=(lep_E->at(i))*MeV2GeV;    //lepton E 
			   mu_v[n_mu][5]=lep_charge->at(i); //lepton charge
			   mu_v[n_mu][6]=(lep_ptcone30->at(i))*MeV2GeV;   //track isolation
			   mu_v[n_mu][7]=(lep_etcone20->at(i))*MeV2GeV;   //calorimeter isolation
			   mu_v[n_mu][8]=lep_tracksigd0pvunbiased->at(i);  // Impact parameter significance in the R-phi plane
			   mu_v[n_mu][9]=lep_z0->at(i);  // lepton extrpolated position on the Z axis
			   n_mu++; 			   

		   }
	    }

//------------------------------------------------------------------------------



 

   }


// At this point we write our histograms to the file and we close the file to finish
   h_pt_mu->Write();
   h_eta_mu->Write();
   h_phi_mu->Write();
   h_pt_el->Write();
   h_eta_el->Write();
   h_phi_el->Write();
   rfile->Close();
//------------------------------------------------------------------------------

}
