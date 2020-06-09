#define Analyzetw_cxx
#include "../include/Analyzetw.hpp"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void Analyzetw::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Analyzetw.C
//      root> Analyzetw t
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
   if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

   pMother1 = new TLorentzVector[nentries];
   b        = new TLorentzVector[nentries];
   pt       = new TLorentzVector[nentries];
   p3       = new TLorentzVector[nentries];
   pTmu    = new double[nentries];
   etamu    = new double[nentries];
   pTelec    = new double [nentries];
   etaElec   = new double [nentries];
   pTb    = new double [nentries];
   eta_b   = new double [nentries];
   pTnon_b    = new double [nentries];
   eta_non_b   = new double [nentries];
   nature   = new int[nentries];
   p3_PID   = new int[nentries];
   pMother1_PID = new int[nentries];
   b_PID = new int[nentries];


  for (int i=0; i<nentries; i++){
    nature[i] = 0;
    pTmu[i] = 0;
    pTelec[i] = 0;
    pTb[i] = 0;
    pTnon_b[i] = 0;
  }



   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

    for(int part = 0; part<kMaxParticle; part++)
    {

      if(Particle_PID[part]>0 && abs(Particle_PID[part])!=5 && Particle_Status[part]==-1)
      {
        nature[jentry] = 1; // One mother is a quark which isn't a b ==> bq t-channel
      }
      if(Particle_PID[part]<0 && abs(Particle_PID[part])!=5 && Particle_Status[part]==-1)
      {
        nature[jentry] = 2; // One mother is an antiquark which isn't a b ==> bqbar t-channel
      }
      if(abs(Particle_PID[part])==21 && Particle_Status[part]==-1)
      {
        nature[jentry] = 3; // One mother is a gluon ==> bg
      }

      //Fill the c particle 4-vector
      if(abs(Particle_PID[part])==4 && Particle_Status[part]==2)
      {
        pt[jentry].SetPx(Particle_Px[part]);
        pt[jentry].SetPy(Particle_Py[part]);
        pt[jentry].SetPz(Particle_Pz[part]);
        pt[jentry].SetE(Particle_E[part]);
      }
      //Fill the s,d particle 4-vector
      if( abs(Particle_PID[part]) < 4 && Particle_Status[part]==-1)
      {
        b[jentry].SetPx(Particle_Px[part]);
        b[jentry].SetPy(Particle_Py[part]);
        b[jentry].SetPz(Particle_Pz[part]);
        b[jentry].SetE(Particle_E[part]);
        b_PID[jentry] = abs(Particle_PID[part]);
      }
      //Fill the random mother 4-vector
      if(abs(Particle_PID[part])!=1 && abs(Particle_PID[part])!=3 && Particle_Status[part]==-1)
      {
        pMother1[jentry].SetPx(Particle_Px[part]);
        pMother1[jentry].SetPy(Particle_Py[part]);
        pMother1[jentry].SetPz(Particle_Pz[part]);
        pMother1[jentry].SetE(Particle_E[part]);
      }
      //Fill the random daughter 4-vector
      if(abs(Particle_PID[part])!=5 && Particle_Status[part]==1 && abs(Particle_PID[part])!=13 && abs(Particle_PID[part])!=11 && abs(Particle_PID[part])!=14 && abs(Particle_PID[part])!=12)
      {
        p3[jentry].SetPx(Particle_Px[part]);
        p3[jentry].SetPy(Particle_Py[part]);
        p3[jentry].SetPz(Particle_Pz[part]);
        p3[jentry].SetE(Particle_E[part]);
      }

      //Get transverse momentum and rapidity for the charged lepton
      if(abs(Particle_PID[part])==11 && Particle_Status[part]== 1 )
      {
        if (pTelec[jentry]>Particle_PT[part] || etaElec[jentry]<abs(Particle_Eta[part]))
        {
          pTelec[jentry] = Particle_PT[part];
          etaElec[jentry]= abs(Particle_Eta[part]);

        }
      }
      if(abs(Particle_PID[part])==13 && Particle_Status[part]== 1)
      {
        if (pTmu[jentry]>Particle_PT[part] || etamu[jentry]<abs(Particle_Eta[part]))
        {
          pTmu[jentry] = Particle_PT[part];
          etamu[jentry]= abs(Particle_Eta[part]);
        }

      }
        ///-------------A remplir--------------//
      if(abs(Particle_PID[part])==5  && (Particle_Status[part]==1 or Particle_Status[part]==-1))
      {
        if (pTb[jentry]>Particle_PT[part] || eta_b[jentry]<abs(Particle_Eta[part]))
        {
          pTb[jentry]    = Particle_PT[part];
          eta_b[jentry]  = abs(Particle_Eta[part]);
        }
      }
      else if(abs(Particle_PID[part])<=6  && Particle_Status[part]==1)
      {
        if (pTnon_b[jentry]>Particle_PT[part] || eta_non_b[jentry]<abs(Particle_Eta[part]))
        {
          pTnon_b[jentry] = Particle_PT[part];
          eta_non_b[jentry] = abs(Particle_Eta[part]);
        }
      }

    }


}}
