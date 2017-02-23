#define findRange_cxx
#include "findRange.h"
#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

void findRange::Run()
{
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TCanvas *c6 = new TCanvas("c6", "hEnergyAtInterface", 800, 600);

   printf("RUNNING WITH ENERGY %.2f.\n", run_energy);

   TH1F *hEnergyAtInterface = new TH1F("hEnergyAtInterface", "Remaining energy after degrader;Energy [MeV];Entries", 500, 2, 3);

   gStyle->SetOptStat(0);

   Int_t lastEvent = -1;

   Int_t lastID = -1;
   Float_t dE = 0;

   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (lastID < 0) lastID = eventID;

      if (parentID == 0) {
      
         if (eventID == lastID) {
            dE += edep;
         }
         
         else {
            hEnergyAtInterface->Fill(dE);
            dE = edep;
         }
         
         lastID = eventID;
      }
   }
   
   c6->cd();
      hEnergyAtInterface->SetFillColor(kBlue-7);
      hEnergyAtInterface->SetLineColor(kBlack);
      hEnergyAtInterface->Draw();
   
}
