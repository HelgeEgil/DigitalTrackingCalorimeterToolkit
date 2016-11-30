#define findFilterEnergyLoss_cxx

#include <iostream>
#include <string.h>

#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TList.h>
#include <TSpectrum.h>
#include <TPolyMarker.h>
#include <TMath.h>

#include "include/findFilterEnergyLoss.h"

using namespace std;

void findFilterEnergyLoss::BinLogY(TH2 *h) {
   TAxis *axis = h->GetYaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;

   Axis_t *new_bins = new Axis_t[bins+1];

   for (int i=0; i <= bins; i++) {
      new_bins[i] = TMath::Power(10, from + i * width);
   }

   axis->Set(bins, new_bins);
   delete new_bins;
}


void findFilterEnergyLoss::Loop()
{
   Int_t energy = run_energy;
   
   Double_t alpha = 0.0014467;
   Double_t alpha_prime = 0.203815;
   Double_t p = 1.707283;
   
   Float_t run_thickness = 1;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
   TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
   TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
   
   TH1F *hEnergyLoss = new TH1F("hEnergyLoss", Form("Energy loss for a %d MeV proton beam on a %.1f mm aluminum foil", energy, run_thickness), 400, 0, 20);
   TH1F *hEnergy= new TH1F("hEnergy", Form("Final energy for a %d MeV proton beam on a %.1f mm aluminum foil", energy, run_thickness), 5000, 0, 250);
   TH2F *hEnergy2D = new TH2F("hEnergy2D", Form("Final energy for a %d MeV proton beam on focal vs x", energy), 128, -20, 20, 200, 0, 200);
   
   Int_t lastID = -1;
   Float_t sum_edep = 0;

   Float_t lastX = 0;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
   
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;

//       if (abs(posX) < 1 && abs(posY) < 1) {
//          continue;
//       }     

      if (eventID != lastID) {
         hEnergyLoss->Fill(sum_edep);
         hEnergy->Fill(energy - sum_edep);
         hEnergy2D->Fill(lastX, energy - sum_edep);
         sum_edep = edep;
      }

      else
         sum_edep += edep;

      lastID = eventID;
      lastX = posX;
   }
   
   TF1 *fEnergy = new TF1("f_Energy", "gaus(0)", 0, energy+10);
   TF1 *fEnergyLoss = new TF1("f_EnergyLoss", "gaus(0)", 0, 30);

// fEnergy->SetParameters(100, energy-10, 0.3);
// fEnergyLoss->SetParameters(160, 10, 0.3);
   
   hEnergy->Fit("f_Energy");
   hEnergyLoss->Fit("gaus", "M");

// fEnergyLoss->SetParLimits(1, finalLoss*0.8, finalLoss*1.2);
// fEnergyLoss->SetParLimits(2, 0.05, 0.5);

   c1->cd();
   hEnergy->SetXTitle("Final energy");
   hEnergy->SetYTitle("Number of particles");
   hEnergy->SetFillColor(kYellow - 10);
   hEnergy->SetLineColor(kBlack);
   hEnergy->Draw();
   
   c2->cd();
   hEnergyLoss->SetXTitle("Total energy loss");
   hEnergyLoss->SetYTitle("Number of particles");
   hEnergyLoss->Draw();
   
   c3->cd();
   hEnergy2D->SetXTitle("X position [mm]");
   hEnergy2D->SetYTitle("Final energy [MeV]");
// hEnergy2D->Draw("COLZ");
   

   cout << "For a " << energy << " MeV beam, the energy loss is " << Form("%.2f", fEnergyLoss->GetParameter(1)) << " +- " << Form("%.2f", fEnergyLoss->GetParameter(2)) << " MeV.\n";
}
