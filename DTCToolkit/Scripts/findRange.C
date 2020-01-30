#define findRange_cxx
#include "findRange.h"
#include "../GlobalConstants/Constants.h"
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
#include <TSpline.h>
#include <iostream>
#include <fstream>
#include <TPad.h>
#include <TMath.h>

using namespace std;

vector<Float_t> findRange::Run(Double_t energy, Double_t sigma_mev, Int_t mm, Int_t degrader, Int_t fileIdx)
{
   vector<Float_t> returnValues;
   returnValues.reserve(100);
   if (fChain == 0) return returnValues;
   if (fChain->GetEntries() == 0) return returnValues;


   Float_t energyBeforeTrackerNoDegrader = 229.88; // proton
   // Float_t energyBeforeTrackerNoDegrader = 916.46 ; // Helium

   Long64_t nentries = fChain->GetEntriesFast();
   Bool_t useDegrader = true; 

   TFile *fEkin = new TFile(Form("Data/MonteCarlo/DTC_Full_Final_Degrader%03dmm_230MeV_psa.root", degrader));
   TTree *tEkin = (TTree*) fEkin->Get("PhaseSpace");

   Float_t mean = 0;
   Float_t sigma = 0;
   Float_t N = 0;
   Float_t Ekine;
   tEkin->SetBranchAddress("Ekine", &Ekine);

   for (Int_t i=0; i<tEkin->GetEntriesFast(); i++) {
      tEkin->GetEntry(i);
      if (Ekine<1) continue;
      mean += Ekine;
      N++;
   }
   mean /= N;

   for (Int_t i=0; i<tEkin->GetEntriesFast(); i++) {
      tEkin->GetEntry(i);
      if (Ekine<1) continue;
      sigma += pow(mean-Ekine,2);
   }
   sigma = sqrt(sigma/N);

   TH1I *hEkine = new TH1I("hEkine", "Energy spectrum before first tracker;Energy [MeV];Entries", 200, fmax(0,mean-4*sigma), mean+4*sigma);
   for (Int_t i=0; i<tEkin->GetEntriesFast(); i++) {
      tEkin->GetEntry(i);
      if (Ekine<1) continue;
      hEkine->Fill(Ekine);
   }

   TF1 *fitEnergy = new TF1("fitEnergy", "gaus");
   fitEnergy->SetParameters(100, mean, sigma);
   hEkine->Fit("fitEnergy");

   Float_t energyBeforeTracker = fitEnergy->GetParameter(1);
   Float_t energyBeforeTrackerSigma = fitEnergy->GetParameter(2);

   printf("From phase space: Energy before first tracker is %.2f MeV +- %.2f Mev\n", energyBeforeTracker, energyBeforeTrackerSigma);

   TCanvas *c2 = new TCanvas("c2", "Ranges and energies", 500, 500);
   c2->Divide(2,1);
   gStyle->SetOptStat(0);

   Int_t    lastEvent = -1;
   Float_t  lastRange = 0;
   Int_t    lastID = -1;
   Float_t  dE = 0, dE2 = 0;

   TH1F *hFirstRange = new TH1F("firstRange", "firstRange", 1000, 0, 400);
//   TH1F *hFirstEnergy = new TH1F("firstEnergy", "firstEnergy", 1000, 0, 1000);
//   TH1F *hFirstEnergy2 = new TH1F("firstEnergy2", "Energy loss in trackers", 1000, 0, 5);

   Float_t allRanges[150000];
//   Float_t allEnergies[150000];
   Int_t rangeIdx = 0, energyIdx = 0;

   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   if (nentries == 0) return returnValues;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      fChain->GetEntry(jentry);
      
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }

      if (lastID < 0) lastID = eventID;
      if (parentID == 0) {
/*      
         if (baseID == 0) { // inside degrader + first layers
//            printf("dE = %.2f + %.2f\n", dE, edep);
            dE += edep;
         }
   
         else if (volumeID[1] == 1) {
            dE2 += edep;
         }

         else if (dE > 0) {
            hFirstEnergy->Fill(energy-dE);
            hFirstEnergy2->Fill(dE2);
            allEnergies[energyIdx++] = energy - dE;
            dE = 0;
            dE2 = 0;
         }
*/
         if (eventID != lastID) {
            hFirstRange->Fill(lastRange - 225);
            allRanges[rangeIdx++] = lastRange - 225;
         }

         lastRange = posZ;
         lastID = eventID;
      }
   }

   float firstRange = hFirstRange->GetBinCenter(hFirstRange->GetMaximumBin());
//   float firstEnergy = hFirstEnergy->GetBinCenter(hFirstEnergy->GetMaximumBin());

   printf("Found %d proton histories in file.\n", lastID);
//   printf("From first search, range = %.2f mm and energy = %.2f MeV\n", firstRange, firstEnergy);
   printf("From first search, range = %.2f mm\n", firstRange);
  
   // WHAT about WEPL range = Residual WET ?? Absolutely no inverse calculations of the spline

   Float_t xfrom = firstRange - 10, xto = firstRange + 10;
   TH1F *hRange = new TH1F("hRange", "Projected range in DTC;Range [mm];Entries", 500, xfrom,xto);
   for (Int_t i=0; i<rangeIdx; i++) hRange->Fill(allRanges[i]);

//   Float_t energyfrom = firstEnergy - 15, energyto = firstEnergy + 15;
//   TH1F *hEnergyAtInterface = new TH1F("hEnergyAtInterface", "Remaining energy after degrader;Energy [MeV];Entries", 500, energyfrom, energyto);
//   for (Int_t i=0; i<energyIdx; i++) hEnergyAtInterface->Fill(allEnergies[i]);

   c2->cd(1);
   // The fitting parameters below
   TF1 *fRange = new TF1("fit_range", "gaus", xfrom, xto);
   fRange->SetLineWidth(3);
   fRange->SetParameters(hRange->GetMaximum(), firstRange, 1.5);
   hRange->Fit("fit_range", "B");//, "M,W,B", "", xfrom, xto);
   Float_t fR = fRange->GetParameter(1);
   Float_t fRS = fRange->GetParameter(2);
   printf("Estimated range from Gaussian fitting = %.3f +- %.3f\n", fR, fRS); 
   
   Float_t hR = hRange->GetMean();
   Float_t hRS = hRange->GetStdDev();
   printf("Estimated range from bin counting = %.3f +- %.3f\n", hR, hRS);
   
   c2->cd(1);
   hRange->SetFillColor(kBlue-7);
   hRange->SetLineColor(kBlack);
   hRange->Draw();
  
   c2->cd(2); 
   /*
   TF1 *fRemainingEnergy = new TF1("fRemainingEnergy", "gaus");
   fRemainingEnergy->SetLineWidth(3);
   fRemainingEnergy->SetParameters(hEnergyAtInterface->GetMaximum(), firstEnergy, 3);
   hEnergyAtInterface->Fit("fRemainingEnergy", "B");
   Float_t fE = fRemainingEnergy->GetParameter(1);
   Float_t fES = fRemainingEnergy->GetParameter(2);

   printf("Estimated remaining energy from Gaussian fitting: %.2f +- %.2f MeV.\n", fE, fES);

   returnValues.push_back(hR);
   returnValues.push_back(hRS);
   returnValues.push_back(fE);
   returnValues.push_back(fES);

   c2->SaveAs(Form("OutputFiles/straggling/straggling_absorber%dmm_degrader%.0fmm.png", mm, run_degraderThickness));

   c2->cd(2);
   hEnergyAtInterface->SetFillColor(kBlue-7);
   hEnergyAtInterface->Draw();
   
   c2->cd(3);
   hFirstEnergy2->SetFillColor(kBlue-7);
   hFirstEnergy2->Draw();
*/
   c2->cd(2);
   hEkine->SetFillColor(kBlue-7);
   hEkine->SetLineColor(kBlack);
   hEkine->Draw();
   
   // helium: 331.7; proton: 330.9
   // Find new numbers after accounting for energy loss due to air!
   // New calibration: pol2 where E=0
   // Helium: 332.3 (+0.6 mm diff)
   // Proton: 
   
   std::ofstream filename(Form("OutputFiles/findManyRangesDegrader_final_idx%d.csv", fileIdx));// , std::ofstream::out | std::ofstream::app);
   filename << degrader << " " << 330.9 - degrader << " " << hR << " " << hRS << " " << energyBeforeTracker << " " << energyBeforeTrackerSigma << endl; 
   
   delete c2;
   delete hRange;
   delete hFirstRange;
//   delete hFirstEnergy;
//   delete fRemainingEnergy;
   delete fRange;
//   delete hEnergyAtInterface;
   
   return returnValues;
}
