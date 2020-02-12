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

   Float_t allRanges[150000];
   Int_t rangeIdx = 0;

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
         if (eventID != lastID) {
            hFirstRange->Fill(lastRange - 225);
            allRanges[rangeIdx++] = lastRange - 225;
         }

         lastRange = posZ;
         lastID = eventID;
      }
   }

   float firstRange = hFirstRange->GetBinCenter(hFirstRange->GetMaximumBin());

   printf("Found %d proton histories in file.\n", lastID);
   printf("From first search, range = %.2f mm\n", firstRange);
  
   Float_t xfrom = firstRange - 10, xto = firstRange + 10;
   TH1F *hRange = new TH1F("hRange", "Projected range in DTC;Range [mm];Entries", 500, xfrom,xto);
   for (Int_t i=0; i<rangeIdx; i++) hRange->Fill(allRanges[i]);

   c2->cd(1);
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
   hEkine->SetFillColor(kBlue-7);
   hEkine->SetLineColor(kBlack);
   hEkine->Draw();
   
   // helium: 331.7; proton: 330.9
   // Find new numbers after accounting for energy loss due to air!
   // proton: 333.7 -> DELTA: 2.8
   // Helium: 332.3 -> DELTA: 0.6
   
   std::ofstream filename(Form("OutputFiles/findManyRangesDegrader_final_Proton_idx%d.csv", fileIdx));// , std::ofstream::out | std::ofstream::app);
   filename << degrader << " " << 330.9 - degrader << " " << hR << " " << hRS << " " << energyBeforeTracker << " " << energyBeforeTrackerSigma << endl; 
   
   delete c2;
   delete hRange;
   delete hFirstRange;
   delete fRange;
   
   return returnValues;
}
