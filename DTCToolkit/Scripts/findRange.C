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

vector<Float_t> findRange::Run(Double_t energy, Double_t sigma_mev)
{
   vector<Float_t> returnValues;
   if (fChain == 0) return returnValues;
   if (fChain->GetEntries() == 0) return returnValues;

   Float_t degraderThickness = run_degraderThickness;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Bool_t useDegrader = (degraderThickness > 0) ? true : false;

   // Load phase space spline
   Double_t phaseSpaceDegraderthickness[300];
   Double_t phaseSpaceEnergy[300];
   Double_t dt, e, es;
   Int_t idx = 0;
   ifstream in;
   in.open("../Data/Ranges/EnergyAfterDegraderPSTAR.csv");

   while (1) {
      in >> dt >> e;
      if (!in.good()) break;
      phaseSpaceDegraderthickness[idx] = dt;
      phaseSpaceEnergy[idx] = e;
   }
   in.close();

   TSpline3 *phaseSpaceSpline = new TSpline3("phaseSpaceSpline", phaseSpaceDegraderthickness, phaseSpaceEnergy, idx);

   run_energy = phaseSpaceSpline->Eval(run_degraderThickness);

   TCanvas *c2 = new TCanvas("c2", "Ranges and energies", 1200, 900);
//   c2->Divide(2, 1, 0.001, 0.001);i

   Int_t nbinsx = 500;
   // 2 mm: 0.0096, 1.784
   // 3 mm: 0.0097, 1.7825
   // 4 mm: 0.0098, 1.7806
   // H20:  0.0239, 1.7548

   Float_t a = 0.0098, p = 1.7806;
   Float_t aw = 0.0239, pw = 1.7548;

   Float_t expectedRange = a * pow(run_energy, p);
   Float_t xfrom = expectedRange - 15;
   if (xfrom < 0) xfrom = 0;
   Float_t xto = expectedRange + 15;

   Float_t x_compensate = 0;

   Int_t energyFrom = run_energy - 15;
   Int_t energyTo = run_energy + 15;
   if (run_energy > 70) {
      energyFrom = run_energy - 15;
      energyTo = run_energy + 25;
   }

   printf("RUNNING WITH ENERGY %.2f.\n", run_energy);

   TH1F *hZ = new TH1F("hZ", "Z profile", nbinsx/3, xfrom + x_compensate, xto + x_compensate);
   TH1F *hRange = new TH1F("hRange", "Projected range in DTC", nbinsx, xfrom + x_compensate, xto + x_compensate);
   TH1F *hRangeWEPL = new TH1F("hRangeWEPL", "Primary ranges in WEPL", nbinsx, xfrom*1.9 + x_compensate, xto*2.5 + x_compensate);
   TH1F *hTracklength = new TH1F("hTracklength", "Beam range straggling", nbinsx, -20 , 20);
   TH1F *hActualTracklength = new TH1F("hActualTracklength", "Ranges in DTC deconvoluted from beam energy straggling", nbinsx, xfrom + x_compensate, xto + x_compensate);
   TH1F *hStepLength = new TH1F("hStepLength", "Steplenghths", 1000, 0, 1);
   TH1F *hEnergyAtInterface = new TH1F("hEnergyAtInterface", "Remaining energy after degrader;Energy [MeV];Entries", 150, energyFrom, energyTo);

   gStyle->SetOptStat(0);

   Int_t lastEvent = -1;

   Float_t lastRange = 0;
   Int_t lastID = -1;
   Float_t lastX = 0;
   Float_t lastY = 0;
   Float_t lastZ = 0;
   Float_t firstX, firstY, firstZ;
   Float_t dE = 0;
   Float_t dTL = 0;
   Float_t dE_random = 0;
   Float_t thisEnergy = 0;
   Float_t thisRange = 0;
   Int_t ignoreID = -5;

   Float_t tl = 0;
   Int_t n = 0;
   Char_t lastProcessName[17];
   
   TRandom3 *gRandom = new TRandom3();

   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   
   firstX = posX;
   firstY = posY;
   firstZ = posZ;

   Int_t lastP = 0;
   
   if (nentries == 0) return returnValues;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (lastID < 0) lastID = eventID;

      if (parentID == 0) {
      
//         hStepLength->Fill(stepLength);

         Float_t z = posZ;
         Float_t y = posY;
         Float_t x = posX;
      
//         hZ->Fill(z + x_compensate);
         n++;

         if (useDegrader) {
            if (baseID == 0) { // inside degrader
               dE += edep;
            }

            else if (dE > 0) {
               thisEnergy = 250 - dE;
               thisRange = aw * pow(run_energy, pw) - aw * pow(thisEnergy, pw);
               hEnergyAtInterface->Fill(thisEnergy);
               dE = 0;
            }
         }

         if (eventID != lastID) {
            n = 0;
            
            Float_t diff = sqrt( pow(firstX - lastX, 2) + pow(firstY - lastY, 2) + pow(0 - lastZ, 2));

            hRange->Fill(lastRange);
            Float_t weplRange = aw / a * pow(lastRange / aw, 1 - pw/p) * lastRange;
            Float_t dtcRange = degraderThickness - thisRange + weplRange;
            hActualTracklength->Fill(dtcRange);
            hRangeWEPL->Fill(weplRange);
            hTracklength->Fill(degraderThickness - thisRange);
            if (lastProcessName[0] == 'P') { lastP++; }

            firstX = posX;
            firstY = posY;
            firstZ = posZ;
            tl = 0;
         }

         else if (jentry>0) {
            Float_t diff = sqrt( pow(x - lastX, 2) + pow(y - lastY, 2) + pow(z - lastZ, 2));
            tl += diff;
         }

         lastRange = posZ;
         lastX = posX;
         lastID = eventID;
         lastY = posY;
         lastZ = posZ;

         for (Int_t j=0; j<17; j++) {
            lastProcessName[j] = processName[j];
         }
      }
   }
   
   c2->cd(1);
  
   // The fitting parameters below
   TF1 *fRange = new TF1("fit_range", "gaus", xfrom, xto);
   fRange->SetLineWidth(3);
   fRange->SetParameters(hRange->GetMaximum(), expectedRange, 2);
   hRange->Fit("fit_range");//, "M,W,B", "", xfrom, xto);

   Float_t cutoff = fRange->GetParameter(1) - 6*fabs(fRange->GetParameter(2));
   if (cutoff < 0) cutoff = 0;
   Float_t cutoffHigh = fRange->GetParameter(1) + 6*fabs(fRange->GetParameter(2));
   Int_t binLow = hRange->GetXaxis()->FindBin(cutoff);
   Int_t binHigh = hRange->GetXaxis()->FindBin(cutoffHigh);
   Float_t total = hRange->Integral();
   Float_t attenuation = hRange->Integral(0, binLow);
   Float_t totalUnderCurve = 0;

   Float_t sumRangeWeight = 0;
   Float_t sigmaRangeWeight = 0;
   for (Int_t j=binLow; j<=binHigh; j++) {
      sumRangeWeight += hRange->GetXaxis()->GetBinCenter(j) * hRange->GetBinContent(j);
      totalUnderCurve += hRange->GetBinContent(j);
   }
   sumRangeWeight /= totalUnderCurve;

   for (Int_t j=binLow; j<binHigh; j++) {
      sigmaRangeWeight += hRange->GetBinContent(j) * pow(hRange->GetXaxis()->GetBinCenter(j) - sumRangeWeight, 2);
   }

   sigmaRangeWeight /= totalUnderCurve-1;
   sigmaRangeWeight = sqrt(sigmaRangeWeight);

   printf("Expecting mean range = %.2f mm. Searching from %.2f mm to %.2f mm.\n", expectedRange, cutoff, cutoffHigh);

   cout << "Mean = " << fRange->GetParameter(1) << endl;
   cout << "3 sigma = " << cutoff << " to " << cutoffHigh << endl;
   cout << "Number of protons attenuated (more than 4 sigma below) = \033[1m" << 100 * attenuation / total << " %\033[0m.\n";
   printf("Estimated range from histogram weighing = %.3f +- %.3f\n",sumRangeWeight, sigmaRangeWeight);
   printf("Estimated range from Gaussian fitting = %.3f +- %.3f\n", fRange->GetParameter(1), fabs(fRange->GetParameter(2)));
   
   TF1 *fRemainingEnergy = new TF1("fRemainingEnergy", "gaus");
   fRemainingEnergy->SetLineWidth(3);
   hEnergyAtInterface->Fit("fRemainingEnergy", "Q");
   printf("Estimated remaining energy: %.2f MeV.\n", phaseSpaceSpline->Eval(run_degraderThickness),);

   printf("Total range and straggling: %.2f +- %.2f mm.\n", run_degraderThickness + aw * pow(fRemainingEnergy->GetParameter(1), pw), sqrt(pow(fRange->GetParameter(2), 2) + pow(fRemainingEnergy->GetParameter(2) * a * p * pow(fRemainingEnergy->GetParameter(1), p-1), 2)));
   printf("Total WEPL straggling: %.2f mm.\n", sqrt(pow(fRemainingEnergy->GetParameter(2) * aw * pw * pow(fRemainingEnergy->GetParameter(1), p-1), 2) + pow(aw/a * pow(fRange->GetParameter(1) / aw, 1-pw/p) * fRange->GetParameter(2), 2)));

//   returnValues.push_back(sumRangeWeight);
//   returnValues.push_back(sigmaRangeWeight);
   returnValues.push_back(fRange->GetParameter(1));
   returnValues.push_back(fabs(fRange->GetParameter(2)));
   returnValues.push_back(100 * attenuation / total);
   returnValues.push_back(fRemainingEnergy->GetParameter(1));
   returnValues.push_back(fRemainingEnergy->GetParameter(2));

/*
   c1->cd();
      hZ->SetXTitle("Z [mm]");
      hZ->SetYTitle("Edep [MeV]");
      hZ->SetFillColor(kBlue-7);
      hZ->SetLineColor(kBlack);
      hZ->Draw();
*/

   c2->cd(1);
      hRange->SetXTitle("Range [mm]");
      hRange->SetYTitle("Number of primaries");
      hRange->SetFillColor(kBlue-7);
      hRange->SetLineColor(kBlack);
      hRange->Draw();
      gPad->Update();
      TLine *l = new TLine(cutoff, 0, cutoff, gPad->GetUymax());
      TLine *l2 = new TLine(cutoffHigh, 0, cutoffHigh, gPad->GetUymax());
      l->SetLineWidth(2);
      l->SetLineStyle(9);
      l->Draw("same");
      l2->SetLineWidth(2);
      l2->SetLineStyle(9);
      l2->Draw("same");


/*
   c7->cd();
      hRangeWEPL->SetXTitle("Range WEPL [mm]");
      hRangeWEPL->SetYTitle("Number of primaries");
      hRangeWEPL->SetFillColor(kBlue-7);
      hRangeWEPL->SetLineColor(kBlack);
      hRangeWEPL->Draw();

   c3->cd();
      hTracklength->SetXTitle("Range [mm]");
      hTracklength->SetYTitle("Number of primaries");
      hTracklength->SetFillColor(kBlue-7);
      hTracklength->SetLineColor(kBlack);
      hTracklength->Draw();
      
   c4->cd();
      hActualTracklength->SetXTitle("Range [mm]");
      hActualTracklength->SetYTitle("Number of primaries");
      hActualTracklength->SetFillColor(kBlue-7);
      hActualTracklength->SetLineColor(kBlack);
      hActualTracklength->Draw();

      
   c5->cd();
      hStepLength->SetXTitle("Steplength [mm]");
      hStepLength->SetYTitle("Number of steps");
      hStepLength->SetFillColor(kBlue-7);
      hStepLength->SetLineColor(kBlack);
      hStepLength->Draw();

   c2->cd(2);
      hEnergyAtInterface->SetFillColor(kBlue-7);
      hEnergyAtInterface->SetLineColor(kBlack);
      hEnergyAtInterface->Draw();
*/

      c2->SaveAs(Form("../OutputFiles/straggling/straggling_%.0f_mm_degrader.png", degraderThickness));
   return returnValues;
}
