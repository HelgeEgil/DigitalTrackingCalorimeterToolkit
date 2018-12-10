// #include "../GlobalConstants/Constants.h"
#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TTree.h>
#include <vector>
#include <fstream>
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

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement = 1, Int_t mmTo = -1);
vector<Float_t> findRange(Float_t energy, Float_t mm, Float_t degrader, TSpline3 *);

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement, Int_t mmTo) {
   if (mmTo < 0) mmTo = mmFrom;

   // Load phase space spline
   Double_t phaseSpaceDegraderthickness[500];
   Double_t phaseSpaceEnergy[500];
   Double_t e, es;
   Int_t dt;
   Int_t idx = 0;
   ifstream in;
   in.open("../Data/Ranges/EnergyAfterDegrader230MeV.csv");
   while (1) {
      in >> dt >> e;
      if (!in.good()) break;
      phaseSpaceDegraderthickness[idx] = dt;
      phaseSpaceEnergy[idx++] = e;
   }
   in.close();

   TSpline3 *phaseSpaceSpline = new TSpline3("phaseSpaceSpline", phaseSpaceDegraderthickness, phaseSpaceEnergy, idx);

   vector<Float_t> resultVector;
   ofstream filename("../OutputFiles/findManyRangesDegrader_230MeV.csv", ofstream::out | ofstream::app);

   for (Int_t degrader=degraderFrom; degrader<=degraderTo; degrader += degraderIncrement) {
      for (Int_t mm=mmFrom; mm<=mmTo; mm += mmIncrement) {
      	resultVector = findRange(230, mm, degrader, phaseSpaceSpline);
         if (resultVector.size() > 1) {
            Float_t expectedRange = resultVector.at(0);
            Float_t expectedSigma = resultVector.at(1);
            Float_t attenuation   = resultVector.at(2);
            printf("out ... degrader = %d\n", degrader);
            Float_t expectedEnergy = phaseSpaceSpline->Eval(degrader);
            printf("out2 ... degrader = %d\n", degrader);
            Float_t expectedEnergySpread = 6.11e-14*pow(mm,6) - 5.59e-11*pow(mm,5) + 1.90e-8*pow(mm,4) - 2.84e-6*pow(mm,3) + 1.57e-4*pow(mm,2) + 6.88e-3*mm + 2.23e-1; // parameterized from energy spectrum EnergyAfterDegrader.csv, might be a bit wrong but not much
            printf("save\n");
            filename << degrader << " " << mm << " " << expectedRange << " " << expectedSigma << " " << attenuation << " " <<  expectedEnergy << " " << expectedEnergySpread  << endl;
         }
      }
   }
   filename.close();

   delete phaseSpaceSpline;
}


vector<Float_t> findRange(Float_t energy, Float_t mm, Float_t degrader, TSpline3 *phaseSpaceSpline) {
   vector<Float_t> returnValues;
   returnValues.reserve(100);
   Bool_t useDegrader = (degrader > 0) ? true : false;
   float run_energy = phaseSpaceSpline->Eval(degrader);

   TFile *f = new TFile(Form("../Data/MonteCarlo/DTC_Full_Aluminium_Absorber%.0fmm_Degrader%03.0fmm_%dMeV.root", mm, degrader, 230));
   TTree *tree = (TTree*) f->Get("Hits");
   if (!tree) { cout << "Could not find file...\n"; return returnValues; }
 
   Long64_t nentries = tree->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   Float_t posX, posY, posZ, edep;
   Int_t baseID, eventID, parentID;
   tree->SetBranchAddress("posX", &posX);
   tree->SetBranchAddress("posY", &posY);
   tree->SetBranchAddress("posZ", &posZ);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("baseID", &baseID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("eventID", &eventID);

   TCanvas *c2 = new TCanvas("c2", "Ranges and energies", 500, 500);

   Int_t nbinsx = 500;
   // 2 mm: 0.0096, 1.784
   // 3 mm: 0.0097, 1.7825
   // 4 mm: 0.0098, 1.7806
   // Focal: 0.0004461, 1.6677
   // H20:  0.0239, 1.7548

   Float_t  a = 0.011, p = 1.7806;
   Float_t aw = 0.0239, pw = 1.7548;

   Float_t expectedRange = a * pow(run_energy, p);
   Float_t xfrom = expectedRange - 15; // 15 for Al case

   printf("Expected range from aEp = %.2f mm.\n", expectedRange);

   if (xfrom < 0) xfrom = 0;
   Float_t xto = expectedRange + 15; // 15 for Al case

   Float_t x_compensate = 0;

   Int_t energyFrom = run_energy - 15;
   Int_t energyTo = run_energy + 15;
   if (run_energy > 70) {
      energyFrom = run_energy - 15;
      energyTo = run_energy + 25;
   }

   printf("RUNNING WITH ENERGY %.2f.\n", run_energy);

   TH1F *hRange = new TH1F("hRange", "Projected range in DTC", nbinsx, xfrom + x_compensate, xto + x_compensate);
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
   
   TRandom3 *gRandom = new TRandom3();

   tree->GetEntry(0);
   firstX = posX;
   firstY = posY;
   firstZ = posZ;

   Int_t lastP = 0;
   
   if (nentries == 0) return returnValues;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = tree->GetEntry(jentry);
      
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }

      if (lastID < 0) lastID = eventID;

      if (parentID == 0) {
      
         Float_t z = posZ;
         Float_t y = posY;
         Float_t x = posX;
      
         n++;

         if (useDegrader) {
            if (baseID == 0) { // inside degrader
               dE += edep;
            }

            else if (dE > 0) {
               thisEnergy = energy - dE;
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
            Float_t dtcRange = degrader - thisRange + weplRange;

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

      }
   }
   printf("Found %d proton histories in file.\n", lastID);
   
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
//   printf("Estimated remaining energy: %.2f MeV.\n", phaseSpaceSpline->Eval(run_degraderThickness),);

   printf("Total range and straggling: %.2f +- %.2f mm.\n", degrader + aw * pow(fRemainingEnergy->GetParameter(1), pw), sqrt(pow(fRange->GetParameter(2), 2) + pow(fRemainingEnergy->GetParameter(2) * a * p * pow(fRemainingEnergy->GetParameter(1), p-1), 2)));
   printf("Total WEPL straggling: %.2f mm.\n", sqrt(pow(fRemainingEnergy->GetParameter(2) * aw * pw * pow(fRemainingEnergy->GetParameter(1), p-1), 2) + pow(aw/a * pow(fRange->GetParameter(1) / aw, 1-pw/p) * fRange->GetParameter(2), 2)));

   Float_t fR = fRange->GetParameter(1);
   Float_t fRS = fRange->GetParameter(2);

   returnValues.push_back(fR);
   returnValues.push_back(fRS);
   returnValues.push_back(100 * attenuation / total);
   returnValues.push_back(fRemainingEnergy->GetParameter(1));
   returnValues.push_back(fRemainingEnergy->GetParameter(2));

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

   c2->SaveAs(Form("../OutputFiles/straggling/straggling_%.0f_mm_degrader.png", degrader));

   delete c2;
   delete hRange;
   delete fRemainingEnergy;
   delete hEnergyAtInterface;

   delete f;

   return returnValues;
}
