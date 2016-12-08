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
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

Float_t p = 1.7813;
Float_t alpha = 0.02073;
Float_t alphaprime = 0.0087;

Double_t a0 = -3.7279603175;
Double_t a1 =  0.1958689129;
Double_t a2 =  0.0066110840;
Double_t a3 = -0.0000051372;

Float_t getTLFromEnergyQuadratic(Float_t energy) {
   return a0 + a1 * energy + a2 * pow(energy,2) + a3 * pow(energy,3);
}

Double_t getEnergyFromTLQuadratic(Float_t tl) {
   Double_t a = a3;
   Double_t b = a2;
   Double_t c = a1;
   Double_t d = a0;

   Double_t qq = (2*pow(b,3) - 9*a*b*c + 27*pow(a,2) * (d - tl)) / (27*pow(a,3));
   Double_t pp = (3*a*c - pow(b,2)) / (3 * pow(a,2));

   Int_t nRoot = 1;
   Double_t t = 2 * sqrt(-pp/3.) * cos(1/3. * acos((3*qq)/(2*pp) * sqrt(-3/pp)) - nRoot * 2 * 3.14159265/3.);
   Double_t E = t - b/(3*a);

   return E;
}

vector<Float_t> findRange::Run(Double_t energy, Double_t sigma_mev)
{
   vector<Float_t> returnValues;
   if (fChain == 0) return returnValues;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TCanvas *c1 = new TCanvas("c1", "hZ", 800, 600);
   TCanvas *c2 = new TCanvas("c2", "hRange", 800, 600);
   TCanvas *c3 = new TCanvas("c3", "hTracklength", 800, 600);
   TCanvas *c4 = new TCanvas("c4", "hActualTracklength", 800, 600);
   TCanvas *c5 = new TCanvas("c5", "hSteplength", 800, 600);

   Int_t nbinsx = 750;
   // 2 mm: 0.0096, 1.784
   // 3 mm: 0.0097, 1.7825
   // 4 mm: 0.0098, 1.7806
   // H20:  0.0239, 1.7548 
   Float_t expectedRange = 0.0096 * pow(run_energy, 1.784);
   Float_t xfrom = expectedRange * 0.5;
   Float_t xto = expectedRange * 2;

   Float_t x_compensate = 0;

   printf("RUNNING WITH ENERGY %d.\n", run_energy);

   TH1F *hZ = new TH1F("hZ", "Z profile", nbinsx/3, xfrom + x_compensate, xto + x_compensate);
   TH1F *hRange = new TH1F("hRange", "Primary ranges", nbinsx, xfrom + x_compensate, xto + x_compensate);
   TH1F *hTracklength = new TH1F("hTracklength", "Straight tracklengths", nbinsx, xfrom + x_compensate, xto + x_compensate);
   TH1F *hActualTracklength = new TH1F("hActualTracklength", "Actual tracklengths", nbinsx, xfrom + x_compensate, xto + x_compensate);
   TH1F *hStepLength = new TH1F("hStepLength", "Steplenghths", 1000, 0, 1);

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
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   Long64_t ientry = LoadTree(jentry);
   
   if (ientry < 0) {
      cout << "Aborting run at jentry = " << jentry << endl;
      break;
   }

   nb = fChain->GetEntry(jentry);   nbytes += nb;


   if (lastID < 0) lastID = eventID;

   if (parentID == 0) {
   
      hStepLength->Fill(stepLength);

      Float_t z = posZ;
      Float_t y = posY;
      Float_t x = posX;
   
      hZ->Fill(z + x_compensate);
      n++;

      if (processName[0] == 'P') {
         hTracklength->Fill(z + firstZ);
      }

      if (eventID != lastID) {
         n = 0;
         
         Float_t diff = sqrt( pow(firstX - lastX, 2) + pow(firstY - lastY, 2) + pow(0 - lastZ, 2));

         hRange->Fill(lastRange);
         hActualTracklength->Fill(tl + firstZ);
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
      for (Int_t j=0; j<17; j++) lastProcessName[j] = processName[j];
      }
   }
   
   c2->cd();
   
   TF1 *fRange = new TF1("fit_range", "gaus", xfrom, xto);
   fRange->SetParameters(100, expectedRange, 4);
   hRange->Fit("fit_range", "Q,M,W,B", "", xfrom, xto);

   Float_t cutoff = fRange->GetParameter(1) - 6*fabs(fRange->GetParameter(2));
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

   cout << "Mean = " << fRange->GetParameter(1) << endl;
   cout << "3 sigma = " << cutoff << " to " << cutoffHigh << endl;
   cout << "Number of protons attenuated (more than 4 sigma below) = \033[1m" << 100 * attenuation / total << " %\033[0m.\n";
   printf("Estimated range from histogram weighing = %.3f +- %.3f\n",sumRangeWeight, sigmaRangeWeight);

   returnValues.push_back(sumRangeWeight);
   returnValues.push_back(sigmaRangeWeight);
   returnValues.push_back(100 * attenuation / total);

   c1->cd();
      hZ->SetXTitle("Z [mm]");
      hZ->SetYTitle("Edep [MeV]");
      hZ->SetFillColor(kBlue-7);
      hZ->SetLineColor(kBlack);
      hZ->Draw();

   c2->cd();
      hRange->SetXTitle("Range [mm]");
      hRange->SetYTitle("Number of primaries");
      hRange->SetFillColor(kBlue-7);
      hRange->SetLineColor(kBlack);
      hRange->Draw();
      TLine *l = new TLine(cutoff, 0, cutoff, hActualTracklength->GetMaximum()*1.95);
      TLine *l2 = new TLine(cutoffHigh, 0, cutoffHigh, hActualTracklength->GetMaximum()*1.95);
      l->SetLineWidth(2);
      l->SetLineStyle(9);
      l->Draw("same");
      l2->SetLineWidth(2);
      l2->SetLineStyle(9);
      l2->Draw("same");
      
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
   
   return returnValues;
}
