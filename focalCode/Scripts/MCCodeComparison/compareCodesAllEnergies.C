#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using namespace std;
enum phantomType {kAluminium, kWater, kComplex};

Float_t findFWHM(TH1F *h) {
   Int_t bin1, bin2, n=0;
   Float_t plateau, fwhm, avg=0;
  
   bin1 = h->FindBin(-2);
   bin2 = h->FindBin(2);
   for (int i=bin1; i<=bin2; i++) {
      avg += h->GetBinContent(i);
      n++;
   }
   avg /= n;

   bin1 = h->FindFirstBinAbove(avg/10.);
   bin2 = h->FindLastBinAbove(avg/10.);
   fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
   printf("FWHM = %.2f\n", fwhm);
 
   return fwhm;
}

void Run()
{
   // RUN SETTINGS
   Bool_t   activateGATE = true;
   Bool_t   activateMCNP = true;
   Bool_t   activateFLUKA = true;
   Bool_t   activateTRIM = false;
   Int_t    phantom = kAluminium;

   Int_t    nbinsxy = 400;
   Int_t    xyfrom = -60;
   Int_t    xyto = 60;
   Int_t    nbinsz = 1000;
   Int_t    zfrom = 0;
   Int_t    zto = 40;
   Int_t    nMax = 500000;
   Float_t  mu, sigma, nominalRange, distributionCutoff, fraction, fwhm;
   Int_t    distributionCutoffBin, totalIntegral, cutoffIntegral, stoppedIntegral, nominalEnergy, bin1, bin2;
   TF1    * fitFunction = nullptr;
      
   Float_t  energies[19]; 
   Float_t  energiesMCNP[19];
   Float_t  energiesFLUKA[19];
   Float_t  energiesGATE[19];
   Float_t  energiesPSTAR[19];
   Float_t  fractionNIMCNP[19];
   Float_t  fractionNIFLUKA[19];
   Float_t  fractionNIGATE[19];
   Float_t  fwhmGATE[19];
   Float_t  fwhmMCNP[19];
   Float_t  fwhmFLUKA[19];
   Float_t  energyError[19] = {};
   Float_t  rangesMCNP[19] = {};
   Float_t  rangesMCNPdiff[19] = {};
   Float_t  sigmaMCNP[19] = {};
   Float_t  rangesFLUKA[19] = {};
   Float_t  rangesFLUKAdiff[19] = {};
   Float_t  sigmaFLUKA[19] = {};
   Float_t  rangesGATE[19] = {};
   Float_t  rangesGATEdiff[19] = {};
   Float_t  sigmaGATE[19] = {};
   Float_t  rangesPSTARWater[19] = {2.224, 3.089, 4.075, 5.176, 6.389, 7.707, 9.128, 10.65, 12.26, 13.96, 15.76, 17.63, 19.59, 21.63, 23.74, 25.93, 28.19, 30.52, 32.91};
   Float_t  rangesPSTARAl[19] = {1.08, 1.5, 1.97, 2.49, 3.07, 3.7, 4.37, 5.09, 5.85, 6.66, 7.51, 8.4, 9.32, 10.28, 11.28, 12.31, 13.38, 14.47, 15.60};
   Float_t  rangesPSTAR[19] = {};
   Float_t  energiesJanni[12] = {50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250};
   Float_t  rangesJanni[12] = {};
   Float_t  stragglingJanni[12] = {}; // %
   Float_t  stragglingJanniError[12] = {}; // % of %
   Float_t  fractionNIJanni[12] = {}; // % 
   Float_t  fractionNIJanniError[12] = {}; // % of %

   Float_t  rangesJanniAl[12] = {1.0733, 1.4848, 1.9527, 2.4744, 3.0475, 3.6698, 5.4282, 7.4545, 9.7252, 12.220, 14.920, 17.811};
   Float_t  stragglingJanniAl[12] = {1.302, 1.276, 1.254, 1.236, 1.220, 1.207, 1.176, 1.153, 1.133, 1.115, 1.100, 1.087};
   Float_t  stragglingJanniErrorAl[12] = {0.019, 0.017, 0.015, 0.014, 0.013, 0.012, 0.01, 0.0088, 0.0078, 0.0071, 0.0066, 0.0061};
   Float_t  fractionNIJanniAl[12] = {4.023, 5.375, 6.817, 8.324, 9.892, 11.51, 15.70, 20.00, 24.36, 28.75, 33.13, 37.46};
   Float_t  fractionNIJanniErrorAl[12] = {0.24, 0.24, 0.24, 0.24, 0.23, 0.23, 0.22, 0.22, 0.21, 0.20, 0.19, 0.18};

   Float_t  rangesJanniWater[12] = {2.252, 3.1259, 4.1226, 5.2360, 6.4610, 7.7929, 11.564, 15.919, 20.807, 26.185, 32.013, 38.257};
   Float_t  stragglingJanniWater[12] = {1.249, 1.227, 1.210, 1.195, 1.181, 1.169, 1.144, 1.123, 1.105, 1.090, 1.076, 1.063};
   Float_t  stragglingJanniErrorWater[12] = {0.017, 0.017, 0.016, 0.015, 0.015, 0.015, 0.014, 0.014, 0.013, 0.013, 0.013, 0.013};
   Float_t  fractionNIJanniWater[12] = {3.487, 4.636, 5.860, 7.136, 8.465, 9.833, 13.29, 16.63, 19.92, 23.30, 26.75, 30.22};
   Float_t  fractionNIJanniErrorWater[12] = {0.24, 0.24, 0.24, 0.24, 0.24, 0.23, 0.23, 0.22, 0.22, 0.21, 0.21, 0.20};

   if (phantom == kAluminium) {
      for (Int_t i=0; i<19; i++) {
         rangesPSTAR[i] = rangesPSTARAl[i];
      }
      
      for (Int_t i=0; i<12; i++) {
         rangesJanni[i] = rangesJanniAl[i];
         stragglingJanni[i] = stragglingJanniAl[i] * rangesJanni[i] / 100.;
         stragglingJanniError[i] = stragglingJanniErrorAl[i];
         fractionNIJanni[i] = fractionNIJanniAl[i] / 100.;
         fractionNIJanniError[i] = fractionNIJanniErrorAl[i] * fractionNIJanni[i];
      }
   }
   else if (phantom == kWater) {
      for (Int_t i=0; i<19; i++) {
         rangesPSTAR[i] = rangesPSTARWater[i];
      }
      for (Int_t i=0; i<12; i++) {
         rangesJanni[i] = rangesJanniWater[i];
         stragglingJanni[i] = stragglingJanniWater[i] * rangesJanni[i] / 100.;
         stragglingJanniError[i] = stragglingJanniErrorWater[i];
         fractionNIJanni[i] = fractionNIJanniWater[i] / 100.;
         fractionNIJanniError[i] = fractionNIJanniErrorWater[i] * fractionNIJanni[i];
      }
   }
   else if (phantom == kComplex) {
      for (Int_t i=0; i<19; i++) {
         rangesPSTAR[i] = 0; // no experimental values for complex geometry
      }
   }

   Float_t  sigmaPSTAR[19] = {};

   Float_t separation = 1.15;

   for (Int_t i=0; i<19; i++) {
      energies[i] = (i+5)*10;
      sigmaPSTAR[i] = rangesPSTAR[i] * pow(energies[i], -0.104) * 0.0188; // FIT TO YANNI
   }

   for (Int_t i=0; i<19; i++) {
      energiesFLUKA[i] = energies[i] + separation*3;
      energiesMCNP[i] = energies[i] + separation;
      energiesGATE[i] = energies[i] - separation;
      energiesPSTAR[i] = energies[i] - separation*3;
   }

   if (phantom == kComplex) {
      for (Int_t i=0; i<19; i++) {
         energiesFLUKA[i] = energies[i] + separation*2;
         energiesMCNP[i] = energies[i];
         energiesGATE[i] = energies[i] - separation*2;
      }
   }
   
   TCanvas *c5 = new TCanvas("c5", "Straggling distribution", 1200, 1000);
   TCanvas *c4 = new TCanvas("c4", "Lateral BP distribution", 1200, 1000);
   TCanvas *c3 = new TCanvas("c3", "Fraction of nuclear interactions", 1200, 1000);
   TCanvas *c2 = new TCanvas("c2", "test", 800, 800);
   TCanvas *c1 = new TCanvas("c1", "Ranges for all codes", 1200, 1000);
   TPad *pad1 = nullptr;
   TPad *pad2 = nullptr;

   if (phantom != kComplex) {
      pad1 = new TPad("pad1", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
      pad2 = new TPad("pad2", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   }
   else {
      pad1 = new TPad("pad1", "The pad 80% of the height", 0.0, 0.0, 1.0, 1.0, 0);
      pad2 = new TPad("pad2", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.0, 0);
   }

   pad1->Draw();
   pad2->Draw();

   gPad->SetFillColor(50);
   gPad->Modified();

   if (activateGATE) {
      Float_t  x,y,z,edep;
      Int_t    parentID, eventID;
      Bool_t   isInelastic;
      TFile   *f1 = nullptr;
      Float_t    nTotal = 0, nInelastic = 0;
     
      for (Int_t i=0; i<19; i++) {
         printf("i = %d\n", i);
         nominalEnergy = (i+5) * 10;
         if (phantom == kAluminium) {
            nominalRange = 0.0012 * pow(nominalEnergy, 1.7483);
         }
         else if (phantom == kWater) {
            nominalRange = 0.0022 * pow(nominalEnergy, 1.77);
         }
         else if (phantom == kComplex) {
            nominalRange = 0.0013 * pow(nominalEnergy, 1.7447); // GUESS!!!!
         }

         TH1F *hGATE = new TH1F("hGATE", "Proton ranges in single GATE dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange - 2), nominalRange + 2);
         TH1F *hGATEAll = new TH1F("hGATEAll", "Proton ranges in single GATE dataset including nuclear interactions;Range [cm];Number of primaries", 200, 0, nominalRange+2);
         TH1F *hGATELateral = new TH1F("hGATELateral", "Lateral BP distribution;X position [cm];Number of end points", 500, -6,6);

         if (phantom == kAluminium) {
            f1 = new TFile(Form("Data/GATE/Aluminium/compressed_aluminium_%dMeV.root", nominalEnergy));
         }
         else if (phantom == kWater) {
            f1 = new TFile(Form("Data/GATE/Water/compressed_water_%dMeV.root", nominalEnergy));
         }
         else if (phantom == kComplex) {
            f1 = new TFile(Form("Data/GATE/ComplexGeometry/compressed_complex_%dMeV.root", nominalEnergy));
         }

         TTree   *treeBic = (TTree*) f1->Get("treeOut");

         cout << "Setting addresses...\n";
         treeBic->SetBranchAddress("posX",&x);
         treeBic->SetBranchAddress("posY",&y);
         treeBic->SetBranchAddress("posZ",&z);
         treeBic->SetBranchAddress("edep",&edep);
         treeBic->SetBranchAddress("eventID",&eventID);
         treeBic->SetBranchAddress("parentID",&parentID);
         treeBic->SetBranchAddress("isInelastic",&isInelastic);

         for (Int_t j=0, N = treeBic->GetEntries(); j<N; ++j) {
            treeBic->GetEntry(j);
            hGATEAll->Fill(z/10.);
            nTotal++;
            if (!isInelastic) {
               hGATE->Fill(z/10.);
               hGATELateral->Fill(y/10.);
            }
            else nInelastic++;
         }
      
         mu = 0; sigma = 0;
         for (Int_t j=0; j<hGATE->GetNbinsX(); j++) {
            mu += hGATE->GetBinContent(j) * hGATE->GetBinCenter(j);
         }
         
         mu /= hGATE->Integral();
         
         for (Int_t j=0; j<hGATE->GetNbinsX(); j++) {
            sigma += hGATE->GetBinContent(j) * pow(hGATE->GetBinCenter(j) - mu, 2);
         }

         sigma /= hGATE->Integral();
         sigma = sqrt(sigma);

         printf("GATE, estimated parameters are %.2f.\n", mu);

         TF1 *fit = new TF1("fit", "gaus");
         hGATE->Fit("fit", "N,Q,WW");
         
         mu = fit->GetParameter(1);
         sigma = fabs(fit->GetParameter(2));

         printf("GATE, %d MeV. nominal Range is %.1f, fit range is %.1f.\n", nominalEnergy, nominalRange, mu);

         rangesGATE[i] = mu;
         sigmaGATE[i] = sigma;
         
         // Calculate fraction of nuclear interactions
         distributionCutoff = mu - 3*sigma;
         distributionCutoffBin = hGATEAll->GetXaxis()->FindBin(distributionCutoff);
         totalIntegral = hGATEAll->Integral();
         cutoffIntegral = hGATEAll->Integral(0, distributionCutoffBin);
         fraction = cutoffIntegral / float(totalIntegral);
         fractionNIGATE[i] = fraction;

         // Calculate FWHM of lateral distribution
//         fwhmGATE[i] = findFWHM(hGATELateral);
         
         Int_t bin1 = hGATELateral->FindFirstBinAbove(hGATELateral->GetMaximum()/10.);
         Int_t bin2 = hGATELateral->FindLastBinAbove(hGATELateral->GetMaximum()/10.);
         Float_t fwhm = hGATELateral->GetBinCenter(bin2) - hGATELateral->GetBinCenter(bin1);
         fwhmGATE[i] = fwhm;

         delete f1;
         delete hGATEAll;
         delete hGATE;
         if (nominalEnergy == 200) {
            c2->cd();
            hGATELateral->SetLineColor(kBlue);
            hGATELateral->Draw();
         }
         else {
            delete hGATELateral;
         }
      }

   }
   
   if (activateMCNP) {
      //
      // MCNP PART
      // 

      ifstream in, in2, in3;
      Int_t    terminationType;
      Int_t    EOLIdentifier;
      Int_t    historyType;
      Int_t    branchNumber;
      Float_t  startZ = 0;
      string   line;
      Float_t  x, y, z;
      Float_t  u, v, w; // vectors
      Float_t  energy, weight, time;
      Float_t  nFilling = 0;

      Int_t    nTotal = 0, nInelastic = 0;

      cout << "READING MCNP FILES...\n";
      for (Int_t i=0; i<19; i++) {
         nominalEnergy = (i+5) * 10;
         cout << "Nominal energy " << nominalEnergy << endl;
         if (phantom == kAluminium) {
            nominalRange = 0.0012 * pow(nominalEnergy, 1.7483);
         }
         else if (phantom == kWater) {
            nominalRange = 0.0022 * pow(nominalEnergy, 1.77);
         }
         else if (phantom == kComplex) {
            nominalRange = 0.0013 * pow(nominalEnergy, 1.7447); // GUESS!!!!
         }

         cout << nominalEnergy << "... ";
         TH1F * hMCNP = new TH1F("hMCNP", "Proton ranges in single MCNP dataset;Range [cm];Number of primaries", 400, fmax(0, nominalRange - 5), nominalRange + 5);
         TH1F * hMCNPAll = new TH1F("hMCNPAll", "Proton ranges in single MCNP dataset including nuclear interactions;Range [cm];Number of primaries", 200, 0, nominalRange+2);
         TH1F * hMCNPLateral = new TH1F("hMCNPLateral", "Lateral BP distribution;X position [cm];Number of end points", 300, -6,6);
         
         if (phantom == kAluminium)  {
            in.open(Form("Data/MCNP/Aluminium/%dMeV_Alp", nominalEnergy));
         }
         else if (phantom == kWater) {
            in.open(Form("Data/MCNP/Water/%dMeV_vannp", nominalEnergy));
         }
         else if (phantom == kComplex) {
            in.open(Form("Data/MCNP/ComplexGeometry/%dMeV_sdetp", nominalEnergy));
         }

         while (! in.eof() ) {
            getline(in, line);

            if (line.size() < 11) continue;
            historyType = atoi(line.substr(6,5).c_str());

            if (historyType == 5000 || historyType == 9000) {
            
               EOLIdentifier = atoi(line.substr(15,7).c_str());
               if (EOLIdentifier == 5000) {
                  continue;
               }

               branchNumber = atoi(line.substr(40,3).c_str());
               terminationType = atoi(line.substr(29,2).c_str());

               if (branchNumber == 1) { // primary particle
                  in >> x >> y >> z >> u >> v >> w >> energy >> weight >> time;

                  hMCNPAll->Fill(z - startZ);
                  nTotal++;
                  if (terminationType != 13 && terminationType != 16) { // not a nuclear interaction
                     hMCNP->Fill((z - startZ));
                     hMCNPLateral->Fill(x);
                  }
                  else nInelastic++;
               }
            }
         }

         in.close();
     
         mu = 0; sigma = 0;
         for (Int_t j=0; j<hMCNP->GetNbinsX(); j++) {
            mu += hMCNP->GetBinContent(j) * hMCNP->GetBinCenter(j);
         }
         
         mu /= hMCNP->Integral();
         
         for (Int_t j=0; j<hMCNP->GetNbinsX(); j++) {
            sigma += hMCNP->GetBinContent(j) * pow(hMCNP->GetBinCenter(j) - mu, 2);
         }

         sigma /= hMCNP->Integral();
         sigma = sqrt(sigma);

         // Find mu, sigma parameters from fit
         TF1 *fit = new TF1("fit", "gaus");
         fit->SetParameter(1, mu);
         fit->SetParameter(2, sigma);
         fit->SetParLimits(1, mu-2, mu+2);
         hMCNP->Fit("fit", "N,Q,M,B");
         mu = fit->GetParameter(1);
         sigma = fabs(fit->GetParameter(2));
         rangesMCNP[i] = mu;
         sigmaMCNP[i] = sigma;

         // Find fraction of nuclear interactions
         distributionCutoff = mu - 3*sigma;
         distributionCutoffBin = hMCNPAll->GetXaxis()->FindBin(distributionCutoff);
         totalIntegral = hMCNPAll->Integral();
         cutoffIntegral = hMCNPAll->Integral(0, distributionCutoffBin);
         fraction = cutoffIntegral / float(totalIntegral);
         fractionNIMCNP[i] = fraction;
         
         // Calculate FWHM of lateral distribution
         fwhmMCNP[i] = findFWHM(hMCNPLateral);
         
         if (nominalEnergy == 200) {
            c2->cd();
            hMCNPLateral->SetLineColor(kRed);
            hMCNPLateral->Draw("same");
         }
         else {
            delete hMCNPLateral;
         }
         delete hMCNP;
         delete hMCNPAll;
         delete fit;
      }
   }
   if (activateFLUKA) {
      //
      // FLUKA PART
      // 

      ifstream in, in2, in3;
      Int_t    historyNumber;
      Int_t    lastHistoryNumber = -1;
      Int_t    terminationType;
      Float_t  x, y, z;
      Float_t  mu;
      Float_t  sigma; 
      Int_t    nTotal = 0, nInelastic = 0;

      for (Int_t i=0; i<19; i++) {
         nominalEnergy = (i+5)*10;
         cout << "Nominal energy " << nominalEnergy << endl;
         if (phantom == kAluminium) { 
            nominalRange = 0.0012 * pow(nominalEnergy, 1.7483);
         }
         else if (phantom == kWater) {
            nominalRange = 0.0022 * pow(nominalEnergy, 1.77);
         }
         else if (phantom == kComplex) {
            nominalRange = 0.0013 * pow(nominalEnergy, 1.7447) * 0.97;
         }

         TH1F * hFLUKA = new TH1F("hFLUKA", "All protons in FLUKA dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange-2), nominalRange+2);
         TH1F * hFLUKAAll = new TH1F("hFLUKAAll", "Proton ranges in single FLUKA dataset including nuclear interactions;Range [cm];Number of primaries", 200, 0, nominalRange+2);
         TH1F * hFLUKALateral = new TH1F("hFLUKALateral", "Lateral BP distribution;X position [cm];Number of end points", 300, -6,6);

         if (phantom == kAluminium) { 
            in.open(Form("Data/FLUKA/Aluminium/%dMevProtons.txt", nominalEnergy));
         }
         else if (phantom == kWater) {
            in.open(Form("Data/FLUKA/Water/%dMeVProtons.txt", nominalEnergy));
         }
         else if (phantom == kComplex) {
            in.open(Form("Data/FLUKA/ComplexGeometry/%dMeVProtons.txt", nominalEnergy));
         }

         while (! in.eof() ) {
            in >> historyNumber >> x >> y >> z >> terminationType;

            if (!in.good()) break;

            hFLUKAAll->Fill(z);
            nTotal++;
            if (terminationType != 11) {
               hFLUKA->Fill(z);
               hFLUKALateral->Fill(x);
            }
            else nInelastic++;
         }

         in.close();

         mu = nominalRange;
         sigma = 0.1;

         // Find mu, sigma parameters
         TF1 *fit = new TF1("fit", "gaus");
         hFLUKA->Fit("fit", "N,Q,M");
         mu = fit->GetParameter(1);
         sigma = fabs(fit->GetParameter(2));
         rangesFLUKA[i] = mu;
         sigmaFLUKA[i] = sigma;
        
         // Find fraction of nuclear interactions
         distributionCutoff = mu - 3*sigma;
         distributionCutoffBin = hFLUKAAll->GetXaxis()->FindBin(distributionCutoff);
         totalIntegral = hFLUKAAll->Integral();
         cutoffIntegral = hFLUKAAll->Integral(0, distributionCutoffBin);
         fraction = cutoffIntegral / float(totalIntegral);
         fractionNIFLUKA[i] = fraction;
         
         // Calculate FWHM of lateral distribution
         fwhmFLUKA[i] = findFWHM(hFLUKALateral);

         delete hFLUKA;
         delete hFLUKAAll;
         if (nominalEnergy == 200) {
            c2->cd();
            hFLUKALateral->SetLineColor(kBlack);
            hFLUKALateral->Draw("same");
         }
         else {
            delete hFLUKALateral;
         }
      }
   }

   cout << "Drawing in pad1\n";
   pad1->cd();
   cout << "Opened pad1\n";

   Float_t errorScalingFactor = 1;
   for (int i=0; i<19; i++) {
      sigmaMCNP[i] *= errorScalingFactor;
      sigmaGATE[i] *= errorScalingFactor;
      sigmaPSTAR[i] *= errorScalingFactor;
      sigmaFLUKA[i] *= errorScalingFactor;

      rangesMCNPdiff[i] = (rangesMCNP[i] - rangesPSTAR[i]) * 10;
      rangesGATEdiff[i] = (rangesGATE[i] - rangesPSTAR[i]) * 10;
      rangesFLUKAdiff[i] = (rangesFLUKA[i] - rangesPSTAR[i]) * 10;
   }

   if (phantom != kComplex) {
      Float_t maxErrorGate = -100;
      Float_t maxErrorMCNP = -100;
      Float_t maxErrorFLUKA = -100;
      Float_t errorGate, errorMCNP, errorFLUKA;

      for (Int_t i=0; i<19; i++) {
         errorGate = rangesGATEdiff[i] / rangesPSTAR[i] * 10;
         errorMCNP = rangesMCNPdiff[i] / rangesPSTAR[i] * 10;
         errorFLUKA = rangesFLUKAdiff[i] / rangesPSTAR[i] * 10;
         maxErrorGate = max(fabs(errorGate), maxErrorGate);
         maxErrorMCNP = max(fabs(errorMCNP), maxErrorMCNP);
         maxErrorFLUKA = max(fabs(errorFLUKA), maxErrorFLUKA);

         cout << "At " << energies[i] << " MeV, the difference between PSTAR and MCNP is " << rangesMCNPdiff[i] / rangesPSTAR[i] * 10 << " %.\n";
         cout << "At " << energies[i] << " MeV, the difference between PSTAR and FLUKA is " << rangesFLUKAdiff[i] / rangesPSTAR[i] * 10 << " %.\n";
         cout << "At " << energies[i] << " MeV, the difference between PSTAR and GATE is " << rangesGATEdiff[i] / rangesPSTAR[i] * 10 << " %.\n";
      }

      printf("Max error GATE %.2f %, max error MCNP %.2f %, max error FLUKA %.2f %.\n", maxErrorGate, maxErrorMCNP, maxErrorFLUKA);
   }
   else {
      Float_t avgDeviationMCNP = 0;
      Float_t maxDeviationMCNP = -100;
      Float_t avgDeviationFLUKA = 0;
      Float_t maxDeviationFLUKA = -100;
      Float_t errorMCNP, errorFLUKA;

      for  (Int_t i=0; i<18; i++) {
         errorMCNP = (rangesMCNP[i] - rangesGATE[i] ) / rangesGATE[i] * 100;
         errorFLUKA = (rangesFLUKA[i] - rangesGATE[i] ) / rangesGATE[i] * 100;

         avgDeviationFLUKA += errorFLUKA;
         avgDeviationMCNP += errorMCNP;
         maxDeviationFLUKA = max(maxDeviationFLUKA, fabs(errorFLUKA));
         maxDeviationMCNP = max(maxDeviationMCNP, fabs(errorMCNP));
      }
      avgDeviationFLUKA /= 18;
      avgDeviationMCNP /= 18;

      printf("The average error between FLUKA and GATE is %.2f %. The max error is %.2f %.", avgDeviationFLUKA, maxDeviationFLUKA);
      printf("The average error between MCNP and GATE is %.2f %. The max error is %.2f %.", avgDeviationMCNP, maxDeviationMCNP);
   }

   Int_t numberOfPoints = 19;
   if (phantom == kComplex) numberOfPoints = 18;

   TGraphErrors *gMCNP = new TGraphErrors(numberOfPoints-1, energiesMCNP, rangesMCNP, energyError, sigmaMCNP);
   if (phantom == kWater) {
      gMCNP->SetTitle(Form("Water range comparison between different codes; Energy [MeV];Range [cm]", errorScalingFactor));
   }
   else if (phantom == kAluminium) {
      gMCNP->SetTitle(Form("Aluminium range comparison between different codes; Energy [MeV];Range [cm]", errorScalingFactor));
   }
   else if (phantom == kComplex) {
      gMCNP->SetTitle(Form("Complex detector geometry range comparison between different codes; Energy [MeV];Range [cm]", errorScalingFactor));
   }

   gMCNP->SetMarkerColor(kRed);
   gMCNP->SetMarkerStyle(7);
   gMCNP->GetXaxis()->SetNdivisions(30);
   gMCNP->Draw("AP");
   

   TGraphErrors *gGATE = new TGraphErrors(numberOfPoints, energiesGATE, rangesGATE, energyError, sigmaGATE);
   gGATE->SetMarkerColor(kBlue);
   gGATE->SetMarkerStyle(7);
   gGATE->Draw("same, P");
   
   TGraphErrors *gPSTAR = new TGraphErrors(numberOfPoints, energiesPSTAR, rangesPSTAR, energyError, sigmaPSTAR);
   gPSTAR->SetMarkerColor(kBlack);
   gPSTAR->SetMarkerStyle(7);
   if (phantom != kComplex) {
      gPSTAR->Draw("same, P");
   }

   TGraphErrors *gFLUKA = new TGraphErrors(numberOfPoints, energiesFLUKA, rangesFLUKA, energyError, sigmaFLUKA);
   gFLUKA->SetMarkerColor(kBlack);
   gFLUKA->SetMarkerStyle(7);
   gFLUKA->Draw("same, P");
   
   pad1->Update();
   for (Int_t i=0; i<20; i++) {
      nominalEnergy = (i+5)*10 - 5;
      TLine *l = new TLine(nominalEnergy, pad1->GetUymin(), nominalEnergy, pad1->GetUymax());
      l->Draw();
   }

   TLegend *leg = new TLegend(0.17, 0.72, 0.35, 0.88);
   if (phantom != kComplex) leg->AddEntry(gPSTAR, "PSTAR", "Ple");
   leg->AddEntry(gGATE, "GATE", "Pel");
   leg->AddEntry(gMCNP, "MCNP", "Pel");
   leg->AddEntry(gFLUKA, "FLUKA", "Pel");
   leg->Draw();

   if (phantom != kComplex) {
      pad2->cd();
      cout << "Drawing in pad2\n";

      TGraph *gMCNPdiff = new TGraph(numberOfPoints-1, energiesMCNP, rangesMCNPdiff);
      gMCNPdiff->SetTitle("; Energy [MeV];Range error [mm]");
      gMCNPdiff->SetLineColor(kRed);
      gMCNPdiff->SetLineWidth(3);
      gMCNPdiff->GetXaxis()->SetTitleSize(0.09);
      gMCNPdiff->GetYaxis()->SetTitleSize(0.09);
      gMCNPdiff->GetYaxis()->SetTitleOffset(0.35);
      gMCNPdiff->GetXaxis()->SetTitleOffset(0.8);
      gMCNPdiff->GetXaxis()->SetLabelSize(0.08);
      gMCNPdiff->GetYaxis()->SetLabelSize(0.08);
      if (phantom == kComplex) gMCNPdiff->GetYaxis()->SetRangeUser(-0.5, 1);
      else if (phantom == kAluminium) gMCNPdiff->GetYaxis()->SetRangeUser(-0.5, 0.8);
      gMCNPdiff->GetXaxis()->SetNdivisions(30);
      gMCNPdiff->Draw("AL");
      
      TGraph *gGATEdiff = new TGraph(numberOfPoints, energiesGATE, rangesGATEdiff);
      gGATEdiff->SetLineColor(kBlue);
      gGATEdiff->SetLineWidth(3);
      gGATEdiff->Draw("same, L");

      TGraph *gFLUKAdiff = new TGraph(numberOfPoints, energiesFLUKA, rangesFLUKAdiff);
      gFLUKAdiff->SetLineColor(kBlack);
      gFLUKAdiff->SetLineWidth(3);
      gFLUKAdiff->Draw("same,L");

      pad2->Update();

      TLine *l = new TLine(pad2->GetUxmin(), 0, pad2->GetUxmax(), 0);
      l->Draw();
   }

   // Draw fraction of nuclear interactions
   c3->cd();
   Float_t dummyArray[12] = {};

   TGraph         *gGATEFractionNI  = new TGraph(numberOfPoints, rangesGATE, fractionNIGATE);
   TGraph         *gMCNPFractionNI  = new TGraph(numberOfPoints-1, rangesMCNP, fractionNIMCNP);
   TGraph         *gFLUKAFractionNI = new TGraph(numberOfPoints, rangesFLUKA, fractionNIFLUKA);
   TGraphErrors   *gJanniFractionNI = new TGraphErrors(12, rangesJanni, fractionNIJanni, dummyArray, fractionNIJanniError);

   gGATEFractionNI->SetTitle("Fraction of nuclear interactions;Range [cm];Fraction of Nuclear Interactions (endpoints below 3#sigma)");
   gGATEFractionNI->SetLineColor(kBlue);
   gGATEFractionNI->SetLineWidth(3);
   gFLUKAFractionNI->SetLineColor(kGreen);
   gFLUKAFractionNI->SetLineWidth(3);
   gMCNPFractionNI->SetLineColor(kRed);
   gMCNPFractionNI->SetLineWidth(3);
   gJanniFractionNI->SetMarkerColor(kBlack);
   gJanniFractionNI->SetMarkerStyle(7);

   gGATEFractionNI->Draw("AL");
   gFLUKAFractionNI->Draw("same, L");
   gMCNPFractionNI->Draw("same, L");
   gJanniFractionNI->Draw("same, P");

   TLegend *legFraction = new TLegend(0.17, 0.72, 0.35, 0.88);
   legFraction->AddEntry(gGATEFractionNI, "GATE", "l");
   legFraction->AddEntry(gMCNPFractionNI, "MCNP", "l");
   legFraction->AddEntry(gFLUKAFractionNI, "FLUKA", "l");
   if (phantom != kComplex) legFraction->AddEntry(gJanniFractionNI, "PSTAR", "Ple");
   legFraction->Draw();

   // Draw lateral BP distribution FWHM
   c4->cd();
   TGraph *gGATELateralFWHM = new TGraph(numberOfPoints, rangesGATE, fwhmGATE);
   TGraph *gFLUKALateralFWHM = new TGraph(numberOfPoints, rangesFLUKA, fwhmFLUKA);
   TGraph *gMCNPLateralFWHM = new TGraph(numberOfPoints-1, rangesMCNP, fwhmMCNP);

   gGATELateralFWHM->SetTitle("Lateral Bragg Peak distribution FWHM;Range [cm];Lateral Bragg Peak size [FWHM cm]");
   gGATELateralFWHM->SetLineColor(kBlue);
   gGATELateralFWHM->SetLineWidth(3);
   gFLUKALateralFWHM->SetLineColor(kGreen);
   gFLUKALateralFWHM->SetLineWidth(3);
   gMCNPLateralFWHM->SetLineColor(kRed);
   gMCNPLateralFWHM->SetLineWidth(3);

   gGATELateralFWHM->Draw("AL");
   gFLUKALateralFWHM->Draw("same, L");
   gMCNPLateralFWHM->Draw("same, L");

   TLegend *legFWHM = new TLegend(0.17, 0.72, 0.35, 0.88);
   legFWHM->AddEntry(gGATELateralFWHM, "GATE", "l");
   legFWHM->AddEntry(gMCNPLateralFWHM, "MCNP", "l");
   legFWHM->AddEntry(gFLUKALateralFWHM, "FLUKA", "l");
   legFWHM->Draw();

   for (int i=0; i<19; i++) {
      printf("Lateral FWHM @ %.0f MeV: %.2f cm (FLUKA) - %.2f cm (GATE) - %.2f cm (MCNP).\n", energies[i], fwhmFLUKA[i], fwhmGATE[i], fwhmMCNP[i]);
   }

   // Draw straggling distribution
   c5->cd();
   TGraph *gGATEStraggling = new TGraph(numberOfPoints, rangesGATE, sigmaGATE);
   TGraph *gMCNPStraggling = new TGraph(numberOfPoints-1, rangesMCNP, sigmaMCNP);
   TGraph *gFLUKAStraggling = new TGraph(numberOfPoints, rangesFLUKA, sigmaFLUKA);
   TGraphErrors *gJanniStraggling = new TGraphErrors(12, rangesJanni, stragglingJanni, dummyArray, stragglingJanniError);

   gGATEStraggling->SetTitle("Range straggling for different ranges;Range [cm];Range straggling [cm]");
   gGATEStraggling->SetLineColor(kBlue);
   gGATEStraggling->SetLineWidth(3);
   gFLUKAStraggling->SetLineColor(kGreen);
   gFLUKAStraggling->SetLineWidth(3);
   gMCNPStraggling->SetLineColor(kRed);
   gMCNPStraggling->SetLineWidth(3);
   gJanniStraggling->SetMarkerColor(kBlack);
   gJanniStraggling->SetMarkerStyle(21);

   gGATEStraggling->Draw("AL");
   gFLUKAStraggling->Draw("same, L");
   gMCNPStraggling->Draw("same, L");
   gJanniStraggling->Draw("same, P");

   TLegend *legStraggling = new TLegend(0.17, 0.72, 0.35, 0.88);
   legStraggling->AddEntry(gGATEStraggling, "GATE", "l");
   legStraggling->AddEntry(gMCNPStraggling, "MCNP", "l");
   legStraggling->AddEntry(gFLUKAStraggling, "FLUKA", "l");
   legStraggling->AddEntry(gJanniStraggling, "Janni", "PLE");
   legStraggling->Draw();

}
