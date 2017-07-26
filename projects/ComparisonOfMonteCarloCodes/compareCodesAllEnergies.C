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
#include <TLatex.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using namespace std;
enum phantomType {kAluminium, kWater, kComplex};

Float_t getBeamSpreading(Float_t initialEnergy, Int_t phantom) {   
   Float_t  X0, angle;
   Float_t  gamma, proton_mass, beta, momentum, mcs, p, a, pv, energy;
   proton_mass = 938.27;

   if       (phantom == kAluminium) {
      X0 = 8.897; // cm
      p = 1.7806;
      a = 0.00098;
   }
   else if  (phantom == kComplex) {
      X0 = 7.241; // cm
      p = 1.7806;
      a = 0.00098;
   }
   else if  (phantom == kWater) {
      X0 = 36.08; // cm
      p = 1.7548;
      a = 0.00239;
   }
   else {
      printf("X0 undefined! Using water values.\n");
      X0 = 36.08;
      p = 1.7548;
      a = 0.00239;
   }
   
   Float_t range = a * pow(initialEnergy, p);
   Float_t halfRange = range / 2;
   Float_t halfRangeEnergy = pow(halfRange / a, 1/p);

   gamma = (halfRangeEnergy + proton_mass) / proton_mass;
   beta = sqrt(1 - pow(gamma, -2));
   momentum = gamma * beta * proton_mass;
   pv = beta * momentum;

   Float_t beamSpreading = 15 / (2 * pv) * sqrt(range/X0);
   return beamSpreading;
}

Float_t getMCSCAngle(Float_t initialEnergy, Float_t range, Int_t phantom) {
   Float_t  X0, angle;
   Float_t  gamma, proton_mass, beta, momentum, mcs, p, a, pv, energy;
   proton_mass = 938.27;

   if       (phantom == kAluminium) {
      X0 = 8.897; // cm
      p = 1.7806;
      a = 0.00098;
   }
   else if  (phantom == kComplex) {
      X0 = 7.241; // cm
      p = 1.7806;
      a = 0.00098;
   }
   else if  (phantom == kWater) {
      X0 = 36.08; // cm
      p = 1.7548;
      a = 0.00239;
   }
   else {
      printf("X0 undefined! Using water values.\n");
      X0 = 36.08;
      p = 1.7548;
      a = 0.00239;
   }

   Float_t sumIntegral = 0, depth = 0;
   Int_t nTerms = 50000;
   Float_t binThickness = range / float(nTerms);
   for (Int_t i=1; i<nTerms; i++) {
      // find energy
      // R = a E^p
      // R - z = a (E)^p
      // E^p = (R - z) / a
      // E = [(R - z) / a] ^ 1/p
      // z = 0 -> E = E0
      // z = R -> E = 0
      
      depth = i * range / nTerms;
      energy = pow((range - depth) / a, 1/p);

      gamma = (initialEnergy + proton_mass) / proton_mass;
      beta = sqrt(1 - pow(gamma, -2));
      momentum = gamma * beta * proton_mass;
      pv = beta * momentum;

      sumIntegral += pow(1/pv, 2) * binThickness / X0;
   }

   mcs = 19.2 * ( 1 + 0.038 * log(range / X0)) * sqrt(sumIntegral);
   mcs *= 180 / 3.14159265358979;
      
   gamma = (initialEnergy + proton_mass) / proton_mass;
   beta = sqrt(1 - pow(gamma, -2));
   momentum = gamma * beta * proton_mass;
   pv = beta * momentum;

   Float_t mcsSimple = 19.2 / pv * sqrt(range/X0) * (1 + 0.038 * log(range / X0)) * 180 / 3.14159265;
   printf("MCS with a simple calculation is %.2f deg. MCS with the integral calculation is %.2f deg.\n", mcsSimple, mcs);

   return mcs;
}


Float_t findFWHM(TH1F *h) {
   Int_t bin1, bin2, n=0;
   Float_t fwhm, maximum;

   maximum = h->GetMaximum();

   bin1 = h->FindFirstBinAbove(maximum/2.);
   bin2 = h->FindLastBinAbove(maximum/2.);
   fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
 
   return fwhm;
}

void Run()
{
   // RUN SETTINGS
   Bool_t   activateGATE = true;
   Bool_t   activateMCNP = true;
   Bool_t   activateFLUKA = true;
   Bool_t   activateTRIM = false;
   Int_t    phantom = kComplex;

   gStyle->SetLabelSize(0.04);
   gStyle->SetLabelSize(0.04, "Y");
   gStyle->SetTitleSize(0.04);
   gStyle->SetTitleSize(0.04, "Y");
   gStyle->SetTitleOffset(1.1);
//   gStyle->SetTitleOffset(1.35, "Y");
   gStyle->SetOptStat(0);
   
   Float_t  X0, angle;
   Float_t  gamma, proton_mass, beta, momentum, mcs, p, a, ap, pv, energy;
   proton_mass = 938.27;

   if       (phantom == kAluminium) {
      X0 = 8.897; // cm
      p = 1.7806;
      a = 0.00098;
      ap = 0.087 * 2.36; // From GAMMEX manual

   }
   else if  (phantom == kComplex) {
      X0 = 7.241; // cm
      p = 1.7806;
      a = 0.00098;
      ap = 0.087 * 2.27; // Calculated from GAMMEX manual parameterization applied on DTC materials, and then summed each material's contribution to the fractional ED relative to water... Maybe not correct, but as the Catphan manufacturer is named: D. Goodenough.
   }
   else if  (phantom == kWater) {
      X0 = 36.08; // cm
      p = 1.7548;
      a = 0.00239;
      ap = 0.087;
   }
   else {
      printf("X0 undefined! Using water values.\n");
      X0 = 36.08;
      p = 1.7548;
      a = 0.00239;
   }

   Int_t    nbinsxy = 400;
   Int_t    xyfrom = -60;
   Int_t    xyto = 60;
   Int_t    nbinsz = 1000;
   Int_t    zfrom = 0;
   Int_t    zto = 40;
   Int_t    nMax = 500000;
   Int_t    lateralfrom = -2;
   Int_t    lateralto = 2;
   Int_t    comparisonEnergy = 120;
   Float_t  mu, sigma, nominalRange, distributionCutoff, fraction, fwhm, mutl;
   Int_t    distributionCutoffBin, totalIntegral, cutoffIntegral, stoppedIntegral, nominalEnergy, bin1, bin2;
   TF1    * fitFunction = nullptr;

   Int_t    flukaColor = kGreen-3;
   Int_t    gateColor = kBlue-7;
   Int_t    expColor = kBlack;
   Int_t    mcnpColor = kRed;

   Int_t    flukaStyle = 1;
   Int_t    gateStyle = 9;
   Int_t    mcnpStyle = 2;
      
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
   Float_t  beamSpreadingTheory[19] = {};
   Float_t  sigmaFLUKA[19] = {};
   Float_t  rangesGATE[19] = {};
   Float_t  rangesGATE75[19] = {};
   Float_t  rangesGATE74[19] = {};
   Float_t  rangesGATE73[19] = {};
   Float_t  rangesGATE72[19] = {};
   Float_t  rangesGATE71[19] = {};
   Float_t  rangesGATEdiff[19] = {};
   Float_t  rangesGATEdiff75[19] = {};
   Float_t  rangesGATEdiff74[19] = {};
   Float_t  rangesGATEdiff73[19] = {};
   Float_t  rangesGATEdiff72[19] = {};
   Float_t  rangesGATEdiff71[19] = {};
   Float_t  detourGATE[19] = {};
   Float_t  detourFLUKA[19] = {};
   Float_t  detourMCNP[19] = {};
   Float_t  sigmaGATE[19] = {};
   Float_t  rangesPSTARWater[19] = {2.224, 3.089, 4.075, 5.176, 6.389, 7.707, 9.128, 10.65, 12.26, 13.96, 15.76, 17.63, 19.59, 21.63, 23.74, 25.93, 28.19, 30.52, 32.91};
   Float_t  rangesPSTARAl[19] = {1.08, 1.5, 1.97, 2.49, 3.07, 3.7, 4.37, 5.09, 5.85, 6.66, 7.51, 8.4, 9.32, 10.28, 11.28, 12.31, 13.38, 14.47, 15.60};
   Float_t  rangesPSTAR[19] = {};
   Float_t  energiesJanni[12] = {50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 225, 250};
   Float_t  rangesJanni[12] = {};
   Float_t  stragglingJanni[25] = {}; // %
   Float_t  stragglingJanniError[25] = {}; // % of %
   Float_t  fractionNIJanni[25] = {}; // % 
   Float_t  fractionNIJanniError[25] = {}; // % of %

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
  
   Float_t  detourPSTARWater[19] = {0.9985, 0.9986, 0.9986, 0.9986, 0.9986, 0.9987, 0.9987, 0.9987, 0.9987, 0.9987, 0.9987, 0.9988, 0.9988, 0.9988, 0.9988, 0.9988, 0.9988, 0.9988, 0.9988};
   Float_t  detourPSTARAl[19] =    {0.9967, 0.9968, 0.9969, 0.9970, 0.9970, 0.9971, 0.9971, 0.9972, 0.9972, 0.9972, 0.9973, 0.9973, 0.9973, 0.9974, 0.9974, 0.9974, 0.9974, 0.9975, 0.9975};

   Float_t  highlandMCSAl[19] = {93.96, 93.66, 94.46, 93.33, 93.25, 93.18, 93.16, 93.14, 93.12, 92.12, 93.12, 93.12, 93.13, 93.14, 93.15, 93.17, 93.18, 93.19, 93.21}; // mrad
   Float_t  highlandMCSWater[19] = {63.16, 63.25, 63.34, 63.45, 63.54, 63.64, 63.73, 63.81, 63.89, 63.96, 64.03, 64.10, 64.17, 64.23, 64.29, 64.34, 64.40, 64.45, 64.49};

   Float_t  sigmaPSTAR[19] = {};

   for (Int_t i=0; i<19; i++) {
      energies[i] = (i+5)*10;
      sigmaPSTAR[i] = rangesPSTAR[i] * pow(energies[i], -0.104) * 0.0188; // FIT TO YANNI
   }

   if (phantom == kAluminium) {
      for (Int_t i=0; i<19; i++) {
         rangesPSTAR[i] = rangesPSTARAl[i];
      }
      
      for (Int_t i=0; i<12; i++) {
         rangesJanni[i] = rangesJanniAl[i];
         stragglingJanni[i] = stragglingJanniAl[i] * rangesJanni[i] / 100.*10.;
         stragglingJanniError[i] = stragglingJanniErrorAl[i] * 10.;
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
         stragglingJanni[i] = stragglingJanniWater[i] * rangesJanni[i] / 100. *10.;
         stragglingJanniError[i] = stragglingJanniErrorWater[i] * 10.;
         fractionNIJanni[i] = fractionNIJanniWater[i] / 100.;
         fractionNIJanniError[i] = fractionNIJanniErrorWater[i] * fractionNIJanni[i];
      }
   }
   else if (phantom == kComplex) {
      for (Int_t i=0; i<18; i++) {
         rangesPSTAR[i] = 0; // no experimental values for complex geometry
      }
   }


   Float_t separation = 1.15;


   for (Int_t i=0; i<19; i++) {
      energiesFLUKA[i] = energies[i] + separation*3;
      energiesMCNP[i] = energies[i] + separation;
      energiesGATE[i] = energies[i] - separation;
      energiesPSTAR[i] = energies[i] - separation*3;

      beamSpreadingTheory[i] = getBeamSpreading(energies[i], phantom);
   }

   if (phantom == kComplex) {
      for (Int_t i=0; i<19; i++) {
         energiesFLUKA[i] = energies[i] + separation*2;
         energiesMCNP[i] = energies[i];
         energiesGATE[i] = energies[i] - separation*2;
      }
   }
   
   TCanvas *c7 = new TCanvas("c7", "Range comparisons", 1000, 1000);
   TCanvas *c6 = new TCanvas("c6", "GATE different IP vs MCNP water", 1200, 1000);
   TCanvas *c5 = new TCanvas("c5", "Straggling distribution", 1000, 1000);
   TCanvas *c4 = new TCanvas("c4", "Lateral BP distribution", 1000, 1000);
   TCanvas *c3 = new TCanvas("c3", "Fraction of nuclear interactions", 1000, 1000);
   TCanvas *c2 = new TCanvas("c2", "test", 800, 800);
   TCanvas *c1 = new TCanvas("c1", "Ranges for all codes", 1200, 1000);
   TCanvas *c8 = new TCanvas("c8", "Detour factor", 1200, 1000);
         
   
   TH1F *hGATELateral = new TH1F("hGATELateral", ";Bragg peak displacement in x [mm];Number of end points", 500, lateralfrom, lateralto);
   TH1F * hFLUKALateral = new TH1F("hFLUKALateral", ";Bragg peak displacement in x [mm];Number of end points", 500, lateralfrom,lateralto);
   TH1F * hMCNPLateral = new TH1F("hMCNPLateral", "Angular BP distribution;#theta [mrad];Number of end points", 500, lateralfrom,lateralto);

   c1->cd();
   TPad *pad1 = nullptr;
   TPad *pad2 = nullptr;

   if (true) {
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
         TH1F *hGATETL = new TH1F("hGATE", "Proton tracklengths in single GATE dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange - 2), nominalRange + 2);
         TH1F *hGATEAll = new TH1F("hGATEAll", "Proton ranges in single GATE dataset including nuclear interactions;Range [cm];Number of primaries", 200, 0, nominalRange+2);
         
         TH1F *hGATE74 = new TH1F("hGATE74", "Proton ranges in single GATE dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange - 2), nominalRange + 2);
         TH1F *hGATE73 = new TH1F("hGATE73", "Proton ranges in single GATE dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange - 2), nominalRange + 2);
         TH1F *hGATE72 = new TH1F("hGATE72", "Proton ranges in single GATE dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange - 2), nominalRange + 2);
         TH1F *hGATE71 = new TH1F("hGATE71", "Proton ranges in single GATE dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange - 2), nominalRange + 2);

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

         Float_t lateralDeviation, angle;
         for (Int_t j=0, N = treeBic->GetEntries(); j<N; ++j) {
            treeBic->GetEntry(j);
            hGATEAll->Fill(z/10.);
            nTotal++;
            if (!isInelastic) {
               hGATE->Fill(z/10.);
               hGATETL->Fill(sqrt(pow(z/10., 2) + pow(x/10., 2) + pow(y/10., 2)));
               angle = atan(x/10. / (z/10.)) * 1000;
//               hGATELateral->Fill(angle);
               hGATELateral->Fill(x/10.);
            }
            else nInelastic++;
         }
         delete f1;
      
         for (Int_t ip = 74; ip>70; ip--) {
            f1 = new TFile(Form("Data/GATE/Water/compressed_water_%d_%dMeV.root", ip, nominalEnergy));
            TTree *tree = (TTree*) f1->Get("treeOut");
            tree->SetBranchAddress("posX",&x);
            tree->SetBranchAddress("posY",&y);
            tree->SetBranchAddress("posZ",&z);
            tree->SetBranchAddress("edep",&edep);
            tree->SetBranchAddress("eventID",&eventID);
            tree->SetBranchAddress("parentID",&parentID);
            tree->SetBranchAddress("isInelastic",&isInelastic);
            for (Int_t j=0, N = tree->GetEntries(); j<N; ++j) {
               tree->GetEntry(j);
               if (!isInelastic) {
                  if      (ip == 74)   hGATE74->Fill(z/10.);
                  else if (ip == 73)   hGATE73->Fill(z/10.);
                  else if (ip == 72)   hGATE72->Fill(z/10.);
                  else if (ip == 71)   hGATE71->Fill(z/10.);
               }
            }
         }
       
         printf("Mean values at %d MeV: 75 = %.2f, 74 = %.2f, 73 = %.2f, 72 = %.2f, 71 = %.2f\n", nominalEnergy, hGATE->GetMean(), hGATE74->GetMean(), hGATE73->GetMean(), hGATE72->GetMean(), hGATE71->GetMean());

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
         hGATE->Fit("fit", "N,WW");
         
         mu = fit->GetParameter(1);
         sigma = fabs(fit->GetParameter(2));
         
         hGATETL->Fit("fit", "N,Q,M,B");
         mutl = fit->GetParameter(1);

         detourGATE[i] = mu / mutl;
         printf("Detour factor GATE = %.4f.\n", detourGATE[i]);

         printf("GATE, %d MeV. nominal Range is %.1f, fit range is %.1f.\n", nominalEnergy, nominalRange, mu);

         rangesGATE[i] = mu;
         sigmaGATE[i] = sigma*10;

         // FIND RANGES FOR DIFFERENT IONIZATION POTENTIALS
         hGATE74->Fit("fit", "N,Q,WW");
         rangesGATE74[i] = fit->GetParameter(1);
         hGATE73->Fit("fit", "N,Q,WW");
         rangesGATE73[i] = fit->GetParameter(1);
         hGATE72->Fit("fit", "N,Q,WW");
         rangesGATE72[i] = fit->GetParameter(1);
         hGATE71->Fit("fit", "N,Q,WW");
         rangesGATE71[i] = fit->GetParameter(1);
         
         printf("Different IPs: %.2f (75), %.2f, %.2f, %.2f, %.2f.\n", mu, rangesGATE74[i], rangesGATE73[i], rangesGATE72[i], rangesGATE71[i]);
         // Calculate fraction of nuclear interactions
         distributionCutoff = mu - 3*sigma;
         distributionCutoffBin = hGATEAll->GetXaxis()->FindBin(distributionCutoff);
         totalIntegral = hGATEAll->Integral();
         cutoffIntegral = hGATEAll->Integral(0, distributionCutoffBin);
         fraction = cutoffIntegral / float(totalIntegral);
//         fractionNIGATE[i] = fraction;
         fractionNIGATE[i] = 1 - hGATE->Integral() / hGATEAll->Integral();


         // Calculate FWHM of lateral distribution
//         fwhmGATE[i] = findFWHM(hGATELateral);
         fwhmGATE[i] = hGATELateral->GetRMS() / mu;
//         hGATELateral->Fit("fit", "N,Q,WW");
//         fwhmGATE[i] = fit->GetParameter(2);
         
         delete f1;
         delete hGATEAll;
         delete hGATE;
         if (nominalEnergy == comparisonEnergy) {
            c2->cd();
//            hGATELateral->Scale(1/hGATELateral->Integral());
            hGATELateral->SetLineColor(gateColor);
            hGATELateral->SetLineStyle(gateStyle);
            hGATELateral->SetLineWidth(2);
            hGATELateral->GetXaxis()->SetTitleFont(22);
            hGATELateral->GetXaxis()->SetLabelFont(22);
            hGATELateral->GetYaxis()->SetTitleFont(22);
            hGATELateral->GetYaxis()->SetLabelFont(22);
            hGATELateral->GetYaxis()->SetTitleOffset(1.4);
            gPad->SetLogy();
            hGATELateral->Draw();
         }
         hGATELateral->Reset();
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
         TH1F * hMCNPTL = new TH1F("hMCNP", "Proton tracklengths in single MCNP dataset;Range [cm];Number of primaries", 400, fmax(0, nominalRange - 5), nominalRange + 5);
         TH1F * hMCNPAll = new TH1F("hMCNPAll", "Proton ranges in single MCNP dataset including nuclear interactions;Range [cm];Number of primaries", 200, 0, nominalRange+2);
         
         if (phantom == kAluminium)  {
            in.open(Form("Data/MCNP/Aluminium/%dMeV_Alp", nominalEnergy));
         }
         else if (phantom == kWater) {
            in.open(Form("Data/MCNP/Water/%dMeV_vannp", nominalEnergy));
         }
         else if (phantom == kComplex) {
            in.open(Form("Data/MCNP/ComplexGeometry/%dMeV_sdetp", nominalEnergy));
         }

         Float_t lateralDeviation, angle;
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
                     hMCNPTL->Fill(sqrt(pow(z - startZ, 2) + pow(x, 2) + pow(y, 2)));
                     angle = atan(x / z) * 1000;
//                     hMCNPLateral->Fill(angle);
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
         sigmaMCNP[i] = sigma*10;

         hMCNPTL->Fit("fit", "N,Q,M,B");
         mutl = fit->GetParameter(1);

         detourMCNP[i] = mu / mutl;
         printf("Detour factor MCNP = %.4f.\n", detourMCNP[i]);

         // Find fraction of nuclear interactions
         distributionCutoff = mu - 3*sigma;
         distributionCutoffBin = hMCNPAll->GetXaxis()->FindBin(distributionCutoff);
         totalIntegral = hMCNPAll->Integral();
         cutoffIntegral = hMCNPAll->Integral(0, distributionCutoffBin);
         fraction = cutoffIntegral / float(totalIntegral);
//         fractionNIMCNP[i] = fraction;
         fractionNIMCNP[i] = 1 - hMCNP->Integral() / hMCNPAll->Integral();
         
         // Calculate FWHM of lateral distribution
//         fwhmMCNP[i] = findFWHM(hMCNPLateral);
         fwhmMCNP[i] = hMCNPLateral->GetRMS() / mu;
         
         if (nominalEnergy == comparisonEnergy) {
            c2->cd();
//            hMCNPLateral->Scale(1/hMCNPLateral->Integral());
            hMCNPLateral->SetLineColor(mcnpColor);
            hMCNPLateral->SetLineStyle(mcnpStyle);
            hMCNPLateral->SetLineWidth(2);
            hMCNPLateral->Draw("same");
         }
         hMCNPLateral->Reset();

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
         TH1F * hFLUKATL = new TH1F("hFLUKATL", "All protons tracklengths in FLUKA dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange-2), nominalRange+2);
         TH1F * hFLUKAAll = new TH1F("hFLUKAAll", "Proton ranges in single FLUKA dataset including nuclear interactions;Range [cm];Number of primaries", 200, 0, nominalRange+2);

         if (phantom == kAluminium) { 
            in.open(Form("Data/FLUKA/Aluminium/%dMevProtonsPencil.txt", nominalEnergy));
         }
         else if (phantom == kWater) {
            in.open(Form("Data/FLUKA/Water/%dMeVProtonsPencil.txt", nominalEnergy));
         }
         else if (phantom == kComplex) {
            in.open(Form("Data/FLUKA/ComplexGeometry/%dMeVProtonsPencil.txt", nominalEnergy));
         }

         Float_t lateralDeviation, angle;
         while (! in.eof() ) {
            in >> historyNumber >> x >> y >> z >> terminationType;

            if (!in.good()) break;

            hFLUKAAll->Fill(z);
            nTotal++;
            if (terminationType != 11) {
               hFLUKA->Fill(z);
               hFLUKATL->Fill(sqrt(pow(z, 2) + pow(x,2) + pow(y,2)));
               angle = atan(x / z) * 1000;
              // hFLUKALateral->Fill(angle);
               hFLUKALateral->Fill(x);
            }
            else {
               nInelastic++;
            }
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
         sigmaFLUKA[i] = sigma*10;

         hFLUKATL->Fit("fit", "N,Q,M");
         mutl = fit->GetParameter(1);
         detourFLUKA[i] = mu / mutl;
         printf("FLUKA detour factor is %.4f\n", detourFLUKA[i]);
        
         // Find fraction of nuclear interactions
         distributionCutoff = mu - 3*sigma;
         distributionCutoffBin = hFLUKAAll->GetXaxis()->FindBin(distributionCutoff);
         totalIntegral = hFLUKAAll->Integral();
         cutoffIntegral = hFLUKAAll->Integral(0, distributionCutoffBin);
         fraction = cutoffIntegral / float(totalIntegral);
//         fractionNIFLUKA[i] = fraction;
         fractionNIFLUKA[i] = 1 - hFLUKA->Integral() / hFLUKAAll->Integral();
         
         // Calculate FWHM of lateral distribution
//         fwhmFLUKA[i] = findFWHM(hFLUKALateral);
         fwhmFLUKA[i] = hFLUKALateral->GetRMS() / mu;

         delete hFLUKA;
         delete hFLUKAAll;
         if (nominalEnergy == comparisonEnergy) {
            c2->cd();
//            hFLUKALateral->Scale(1/hFLUKALateral->Integral());
            hFLUKALateral->SetLineColor(flukaColor);
            hFLUKALateral->SetLineStyle(flukaStyle);
            hFLUKALateral->SetLineWidth(2);
            hFLUKALateral->Draw("same");
//            TF1 *mcsGauss = new TF1("mcsGauss", "gaus(0)", 0, 40);
            // energy, depth, phantom
//            mcsGauss->SetParameters(hFLUKALateral->GetMaximum(), 0, getMCSCAngle(nominalEnergy, mu, phantom));
//            mcsGauss->Draw("same");
            TF1 *highland = new TF1("highland", "gaus(0)", -200, 200);
            if (phantom == kWater) {
               highland->SetParameters(hFLUKALateral->GetMaximum(), 0, 63.81);
            }
         
            else if (phantom == kAluminium) {
               highland->SetParameters(hFLUKALateral->GetMaximum(), 0, 93.14);
            }
//            highland->Draw("same");
            printf("FLUKA rms is %.2f mrad\n", hFLUKALateral->GetRMS());
         }
         hFLUKALateral->Reset();
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
      rangesGATEdiff75[i] = (rangesGATE[i] - rangesMCNP[i]) * 10;
      rangesGATEdiff74[i] = (rangesGATE74[i] - rangesMCNP[i]) * 10;
      rangesGATEdiff73[i] = (rangesGATE73[i] - rangesMCNP[i]) * 10;
      rangesGATEdiff72[i] = (rangesGATE72[i] - rangesMCNP[i]) * 10;
      rangesGATEdiff71[i] = (rangesGATE71[i] - rangesMCNP[i]) * 10;
      rangesFLUKAdiff[i] = (rangesFLUKA[i] - rangesPSTAR[i]) * 10;
   }

   if (phantom == kComplex) {
      for (int i=0; i<19; i++) {
         sigmaMCNP[i] *= errorScalingFactor;
         sigmaGATE[i] *= errorScalingFactor;
         sigmaPSTAR[i] *= errorScalingFactor;
         sigmaFLUKA[i] *= errorScalingFactor;

         rangesMCNPdiff[i] = 10 * (rangesMCNP[i] - (rangesMCNP[i] + rangesGATE[i] + rangesFLUKA[i]) / 3);
         rangesGATEdiff[i] = 10 * (rangesGATE[i] - (rangesMCNP[i] + rangesGATE[i] + rangesFLUKA[i]) / 3);
         rangesFLUKAdiff[i] = 10 * (rangesFLUKA[i] - (rangesMCNP[i] + rangesGATE[i] + rangesFLUKA[i]) / 3);
      }
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
   if (phantom == kComplex) numberOfPoints = 17;

   TGraphErrors *gMCNP = new TGraphErrors(numberOfPoints, energiesMCNP, rangesMCNP, energyError, sigmaMCNP);
   if (phantom == kWater) {
      gMCNP->SetTitle(Form("Range comparison between MC codes, in water phantom; Energy bin [MeV];Range [cm]", errorScalingFactor));
   }
   else if (phantom == kAluminium) {
      gMCNP->SetTitle(Form("Range comparison between MC codes, in aluminium phantom; Initial energy [MeV];Range [cm]", errorScalingFactor));
   }
   else if (phantom == kComplex) {
      gMCNP->SetTitle(Form("Range comparison between MC codes, in detector geometry; Initial energy [MeV];Range [cm]", errorScalingFactor));
   }

   gMCNP->SetMarkerColor(mcnpColor);
   gMCNP->SetMarkerStyle(7);
   gMCNP->SetLineWidth(2);
   gMCNP->SetLineColor(mcnpColor);
   gMCNP->SetFillColor(mcnpColor);
   gMCNP->GetXaxis()->SetNdivisions(30);
   gMCNP->Draw("AP");
   gMCNP->GetXaxis()->SetRangeUser(45, 235);

   TGraphErrors *gGATE = new TGraphErrors(numberOfPoints, energiesGATE, rangesGATE, energyError, sigmaGATE);
   gGATE->SetMarkerColor(gateColor);
   gGATE->SetFillColor(gateColor);
   gGATE->SetMarkerStyle(7);
   gGATE->SetLineWidth(2);
   gGATE->SetLineColor(gateColor);
   gGATE->Draw("same, P");
   
   TGraphErrors *gPSTAR = new TGraphErrors(numberOfPoints, energiesPSTAR, rangesPSTAR, energyError, sigmaPSTAR);
   gPSTAR->SetMarkerColor(expColor);
   gPSTAR->SetFillColor(expColor);
   gPSTAR->SetMarkerStyle(7);
   gPSTAR->SetLineWidth(2);
   gPSTAR->SetLineColor(expColor);
   if (phantom != kComplex) {
      gPSTAR->Draw("same, P");
   }

   TGraphErrors *gFLUKA = new TGraphErrors(numberOfPoints, energiesFLUKA, rangesFLUKA, energyError, sigmaFLUKA);
   gFLUKA->SetMarkerColor(flukaColor);
   gFLUKA->SetMarkerStyle(7);
   gFLUKA->SetFillColor(flukaColor);
   gFLUKA->SetLineColor(flukaColor);
   gFLUKA->SetLineWidth(2);
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

   if (true) {
      pad2->cd();
      cout << "Drawing in pad2\n";

      TGraph *gMCNPdiff = new TGraph(numberOfPoints, energies, rangesMCNPdiff);
      gMCNPdiff->SetTitle("; Initial energy [MeV];Range error [mm]");
      gMCNPdiff->SetLineColor(mcnpColor);
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
      
      TGraph *gGATEdiff = new TGraph(numberOfPoints, energies, rangesGATEdiff);
      gGATEdiff->SetLineColor(gateColor);
      gGATEdiff->SetLineWidth(3);
      gGATEdiff->Draw("same, L");

      TGraph *gFLUKAdiff = new TGraph(numberOfPoints, energies, rangesFLUKAdiff);
      gFLUKAdiff->SetLineColor(flukaColor);
      gFLUKAdiff->SetLineWidth(3);
      gFLUKAdiff->Draw("same,L");

      pad2->Update();

      TLine *l = new TLine(pad2->GetUxmin(), 0, pad2->GetUxmax(), 0);
      l->Draw();
   }

   // Test MCNP vs GATE (different IPs)
   c6->cd();
//   TGraph *gMCNP75diff = new TGraph(numberOfPoints, energies, rangesMCNPdiff);
   TGraph *gGATE75diff = new TGraph(numberOfPoints, energies, rangesGATEdiff75);
   TGraph *gGATE74diff = new TGraph(numberOfPoints, energies, rangesGATEdiff74);
   TGraph *gGATE73diff = new TGraph(numberOfPoints, energies, rangesGATEdiff73);
   TGraph *gGATE72diff = new TGraph(numberOfPoints, energies, rangesGATEdiff72);
   TGraph *gGATE71diff = new TGraph(numberOfPoints, energies, rangesGATEdiff71);

//   gMCNP75diff->SetLineColor(mcnpColor);
//   gMCNP75diff->SetLineWidth(3);
//   gMCNP75diff->Draw("LA");
   
   gGATE75diff->SetTitle("; Initial energy [MeV];Range deviation MCNP - GATE [mm]");
   gGATE75diff->SetLineColor(kBlue-10);
   gGATE74diff->SetLineColor(kBlue-8);
   gGATE73diff->SetLineColor(kBlue-5);
   gGATE72diff->SetLineColor(kBlue-1);
   gGATE71diff->SetLineColor(kBlue+4);
   gGATE75diff->SetLineWidth(3);
   gGATE74diff->SetLineWidth(3);
   gGATE73diff->SetLineWidth(3);
   gGATE72diff->SetLineWidth(3);
   gGATE71diff->SetLineWidth(3);

   gGATE75diff->GetXaxis()->SetTitleFont(22);
   gGATE75diff->GetYaxis()->SetTitleFont(22);
   gGATE75diff->GetXaxis()->SetLabelFont(22);
   gGATE75diff->GetYaxis()->SetLabelFont(22);
   gGATE75diff->GetYaxis()->SetRangeUser(-1, 1.5);
   

   
   for (int i=0; i<numberOfPoints; i++) {
      printf("At %.0f MeV, GATE75 = %.2f.\n", energies[i], rangesGATE[i]);
      printf("At %.0f MeV, GATE74 = %.2f.\n", energies[i], rangesGATE74[i]);
      printf("At %.0f MeV, GATE73 = %.2f.\n", energies[i], rangesGATE73[i]);
      printf("At %.0f MeV, GATE72 = %.2f.\n", energies[i], rangesGATE72[i]);
      printf("At %.0f MeV, GATE71 = %.2f.\n", energies[i], rangesGATE71[i]);
   }

   gGATE75diff->Draw("LA");
   gGATE74diff->Draw("L, same");
   gGATE73diff->Draw("L, same");
   gGATE72diff->Draw("L, same");
   gGATE71diff->Draw("L, same");
   
   // make text labels
   Int_t energy_pos = 233;
   Float_t subvalue = 0.02;
   TLatex *text75 = new TLatex(energy_pos, gGATE75diff->Eval(230) - subvalue, "75 eV");
   TLatex *text74 = new TLatex(energy_pos, gGATE74diff->Eval(230) - subvalue, "74 eV");
   TLatex *text73 = new TLatex(energy_pos, gGATE73diff->Eval(230) - subvalue, "73 eV");
   TLatex *text72 = new TLatex(energy_pos, gGATE72diff->Eval(230) - subvalue, "72 eV");
   TLatex *text71 = new TLatex(energy_pos, gGATE71diff->Eval(230) - subvalue, "71 eV");

   Float_t latexTextSize = 0.03;
   text75->SetTextSize(latexTextSize);
   text75->SetTextFont(22);
   text74->SetTextSize(latexTextSize);
   text74->SetTextFont(22);
   text73->SetTextSize(latexTextSize);
   text73->SetTextFont(22);
   text72->SetTextSize(latexTextSize);
   text72->SetTextFont(22);
   text71->SetTextSize(latexTextSize);
   text71->SetTextFont(22);

   text75->Draw();
   text74->Draw();
   text73->Draw();
   text72->Draw();
   text71->Draw();


/*
   TLegend *legIP = new TLegend(0.17, 0.72, 0.35, 0.88);
   legIP->AddEntry(gGATE75diff, "GATE I = 75 eV", "l");
   legIP->AddEntry(gGATE74diff, "GATE I = 74 eV", "l");
   legIP->AddEntry(gGATE73diff, "GATE I = 73 eV", "l");
   legIP->AddEntry(gGATE72diff, "GATE I = 72 eV", "l");
   legIP->AddEntry(gGATE71diff, "GATE I = 71 eV", "l");
//   legIP->AddEntry(gMCNP75diff, "MCNP", "L");
   legIP->Draw();
*/

   // Draw fraction of nuclear interactions
   c3->cd();
   Float_t dummyArray[50] = {};
   
   TGraph         *gGATEFractionNI  = new TGraph(numberOfPoints, energies, fractionNIGATE);
   TGraph         *gMCNPFractionNI  = new TGraph(numberOfPoints, energies, fractionNIMCNP);
   TGraph         *gFLUKAFractionNI = new TGraph(numberOfPoints, energies, fractionNIFLUKA);
   TGraphErrors   *gJanniFractionNI = new TGraphErrors(12, energiesJanni, fractionNIJanni, dummyArray, fractionNIJanniError);

   gPad->SetGridx();
   gPad->SetGridy();

   if (phantom == kAluminium) {
      gGATEFractionNI->SetTitle(";Initial energy [MeV];Fraction of Nuclear Interactions");
   }
   else if (phantom == kWater) {
      gGATEFractionNI->SetTitle(";Initial energy [MeV];Fraction of Nuclear Interactions");
   }
   else if (phantom == kComplex) {
      gGATEFractionNI->SetTitle(";Initial energy [MeV];Fraction of Nuclear Interactions");
   }
   gGATEFractionNI->GetYaxis()->SetRangeUser(0, 0.45);
   gGATEFractionNI->SetLineColor(gateColor);
   gGATEFractionNI->SetLineStyle(gateStyle);
   gGATEFractionNI->SetLineWidth(3);
   gGATEFractionNI->GetYaxis()->SetTitleOffset(1.35);
   gFLUKAFractionNI->SetLineColor(flukaColor);
   gFLUKAFractionNI->SetLineStyle(flukaStyle);
   gFLUKAFractionNI->SetLineWidth(3);
   gMCNPFractionNI->SetLineColor(mcnpColor);
   gMCNPFractionNI->SetLineStyle(mcnpStyle);
   gMCNPFractionNI->SetLineWidth(3);
   gJanniFractionNI->SetMarkerColor(expColor);
   gJanniFractionNI->SetMarkerStyle(21);

   gGATEFractionNI->GetXaxis()->SetTitleFont(22);
   gGATEFractionNI->GetXaxis()->SetLabelFont(22);
   gGATEFractionNI->GetYaxis()->SetTitleFont(22);
   gGATEFractionNI->GetYaxis()->SetLabelFont(22);

   gGATEFractionNI->Draw("AL");
   gFLUKAFractionNI->Draw("same, L");
   gMCNPFractionNI->Draw("same, L");
   gJanniFractionNI->Draw("same, P");

   TLegend *legFraction = new TLegend(0.17, 0.72, 0.35, 0.88);
   legFraction->AddEntry(gGATEFractionNI, "GATE", "l");
   legFraction->AddEntry(gMCNPFractionNI, "MCNP", "l");
   legFraction->AddEntry(gFLUKAFractionNI, "FLUKA", "l");
   if (phantom != kComplex) legFraction->AddEntry(gJanniFractionNI, "Janni", "Ple");
   legFraction->SetTextFont(22);
   legFraction->Draw();

   // Draw lateral BP distribution FWHM
   c4->cd();
   TGraph *gGATELateralFWHM = new TGraph(numberOfPoints, energies, fwhmGATE);
   TGraph *gFLUKALateralFWHM = new TGraph(numberOfPoints, energies, fwhmFLUKA);
   TGraph *gMCNPLateralFWHM = new TGraph(numberOfPoints, energies, fwhmMCNP);
   TGraph *gTheoryLateralFWHM = new TGraph(numberOfPoints, energies, beamSpreadingTheory);

   if (phantom == kWater) {
      gGATELateralFWHM->SetTitle(";Initial Energy [MeV];Beam Spreading #sigma_{RMS} / R");
   }
   else if (phantom == kAluminium) {
      gGATELateralFWHM->SetTitle(";Initial Energy [MeV];Beam Spreading #sigma_{RMS} / R");
   }
   else if (phantom == kComplex) {
      gGATELateralFWHM->SetTitle(";Initial Energy [MeV];Beam Spreading #sigma_{RMS} / R");
   }

   gPad->SetGridx();
   gPad->SetGridy();

   gGATELateralFWHM->GetYaxis()->SetTitleOffset(1.4);
   gGATELateralFWHM->GetYaxis()->SetRangeUser(0, 0.07);

   gGATELateralFWHM->GetXaxis()->SetTitleFont(22);
   gGATELateralFWHM->GetXaxis()->SetLabelFont(22);
   gGATELateralFWHM->GetYaxis()->SetTitleFont(22);
   gGATELateralFWHM->GetYaxis()->SetLabelFont(22);

   gGATELateralFWHM->SetLineColor(gateColor);
   gGATELateralFWHM->SetLineStyle(gateStyle);
   gGATELateralFWHM->SetLineWidth(3);
   gFLUKALateralFWHM->SetLineColor(flukaColor);
   gFLUKALateralFWHM->SetLineStyle(flukaStyle);
   gFLUKALateralFWHM->SetLineWidth(3);
   gMCNPLateralFWHM->SetLineColor(mcnpColor);
   gMCNPLateralFWHM->SetLineStyle(mcnpStyle);
   gMCNPLateralFWHM->SetLineWidth(3);
   gTheoryLateralFWHM->SetLineColor(kBlack);
   gTheoryLateralFWHM->SetLineWidth(3);

   gGATELateralFWHM->Draw("AL");
   gFLUKALateralFWHM->Draw("same, L");
   gMCNPLateralFWHM->Draw("same, L");
//   gTheoryLateralFWHM->Draw("same, L");

   TLegend *legFWHM = new TLegend(0.17, 0.72, 0.35, 0.88);
   legFWHM->AddEntry(gGATELateralFWHM, "GATE", "l");
   legFWHM->AddEntry(gMCNPLateralFWHM, "MCNP", "l");
   legFWHM->AddEntry(gFLUKALateralFWHM, "FLUKA", "l");
//   legFWHM->AddEntry(gTheoryLateralFWHM, "Preston and Koehler", "l");
   legFWHM->SetTextFont(22);
   legFWHM->Draw();

   for (int i=0; i<19; i++) {
      printf("Lateral FWHM @ %.0f MeV: %.2f cm (FLUKA) - %.2f cm (GATE) - %.2f cm (MCNP).\n", energies[i], fwhmFLUKA[i], fwhmGATE[i], fwhmMCNP[i]);
   }

   // Draw straggling distribution
   c5->cd();
   TGraph *gGATEStraggling = new TGraph(numberOfPoints, energies, sigmaGATE);
   TGraph *gMCNPStraggling = new TGraph(numberOfPoints, energies, sigmaMCNP);
   TGraph *gFLUKAStraggling = new TGraph(numberOfPoints, energies, sigmaFLUKA);
   TGraphErrors *gJanniStraggling = nullptr;

   // find A and P
   printf("Before fitting. a = %.3f, p = %.3f.\n", a, p);
   Float_t a_fit, p_fit;
   TF1 *BKfunction = new TF1("BKfunction", "[0] * pow(x, [1])");
   BKfunction->SetParameters(a, p);
   a_fit = BKfunction->GetParameter(0)*10;
   p_fit = BKfunction->GetParameter(1);
   gGATE->Fit("BKfunction","N,Q,WW");
   printf("After fitting. a = %.3f, p = %.3f.\n", a_fit, p_fit);

   cout << "a" << endl;
   cout << "a = " << a << ", p = " << p << endl;
   if (phantom == kComplex) {
      for (Int_t i=0; i<numberOfPoints; i++) {
         cout << i << endl;
         Float_t first = pow(p_fit, 2) * pow(a_fit, 2/p_fit);
         Float_t second = 3-2/p_fit;
         Float_t third = pow(rangesFLUKA[i], 3-2/p_fit);
         stragglingJanni[i] = sqrt(ap * first / second * third);
      }
      cout << "Making TGE\n";
      gJanniStraggling = new TGraphErrors(numberOfPoints, rangesFLUKA, stragglingJanni, dummyArray, dummyArray);
   }
   else {
      gJanniStraggling = new TGraphErrors(12, energiesJanni, stragglingJanni, dummyArray, stragglingJanniError);
   }

   gPad->SetGridx();
   gPad->SetGridy();
   if (phantom == kComplex) {
      gGATEStraggling->SetTitle(";Initial energy [MeV];Range straggling [mm]");
   }
   else if (phantom == kAluminium) {
      gGATEStraggling->SetTitle(";Initial energy [MeV];Range straggling [mm]");
   }
   else if (phantom == kWater) {
      gGATEStraggling->SetTitle(";Initial energy [MeV];Range straggling [mm]");
   }

   gGATEStraggling->SetLineColor(gateColor);
   gGATEStraggling->SetLineStyle(gateStyle);
   gGATEStraggling->SetLineWidth(3);
   gFLUKAStraggling->SetLineColor(flukaColor);
   gFLUKAStraggling->SetLineStyle(flukaStyle);
   gFLUKAStraggling->SetLineWidth(3);
   gMCNPStraggling->SetLineColor(mcnpColor);
   gMCNPStraggling->SetLineStyle(mcnpStyle);
   gMCNPStraggling->SetLineWidth(3);
   gJanniStraggling->SetMarkerColor(expColor);
   gJanniStraggling->SetMarkerStyle(21);

   gGATEStraggling->GetXaxis()->SetTitleFont(22);
   gGATEStraggling->GetYaxis()->SetTitleFont(22);
   gGATEStraggling->GetXaxis()->SetLabelFont(22);
   gGATEStraggling->GetYaxis()->SetLabelFont(22);

   gGATEStraggling->Draw("AL");
   gFLUKAStraggling->Draw("same, L");
   gMCNPStraggling->Draw("same, L");
   gJanniStraggling->Draw("same, P");

   TLegend *legStraggling = new TLegend(0.17, 0.72, 0.35, 0.88);
   legStraggling->AddEntry(gGATEStraggling, "GATE", "l");
   legStraggling->AddEntry(gMCNPStraggling, "MCNP", "l");
   legStraggling->AddEntry(gFLUKAStraggling, "FLUKA", "l");
   if (phantom != kComplex) legStraggling->AddEntry(gJanniStraggling, "Janni", "PLE");
   legStraggling->SetTextFont(22);
   legStraggling->Draw();

   // Direct range comparisons
   c7->cd();
   gPad->SetGridy();
   gPad->SetGridx();
   TGraph *gGATEvsPSTAR = new TGraph(numberOfPoints, energies, rangesGATEdiff);
   TGraph *gMCNPvsPSTAR = new TGraph(numberOfPoints, energies, rangesMCNPdiff);
   TGraph *gFLUKAvsPSTAR = new TGraph(numberOfPoints, energies, rangesFLUKAdiff);

   if (phantom == kComplex) {
      gMCNPvsPSTAR->SetTitle(";Initial energy [MeV];Range deviation* [mm]");
   }
   else if (phantom == kAluminium) {
      gMCNPvsPSTAR->SetTitle(";Initial energy [MeV];Range deviation* [mm]");
   }
   else if (phantom == kWater) {
      gMCNPvsPSTAR->SetTitle(";Initial energy [MeV];Range deviation* [mm]");
   }

   gMCNPvsPSTAR->GetYaxis()->SetTitleOffset(1.3);
   gMCNPvsPSTAR->SetLineColor(mcnpColor);
   gMCNPvsPSTAR->SetLineStyle(mcnpStyle);
   gMCNPvsPSTAR->SetLineWidth(3);
   gMCNPvsPSTAR->GetYaxis()->SetRangeUser(-2, 0.6);
   gMCNPvsPSTAR->GetXaxis()->SetNdivisions(12);
   gMCNPvsPSTAR->GetXaxis()->SetLabelFont(22);
   gMCNPvsPSTAR->GetXaxis()->SetTitleFont(22);
   gMCNPvsPSTAR->GetYaxis()->SetLabelFont(22);
   gMCNPvsPSTAR->GetYaxis()->SetTitleFont(22);
   gMCNPvsPSTAR->Draw("AL");
   
   gGATEvsPSTAR->SetLineColor(gateColor);
   gGATEvsPSTAR->SetLineStyle(gateStyle);
   gGATEvsPSTAR->SetLineWidth(3);
   gGATEvsPSTAR->Draw("same, L");

   gFLUKAvsPSTAR->SetLineColor(flukaColor);
   gFLUKAvsPSTAR->SetLineStyle(flukaStyle);
   gFLUKAvsPSTAR->SetLineWidth(3);
   gFLUKAvsPSTAR->Draw("same, L");

   TLegend *leg3 = new TLegend(0.17, 0.72, 0.35, 0.88);
   leg3->AddEntry(gGATEvsPSTAR, "GATE", "l");
   leg3->AddEntry(gMCNPvsPSTAR, "MCNP", "l");
   leg3->AddEntry(gFLUKAvsPSTAR, "FLUKA", "l");
   leg3->SetTextFont(22);
   leg3->Draw();

   // MAKE DATA TABLES
   printf("Geometry; Monte Carlo code; Energy [MeV]; Range [cm]; Range error [mm]; Range straggling [mm]; Beam spreading [1]; Fraction of inelastic collisions [%]\n");
   for (Int_t i=0; i<numberOfPoints; i++) {
      if (phantom == kWater) {
         printf("Water; GATE 71 eV; %.0f; %.5f; %.5f; 0; 0; 0\n", energies[i], rangesGATE71[i], rangesGATEdiff71[i]);
         printf("Water; GATE 72 eV; %.0f; %.5f; %.5f; 0; 0; 0\n", energies[i], rangesGATE72[i], rangesGATEdiff72[i]);
         printf("Water; GATE 73 eV; %.0f; %.5f; %.5f; 0; 0; 0\n", energies[i], rangesGATE73[i], rangesGATEdiff73[i]);
         printf("Water; GATE 74 eV; %.0f; %.5f; %.5f; 0; 0; 0\n", energies[i], rangesGATE74[i], rangesGATEdiff74[i]);
         printf("Water; GATE 75 eV; %.0f; %.5f; %.5f; 0; 0; 0\n", energies[i], rangesGATE[i], rangesGATEdiff75[i]);
         printf("Water; GATE (standard); %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesGATE[i], rangesGATEdiff[i], sigmaGATE[i], fwhmGATE[i], fractionNIGATE[i]);
         printf("Water; FLUKA; %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesFLUKA[i], rangesFLUKAdiff[i], sigmaFLUKA[i], fwhmFLUKA[i], fractionNIFLUKA[i]);
         printf("Water; MCNP; %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesMCNP[i], rangesMCNPdiff[i], sigmaMCNP[i], fwhmMCNP[i], fractionNIMCNP[i]);
         printf("Water; PSTAR; %.0f; %.5f; 0; 0; 0\n", energies[i], rangesPSTAR[i]);
         if (i < 13) printf("Water; JANNI; %.0f; %.5f; %.5f; 0; %.5f\n", energiesJanni[i], rangesJanni[i], stragglingJanni[i], fractionNIJanni[i]);
      }
      if (phantom == kComplex) {
         printf("DTC; GATE; %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesGATE[i], rangesGATEdiff[i], sigmaGATE[i], fwhmGATE[i], fractionNIGATE[i]);
         printf("DTC; FLUKA; %.0f; %.5f; %.5f; %.5f;  %.5f; %.5f\n", energies[i], rangesFLUKA[i], rangesFLUKAdiff[i], sigmaFLUKA[i], fwhmFLUKA[i], fractionNIFLUKA[i]);
         printf("DTC; MCNP; %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesMCNP[i], rangesMCNPdiff[i], sigmaMCNP[i], fwhmMCNP[i], fractionNIMCNP[i]);
         printf("DTC; PSTAR; %.0f; %.5f; 0; 0; 0\n", energies[i], rangesPSTAR[i]);
         if (i < 13) printf("DTC; Bortfeld; %.0f; %.5f\n", energiesJanni[i], stragglingJanni[i]);
      }
      if (phantom == kAluminium) {
         printf("Aluminium; GATE; %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesGATE[i], rangesGATEdiff[i], sigmaGATE[i], fwhmGATE[i], fractionNIGATE[i]);
         printf("Aluminium; FLUKA; %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesFLUKA[i], rangesFLUKAdiff[i], sigmaFLUKA[i], fwhmFLUKA[i], fractionNIFLUKA[i]);
         printf("Aluminium; MCNP; %.0f; %.5f; %.5f; %.5f; %.5f; %.5f\n", energies[i], rangesMCNP[i], rangesMCNPdiff[i], sigmaMCNP[i], fwhmMCNP[i], fractionNIMCNP[i]);
         printf("Aluminium; PSTAR; %.0f; %.5f; 0; 0; 0\n", energies[i], rangesPSTAR[i]);
         if (i < 13) printf("Aluminium; JANNI; %.0f; %.5f; %.5f; 0; %.5f\n", energiesJanni[i], rangesJanni[i], stragglingJanni[i], fractionNIJanni[i]);
      }
   }

   c8->cd();
   TGraph *gDetourGATE = new TGraph(numberOfPoints, energies, detourGATE);
   TGraph *gDetourFLUKA = new TGraph(numberOfPoints, energies, detourFLUKA);
   TGraph *gDetourMCNP = new TGraph(numberOfPoints, energies, detourMCNP);
   TGraph *gDetourPSTARWater = new TGraph(numberOfPoints, energies, detourPSTARWater);
   TGraph *gDetourPSTARAl = new TGraph(numberOfPoints, energies, detourPSTARAl);
   
   gDetourGATE->SetTitle("Detour factors;Energy [MeV];Detour Factor");

   gDetourGATE->SetLineColor(gateColor);
   gDetourGATE->SetLineWidth(3);
   gDetourFLUKA->SetLineColor(flukaColor);
   gDetourFLUKA->SetLineWidth(3);
   gDetourMCNP->SetLineColor(mcnpColor);
   gDetourMCNP->SetLineWidth(3);
   gDetourPSTARAl->SetMarkerColor(expColor);
   gDetourPSTARAl->SetMarkerStyle(21);
   gDetourPSTARWater->SetMarkerColor(expColor);
   gDetourPSTARWater->SetMarkerStyle(21);

   gDetourGATE->Draw("AL");
   gDetourFLUKA->Draw("same, L");
   gDetourMCNP->Draw("same, L");
   
   if (phantom == kWater) {
      gDetourPSTARWater->Draw("same, P");
   }

   else if (phantom == kAluminium) {
      gDetourPSTARAl->Draw("same, P");
   }
   
   TLegend *legDetour = new TLegend(0.17, 0.72, 0.35, 0.88);
   legDetour->AddEntry(gDetourGATE, "GATE", "l");
   legDetour->AddEntry(gDetourMCNP, "MCNP", "l");
   legDetour->AddEntry(gDetourFLUKA, "FLUKA", "l");
   if (phantom == kWater) {
      legDetour->AddEntry(gDetourPSTARWater, "PSTAR", "PLE");
   }
   else if (phantom == kAluminium) {
      legDetour->AddEntry(gDetourPSTARAl, "PSTAR", "PLE");
   }

   legDetour->Draw();

   c2->cd();
   TLegend *legLateral = new TLegend(0.17, 0.72, 0.35, 0.88);
   legLateral->AddEntry(hGATELateral, "GATE", "l");
   legLateral->AddEntry(hMCNPLateral, "MCNP", "l");
   legLateral->AddEntry(hFLUKALateral, "FLUKA", "l");
   legLateral->SetTextFont(22);
   legLateral->Draw();

}
