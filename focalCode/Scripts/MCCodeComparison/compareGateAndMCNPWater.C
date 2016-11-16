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

void Run()
{
   // GATE PART
   Bool_t activateGATE = true;
   Bool_t activateMCNP = true;
   Bool_t activateFLUKA = true;
   Bool_t activateTRIM = false;

   enum phantomType {kAluminium, kWater, kComplex};

   Int_t phantom = kWater;

   Int_t    nbinsxy = 400;
   Int_t    xyfrom = -60;
   Int_t    xyto = 60;
   Int_t    nbinsz = 1000;
   Int_t    zfrom = 0;
   Int_t    zto = 40;
   Int_t    nMax = 500000;
   TF1     *fitFunction = nullptr;
   Float_t  mu;
   Float_t  sigma; 
   Int_t    nominalEnergy;
   Float_t  nominalRange;
      
   Float_t  energies[19]; 
   Float_t  energiesMCNP[19];
   Float_t  energiesFLUKA[19];
   Float_t  energiesGATE[19];
   Float_t  energiesPSTAR[19];
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
//   Float_t  rangesPSTAR[19] = {2.224, 4.075, 6.389, 9.128, 12.26, 15.76, 19.59, 23.74, 28.19, 32.91, 37.90, 43.12, 48.58, 54.26};
   Float_t  rangesPSTARWater[19] = {2.224, 3.089, 4.075, 5.176, 6.389, 7.707, 9.128, 10.65, 12.26, 13.96, 15.76, 17.63, 19.59, 21.63, 23.74, 25.93, 28.19, 30.52, 32.91};
   Float_t  rangesPSTARAl[19] = {1.08, 1.5, 1.97, 2.49, 3.07, 3.7, 4.37, 5.09, 5.85, 6.66, 7.51, 8.4, 9.32, 10.28, 11.28, 12.31, 13.38, 14.47, 15.60};
   Float_t  rangesPSTAR[19] = {};

   if (phantom == kAluminium) {
      for (Int_t i=0; i<19; i++) {
         rangesPSTAR[i] = rangesPSTARAl[i];
      }
   }
   else if (phantom == kWater) {
      for (Int_t i=0; i<19; i++) {
         rangesPSTAR[i] = rangesPSTARWater[i];
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
            if (!isInelastic) hGATE->Fill(z/10.);
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
//         fit->SetParameter(1, mu);
 //        fit->SetParameter(2, sigma);

  //       fit->SetParLimits(1, mu*0.9, mu*1.1);
    //     fit->SetParLimits(2, 0, 3);

      //   fit->SetNpx(1000);

      //   hGATE->Fit("fit", "N,Q,B,W", "", mu*0.9, mu*1.1);
         hGATE->Fit("fit", "N,Q,WW");
         mu = fit->GetParameter(1);
         sigma = fabs(fit->GetParameter(2));

         printf("GATE, %d MeV. nominal Range is %.1f, fit range is %.1f.\n", nominalEnergy, nominalRange, mu);

         rangesGATE[i] = mu;
         sigmaGATE[i] = sigma;
            
         delete f1;
         if (nominalEnergy == 190) {
            c2->cd();
            hGATE->Draw();
         }
         else {
            delete hGATE;
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
         TH1F    *hMCNP = new TH1F("hMCNP", "Proton ranges in single MCNP dataset;Range [cm];Number of primaries", 400, fmax(0, nominalRange - 5), nominalRange + 5);
         
         if (phantom == kAluminium)  {
            in.open(Form("Data/MCNP/Aluminium/%dMeV_Alp", nominalEnergy));
         }
         else if (phantom == kWater) {
            in.open(Form("Data/MCNP/Water/%dMeV_vannp", nominalEnergy));
         }
         else if (phantom == kComplex) {
            in.open(Form("Data/MCNP/ComplexGeometry/cem03.03/10k_protons/%dMeV_sdetp", nominalEnergy));
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

               if (branchNumber == 1 && terminationType != 13 && terminationType != 16) { // primary particle
                  in >> x >> y >> z >> u >> v >> w >> energy >> weight >> time;

                  hMCNP->Fill((z - startZ));
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


         TF1 *fit = new TF1("fit", "gaus");
         fit->SetParameter(1, mu);
         fit->SetParameter(2, sigma);

         fit->SetParLimits(1, mu-2, mu+2);
         hMCNP->Fit("fit", "N,Q,M,B");
         mu = fit->GetParameter(1);
         sigma = fabs(fit->GetParameter(2));

         cout << "IN " << nominalEnergy << " DATASET, MU = " << mu << ", SIGMA = " << sigma << endl;
         rangesMCNP[i] = mu;
         sigmaMCNP[i] = sigma;
         if (i == 17-5) {
            c2->cd();
//            hMCNP->Draw();
         }
         else {
            delete hMCNP;
         }
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
            nominalRange = 0.0013 * pow(nominalEnergy, 1.7447);
         }

         TH1F    *hFLUKA= new TH1F("hFLUKA", "All protons in FLUKA dataset;Range [cm];Number of primaries", 100, fmax(0, nominalRange-2), nominalRange+2);

         if (phantom == kAluminium) { 
            in.open(Form("Data/FLUKA/Aluminium/aluminium%dMeV.txt", nominalEnergy));
         }
         else if (phantom == kWater) {
            in.open(Form("Data/FLUKA/Water/water%dMeV.txt", nominalEnergy));
         }
         else if (phantom == kComplex) {
            in.open(Form("Data/FLUKA/ComplexGeometry/detector%dMeV.txt", nominalEnergy));
         }

         while (! in.eof() ) {
            in >> historyNumber >> x >> y >> z >> terminationType;

            if (!in.good()) break;

            if (terminationType != 11) {
               hFLUKA->Fill(z);
            }
         }

         in.close();

         mu = nominalRange;
         sigma = 0.1;

         TF1 *fit = new TF1("fit", "gaus");
/*         fit->SetParameter(1, mu);
         fit->SetParameter(2, sigma);

         fit->SetParLimits(1, mu-2, mu+2);
         fit->SetParLimits(2, 0, 3); */
         hFLUKA->Fit("fit", "N,Q,M");
         mu = fit->GetParameter(1);
         sigma = fabs(fit->GetParameter(2));

         cout << "IN " << nominalEnergy << " DATASET, MU = " << mu << ", SIGMA = " << sigma << endl;
         rangesFLUKA[i] = mu;
         sigmaFLUKA[i] = sigma;
        
         if (nominalEnergy == 220) {
            c2->cd();
            hFLUKA->Draw();
         }
         else {
            delete hFLUKA;
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

   TGraphErrors *gMCNP = new TGraphErrors(numberOfPoints, energiesMCNP, rangesMCNP, energyError, sigmaMCNP);
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
   gFLUKA->SetMarkerColor(kOrange);
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

      TGraph *gMCNPdiff = new TGraph(19, energiesMCNP, rangesMCNPdiff);
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
      
      TGraph *gGATEdiff = new TGraph(19, energiesGATE, rangesGATEdiff);
      gGATEdiff->SetLineColor(kBlue);
      gGATEdiff->SetLineWidth(3);
      gGATEdiff->Draw("same, L");

      TGraph *gFLUKAdiff = new TGraph(19, energiesFLUKA, rangesFLUKAdiff);
      gFLUKAdiff->SetLineColor(kOrange);
      gFLUKAdiff->SetLineWidth(3);
      gFLUKAdiff->Draw("same,L");

      pad2->Update();

      TLine *l = new TLine(pad2->GetUxmin(), 0, pad2->GetUxmax(), 0);
      l->Draw();
   }
}
