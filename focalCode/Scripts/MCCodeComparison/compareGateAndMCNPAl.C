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
   Bool_t activateFLUKA = false;
   Bool_t activateTRIM = false;

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
   Float_t  rangesPSTAR[19] = {1.08, 1.5, 1.97, 2.49, 3.07, 3.7, 4.37, 5.09, 5.85, 6.66, 7.51, 8.4, 9.32, 10.28, 11.28, 12.31, 13.38, 14.47, 15.60};
//   Float_t  rangesPSTAR[19] = {1.08, 1.97, 3.07, 4.37, 5.85, 7.51, 9.32, 11.28, 13.38, 15.60, 17.94, 20.40, 22.97, 25.64};
   Float_t  sigmaPSTAR[19] = {};

   for (Int_t i=0; i<19; i++) {
      energies[i] = (i+5)*10;
//      rangesPSTAR[i] *= 10;
      sigmaPSTAR[i] = rangesPSTAR[i] * pow(energies[i], -0.104) * 0.0188; // FIT TO YANNI
   }

   for (Int_t i=0; i<19; i++) {
      energiesMCNP[i] = energies[i] - 3;
      energiesGATE[i] = energies[i];
      energiesFLUKA[i] = energies[i] + 3;
      energiesPSTAR[i] = energies[i] + 6;
   }
   
   TCanvas *c2 = new TCanvas("c2", "test", 800, 800);
   TCanvas *c1 = new TCanvas("c1", "Ranges for all codes", 1200, 1000);
   TPad *pad1 = new TPad("pad1", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad2 = new TPad("pad2", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad1->Draw();
   pad2->Draw();

   gPad->SetFillColor(50);
   gPad->Modified();

   if (activateGATE) {
      Float_t  x,y,z,edep;
      Int_t    parentID, eventID;
      Bool_t   isInelastic;
      
      for (Int_t i=0; i<18; i++) {
         nominalEnergy = (i+5)*10;
         cout << "Nominal energy " << nominalEnergy << endl;
         nominalRange = 0.0012 * pow(nominalEnergy, 1.7483);
         TH1F *hGATE = new TH1F("hGATE", "Proton ranges in single GATE dataset;Range [cm];Number of primaries", 400, fmin(0, nominalRange - 5), nominalRange + 5);
         
         cout << "Reading file " << Form("Data/GATE/Aluminium/compressed_aluminium_%dMeV.root\n", nominalEnergy);
         TFile   *f1 = new TFile(Form("Data/GATE/Aluminium/compressed_aluminium_%dMeV.root", nominalEnergy));
         cout << "Opening tree...\n";

         TTree   *treeBic = (TTree*) f1->Get("treeOut");

         cout << "Setting addresses...\n";
         treeBic->SetBranchAddress("posX",&x);
         treeBic->SetBranchAddress("posY",&y);
         treeBic->SetBranchAddress("posZ",&z);
         treeBic->SetBranchAddress("edep",&edep);
         treeBic->SetBranchAddress("eventID",&eventID);
         treeBic->SetBranchAddress("parentID",&parentID);
         treeBic->SetBranchAddress("isInelastic",&isInelastic);

         printf("Reading %d entries from GATE file.\n", treeBic->GetEntries());
         for (Int_t j=0, N = treeBic->GetEntries(); j<N; ++j) {
            treeBic->GetEntry(j);
            if (!isInelastic) hGATE->Fill(z/10.0);
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

         printf("before fit GATE, mu = %.2f, sigma = %.2f.\n", mu, sigma);
      
         TF1 *fit = new TF1("fit", "gaus");
         fit->SetParameter(1, mu); fit->SetParameter(2, sigma);
         fit->SetParameter(0, 1000);
         fit->SetParLimits(1, mu*0.75, mu*1.5);
         fit->SetParLimits(2, sigma*0.75, sigma*1.5);
         fit->SetNpx(1000);

         hGATE->Fit("fit", "N,M,B");
         mu = fit->GetParameter(1);
         sigma = fit->GetParameter(2);
         printf("after fit GATE, mu = %.2f, sigma = %.2f.\n", mu, sigma);

         cout << "IN " << nominalEnergy << " DATASET, MU = " << mu << ", SIGMA = " << sigma << endl;
         rangesGATE[i] = mu;
         sigmaGATE[i] = sigma;
            
         delete f1;
         if (nominalEnergy == 110) {
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
      Float_t  startZ = -40;
      string   line;
      Float_t  x, y, z;
      Float_t  u, v, w; // vectors
      Float_t  energy, weight, time;
      Float_t  nFilling = 0;

      cout << "READING MCNP FILES...\n";
      for (Int_t i=0; i<19; i++) {
         nominalEnergy = (i+5) * 10;
         cout << "Nominal energy " << nominalEnergy << endl;
         nominalRange = 0.0012 * pow(nominalEnergy, 1.7341);
         cout << nominalEnergy << "... ";
         TH1F    *hMCNP = new TH1F("hMCNP", "Proton ranges in single MCNP dataset;Range [cm];Number of primaries", 400, fmin(0, nominalRange - 5), nominalRange + 5);

         in.open(Form("Data/MCNP/Aluminium/%dMeV_Alp", nominalEnergy));

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

                  hMCNP->Fill((z - startZ)/10.);
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

         printf("before fit MCNP, mu = %.2f, sigma = %.2f.\n", mu, sigma);

         TF1 *fit = new TF1("fit", "gaus");
         hMCNP->Fit("fit", "N");
         mu = fit->GetParameter(1);
         sigma = fit->GetParameter(2);
         printf("after fit MCNP, mu = %.2f, sigma = %.2f.\n", mu, sigma);

         rangesMCNP[i] = mu;
         sigmaMCNP[i] = sigma;
         
         delete hMCNP;
      }
   }
/*
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

      TH1F    *hFLUKAPrecisio = new TH1F("hFLUKAPrecisio", "All protons in FLUKA PRECISIO dataset;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH1F    *hFLUKAPrecisioStop = new TH1F("hFLUKAPrecisioStop", "Protons in FLUKA PRECISIO dataset with termination type 11;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH2F    *hFLUKAPrecisioLateral = new TH2F("hFLUKAPrecisioLateral", "2D distribution of BP position in FLUKA PRECISIO dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);
      TH1F    *hFLUKAHadrothe = new TH1F("hFLUKAHadrothe", "All protons in FLUKA HADROTHE dataset;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH1F    *hFLUKAHadrotheStop = new TH1F("hFLUKAHadrotheStop", "Protons in FLUKA HADROTHE dataset with termination type 11;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH2F    *hFLUKAHadrotheLateral = new TH2F("hFLUKAHadrotheLateral", "2D distribution of BP position in FLUKA HADROTHE dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);
      TH1F    *hFLUKACalorime = new TH1F("hFLUKACalorime", "All protons in FLUKA CALORIME dataset;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH1F    *hFLUKACalorimeStop = new TH1F("hFLUKACalorimeStop", "Protons in FLUKA CALORIME dataset with termination type 11;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH2F    *hFLUKACalorimeLateral = new TH2F("hFLUKACalorimeLateral", "2D distribution of BP position in FLUKA CALORIME dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);

      cout << "READING PRECISIO FILE\n";
      in.open("Data/FLUKA/FLUKA_PRECISIO.txt");

      while (! in.eof() ) {
         in >> historyNumber >> x >> y >> z >> terminationType;

         if (lastHistoryNumber < 0) lastHistoryNumber = historyNumber;

         if (historyNumber != lastHistoryNumber && terminationType != 23) { // primary particle
            hFLUKAPrecisio->Fill(z);

            if (terminationType == 11) {
               hFLUKAPrecisioStop->Fill(z);
            }

            else {
               hFLUKAPrecisioLateral->Fill(10*x,10*y);
            }
         }

         lastHistoryNumber = historyNumber; // all following entries with same history number are 2ndies
      }

      in.close();

      lastHistoryNumber = -1;
      in.open("Data/FLUKA/FLUKA_HADROTHE.txt");
      while (! in.eof() ) {
         in >> historyNumber >> x >> y >> z >> terminationType;

         if (lastHistoryNumber < 0) lastHistoryNumber = historyNumber;

         if (historyNumber != lastHistoryNumber && terminationType != 23) { // primary particle
            hFLUKAHadrothe->Fill(z);

            if (terminationType == 11) {
               hFLUKAHadrotheStop->Fill(z);
            }

            else {
               hFLUKAHadrotheLateral->Fill(10*x,10*y);
            }
         }

         lastHistoryNumber = historyNumber; // all following entries with same history number are 2ndies
      }
      in.close();

      lastHistoryNumber = -1;
      in.open("Data/FLUKA/FLUKA_CALORIME.txt");
//      in.open("Data/FLUKA/WATER_PRECISIO.txt");
      while (! in.eof() ) {
         in >> historyNumber >> x >> y >> z >> terminationType;

         if (lastHistoryNumber < 0) lastHistoryNumber = historyNumber;

         if (historyNumber != lastHistoryNumber && terminationType != 23) { // primary particle
            hFLUKACalorime->Fill(z);

            if (terminationType == 11) {
               hFLUKACalorimeStop->Fill(z);
            }

            else {
               hFLUKACalorimeLateral->Fill(10*x,10*y);
            }
         }

         lastHistoryNumber = historyNumber; // all following entries with same history number are 2ndies
      }
      in.close();

      cout << "PLOTTING\n";

      c3->cd(1);
      hFLUKAPrecisio->SetFillColor(kGreen-3);
      hFLUKAPrecisio->Draw();

      fitFunction = new TF1("fitFunction", "gaus");
      fitFunction->SetParameter(1, 11.5);
      fitFunction->SetParameter(2, 0.15);
      fitFunction->SetParLimits(1, 9, 13);
      fitFunction->SetParLimits(2, 0.05, 0.4);

      hFLUKAPrecisio->Fit("fitFunction", "B, M,Q");

      mu = fitFunction->GetParameter(1);
      sigma = fabs(fitFunction->GetParameter(2));

      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hFLUKAPrecisio->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hFLUKAPrecisio->Integral();
      cutoffIntegral = hFLUKAPrecisio->Integral(0, distributionCutoffBin);
      stoppedIntegral = hFLUKAPrecisioStop->Integral();
      fraction = cutoffIntegral / float(totalIntegral);

      TLegend *leg4 = new TLegend(.18, .7, .68, .84);
      leg4->AddEntry(hFLUKAPrecisio, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      leg4->AddEntry(fitFunction, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      leg4->Draw();

      cout << Form("Total number of events: %d. Number of events with terminationType == (13, 16): %.0f. Fraction: %.2f.", totalIntegral, hFLUKAPrecisioStop->Integral(), fraction) << endl;
      
      c3->cd(2);
      hFLUKAPrecisioStop->SetFillColor(kGreen-3);
      hFLUKAPrecisioStop->Draw();

      c3->cd(3);
      hFLUKAPrecisioLateral->Draw("COLZ");

      c3->cd(4);
      hFLUKAHadrothe->SetFillColor(kGreen-3);
      hFLUKAHadrothe->Draw();

      fitFunction = new TF1("fitFunction", "gaus");
      fitFunction->SetParameter(1, 11.5);
      fitFunction->SetParameter(2, 0.15);
      fitFunction->SetParLimits(1, 9, 13);
      fitFunction->SetParLimits(2, 0.05, 0.4);

      hFLUKAHadrothe->Fit("fitFunction", "B, M,Q");

      mu = fitFunction->GetParameter(1);
      sigma = fabs(fitFunction->GetParameter(2));

      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hFLUKAHadrothe->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hFLUKAHadrothe->Integral();
      cutoffIntegral = hFLUKAHadrothe->Integral(0, distributionCutoffBin);
      stoppedIntegral = hFLUKAHadrotheStop->Integral();
      fraction = cutoffIntegral / float(totalIntegral);

      TLegend *leg5 = new TLegend(.18, .7, .68, .84);
      leg5->AddEntry(hFLUKAHadrothe, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      leg5->AddEntry(fitFunction, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      leg5->Draw();

      cout << Form("Total number of events: %d. Number of events with terminationType == (13, 16): %.0f. Fraction: %.2f.", totalIntegral, hFLUKAHadrotheStop->Integral(), fraction) << endl;
      
      c3->cd(5);
      hFLUKAHadrotheStop->SetFillColor(kGreen-3);
      hFLUKAHadrotheStop->Draw();

      c3->cd(6);
      hFLUKAHadrotheLateral->Draw("COLZ");

      c3->cd(7);
      hFLUKACalorime->SetFillColor(kGreen-3);
      hFLUKACalorime->Draw();

      fitFunction = new TF1("fitFunction", "gaus");
      fitFunction->SetParameter(1, 11.5);
      fitFunction->SetParameter(2, 0.15);
      fitFunction->SetParLimits(1, 9, 13);
      fitFunction->SetParLimits(2, 0.05, 0.4);

      hFLUKACalorime->Fit("fitFunction", "B, M,Q");

      mu = fitFunction->GetParameter(1);
      sigma = fabs(fitFunction->GetParameter(2));

      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hFLUKACalorime->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hFLUKACalorime->Integral();
      cutoffIntegral = hFLUKACalorime->Integral(0, distributionCutoffBin);
      stoppedIntegral = hFLUKACalorimeStop->Integral();
      fraction = cutoffIntegral / float(totalIntegral);

      TLegend *leg6 = new TLegend(.18, .7, .68, .84);
      leg6->AddEntry(hFLUKACalorime, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      leg6->AddEntry(fitFunction, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      leg6->Draw();

      cout << Form("Total number of events: %d. Number of events with terminationType == (13, 16): %.0f. Fraction: %.2f.", totalIntegral, hFLUKACalorimeStop->Integral(), fraction) << endl;
      
      c3->cd(8);
      hFLUKACalorimeStop->SetFillColor(kGreen-3);
      hFLUKACalorimeStop->Draw();

      c3->cd(9);
      hFLUKACalorimeLateral->Draw("COLZ");

      TH1D *hProfile = hFLUKACalorimeLateral->ProjectionX();
      int bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      int bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      double fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM FLUKA CALORIMETER 2D = " << fwhm << " mm." << endl;
      hProfile = hFLUKAPrecisioLateral->ProjectionX();
      bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM FLUKA PRECISION 2D = " << fwhm << " mm." << endl;
      hProfile = hFLUKAHadrotheLateral->ProjectionX();
      bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      cout << "Maximum/2 for HADRO = " << hProfile->GetMaximum()/2 << endl;
      fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM FLUKA HADROTHERAPY 2D = " << fwhm << " mm." << endl;
      TCanvas *c10 = new TCanvas("c10", "c10", 800, 800);
      c10->cd(); hProfile->Draw();
   }
   */

   cout << "Drawing in pad1\n";
   pad1->cd();
   cout << "Opened pad1\n";

   Float_t errorScalingFactor = 10;
   for (int i=0; i<19; i++) {
      sigmaMCNP[i] *= errorScalingFactor;
      sigmaGATE[i] *= errorScalingFactor;
      sigmaPSTAR[i] *= errorScalingFactor;
      sigmaFLUKA[i] *= errorScalingFactor;

      rangesMCNPdiff[i] = (rangesMCNP[i] - rangesPSTAR[i]) * 10;
      rangesGATEdiff[i] = (rangesGATE[i] - rangesPSTAR[i]) * 10;
      rangesFLUKAdiff[i] = (rangesFLUKA[i] - rangesPSTAR[i] * 10);
   }


   TGraphErrors *gMCNP = new TGraphErrors(19, energiesMCNP, rangesMCNP, energyError, sigmaMCNP);
   gMCNP->SetTitle(Form("Aluminum range comparison between different codes (%.1fx error bars); Energy [MeV];Range [cm]", errorScalingFactor));
   gMCNP->SetMarkerColor(kRed);
   gMCNP->SetMarkerStyle(22);
   gMCNP->Draw("AP");
   
   TGraphErrors *gGATE = new TGraphErrors(19, energiesGATE, rangesGATE, energyError, sigmaGATE);
   gGATE->SetMarkerColor(kBlue);
   gGATE->SetMarkerStyle(22);
   gGATE->Draw("same, P");
   
   TGraphErrors *gPSTAR = new TGraphErrors(19, energiesPSTAR, rangesPSTAR, energyError, sigmaPSTAR);
   gPSTAR->SetMarkerColor(kBlack);
   gPSTAR->SetMarkerStyle(22);
   gPSTAR->Draw("same, P");

   TLegend *leg = new TLegend(0.17, 0.72, 0.35, 0.88);
   leg->AddEntry(gMCNP, "MCNP", "Pel");
   leg->AddEntry(gGATE, "GATE", "Pel");
   leg->AddEntry(gPSTAR, "PSTAR", "Ple");
   leg->Draw();

   pad2->cd();
   cout << "Drawing in pad2\n";

   TGraph *gMCNPdiff = new TGraph(19, energiesMCNP, rangesMCNPdiff);
   gMCNPdiff->SetTitle("Difference between MC codes and PSTAR range; Energy [MeV];Range [mm]");
   gMCNPdiff->SetMarkerColor(kRed);
   gMCNPdiff->SetMarkerStyle(22);
   gMCNPdiff->GetXaxis()->SetTitleSize(0.09);
   gMCNPdiff->GetYaxis()->SetTitleSize(0.09);
   gMCNPdiff->GetYaxis()->SetTitleOffset(0.35);
   gMCNPdiff->GetXaxis()->SetTitleOffset(0.8);
   gMCNPdiff->GetXaxis()->SetLabelSize(0.08);
   gMCNPdiff->GetYaxis()->SetLabelSize(0.08);
   gMCNPdiff->Draw("AL");
   
   TGraph *gGATEdiff = new TGraph(19, energiesGATE, rangesGATEdiff);
   gGATEdiff->SetMarkerColor(kBlue);
   gGATEdiff->SetMarkerStyle(22);
   gGATEdiff->Draw("same, L");
}
