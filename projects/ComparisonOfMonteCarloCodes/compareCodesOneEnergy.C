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

using namespace std;

void Run()
{
   // GATE PART
   Bool_t activateGATE = true;
   Bool_t activateMCNP = true;
   Bool_t activateFLUKA = true;
   Bool_t activateTRIM = false;

   Int_t    nbinsxy = 400;
   Int_t    xyfrom = -60;
   Int_t    xyto = 60;
   Int_t    nbinsz = 400;
   Int_t    zfrom = 0;
   Int_t    zto = 13;
   Int_t    nMax = 500000;
   TF1     *fitFunction = nullptr;
   Float_t  mu;
   Float_t  sigma; 

   Float_t  distributionCutoff;
   Int_t    distributionCutoffBin;
   Int_t    totalIntegral;
   Int_t    cutoffIntegral;
   Int_t    stoppedIntegral;
   Float_t  fraction;
   
   TCanvas *c1 = new TCanvas("c1", "GATE", 1200 , 1050);
   c1->Divide(3,3,0.001,0.001);
   c1->SetLogy();
   
   TCanvas *c2 = new TCanvas("c2", "MCNP", 1200, 1050);
   c2->Divide(3,3,0.001,0.001);
   c2->SetLogy();

   TCanvas *c3 = new TCanvas("c3", "FLUKA", 1200, 1050);
   c3->Divide(3,3,0.001,0.001);

   // Set logy in all pads
   TPad *p1 = (TPad *) c1->cd(1);
   TPad *p2 = (TPad *) c1->cd(2);
   TPad *p4 = (TPad *) c1->cd(4);
   TPad *p5 = (TPad *) c1->cd(5);
   TPad *p7 = (TPad *) c2->cd(1);
   TPad *p8 = (TPad *) c2->cd(2);
   TPad *p10 = (TPad *) c2->cd(4);
   TPad *p11 = (TPad *) c2->cd(5);
   TPad *p13 = (TPad *) c2->cd(7);
   TPad *p14 = (TPad *) c2->cd(8);
   TPad *p16 = (TPad *) c3->cd(1);
   TPad *p17 = (TPad *) c3->cd(2);
   TPad *p19 = (TPad *) c3->cd(4);
   TPad *p20 = (TPad *) c3->cd(5);
   TPad *p22 = (TPad *) c3->cd(7);
   TPad *p23 = (TPad *) c3->cd(8);
   TPad *p25 = (TPad *) c1->cd(7);
   TPad *p26 = (TPad *) c1->cd(8);
   p1->SetLogy();
   p2->SetLogy();
   p4->SetLogy();
   p5->SetLogy();
   p7->SetLogy();
   p8->SetLogy();
   p10->SetLogy();
   p11->SetLogy();
   p13->SetLogy();
   p14->SetLogy();
   p16->SetLogy();
   p17->SetLogy();
   p19->SetLogy();
   p20->SetLogy();
   p22->SetLogy();
   p23->SetLogy();
   p25->SetLogy();
   p26->SetLogy();

   if (activateGATE) {
      TH1F    *hGATE = new TH1F("hGATE", "All protons in GATE QGSP-BIC-EMY (LT) dataset", nbinsz, zfrom, zto);
      TH1F    *hGATEStop = new TH1F("hGATEStop", "Protons in GATE QGSP-BIC-EMY (LT) dataset with processName == ProtonInelastic", nbinsz, zfrom, zto);
      TH2F    *hGATELateral = new TH2F("hGATELateral", "2D distribution of BP position in GATE QGSP-BIC-EMY (LT) dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);
      TH1F    *hGATEDefaultThreshold = new TH1F("hGATEDefaultThreshold", "All protons in GATE QGSP-BIC-EMY (DT) dataset", nbinsz, zfrom, zto);
      TH1F    *hGATEDefaultThresholdStop = new TH1F("hGATEDefaultThresholdStop", "Protons in GATE QGSP-BIC-EMY (DT) dataset with processName == ProtonInelastic", nbinsz, zfrom, zto);
      TH2F    *hGATEDefaultThresholdLateral = new TH2F("hGATEDefaultThresholdLateral", "2D distribution of BP position in GATE QGSP-BIC-EMY (DT) dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);
      TH1F    *hGATESingleScatter = new TH1F("hGATESingleScatter", "All protons in GATE QGSP-BIC-EMY (SS) dataset", nbinsz, zfrom, zto);
      TH1F    *hGATESingleScatterStop = new TH1F("hGATESingleScatterStop", "Protons in GATE QGSP-BIC-EMY (SS) dataset with processName == ProtonInelastic", nbinsz, zfrom, zto);
      TH2F    *hGATESingleScatterLateral = new TH2F("hGATESingleScatterLateral", "2D distribution of BP position in GATE QGSP-BIC-EMY (SS) dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);

      Float_t  x,y,z,edep;
      Int_t    parentID, eventID;
      Bool_t   isInelastic;
      
      TFile   *f1 = new TFile("Data/GATE/compressed_complex_geometry_bic.root");
      TTree   *treeBic = (TTree*) f1->Get("treeOut");
      treeBic->SetBranchAddress("posX",&x);
      treeBic->SetBranchAddress("posY",&y);
      treeBic->SetBranchAddress("posZ",&z);
      treeBic->SetBranchAddress("edep",&edep);
      treeBic->SetBranchAddress("eventID",&eventID);
      treeBic->SetBranchAddress("parentID",&parentID);
      treeBic->SetBranchAddress("isInelastic",&isInelastic);

      cout << "Analysing GATE-BIC-EMY tree. Found " << treeBic->GetEntries() << " entries.\n";
      for (Int_t i=0, N = treeBic->GetEntries(); i<N; ++i) {
         treeBic->GetEntry(i);

         hGATE->Fill(z/10.);
         if (!isInelastic)     hGATELateral->Fill(x,y);
         else                  hGATEStop->Fill(z/10.);
      }
      
      TFile   *f2 = new TFile("Data/GATE/compressed_complex_geometry_emy_nolimit.root");
      TTree   *treeDefaultThreshold = (TTree*) f2->Get("treeOut");
      treeDefaultThreshold->SetBranchAddress("posX",&x);
      treeDefaultThreshold->SetBranchAddress("posY",&y);
      treeDefaultThreshold->SetBranchAddress("posZ",&z);
      treeDefaultThreshold->SetBranchAddress("edep",&edep);
      treeDefaultThreshold->SetBranchAddress("eventID",&eventID);
      treeDefaultThreshold->SetBranchAddress("parentID",&parentID);
      treeDefaultThreshold->SetBranchAddress("isInelastic",&isInelastic);

      cout << "Analysing GATE-DefaultThreshold tree. Found " << treeDefaultThreshold->GetEntries() << " entries.\n";

      for (Int_t i=0, N = treeDefaultThreshold->GetEntries(); i<N; ++i) {
         treeDefaultThreshold->GetEntry(i);

         hGATEDefaultThreshold->Fill(z/10.);
         if (!isInelastic)    hGATEDefaultThresholdLateral->Fill(x,y);
         else                 hGATEDefaultThresholdStop->Fill(z/10.);
      }

      TFile   *f3 = new TFile("Data/GATE/compressed_complex_geometry_emy_nolimit_singlescatter_1k.root");
      TTree   *treeSingleScatter = (TTree*) f3->Get("treeOut");
      treeSingleScatter->SetBranchAddress("posX",&x);
      treeSingleScatter->SetBranchAddress("posY",&y);
      treeSingleScatter->SetBranchAddress("posZ",&z);
      treeSingleScatter->SetBranchAddress("edep",&edep);
      treeSingleScatter->SetBranchAddress("eventID",&eventID);
      treeSingleScatter->SetBranchAddress("parentID",&parentID);
      treeSingleScatter->SetBranchAddress("isInelastic",&isInelastic);

      cout << "Analysing GATE-SingleScatter tree. Found " << treeSingleScatter->GetEntries() << " entries.\n";

      for (Int_t i=0, N = treeSingleScatter->GetEntries(); i<N; ++i) {
         treeSingleScatter->GetEntry(i);

         hGATESingleScatter->Fill(z/10.);
         if (!isInelastic)    hGATESingleScatterLateral->Fill(x,y);
         else                 hGATESingleScatterStop->Fill(z/10.);
      }

      
      c1->cd(1);

      fitFunction = new TF1("fitFunction", "gaus", 20, 30);
      fitFunction->SetParameter(1, 11.61);
      fitFunction->SetParameter(2, 0.13);
      fitFunction->SetParLimits(1, 9, 13);
      fitFunction->SetParLimits(2, 0.05, 0.4);

      hGATE->SetXTitle("Range [cm]");
      hGATE->SetYTitle("Number of primaries");
      hGATE->SetFillColor(kBlue-7);
      hGATE->SetLineColor(kBlack);
      hGATE->Draw();
      hGATE->Fit("fitFunction", "B,M,Q");
      
      mu = fitFunction->GetParameter(1);
      sigma = fabs(fitFunction->GetParameter(2));
      
      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hGATE->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hGATE->Integral();
      cutoffIntegral = hGATE->Integral(0, distributionCutoffBin);
      Float_t BPIntegral = hGATE->Integral(distributionCutoffBin, hGATE->GetNbinsX());
      fraction = cutoffIntegral / float(totalIntegral);

      cout << "Cutoff at " << distributionCutoff << endl;
      cout << "Total integral = " << totalIntegral << endl;
      cout << "BP integral = " << BPIntegral << endl;
      cout << "cutoff integral = " << cutoffIntegral << endl;

      TLegend *legG = new TLegend(.18, .7, .68, .84);
      legG->AddEntry(hGATE, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      legG->AddEntry(fitFunction, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      legG->Draw();

      c1->cd(2);
      hGATEStop->SetXTitle("Range [cm]");
      hGATEStop->SetYTitle("Number of primaries");
      hGATEStop->SetFillColor(kBlue-7);
      hGATEStop->SetLineColor(kBlack);
      hGATEStop->Draw();

      c1->cd(3);
      hGATELateral->Draw("COLZ");

      ///

      TF1 *fitFunctionG2 = new TF1("fitFunctionG2", "gaus", 20, 30);
      fitFunctionG2->SetParameter(1, 11.61);
      fitFunctionG2->SetParameter(2, 0.13);
      fitFunctionG2->SetParLimits(1, 9, 13);
      fitFunctionG2->SetParLimits(2, 0.05, 0.4);

      c1->cd(4);
      hGATEDefaultThreshold->SetXTitle("Range [cm]");
      hGATEDefaultThreshold->SetYTitle("Number of primaries");
      hGATEDefaultThreshold->SetFillColor(kBlue-7);
      hGATEDefaultThreshold->SetLineColor(kBlack);
      hGATEDefaultThreshold->Draw();
      hGATEDefaultThreshold->Fit("fitFunctionG2", "B,M,Q");
      
      mu = fitFunctionG2->GetParameter(1);
      sigma = fabs(fitFunctionG2->GetParameter(2));
      
      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hGATEDefaultThreshold->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hGATEDefaultThreshold->Integral();
      cutoffIntegral = hGATEDefaultThreshold->Integral(1, distributionCutoffBin);
      fraction = cutoffIntegral / float(totalIntegral);

      cout << "Cutoff at " << distributionCutoff << endl;
      cout << "Total integral = " << totalIntegral << endl;
      cout << "cutoff integral = " << cutoffIntegral << endl;

      TLegend *legG2 = new TLegend(.18, .7, .68, .84);
      legG2->AddEntry(hGATEDefaultThreshold, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      legG2->AddEntry(fitFunctionG2, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      legG2->Draw();

      c1->cd(5);
      hGATEDefaultThresholdStop->SetXTitle("Range [cm]");
      hGATEDefaultThresholdStop->SetYTitle("Number of primaries");
      hGATEDefaultThresholdStop->SetFillColor(kBlue-7);
      hGATEDefaultThresholdStop->SetLineColor(kBlack);
      hGATEDefaultThresholdStop->Draw();

      c1->cd(6);
      hGATEDefaultThresholdLateral->Draw("COLZ");

      fitFunctionG2->SetParameter(1, 11.61);
      fitFunctionG2->SetParameter(2, 0.13);
      fitFunctionG2->SetParLimits(1, 9, 13);
      fitFunctionG2->SetParLimits(2, 0.05, 0.4);

      c1->cd(7);
      hGATESingleScatter->SetXTitle("Range [cm]");
      hGATESingleScatter->SetYTitle("Number of primaries");
      hGATESingleScatter->SetFillColor(kBlue-7);
      hGATESingleScatter->SetLineColor(kBlack);
      hGATESingleScatter->Draw();
      hGATESingleScatter->Fit("fitFunctionG2", "B,M,Q");
      
      mu = fitFunctionG2->GetParameter(1);
      sigma = fabs(fitFunctionG2->GetParameter(2));
      
      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hGATESingleScatter->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hGATESingleScatter->Integral();
      cutoffIntegral = hGATESingleScatter->Integral(1, distributionCutoffBin);
      fraction = cutoffIntegral / float(totalIntegral);

      cout << "Cutoff at " << distributionCutoff << endl;
      cout << "Total integral = " << totalIntegral << endl;
      cout << "cutoff integral = " << cutoffIntegral << endl;

      TLegend *legG3 = new TLegend(.18, .7, .68, .84);
      legG3->AddEntry(hGATESingleScatter, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      legG3->AddEntry(fitFunctionG2, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      legG3->Draw();

      c1->cd(8);
      hGATESingleScatterStop->SetXTitle("Range [cm]");
      hGATESingleScatterStop->SetYTitle("Number of primaries");
      hGATESingleScatterStop->SetFillColor(kBlue-7);
      hGATESingleScatterStop->SetLineColor(kBlack);
      hGATESingleScatterStop->Draw();

      c1->cd(9);
      hGATESingleScatterLateral->Draw("COLZ");
      
      TH1D *hProfile = hGATELateral->ProjectionX();
      int bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      int bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      double fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM GATE Low Threshold 2D = " << fwhm << " mm." << endl;
      hProfile = hGATEDefaultThresholdLateral->ProjectionX();
      bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM GATE Default Threshold 2D = " << fwhm << " mm." << endl;
      hProfile = hGATESingleScatterLateral->ProjectionX();
      bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM GATE Single Scatter 2D = " << fwhm << " mm." << endl;
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
      
      TH1F    *hMCNPCEM = new TH1F("hMCNPCEM", "All protons in MCNP CEM dataset;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH1F    *hMCNPCEMStop = new TH1F("hMCNPCEMStop", "Protons in MCNP CEM dataset with termination type 13 or 16;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH2F    *hMCNPCEMLateral = new TH2F("hMCNPLateral", "2D distribution of BP position in MCNP CEM dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);
      TH1F    *hMCNPBERT = new TH1F("hMCNPBERT", "All protons in MCNP BERT dataset;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH1F    *hMCNPBERTStop = new TH1F("hMCNPBERTStop", "Protons in MCNP BERT dataset with termination type 13 or 16;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH2F    *hMCNPBERTLateral = new TH2F("hMCNPBERTLateral", "2D distribution of BP position in MCNP BERT dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);
      TH1F    *hMCNPINC = new TH1F("hMCNPINC", "All protons in MCNP INC dataset;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH1F    *hMCNPINCStop = new TH1F("hMCNPINCStop", "Protons in MCNP INC dataset with termination type 13 or 16;Range [cm];Number of primaries", nbinsz, zfrom, zto);
      TH2F    *hMCNPINCLateral = new TH2F("hMCNPINCLateral", "2D distribution of BP position in MCNP INC dataset;X position [mm];Y position [mm]", nbinsxy, xyfrom, xyto, nbinsxy, xyfrom, xyto);


      cout << "READING CEM FILE\n";
      in.open("Data/MCNP/190MeV_sdet_p_cem");

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

               hMCNPCEM->Fill(z - startZ);
               if (terminationType == 13 || terminationType == 16) {
                  hMCNPCEMStop->Fill(z - startZ);
               }
               else {
                  hMCNPCEMLateral->Fill(10*x,10*y);
               }
            }
         }
      }

      in.close();
      
      cout << "READING BERT FILE\n";
      
      in2.open("Data/MCNP/190MeV_sdet_p_bert");

      while (! in2.eof() ) {
         getline(in2, line);

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
               in2 >> x >> y >> z >> u >> v >> w >> energy >> weight >> time;

               hMCNPBERT->Fill(z - startZ);
               if (terminationType == 13 || terminationType == 16) {
                  hMCNPBERTStop->Fill(z - startZ);
               }
               else {
                  hMCNPBERTLateral->Fill(10*x,10*y);
               }
            }
         }
      }
      in2.close();

      cout << "READING INC FILE\n";
      in3.open("Data/MCNP/190MeV_sdet_p_inc");

      while (! in3.eof() ) {
         getline(in3, line);

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
               in3 >> x >> y >> z >> u >> v >> w >> energy >> weight >> time;

               hMCNPINC->Fill(z - startZ);
               if (terminationType == 13 || terminationType == 16) {
                  hMCNPINCStop->Fill(z - startZ);
               }
               else {
                  hMCNPINCLateral->Fill(10*x,10*y);
               }
            }
         }
      }
      in3.close();

      cout << "PLOTTING\n";

      c2->cd(1);
      hMCNPCEM->SetFillColor(kGreen-3);
      hMCNPCEM->Draw();

      fitFunction = new TF1("fitFunction", "gaus");
      fitFunction->SetParameter(1, 11.5);
      fitFunction->SetParameter(2, 0.15);
      fitFunction->SetParLimits(1, 9, 13);
      fitFunction->SetParLimits(2, 0.05, 0.4);

      hMCNPCEM->Fit("fitFunction", "B, M,Q");

      mu = fitFunction->GetParameter(1);
      sigma = fabs(fitFunction->GetParameter(2));

      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hMCNPCEM->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hMCNPCEM->Integral();
      cutoffIntegral = hMCNPCEM->Integral(0, distributionCutoffBin);
      stoppedIntegral = hMCNPCEMStop->Integral();
      fraction = cutoffIntegral / float(totalIntegral);

      TLegend *leg = new TLegend(.18, .7, .68, .84);
      leg->AddEntry(hMCNPCEM, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      leg->AddEntry(fitFunction, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      leg->Draw();

      cout << Form("Total number of events: %d. Number of events with terminationType == (13, 16): %.0f. Fraction: %.2f.", totalIntegral, hMCNPCEMStop->Integral(), fraction) << endl;
      
      c2->cd(2);
      hMCNPCEMStop->SetFillColor(kGreen-3);
      hMCNPCEMStop->Draw();

      c2->cd(3);
      hMCNPCEMLateral->Draw("COLZ");
////////////////
      c2->cd(4);
      hMCNPBERT->SetFillColor(kGreen-3);
      hMCNPBERT->Draw();

      TF1 * fitFunction = new TF1("fitFunction", "gaus");
      fitFunction->SetParameter(1, 11.5);
      fitFunction->SetParameter(2, 0.15);
      fitFunction->SetParLimits(1, 9, 13);
      fitFunction->SetParLimits(2, 0.05, 0.4);

      hMCNPBERT->Fit("fitFunction", "B, M,Q");

      mu = fitFunction->GetParameter(1);
      sigma = fabs(fitFunction->GetParameter(2));

      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hMCNPBERT->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hMCNPBERT->Integral();
      cutoffIntegral = hMCNPBERT->Integral(0, distributionCutoffBin);
      stoppedIntegral = hMCNPBERTStop->Integral();
      fraction = cutoffIntegral / float(totalIntegral);

      TLegend *leg2 = new TLegend(.18, .7, .68, .84);
      leg2->AddEntry(hMCNPBERT, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      leg2->AddEntry(fitFunction, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      leg2->Draw();

      cout << Form("Total number of events: %d. Number of events with terminationType == (13, 16): %.0f. Fraction: %.2f.", totalIntegral, hMCNPBERTStop->Integral(), fraction) << endl;
      
      c2->cd(5);
      hMCNPBERTStop->SetFillColor(kGreen-3);
      hMCNPBERTStop->Draw();

      c2->cd(6);
      hMCNPBERTLateral->Draw("COLZ");
//////////////
      c2->cd(7);
      hMCNPINC->SetFillColor(kGreen-3);
      hMCNPINC->Draw();

      TF1 * fitFunction3 = new TF1("fitFunction3", "gaus");
      fitFunction3->SetParameter(1, 11.5);
      fitFunction3->SetParameter(2, 0.15);
      fitFunction3->SetParLimits(1, 9, 13);
      fitFunction3->SetParLimits(2, 0.05, 0.4);

      hMCNPINC->Fit("fitFunction3", "B, M,Q");

      mu = fitFunction3->GetParameter(1);
      sigma = fabs(fitFunction3->GetParameter(2));

      distributionCutoff = mu - 3*sigma;
      distributionCutoffBin = hMCNPINC->GetXaxis()->FindBin(distributionCutoff);
      totalIntegral = hMCNPINC->Integral();
      cutoffIntegral = hMCNPINC->Integral(0, distributionCutoffBin);
      stoppedIntegral = hMCNPINCStop->Integral();
      fraction = cutoffIntegral / float(totalIntegral);

      TLegend *leg3 = new TLegend(.18, .7, .68, .84);
      leg3->AddEntry(hMCNPINC, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
      leg3->AddEntry(fitFunction3, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
      leg3->Draw();

      cout << Form("Total number of events: %d. Number of events with terminationType == (13, 16): %.0f. Fraction: %.2f.", totalIntegral, hMCNPINCStop->Integral(), fraction) << endl;
      
      c2->cd(8);
      hMCNPINCStop->SetFillColor(kGreen-3);
      hMCNPINCStop->Draw();

      c2->cd(9);
      hMCNPINCLateral->Draw("COLZ");
      
      TH1D *hProfile = hMCNPINCLateral->ProjectionX();
      int bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      int bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      double fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM MCNP INC 2D = " << fwhm << " mm." << endl;
      hProfile = hMCNPBERTLateral->ProjectionX();
      bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM MCNP BERT 2D = " << fwhm << " mm." << endl;
      hProfile = hMCNPCEMLateral->ProjectionX();
      bin1 = hProfile->FindFirstBinAbove(hProfile->GetMaximum()/2);
      bin2 = hProfile->FindLastBinAbove(hProfile->GetMaximum()/2);
      fwhm = hProfile->GetBinCenter(bin2) - hProfile->GetBinCenter(bin1);
      cout << "FWHM MCNP CEM 2D = " << fwhm << " mm." << endl;
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
}
