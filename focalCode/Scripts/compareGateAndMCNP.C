#define compareGateAndMCNP_cxx
#include "compareGateAndMCNP.h"
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

using namespace std;

void compareGateAndMCNP::Run()
{
   // GATE PART

   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
   c1->Divide(2,2,0.001,0.001);
   c1->SetLogy();

   // Set logy in all pads
   TPad *p1 = (TPad *) c1->cd(1);
   TPad *p2 = (TPad *) c1->cd(2);
   TPad *p3 = (TPad *) c1->cd(3);
   TPad *p4 = (TPad *) c1->cd(4);
   p1->SetLogy();
   p2->SetLogy();
   p3->SetLogy();
   p4->SetLogy();

   Int_t nbinsx = 500;
   Int_t xfrom = 0;
   Int_t xto = 27;

	TH1F *hGATE = new TH1F("hGATE", "All protons in GATE dataset", nbinsx, xfrom, xto);
	TH1F *hGATEStop = new TH1F("hGATEStop", "Protons in GATE dataset with processName == ProtonInelastic", nbinsx, xfrom, xto);

   Int_t    lastEvent = -1;
	Float_t  lastRange = 0;
	Int_t    lastID = -1;

   // GATE READOUT

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
	
   Long64_t ientry = LoadTree(0);
	fChain->GetEntry(0);
	
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	
	if (ientry < 0) {
		break;
	}

	nb = fChain->GetEntry(jentry);   nbytes += nb;


	if (lastID < 0) lastID = eventID;

	if (parentID == 0) {

		if (eventID != lastID) {
			hGATE->Fill(lastRange/10.);
   	}

      else {
         if (processName[0] == 'P') hGATEStop->Fill(lastRange/10.);
      }

		lastRange = posZ;
      lastID = eventID;
		}
   }
   
   TF1 *fitFunctionG = new TF1("fitFunctionG", "gaus", 20, 30);
   fitFunctionG->SetParameter(1, 24);
   fitFunctionG->SetParameter(2, 0.5);

   c1->cd(1);
   hGATE->SetXTitle("Range [cm]");
   hGATE->SetYTitle("Number of primaries");
   hGATE->SetFillColor(kBlue-7);
   hGATE->SetLineColor(kBlack);
   hGATE->Draw();
   hGATE->Fit("fitFunctionG", "B,M");
   
   Float_t muG = fitFunctionG->GetParameter(1);
   Float_t sigmaG = fabs(fitFunctionG->GetParameter(2));
   
   Float_t  distributionCutoffG = muG - 3*sigmaG;
   Int_t    distributionCutoffBinG = hGATE->GetXaxis()->FindBin(distributionCutoffG);
   Int_t    totalIntegralG = hGATE->Integral();
   Int_t    cutoffIntegralG = hGATE->Integral(0, distributionCutoffBinG);
   Float_t  fractionG = cutoffIntegralG / float(totalIntegralG);

   TLegend *legG = new TLegend(.18, .7, .68, .84);
   legG->AddEntry(hGATE, Form("Fraction below 3#sigma = %.2f%%", 100*fractionG), "F");
   legG->AddEntry(fitFunctionG, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", muG, sigmaG), "L");
   legG->Draw();

   c1->cd(2);
   hGATEStop->SetXTitle("Range [cm]");
   hGATEStop->SetYTitle("Number of primaries");
   hGATEStop->SetFillColor(kBlue-7);
   hGATEStop->SetLineColor(kBlack);
   hGATEStop->Draw();

   
   //
   // MCNP PART
   // 

   ifstream in;
   Int_t    terminationType;
   Int_t    EOLIdentifier;
   Int_t    historyType;
   Int_t    branchNumber;
   Float_t  startZ = -40;
   string   line;
   Float_t  x, y, z;
   Float_t  u, v, w; // vectors
   Float_t  energy, weight, time;
   
   TH1F    *hMCNP = new TH1F("hMCNP", "All protons in MCNP dataset;Range [cm];Number of primaries", nbinsx, xfrom, xto);
   TH1F    *hMCNPStop = new TH1F("hMCNPStop", "Protons in MCNP dataset with termination type 13 or 16;Range [cm];Number of primaries", nbinsx, xfrom, xto);

   in.open("hp_proton_190MeVr");

   while (! in.eof() ) {
      getline(in, line);

      if (line.size() < 11) continue;
      historyType = atoi(line.substr(6,5).c_str());

      if (historyType == 5000 || historyType == 9000) {
      
         EOLIdentifier = atoi(line.substr(15,7).c_str());
         if (EOLIdentifier == 5000) continue;

         branchNumber = atoi(line.substr(40,3).c_str());
         terminationType = atoi(line.substr(29,2).c_str());

         if (branchNumber == 1) { // primary particle
            in >> x >> y >> z >> u >> v >> w >> energy >> weight >> time;

            hMCNP->Fill(z - startZ);
            if (terminationType == 13 || terminationType == 16) {
               hMCNPStop->Fill(z - startZ);
            }
         }
      }
   }
   in.close();

   c1->cd(3);
   hMCNP->SetFillColor(kGreen-3);
   hMCNP->Draw();

   TF1 * fitFunction = new TF1("fitFunction", "gaus", 20, 30);
   fitFunction->SetParameter(1, 24);
   fitFunction->SetParameter(2, 0.5);

   hMCNP->Fit("fitFunction", "V, B, M");

   Float_t mu = fitFunction->GetParameter(1);
   Float_t sigma = fabs(fitFunction->GetParameter(2));

   Float_t  distributionCutoff = mu - 3*sigma;
   Int_t    distributionCutoffBin = hMCNP->GetXaxis()->FindBin(distributionCutoff);
   Int_t    totalIntegral = hMCNP->Integral();
   Int_t    cutoffIntegral = hMCNP->Integral(0, distributionCutoffBin);
   Int_t    stoppedIntegral = hMCNPStop->Integral();
   Float_t  fraction = cutoffIntegral / float(totalIntegral);

   TLegend *leg = new TLegend(.18, .7, .68, .84);
   leg->AddEntry(hMCNP, Form("Fraction below 3#sigma = %.2f%%", 100*fraction), "F");
   leg->AddEntry(fitFunction, Form("(#mu, #sigma) = (%.2f cm, %.2f cm)", mu, sigma), "L");
   leg->Draw();

   cout << Form("Total number of events: %d. Number of events with terminationType == (13, 16): %.0f. Fraction: %.2f.", totalIntegral, hMCNPStop->Integral(), fraction) << endl;
   
   c1->cd(4);
   hMCNPStop->SetFillColor(kGreen-3);
   hMCNPStop->Draw();
}
