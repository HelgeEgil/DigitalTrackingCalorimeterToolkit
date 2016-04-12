#define findRange_cxx
#include "findRange.h"
#include "Constants.h"
#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

Float_t p= 1.658973;
Float_t alpha = 0.004663;

Float_t getEnergyFromTL(Float_t tl) {
	return pow(tl / alpha, 1/p);
}

Float_t getTLFromEnergy(Float_t energy) {
	return alpha * pow(energy, p);
}

void findRange::Loop(Double_t energy, Double_t sigma_mev)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
   TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
   TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
   TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
//    TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
//    TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);
//    TCanvas *c7 = new TCanvas("c7", "c7", 800, 600);
//    TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);

   Int_t nbinsx = 10000;
   Int_t xfrom = 0;
   Int_t xto = 50;

	Float_t x_compensate = 0;

	TH1F *hZ = new TH1F("hZ", "Z profile", nbinsx/3, xfrom + x_compensate, xto + x_compensate);
	TH1F *hRange = new TH1F("hRange", "Primary ranges", nbinsx, xfrom + x_compensate, xto + x_compensate);
	TH1F *hTracklength = new TH1F("hTracklength", "Straight tracklengths", nbinsx, xfrom + x_compensate, xto + x_compensate);
	TH1F *hActualTracklength = new TH1F("hActualTracklength", "Actual tracklengths", nbinsx, xfrom + x_compensate, xto + x_compensate);
// 	TH1F *hRangec = new TH1F("hRangec", "Corrected primary ranges", nbinsx, xfrom + x_compensate, xto + x_compensate);
// 	TH2F *hRange2D = new TH2F("hRange2D", "Proton ranges in scintillators + focal", 300, 18, 30, 128, -25, 25);
// 	TH2F *hEnergy2D = new TH2F("hEnergy2D", "Estimated initial energy", 300, 150, 200, 128, -25, 25);
// 	TH2F *hEnergy2Dc = new TH2F("hEnergy2Dc", "Estimated initial energy, corrected for scintillators", 300, 150, 200, 128, -25, 25);
// 	TH2F *hEnergy2Dcg = new TH2F("hEnergy2Dcg", "Estimated initial energy, Gauss-corrected for scintillators", 300, 150, 200, 128, -25, 25);
// 	TH2F *hRange2Dtrue = new TH2F("hRange2Dtrue", "True proton ranges in scintillators + focal", 300, 18, 30, 128, -25, 25);
// 	TH2F *hEnergy2Dtrue = new TH2F("hEnergy2Dtrue", "True Estimated initial energy", 300, 150, 200, 128, -25, 25);
	

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
	
	Float_t tl = 0;
	
	TRandom3 *gRandom = new TRandom3();

	Long64_t ientry = LoadTree(0);
	fChain->GetEntry(0);
	
	firstX = posX;
	firstY = posY;
	firstZ = posZ;
	
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	
	if (ientry < 0) {
		cout << "Aborting run at jentry = " << jentry << endl;
		break;
	}

	nb = fChain->GetEntry(jentry);   nbytes += nb;


	if (lastID < 0) lastID = eventID;


	if (parentID == 0) {
	
		Float_t z = posZ;
		Float_t y = posY;
		Float_t x = posX;
		
		if (fabs(x) < 20 && fabs(y) < 20) hZ->Fill(z + x_compensate, edep);
			
		if (eventID != lastID) {
			Float_t diff = sqrt( pow(firstX - lastX, 2) + pow(firstY - lastY, 2) + pow(firstZ - lastZ, 2));
			hRange->Fill(lastRange - firstZ);
			hTracklength->Fill(diff);
			hActualTracklength->Fill(tl);
			
	// 				hRange2D->Fill(lastRange + x_compensate, lastX);
	// 				hEnergy2D->Fill(getEnergyFromTL(lastRange + x_compensate), lastX);
	// 				
	// 				if (abs(lastX)>5 && abs(lastY)>5) {
	// 					dE = 6.65;
	// 					dE_random = gRandom->Gaus(dE, 0.35);
	// 				}
	// 				else if (abs(lastX)<5 && abs(lastY)>5) {
	// 					dE = 9.23;
	// 					dE_random = gRandom->Gaus(dE, 0.43);
	// 				}
	// 				else if (abs(lastX)>5 && abs(lastY)<5) {
	// 					dE = 9.23;
	// 					dE_random = gRandom->Gaus(dE, 0.43);
	// 				}
	// 				else if (abs(lastX)<5 && abs(lastY)<5) {
	// 					dE = 11.79;
	// 					dE_random = gRandom->Gaus(dE, 0.55);
	// 				}
	// 				
	// 				dTL = getTLFromEnergy(energy) - getTLFromEnergy(energy - dE);
	// 				hEnergy2Dc->Fill(getEnergyFromTL(lastRange + x_compensate) + dE, lastX);
	// 				hEnergy2Dcg->Fill(getEnergyFromTL(lastRange + x_compensate) + dE_random, lastX);
	// 				hRangec->Fill(lastRange + x_compensate + dTL);
					
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

	// range 1: 80 % of maximum for bragg peak on distal edge
	Double_t range_1 = hZ->GetXaxis()->GetBinCenter(hZ->FindLastBinAbove(hZ->GetMaximum() * 0.8));
	Double_t range_3 = hRange->GetXaxis()->GetBinCenter(hRange->GetMaximumBin());
	
 	cout << "Maximum from bragg peak plot, 80\% of maximum on distal edge: " << range_1 << " mm.\n";
// 	cout << "Maximum from hRange plot: Max bin " << range_3 << " mm.\n";

	c2->cd();
	
	TF1 *fRange = new TF1("fit_range", "gaus", 0, 900);
	hRange->Fit("fit_range", "Q, W", "", 0, 900);
	cout << Form("Mean range: %.2f mm +- %.2f mm.\n", fRange->GetParameter(1), fRange->GetParameter(2));
	
 	hTracklength->Fit("fit_range", "Q, W", "", 0, 900);
 	cout << Form("Mean tracklength: %.2f mm +- %.2f mm.\n", fRange->GetParameter(1), fRange->GetParameter(2));
	
 	hActualTracklength->Fit("fit_range", "Q, W", "", 0, 900);
 	cout << Form("Mean actual tracklengths: %.2f mm +- %.2f mm.\n", fRange->GetParameter(1), fRange->GetParameter(2));
	
	
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
		
// 		hRange2D->SetYTitle("X position [mm]");
// 		hRange2D->SetXTitle("Range [mm]");
// 		hRange2D->SetFillColor(kBlue-7);
// 		hRange2D->SetLineColor(kBlack);
// 		hRange2D->Draw("COLZ");
		
// 	c4->cd();
// 		hEnergy2D->SetYTitle("X position [mm]");
// 		hEnergy2D->SetXTitle("Energy [MeV]");
// 		hEnergy2D->SetFillColor(kBlue-7);
// 		hEnergy2D->SetLineColor(kBlack);
// 		hEnergy2D->Draw("COLZ");
// 		
// 	c5->cd();
// 		hEnergy2Dc->SetYTitle("X position [mm]");
// 		hEnergy2Dc->SetXTitle("Corrected energy [MeV]");
// 		hEnergy2Dc->SetFillColor(kBlue-7);
// 		hEnergy2Dc->SetLineColor(kBlack);
// 		hEnergy2Dc->Draw("COLZ");
// 		
// 	c6->cd();
// 		hEnergy2Dcg->SetYTitle("X position [mm]");
// 		hEnergy2Dcg->SetXTitle("Corrected energy w/Gauss [MeV]");
// 		hEnergy2Dcg->SetFillColor(kBlue-7);
// 		hEnergy2Dcg->SetLineColor(kBlack);
// 		hEnergy2Dcg->Draw("COLZ");
// 
// 	c7->cd();
// 		hRangec->SetXTitle("Range [mm]");
// 		hRangec->SetYTitle("Number of primaries");
// 		hRangec->SetFillColor(kBlue-7);
// 		hRangec->SetLineColor(kBlack);
// 		hRangec->Draw();
		
}