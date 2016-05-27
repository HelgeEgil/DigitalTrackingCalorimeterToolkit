#define findRange_cxx
#include "include/findRange.h"
#include "../GlobalConstants/Constants.h"
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

void findRange::Loop(Double_t energy, Double_t sigma_mev)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
   TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
   TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
   TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
   TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);

   Int_t nbinsx = 20000;
   Int_t xfrom = 0;
   Int_t xto = 50;

	Float_t x_compensate = 0;

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
	
	Float_t tl = 0;
	Int_t n = 0;
	
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

//	if (posY<0) continue;

	if (parentID == 0) {
	
		hStepLength->Fill(stepLength);

		Float_t z = posZ;
		Float_t y = posY;
		Float_t x = posX;
		
		
		if (fabs(x) < 20 && fabs(y) < 20 && volumeID[4] == 4) hZ->Fill(z + x_compensate, edep);
		n++;
			
		if (eventID != lastID) {
			n = 0;
			
			Float_t diff = sqrt( pow(firstX - lastX, 2) + pow(firstY - lastY, 2) + pow(0 - lastZ, 2));
			hRange->Fill(lastRange);
			hTracklength->Fill(diff);
			hActualTracklength->Fill(tl + firstZ);

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

	c2->cd();
	
	Float_t fit_tl = 0, fit_range = 0;

	TF1 *fRange = new TF1("fit_range", "gaus", xfrom, xto);
//	fRange->SetParameters(100, 157, 0.5);
	hRange->Fit("fit_range", "Q,M,WW", "", xfrom, xto);
 	fit_range = fRange->GetParameter(1);
	cout << Form("Range: %.3f mm +- %.3f mm.\n", fRange->GetParameter(1), fRange->GetParameter(2));
	
 	hTracklength->Fit("fit_range", "Q,M,WW", "", xfrom, xto);
 	cout << Form("Straight line: %.3f mm +- %.3f mm.\n", fRange->GetParameter(1), fRange->GetParameter(2));
	
 	hActualTracklength->Fit("fit_range", "Q,M,WW", "", xfrom, xto);
 	fit_tl = fRange->GetParameter(1);
 	cout << Form("Tracklength: %.3f mm +- %.3f mm.\n", fit_tl, fRange->GetParameter(2));
	
 	cout << "Detour factor = " << fit_range / fit_tl << endl;
 	Float_t straggling = sqrt( alphaprime * (pow(p, 2) * pow(alpha, 2/p) / (3-2/p) * pow(fit_range, 3-2/p)) );
 	cout << "Expected Energy from Bortfeld = " << getEnergyFromTLQuadratic(fit_range) << endl;
 	cout << "Expected straggling from Bortfeld = " << straggling << endl;

 	Float_t cutoff = fit_range - 4*straggling;
 	Float_t total = hRange->Integral();
 	Float_t attenuation = hRange->Integral(0, hRange->GetXaxis()->FindBin(cutoff));

 	cout << "Number of protons attenuated (more than 4 sigma below) = " << 100 * attenuation / total << " %.\n";
	
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
		
	c5->cd();
		hStepLength->SetXTitle("Steplength [mm]");
		hStepLength->SetYTitle("Number of steps");
		hStepLength->SetFillColor(kBlue-7);
		hStepLength->SetLineColor(kBlack);
		hStepLength->Draw();

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
