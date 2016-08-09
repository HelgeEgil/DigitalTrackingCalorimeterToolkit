#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TText.h>
#include <iostream>
#include <fstream>
#include <TLatex.h>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

void makePlots() {
   TCanvas *c1	= new TCanvas("c1", "c1", 1200, 800);
	TCanvas *c2 = new TCanvas("c2", "Correct Tracks fraction", 1200, 800);
	TCanvas *c3 = new TCanvas("c3", "Reconstruction efficiency", 1200, 800);

	Float_t arrayE[200] = {0}; // energy MC
	Float_t arrayEE[200] = {0}; // error on energy MC
	Float_t arrayMC[200] = {0}; // range MC
	Float_t arrayEMC[200] = {0}; // error on range MC
	Float_t arrayPSTAR[200] = {0};
	Float_t arrayEPstar[200] = {0};
	Float_t arrayE2[200] = {0}; // energy data 
	Float_t arrayEE2[200] = {0}; // error on energy data
	Float_t arrayEData[200] = {0};  // range data
	Float_t arrayData[200] = {0}; // error on range data
	Float_t arrayEfficiencyEnergy[300];
	Float_t arrayEfficiencyFinal[300];
	Float_t arrayEfficiencyFinal2[300];

	Int_t	nThisEnergy = 0, lastEnergy = 0;
	
   gStyle->SetOptStat(0);

	ifstream in;
	in.open("OutputFiles/result_makebraggpeakfit.csv");

	cout << "Opened file.\n";

	Float_t nomrange_, estrange_, sigmaRange_, lastRange_;
	Int_t energy_;

	Int_t nlines = 0;
	TNtuple *ntuple = new TNtuple("ntuple", "data from file", "energy_:nomrange_:estrange_:sigmaRange:lastRange_");
	
	Int_t MC2Data = 61;

	while (1) {
		in >> energy_ >> nomrange_ >> estrange_ >> sigmaRange_ >> lastRange_ ;

		if (!in.good()) {
			cout << energy_ << ", " << nomrange_ << ", " << estrange_ << ", " << sigmaRange_ << ", " << lastRange_ << endl;
			break;
		}

		if (nlines < MC2Data) {
			cout << "Line " << nlines << ", energy " << energy_ << ",  MC" << endl;
			arrayE[nlines] = energy_;
			arrayEE[nlines] = 0;
			arrayMC[nlines] = estrange_;
			arrayEMC[nlines] = sigmaRange_;
			arrayPSTAR[nlines] = nomrange_;
			arrayEPstar[nlines] = 0;
		}

		else {
			cout << "Line " << nlines << ", (or " << nlines-MC2Data << "), energy " << energy_ << ", data" << endl;
			arrayE2[nlines-MC2Data] = energy_;
			arrayEE2[nlines-MC2Data] = 0;
			arrayData[nlines-MC2Data] = estrange_;
			arrayEData[nlines-MC2Data] = sigmaRange_;
		}

		nlines++;
	}

	in.close();
	
	ifstream in2;
	in2.open("OutputFiles/lastLayerCorrect_different_nRuns.csv");
	Float_t factor, np, correctLast, correctWhole, lastIsFirst, lastIsAlmostFirst;
	Float_t arrayFractionX[200] = {0};
	Float_t arrayFractionY[200] = {0};
	Float_t arrayFractionY2[200] = {0};
	Float_t arrayFractionY3[200] = {0};
	Int_t nlines2 = 0;
	while (1) {
		in2 >> factor >> np >> correctWhole >> lastIsFirst >> lastIsAlmostFirst;

		if (!in2.good()) break;
		
		arrayFractionX[nlines2] = np;
		arrayFractionY[nlines2] = correctWhole * 100;
		arrayFractionY2[nlines2] = lastIsFirst * 100;
		arrayFractionY3[nlines2] = lastIsAlmostFirst * 100;

		nlines2++;
	}
	cout << "Found " << nlines2 << " lines in lastLayerCorrect.\n";
	
	in2.close();
	
	ifstream in3;
	in3.open("OutputFiles/efficiency_500.csv");
	Int_t energy, nRecon, nNotLeaving, nFinal;
	Int_t nEnergies = 0;

	while (1) {
		in3 >> energy >> np >> nRecon >> nNotLeaving >> nFinal;

		if (!in3.good()) break;
		if (!lastEnergy) {
			arrayEfficiencyEnergy[0] = energy;
			arrayEfficiencyFinal[0] = nFinal;
			nThisEnergy = np;
			lastEnergy = energy;
		}

		else if (lastEnergy == energy) {
			arrayEfficiencyFinal[nEnergies] += nFinal;
			nThisEnergy += np;
			lastEnergy = energy;
		}

		else if (lastEnergy != energy) {
			arrayEfficiencyFinal[nEnergies] /= nThisEnergy;
			nEnergies++;
			arrayEfficiencyEnergy[nEnergies] = energy;
			arrayEfficiencyFinal[nEnergies] = nFinal;
			nThisEnergy = np;
			lastEnergy = energy;
		}
	}

	in3.close();

	ifstream in4;
	in4.open("OutputFiles/efficiency_200.csv");
	nEnergies = 0;
	lastEnergy = 0;
	while (1) {
		in4 >> energy >> np >> nRecon >> nNotLeaving >> nFinal;

		if (!in4.good()) break;
		if (!lastEnergy) {
			arrayEfficiencyFinal2[0] = nFinal;
			nThisEnergy = np;
			lastEnergy = energy;
		}

		else if (lastEnergy == energy) {
			arrayEfficiencyFinal2[nEnergies] += nFinal;
			nThisEnergy += np;
			lastEnergy = energy;
		}

		else if (lastEnergy != energy) {
			arrayEfficiencyFinal2[nEnergies] /= nThisEnergy;
			nEnergies++;
			arrayEfficiencyFinal2[nEnergies] = nFinal;
			nThisEnergy = np;
			lastEnergy = energy;
		}
	}

	c1->cd();
	
	TGraphErrors *hMC = new TGraphErrors(MC2Data, arrayE, arrayMC, arrayEE, arrayEMC);
	TGraphErrors *hData = new TGraphErrors(nlines-MC2Data, arrayE2, arrayData, arrayEE2, arrayEData);
	TGraphErrors *pstar = new TGraphErrors(MC2Data, arrayE, arrayPSTAR, arrayEE, arrayEPstar);
	
	hMC->GetXaxis()->SetTitleFont(22);
	hMC->GetYaxis()->SetTitleFont(22);
	hMC->GetXaxis()->SetTitleOffset(1.2);
	hMC->GetYaxis()->SetTitleOffset(1.2);
	hMC->GetXaxis()->SetLabelFont(22);
	hMC->GetYaxis()->SetLabelFont(22);

	hMC->GetXaxis()->SetRangeUser(145, 200);
	hMC->GetYaxis()->SetRangeUser(145, 270);

	hMC->SetMarkerColor(kBlue);
	hMC->SetMarkerStyle(21);
	hMC->SetMarkerSize(1.5);

	hData->SetMarkerColor(kRed);
	hData->SetMarkerStyle(22);
	hData->SetMarkerSize(2);
	
	pstar->SetLineWidth(3);
	pstar->SetLineColor(kRed);

	gStyle->SetPadTickY(1);
	hMC->SetTitle("Range estimation of proton tracks using weighted Gaussian approach;Energy [MeV];Projected range [mm]");
	hMC->Draw("AP");
	hData->Draw("P");
	pstar->Draw("L");

	gPad->Update();
	TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
	title->SetTextFont(22);
	gPad->Modified();

	TLegend *leg = new TLegend(0.15, 0.68, 0.40, 0.85);
	leg->SetTextSize(0.035);
	leg->SetTextFont(22);
	leg->AddEntry(pstar, "PSTAR range", "L");
	leg->AddEntry(hMC, "Monte Carlo", "PE");
	leg->AddEntry(hData, "Experimental data", "PE");
	leg->Draw();

	c1->Update();

	c2->cd();
	TGraph *gFraction = new TGraph(nlines2, arrayFractionX, arrayFractionY);
	TGraph *gFraction2 = new TGraph(nlines2, arrayFractionX, arrayFractionY2);
	TGraph *gFraction3 = new TGraph(nlines2, arrayFractionX, arrayFractionY3);
	gFraction->GetXaxis()->SetRangeUser(10, 6000);
	gFraction->SetMaximum(100);
	gFraction->SetMinimum(0);
	gFraction->GetXaxis()->SetTitleFont(22);
	gFraction->GetYaxis()->SetTitleFont(22);
	gFraction->GetXaxis()->SetTitleOffset(1.2);
	gFraction->GetYaxis()->SetTitleOffset(1.2);
	gFraction->GetXaxis()->SetLabelFont(22);
	gFraction->GetYaxis()->SetLabelFont(22);
	gFraction->SetTitle("Fraction of correctly reconstructed tracks;Number of protons in frame;Fraction of correctly reconstructed tracks");
	gFraction->SetLineColor(kGreen+2);
	gFraction->SetLineWidth(3);
	gFraction2->SetLineColor(kAzure+4);
	gFraction2->SetLineWidth(3);
	gFraction3->SetLineColor(kPink+4);
	gFraction3->SetLineWidth(3);

	gFraction->Draw("AL");
	gFraction2->Draw("L");
	gFraction3->Draw("L");
	
	TText *t = new TText();
	gFraction->GetYaxis()->SetLabelOffset(5);
	t->SetTextAlign(32);
	t->SetTextSize(0.035);
	t->SetTextFont(22);
	for (Int_t i=0; i<6;i++) {
		cout << "Drawing text at " << -0.42 << ", " << i*20 << endl;
		t->DrawText(-0.42, i*20, Form("%d%%", i*20));
	}

	gPad->SetLogx();
	gFraction->GetXaxis()->SetNoExponent();

	gPad->Update();
	title = (TPaveText*) gPad->GetPrimitive("title");
	title->SetTextFont(22);
	gPad->Modified();
	
	TLegend * leg2 = new TLegend(0.67, 0.76, 0.97, 0.93);
	leg2->SetTextSize(0.025);
	leg2->SetTextFont(22);
	leg2->AddEntry(gFraction, "Whole track correct", "L");
	leg2->AddEntry(gFraction2, "Correct endpoints (ID_{first} = ID_{last})", "L");
	leg2->AddEntry(gFraction3, "Close endpoints (#pm 0.5 mm, #pm 0.5#circ)", "L");
	leg2->Draw();

	c2->Update();

	c3->cd();
	TGraph *gEfficiency = new TGraph(nEnergies, arrayEfficiencyEnergy, arrayEfficiencyFinal);
	TGraph *gEfficiency2 = new TGraph(nEnergies, arrayEfficiencyEnergy, arrayEfficiencyFinal2);
	gEfficiency->SetTitle("Efficiency of tracking algorithm at n_{p} = 500;Energy [MeV];Tracks reconstructed / n_{p}");
	gEfficiency->SetMinimum(0.7);
	gEfficiency->SetMaximum(1);
	gEfficiency->GetXaxis()->SetTitleFont(22);
	gEfficiency->GetYaxis()->SetTitleFont(22);
	gEfficiency->GetXaxis()->SetTitleOffset(1.2);
	gEfficiency->GetYaxis()->SetTitleOffset(1.2);
	gEfficiency->GetXaxis()->SetLabelFont(22);
	gEfficiency->GetYaxis()->SetLabelFont(22);
	gEfficiency->SetLineColor(kGreen+3);
	gEfficiency->SetLineWidth(3);
	gEfficiency2->SetLineColor(kGreen-3);
	gEfficiency2->SetLineWidth(3);
	gEfficiency->Draw("AL");
	gEfficiency2->Draw("L");

	c1->SaveAs("OutputFiles/figures/finalPlotsForArticle/estimated_ranges_all_energies.eps");
	c1->SaveAs("OutputFiles/figures/finalPlotsForArticle/estimated_ranges_all_energies.root");
	c2->SaveAs("OutputFiles/figures/finalPlotsForArticle/fraction_of_correct_tracks.eps");
	c2->SaveAs("OutputFiles/figures/finalPlotsForArticle/fraction_of_correct_tracks.root");
}
