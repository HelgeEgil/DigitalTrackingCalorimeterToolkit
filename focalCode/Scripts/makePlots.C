#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
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
   TCanvas *c1	= new TCanvas("c1", "c1", 2000, 1500);
	TCanvas *c2 = new TCanvas("c2", "Correct Tracks fraction", 2000, 1500);

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
	

   gStyle->SetOptStat(0);

	ifstream in;
	in.open("OutputFiles/result_makebraggpeakfit.csv");

	cout << "Opened file.\n";

	Float_t nomrange_, estrange_, sigmaRange_, lastRange_;
	Int_t energy_;

	Int_t nlines = 0;
	TNtuple *ntuple = new TNtuple("ntuple", "data from file", "energy_:nomrange_:estrange_:sigmaRange:lastRange_");
	
	Int_t MC2Data = 88;

	while (1) {
		in >> energy_ >> nomrange_ >> estrange_ >> sigmaRange_ >> lastRange_ ;

		if (!in.good()) {
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
	Float_t factor, np, correctLast, correctWhole, lastIsFirst;
	Float_t arrayFractionX[200] = {0};
	Float_t arrayFractionY[200] = {0};
	Int_t nlines2 = 0;
	while (1) {
		in2 >> factor >> np >> correctLast >> correctWhole >> lastIsFirst;

		if (!in2.good()) break;
		
		arrayFractionX[nlines2] = np;
		arrayFractionY[nlines2] = lastIsFirst * 100;
		cout << "line " << nlines2 << ", np = " << np << ", lastIsFirst = " << lastIsFirst << endl;

		nlines2++;
	}
	cout << "Found " << nlines2 << " lines in lastLayerCorrect.\n";
	
	in2.close();

	c1->cd();
	
	TGraphErrors *hMC = new TGraphErrors(MC2Data, arrayE, arrayMC, arrayEE, arrayEMC);
	TGraphErrors *hData = new TGraphErrors(nlines-MC2Data, arrayE2, arrayData, arrayEE2, arrayEData);
	TGraphErrors *pstar = new TGraphErrors(MC2Data, arrayE, arrayPSTAR, arrayEE, arrayEPstar);
	
	hMC->SetTitleFont(22);
	hMC->GetXaxis()->SetTitleFont(22);
	hMC->GetYaxis()->SetTitleFont(22);
	hMC->GetXaxis()->SetTitleOffset(1.2);
	hMC->GetYaxis()->SetTitleOffset(1.2);

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

	TLegend *leg = new TLegend(0.15, 0.68, 0.44, 0.85);
	leg->SetTextSize(0.035);
	leg->AddEntry(pstar, "PSTAR range", "L");
	leg->AddEntry(hMC, "Monte Carlo", "PE");
	leg->AddEntry(hData, "Experimental data", "PE");
	leg->Draw();

	c1->Update();

	c2->cd();
	TGraph *gFraction = new TGraph(nlines2, arrayFractionX, arrayFractionY);
	gFraction->GetXaxis()->SetRangeUser(10, 6000);
	gFraction->GetXaxis()->SetTitleOffset(1.2);
	gFraction->GetYaxis()->SetTitleOffset(1.2);
	gFraction->SetTitle("Fraction of correctly reconstructed tracks;Number of protons in frame;Fraction of tracks where first and last track ID is equal [%]");
	gFraction->SetLineColor(kBlue);
	gFraction->SetLineWidth(3);
	gFraction->Draw("AL");

	gPad->SetLogx();
	gFraction->GetXaxis()->SetNoExponent();

	TLine *line80 = new TLine(0, 80, 234, 80);
	line80->SetLineColor(kRed);
	line80->SetLineWidth(2);
	line80->SetLineStyle(7);
	line80->Draw("same");

	TLine *vLine80 = new TLine(234, 80, 234, 0);
	vLine80->SetLineColor(kRed);
	vLine80->SetLineWidth(2);
	vLine80->SetLineStyle(7);
	vLine80->Draw("same");
	
	TLine *line1000 = new TLine(0, 37.95, 1000, 37.95);
	line1000->SetLineColor(kRed);
	line1000->SetLineWidth(2);
	line1000->SetLineStyle(7);
	line1000->Draw("same");

	TLine *vLine1000 = new TLine(1000, 37.95, 1000, 0);
	vLine1000->SetLineColor(kRed);
	vLine1000->SetLineWidth(2);
	vLine1000->SetLineStyle(7);
	vLine1000->Draw("same");

	c2->Update();

	c1->SaveAs("OutputFiles/estimated_ranges_all_energies.pdf");
	c2->SaveAs("OutputFiles/fraction_of_correct_tracks.pdf");
}
