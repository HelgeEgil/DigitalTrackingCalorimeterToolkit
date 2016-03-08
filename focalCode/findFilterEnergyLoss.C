#define findFilterEnergyLoss_cxx
#include "findFilterEnergyLoss.h"
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
#include <TPad.h>
#include <TList.h>
#include <TSpectrum.h>
#include <TPolyMarker.h>
#include <TMath.h>

using namespace std;

void findFilterEnergyLoss::BinLogY(TH2 *h) {
   TAxis *axis = h->GetYaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;

   Axis_t *new_bins = new Axis_t[bins+1];

   for (int i=0; i <= bins; i++) {
      new_bins[i] = TMath::Power(10, from + i * width);
   }

   axis->Set(bins, new_bins);
   delete new_bins;
}


void findFilterEnergyLoss::Loop()
{
	Int_t energy = run_energy;
	
	Double_t alpha = 0.0014467;
	Double_t alpha_prime = 0.203815;
	Double_t p = 1.707283;
	
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	TH1F *hEnergyLoss = new TH1F("hEnergyLoss", Form("Energy loss for a %d MeV proton beam on a 1.5 mm aluminum foil", energy), 5000, 0, 250);
	TH1F *hEnergy= new TH1F("hEnergy", Form("Final energy for a %d MeV proton beam on a 1.5 mm aluminum foil", energy), 5000, 0, 250);

	Int_t lastID = -1;
	Float_t sum_edep = 0;

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
	
		if (ientry < 0) {
			cout << "Aborting run at jentry = " << jentry << endl;
			break;
		}

		nb = fChain->GetEntry(jentry);   nbytes += nb;

		if (abs(posX) < 1 && abs(posY) < 1) {
			continue;
		}
		
		sum_edep += edep;
		if (eventID != lastID) {
			hEnergyLoss->Fill(sum_edep);
			hEnergy->Fill(energy - sum_edep);
			sum_edep = 0;
		}
		lastID = eventID;
	}
	
	Int_t res;
	
	TSpectrum *s = new TSpectrum(6);
	TSpectrum *sLoss = new TSpectrum(6);
	
	res = s->Search(hEnergy);
	TList *functions = hEnergy->GetListOfFunctions();
	TPolyMarker *pm = (TPolyMarker*) functions->FindObject("TPolyMarker");
	
	res = sLoss->Search(hEnergyLoss);
	TList *functionsLoss = hEnergyLoss->GetListOfFunctions();
	TPolyMarker *pmLoss = (TPolyMarker*) functionsLoss->FindObject("TPolyMarker");
	
	Int_t n = pm->GetN();
	Int_t nLoss = pmLoss->GetN();
	
	Double_t *x = pm->GetX();
	Double_t *xLoss = pmLoss->GetX();
	
	for (Int_t i=0; i<n; i++) {
		cout << "Final energy peak " << i + 1 << " at " << x[i] << endl;
	}
	
	for (Int_t i=0; i<nLoss; i++) {
		cout << "Energy loss peak " << i + 1 << " at " << xLoss[i] << endl;
	}
	
	TF1 *fEnergy = new TF1("fEnergy", "gaus(0)", energy-50, energy+15);
	TF1 *fEnergyLoss = new TF1("fEnergyLoss", "gaus(0)", 0, 15);
	
	fEnergy->SetParameters(100, energy-5, 0.3);
	fEnergyLoss->SetParameters(100, 5, 0.3);
	
	c1->cd();
	hEnergy->SetXTitle("Final energy");
	hEnergy->SetYTitle("Number of particles");
	hEnergy->Draw();
	
	c2->cd();
	hEnergyLoss->SetXTitle("Total energy loss");
	hEnergyLoss->SetYTitle("Number of particles");
	hEnergyLoss->Draw();
	
	hEnergy->Fit("fEnergy", "B, Q", "", 0, energy+50);
	hEnergyLoss->Fit("fEnergyLoss", "B, Q", "", 0, 15);
	
	cout << "For a " << energy << " MeV beam, the final energy is " << fEnergy->GetParameter(1) 
		 << " MeV, with a total energy loss of " << fEnergyLoss->GetParameter(1) << " +- " << fEnergy->GetParameter(2) << " MeV.\n";
	
}
