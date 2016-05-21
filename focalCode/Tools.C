/*
 * Tools.C
 *
 *  Created on: Sep 2, 2015
 *      Author: local1
 */

#include "Cluster.h"
#include "Hit.h"
#include "Constants.h"
#include "Track_conversion.h"
#include "MaterialConstants.h"
#include <vector>
#include <algorithm>
#include <TObject.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>


using namespace std;

Bool_t isItemInVector(Int_t i, vector<Int_t> *v) {
	return find(v->begin(), v->end(), i) != v->end();
}

Bool_t existsEnergyFile(Int_t energy) {
	Bool_t res = false;
	for (Int_t i=0; i<nEnergies; i++) {
		if (energies[i] == energy) res = true;
	}
	return res;
}

Float_t diffmmXY(Cluster *p1, Cluster *p2) {
	if (!p1 || !p2)
		return -1;

	Double_t diffx = p2->getXmm() - p1->getXmm();
	Double_t diffy = p2->getYmm() - p1->getYmm();

	return sqrt(pow(diffx,2) + pow(diffy,2));
}

Float_t diffmmXY(Hit *h1, Hit *h2) {
	if (!h1 || !h2)
		return -1;

	Double_t diffx = h2->getXmm() - h1->getXmm();
	Double_t diffy = h2->getYmm() - h1->getYmm();

	return sqrt(pow(diffx,2) + pow(diffy,2));
}

Float_t diffmmXYZ(Cluster *p1, Cluster *p2) {
	return sqrt(pow(p2->getXmm() - p1->getXmm(), 2) + 
				pow(p2->getYmm() - p1->getYmm(), 2) +
				pow(p2->getLayermm() - p1->getLayermm(), 2));
}

Double_t fitfunc_DBP(Double_t *v, Double_t *par) {
	// Based on Bortfeld and Schlegel 1996

	Float_t depth = v[0];
	Float_t energy = par[0];
	Float_t scale = par[1];

	// alpha, p and pinv are values from MaterialConstants.C
	
	Double_t fitval = 0;
	
	if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		Float_t range = getWEPLFromEnergy(energy);
		fitval = scale / ( p_water * pow(alpha_water, pinv_water) * pow((range - depth), 1-pinv_water) );
	}
	else if (kOutputUnit == kPhysical) {
		Float_t range = getTLFromEnergy(energy);
		fitval = scale / ( p * pow(alpha, pinv) * pow((range - depth), 1-pinv) );
	}
	
	if (isnan(fitval)) fitval = 0;

	return fitval;
}

Double_t double_landau(Double_t *v, Double_t *par) {
	Double_t x = v[0];
	Double_t sigma1 = par[1];
	Double_t sigma2 = par[2];
	Double_t mpv1 = par[0];
	Double_t mpv2 = mpv1 + dz;	
	if (kOutputUnit == kWEPL || kOutputUnit == kEnergy)
		mpv2 = mpv1 + getWEPLFactorFromEnergy(run_energy) * dz;
	
	return TMath::Landau(x, mpv1, sigma1, 0) + TMath::Landau(x, mpv2, sigma2, 0);
}

char * getMaterialChar() {
	char *sMaterial;

	if (kMaterial == kTungsten) sMaterial = (char*) "tungsten";
	if (kMaterial == kAluminum) sMaterial = (char*) "aluminum";
	if (kMaterial == kPMMA) sMaterial = (char*) "PMMA";

	return sMaterial;
}

char * getDataTypeChar(Int_t dataType) {
	char * sDataType = (char*) "";
	if (dataType == kMC) sDataType = (char*) "MC";
	else if (dataType == kData) sDataType = (char*) "Exp. data";

	return sDataType;
}

Int_t getFWxMInRange(TH1F* h, Float_t first, Float_t last, Int_t div) {
	// find the full with at maximum / div (2 = FWHM, 3 = FWTM, etc.)
	
	Int_t firstRangeBin = h->GetBin(first);
	Int_t lastRangeBin = h->GetBin(last);
	
//	cout << "firstRangeBin = " << firstRangeBin << ", lastRangeBin = " << lastRangeBin << endl;
	
	Float_t maxBin = 0;
	Int_t maxIdx = 0;
	for (Int_t i=firstRangeBin; i<=lastRangeBin; i++) {
// 		cout << h->GetBinContent(i) << endl;
		if (h->GetBinContent(i) > maxBin ) {
			maxBin = h->GetBinContent(i);
// 			cout << "new max = " << maxBin << endl;
			maxIdx = i;
		}
	}
	
// 	cout << "max value is " << maxBin << " at " << maxIdx << endl;
	
	Float_t firstCenter = 0, lastCenter = 0;
	for (Int_t i=firstRangeBin; i<=lastRangeBin; i++) {
		if (!firstCenter) {
			if (h->GetBinContent(i) > maxBin/div) {
				firstCenter = h->GetBinCenter(i);
// 				cout << "firstCenter at " << firstCenter << endl;
			}
		}
		
		else {
			if (!lastCenter) {
				if (h->GetBinContent(i) < maxBin/div) {
					lastCenter = h->GetBinCenter(i-1);
// 					cout << "lastCenter at " << firstCenter << endl;
				}
			}
		}
	}
	
	if (!firstCenter || !lastCenter) return 0;	
	Float_t width = lastCenter - firstCenter;

	return width * 2 * sqrt( 2 * log(div));
}

Int_t getMinimumTrackLength(Float_t energy) {
	Int_t minTL = 0;

	if (energy < 150) minTL = 5;
	else if (energy < 170) minTL = 10;
	else if (energy < 190) minTL = 15;
	else if (energy < 200) minTL = 20;
	else if (energy < 230) minTL = 23;

	return minTL;
}

Float_t quadratureAdd(Float_t a, Float_t b) {
	return sqrt(pow(a,2) + pow(b, 2));
}

void convertXYToWEPL(Float_t *x_energy, Float_t *y_energy, Int_t eventID) {
	Float_t WEPLFactor = getWEPLFactorFromEnergy(run_energy);
	Long64_t n=0;
	Long64_t j=0;

	for (Long64_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) {
		x_energy[i] *=  WEPLFactor;
	}
}

Float_t getEnergyFromXY(Float_t *x_energy, Float_t *y_energy, Int_t eventID) {
	Float_t WEPLFactor = getWEPLFactorFromEnergy(run_energy);
	Long64_t n=0;
	Long64_t j=0;

	for (Long64_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) {
		if (y_energy[i] == 0) {
			n = j;
			break;
		}
		j++;
	}

	return getEnergyFromWEPL(x_energy[eventID*sizeOfEventID + n-1]);
}

Double_t correctForEnergyParameterisation(Float_t energy) {
	// By using the range = alpha * pow ( energy, p ) parameterisation to find
	// the depth dose curve used in Bragg Peak fitting, a 0.2 % error is introduced
	// This may be fixed by calculating back the WEPL range using the r = a E^p method,
	// and then finding the corresponding energy using the cubic and more accurate method

	Double_t wepl = alpha_water * pow(energy, p_water);
	Double_t corrected_energy = getEnergyFromWEPL(wepl);

	cout << "Energy = " << energy << ", corrected energy = " << corrected_energy << endl;

	return corrected_energy;
}
