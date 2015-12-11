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
#include <vector>
#include <algorithm>
#include <TObject.h>
#include <TAxis.h>
#include <TGraph.h>

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

Float_t diffmm(Cluster *p1, Cluster *p2) {
	if (!p1 || !p2)
		return -1;

	Double_t diffx = p2->getXmm() - p1->getXmm();
	Double_t diffy = p2->getYmm() - p1->getYmm();

	return sqrt(pow(diffx,2) + pow(diffy,2));
}

Float_t diffmm(Hit *h1, Hit *h2) {
	if (!h1 || !h2)
		return -1;

	Double_t diffx = h2->getXmm() - h1->getXmm();
	Double_t diffy = h2->getYmm() - h1->getYmm();

	return sqrt(pow(diffx,2) + pow(diffy,2));
}

Double_t fitfunc_DBP(Double_t *v, Double_t *par) {
	// Based on Bortfeld and Schlegel 1996

	Float_t depth = v[0];
	Float_t energy = par[0];
	Float_t scale = par[1];

// 	Float_t alpha = 0.0019;
	Float_t alpha = 0.033; // from getPValues()
	Float_t p = 1.6891; // from getPValues();
	Float_t pinv = 1/p;
	Float_t alpha_p = p * pow(alpha, pinv);
	
// 	Float_t R_mm = 0.0019 * pow(energy,1.8) * 10;
	Float_t R_mm = getWEPLFromEnergy(energy);
	
	Double_t fitval = scale / ( alpha_p * pow((R_mm - depth), 1-pinv) );

	if (isnan(fitval)) fitval = 0;

	return fitval;
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

Int_t getMinimumTrackLength(Int_t energy) {
	Int_t minTL = 0;

	if (energy < 150) minTL = 5;
	else if (energy < 170) minTL = 10;
	else if (energy < 190) minTL = 15;
	else if (energy < 200) minTL = 20;
	else if (energy < 230) minTL = 23;

	return minTL;
}
