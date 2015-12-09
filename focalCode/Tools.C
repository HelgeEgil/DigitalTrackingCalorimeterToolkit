/*
 * Tools.C
 *
 *  Created on: Sep 2, 2015
 *      Author: local1
 */

#include "TObject.h"
#include "Cluster.h"
#include "Hit.h"
#include "Constants.h"
#include <vector>
#include <algorithm>

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

	Float_t R_mm = 0.0019 * pow(energy,1.8) * 10;
	Double_t fitval = scale / ( 0.0554 * pow((R_mm - depth), 0.444) );

	if (isnan(fitval)) fitval = 0;

	return fitval;
}
