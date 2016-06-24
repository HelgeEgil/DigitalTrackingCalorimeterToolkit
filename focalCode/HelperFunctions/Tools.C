#include <vector>
#include <algorithm>

#include <TObject.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>

#include "../Classes/Cluster/Cluster.h"
#include "../Classes/Hit/Hit.h"
#include "../GlobalConstants/Constants.h"
#include "../GlobalConstants/MaterialConstants.h"
#include "../Classes/Track/conversionFunctions.h"


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

Float_t diffXY(Cluster *p, Hit *h) {
	if (!p || !h) return 1e5;

	Float_t diffx = p->getX() - h->getX();
	Float_t diffy = p->getY() - h->getY();

	return sqrt(pow(diffx, 2) + pow(diffy, 2));
}

Float_t diffXY(Cluster *p1, Cluster *p2) {
	if (!p1 || !p2)
		return -1;

	Double_t diffx = p2->getX() - p1->getX();
	Double_t diffy = p2->getY() - p1->getY();

	return sqrt(pow(diffx,2) + pow(diffy,2));
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

	Float_t depth = v[0];
	Float_t energy = par[0];
	Float_t scale = par[1];

	Double_t fitval = 0;

	Float_t range = getWEPLFromEnergy(energy); // cubic calculation, more accurate
	fitval = scale / ( p_water * pow(alpha_water, 1/p_water) * pow((range - depth), 1-1/p_water) ); // power calculation, less accurate but with the depth dose function

	/*
	 * I tried to use an improved depth-dose curve in order to improve accuracy
	 * Alas it didn't work. Maybe you could give it a go?
	 * Start with: W. Ulmer, Rad. Phys. and Chem. 76 (2007)
	 * and implement the depth-dose curve (28) with maybe a nice Bragg Peak
	 * I had some bumps... Back to Bragg-Kleeman.
	 *
	Double_t rz = range - depth;
	if (rz<0 || depth < 0) return 0;

	Double_t sum_1a = rz * c1_water * exp(-l1_water * rz);
	Double_t sum_1b = rz * c2_water * exp(-l2_water * rz);
	Double_t sum_1c = rz * c3_water * exp(-l3_water * rz);
	Double_t sum_1d = rz * c4_water * exp(-l4_water * rz);
	Double_t sum_1e = rz * c5_water * exp(-l5_water * rz);
	
	Double_t Ez = sum_1a + sum_1b + sum_1c + sum_1d + sum_1e;
	if (Ez < 0) return 0;

	Double_t sum_2a = l1_water * sum_1a;
	Double_t sum_2b = l2_water * sum_1b;
	Double_t sum_2c = l3_water * sum_1c;
	Double_t sum_2d = l4_water * sum_1d;
	Double_t sum_2e = l5_water * sum_1e;

	Double_t Ek = sum_2a + sum_2b + sum_2c + sum_2d + sum_2e;

	Double_t dEdz = scale * (Ez / rz - Ek);
	if (isnan(dEdz)) dEdz = 0;

	cout << "dEdz = " << dEdz << ", fitval = " << fitval << endl;
	*/

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

	if (kMaterial == kTungsten) sMaterial = (char*) "Tungsten";
	if (kMaterial == kAluminum) sMaterial = (char*) "Aluminium";
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

	return corrected_energy;
}

Float_t calculateMCS(Float_t energy, Float_t depth, Bool_t kUseDifferentX0) {
	Float_t gamma = (energy + proton_mass) / proton_mass;
	Float_t beta = sqrt(1 - pow(gamma, -2));
	Float_t momentum = gamma * beta * proton_mass;

	Float_t thisX0 = X0;
	if (kUseDifferentX0) {
		thisX0 = X0_firstlayer;
	}

	Float_t mcs = 13.6 / (beta * momentum) * sqrt(depth / thisX0) * (1 + 0.038 * log(depth / thisX0));

	return mcs; // radians
}

Float_t findMCSPixelRadiusForLayer(Float_t layer, Float_t E0) {
	Float_t currentEnergy = getEnergyAtTL(E0, getLayerPositionmm(layer));
	currentEnergy -= getAverageEnergyLoss(E0);
	
	Bool_t kUseDifferentX0 = false;

	Float_t depth = dz;
	if (layer == 0) {
		// Expected depth of scintillator + preabsorber (and no tungsten absorber before sensor)
		depth = 16.6; 
		kUseDifferentX0 = true;
	}
	
	Float_t restRange = getTLFromEnergy(currentEnergy);
	if (restRange < depth) depth = restRange;

	Float_t mcs = calculateMCS(currentEnergy, depth, kUseDifferentX0);
	Float_t radius = depth * tan (mcs);

	return radius;
}

Float_t findMCSAtLayerRad(Int_t layer, Float_t E0) {
	Float_t energy = E0 - getAverageEnergyLoss(E0);
	Float_t restRange = getTLFromEnergy(energy);

	Float_t mcs = calculateMCS(energy, restRange, false);
	return mcs;
}

void fillMCSRadiusList(Float_t angleFactor) {
	Float_t pixelRadius;

	if (!run_energy) {
		cout << "\033[1mRUN_ENERGY NOT DEFINED! ASSUMING 180 MeV.\033[0m\n";
	}

	for (Int_t i=0; i<nLayers; i++) {
		pixelRadius = findMCSPixelRadiusForLayer(i, run_energy);
		if (!isnan(pixelRadius)) {
			mcs_radius_per_layer[i] = pixelRadius * angleFactor;
		}
	}
}

Float_t getSearchRadiusForLayer(Int_t layer) {
	return mcs_radius_per_layer[layer];
}

Float_t getMCSAngleForLayer(Int_t layer) {
	Float_t radius = getSearchRadiusForLayer(layer);
	Float_t depth = dz;
	if (layer == 0) {
		depth = 16.6;
	}
	
	Float_t angle = atan( radius / depth );
	return angle * 180 / 3.14159265;
}

Float_t getAverageEnergyLoss(Float_t energy) {
	Float_t dE_1 = getEnergyLossFromScintillators(energy, 1);
	Float_t dE_2 = getEnergyLossFromScintillators(energy, 2);
	Float_t dE_3 = getEnergyLossFromScintillators(energy, 3);
	Float_t dE_Al = getEnergyLossFromAluminumAbsorber(energy);

	Float_t probabilityOf1Scintillator = 0.5625;
	Float_t probabilityOf2Scintillators = 0.375;
	Float_t probabilityOf3Scintillators = 0.0675;

	Float_t avgEnergyLoss = dE_1 * probabilityOf1Scintillator +
									dE_2 * probabilityOf2Scintillators +
									dE_3 * probabilityOf3Scintillators +
									dE_Al;

	return avgEnergyLoss;
}

Hit * sumHits(Hit * a, Hit * b) {
	Int_t x, y, layer, eventID;
	Float_t eDep;

	if (!a) cout << "!a\n";
	if (!b) cout << "!b\n";

	eDep = a->getEdep() + b->getEdep(); // sum 
	x = ( (int) a->getX() + (int) b->getX() + 0.5 ) / 2; // average
	y = ( (int) a->getY() + (int) b->getY() + 0.5 ) / 2; // average
	layer = a->getLayer();
	eventID = a->getEventID();

	if (a->getEventID() != b->getEventID() || a->getLayer() != b->getLayer()) {
		cout << "\033[1mASSERTION ERROR IN sumHits: THE EVENT IDS OR LAYER# ARE DIFFERENT!!\033[0m\n";
	}

	Hit * newHit = new Hit(x, y, layer, eventID, eDep);

	return newHit;
}

Cluster * getTrackPropagationToLayer(Track * track, Int_t layer) {
	Int_t last = track->GetEntriesFast() - 1;
	Int_t diffLayer = layer - track->getLayer(last);

	Cluster p1(track->getX(last-1), track->getY(last-1));
	Cluster p2(track->getX(last), track->getY(last));

	Cluster slope(p2.getX() - p1.getX(), p2.getY() - p1.getY());

	Float_t x = p2.getX() + diffLayer * slope.getX();
	Float_t y = p2.getY() + diffLayer * slope.getY();

	return new Cluster(x, y, layer);
}


Bool_t isPointOutOfBounds(Cluster *point, Float_t padding) {
	Bool_t isOutside;
   Float_t x = point->getX();
   Float_t y = point->getY();

   if (!point)
   	isOutside = kTRUE;
   else {
		if (x < 0 - padding || x > 2*nx + padding ||
			 y < 0 - padding || y > 2*ny + padding)
			isOutside = kTRUE;
		else
			isOutside = kFALSE;
   }

   return isOutside;
}
