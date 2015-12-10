#include "Tools.h"
#include <iostream>
#include <cmath>
#include <TGraph.h>
#include "Constants.h"
#include "MaterialConstants.h"
#include <TClonesArray.h>

using namespace std;

Float_t getEnergyFromTL(Float_t tl) {
	// TODO: optimize via stackoverflow.com/questions/11396860
	// (and create bigger table for tungsten....)

	for (Int_t i=1; i<nPLEnergies; i++) {
		if (kPLFocal[i] < tl) { continue; }
		else {
			Float_t ratio = (tl - kPLFocal[i-1]) / (kPLFocal[i] - kPLFocal[i-1]);
			return kPLEnergies[i-1] + ratio * (kPLEnergies[i] - kPLEnergies[i-1]);
		}
	}
}

Float_t getEnergyFromWEPL(Float_t wepl) {

	for (Int_t i=1; i<nPLEnergies; i++) {
		if (kPLFocal[i] * kWEPLRatio[i] < wepl) { continue; }
		else {
			Float_t ratio = (wepl - kPLFocal[i-1] * kWEPLRatio[i-1]) / (kPLFocal[i]* kWEPLRatio[i] - kPLFocal[i-1]* kWEPLRatio[i-1]);
			return kPLEnergies[i-1] + ratio * (kPLEnergies[i] - kPLEnergies[i-1]);
		}
	}
}

Float_t getWEPLFactorFromEnergy(Float_t energy) {
	for (Int_t i=1; i<nPLEnergies; i++) {
		if (kPLEnergies[i] < energy) { continue; }
		else {
			Float_t ratio = (energy - kPLEnergies[i-1]) / (kPLEnergies[i] - kPLEnergies[i-1]);
			return kWEPLRatio[i-1] + ratio * (kWEPLRatio[i] - kWEPLRatio[i-1]);
		}
	}
}

Float_t getWEPLFromTL(Float_t tl) {

	for (Int_t i=1; i<nPLEnergies; i++) {
		if (kPLFocal[i] < tl) { continue; }
		else {
			Float_t ratio = (tl - kPLFocal[i-1]) / (kPLFocal[i] - kPLFocal[i-1]);
			return tl * (kWEPLRatio[i-1] + ratio * (kWEPLRatio[i] - kWEPLRatio[i-1]));
		}
	}
}

Float_t getWEPLFromEnergy(Float_t energy) {
	Float_t tl = getRangeFromEnergy(energy);
	return getWEPLFromTL(tl);
}

Float_t getRangeStraggling(Float_t energy) {
	for (Int_t i=1; i<nPLEnergies; i++) {
		if (kPLEnergies[i] < energy) { continue; }
		else {
			Float_t ratio = (energy - kPLEnergies[i-1]) / (kPLEnergies[i] - kPLEnergies[i-1]);
			return kStraggling[i-1] + ratio * (kStraggling[i] - kStraggling[i-1]);
		}
	}
}

Float_t getRangeFromEnergy(Float_t energy) {
	for (Int_t i=1; i<nPLEnergies; i++) {
		if (kPLEnergies[i] < energy) { continue; }
		else {
			Float_t ratio = (energy - kPLEnergies[i-1]) / (kPLEnergies[i] - kPLEnergies[i-1]);
			return kPLFocal[i-1] + ratio * (kPLFocal[i] - kPLFocal[i-1]);
		}
	}
}

Float_t getEnergyStragglingFromRangeStraggling(Float_t range, Float_t rangeStraggling) {
	Float_t upperEnergy = getEnergyFromTL(range + rangeStraggling / 2);
	Float_t lowerEnergy = getEnergyFromTL(range - rangeStraggling / 2);
	return upperEnergy - lowerEnergy;
}

Float_t getEnergyStragglingFromEnergy(Float_t energy) {
	Float_t range = getRangeFromEnergy(energy);
	if (!range) return 0;
	
	Float_t rangeStraggling = getRangeStraggling(energy);
	if (!rangeStraggling) return 0;
	
	Float_t energyStraggling = getEnergyStragglingFromRangeStraggling(range, rangeStraggling);
	
	return energyStraggling;
}