#include "Tools.h"
#include <iostream>
#include <cmath>
#include <TGraph.h>
#include "Constants.h"
#include "MaterialConstants.h"
#include <TClonesArray.h>

using namespace std;

Float_t getEnergyFromTL(Float_t tl) {
	return pow(tl / alpha, pinv);
}

Float_t getTLFromWEPL(Float_t wepl) {
	return alpha * pow(wepl / alpha_water, p / p_water);
}

Float_t getTLFromEnergy(Float_t energy) {
	return alpha * pow(energy, p);
}

Float_t getEnergyFromWEPL(Float_t wepl) {
	return pow(wepl / alpha_water, 1/(p_water));
}

Float_t getWEPLFactorFromEnergy(Float_t energy) {
	return (alpha_water / alpha) * pow(energy, p_water - p);
}

Float_t getWEPLFromTL(Float_t tl) {
	return alpha_water * pow(tl / alpha, p_water / p);
}

Float_t getWEPLFromEnergy(Float_t energy) {
	return alpha_water * pow(energy, p_water);
}

Float_t getRangeStragglingFromRange(Float_t range) {
	Float_t sigma = sqrt(alpha_prime * (pow(p, 2) * pow(alpha, 2*pinv) / (3-2*pinv) * pow(range, 3-2*pinv))) ;
	return sigma;
}

Float_t getRangeStragglingFromEnergy(Float_t energy) {
	Float_t range = getTLFromEnergy(energy);
	return getRangeStragglingFromRange(range);
}

Float_t getWEPLStragglingFromWEPL(Float_t wepl) {
	Float_t sigma = sqrt(alpha_prime_water * (pow(p_water, 2) * pow(alpha_water, 2*pinv_water) / (3-2*pinv_water) * pow(wepl, 3-2*pinv_water))) ;
	return sigma;
}

Float_t getEnergyStragglingFromEnergy(Float_t energy) {
	Float_t range = getTLFromEnergy(energy);
	return getEnergyStragglingFromRange(range);
}

Float_t getEnergyStragglingFromRange(Float_t range) {
	Float_t straggling = getRangeStragglingFromRange (range);

	Float_t upper = range + straggling / 2;
	Float_t lower = range - straggling / 2;

	Float_t upper_energy = getEnergyFromTL(upper);
	Float_t lower_energy = getEnergyFromTL(lower);

	Float_t energy_straggling = (upper_energy - lower_energy);

	return energy_straggling;
}
