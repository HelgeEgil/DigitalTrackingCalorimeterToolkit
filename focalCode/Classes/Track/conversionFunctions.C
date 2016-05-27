#include <iostream>
#include <cmath>

#include <TGraph.h>
#include <TClonesArray.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;

// New calculation as of 2016-05-16
// Using QSGP-BIC-EMY instead of EMSTANDARD_OPT3

Double_t getQ(Double_t a1, Double_t a2,Double_t a3, Double_t a4) {
		return (2*pow(a2,3) - 9*a1*a2*a3 + 27*pow(a1,2) * a4) / (27*pow(a1,3));
}

Double_t getP(Double_t a1, Double_t a2, Double_t a3) {
	return (3*a1 * a3 - pow(a2, 2)) / (3 * pow(a1, 2));
}

Double_t inverse_cubic(Double_t p, Double_t q) {
	return 2 * sqrt(-p/3) * cos(1/3. * acos((3*q)/(2*p) * sqrt(-3./p)) - 2*3.14159265/3.);
}

Float_t getTLFromEnergy(Float_t energy, Double_t a1, Double_t a2, Double_t a3, Double_t a4) {
	return a1 * pow(energy, 3) + a2 * pow(energy, 2) + a3 * energy + a4;
}

Float_t getEnergyFromTL(Float_t range, Double_t a1, Double_t a2, Double_t a3, Double_t a4) {
	Double_t p = getP(a1, a2, a3);
	Double_t q = getQ(a1, a2, a3, a4 - range);

	Double_t tk = inverse_cubic(p, q);
	Double_t energy = tk - a2 / (3 * a1);

	return energy;
}

Float_t getTLFromEnergy(Float_t energy) {
	return getTLFromEnergy(energy, a1_tungsten, a2_tungsten, a3_tungsten, a4_tungsten);
}

Float_t getEnergyFromTL(Float_t range) {
	return getEnergyFromTL(range, a1_tungsten, a2_tungsten, a3_tungsten, a4_tungsten);
}

Float_t getWEPLFromEnergy(Float_t energy) {
	return getTLFromEnergy(energy, a1_water, a2_water, a3_water, a4_water);
}

Float_t getEnergyFromWEPL(Float_t wepl) {
	return getEnergyFromTL(wepl, a1_water, a2_water, a3_water, a4_water);
}

Float_t getWEPLFactorFromEnergy(Float_t energy) {
	Float_t range = getTLFromEnergy(energy);
	Float_t wepl = getWEPLFromEnergy(energy);

	return wepl / range;
}

Float_t getWEPLFromTL(Float_t tl) {
	Double_t energy = getEnergyFromTL(tl);
	Double_t wepl = getWEPLFromEnergy(energy);

	return wepl;
}

Float_t getTLStragglingFromTL(Float_t tl, Float_t sigma_energy) {
	Float_t energy = getEnergyFromTL(tl);
	
	Float_t sigma_a = alpha_prime * (pow(p, 2) * pow(alpha, 2/p) / (3-2/p) * pow(tl, 3-2/p));
	Float_t sigma_b = pow(sigma_energy * alpha * p, 2) * pow(energy, 2*p-2);
	
	Float_t sigma = sqrt(sigma_a + sigma_b);
	
	return sigma;
}

Float_t getTLStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
	Float_t tl = getTLFromEnergy(energy);
	return getTLStragglingFromTL(tl, sigma_energy);
}

Float_t getWEPLStragglingFromWEPL(Float_t wepl, Float_t sigma_energy) {
	Float_t energy = getEnergyFromWEPL(wepl);
	
	Float_t tlStraggling = getTLStragglingFromEnergy(energy, sigma_energy);
	Float_t tl = getTLFromEnergy(energy);

	Float_t upperTL = tl + tlStraggling / 2;
	Float_t lowerTL = tl - tlStraggling / 2;

	Float_t weplStraggling = getWEPLFromTL(upperTL) - getWEPLFromTL(lowerTL);

	return weplStraggling;
}

Float_t getWEPLStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
	Float_t wepl = getWEPLFromEnergy(energy);
	return getWEPLStragglingFromWEPL(wepl, sigma_energy);
}	

Float_t getEnergyStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
	Float_t tl = getTLFromEnergy(energy);
	return getEnergyStragglingFromTL(tl, sigma_energy);
}

Float_t getEnergyStragglingFromTL(Float_t tl, Float_t sigma_energy) {
	Float_t straggling = getTLStragglingFromTL (tl, sigma_energy);

	Float_t upper = tl + straggling / 2;
	Float_t lower = tl - straggling / 2;

	Float_t upper_energy = getEnergyFromTL(upper);
	Float_t lower_energy = getEnergyFromTL(lower);

	Float_t energy_straggling = (upper_energy - lower_energy);

	return energy_straggling;
}

Float_t getEnergyStragglingFromTLStraggling(Float_t tl, Float_t TLStraggling) {
	Float_t upper = tl + TLStraggling / 2;
	Float_t lower = tl - TLStraggling / 2;

	Float_t upper_energy = getEnergyFromTL(upper);
	Float_t lower_energy = getEnergyFromTL(lower);

	Float_t energy_straggling = (upper_energy - lower_energy);

	return energy_straggling;
}

Float_t getEnergyStragglingFromWEPLStraggling(Float_t wepl, Float_t WEPLStraggling) {
	Float_t upper = wepl + WEPLStraggling / 2;
	Float_t lower = wepl - WEPLStraggling / 2;

	Float_t upper_energy = getEnergyFromWEPL(upper);
	Float_t lower_energy = getEnergyFromWEPL(lower);

	Float_t energy_straggling = (upper_energy - lower_energy);

	return energy_straggling;
}

Float_t getUnitFromEnergy(Float_t energy) {
	Float_t res = 0;
	if		  (kOutputUnit == kEnergy) 	res = energy;
	else if (kOutputUnit == kPhysical) 	res = getTLFromEnergy(energy);
	else if (kOutputUnit == kWEPL) 		res = getWEPLFromEnergy(energy);

	return res;
}

Float_t getEnergyFromUnit(Float_t unit) {
	Float_t res = 0;

	if		(kOutputUnit == kEnergy)	res = unit;
	else if	(kOutputUnit == kPhysical)	res = getEnergyFromTL(unit);
	else if (kOutputUnit == kWEPL)		res = getEnergyFromWEPL(unit);
	
	return res;
}

Float_t getUnitStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
	Float_t res = 0;
	if		  (kOutputUnit == kEnergy)		res = getEnergyStragglingFromEnergy(energy, sigma_energy);
	else if (kOutputUnit == kWEPL) 		res = getWEPLStragglingFromEnergy(energy, sigma_energy);
	else if (kOutputUnit == kPhysical) 	res = getTLStragglingFromEnergy(energy, sigma_energy);

	return res;
}
