#include <iostream>
#include <cmath>

#include <TGraph.h>
#include <TClonesArray.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;

// Full revision of range-energy and depth-dose relationship as of 2016-06
// Using Ulmer Rad Phys and Chem 76 (2007) 1089-1107
// In short:
// R = a1 E0 * (1 + sum_(k=1)^2 (bk - bk * exp(-gk * E0)))
// E = R * sum_(i=1)^5 (ck * exp(-lambdak * R))
//		 ( and  E(z) = (R-z) * sum_(i=1)^5 (ck * exp(-lambdak * (R-z))) )
// -dE/dz = E(z) / (R-z) - sum_(k=1)^5 (lambdak * ck * (R-z) * exp(-lambdak * (R-z)))
// The parameters a1, bk(1-2), gk(1-2), ck(1-5), lambdak(1-5) are found through range-energy fitting on different materials in GATE simulations, and defined below
// All calculations are performed in Classes/Track/conversionFunctions.C


Float_t getTLFromEnergy(Float_t energy, Double_t a1, Double_t b1, Double_t b2, Double_t g1, Double_t g2) {
	Double_t sum1 = b1 * (1 - exp( -g1 * energy ));
	Double_t sum2 = b2 * (1 - exp( -g2 * energy ));
	Float_t tl = a1 * energy * (1 + sum1 + sum2);

	return tl;	
}

Float_t getEnergyFromTL(Float_t range, Double_t c1, Double_t c2, Double_t c3, Double_t c4, Double_t c5,
													Double_t l1, Double_t l2, Double_t l3, Double_t l4, Double_t l5) {
	Double_t sum1 = c1 * exp(-l1 * range);
	Double_t sum2 = c2 * exp(-l2 * range);
	Double_t sum3 = c3 * exp(-l3 * range);
	Double_t sum4 = c4 * exp(-l4 * range);
	Double_t sum5 = c5 * exp(-l5 * range);
	Float_t energy = range * (sum1 + sum2 + sum3 + sum4 + sum5);

	return energy;
}

Float_t getEnergyAtTL(Float_t E0, Float_t depth, Double_t c1, Double_t c2, Double_t c3, Double_t c4, Double_t c5,
													Double_t l1, Double_t l2, Double_t l3, Double_t l4, Double_t l5) {

	Double_t range = getTLFromEnergy(E0);

	if (range < depth) return 0;

	Double_t sum1 = c1 * exp(-l1 * (range - depth));
	Double_t sum2 = c2 * exp(-l2 * (range - depth));
	Double_t sum3 = c3 * exp(-l3 * (range - depth));
	Double_t sum4 = c4 * exp(-l4 * (range - depth));
	Double_t sum5 = c5 * exp(-l5 * (range - depth));
	Float_t energy = (range - depth) * (sum1 + sum2 + sum3 + sum4 + sum5);

	return energy;
}

Float_t getTLFromEnergy(Float_t energy) {
	return getTLFromEnergy(energy, a1_material, b1_material, b2_material, g1_material, g2_material);
}

Float_t getEnergyFromTL(Float_t range) {
	return getEnergyFromTL(range, c1_material, c2_material, c3_material, c4_material, c5_material,
											l1_material, l2_material, l3_material, l4_material, l5_material);
}

Float_t getEnergyAtTL(Float_t E0, Float_t depth) {
	return getEnergyAtTL(E0, depth, c1_material, c2_material, c3_material, c4_material, c5_material,
											l1_material, l2_material, l3_material, l4_material, l5_material);
}

Float_t getWEPLFromEnergy(Float_t energy) {
	return getTLFromEnergy(energy, a1_water, b1_water, b2_water, g1_water, g2_water);
}

Float_t getEnergyFromWEPL(Float_t wepl) {
	return getEnergyFromTL(wepl, c1_water, c2_water, c3_water, c4_water, c5_water,
											l1_water, l2_water, l3_water, l4_water, l5_water);
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


Float_t getEnergyStragglingFromTL(Float_t tl, Float_t sigma_energy) {
	Float_t straggling = getTLStragglingFromTL (tl, sigma_energy);

	Float_t upper = tl + straggling / 2;
	Float_t lower = tl - straggling / 2;

	Float_t upper_energy = getEnergyFromTL(upper);
	Float_t lower_energy = getEnergyFromTL(lower);

	Float_t energy_straggling = (upper_energy - lower_energy);

	return energy_straggling;
}

Float_t getEnergyStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
	Float_t tl = getTLFromEnergy(energy);
	return getEnergyStragglingFromTL(tl, sigma_energy);
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
