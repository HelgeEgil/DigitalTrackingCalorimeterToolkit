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

Float_t getTLFromMaterial(Float_t tl_mat, Int_t material) {
	Float_t tl = 0;
	if (material == kTungsten) {
		tl = alpha * pow(tl_mat / alpha_tungsten, p / p_tungsten);
	}
	else if (material == kAluminum) {
		tl = alpha * pow(tl_mat / alpha_aluminum, p / p_aluminum);
	}
	else {
		cout << "MATERIAL" << material << " NOT FOUND! Reverting to input TL\n";
		tl = tl_mat;
	}

	return tl;
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

Float_t getTLStragglingFromTL(Float_t tl, Float_t sigma_energy) {
	Float_t energy = getEnergyFromTL(tl);
	
	Float_t sigma_a = alpha_prime * (pow(p, 2) * pow(alpha, 2*pinv) / (3-2*pinv) * pow(tl, 3-2*pinv));
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
	
	Float_t sigma_a = alpha_prime_water * (pow(p_water, 2) * pow(alpha_water, 2*pinv_water) / (3-2*pinv_water) * pow(wepl, 3-2*pinv_water)) ;
	Float_t sigma_b = pow(sigma_energy * alpha_water * p_water, 2) * pow(energy, 2*p_water-2);

	Float_t sigma = sqrt(sigma_a + sigma_b);

	return sigma;
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