#include <TObject.h>
#include "MaterialConstants.h"
#include "Constants.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

void MaterialConstants() {

	// New calculation as of 2016-05-16
	// Using QSGP-BIC-EMY instead of EMSTANDARD_OPT3

	a1_water = -0.0000051372;
	a2_water =  0.0066110840;
	a3_water =  0.1958689129;
	a4_water = -3.7279603175;

	a1_tungsten = -0.0000004333;
	a2_tungsten =  0.0006465835;
	a3_tungsten =  0.0435509812;
	a4_tungsten = -1.0445000000;

	firstUpperLayerZ = 0.3;
	firstLowerLayerZ = 0.6;
	
	p_water_high = 1.7547;
	alpha_water_high = 0.02387;
	
	p_water_low = 1.77251; 			// below 150 MeV
	alpha_water_low = 0.02195; 		// below 150 MeV

	p_water = p_water_high;
	alpha_water = alpha_water_high;
	pinv_water = 1/p_water_high;
	alpha_prime_water = 0.0087;
	
	p_aluminum = 1.677;
	alpha_aluminum = 0.01738;
	alpha_prime_aluminum = 0.017102;
	
	p_tungsten = 1.6677;
	alpha_tungsten = 0.004461;
	
	alpha_prime_tungsten = 0.1086784;

	if (kMaterial == kTungsten) {
//		nLayers = 24;
		nLayers = 41;
		p = p_tungsten;
		alpha = alpha_tungsten; // mm
		alpha_prime = alpha_prime_tungsten; // MeV^2 / mm
	}

	else if (kMaterial == kAluminum) {
		nLayers = 41;

		p = p_aluminum;
		alpha = alpha_aluminum;
		alpha_prime = alpha_prime_aluminum;
	}

	else if (kMaterial == kPMMA) {
		nLayers = 65;

		alpha_prime = alpha_prime_water; // MUST VERIFY
		p = p_water; // MUST VERIFY
		alpha = alpha_water; // MUST VERIFY
	}

	pinv = 1./p;
}

void setWaterPValues(Int_t setting) {
	if (setting == kLow) {
		p_water = p_water_low;
		alpha_water = alpha_water_low;
	}
	else if (setting == kHigh) {
		p_water = p_water_high;
		alpha_water = alpha_water_high;
	}
	
	pinv_water = 1/p_water;
}

Float_t getEnergyLossFromScintillators(Float_t energy, Int_t nScintillators) {
	Float_t energyLoss = 0;

	// New calculation as of 2016-05-16
	// Using QSGP-BIC-EMY instead of EMSTANDARD_OPT3
	
 	if 	  (nScintillators == 0) energyLoss = 0;
 	else if (nScintillators == 1) energyLoss = 214.88444 * pow(energy, -0.7264956);
 	else if (nScintillators == 2) energyLoss = 341.79255 * pow(energy, -0.7357198);
 	else if (nScintillators == 3) energyLoss = 482.63531 * pow(energy, -0.7453624);
	
	else {
		cout << "ASSERTION ERROR ON NUMBER OF SCINTILLATORS!\n";
		energyLoss = 0;
	}
	
	return energyLoss;
}

Float_t getEnergyLossErrorFromScintillators(Int_t nScintillators) {
	Float_t energyLossError = 0;
	
	// New calculation as of 2016-05-16
	// Using QSGP-BIC-EMY instead of EMSTANDARD_OPT3

	if 		(nScintillators == 0) energyLossError = 0;
	else if (nScintillators == 1) energyLossError = 0.30;
	else if (nScintillators == 2) energyLossError = 0.38;
	else if (nScintillators == 3) energyLossError = 0.44;
	
	return energyLossError;
}

Float_t getEnergyLossFromAluminumAbsorber(Float_t energy) {
	// New calculation as of 2016-05-16
	// Using QSGP-BIC-EMY instead of EMSTANDARD_OPT3

	return 62.00441 * pow(energy, -0.7158403);
}

Float_t getEnergyLossErrorFromAluminumAbsorber() {
	return 0.135;
}

Double_t getLayerPositionmm(Int_t i) {
	Double_t z = 0;

	if (i>0) {
		z = ( firstUpperLayerZ + firstLowerLayerZ ) / 2 + dz * i;
	}

	return z;
}

Float_t getSigmaEnergy(Int_t energy) { 
	Float_t sigma_energy = 0;
	
		if 	  (energy == 188) sigma_energy = 0.3; // approx
	else if (energy == 180) sigma_energy = 2.3;
	else if (energy == 170) sigma_energy = 4.5; // was 4.5
	else if (energy == 160) sigma_energy = 8.3;
	else if (energy == 150) sigma_energy = 20;
	
	return sigma_energy;
}

void setCalibratedTungstenPValues() {
// 	alpha_tungsten /= 1.028;
	alpha_tungsten /= 1.037;
	alpha = alpha_tungsten; // mm
}
