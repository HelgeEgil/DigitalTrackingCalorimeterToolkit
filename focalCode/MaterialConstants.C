#include <TObject.h>
#include "MaterialConstants.h"
#include "Constants.h"
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

void MaterialConstants() {
	firstUpperLayerZ = 0.3;
	firstLowerLayerZ = 0.6;

	p_water_high = 1.725517; 		// above 150 MeV
	alpha_water_high = 0.027681;	// above 150 Mev
	
	p_water_low = 1.77251; 			// below 150 MeV
	alpha_water_low = 0.02195; 		// below 150 MeV

	p_water = p_water_high;
	alpha_water = alpha_water_high;
	pinv_water = 1/p_water_high;
	alpha_prime_water = 0.0087;
	
	p_aluminum = 1.677;
	alpha_aluminum = 0.01738;
	alpha_prime_aluminum = 0.017102;
	
	p_tungsten = 1.6413;
	alpha_tungsten = 0.00495;
	alpha_prime_tungsten = 0.1086784;

	if (kMaterial == kTungsten) {
		memcpy(kPLFocal, kPLFocal_W, sizeof(kPLFocal_W));
		memcpy(kWEPLRatio, kWEPLRatio_W, sizeof(kWEPLRatio_W));
		memcpy(kStraggling, kStraggling_W, sizeof(kStraggling_W));
//		nLayers = 24;
		nLayers = 41;
		p = p_tungsten;
		alpha = alpha_tungsten; // mm
		alpha_prime = alpha_prime_tungsten; // MeV^2 / mm
	}

	else if (kMaterial == kAluminum) {
		memcpy(kPLFocal, kPLFocal_Al, sizeof(kPLFocal_Al));
		memcpy(kWEPLRatio, kWEPLRatio_Al, sizeof(kWEPLRatio_Al));
		memcpy(kStraggling, kStraggling_Al, sizeof(kStraggling_Al));
		nLayers = 41;

		p = p_aluminum;
		alpha = alpha_aluminum;
		alpha_prime = alpha_prime_aluminum;
	}

	else if (kMaterial == kPMMA) {
		memcpy(kPLFocal, kPLFocal_PMMA, sizeof(kPLFocal_PMMA));
		memcpy(kWEPLRatio, kWEPLRatio_PMMA, sizeof(kWEPLRatio_PMMA));
		memcpy(kStraggling, kStraggling_PMMA, sizeof(kStraggling_PMMA));
		nLayers = 65;

		alpha_prime = alpha_prime_water; // MUST VERIFY
		p = p_water; // MUST VERIFY
		alpha = alpha_water; // MUST VERIFY
	}

	pinv = 1./p;
	firstHalfLayer = 1.65;
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
	
	if (nScintillators == 0) {
		energyLoss = 0;
	}

	else if (nScintillators == 1) {
		// only 4x4 scintillator
		energyLoss = 439.79 * pow(energy, -0.87);
	}
	else if (nScintillators == 2) {
		// 4x4 plus one of the other .. They are equal in thickness and material
		energyLoss = 530.5 * pow(energy, -0.822);
	}
	else if (nScintillators == 3) {
		// 4x4 plus both 1x4 and 4x1 scintillators
		energyLoss = 650.42 * pow(energy, -0.804);
	}

	else {
		cout << "ASSERTION ERROR ON NUMBER OF SCINTILLATORS!\n";
		energyLoss = 0;
	}

	return energyLoss;
}

Float_t getEnergyLossErrorFromScintillators(Int_t nScintillators) {
	Float_t energyLossError = 0;

	if (nScintillators == 0) {
		energyLossError = 0;
	}

	else if (nScintillators == 1) {
		energyLossError = 0.31;
	}

	else if (nScintillators == 2) {
		energyLossError = 0.37;
	}

	else if (nScintillators == 3) {
		energyLossError = 0.44;
	}

	return energyLossError;
}

Float_t getEnergyLossFromAluminumAbsorber(Float_t energy) {
	Float_t energyLoss = 92.135 * pow(energy, -0.8);
	return energyLoss;
}

Float_t getEnergyLossErrorFromAluminumAbsorber() {
	return 0.16;
}
