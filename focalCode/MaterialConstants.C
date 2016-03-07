#include <TObject.h>
#include "MaterialConstants.h"
#include "Constants.h"
#include <vector>

using namespace std;

void MaterialConstants() {
	p_water = 1.725517;
	alpha_water = 0.027681;
	alpha_prime_water = 0.0087;
	pinv_water = 1/p_water;
	
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
