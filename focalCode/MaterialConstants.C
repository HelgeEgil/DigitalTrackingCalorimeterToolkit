#include <TObject.h>
#include "MaterialConstants.h"
#include "Constants.h"
#include <vector>

using namespace std;

// Float_t kPLFocal[nPLEnergies];
// Float_t kWEPLRatio[nPLEnergies];
// Int_t nLayers;

void MaterialConstants() {

	if (kMaterial == kTungsten) {
		memcpy(kPLFocal, kPLFocal_W, sizeof(kPLFocal_W));
		memcpy(kWEPLRatio, kWEPLRatio_W, sizeof(kWEPLRatio_W));
		nLayers = 24;
	}

	else if (kMaterial == kAluminum) {
		memcpy(kPLFocal, kPLFocal_W, sizeof(kPLFocal_Al));
		memcpy(kWEPLRatio, kWEPLRatio_W, sizeof(kWEPLRatio_Al));
		nLayers = 41;
	}

	else if (kMaterial == kPMMA) {
		memcpy(kPLFocal, kPLFocal_W, sizeof(kPLFocal_PMMA));
		memcpy(kWEPLRatio, kWEPLRatio_W, sizeof(kWEPLRatio_PMMA));
		nLayers = 65;
	}
}