#include <TObject.h>
#include "MaterialConstants.h"
#include "Constants.h"
#include <vector>

using namespace std;

void MaterialConstants() {

	if (kMaterial == kTungsten) {
		memcpy(kPLFocal, kPLFocal_W, sizeof(kPLFocal_W));
		memcpy(kWEPLRatio, kWEPLRatio_W, sizeof(kWEPLRatio_W));
		memcpy(kStraggling, kStraggling_W, sizeof(kStraggling_W));
//		nLayers = 24;
		nLayers = 41;
	}

	else if (kMaterial == kAluminum) {
		memcpy(kPLFocal, kPLFocal_Al, sizeof(kPLFocal_Al));
		memcpy(kWEPLRatio, kWEPLRatio_Al, sizeof(kWEPLRatio_Al));
		memcpy(kStraggling, kStraggling_Al, sizeof(kStraggling_Al));
		nLayers = 41;
	}

	else if (kMaterial == kPMMA) {
		memcpy(kPLFocal, kPLFocal_PMMA, sizeof(kPLFocal_PMMA));
		memcpy(kWEPLRatio, kWEPLRatio_PMMA, sizeof(kWEPLRatio_PMMA));
		memcpy(kStraggling, kStraggling_PMMA, sizeof(kStraggling_PMMA));
		nLayers = 65;
	}
}
