#ifndef Track_conversion_h
#define Track_conversion_h
#include "TClonesArray.h"
#include "Cluster.h"

Float_t getEnergyFromWEPL(Float_t wepl);
Float_t getEnergyFromTL(Float_t tl);

Float_t getRangeFromEnergy(Float_t energy);
Float_t getTLFromEnergy(Float_t energy);
Float_t getTLFromWEPL(Float_t wepl);
Float_t getTLFromMaterial(Float_t tl_mat, Int_t material);
Float_t getWEPLFromEnergy(Float_t energy);
Float_t getWEPLFactorFromEnergy ( Float_t energy );
Float_t getWEPLFromTL(Float_t tl);

Float_t getTLStragglingFromEnergy(Float_t energy, Float_t sigma_energy);
Float_t getEnergyStragglingFromTLStraggling(Float_t range, Float_t rangeStraggling);
Float_t getEnergyStragglingFromWEPLStraggling(Float_t range, Float_t rangeStraggling);
Float_t getTLStragglingFromTL ( Float_t range, Float_t sigma_energy);
Float_t getWEPLStragglingFromWEPL(Float_t wepl, Float_t sigma_energy);
Float_t getWEPLStragglingFromEnergy(Float_t energy, Float_t sigma_energy);
Float_t getEnergyStragglingFromTL(Float_t range, Float_t sigma_energy);
Float_t getEnergyStragglingFromEnergy(Float_t energy, Float_t sigma_energy);

Float_t getUnitFromEnergy(Float_t energy);
Float_t getUnitStragglingFromEnergy(Float_t energy, Float_t sigma_energy);

#endif
