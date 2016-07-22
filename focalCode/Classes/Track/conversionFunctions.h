#ifndef Track_conversion_h
#define Track_conversion_h

#include <TClonesArray.h>

#include "Classes/Cluster/Cluster.h"

// Actual calculation functions
Float_t getTLFromEnergy(Float_t energy, Double_t a1, Double_t a2, Double_t a3, Double_t a4);
Float_t getEnergyFromTL(Float_t range, Double_t a1, Double_t a2, Double_t a3, Double_t a4);
Float_t getEnergyAtTL(Float_t E0, Float_t depth, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

// Conversion to energy
Float_t getEnergyFromTL(Float_t tl);
Float_t getEnergyFromWEPL(Float_t wepl);
Float_t getEnergyAtTL(Float_t energy, Float_t depth);

// Conversion to range
Float_t getTLFromEnergy(Float_t energy);

// Conversion to Water Equivalent Path Length
Float_t getWEPLFromEnergy(Float_t energy);
Float_t getWEPLFromTL(Float_t tl);
Float_t getWEPLFactorFromEnergy ( Float_t energy);

// Calculation of Range straggling
Float_t getTLStragglingFromTL( Float_t range, Float_t sigma_energy);
Float_t getTLStragglingFromEnergy(Float_t energy, Float_t sigma_energy);

// Conversion to WEPL straggling
Float_t getWEPLStragglingFromWEPL(Float_t wepl, Float_t sigma_energy);
Float_t getWEPLStragglingFromEnergy(Float_t energy, Float_t sigma_energy);

// Conversion to energy straggling
Float_t getEnergyStragglingFromTL(Float_t range, Float_t sigma_energy);
Float_t getEnergyStragglingFromEnergy(Float_t energy, Float_t sigma_energy);
Float_t getEnergyStragglingFromTLStraggling(Float_t range, Float_t rangeStraggling);
Float_t getEnergyStragglingFromWEPLStraggling(Float_t range, Float_t rangeStraggling);

// Conversion to/from defined kOutputUnit
Float_t getUnitFromEnergy(Float_t energy);
Float_t getEnergyFromUnit(Float_t unit);
Float_t getUnitStragglingFromEnergy(Float_t energy, Float_t sigma_energy);

#endif
