#ifndef Track_conversion_h
#define Track_conversion_h
#include "TClonesArray.h"
#include "Cluster.h"

Float_t getWEPLFactorFromEnergy ( Float_t energy );
Float_t getWEPLFromTL(Float_t tl);
Float_t getWEPLFromEnergy(Float_t energy);
Float_t getEnergyFromWEPL(Float_t wepl);
Float_t getEnergyFromTL(Float_t tl);
Float_t getEnergyStragglingFromEnergy(Float_t energy);
Float_t getRangeStraggling ( Float_t energy );
Float_t getRangeFromEnergy(Float_t energy);
Float_t getEnergyStragglingFromRangeStraggling(Float_t range, Float_t rangeStraggling);

#endif