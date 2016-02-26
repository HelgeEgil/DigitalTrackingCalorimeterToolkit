#ifndef Track_conversion_h
#define Track_conversion_h
#include "TClonesArray.h"
#include "Cluster.h"

Float_t getEnergyFromWEPL(Float_t wepl);
Float_t getEnergyFromTL(Float_t tl);

Float_t getRangeFromEnergy(Float_t energy);
Float_t getTLFromEnergy(Float_t energy);
Float_t getTLFromWEPL(Float_t wepl);
Float_t getWEPLFromEnergy(Float_t energy);
Float_t getWEPLFactorFromEnergy ( Float_t energy );
Float_t getWEPLFromTL(Float_t tl);

Float_t getRangeStragglingFromEnergy(Float_t energy);
Float_t getRangeStragglingFromRange ( Float_t range );
Float_t getWEPLStragglingFromWEPL(Float_t wepl);
Float_t getWEPLStragglingFromEnergy(Float_t energy);
Float_t getEnergyStragglingFromRange(Float_t range);
Float_t getEnergyStragglingFromEnergy(Float_t energy);


#endif