#ifndef Track_conversion_h
#define Track_conversion_h

#include <TClonesArray.h>

#include "Classes/Cluster/Cluster.h"

Float_t getEnergyFromWEPL(Float_t wepl);
Float_t getEnergyFromTL(Float_t tl);

Double_t getQ(Double_t a1, Double_t a2,Double_t a3, Double_t a4);
Double_t getP(Double_t a1, Double_t a2, Double_t a3);

Float_t getTLFromEnergy(Float_t energy, Double_t a1, Double_t a2, Double_t a3, Double_t a4);
Float_t getEnergyFromTL(Float_t range, Double_t a1, Double_t a2, Double_t a3, Double_t a4);

Float_t getTLFromEnergy(Float_t energy);
Float_t getTLFromWEPL(Float_t wepl);
Float_t getWEPLFromEnergy(Float_t energy);
Float_t getWEPLFactorFromEnergy ( Float_t energy);
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
Float_t getEnergyFromUnit(Float_t unit);
Float_t getUnitStragglingFromEnergy(Float_t energy, Float_t sigma_energy);

#endif
