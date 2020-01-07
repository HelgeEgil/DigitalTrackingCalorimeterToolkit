#ifndef RangeAndEnergyCalculations_h
#define RangeAndEnergyCalculations_h

// Actual calculation functions
Float_t getTLFromEnergy(Float_t energy, Double_t a1, Double_t b1, Double_t b2, Double_t g1, Double_t g2);

Float_t getEnergyFromTL(Float_t range, Double_t c1, Double_t c2, Double_t c3, Double_t c4, Double_t c5,
                        Double_t l1, Double_t l2, Double_t l3, Double_t l4, Double_t l5);

Float_t getEnergyAtTL(Float_t E0, Double_t depth, Double_t c1, Double_t c2, Double_t c3, Double_t c4, Double_t c5,
                        Double_t l1, Double_t l2, Double_t l3, Double_t l4, Double_t l5,
                        Double_t a1, Double_t b1, Double_t b2, Double_t g1, Double_t g2);

// Conversion to energy
Float_t getEnergyFromTL(Float_t tl);
Float_t getEnergyFromWEPL(Float_t wepl);
Float_t getEnergyAtTL(Float_t energy, Float_t depth);
Float_t getEnergyAtWEPL(Float_t energy, Float_t depth);
Float_t getEnergyAtTLFromPureAluminum(Float_t energy, Float_t depth);

Float_t getWETFromLayer(Float_t layer);
Float_t getWEPLFromLayer(Float_t layer);
Float_t getWETFromDegrader(Float_t degrader);

// Conversion to range
Float_t getTLFromEnergy(Float_t energy);
Float_t getTLFromWELP(Float_t wepl);

// Conversion to Water Equivalent Path Length
Float_t getWEPLFromEnergy(Float_t energy);
Float_t getWEPLFromTL(Float_t tl);
Float_t getWEPLFactorFromEnergy ( Float_t energy);
Float_t getWEPLFactorFromWEPL(Float_t wepl);

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
Float_t getUnitFromWEPL(Float_t wepl);
Float_t getUnitFromTL(Float_t tl);
Float_t getEnergyFromUnit(Float_t unit);
Float_t getEnergyAtUnit(Float_t unit, Float_t degraderthickness);
Float_t getUnitStragglingFromEnergy(Float_t energy, Float_t sigma_energy);

// Calculation of energy loss
Float_t  getEnergyLossFromScintillators(Float_t energy, Int_t nScintillators);
Float_t  getEnergyLossErrorFromScintillators(Int_t nScintillators);
Float_t  getEnergyLossFromAluminumAbsorber(Float_t energy);
Float_t  getEnergyLossFromTracker(Float_t energy);
Float_t  getEnergyLossErrorFromAluminumAbsorber();

Float_t  getEnergyLossAtWEPL(Float_t wepl, Float_t initialEnegy);
Float_t  getEnergyLossAtTL(Float_t tl, Float_t initialEnergy);

Float_t getEnergyFromDegraderThickness(Double_t degraderThickness);

// Fallback Bragg-Kleeman functions
Float_t  getBKTL(Float_t energy, Float_t aa, Float_t pp);
Float_t  getBKEnergy(Float_t tl, Float_t aa, Float_t pp);
Float_t  getBKEnergyHigh(Float_t tl);
Float_t  getBKTLHigh(Float_t energy);
Float_t  getBKEnergyLow(Float_t tl);
Float_t  getBKTLLow(Float_t energy);

#endif
