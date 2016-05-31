#ifndef MaterialConstants_h
#define MaterialConstants_h

#include <TObject.h>

Int_t nLayers;

// range = alpha * energy ^ p
// range straggling is proportionate to alpha_prime
Double_t p, alpha, alpha_prime;
Double_t p_water, alpha_water, alpha_prime_water;
Double_t p_aluminum, alpha_aluminum, alpha_prime_aluminum;
Double_t p_tungsten, alpha_tungsten, alpha_prime_tungsten;

// range = a1 * E^3 + a2 * E^2 + a3 * e + a4: More accurate
Double_t a1_tungsten, a2_tungsten, a3_tungsten, a4_tungsten, a5_tungsten;
Double_t a1_material, a2_material, a3_material, a4_material, a5_material;
Double_t a1_water, a2_water, a3_water, a4_water, a5_water;

Double_t firstUpperLayerZ, firstLowerLayerZ;

Float_t getEnergyLossFromScintillators(Float_t energy, Int_t nScintillators);
Float_t getEnergyLossErrorFromScintillators(Int_t nScintillators);
Float_t getEnergyLossFromAluminumAbsorber(Float_t energy);
Float_t getEnergyLossErrorFromAluminumAbsorber();
Double_t getLayerPositionmm(Int_t i);
Float_t getSigmaEnergy(Int_t energy);

#endif
