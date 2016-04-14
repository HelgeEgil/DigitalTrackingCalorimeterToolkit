#ifndef MaterialConstants_h
#define MaterialConstants_h

#include <TObject.h>

const Int_t nPLEnergies = 96;

Float_t kPLFocal[nPLEnergies];
Float_t kWEPLRatio[nPLEnergies];
Float_t kStraggling[nPLEnergies];
Int_t nLayers;

enum eWaterSetting {kLow, kHigh};

Double_t p, pinv, alpha, alpha_prime;
Double_t p_water, alpha_water, alpha_prime_water, pinv_water;
Double_t p_water_low, alpha_water_low;
Double_t p_water_high, alpha_water_high;
Double_t p_aluminum, alpha_aluminum, alpha_prime_aluminum;
Double_t p_tungsten, alpha_tungsten, alpha_prime_tungsten;
Double_t firstHalfLayer;

Double_t firstUpperLayerZ;
Double_t firstLowerLayerZ;

Float_t getEnergyLossFromScintillators(Float_t energy, Int_t nScintillators);
Float_t getEnergyLossErrorFromScintillators(Int_t nScintillators);
Float_t getEnergyLossFromAluminumAbsorber(Float_t energy);
Float_t getEnergyLossErrorFromAluminumAbsorber();
Double_t getLayerPositionmm(Int_t i);
Float_t getSigmaEnergy(Int_t energy);
void setWaterPValues(Int_t setting);

#endif
