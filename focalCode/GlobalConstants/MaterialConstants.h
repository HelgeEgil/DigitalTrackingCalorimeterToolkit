#ifndef MaterialConstants_h
#define MaterialConstants_h

#include <TObject.h>

class Hit;
class Cluster;

Int_t nLayers;

// range = alpha * energy ^ p
// range straggling is proportionate to alpha_prime
Double_t p, alpha, alpha_prime;
Double_t p_water, alpha_water, alpha_prime_water;
Double_t p_aluminum, alpha_aluminum, alpha_prime_aluminum;
Double_t p_tungsten, alpha_tungsten, alpha_prime_tungsten;

Double_t a1_water, b1_water, b2_water, g1_water, g2_water;
Double_t c1_water, c2_water, c3_water, c4_water, c5_water;
Double_t l1_water, l2_water, l3_water, l4_water, l5_water;
Double_t a1_material, b1_material, b2_material, g1_material, g2_material;
Double_t c1_material, c2_material, c3_material, c4_material, c5_material;
Double_t l1_material, l2_material, l3_material, l4_material, l5_material;
Double_t a1_tungsten, b1_tungsten, b2_tungsten, g1_tungsten, g2_tungsten;
Double_t c1_tungsten, c2_tungsten, c3_tungsten, c4_tungsten, c5_tungsten;
Double_t l1_tungsten, l2_tungsten, l3_tungsten, l4_tungsten, l5_tungsten;
Double_t a1_aluminium, b1_aluminium, b2_aluminium, g1_aluminium, g2_aluminium;
Double_t c1_aluminium, c2_aluminium, c3_aluminium, c4_aluminium, c5_aluminium;
Double_t l1_aluminium, l2_aluminium, l3_aluminium, l4_aluminium, l5_aluminium;
Double_t a1_pmma, b1_pmma, b2_pmma, g1_pmma, g2_pmma;
Double_t c1_pmma, c2_pmma, c3_pmma, c4_pmma, c5_pmma;
Double_t l1_pmma, l2_pmma, l3_pmma, l4_pmma, l5_pmma;

Double_t X0, X0_tungsten, X0_aluminum, X0_pmma, X0_firstlayer;
Double_t proton_mass;

Float_t mcs_radius_per_layer[100]; // change if max(nLayers) > 100

Double_t firstUpperLayerZ, firstLowerLayerZ;

void     loadRangeValuesForTungsten();
void     loadRangeValuesForAluminium();
Float_t  getEnergyLossFromScintillators(Float_t energy, Int_t nScintillators);
Float_t  getEnergyLossErrorFromScintillators(Int_t nScintillators);
Float_t  getEnergyLossFromAluminumAbsorber(Float_t energy);
Float_t  getEnergyLossErrorFromAluminumAbsorber();
Double_t getLayerPositionmm(Int_t i);
Float_t  getSigmaEnergy(Int_t energy);
Bool_t   isChipLowResistivity(Int_t chipIdx);
Float_t  getChipCalibrationFactor(Int_t chipIdx);
Bool_t   isBadData(Cluster * estimatedPosition);

#endif
