#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <TObject.h>
#include <TSpline.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "GlobalConstants/RangeAndEnergyCalculations.h"

using namespace std;

Float_t getTLFromEnergy(Float_t energy, TSpline3 *spline) {
   return spline->Eval(energy);
}

Float_t getEnergyFromTL(Float_t range, TSpline3 *invSpline) {
   return invSpline->Eval(range);
}

Float_t getEnergyAtTL(Float_t E0, Float_t depth, TSpline3 *spline, TSpline3 *invSpline) {
   Float_t range = spline->Eval(E0);
   Float_t energy = invSpline->Eval(range - depth);
   return energy;
}

/*
Float_t getTLFromEnergy(Float_t energy, Double_t a1, Double_t b1, Double_t b2, Double_t g1, Double_t g2) {
   Double_t sum1 = b1 * (1 - exp( -g1 * energy ));
   Double_t sum2 = b2 * (1 - exp( -g2 * energy ));
   Float_t tl = a1 * energy * (1 + sum1 + sum2);

   return tl;  
}

Float_t getEnergyFromTL(Float_t range, Double_t c1, Double_t c2, Double_t c3, 
                        Double_t c4, Double_t c5, Double_t l1, Double_t l2, 
                        Double_t l3, Double_t l4, Double_t l5) {
   Double_t sum1 = c1 * exp(-l1 * range);
   Double_t sum2 = c2 * exp(-l2 * range);
   Double_t sum3 = c3 * exp(-l3 * range);
   Double_t sum4 = c4 * exp(-l4 * range);
   Double_t sum5 = c5 * exp(-l5 * range);
   Float_t energy = range * (sum1 + sum2 + sum3 + sum4 + sum5);

   return energy;
}

Float_t getEnergyAtTL(Float_t E0, Float_t depth, Double_t c1, Double_t c2, 
                      Double_t c3, Double_t c4, Double_t c5, Double_t l1, 
                      Double_t l2, Double_t l3, Double_t l4, Double_t l5, 
                      Double_t a1, Double_t b1, Double_t b2, Double_t g1, Double_t g2) {
   Double_t range = getTLFromEnergy(E0, a1, b1, b2, g1, g2);
   cout << "RANGE = " << range << endl;
   if (range < depth) return 0;

   Double_t sum1 = c1 * exp(-l1 * (range - depth));
   Double_t sum2 = c2 * exp(-l2 * (range - depth));
   Double_t sum3 = c3 * exp(-l3 * (range - depth));
   Double_t sum4 = c4 * exp(-l4 * (range - depth));
   Double_t sum5 = c5 * exp(-l5 * (range - depth));

   Float_t energy = (range - depth) * (sum1 + sum2 + sum3 + sum4 + sum5);

   return energy;
}
*/

/*
//////////////////////////
// Conversion to energy //
//////////////////////////

Float_t getEnergyFromTL(Float_t range) {
   return getEnergyFromTL(range, c1_material, c2_material, c3_material, c4_material, c5_material,
                                 l1_material, l2_material, l3_material, l4_material, l5_material);
}

Float_t getEnergyFromWEPL(Float_t wepl) {
   return getEnergyFromTL(wepl, c1_water, c2_water, c3_water, c4_water, c5_water,
                                 l1_water, l2_water, l3_water, l4_water, l5_water);
}

Float_t getEnergyAtTL(Float_t E0, Float_t depth) {
   return getEnergyAtTL(E0, depth, c1_material, c2_material, c3_material, c4_material, c5_material,
                                 l1_material, l2_material, l3_material, l4_material, l5_material,
                                 a1_material, b1_material, b2_material, g1_material, g2_material);
}

Float_t getEnergyAtTLFromPureAluminium(Float_t E0, Float_t depth) {
   return getEnergyAtTL(E0, depth, c1_pure_aluminium, c2_pure_aluminium, c3_pure_aluminium, c4_pure_aluminium, c5_pure_aluminium,
                                 l1_pure_aluminium, l2_pure_aluminium, l3_pure_aluminium, l4_pure_aluminium, l5_pure_aluminium,
                                 a1_pure_aluminium, b1_pure_aluminium, b2_pure_aluminium, g1_pure_aluminium, g2_pure_aluminium);
}
*/

Float_t  getEnergyFromTL(Float_t range) {
   return getEnergyFromTL(range, splineMaterialInv);
}

Float_t getEnergyFromWEPL(Float_t wepl) {
   return getEnergyFromTL(wepl, splineWaterInv);
}

Float_t getEnergyAtTL(Float_t E0, Float_t depth) {
   return getEnergyAtTL(E0, depth, splineMaterial, splineMaterialInv);
}

Float_t getEnergyAtTLFromPureAluminium(Float_t E0, Float_t depth) {
   return getEnergyAtTL(E0, depth, splinePureAl, splinePureAlInv);
}

//////////////////////////
// Conversion to range  //
//////////////////////////
/*
Float_t getTLFromEnergy(Float_t energy) {
   return getTLFromEnergy(energy, a1_material, b1_material, b2_material, g1_material, g2_material);
}

Float_t getWEPLFromEnergy(Float_t energy) {
   return getTLFromEnergy(energy, a1_water, b1_water, b2_water, g1_water, g2_water);
}
*/

Float_t getTLFromEnergy(Float_t energy) {
   return getTLFromEnergy(energy, splineMaterial);
}

Float_t getWEPLFromEnergy(Float_t energy) {
   return getTLFromEnergy(energy, splineWater);
}

Float_t getTLFromWEPL(Float_t wepl) {
   Float_t energy = getEnergyFromWEPL(wepl);
   Float_t tl = getTLFromEnergy(energy);

   return tl;
}

//////////////////////////
// Conversion to WEPL   //
//////////////////////////

Float_t getWEPLFactorFromEnergy(Float_t energy) {
   Float_t range = getTLFromEnergy(energy);
   Float_t wepl = getWEPLFromEnergy(energy);

   return wepl / range;
}

Float_t getWEPLFromTL(Float_t tl) {
   Double_t energy = getEnergyFromTL(tl);
   Double_t wepl = getWEPLFromEnergy(energy);

   return wepl;
}


/////////////////////////////////////
// calculation of range straggling //
/////////////////////////////////////

Float_t getTLStragglingFromTL(Float_t tl, Float_t sigma_energy) {
   Float_t energy = getEnergyFromTL(tl);
   
   Float_t sigma_a = alpha_prime * (pow(p, 2) * pow(alpha, 2/p) / (3-2/p) * pow(tl, 3-2/p));
   Float_t sigma_b = pow(sigma_energy * alpha * p, 2) * pow(energy, 2*p-2);
   
   Float_t sigma = sqrt(sigma_a + sigma_b);
   
   return sigma;
}

Float_t getTLStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
   Float_t tl = getTLFromEnergy(energy);
   return getTLStragglingFromTL(tl, sigma_energy);
}

/////////////////////////////////////
// Conversion to WEPL straggling   //
/////////////////////////////////////

Float_t getWEPLStragglingFromWEPL(Float_t wepl, Float_t sigma_energy) {
   Float_t energy = getEnergyFromWEPL(wepl);
   
   Float_t tlStraggling = getTLStragglingFromEnergy(energy, sigma_energy);
   Float_t tl = getTLFromEnergy(energy);

   Float_t upperTL = tl + tlStraggling / 2;
   Float_t lowerTL = tl - tlStraggling / 2;

   Float_t weplStraggling = getWEPLFromTL(upperTL) - getWEPLFromTL(lowerTL);

   return weplStraggling;
}

Float_t getWEPLStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
   Float_t wepl = getWEPLFromEnergy(energy);
   return getWEPLStragglingFromWEPL(wepl, sigma_energy);
}  

///////////////////////////////////////
// Conversion to energy straggling   //
///////////////////////////////////////

Float_t getEnergyStragglingFromTL(Float_t tl, Float_t sigma_energy) {
   Float_t straggling = getTLStragglingFromTL (tl, sigma_energy);

   Float_t upper = tl + straggling / 2;
   Float_t lower = tl - straggling / 2;

   Float_t upper_energy = getEnergyFromTL(upper);
   Float_t lower_energy = getEnergyFromTL(lower);

   Float_t energy_straggling = (upper_energy - lower_energy);

   return energy_straggling;
}

Float_t getEnergyStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
   Float_t tl = getTLFromEnergy(energy);
   return getEnergyStragglingFromTL(tl, sigma_energy);
}

Float_t getEnergyStragglingFromTLStraggling(Float_t tl, Float_t TLStraggling) {
   Float_t upper = tl + TLStraggling / 2;
   Float_t lower = tl - TLStraggling / 2;

   Float_t upper_energy = getEnergyFromTL(upper);
   Float_t lower_energy = getEnergyFromTL(lower);

   Float_t energy_straggling = (upper_energy - lower_energy);

   return energy_straggling;
}

Float_t getEnergyStragglingFromWEPLStraggling(Float_t wepl, Float_t WEPLStraggling) {
   Float_t upper = wepl + WEPLStraggling / 2;
   Float_t lower = wepl - WEPLStraggling / 2;

   Float_t upper_energy = getEnergyFromWEPL(upper);
   Float_t lower_energy = getEnergyFromWEPL(lower);

   Float_t energy_straggling = (upper_energy - lower_energy);

   return energy_straggling;
}

/////////////////////////////////////////////
// Conversion to/from defined kOutputUnit  //
/////////////////////////////////////////////

Float_t getUnitFromEnergy(Float_t energy) {
   Float_t res = 0;
   if       (kOutputUnit == kEnergy)    res = energy;
   else if  (kOutputUnit == kPhysical)  res = getTLFromEnergy(energy);
   else if  (kOutputUnit == kWEPL)      res = getWEPLFromEnergy(energy);

   return res;
}

Float_t getUnitFromWEPL(Float_t wepl) {
   Float_t res = 0;
   if       (kOutputUnit == kEnergy)   res = getEnergyFromWEPL(wepl);
   else if  (kOutputUnit == kPhysical) res = getTLFromWEPL(wepl);
   else if  (kOutputUnit == kWEPL)     res = wepl;
   
   return res;
}

Float_t getEnergyFromUnit(Float_t unit) {
   Float_t res = 0;

   if    (kOutputUnit == kEnergy)   res = unit;
   else if  (kOutputUnit == kPhysical) res = getEnergyFromTL(unit);
   else if (kOutputUnit == kWEPL)      res = getEnergyFromWEPL(unit);
   
   return res;
}

Float_t getUnitStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
   Float_t res = 0;
   if      (kOutputUnit == kEnergy)    res = getEnergyStragglingFromEnergy(energy, sigma_energy);
   else if (kOutputUnit == kWEPL)      res = getWEPLStragglingFromEnergy(energy, sigma_energy);
   else if (kOutputUnit == kPhysical)  res = getTLStragglingFromEnergy(energy, sigma_energy);

   return res;
}

Float_t getEnergyLossFromScintillators(Float_t energy, Int_t nScintillators) {
   Float_t energyLoss = 0;
   
   if      (nScintillators == 0) energyLoss = 0;
   else if (nScintillators == 1) energyLoss = 214.88444 * pow(energy, -0.7264956);
   else if (nScintillators == 2) energyLoss = 341.79255 * pow(energy, -0.7357198);
   else if (nScintillators == 3) energyLoss = 482.63531 * pow(energy, -0.7453624);
   
   else {
      cout << "ASSERTION ERROR ON NUMBER OF SCINTILLATORS!\n";
      energyLoss = 0;
   }
   
   return energyLoss;
}

Float_t getEnergyLossErrorFromScintillators(Int_t nScintillators) {
   Float_t energyLossError = 0;

   if       (nScintillators == 0) energyLossError = 0;
   else if (nScintillators == 1) energyLossError = 0.30;
   else if (nScintillators == 2) energyLossError = 0.38;
   else if (nScintillators == 3) energyLossError = 0.44;
   
   return energyLossError;
}

Float_t getEnergyLossFromAluminumAbsorber(Float_t energy) {
   // We assume that the energy loss function is constant in the high energy range
   // (somewhat true)
   // and calculate the energy loss of kAbsorbatorThickness mm Aluminum by scaling the
   // measured MC energy loss in a 1.5 mm Al absorber

   return  62.00441 * pow(energy, -0.7158403) * (kAbsorbatorThickness/1.5);
}

Float_t getEnergyLossErrorFromAluminumAbsorber() {
   return 0.135;
}
