#ifndef rangeandenergycalculations_cxx
#define rangeandenergycalculations_cxx

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

Float_t getWETFromLayer(Float_t layer) {
   return splineWET->Eval(layer);
}

Float_t getWEPLFromLayer(Float_t layer) {
   return splineWEPL->Eval(layer);
}

Float_t getWETFromDegrader(Float_t degrader) {
   return splineWETFromDegrader->Eval(degrader);
}

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

Float_t  getEnergyFromTL(Float_t range) {
   if (range == 0) return 0;

   else if (range < splineMaterialInv->GetXmin()) {
      return getBKEnergyLow(range);
   }
   
   else if (range > splineMaterialInv->GetXmax()) {
      return getBKEnergyHigh(range);
   }

   return getEnergyFromTL(range, splineMaterialInv);
}

Float_t getEnergyFromWEPL(Float_t wepl) {
   return getEnergyFromTL(wepl, splineWaterInv);
}

Float_t getEnergyAtTL(Float_t E0, Float_t depth) {
   return getEnergyAtTL(E0, depth, splineMaterial, splineMaterialInv);
}

Float_t getEnergyAtWEPL(Float_t E0, Float_t depth) {
   return getEnergyAtTL(E0, depth, splineWater, splineWaterInv);
}

Float_t getEnergyAtTLFromPureAluminum(Float_t E0, Float_t depth) {
   return getEnergyAtTL(E0, depth, splinePureAl, splinePureAlInv);
}

//////////////////////////
// Conversion to range  //
//////////////////////////

Float_t getTLFromEnergy(Float_t energy) {
   if (energy == 0) return 0;

   else if (energy < splineMaterial->GetXmin()) {
      return getBKTLLow(energy);
   }

   else if (energy > splineMaterial->GetXmax()) {
      return getBKTLHigh(energy);
   }

   return getTLFromEnergy(energy, splineMaterial);
}

Float_t getWEPLFromEnergy(Float_t energy) {
   if (energy == 0) return 0;
   Float_t wepl = getTLFromEnergy(energy, splineWater);
/*
   if (kIsFirstLayerAir) {
      // This is an approximation to increase the accuracy of the range determination
      // Due to the lack of absorber in the first layer, and thus dz * tl does not represent the 'average' energy loss path length (but DZ does)
      // dz = DZ + absorberLength

      Float_t tl = getTLFromEnergy(energy, splineMaterial);
      wepl -= wepl/tl * (dz);
   }
  */ 
   return fmax(wepl, 0);
}

Float_t getTLFromWEPL(Float_t wepl) {
   if (wepl == 0) return 0;
   Float_t energy = getEnergyFromWEPL(wepl);
   Float_t tl = getTLFromEnergy(energy);

   return tl;
}

//////////////////////////
// Conversion to WEPL   //
//////////////////////////

Float_t getWEPLFactorFromEnergy(Float_t energy) {
   if (energy == 0) return 0;

   Float_t range = getTLFromEnergy(energy);
   Float_t wepl = getWEPLFromEnergy(energy);

   return wepl / range;
}

Float_t getWELPFactorFromWEPL(Float_t wepl) {
   Float_t range = getTLFromWEPL(wepl);

   return wepl / range;
}

Float_t getWEPLFromTL(Float_t tl) {
   if (tl == 0) return 0;
   Double_t energy = getEnergyFromTL(tl);
   Double_t wepl = getWEPLFromEnergy(energy);

   return wepl;
}


/////////////////////////////////////
// calculation of range straggling //
/////////////////////////////////////

Float_t getTLStragglingFromTL(Float_t tl, Float_t sigma_energy) {
   if (tl == 0) return 0;
   Float_t energy = getEnergyFromTL(tl);
   
   // The straggling is calculated empirically from fit from full MC data
   Float_t sigma_a = pow(straggling_a + tl * straggling_b, 2);
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
   Float_t estimatedStraggling;

   Float_t tl = getTLFromWEPL(wepl);
   Float_t tl_straggling = getTLStragglingFromTL(tl, sigma_energy);

   Float_t tl_high = tl + tl_straggling / 2;
   Float_t tl_low = tl - tl_straggling / 2;

   Float_t wepl_high = getWEPLFromTL(tl_high);
   Float_t wepl_low = getWEPLFromTL(tl_low);

   estimatedStraggling = wepl_high - wepl_low;

   return estimatedStraggling;
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

Float_t getUnitFromTL(Float_t tl) {
   Float_t res = 0;
   if       (kOutputUnit == kEnergy)   res = getEnergyFromTL(tl);
   else if  (kOutputUnit == kPhysical) res = tl;
   else if  (kOutputUnit == kWEPL)     res = getWEPLFromTL(tl);

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

   if       (kOutputUnit == kEnergy)   res = unit;
   else if  (kOutputUnit == kPhysical) res = getEnergyFromTL(unit);
   else if  (kOutputUnit == kWEPL)     res = getEnergyFromWEPL(unit);
   
   return res;
}

Float_t getEnergyAtUnit(Float_t unit, Float_t degraderthickness) {
   Float_t res = 0;

   if       (kOutputUnit == kEnergy)      res = getEnergyAtWEPL(getWEPLFromEnergy(unit), degraderthickness); // Eh?
   else if  (kOutputUnit == kWEPL)        res = getEnergyAtWEPL(unit, degraderthickness);
   else if  (kOutputUnit == kPhysical)    res = getEnergyAtTL(unit, degraderthickness);

   return res;
}

Float_t getUnitStragglingFromEnergy(Float_t energy, Float_t sigma_energy) {
   Float_t res = 0;
   if      (kOutputUnit == kEnergy)    res = getEnergyStragglingFromEnergy(energy, sigma_energy);
   else if (kOutputUnit == kWEPL)      res = getWEPLStragglingFromEnergy(energy, sigma_energy);
   else if (kOutputUnit == kPhysical)  res = getTLStragglingFromEnergy(energy, sigma_energy);

   return res;
}

Float_t getEnergyLossAtWEPL(Float_t wepl, Float_t initialEnergy) {
   Float_t sigma = 1e-3; // 1 um
   Float_t totalRange = getWEPLFromEnergy(initialEnergy);
   
   Float_t energy_low = getEnergyFromWEPL(totalRange - (wepl - sigma / 2));
   Float_t energy_high = getEnergyFromWEPL(totalRange - (wepl + sigma / 2));
   
   Float_t gradient = 1e3 * (energy_low - energy_high); // keV/um

   return gradient;
}

Float_t getEnergyLossAtTL(Float_t tl, Float_t initialEnergy) {
   Float_t sigma = 1e-2; // 0.1 um
   Float_t totalRange = getTLFromEnergy(initialEnergy);

   Float_t energy_low = getEnergyFromTL(totalRange - (tl - sigma / 2));
   Float_t energy_high = getEnergyFromTL(totalRange - (tl + sigma / 2));

   Float_t gradient = 1e2 * (energy_low - energy_high); // keV/um

   return gradient;
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
   else if  (nScintillators == 1) energyLossError = 0.30;
   else if  (nScintillators == 2) energyLossError = 0.38;
   else if  (nScintillators == 3) energyLossError = 0.44;
   
   return energyLossError;
}

Float_t getEnergyLossFromAluminumAbsorber(Float_t energy) {
   // We assume that the energy loss function is constant in the high energy range
   // (somewhat true)
   // and calculate the energy loss of kAbsorberThickness mm Aluminum by scaling the
   // measured MC energy loss in a 1.5 mm Al absorber

   return  62.00441 * pow(energy, -0.7158403) * (kAbsorberThickness/1.5);
}

Float_t getEnergyLossFromTracker(Float_t energy) {
   // NEW CALCULATION, BASED ON STOPPING POWER TABLES
   return 23.775 * pow(energy, -0.737);
}

Float_t getEnergyLossErrorFromAluminumAbsorber() {
   return 0.135;
}

Float_t getEnergyFromDegraderThickness(Double_t degraderThickness) {
   // THIS IS A ONE-TIME OPERATION, SO NO NEED TO WORRY ABOUT SPEED

   Double_t phaseSpaceDegraderthickness[1200];
   Double_t phaseSpaceEnergy[1200];
   Double_t e, es;
   Int_t    dt;
   Int_t idx = 0;
   ifstream in;
   
   if (kHelium) {
      in.open("Data/Ranges/EnergyAfterDegraderHelium.csv");
      printf("Loading EnergyAfterDegrader.csv - Helium 917 MeV\n");
      while (1) {
         in >> dt >> e >> es;
         if (!in.good()) break;
         phaseSpaceDegraderthickness[idx] = double(dt);
         phaseSpaceEnergy[idx++] = e;
      }
   }

   else if (kEnergy == 250) {
      in.open("Data/Ranges/EnergyAfterDegraderPSTAR.csv");
      printf("Loading EnergyAfterDegrader.csv - proton 250 MeV\n");
      while (1) {
         in >> dt >> e >> es;
         if (!in.good()) break;
         phaseSpaceDegraderthickness[idx] = double(dt);
         phaseSpaceEnergy[idx++] = e;
      }
   }
   else if (kEnergy == 230) {
      in.open("Data/Ranges/EnergyAfterDegraderProton.csv");
      printf("Loading EnergyAfterDegraderHelium - proton 230 MeV (78 eV)\n");
      while (1) {
         in >> dt >> e >> es;
         if (!in.good()) break;
         phaseSpaceDegraderthickness[idx] = double(dt);
         phaseSpaceEnergy[idx++] = e;
      }
   }

   in.close();

   TSpline3 *phaseSpaceSpline = new TSpline3("phaseSpaceSpline", phaseSpaceDegraderthickness, phaseSpaceEnergy, idx);

   cout << "Degraderthickness is " << degraderThickness << ". Remaining energy is " << phaseSpaceSpline->Eval(degraderThickness) << endl;

   Double_t result = phaseSpaceSpline->Eval(degraderThickness);

   delete phaseSpaceSpline;
   return (float) result;
}


// BRAGG-KLEEMAN DERIVED FALLBACK FUNCTIONS
// FOR LOWER-QUALITY CALCULATIONS OUTSIDE SPLINE REGIONS


// General functions
Float_t getBKEnergy(Float_t tl, Float_t aa, Float_t pp) {
   return pow(tl / aa, 1/pp);
}

Float_t getBKTL(Float_t energy, Float_t aa, Float_t pp) {
   return aa * pow(energy, pp);
}

// Specific functions
Float_t getBKEnergyHigh(Float_t tl) {
   return getBKEnergy(tl, alpha_material_high, p_material_high);
}

Float_t getBKEnergyLow(Float_t tl) {
   return getBKEnergy(tl, alpha_material_low, p_material_low);
}

Float_t getBKTLHigh(Float_t energy) {
   return getBKTL(energy, alpha_material_high, p_material_high);
}

Float_t getBKTLLow(Float_t energy) {
   return getBKTL(energy, alpha_material_low, p_material_low);
}

#endif
