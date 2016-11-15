#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <TObject.h>

#include "GlobalConstants/MaterialConstants.h"
#include "GlobalConstants/Constants.h"
#include "Classes/Hit/Hit.h"
#include "Classes/Cluster/Cluster.h"
// #include "Classes/Track/conversionFunctions.h"

using namespace std;
   
// Full revision of range-energy and depth-dose relationship as of 2016-06
// Using Ulmer Rad Phys and Chem 76 (2007) 1089-1107
// In short:
// R = a1 E0 * (1 + sum_(k=1)^2 (bk - bk * exp(-gk * E0)))
// E = R * sum_(i=1)^5 (ck * exp(-lambdak * R))
//     ( and  E(z) = (R-z) * sum_(i=1)^5 (ck * exp(-lambdak * (R-z))) )
// -dE/dz = E(z) / (R-z) - sum_(k=1)^5 (lambdak * ck * (R-z) * exp(-lambdak * (R-z)))
// The parameters a1, bk(1-2), gk(1-2), ck(1-5), lambdak(1-5) are found through range-energy fitting on different materials in GATE simulations, and defined below
// All calculations are performed in Classes/Track/conversionFunctions.C

void MaterialConstants() {

   a1_water = 0.0727504;
   b1_water = 44.7718;
   g1_water = 0.00108048;
   b2_water = 13.4883;
   g2_water = 0.00457689;

   c1_water = 10.7425;
   l1_water = 0.500803;
   c2_water = 1.7655;
   l2_water = 0.0734755;
   c3_water = 0.95865;
   l3_water = 0.0236093;
   c4_water = 0.723177;
   l4_water = 0.00662777;
   c5_water = 0.733621;
   l5_water = 0.000523239;

   cout << "LoadRangeValues 1\n";
   loadRangeValuesForTungsten();
   cout << "LoadRangeValues 2\n";
   loadRangeValuesForAluminium();

   cout << "create splines\n";
   createSplines();

   // USING VALUES FOR WATER FOR PMMA FIXME
   
   a1_pmma = 0.0727504;
   b1_pmma = 44.7718;
   g1_pmma = 0.00108048;
   b2_pmma = 13.4883;
   g2_pmma = 0.00457689;
   c1_pmma = 10.7425;
   l1_pmma = 0.500803;
   c2_pmma = 1.7655;
   l2_pmma = 0.0734755;
   c3_pmma = 0.95865;
   l3_pmma = 0.0236093;
   c4_pmma = 0.723177;
   l4_pmma = 0.0662777;
   c5_pmma = 0.733621;
   l5_pmma = 0.000523239;

   firstUpperLayerZ = 0.3;
   firstLowerLayerZ = 0.6;

   p_water = 1.7547;
   alpha_water = 0.02387;
   alpha_prime_water = 0.0087;

   proton_mass = 938.27;
   X0_firstlayer = 33.36;
   X0_pmma = 16.52;

   for (Int_t i=0; i<100; i++) {
      mcs_radius_per_layer[i] = 0;
   }

   if (kMaterial == kTungsten) {
      nLayers = 41;
      p = p_tungsten;
      alpha = alpha_tungsten;
      alpha_prime = alpha_prime_tungsten;
      X0 = X0_tungsten;
      
      a1_material = a1_tungsten;
      b1_material = b1_tungsten;
      g1_material = g1_tungsten;
      b2_material = b2_tungsten;
      g2_material = g2_tungsten;
      c1_material = c1_tungsten;
      l1_material = l1_tungsten;
      c2_material = c2_tungsten;
      l2_material = l2_tungsten;
      c3_material = c3_tungsten;
      l3_material = l3_tungsten;
      c4_material = c4_tungsten;
      l4_material = l4_tungsten;
      c5_material = c5_tungsten;
      l5_material = l5_tungsten;

      splineMaterial = splineWater; // FIX WITH UPDATED W VALUES
      splineMaterialInv = splineWaterInv; // FIX WITH UPDATED W VALUES
   }

   else if (kMaterial == kAluminium) {
      nLayers = 100;
      p = p_aluminum;
      alpha = alpha_aluminum;
      alpha_prime = alpha_prime_aluminum;
      X0 = X0_aluminum;

      a1_material = a1_aluminium;
      b1_material = b1_aluminium;
      g1_material = g1_aluminium;
      b2_material = b2_aluminium;
      g2_material = g2_aluminium;
      c1_material = c1_aluminium;
      l1_material = l1_aluminium;
      c2_material = c2_aluminium;
      l2_material = l2_aluminium;
      c3_material = c3_aluminium;
      l3_material = l3_aluminium;
      c4_material = c4_aluminium;
      l4_material = l4_aluminium;
      c5_material = c5_aluminium;
      l5_material = l5_aluminium;

      splineMaterial = spline2mmAl;
      splineMaterialInv = spline2mmAlInv;
   }

   else if (kMaterial == kPMMA) {
      nLayers = 65;
      alpha_prime = alpha_prime_water; // MUST VERIFY
      p = p_water; // MUST VERIFY
      alpha = alpha_water; // MUST VERIFY
      X0 = X0_pmma;
      
      a1_material = a1_pmma;
      b1_material = b1_pmma;
      g1_material = g1_pmma;
      b2_material = b2_pmma;
      g2_material = g2_pmma;
      c1_material = c1_pmma;
      l1_material = l1_pmma;
      c2_material = c2_pmma;
      l2_material = l2_pmma;
      c3_material = c3_pmma;
      l3_material = l3_pmma;
      c4_material = c4_pmma;
      l4_material = l4_pmma;
      c5_material = c5_pmma;
      l5_material = l5_pmma;

      splineMaterial = splineWater; // FIX
      splineMaterialInv = splineWaterInv; // FIX
   }
}

void  createSplines() {
   // load energies and ranges for water, orig geometry + AL optimized
   // create forward and backwards splines
   // store them in .h

   cout << "Creating SPLINE files\n";
   ifstream in;
   Float_t    energy;
   Int_t    idx2mmAl = 0;
   Int_t    idxWater = 0;
   Int_t    idxPureAl = 0;
   Double_t  range;
   Double_t  ranges2mmAl[250];
   Double_t  energies2mmAl[250];
   Double_t  rangesWater[250];
   Double_t  energiesWater[250];
   Double_t  rangesPureAl[250];
   Double_t  energiesPureAl[250];

   if (kAbsorbatorThickness == 2) {
      in.open("Data/Ranges/2mm_Al.csv");
   }

   else if ( kAbsorbatorThickness == 3) {
      in.open("Data/Ranges/3mm_Al.csv");
   }

   else {
      in.open("Data/Range/2mm_Al.csv");
   }

   while (1) {
      in >> energy >> range;
      if (!in.good()) break;

      ranges2mmAl[idx2mmAl] = range;
      energies2mmAl[idx2mmAl++] = energy;
   }

   in.close();

   in.open("Data/Ranges/Water.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;

      rangesWater[idxWater] = range*10; // cm to mm
      energiesWater[idxWater++] = energy;
   }
   in.close();
   
   in.open("Data/Ranges/PureAl.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;

      rangesPureAl[idxPureAl] = range*10; // cm to mm
      energiesPureAl[idxPureAl++] = energy;
   }
   

   in.close();

   spline2mmAl = new TSpline3("spline2mmAl", energies2mmAl, ranges2mmAl, idx2mmAl);
   splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);
   splinePureAl = new TSpline3("splinePureAl", energiesPureAl, rangesPureAl, idxPureAl);
   spline2mmAlInv = new TSpline3("spline2mmAlInv", ranges2mmAl, energies2mmAl, idx2mmAl);
   splineWaterInv = new TSpline3("splineWaterInv", rangesWater, energiesWater, idxWater);
   splinePureAlInv = new TSpline3("splineWaterInv", rangesPureAl, energiesPureAl, idxPureAl);

   if (kAbsorbatorThickness == 2) {
      alpha_aluminum = 0.0154651;
      p_aluminum = 1.73118;
   }

   else if (kAbsorbatorThickness == 3) {
      alpha_aluminum = 0.0173389;
      p_aluminum = 1.7751;
   }

//   alpha_aluminum = 0.0140203;
//   p_aluminum = 1.72903;
}

void loadRangeValuesForTungsten() {
   X0_tungsten = 4.2;
   p_tungsten = 1.6677;
   alpha_tungsten = 0.004461;
   alpha_prime_tungsten = 0.1086784;

   a1_tungsten = 0.0128998;
   b1_tungsten = 39.5341;
   g1_tungsten = 0.00080173;
   b2_tungsten = 6.67826;
   g2_tungsten = 0.00687635;

   c1_tungsten = 1.01309;
   l1_tungsten = 22.4927;
   c2_tungsten = 12.6053;
   l2_tungsten = 0.472777;
   c3_tungsten = 4.31023;
   l3_tungsten = 0.0361097;
   c4_tungsten = 6.1832;
   l4_tungsten = 0.122736;
   c5_tungsten = 5.42993;
   l5_tungsten = 0.00293302;
}

void loadRangeValuesForAluminium() {
   // Parameterization of X0 calculation
   Int_t    m = kAbsorbatorThickness;
   Float_t  mm_thickness;
   ifstream in;
   
   X0_aluminum = 5.4775 + 1.3109*m - 0.2608 * pow(m,2) + 0.0248 * pow(m,3) - 0.0009 * pow(m,4);
   alpha_prime_aluminum = 0.017102; // HOW DID I FIND THIS!!!

   /*
   in.open("OutputFiles/RangeEnergyParameters.csv");
   while (1) {
      in >> mm_thickness;

      if (!in.good()) continue;

      if (mm_thickness == kAbsorbatorThickness) {
         in >> alpha_aluminum >> p_aluminum;
         in >> a1_aluminium >> b1_aluminium >> g1_aluminium >> b2_aluminium >> g2_aluminium;
         in >> c1_aluminium >> c2_aluminium >> c3_aluminium >> c4_aluminium >> c5_aluminium;
         in >> l1_aluminium >> l2_aluminium >> l3_aluminium >> l4_aluminium >> l5_aluminium;
         continue; 
      }

      else if (mm_thickness == 10) {
         in >> alpha_pure_aluminum >> p_pure_aluminum;
         in >> a1_pure_aluminium >> b1_pure_aluminium >> g1_pure_aluminium >> b2_pure_aluminium >> g2_pure_aluminium;
         in >> c1_pure_aluminium >> c2_pure_aluminium >> c3_pure_aluminium >> c4_pure_aluminium >> c5_pure_aluminium;
         in >> l1_pure_aluminium >> l2_pure_aluminium >> l3_pure_aluminium >> l4_pure_aluminium >> l5_pure_aluminium; 
         break;
      }
   }
   in.close();
   */
}


Double_t getLayerPositionmm(Int_t i) {
   Double_t z = 0;

   if (i>0) {
      z = ( firstUpperLayerZ + firstLowerLayerZ ) / 2 + dz * i;
   }

   return z;
}

Float_t getSigmaEnergy(Int_t energy) { 
   Float_t sigma_energy = 0;
   
   if      (energy == 188) sigma_energy = 0.3; // approx
   else if (energy == 180) sigma_energy = 2.3;
   else if (energy == 170) sigma_energy = 4.5; // was 4.5
   else if (energy == 160) sigma_energy = 8.3;
   else if (energy == 150) sigma_energy = 20;
   
   return sigma_energy;
}

Bool_t isChipLowResistivity(Int_t chipIdx) {
   Bool_t isHigh = false;

   if (chipIdx == 2  || chipIdx == 3  ||
       chipIdx == 8  || chipIdx == 9  ||
       chipIdx == 16 || chipIdx == 17 ||
       chipIdx == 18 || chipIdx == 21 ||
       chipIdx == 23 || chipIdx == 24 ||
       chipIdx == 27 || chipIdx == 25) {
      isHigh = true;
   }

   return isHigh;
}

Float_t getChipCalibrationFactor(Int_t chip) {
   Float_t f[28] = {1.273, 1.253, 0.836, 0.834,
                    1.267, 1.267, 1.656, 1.417,
                    0.813, 0.835, 1.180, 0.000,
                    1.403, 1.412, 1.328, 1.366,
                    0.874, 0.829, 1.020, 1.221,
                    1.188, 1.510, 1.261, 1.316,
                    1.112, 0.976, 1.097, 1.013};
   
   if (chip > 27) {
      cout << "CHIP NUMBER TOO HIGH, CANNOT CALIBRATE!\n";
      return 0;
   }

   return f[chip];
}

Bool_t isBadData(Cluster *estimatedPosition) {
   Float_t   x = estimatedPosition->getXmm();
   Float_t   y = estimatedPosition->getYmm();
   Int_t     layer = estimatedPosition->getLayer();
   Bool_t    isBad = false;

   // gap between the chips
   if (fabs(y) < 0.05) isBad = true;

   switch (layer) {
      case 1:
         if (x<0) {
            if (y > 4.7 && y < 9.7) isBad = true;
            else if (y > -2.9 && y < 0.0) isBad = true;
            else if (y < -15) isBad = true;
         }
         else {
            if (y > 0 && y < 5.2) isBad = true;
         }
         break;

      case 2:
         if (x < 0 && y >  0 ) isBad = true;
         else if (x > 0 && y < -14.7) isBad = true;
         break;

      case 3:
         if (x < 0 && y < -14.7) isBad = true;
         break;

      case 5:
         if (x > 0 && y > 4.7 && y < 9.7) isBad = true;
         else if (x > 0 && y < -4.7 && y > -9.7) isBad = true;
         break;

      case 6:
         if (x < 0 && y > 14.7) isBad = true;
         else if (x > 0 && y > -9.7 && y < -4.7) isBad = true;
         break;

      case 7:
         if (x < 0) {
            if (y < -4.8 && y > -9.8) isBad = true;
            else if (y < -14.8) isBad = true;
         }
         break;
   }
   
   return isBad;
}
