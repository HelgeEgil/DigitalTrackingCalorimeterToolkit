#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

#include <TObject.h>

#include "GlobalConstants/MaterialConstants.h"
#include "GlobalConstants/Constants.h"
#include "Classes/Cluster/Cluster.h"

using namespace std;
using namespace DTC;

void MaterialConstants() {

   X0_tungsten = 4.2;
   p_tungsten = 1.6677;
   alpha_tungsten = 0.004461;
   alpha_prime_tungsten = 0.1086784;
   
   Float_t m = kAbsorbatorThickness;
   X0_aluminum = 5.4775 + 1.3109*m - 0.2608 * pow(m,2) + 0.0248 * pow(m,3) - 0.0009 * pow(m,4);
   alpha_aluminum = 0.0140203;
   p_aluminum = 1.72903;
   alpha_prime_aluminum = 0.017102; // HOW DID I FIND THIS!!!

   p_water = 1.7547;
   alpha_water = 0.02387;
   alpha_prime_water = 0.0087;
   X0_pmma = 16.52;
   X0_firstlayer = 33.36;
   
   straggling_a = 1.8568;
   straggling_b = 0.000856;
   
   createSplines();

   firstUpperLayerZ = 0.3;
   firstLowerLayerZ = 0.6;

   proton_mass = 938.27;

   for (Int_t i=0; i<100; i++) {
      mcs_radius_per_layer[i] = 0;
   }

   if (kMaterial == kTungsten) {
      nLayers = 41;
      p = p_tungsten;
      alpha = alpha_tungsten;
      alpha_prime = alpha_prime_tungsten;
      X0 = X0_tungsten;
      
      splineMaterial = splineW;
      splineMaterialInv = splineWInv;
   }

   else if (kMaterial == kAluminum) {
      nLayers = 100;
      p = p_aluminum;
      alpha = alpha_aluminum;
      alpha_prime = alpha_prime_aluminum;
      X0 = X0_aluminum;

      splineMaterial = spline2mmAl;
      splineMaterialInv = spline2mmAlInv;
   }

   else if (kMaterial == kPMMA) {
      nLayers = 65;
      alpha_prime = alpha_prime_water; // MUST VERIFY
      p = p_water; // MUST VERIFY
      alpha = alpha_water; // MUST VERIFY
      X0 = X0_pmma;
      
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
   Int_t    idxW = 0;
   Double_t  range;
   Double_t  ranges2mmAl[500];
   Double_t  energies2mmAl[500];
   Double_t  rangesW[500];
   Double_t  energiesW[500];
   Double_t  rangesWater[500];
   Double_t  energiesWater[500];
   Double_t  rangesPureAl[500];
   Double_t  energiesPureAl[500];
   
   if (kAbsorbatorThickness == 2) {
      in.open("Data/Ranges/2mm_Al.csv");
   }

   else if (kAbsorbatorThickness == 3) {
      in.open("Data/Ranges/3mm_Al.csv");
   }

   else if (kAbsorbatorThickness == 4) {
      in.open("Data/Ranges/4mm_Al.csv");
   }
   
   else if (kAbsorbatorThickness == 5) {
      in.open("Data/Ranges/5mm_Al.csv");
   }
   
   else if (kAbsorbatorThickness == 6) {
      in.open("Data/Ranges/6mm_Al.csv");
   }
   
   else if (kAbsorbatorThickness == 7) {
      in.open("Data/Ranges/7mm_Al.csv");
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

      rangesWater[idxWater] = range;
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

   in.open("Data/Ranges/3mm_W.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;

      rangesW[idxW] = range;
      energiesW[idxW++] = energy;
   }

   in.close();

   spline2mmAl = new TSpline3("spline2mmAl", energies2mmAl, ranges2mmAl, idx2mmAl);
   splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);
   splinePureAl = new TSpline3("splinePureAl", energiesPureAl, rangesPureAl, idxPureAl);
   splineW = new TSpline3("splineW", energiesW, rangesW, idxW);
   spline2mmAlInv = new TSpline3("spline2mmAlInv", ranges2mmAl, energies2mmAl, idx2mmAl);
   splineWaterInv = new TSpline3("splineWaterInv", rangesWater, energiesWater, idxWater);
   splinePureAlInv = new TSpline3("splineWaterInv", rangesPureAl, energiesPureAl, idxPureAl);
   splineWInv = new TSpline3("splineWInv", rangesW, energiesW, idxW);

   if (kAbsorbatorThickness == 2) { // updated 2017-03-15
      alpha_aluminum = 0.012626;
      p_aluminum = 1.728743;
      straggling_a = 1.72323;
      straggling_b = 0.00100124;
   }

   else if (kAbsorbatorThickness == 3) { // updated 2017-03-29
      alpha_aluminum = 0.011154;
      p_aluminum = 1.751128;
      straggling_a = 1.80564;
      straggling_b = 0.001971;
   }

   else if (kAbsorbatorThickness == 4) { // updated 2017-03-22
      alpha_aluminum = 0.016493;
      p_aluminum = 1.683775;
      straggling_a = 1.73913;
      straggling_b = 0.00109215;
   }
   
   else if (kAbsorbatorThickness == 5) { // updated 2017-03-23
      alpha_aluminum = 0.018807;
      p_aluminum = 1.661016;
      straggling_a = 1.71864;
      straggling_b = 0.001381;
   }

   else if (kAbsorbatorThickness == 6) { // dummy values
      alpha_aluminum = 0.021480;
      p_aluminum = 1.637919;
      straggling_a = 1.71891;
      straggling_b = 0.001397;
   }
   
   else if (kAbsorbatorThickness == 7) { // dummy values
      alpha_aluminum = 0.024070;
      p_aluminum = 1.618193;
      straggling_a = 1.73329;
      straggling_b = 0.001268;
   }
}

Double_t getLayerPositionmm(Int_t i) {
   Double_t z = 0;

   if (!kUseAlpide) {
      if (i>0) {
         z = ( firstUpperLayerZ + firstLowerLayerZ ) / 2 + dz * i;
      }
   }
   else {
      z = dz * i;
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
