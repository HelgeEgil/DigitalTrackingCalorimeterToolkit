#ifndef materialconstants_cxx
#define materialconstants_cxx

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <list>

#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeIndex.h>
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
   
   Float_t m = kAbsorberThickness;
   X0_aluminum = 5.4775 + 1.3109*m - 0.2608 * pow(m,2) + 0.0248 * pow(m,3) - 0.0009 * pow(m,4);
   alpha_aluminum = 0.0140203;
   p_aluminum = 1.72903;
   alpha_prime_aluminum = 0.017102; // HOW DID I FIND THIS!!!

   p_water = 1.7547;
   alpha_water = 0.02387;
   alpha_prime_water = 0.0087;
   X0_pmma = 16.52;
   X0_firstlayer = 33.36;
   
   createSplines();

   firstUpperLayerZ = 0.3;
   firstLowerLayerZ = 0.6;

   proton_mass = 938.27;

   for (Int_t i=0; i<100; i++) {
      mcs_radius_per_layer[i] = 0;
   }

   if (kMaterial == kTungsten) {
      nLayers = 24;
      p = p_tungsten;
      alpha = alpha_tungsten;
      alpha_prime = alpha_prime_tungsten;
      X0 = X0_tungsten;
      
      splineMaterial = splineW;
      splineMaterialInv = splineWInv;
   }

   else if (kMaterial == kAluminum) {
      if (kAbsorberThickness == 3) nLayers = 75;
      if (kAbsorberThickness == 2) nLayers = 75;
      else { nLayers = 70; };
      p = p_aluminum;
      alpha = alpha_aluminum;
      alpha_prime = alpha_prime_aluminum;
      X0 = X0_aluminum;

      splineMaterial = splineDTC;
      splineMaterialInv = splineDTCInv;
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

   else if (kMaterial == kCarbon) {
      nLayers = 100;
      p = p_aluminum;
      alpha = alpha_aluminum;
      alpha_prime = alpha_prime_aluminum;
      X0 = X0_aluminum;

      splineMaterial = splineDTC;
      splineMaterialInv = splineDTCInv;
   }

}

void  createSplines() {
   // load energies and ranges for water, orig geometry + AL optimized
   // create forward and backwards splines
   // store them in .h

   ifstream in;
   Double_t  energy;
   Int_t    idxDTC = 0;
   Int_t    idxWater = 0;
   Int_t    idxPureAl = 0;
   Int_t    idxW = 0;
   Int_t    layer = 0;
   Double_t range, wet, wepl;
   Double_t mu, sigma;
   Double_t rangesDTC[1200];
   Double_t energiesDTC[1200];
   Double_t rangesW[500];
   Double_t energiesW[500];
   Double_t rangesWater[1200];
   Double_t energiesWater[1200];
   Double_t rangesPureAl[500];
   Double_t energiesPureAl[500];

   Double_t layerNumber[50];
   Double_t layerWET[50];
   Double_t layerWEPL[50];

   Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;

   if (kHelium) {
      if (!kFinalDesign) {
         in.open("Data/Ranges/35mm_Al_Helium.csv"); // for 917 MeV
         if (readoutAbsorber != 35) {
            printf("REMEMBER TO ADD CALIBRATION FOR THIS ABSORBER THICKNESS (using 3.5 mm)!!! (%.1f)\n", kAbsorberThickness);
         }
         if (kEnergy != 917) printf("SYSTEM NOT CALIBRATED FOR %d MeV!!\n", kEnergy);
      }
      else {
         in.open("Data/Ranges/Final_Al_Helium.csv"); // for 917 MeV
         if (kEnergy != 917) printf("SYSTEM NOT CALIBRATED FOR %d MeV!!\n", kEnergy);
        
      }
   }

   else if (kEnergy == 250) { // the energy is off anywho
      if       (kMaterial == kTungsten) {
         in.open(Form("Data/Ranges/%.0fmm_W.csv", readoutAbsorber));
      }
      else if  (kMaterial == kAluminum) {
         in.open(Form("Data/Ranges/%.0fmm_Al.csv", readoutAbsorber));
      }
      else if  (kMaterial == kCarbon) {
         in.open(Form("Data/Ranges/%.0fmm_C.csv", readoutAbsorber));
      }
   }
   else if (kEnergy == 230) {
      if       (kMaterial == kTungsten) {
         in.open(Form("Data/Ranges/%.0fmm_W_230MeV.csv", readoutAbsorber));
      }
      else if  (kMaterial == kAluminum) {
         if (!kFinalDesign) {
            in.open(Form("Data/Ranges/%.0fmm_Al_230MeV.csv", readoutAbsorber));
         }
         else {
            in.open("Data/Ranges/Final_Al_Proton.csv"); // 230 MeV
         }
      }
      else if  (kMaterial == kCarbon) {
         in.open(Form("Data/Ranges/%.0fmm_C_230MeV.csv", readoutAbsorber));
      }
   }

   while (1) {
      in >> energy >> range;
      if (!in.good()) break;

      rangesDTC[idxDTC] = range;
      energiesDTC[idxDTC++] = energy;
   }

   in.close();

   if (!kHelium) {
      in.open("Data/Ranges/WaterPSTAR.csv");
      while (1) {
         in >> energy >> range;
         if (!in.good()) break;
         
         rangesWater[idxWater] = range*10;
         energiesWater[idxWater++] = energy;
      }
      in.close();
   }

   else {
      in.open("Data/Ranges/WaterASTAR.csv");
      while (1) {
         in >> energy >> range;
         if (!in.good()) break;

         rangesWater[idxWater] = range*10;
         energiesWater[idxWater++] = energy;
      }
      in.close();
   }

   in.open("Data/Ranges/PureAl.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;
      
      rangesPureAl[idxPureAl] = range*10; // cm to mm
      energiesPureAl[idxPureAl++] = energy;
   }

   in.close();

   in.open("Data/Ranges/33mm_W.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;
      
      rangesW[idxW] = range;
      energiesW[idxW++] = energy;
   }

   in.close();

   in.open(Form("Data/Scattering/%.0fmmAl_50mm.csv", kAbsorberThickness));
   while (1) {
      in >> layer >> mu >> sigma;
      if (!in.good()) break;
      mcs_radius_per_layer_empirical[layer] = mu + 3 * sigma;
   }
   in.close();

   Int_t layerIdx = 0;
   if (!kHelium) {
      in.open("Data/Ranges/layer_wet_wepl_Proton.csv");
   }
   else {
      in.open("Data/Ranges/layer_wet_wepl_Helium.csv");
   }

   while (1) {
      in >> layer >> wet >> wepl;
      if (!in.good()) break;
      layerNumber[layerIdx] = layer;
      layerWET[layerIdx] = wet;
      layerWEPL[layerIdx++] = wepl;
   }

   splineDTC = new TSpline3("splineDTC", energiesDTC, rangesDTC, idxDTC);
   splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);
   splinePureAl = new TSpline3("splinePureAl", energiesPureAl, rangesPureAl, idxPureAl);
   splineW = new TSpline3("splineW", energiesW, rangesW, idxW);
   splineDTCInv = new TSpline3("splineDTCInv", rangesDTC, energiesDTC, idxDTC);
   splineWaterInv = new TSpline3("splineWaterInv", rangesWater, energiesWater, idxWater);
   splinePureAlInv = new TSpline3("splineWaterInv", rangesPureAl, energiesPureAl, idxPureAl);
   splineWInv = new TSpline3("splineWInv", rangesW, energiesW, idxW);

   splineWET = new TSpline3("splineWET", layerNumber, layerWET, layerIdx);
   splineWEPL = new TSpline3("splineWET", layerNumber, layerWEPL, layerIdx);
   Double_t layerWEPLinv[50], layerWETinv[50];
   for (Int_t i=0; i<50; i++) {
      if (i>=layerIdx) continue;
      layerWEPLinv[i] = layerWEPL[layerIdx-i-1];
      layerWETinv[i] = layerWET[layerIdx-i-1];
   }

   splineWETFromDegrader = new TSpline3("splineWETFromDegrader", layerWEPLinv, layerWETinv, layerIdx);

   // FIND BRAGG-KLEEMAN PARAMETERS
   TGraph * range_energy = new TGraph(idxDTC, energiesDTC, rangesDTC);
   TF1    * range_energy_fit = new TF1("range_energy_fit", "[0] * pow(x, [1])");
   range_energy_fit->SetParameters(0.02, 1.6);
   range_energy->Fit("range_energy_fit", "Q,M");
   alpha_aluminum = range_energy_fit->GetParameter(0);
   p_aluminum = range_energy_fit->GetParameter(1);

//   printf("alpha = %.4f, p = %.4f.\n", alpha_aluminum, p_aluminum);
   
   // FIND BRAGG-KLEEMAN PARAMETERS
   TGraph * range_energyW = new TGraph(idxW, energiesW, rangesW);
   TF1    * range_energyW_fit = new TF1("range_energyW_fit", "[0] * pow(x, [1])");
   range_energyW_fit->SetParameters(0.002, 1.7);
   range_energyW->Fit("range_energyW_fit", "Q,M,N");
   alpha_tungsten = range_energyW_fit->GetParameter(0);
   p_tungsten = range_energyW_fit->GetParameter(1);

   // FIND BRAGG-KLEEMAN PARAMETERS HIGH / LOW
   if (kMaterial == kAluminum) {
      range_energy->Fit("range_energy_fit", "Q,M,N", "", 0, 50); // fit 0 - 40 MeV
      alpha_material_low = range_energy_fit->GetParameter(0);
      p_material_low = range_energy_fit->GetParameter(1);

      range_energy->Fit("range_energy_fit", "Q,M,N", "", 220, 250); // fit 220 - 250 MeV
      alpha_material_high = range_energy_fit->GetParameter(0);
      p_material_high = range_energy_fit->GetParameter(1);
   }

   else if (kMaterial == kTungsten) {
      // FIND BRAGG-KLEEMAN PARAMETERS HIGH / LOW
      range_energyW->Fit("range_energyW_fit", "Q,M,N", "", 0, 40); // fit 0 - 40 MeV
      alpha_material_low = range_energyW_fit->GetParameter(0);
      p_material_low = range_energyW_fit->GetParameter(1);

      range_energyW->Fit("range_energyW_fit", "Q,M,N", "", 220, 250); // fit 220 - 250 MeV
      alpha_material_high = range_energyW_fit->GetParameter(0);
      p_material_high = range_energyW_fit->GetParameter(1);
   }

   else { printf("COULD NOT FIND HIGH / LOW PARAMETERS FOR MATERIAL %d!!", kMaterial); }

   
   // About the same in all geometries
   straggling_a = 1.76;
   straggling_b = 0.0012;

   delete range_energy;
   delete range_energy_fit;
   delete range_energyW;
   delete range_energyW_fit;

   if (kDoDiffusion) {
      Int_t lastClusterSize = -1;
      CDB_fCluster = new TFile("Data/ClusterSizes/ALPIDE/database_final_reduced.root", "READ");
      CDB_treeCluster = (TTree*) CDB_fCluster->Get("database");
      CDB_treeCluster->SetBranchAddress("size", &CDB_clusterSize);
      CDB_treeCluster->SetBranchAddress("x_mean", &CDB_x_mean);
      CDB_treeCluster->SetBranchAddress("y_mean", &CDB_y_mean);
      CDB_treeCluster->SetBranchAddress("hit_array", &CDB_hit_array, &CDB_b_hit_array);

      for (int i=0; i<50; i++) CDB_sortIndex[i] = -1;
      
      // Load pre-sorted size index 
      in.open("Data/ClusterSizes/ALPIDE/sortIndex.csv");
      Int_t clustersize_, index_;
      while (1) {
         in >> clustersize_ >> index_;
         if (!in.good()) break;
         CDB_sortIndex[clustersize_] = index_;
      }
      in.close();
   }
}

Double_t getLayerPositionmm(Int_t i) {
   Double_t z = 0;
   if (kFinalDesign) {
      if (i < 2) {
         z = 0.2315 + i * dz2; // ALSO OK if i < 0
      }
      else {
         z = 2 * dz2 + 3.749 + (i-2) * dz;
      }

   }
   else {
      z = i * dz;
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

#endif
