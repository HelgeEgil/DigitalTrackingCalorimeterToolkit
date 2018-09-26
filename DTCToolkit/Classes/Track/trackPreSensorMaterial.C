#ifndef trackPreSensorMaterial_cxx
#define trackPreSensorMaterial_cxx

#include <iostream>
#include <cmath>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TClonesArray.h>
#include <TF1.h>

#include "Classes/Track/Track.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Hit/Hit.h"
#include "HelperFunctions/Tools.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"

Bool_t Track::isHitOnScintillatorH() {
   // (dx, dy, dz) = (4, 1, 0.5) cm @ -174 mm
   // ...but we get the most accurate results (confirmed with MC tagging),
   // see chapter in tracking documentation) using extrapolation - > 0
   
   Bool_t      isOnScintillator = false;
   Cluster    * extrapolatedCluster = getExtrapolatedClusterAt(0);
   Float_t     y = extrapolatedCluster->getYmm();

   if (fabs(y)<5) {
      isOnScintillator = true;
   }

   delete extrapolatedCluster;

   return isOnScintillator;
}

Bool_t Track::isHitOnScintillatorV() {
   // (dx, dy, dz) = (1, 4, 0.5) cm @ -166 mm
   // // ...but we get the most accurate results (confirmed with MC tagging),
   // see chapter in tracking documentation) using extrapolation - > 0
   
   Bool_t      isOnScintillator = false;
   Cluster    * extrapolatedCluster = getExtrapolatedClusterAt(0);
   Float_t     x = extrapolatedCluster->getXmm();

   if (fabs(x)<5) {
      isOnScintillator = true;
   }

   delete extrapolatedCluster;

   return isOnScintillator;
}

Int_t Track::getNScintillators() {
   return 1 + isHitOnScintillatorH() + isHitOnScintillatorV();

}
Bool_t Track::isHitOnScintillators() {
   return (isHitOnScintillatorH() || isHitOnScintillatorV());
}
Float_t Track::getPreEnergyLoss() {
   // Functions getEnergyLossFromScintillators and getEnergyLossFromAluminumAbsorber
   // are located in GlobalConstants/MaterialConstants.C
   
   Float_t  energyLoss = 0;
   Int_t    nScintillators = 0;

   if (kIsScintillator) {
      nScintillators = 1 + isHitOnScintillatorH() + isHitOnScintillatorV();
      energyLoss = getEnergyLossFromScintillators(run_energy, nScintillators);
   }

   if (kIsAluminumPlate) {
      energyLoss += getEnergyLossFromAluminumAbsorber(run_energy);
   }
   
   return energyLoss;
}

Float_t Track::getPreEnergyLossError() {
   Float_t  energyLossError = 0;
   Bool_t   scintillatorH = isHitOnScintillatorH();
   Bool_t   scintillatorV = isHitOnScintillatorV();
   Bool_t   scintillatorF = true; // 4x4 cm scintillator, full field

   Int_t nScintillators = scintillatorF + scintillatorH + scintillatorV;

   energyLossError = quadratureAdd(getEnergyLossErrorFromScintillators(nScintillators),
                                   getEnergyLossErrorFromAluminumAbsorber());

   return energyLossError;
}

Float_t Track::getPreTL() {
   Float_t  energyLoss = 0;
   Float_t  tl = 0;
   Int_t    nScintillators = 0;

   if (kIsScintillator) {
      nScintillators = getNScintillators();
      energyLoss = getEnergyLossFromScintillators(run_energy, nScintillators);
   }

   if (kIsAluminumPlate) {
      energyLoss += getEnergyLossFromAluminumAbsorber(run_energy);
   }

   if (kIsFirstLayerAir) {
      energyLoss += getEnergyLossFromTracker(run_energy);
   }

   if (energyLoss > 0) {   
      tl = getTLFromEnergy(run_energy) - getTLFromEnergy(run_energy - energyLoss);
   }

   return tl;
}

Float_t Track::getPreWEPL() {
   Float_t  energyLoss = 0;
   Float_t  wepl;
   Int_t    nScintillators = 0;

   if (kIsScintillator) {
      nScintillators = 1 + isHitOnScintillatorH() + isHitOnScintillatorV();
      energyLoss = getEnergyLossFromScintillators(run_energy, nScintillators);
   }

   if (kIsAluminumPlate) {
      energyLoss += getEnergyLossFromAluminumAbsorber(run_energy);
   }
   
   wepl = getWEPLFromEnergy(run_energy) - getWEPLFromEnergy(run_energy - energyLoss);
   return wepl;
}

#endif
