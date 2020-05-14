#ifndef trackProperties_cxx
#define trackProperties_cxx

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

Int_t Track::getFirstLayer() { 
   for (Int_t i = 0; i < GetEntriesFast(); i++) {
      if (!At(i)) continue;
      return getLayer(i);
   }
   return -1;
}

Int_t Track::getLastLayer() {
   for (Int_t i = GetEntriesFast() - 1; i >= 0; i--) {
      if (!At(i)) continue;
      return getLayer(i);
   }
   return -1;
}

Bool_t Track::hasLayer(Int_t layer) {
   for (Int_t i = 0; i < GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (getLayer(i) == layer) return true;
   }
   return false;
}

Int_t Track::getIdxFromLayer(Int_t layer) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (getLayer(i) == layer) return i;
   }

   return -1;
}

Bool_t Track::doesTrackEndAbruptly() {
//   This is limit has been found through MC truths
//   See full information in PhD thesis 
   if (GetEntriesFast() < 3) return 1;

   Bool_t   endsAbruptly = (Last()->getDepositedEnergy() < 3.5);
   if (kDataType == kData) endsAbruptly = (Last()->getDepositedEnergy() < 4);
   if (kHelium) endsAbruptly = (Last()->getDepositedEnergy() < 8);

   return endsAbruptly; 
}

Float_t Track::getRiseFactor(Int_t n) {
   Float_t  riseFactor = 0;

   Int_t last = GetEntriesFast();
   Cluster *c = nullptr;

   Float_t sumRise = 0;
   Float_t sumEdep = 0;
   Float_t nEdep = 0;
   Float_t lastEdep = 0;

   for (Int_t i=1; i<=n; i++) {
      c = At(last-i);
      if (c) {
         sumEdep = c->getDepositedEnergy();
         nEdep++;
      }
   }

   Float_t mean = sumEdep / nEdep;

   for (Int_t i=1; i<=n; i++) {
      c = At(last-i);
      if (c) {
         if (lastEdep > 0) {
            sumRise += pow(mean - c->getDepositedEnergy(), 2);
         }
         lastEdep = c->getDepositedEnergy();
      }
   }

   Float_t sigma = sqrt(sumRise / nEdep);

//   printf("SIGMA POL0 = %.2f\n", sigma);

   return mean;
}

Float_t Track::getAverageDepositedEnergy(Int_t fromIdx, Int_t toIdxExclusive) { // exclusive
   Int_t nClusters = 0;
   Float_t averageEdep = 0;
   Cluster * thisCluster = nullptr;

   if (fromIdx < 0) fromIdx = 0;
   if (toIdxExclusive > GetEntriesFast() || toIdxExclusive < 0) toIdxExclusive = GetEntriesFast();

   for (Int_t i=fromIdx; i<toIdxExclusive; i++) {
      if (At(i)) {
         averageEdep += At(i)->getDepositedEnergy();
         nClusters++;
      }
   }

   averageEdep /= nClusters;

   return averageEdep;
}


Int_t Track::getNMissingLayers() {
   Int_t missingLayers = 0;
   Int_t lastLayer = 0;
   Int_t firstLayer = 1e6;
   Int_t previousLayer = -1;
   Int_t diff = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (At(i)) {
         firstLayer = i;
         break;
      }
   }

   for (Int_t i=GetEntriesFast()-1; i>=0; i--) {
      if (At(i)) {
         lastLayer = i;
         break;
      }   
   }

   previousLayer = firstLayer - 1;
   diff = 0;
  
   for (Int_t i=firstLayer; i<=lastLayer; i++) {
      if (!At(i)) continue;
    
      diff = getLayer(i) - previousLayer - 1;
      if (diff>0) missingLayers += diff;
      previousLayer = getLayer(i);
  }
  
  return missingLayers;
}

#endif
