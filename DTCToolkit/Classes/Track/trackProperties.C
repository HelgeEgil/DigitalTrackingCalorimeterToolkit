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
   Float_t  expectedTL = getTLFromEnergy(run_energy);
   Float_t  actualTL = getTLFromEnergy(getEnergy());
   Float_t  riseFactor = getRiseFactor();
   Float_t edepLimit = 2;
   if (kMaterial == kAluminum) edepLimit = 3;
   Bool_t   endsAbruptly = (Last()->getDepositedEnergy() < edepLimit);
   
   if (endsAbruptly) return true;
   else              return false;
}

Float_t Track::getRiseFactor() {
   Float_t  normalization = 0;
   Float_t  riseFactor = 0;
   Int_t    nNorm = 0;
   Int_t    nNormActual = 0;
   Int_t    n = GetEntriesFast();

   if (!n) return 0;
   if (n-2<0) return 0;

   nNorm = n/2;

   for (Int_t i=0; i<nNorm; i++) {
      if (!At(i)) continue;
      normalization += getSize(i);
      nNormActual++;
   }

   normalization /= nNormActual;

   if (At(n-1) && At(n-2)) {
      riseFactor = (getSize(n-1) + getSize(n-2)) / 2 / normalization;
   }

   else if (At(n-1) && !At(n-1)) {
      riseFactor = getSize(n-1) / normalization;
   }

   else if (!At(n-1) && At(n-1)) {
      riseFactor = getSize(n-2) / normalization;
   }

   else return 0;

   return riseFactor;
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

