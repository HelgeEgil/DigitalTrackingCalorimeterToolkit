#ifndef trackClusterProperties_cxx
#define trackClusterProperties_cxx

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

Int_t Track::getClusterIdx(Float_t x, Float_t y, Int_t layer) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (getLayer(i) < layer) continue;
      if (x == getX(i)) {
         if (y == getY(i)) {
            if (layer == getLayer(i)) { // found it
               return i;
            }
         }
      }
   }

   return -1;
}

Int_t Track::getClusterIdx(Cluster * cluster) {
   if (!cluster) return -1;

   return getClusterIdx(cluster->getX(), cluster->getY(), cluster->getLayer());
}

Bool_t Track::isClusterInTrack(Cluster * cluster) {
   Int_t idx = getClusterIdx(cluster);
   
   return (idx >= 0);
}

Int_t Track::getClusterFromLayer(Int_t layer) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      if (getLayer(i) == layer) {
         return i;
      }
   }

   return -1;
}

Clusters * Track::getConflictClusters() { 
   Clusters *conflictClusters = new Clusters(GetEntriesFast());

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      if (At(i)->isUsed()) {
         conflictClusters->appendCluster(At(i));
      }
   }

   return conflictClusters;
}

Int_t Track::getNumberOfConflictClusters() { 
   Int_t n = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      if (At(i)->isUsed()) n++;
   }

   return n;
}

Bool_t Track::isUsedClustersInTrack() {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      if (isUsed(i)) return true;
   }

   return false;
}

Float_t Track::getAverageCS() {
   Int_t    n = GetEntries();
   Float_t  sum = 0;
   Float_t  avg = 0;

   if (!n) return 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      sum += getSize(i);
   }

   avg = sum / n;

   return avg;
}

Float_t Track::getAverageCSLastN(Int_t last_n) {
   Int_t    n = 0;
   Int_t    nFound = 0;
   Int_t    lastIdx = 0;
   Float_t  sum = 0;
   Float_t  avg = 0;

   for (Int_t i=GetEntriesFast()-1; i>=0; i--) {
      if (!At(i)) continue;
      nFound++;
      if (nFound == last_n) {
         lastIdx = i;
         break;
      }
   }

   for (Int_t i=lastIdx; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      n++;
      sum += getSize(i);
   }

   avg = sum / n;

   return avg;
}

Float_t Track::getMeanSizeToIdx(Int_t toIdx) {
   Int_t tempSum = 0;

   if (toIdx<0 || toIdx >= GetEntriesFast())
      toIdx = GetEntriesFast() - 1;

   for (Int_t i=0; i<toIdx; i++) {
      tempSum += getSize(i);
   }

   return (float) tempSum / (toIdx + 1 - getNMissingLayers());
}

Float_t Track::getStdSizeToIdx(Int_t toIdx) {
   Float_t  tempSum = 0;
   Float_t  std;
   Int_t    mean;

   if (toIdx<0 || toIdx >= GetEntriesFast())
      toIdx = GetEntriesFast() - 1;

   mean = getMeanSizeToIdx(toIdx);

   for (Int_t i=0; i<toIdx; i++) {
      tempSum += pow(getSize(i) - mean, 2);
   }

   std = sqrt((float) tempSum / (toIdx + 1));
   return std;
}

#endif
