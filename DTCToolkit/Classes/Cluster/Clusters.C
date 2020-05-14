#ifndef clusters_c
#define clusters_c

#include <iostream>
#include <vector>

#include <TClonesArray.h>
#include <TF1.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Track/Track.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;
using namespace DTC;

Clusters::Clusters(Bool_t frameType) : clusters_("DTC::Cluster", kEventsPerRun*10),
                                       clustersWithoutTrack_("DTC::Cluster", kEventsPerRun*10) {
   frameType_ = frameType;

   clusters_.SetOwner(kTRUE);
   clustersWithoutTrack_.SetOwner(kTRUE);
   clusters_.SetBit(kCanDelete);
   clustersWithoutTrack_.SetBit(kCanDelete);
   clusters_.SetBit(kMustCleanup);
   clustersWithoutTrack_.SetBit(kMustCleanup);


   kMCSFactorFirstPass = 2; // 1
   kMCSFactorSecondPass = 3; // 4
   kMCSFactorLastPass1 = 3; // 6
   kMCSFactorLastPass2 = 5; // 6
   kMCSFactorLastPass3 = 7; // 6
}

Clusters::~Clusters() {
   // Destructor
   Clear("C+C");
   clustersWithoutTrack_.Clear("C+C");
   layerIndex_.clear();   
   clustersPerEventID_.clear();
   clusters_.Delete();
   clustersWithoutTrack_.Delete();
}

Int_t Clusters::GetEntriesFastLastLayer() {
   Int_t nInLayer[nLayers];
   Int_t lastLayer = 0;

   for (Int_t i=0; i<nLayers; i++) nInLayer[i] = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      nInLayer[getLayer(i)]++;
   }

   for (Int_t i=0; i<nLayers; i++) {
      if (nInLayer[i] > 0) lastLayer = i;
   }
   
   if (nInLayer[lastLayer] < 0.5 * nInLayer[lastLayer - 1]) {
      lastLayer--;
   }

   return nInLayer[lastLayer];
}

Int_t Clusters::GetEntriesInLayer(Int_t layer) {
   Int_t nInLayer = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (getLayer(i) == layer) nInLayer++;
   }

   return nInLayer;
}

void Clusters::Clear(Option_t *option) {
   clusters_.Clear(option);
   clustersWithoutTrack_.Clear(option);

   for (UInt_t i=0; i<layerIndex_.size(); i++) {
      layerIndex_.at(i) = - 1;
   }

   frameType_ = 0;
}

TObject * Clusters::removeClusterAt(Int_t i) {
   appendClusterWithoutTrack(At(i));
   return clusters_.RemoveAt(i);
}

void Clusters::removeCluster(Cluster * cluster) {
   appendClusterWithoutTrack(cluster);
   clusters_.Remove((TObject*) cluster);
}

void Clusters::removeAllClustersInTrack(Track *track) {
   Int_t    lastIndex = 0;
   Float_t  x, y;
   Int_t    layer;

   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      layer = track->getLayer(i);
      x = track->getX(i);
      y = track->getY(i);

      for (Int_t j=lastIndex; j<GetEntriesFast(); j++) {
         if (!At(j))
            continue;

         if (getLayer(j) < layer)
            continue;

         if (getX(j) == x) {
            if (getY(j) == y && getLayer(j) == layer) {
               removeClusterAt(j);
               lastIndex = j+1;
            }
         }
      }
   }
}

void Clusters::removeTrackFromClustersWithoutTrack(Track *track) {
   // FIXME TClonesArray error
   
   Int_t    lastIndex = 0;
   Int_t    layer;
   Float_t  x, y;
   

   showDebug("RemoveTrackFromCWT:\n");
   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      layer = track->getLayer(i);
      x = track->getX(i);
      y = track->getY(i);
      
      for (Int_t j=0; j<GetEntriesFastCWT(); j++) { // j = lastIndex -> 0
         Cluster *cluster = (Cluster *) clustersWithoutTrack_.At(j);

         if (!cluster)
            continue;
         
         if (cluster->getLayer() < layer)
            continue;

         if (cluster->getX() == x) {
            if (cluster->getY() == y && cluster->getLayer() == layer) {
               removeClusterWithoutTrackAt(j);
               lastIndex = j+1;
               showDebug("OK!\n");
               break;
            }
         }
      }
   }
}

void Clusters::removeSmallClusters(Int_t size) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (getSize(i) <= size) {
         At(i)->markUsed();
         removeClusterAt(i);
      }
   }
   Compress();
}

void Clusters::removeAllClustersAfterLayer(Int_t afterLayer) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (getLayer(i) > afterLayer) {
         removeClusterAt(i);
      }
   }
}

void Clusters::appendCluster(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID, Bool_t isSecondary, Int_t PDG) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   c->set(x,y,layer,size, eventID, isSecondary, PDG);
}

void Clusters::appendClusterEdep(Float_t x, Float_t y, Int_t layer, Float_t edep, Int_t eventID, Bool_t isSecondary, Int_t PDG) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   
   // calculate size from edep
   Int_t size = getCSFromEdep(edep);

   c->set(x,y,layer,size, eventID, isSecondary, PDG);
}

void Clusters::appendCluster(Cluster *cluster) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());

   if (cluster->isUsed()) { c->markUsed(); }
   c->setEventID(cluster->getEventID());
   c->setSecondary(cluster->isSecondary());
   c->setPDG(cluster->getPDG());
}

void Clusters::appendClusterWithoutTrack(Cluster *cluster) {
   Int_t i = clustersWithoutTrack_.GetEntriesFast();
   Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());
   c->setEventID(cluster->getEventID());
   c->setSecondary(cluster->isSecondary());
   c->setPDG(cluster->getPDG());
}


void Clusters::appendTrackToClustersWithoutTrack(Track *track) {
   Cluster * cluster = nullptr;
   for (Int_t j=0; j<=track->GetEntriesFast(); j++) {
      Int_t i = clustersWithoutTrack_.GetEntriesFast();
      Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i);
      cluster = track->At(j);

      c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());
      c->setEventID(cluster->getEventID());
      c->setSecondary(cluster->isSecondary());
      c->setPDG(cluster->getPDG());
   }
}

void Clusters::markUsedClusters(Track *track) {
   Int_t    layer, idx;
   Float_t  x, y;

   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      layer = track->getLayer(i);
      x = track->getX(i);
      y = track->getY(i);

      idx = getClusterIdx(x,y,layer);
      if (idx > -1) {
         markUsed(idx);
      }
   }
}

Int_t Clusters::getFirstIndexOfLayer(UInt_t layer) {
   if (layerIndex_.size() < layer) return -1;
   return layerIndex_.at(layer);
}

Int_t Clusters::getLastIndexOfLayer(UInt_t layer) {
   for (Int_t i=layer+1; i<nLayers; i++) {
      if (getFirstIndexOfLayer(i) != -1 ) {
         return getFirstIndexOfLayer(i);
      }
   }
   return GetEntriesFast(); // no hits beyond layer
}

Int_t Clusters::getClusterIdx(Float_t x, Float_t y, Int_t layer) {
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

Int_t Clusters::getLastActiveLayer() {
   Int_t lastActiveLayer = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      if (getLayer(i) > lastActiveLayer)
         lastActiveLayer = getLayer(i);
   }
   return lastActiveLayer;
}

void Clusters::makeLayerIndex() {
   Int_t    lastLayer = -1;
   Int_t    startOffset = 0;
   Bool_t   kStarted = false;

   if (layerIndex_.size() == 0) {
      for (Int_t i=0; i<nLayers; i++) {
         layerIndex_.push_back(-1);
      }
   }

   else {
      for (Int_t i=0; i<nLayers; i++) {
         layerIndex_.at(i) = -1;
      }
   }

   if (!GetEntriesFast()) return;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
         
      if (lastLayer != getLayer(i)) {
         layerIndex_.at(getLayer(i)) = i;
         lastLayer = getLayer(i);
      }
   }

   // set first indices to 0 if getLayer(0)>0
   for (UInt_t i=0; i<layerIndex_.size(); i++) {
      if (layerIndex_.at(i) == -1 && !kStarted) {
         startOffset = i;
      }

      else if (layerIndex_.at(i) != -1) {
         kStarted = true;
      }
   }

   if (startOffset>0) {
      for (Int_t i=0; i<=startOffset; i++) {
         layerIndex_.at(i) = 0;
      }
   }

}

void Clusters::matchWithEventIDs(Hits * eventIDs) {   
   Float_t     minDist = 1e5; // pixels
   Float_t     thisDist = 0;
   Cluster   * thisCluster = nullptr;
   Hit       * thisHit = nullptr;
   Int_t       layer = -1;
   Int_t       minIdx = -1;
   Bool_t      doLoop = true;
   Int_t       cWithoutEventID = 0;
   Float_t     cX, cY;
   Int_t       nClusters = GetEntries();
   Int_t       nHits = eventIDs->GetEntriesFast();

   for (Int_t c=0; c<GetEntriesFast(); c++) {
      thisCluster = At(c);
      if (!thisCluster) continue;

      layer = thisCluster->getLayer();

      cX = thisCluster->getX();
      cY = thisCluster->getY();

      minDist = 1e5;
      minIdx = -1;

      for (Int_t h=0; h<eventIDs->GetEntriesFast(); h++) {
         thisHit = eventIDs->At(h);
         
         if (!thisHit) continue;
         if (thisHit->getLayer() != layer) continue;

         if (fabs(cX - thisHit->getX()) < 10) {
            if (fabs(cY - thisHit->getY()) < 10) {
               thisDist = diffXY(thisCluster, thisHit);
               if (thisDist < minDist) {
                  minDist = thisDist;
                  minIdx = h;
               }
            }
         }
      }

      if (minIdx >= 0 && minDist < 10) {
         thisCluster->setEventID(eventIDs->getEventID(minIdx));
         eventIDs->removeHitAt(minIdx);
      }
   }

   for (Int_t c=0; c<GetEntriesFast(); c++) {
      if (!At(c)) continue;
      if (getEventID(c) < 0) cWithoutEventID++;
   }
}

Int_t Clusters::getClustersForEventID(Int_t eventID) {
   if (clustersPerEventID_.size() == 0) {
//      throw logic_error("Cannot use getClustersForEventID before running findNumberOfClustersForEachEventID");
      cout << "ERROR\n";
      return -1;
   }
   return clustersPerEventID_.at(eventID);
}

void Clusters::findNumberOfClustersForEachEventID() {
   // Should be 0 for eventIDs not in vector
   // DO THIS BEFORE SORTING BY LAYER!!!!

   Int_t highestEventID = Last()->getEventID();
   Int_t cIdx = 0;
   Int_t n;
   
   for (Int_t i=0; i <= highestEventID; i++) {
      n = 0;
      while (cIdx < GetEntriesFast() && i == getEventID(cIdx++)) n++;
      clustersPerEventID_.push_back(n);
   }

   Int_t sumEID = 0;
   for (Int_t i=0; i <= highestEventID; i++) {
      sumEID += clustersPerEventID_.at(i);
      printf("For history with EID %d, there are %d clusters.\n", i, clustersPerEventID_.at(i));
   }
   printf("...In total %d clusters.\n", sumEID);
}

Float_t Clusters::removeClustersInGap(Float_t gapSizemm, Float_t gapPosmm) {
   // Remove all clusters inside Y gap area

   Float_t  nClustersRemoved = 0;
   Float_t  nClustersTotal = GetEntries();
   Cluster *thisCluster = nullptr;
   Float_t  yFrom = gapPosmm - gapSizemm/2;
   Float_t  yTo = gapPosmm + gapSizemm/2;
   Float_t  y;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisCluster = At(i);
      if (!thisCluster) continue;

      y = thisCluster->getYmm();

      if (y > yFrom && y < yTo) {
         // INSIDE GAP AREA
         removeClusterAt(i);
         nClustersRemoved++;
         continue;
      }
   }

   Compress();
   return nClustersRemoved;
}


void Clusters::propagateSecondaryStatusFromTop(Int_t eventID) {
   if (eventID < 0) eventID = Last()->getEventID();
   Int_t n = 0;

   for (Int_t i=GetEntriesFast()-1; i>0; i--) {
      if (!At(i)) continue;

      if (eventID == At(i)->getEventID()) {
         At(i)->setSecondary(true);
         n++;
      }
      else break;
   }
   printf("Back-converted %d particles with eventID %d as secondaries\n", n, eventID);
}

void Clusters::removeHaloAtRadius(Float_t radius) {
   Cluster * cluster = nullptr;
   float x,y;
   Int_t lastRemovedEventID = -1;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      cluster = At(i);
      if (!cluster) continue;

      x = cluster->getXmm();
      y = cluster->getYmm();

      if (sqrt(x*x+y*y) > radius) {
         
         // Propagate any found 2nd status in this track to i-1
         if (i>0) {
            if (lastRemovedEventID != getEventID(i)) {
               for (Int_t j=i; j<i+200; j++) {
                  if (!At(j) || !At(i-1)) break;
                  if (getEventID(j) != getEventID(i)) break;
                  if (At(j)->isSecondary()) {
                     At(i-1)->setSecondary(true);
                     break;
                  }
               }
            }
         }
         lastRemovedEventID = At(i)->getEventID();
         removeClusterAt(i);
      }
   }

   Compress();
}


#endif
