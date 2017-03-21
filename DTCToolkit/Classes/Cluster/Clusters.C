#include <iostream>
#include <vector>

#include <TClonesArray.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Track/Track.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;
using namespace DTC;

Clusters::Clusters(Bool_t frameType) : clusters_("Cluster", kEventsPerRun*2),
                                       clustersWithoutTrack_("Cluster", kEventsPerRun*2) {
   frameType_ = frameType;
}

Clusters::~Clusters() {
   // Destructor
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

void Clusters::clearClusters() {
    clusters_.Clear("C");
    clustersWithoutTrack_.Clear("C");

    for (UInt_t i=0; i<layerIndex_.size(); i++) {
       layerIndex_.at(i) = - 1;
    }

    frameType_ = 0;
}

void Clusters::Clear(Option_t *) {
   clusters_.Clear("C");
   clustersWithoutTrack_.Clear("C");
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

   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      layer = track->getLayer(i);
      x = track->getX(i);
      y = track->getY(i);

      for (Int_t j=lastIndex; j<GetEntriesFastCWT(); j++) {
         Cluster *cluster = (Cluster *) clustersWithoutTrack_.At(j);

         if (!cluster)
            continue;
         
         if (cluster->getLayer() < layer)
            continue;

         if (cluster->getX() == x) {
            if (cluster->getY() == y && cluster->getLayer() == layer) {
               removeClusterWithoutTrackAt(j);
               lastIndex = j+1;
               break;
            }
         }
      }
   }
}

void Clusters::removeSmallClusters(Int_t size) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (getSize(i) <= size) {
         removeClusterAt(i);
      }
   }
}

void Clusters::removeAllClustersAfterLayer(Int_t afterLayer) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (getLayer(i) > afterLayer) {
         removeClusterAt(i);
      }
   }
}

void Clusters::appendCluster(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   c->set(x,y,layer,size, eventID);
}

void Clusters::appendClusterEdep(Float_t x, Float_t y, Int_t layer, Float_t edep, Int_t eventID) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   
   // calculate size from edep
   
   Int_t size = getCSFromEdep(edep);

   c->set(x,y,layer,size, eventID);
}

void Clusters::appendCluster(Cluster *cluster) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());

   if (cluster->isUsed()) { c->markUsed(); }
   c->setEventID(cluster->getEventID());
}

void Clusters::appendClusterWithoutTrack(Cluster *cluster) {
   Int_t i = clustersWithoutTrack_.GetEntriesFast();
   Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());
   c->setEventID(cluster->getEventID());
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

#ifdef USEDEBUG
   printf("Layer Index -------\n");
   for (UInt_t i=0; i<layerIndex_.size(); i++) {
      printf("layer %d: Layer Index = %d.\n", i, layerIndex_.at(i));
   }
#endif

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

// cout << "Number of clusters without eventID: " << cWithoutEventID << " (" << (float) cWithoutEventID / nClusters * 100 << "%)" << endl;
}
