#include "Clusters.h"
#include "Cluster.h"
#include "Track.h"
#include "Constants.h"
#include "MaterialConstants.h"
#include "Tools.h"
#include <iostream>
#include <TClonesArray.h>
#include <vector>

using namespace std;

// ClassImp(ClusterCollection)

Clusters::Clusters(Bool_t frameType) : clusters_("Cluster", kEventsPerRun*2),
											 clustersWithoutTrack_("Cluster", kEventsPerRun*2) {
	frameType_ = frameType;
}

Clusters::~Clusters() {
   // Destructor
   clusters_.Delete();
   clustersWithoutTrack_.Delete();
}

void Clusters::appendCluster(Float_t x, Float_t y, Int_t layer, Int_t size) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   c->set(x,y,layer,size);
}

void Clusters::appendCluster(Cluster *cluster) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clusters_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());

   if (kEventsPerRun == 1) {
   	c->setEventID(cluster->getEventID());
   }

}

void Clusters::appendClusterWithoutTrack(Cluster *cluster) {
   Int_t i = clustersWithoutTrack_.GetEntriesFast();
   Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());
}

void Clusters::Clear(Option_t *) {
	clusters_.Clear("C");
	clustersWithoutTrack_.Clear("C");
}

void Clusters::clearClusters() {
	 clusters_.Clear("C");
	 clustersWithoutTrack_.Clear("C");

	 for (UInt_t i=0; i<layerIndex_.size(); i++) {
		 layerIndex_.at(i) = - 1;
	 }
	 frameType_ = 0;
}

void Clusters::removeAllClustersInTrack(Track *track) {
   Int_t lastIndex = 0;

   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      Int_t layer = track->getLayer(i);
      Float_t x = track->getX(i);
      Float_t y = track->getY(i);

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

Bool_t Clusters::removeClusterFromCoords(Float_t x, Float_t y, Int_t layer) {
   for (Int_t i=0; GetEntriesFast(); i++) {
      if (getLayer(i) < layer) continue;
      if (x == getX(i)) {
         if (y == getY(i)) {
            if (layer == getLayer(i)) { // found it
               removeClusterAt(i);
               return true;
            }
         }
      }
   }
   return false;
}

void Clusters::makeLayerIndex() {
   if (layerIndex_.size() == 0) {
      for (Int_t i=0; i<nLayers; i++) layerIndex_.push_back(-1);
   }
   else {
      for (Int_t i=0; i<nLayers; i++) layerIndex_.at(i) = -1;
   }

   if (!GetEntriesFast()) return;

	Int_t lastLayer = -1;
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i))
			continue;

		if (lastLayer != getLayer(i)) {
			layerIndex_.at(getLayer(i)) = i;
			lastLayer = getLayer(i);
		}
	}

	// set first indices to  0 if getLayer(0)>0
	Int_t startOffset = 0;
	Bool_t kStarted = kFALSE;
	for (UInt_t i=0; i<layerIndex_.size(); i++) {
		if (layerIndex_.at(i) == -1 && !kStarted) {
			startOffset = i;
		}
		else if (layerIndex_.at(i) != -1) {
			kStarted = kTRUE;
		}
	}

	if (startOffset>0) {
		for (Int_t i=0; i<=startOffset; i++)
			layerIndex_.at(i) = 0;
	}
}

Int_t Clusters::getFirstIndexOfLayer(UInt_t layer) {
   // return index in clusters based on fLayerIndex
   // Return value -1 == no cluster in that layer

   if (layerIndex_.size() == 0) return -1;
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

Int_t Clusters::getLastActiveLayer() {
	Int_t lastActiveLayer = 0;
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i))
			continue;
		if (getLayer(i) > lastActiveLayer)
			lastActiveLayer = getLayer(i);
	}
	return lastActiveLayer;
}

Bool_t Clusters::isPointOutOfBounds(Cluster *point) {
	Bool_t isOutside;
   Float_t x = point->getX();
   Float_t y = point->getY();

   if (!point)
   	isOutside = kTRUE;
   else {
		if (x < 0 || x > 2*nx || y < 0 || y > 2*ny)
			isOutside = kTRUE;
		else
			isOutside = kFALSE;
   }

   return isOutside;
}

void Clusters::removeSmallClusters(Int_t size) {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (getSize(i) <= size) {
			removeClusterAt(i);
		}
	}
}
