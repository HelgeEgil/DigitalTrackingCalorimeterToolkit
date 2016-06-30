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

Clusters::Clusters(Bool_t frameType) : clusters_("Cluster", kEventsPerRun*2),
											 clustersWithoutTrack_("Cluster", kEventsPerRun*2),
											 conflictClusters_("Cluster", kEventsPerRun*2) {
	frameType_ = frameType;
}

Clusters::~Clusters() {
   // Destructor
   clusters_.Delete();
   clustersWithoutTrack_.Delete();
	conflictClusters_.Delete();
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

	if (cluster->isUsed()) { c->markUsed(); }
	c->setEventID(cluster->getEventID());

}

void Clusters::appendClusterWithoutTrack(Cluster *cluster) {
   Int_t i = clustersWithoutTrack_.GetEntriesFast();
   Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());
	c->setEventID(cluster->getEventID());
}

void Clusters::appendConflictClusters(Clusters *clusters) {
	Int_t i=conflictClusters_.GetEntriesFast();
	Int_t n=clusters->GetEntriesFast();

	for (Int_t j=0; j<n; j++) {
		Cluster *c = (Cluster*) conflictClusters_.ConstructedAt(i+j);
		c->set(clusters->getX(j), clusters->getY(j), clusters->getLayer(j), clusters->getSize(j));
		c->setEventID(clusters->getEventID(j));
	}
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

void Clusters::removeAllClustersInTrackFromClustersWithoutTrack(Track *track) {
   Int_t lastIndex = 0;

   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      Int_t layer = track->getLayer(i);
      Float_t x = track->getX(i);
      Float_t y = track->getY(i);

      for (Int_t j=lastIndex; j<GetEntriesFast(); j++) {

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
void Clusters::markUsedClusters(Track *track) {
   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      Int_t layer = track->getLayer(i);
      Float_t x = track->getX(i);
      Float_t y = track->getY(i);

		Int_t idx = getClusterIdx(x,y,layer);
		if (idx > -1) {
			markUsed(idx);
		}
	}
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

Bool_t Clusters::removeClusterFromCoords(Float_t x, Float_t y, Int_t layer) {
	Int_t clusterIndex = getClusterIdx(x, y, layer);
	if (clusterIndex>-1) {
		removeClusterAt(clusterIndex);
		return true;
	}
	else return false;
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

void Clusters::removeSmallClusters(Int_t size) {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (getSize(i) <= size) {
			removeClusterAt(i);
		}
	}
}


void Clusters::matchWithEventIDs(Hits * eventIDs) {
	
	TStopwatch t1;
	Float_t minDist = 1e5; // px
	Float_t thisDist = 0;
	Cluster *thisCluster = nullptr;
	Hit *thisHit = nullptr;
	Int_t layer = -1;
	Int_t minIdx = 0;
	Bool_t doLoop = true;
	Int_t nHits = eventIDs->GetEntriesFast();

	Float_t cX, cY;

	t1.Start();

	Int_t nClusters = 0;
	nClusters = GetEntries();

	cout << "Number of clusters is " << nClusters << ", number of Hits is " << eventIDs->GetEntriesFast() << endl;

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

//		else { cout << "Could not find hit match for cluster @ " << *thisCluster << endl; }
	}
	t1.Stop();

	cout << "Total time for function: " << t1.RealTime() << " seconds.\n";
	cout << "Number of Hits not matched (out of " << nHits << "): " << eventIDs->GetEntries() << " ( " << eventIDs->GetEntries() / (float) nHits * 100 << "% )" << endl;

	Int_t cWithoutEventID = 0;
	for (Int_t c=0; c<GetEntriesFast(); c++) {
		if (!At(c)) continue;
		if (getEventID(c) < 0) {
			cWithoutEventID++;
			minDist = 1e5;
			for (Int_t h=0; h<eventIDs->GetEntriesFast(); h++) {
				thisHit = eventIDs->At(h);
				if (!thisHit) continue;
				thisDist = diffXY(At(c), thisHit);
				if (thisDist < minDist) minDist = thisDist;
			}
		}
	}

	cout << "Number of clusters without eventID: " << cWithoutEventID << " (" << (float) cWithoutEventID / nClusters * 100 << "%)" << endl;
}

Int_t Clusters::GetEntriesFastLastLayer() {
	Int_t nInLayer[nLayers];
	for (Int_t i=0; i<nLayers; i++) nInLayer[i] = 0;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		nInLayer[getLayer(i)]++;
	}

	Int_t lastLayer = 0;
	for (Int_t i=0; i<nLayers; i++) {
		if (nInLayer[i] > 0) lastLayer = i;
	}
	
	if (nInLayer[lastLayer] < 0.5 * nInLayer[lastLayer - 1]) {
		lastLayer--;
	}

	return nInLayer[lastLayer];
}



