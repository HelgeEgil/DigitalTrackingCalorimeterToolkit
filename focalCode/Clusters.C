#include "Clusters.h"

#include "Cluster.h"
#include "Track.h"
#include "Constants.h"
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
   clearClusters();
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
}

void Clusters::appendClusterWithoutTrack(Cluster *cluster) {
   Int_t i = clustersWithoutTrack_.GetEntriesFast();
   Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());
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
      if (getFirstIndexOfLayer(i) != -999 ) {
         return getFirstIndexOfLayer(i);
      }
   }
   return GetEntriesFast(); // no hits beyond layer
}

Tracks * Clusters::findTracks() {
   Tracks *tracks = new Tracks(kEventsPerRun * 5);
   Track *bestTrack = new Track();
   Int_t startOffset = 0;

   makeLayerIndex();

   findTracksFromLayer(tracks, 0);
   findTracksFromLayer(tracks, 1);

   cout << "Found " << tracks->GetEntriesFast() << " tracks.\n";

	Int_t clustersLeft = 0;
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (At(i)) {
		   clustersLeft++;
		   appendClusterWithoutTrack(At(i));
      }
	}
	Float_t factor = 100 * (1 - (Float_t) clustersLeft / (Float_t) GetEntriesFast());
   cout << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track! (" << factor << " %)\n";

   return tracks;
}

void Clusters::findTracksFromLayer(Tracks * tracks, Int_t layer) {
	Int_t startOffset = 0;

	Clusters * seeds = findSeeds(layer);
	   for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
			if (!seeds->At(i))
				continue;

			Track * bestTrack = growTracksFromSeed(seeds->At(i));

	      if (bestTrack->GetEntriesFast() > 0) {
	      	if (!bestTrack->At(0))
	      		startOffset = 1;

	         tracks->appendTrack(bestTrack, startOffset);
	         removeAllClustersInTrack(bestTrack);
	      }
	      delete bestTrack;
	   }
	   delete seeds;
}

Clusters * Clusters::findSeeds(Int_t layer) {
	Clusters *seeds = new Clusters(1000);
   Int_t layerIdxFrom = getFirstIndexOfLayer(layer);
   Int_t layerIdxTo = getLastIndexOfLayer(layer);

   if (layerIdxFrom<0)
   	return seeds;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;
      seeds->appendCluster(At(i));
   }
   return seeds;
}

Track * Clusters::growTracksFromSeed(Cluster *seed) {
   Tracks *seedTracks = new Tracks(100);
   Track *currentTrack = new Track();
   Track *longestTrack = new Track();

	Clusters * nextClusters = findNearestClustersInNextLayer(seed);

   for (Int_t i=0; i<nextClusters->GetEntriesFast(); i++) {
      currentTrack->appendCluster(seed);
      currentTrack->appendCluster(nextClusters->At(i));

      doTrackPropagation(currentTrack, nextClusters->getLayer(i));
      if (currentTrack->GetEntriesFast())
         seedTracks->appendTrack(currentTrack);

      currentTrack->clearTrack();
   }

   if (seedTracks->GetEntriesFast())
   	longestTrack = findLongestTrack(seedTracks);

   return longestTrack;
}

Clusters * Clusters::findNearestClustersInNextLayer(Cluster *seed) {
	Clusters *nextClusters = new Clusters(50);
	Clusters *clustersFromThisLayer = 0;

	Int_t layerCounter = 1;

	for (Int_t skipLayers=0; skipLayers<2; skipLayers++) {
		Int_t nextLayer = seed->getLayer() + 1 + skipLayers;
		clustersFromThisLayer = findClustersFromSeedInLayer(seed, nextLayer);

		if (clustersFromThisLayer->GetEntriesFast())
			break;
	}

	for (Int_t i=0; i<clustersFromThisLayer->GetEntriesFast(); i++) {
		nextClusters->appendCluster(clustersFromThisLayer->At(i));
	}

	delete clustersFromThisLayer;
	return nextClusters;
}

Clusters * Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer) {
   Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
   Clusters *clustersFromThisLayer = new Clusters(50);

	if (layerIdxFrom < 0)
		return 0;

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

		if (diffmm(seed, At(i)) < searchRadius)
			clustersFromThisLayer->appendCluster(At(i));
	}
	return clustersFromThisLayer;
}

void Clusters::doTrackPropagation(Track *track, Int_t lastHitLayer) {
	Int_t nSearchLayers = getLastActiveLayer();
	Cluster * projectedPoint = 0;
	Cluster * nearestNeighbour = 0;
	Cluster * skipNearestNeighbour = 0;
	Float_t distance, skipDistance;

   for (Int_t layer = lastHitLayer + 1; layer<nSearchLayers-1; layer++) {
   	projectedPoint = getTrackPropagationToLayer(track, layer);

   	if (isPointOutOfBounds(projectedPoint))
   		break;

      nearestNeighbour = findNearestNeighbour(projectedPoint);
      distance = diffmm(projectedPoint, nearestNeighbour);

      if (distance > secondSearchRadius/2) {
      	projectedPoint = getTrackPropagationToLayer(track, layer+1);
      	if (!isPointOutOfBounds(projectedPoint)) {
      		skipNearestNeighbour = findNearestNeighbour(projectedPoint);
      		skipDistance = diffmm(projectedPoint, skipNearestNeighbour);
      		if (skipDistance * 1.2 < distance  && skipDistance > 0) {
      			cout << "A better match was found in layer+1! (" << skipDistance << " vs " << distance << ")\n";
      			track->appendCluster(skipNearestNeighbour);
      			lastHitLayer = ++layer+1; // don't search next layer...
      		}
      		else if (distance > 0) {
      			cout << "No better matches were found (distance = " << distance << ")\n";
      			track->appendCluster(nearestNeighbour);
      			lastHitLayer = layer+1;
      		}
      	}
      }

      else if (distance > 0) {
         track->appendCluster(nearestNeighbour);
         lastHitLayer = layer;
      }

      if (layer>lastHitLayer+2)
      	break;
   }
   delete projectedPoint;
   delete nearestNeighbour;
   delete skipNearestNeighbour;
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

Cluster * Clusters::getTrackPropagationToLayer(Track *track, Int_t layer) {
   Cluster *nextProjectedPoint = new Cluster();
   Int_t last = track->GetEntriesFast() - 1;
   Int_t diffLayer = layer - track->getLayer(last);

   Cluster p1(track->getX(last-1), track->getY(last-1));
   Cluster p2(track->getX(last), track->getY(last));

   Cluster slope(p2.getX() - p1.getX(), p2.getY() - p1.getY());

   Float_t x = p2.getX() + diffLayer * slope.getX();
   Float_t y = p2.getY() + diffLayer * slope.getY();

   nextProjectedPoint->set(x, y, layer);

   return nextProjectedPoint;
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

Cluster * Clusters::findNearestNeighbour(Cluster *projectedPoint) {
   Cluster *nearestNeighbour = new Cluster();
   Float_t distance = secondSearchRadius; // maximum distance for nearest neighbour
   Float_t delta;
   Bool_t kFoundNeighbour = kFALSE;

   Int_t searchLayer = projectedPoint->getLayer();
   Int_t layerIdxFrom = getFirstIndexOfLayer(searchLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(searchLayer);

   if (layerIdxFrom < 0)
   		return 0;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

      delta = diffmm(projectedPoint, At(i));
      if (delta < distance) {
			nearestNeighbour->set(At(i));
			distance = delta;
			kFoundNeighbour = kTRUE;
      }
   }

   if (!kFoundNeighbour)
   	return 0;

	return nearestNeighbour;
}

Track * Clusters::findLongestTrack(Tracks *seedTracks) {
   Int_t bestIdx = -1;
   Float_t bestLen = -1;

   Track * longestTrack = new Track();

   for (Int_t i=0; i<seedTracks->GetEntriesFast(); i++) {
   	Float_t trackLength = seedTracks->getTrackLengthmm(i);
   	if (trackLength > bestLen) {
         bestLen = trackLength;
         bestIdx = i;
   	}
   }

   Int_t startOffset = 0;
   Int_t trackStartingLayer = seedTracks->At(bestIdx)->getLayer(0);

   if (trackStartingLayer > 0)
   	startOffset = 1;

   longestTrack->setTrack((Track*) seedTracks->At(bestIdx), startOffset);
   return longestTrack;
}

void Clusters::removeSmallClusters(Int_t size) {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (getSize(i) <= size) {
			removeClusterAt(i);
		}
	}
}

Float_t Clusters::diffmm(Cluster *p1, Cluster *p2) {
	if (!p1 || !p2)
		return -1;

	Float_t diffx = p2->getXmm() - p1->getXmm();
	Float_t diffy = p2->getYmm() - p1->getYmm();
	return sqrt(pow(diffx,2) + pow(diffy,2));
}
