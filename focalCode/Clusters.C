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
   Clear();
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
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i);
   c->set(cluster->getX(), cluster->getY(), cluster->getLayer(), cluster->getSize());
}

void Clusters::removeAllClustersInTrack(Track *track) {
   // Remove all points in track from clusters_
   //
   // since both fClusterCollection and track are sorted
   // by increasing layer
   // we keep searching from the index where the last point
   // in track was found
   // Instead of O(100) passes through fClusterCollection,
   // it's only 1

   Int_t lastIndex = 0;

   for (Int_t i = 0; i < track->GetEntriesFast(); i++) {
      if (!track->At(i)) continue;
      Int_t layer = track->getLayer(i);
      Float_t x = track->getX(i);
      Float_t y = track->getY(i);

      for (Int_t j=lastIndex; j<GetEntriesFast(); j++) {

			if (!At(j)) continue; // NaC
         if (getLayer(j) < layer) continue;
         if (getX(j) == x) {
            if (getY(j) == y && getLayer(j) == layer) {
               // it's a match, remove from fClusterCollection
               removeClusterAt(j);
               lastIndex = j+1;
            } // end check for point validity (y,z)
         } // end check for point validity (x)
      } // end loop over this track clusters
   } // end loop over foreign track clusters
} // end function RemoveTrack

Bool_t Clusters::removeClusterFromCoords(Float_t x, Float_t y, Int_t layer) {
   // loop through fClusterCollection to find matching point
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
   // update or make layerIndex_

	Int_t layerOffset = 0;

   // Purge layerIndex_
   if (layerIndex_.size() == 0) {
      for (Int_t i=0; i<nLayers; i++) layerIndex_.push_back(-999);
   }
   else {
      for (Int_t i=0; i<nLayers; i++) layerIndex_.at(i) = -999;
   }

   // Make new list
   if (GetEntriesFast() > 0) {
      Int_t lastLayer = -999;
      for (Int_t i=0; i<GetEntriesFast(); i++) {
			
			if (!At(i)) continue;
         if (lastLayer != (getLayer(i) + layerOffset)) {
				if (getLayer(i)<0 && getLayer(i)>-100) {
					layerOffset = nTrackers;
				}

				cout << "Setting idx " << i << " to layer " << getLayer(i) << " + " << layerOffset << endl;
            layerIndex_.at(getLayer(i) + layerOffset) = i;
            lastLayer = getLayer(i) + layerOffset;
         } // end check for new layer
      } // end loop over all tracks

		Int_t startOffset = 0;
		Bool_t kStarted = kFALSE;
		for (UInt_t i=0; i<layerIndex_.size(); i++) {
			if (layerIndex_.at(i) == -999 && !kStarted) {
				startOffset = i;
			}
			else if (layerIndex_.at(i) != -999) {
				kStarted = kTRUE;
			}
		}

		if (startOffset>0) {
			for (Int_t i=0; i<=startOffset; i++)
				layerIndex_.at(i) = 0;
		}


   }
	for (UInt_t i=0; i<layerIndex_.size(); i++) cout << i << ": " << layerIndex_.at(i) << ". ";
	cout << endl;
}

Int_t Clusters::getFirstIndexOfLayer(Int_t layer) {
   // return index in clusters based on fLayerIndex
   // Return value -1 == no cluster in that layer

   if (layerIndex_.size() == 0) return -1;
	if (layerIndex_.size() < layer) return -1;

   return layerIndex_.at(layer);
}

Int_t Clusters::getLastIndexOfLayer(Int_t layer) {
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
   cout << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track!\n";

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

	Clusters * nextClusters = findNearestClustersInNextLayer(seed);

   for (Int_t i=0; i<nextClusters->GetEntriesFast(); i++) {
      currentTrack->appendCluster(seed);
      currentTrack->appendCluster(nextClusters->At(i));

      doTrackPropagation(currentTrack, nextClusters->getLayer(i));
      if (currentTrack->GetEntriesFast())
         seedTracks->appendTrack(currentTrack);

      currentTrack->Clear();
   }

   Track * longestTrack = findLongestTrack(seedTracks);
   return longestTrack;
}

Clusters * Clusters::findNearestClustersInNextLayer(Cluster *seed) {
	Clusters *nextClusters = new Clusters(50);

	Int_t layerOffset = 0;
	if (frameType_ == kTracker)
	   layerOffset = nTrackers;

	Int_t layerCounter = 1;
	while (!nextClusters->GetEntriesFast()) {
		Int_t nextLayer = seed->getLayer() + layerOffset + layerCounter++;
		Clusters *clustersFromThisLayer = findClustersFromSeedInLayer(seed, nextLayer);

		for (Int_t i=0; i<clustersFromThisLayer->GetEntriesFast(); i++)
			nextClusters->appendCluster(clustersFromThisLayer->At(i));
	}
	return nextClusters;
}

Clusters * Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer) {
   Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
   Clusters *clustersFromThisLayer = new Clusters(50);

	if (layerIdxFrom<0)
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
   for (Int_t layer = lastHitLayer + 1; layer<nSearchLayers-1; layer++) {
   	Cluster * projectedPoint = getTrackPropagationToLayer(track, layer);

   	if (isPointOutOfBounds(projectedPoint))
   		break;

      Cluster * nearestNeighbour = findNearestNeighbour(projectedPoint);
      Float_t distance = diffmm(projectedPoint, nearestNeighbour);

      if (distance > 0) {
         track->appendCluster(nearestNeighbour);
         lastHitLayer = layer;
      }

      if (layer>lastHitLayer+2)
      	break;
   }
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
	Bool_t isInside;
   Float_t x = point->getX();
   Float_t y = point->getY();

   if (x < 0 || x > 2*nx || y < 0 || y > 2*ny)
   	isInside = kFALSE;
   else
   	isInside = kTRUE;

   return isInside;
}

Cluster * Clusters::findNearestNeighbour(Cluster *projectedPoint) {
   Cluster *nearestNeighbour = new Cluster();
   Float_t distance = secondSearchRadius; // maximum distance for nearest neighbour
   Float_t delta;

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
      }
   }
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

   // if first layer is >0 add empty track to make space for extrapolation
   Int_t startOffset = 0;
   if (seedTracks->At(bestIdx)->getLayer(0)>0)
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
