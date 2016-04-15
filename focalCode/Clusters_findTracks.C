#include "Clusters.h"
#include "Cluster.h"
#include "Track.h"
#include "Tracks.h"
#include "Constants.h"
#include "MaterialConstants.h"
#include "Tools.h"
#include <iostream>
#include <TClonesArray.h>
#include <vector>

using namespace std;

Tracks * Clusters::findTrackerTracks() {
	Tracks * tracks = new Tracks(kEventsPerRun * 5);

	makeLayerIndex();

	findTracksFromLayer(tracks, 0);

	Int_t nComplete = 0;
	for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
		Int_t nClusters = 0;
		if (!At(i)) continue;

		for (Int_t j=0; j<tracks->GetEntriesFast(i); j++) {
			if (tracks->At(i)->At(j))
				nClusters++;
		}

		if (nClusters == 4) nComplete++;
	}

	Int_t clustersLeft = 0;
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (At(i)) {
			clustersLeft++;
			appendClusterWithoutTrack(At(i));
		}
	}

	Float_t factor = 100 * (1 - (Float_t) clustersLeft / GetEntriesFast());
	cout << clustersLeft << " of total " << GetEntriesFast() << " tracker clusters were not assigned to track (" << factor << "%)\n";

	return tracks;
}

Tracks * Clusters::findCalorimeterTracks() {
	Tracks *tracks = new Tracks(kEventsPerRun * 5);
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
	Float_t factor = 100 * (1 - (Float_t) clustersLeft / GetEntriesFast());
	cout << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track! (" << factor << " %)\n";

	return tracks;
}

void Clusters::findTracksFromLayer(Tracks * tracks, Int_t layer, Int_t trackFindingAlgorithm) {
	Int_t startOffset = 0;
	Track *bestTrack = 0;
	
	Clusters * seeds = findSeeds(layer);
	
	if (kDebug) {
		cout << "Number of seeds in layer " << layer << ": " << seeds->GetEntriesFast() << endl;
	}
	
	for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
		if (!seeds->At(i))
			continue;

		if (trackFindingAlgorithm < 0)
			trackFindingAlgorithm = kTrackFindingAlgorithm; // in Constants.h

		if (trackFindingAlgorithm == kRecursive)
			bestTrack = recursiveTrackPropagation(seeds->At(i), Track(seeds->At(i)));
		else if (trackFindingAlgorithm == kNearestCluster)
			bestTrack = nearestClusterTrackPropagation(seeds->At(i));

		if (bestTrack->GetEntriesFast() > 0) {
			if (!bestTrack->At(0))
				startOffset = 1;

			tracks->appendTrack(bestTrack, startOffset);
			removeAllClustersInTrack(bestTrack);
		}
	}

	delete bestTrack;
	delete seeds;
}

Clusters * Clusters::findSeeds(Int_t layer) {
	Clusters *seeds = new Clusters(1000);
	
	Int_t layerIdxFrom = getFirstIndexOfLayer(layer);
	Int_t layerIdxTo = getLastIndexOfLayer(layer);

	if (kDebug) {
		cout << "Layer Index From , To: (" << layerIdxFrom << "," << layerIdxTo << ")\n";
	}
	
	if (layerIdxFrom<0)
		return seeds;

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;
		seeds->appendCluster(At(i));
	}
	return seeds;
}

Track * Clusters::recursiveTrackPropagation(Cluster *cluster, Track currentTrack) {
	currentTrack.appendCluster(cluster);

	Tracks *tracksFromThisCluster = new Tracks(100);
	Track *bestTrack = new Track();
	Track *bestNextTrack = 0;
	Cluster *projectedPoint = 0;
	Cluster *nextClosePoint = 0;

	Clusters *closePointsInNextLayer = new Clusters();

	Int_t layer = cluster->getLayer() + 1;

	if (layer >= nLayers) {
		return new Track();
	}

	while (!closePointsInNextLayer->GetEntriesFast()) {
		projectedPoint = getTrackPropagationToLayer(&currentTrack, layer++);
		if (isPointOutOfBounds(projectedPoint))
			return new Track();

		if (layer > cluster->getLayer() + 3)
			break;

		closePointsInNextLayer = findAllClosePointsInNextLayer(projectedPoint);
	}

	for (Int_t i=0; i<closePointsInNextLayer->GetEntriesFast(); i++) {
		nextClosePoint = closePointsInNextLayer->At(i);
		bestNextTrack = recursiveTrackPropagation(nextClosePoint, currentTrack);

		if (bestNextTrack->GetEntriesFast())
			tracksFromThisCluster->appendTrack(bestNextTrack);
	}

	bestTrack->appendCluster(cluster);

	Track *track = findLongestTrack(tracksFromThisCluster);
	if (track->GetEntriesFast()) {
		for (Int_t i=0; i<track->GetEntriesFast(); i++) {
			bestTrack->appendCluster(track->At(i));
		}
	}

	delete tracksFromThisCluster;
	delete projectedPoint;
	delete closePointsInNextLayer;
	delete nextClosePoint;
	delete bestNextTrack;
	delete track;

	return bestTrack;
}

Track * Clusters::nearestClusterTrackPropagation(Cluster *seed) {
	Tracks *seedTracks = new Tracks(100);
	Track *currentTrack = new Track();
	Track *longestTrack = new Track();

	Clusters * nextClusters = findNearestClustersInNextLayer(seed);

	for (Int_t i=0; i<nextClusters->GetEntriesFast(); i++) {
		currentTrack->appendCluster(seed);
		currentTrack->appendCluster(nextClusters->At(i));

		doNearestClusterTrackPropagation(currentTrack, nextClusters->getLayer(i));
		if (currentTrack->GetEntriesFast())
			seedTracks->appendTrack(currentTrack);

		currentTrack->clearTrack();
	}
	
	if (seedTracks->GetEntriesFast()) {
		longestTrack = findLongestTrack(seedTracks);
		if (kDebug) cout << Form("Longest track is %d layers deep with an energy of about %.1f MeV.\n", longestTrack->GetEntriesFast(), longestTrack->getEnergy());
	}

	delete seedTracks;
	delete currentTrack;
	delete nextClusters;

	return longestTrack;
}

Clusters * Clusters::findNearestClustersInNextLayer(Cluster *seed) {
	Clusters *nextClusters = new Clusters(50);
	Clusters *clustersFromThisLayer = 0;

	Int_t layerCounter = 1;

	for (Int_t skipLayers=0; skipLayers<3; skipLayers++) {
		Int_t nextLayer = seed->getLayer() + 1 + skipLayers;
		clustersFromThisLayer = findClustersFromSeedInLayer(seed, nextLayer);

		if (clustersFromThisLayer->GetEntriesFast())
			break;
	}

	for (Int_t i=0; i<clustersFromThisLayer->GetEntriesFast(); i++) {
		if (!clustersFromThisLayer->At(i)) continue;
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

		if (diffmm(seed, At(i)) < initialSearchRadius)
			clustersFromThisLayer->appendCluster(At(i));
	}
	return clustersFromThisLayer;
}

void Clusters::doNearestClusterTrackPropagation(Track *track, Int_t lastHitLayer) {
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

		if (distance > searchRadius/2) {
			projectedPoint = getTrackPropagationToLayer(track, layer+1);

			if (!isPointOutOfBounds(projectedPoint)) {
				skipNearestNeighbour = findNearestNeighbour(projectedPoint);
				skipDistance = diffmm(projectedPoint, skipNearestNeighbour);

				if (skipDistance * 1.2 < distance  && skipDistance > 0) {
					track->appendCluster(skipNearestNeighbour);
					lastHitLayer = ++layer+1; // don't search next layer...
				}
				else if (distance > 0) {
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

Cluster * Clusters::getTrackPropagationToLayer(Track *track, Int_t layer) {
	Int_t last = track->GetEntriesFast() - 1;
	Int_t diffLayer = layer - track->getLayer(last);

	Cluster p1(track->getX(last-1), track->getY(last-1));
	Cluster p2(track->getX(last), track->getY(last));

	Cluster slope(p2.getX() - p1.getX(), p2.getY() - p1.getY());

	Float_t x = p2.getX() + diffLayer * slope.getX();
	Float_t y = p2.getY() + diffLayer * slope.getY();

	return new Cluster(x, y, layer);
}

Cluster * Clusters::findNearestNeighbour(Cluster *projectedPoint) {
	Cluster *nearestNeighbour = new Cluster();
	Float_t distance = searchRadius; // maximum distance for nearest neighbour
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

Clusters * Clusters::findAllClosePointsInNextLayer(Cluster *projectedPoint) {
	Clusters *nearestNeighbours = new Clusters();
	Float_t delta;
	Int_t searchLayer = projectedPoint->getLayer();

	Float_t distance;
	if (searchLayer == 1)
		distance = initialSearchRadius;
	else
		distance = searchRadius;

	Int_t layerIdxFrom = getFirstIndexOfLayer(searchLayer);
	Int_t layerIdxTo = getLastIndexOfLayer(searchLayer);

	if (layerIdxFrom < 0) {
		return nearestNeighbours;
	}

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

		delta = diffmm(projectedPoint, At(i));
		if (delta < distance) {
			nearestNeighbours->appendCluster(At(i));
		}
	}
	return nearestNeighbours;
}

Track * Clusters::findLongestTrack(Tracks *seedTracks) {
	Float_t bestScore = -1;

	if (!seedTracks->GetEntriesFast())
		return new Track();

	Track * longestTrack = new Track();
	Track * track = 0;

	for (Int_t i=0; i<seedTracks->GetEntriesFast(); i++) {
		Float_t score = seedTracks->getTrackScore(i);
		if (score > bestScore) {
			bestScore = score;
			track = seedTracks->At(i);
		}
	}

	Int_t startOffset = 0;
	Int_t trackStartingLayer = track->getLayer(0);

	if (trackStartingLayer > 0)
		startOffset = 1;

	longestTrack->setTrack(track, startOffset);

	return longestTrack;
}
