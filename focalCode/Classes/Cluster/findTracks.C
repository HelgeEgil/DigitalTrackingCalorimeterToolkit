#include <iostream>
#include <vector>

#include <TClonesArray.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;

Tracks * Clusters::findCalorimeterTracks() {
	Tracks *tracks = new Tracks(kEventsPerRun * 5);
	Int_t startOffset = 0;

	makeLayerIndex();

	fillMCSRadiusList(2);
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
	if (clustersLeft>0) {
		cout << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track! (" << factor << " %)\n";
	}

	return tracks;
}

void Clusters::findTracksFromLayer(Tracks * tracks, Int_t layer, Int_t trackFindingAlgorithm) {
	Int_t startOffset = 0;
	Track *bestTrack = 0;
	
	Clusters * seeds = findSeeds(layer);
	
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

		if (clustersFromThisLayer->GetEntriesFast()) {
			break;
		}
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

	Float_t maxDelta = getSearchRadiusForLayer(nextLayer);

	if (layerIdxFrom < 0)
		return clustersFromThisLayer; // empty

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

		if (diffmmXY(seed, At(i)) < maxDelta) // was initialSearchRadius
			clustersFromThisLayer->appendCluster(At(i));
	}
	return clustersFromThisLayer;
}

void Clusters::doNearestClusterTrackPropagation(Track *track, Int_t lastHitLayer) {
	Int_t nSearchLayers = getLastActiveLayer();
	Cluster * projectedPoint = 0;
	Cluster * nearestNeighbour = 0;
	Cluster * skipNearestNeighbour = 0;
	Float_t delta, skipDistance;

	// Used to be layer < nSearchlayers - 1
	// I don't know why... 
	
	for (Int_t layer = lastHitLayer + 1; layer <= nSearchLayers+1; layer++) {
		searchRadius = getSearchRadiusForLayer(layer);
		projectedPoint = getTrackPropagationToLayer(track, layer);

		if (isPointOutOfBounds(projectedPoint)) {
			break;
		}

		nearestNeighbour = findNearestNeighbour(projectedPoint);
		delta = diffmmXY(projectedPoint, nearestNeighbour);

		if (kDebug && delta > 0) {
			cout << "Searching layer " << layer << ", and found nearest neighbour " << * nearestNeighbour << " with delta " << delta << " mm.\n";
		}
		if (kDebug && delta < 0) {
			cout << "Searching layer " << layer << " and did not find any neighbours!\n";
		}
		
		Bool_t skipCurrentLayer = (delta > searchRadius / 2);

		if (skipCurrentLayer) {
			projectedPoint = getTrackPropagationToLayer(track, layer+1);

			if (!isPointOutOfBounds(projectedPoint)) { // repeat process in layer+1
				skipNearestNeighbour = findNearestNeighbour(projectedPoint);
				skipDistance = diffmmXY(projectedPoint, skipNearestNeighbour);

				if (skipDistance * 1.2 < delta  && skipDistance > 0) {
					track->appendCluster(skipNearestNeighbour);
					lastHitLayer = ++layer+1; // don't search next layer...
				}
				
				else if (delta > 0) {
					track->appendCluster(nearestNeighbour);
					lastHitLayer = layer+1;
				}
			}
		}

		else if (delta > 0) { // and <searchRadius/2
			track->appendCluster(nearestNeighbour);
			lastHitLayer = layer;
		}

		if (layer>lastHitLayer+2) {
			break;
		}
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
	Float_t delta;
	Bool_t kFoundNeighbour = kFALSE;

	Int_t searchLayer = projectedPoint->getLayer();
	Int_t layerIdxFrom = getFirstIndexOfLayer(searchLayer);
	Int_t layerIdxTo = getLastIndexOfLayer(searchLayer);
	
	Float_t maxDelta = getSearchRadiusForLayer(searchLayer); // maximum distance for nearest neighbour

	if (layerIdxFrom < 0)
		return 0;

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

		delta = diffmmXY(projectedPoint, At(i));
		if (delta < maxDelta) {
			nearestNeighbour->set(At(i));
			maxDelta = delta;
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
	Float_t maxDelta = getSearchRadiusForLayer(searchLayer);

	Int_t layerIdxFrom = getFirstIndexOfLayer(searchLayer);
	Int_t layerIdxTo = getLastIndexOfLayer(searchLayer);

	if (layerIdxFrom < 0) {
		return nearestNeighbours;
	}

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

		delta = diffmmXY(projectedPoint, At(i));
		if (delta < maxDelta) {
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
