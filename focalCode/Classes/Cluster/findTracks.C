#include <vector>

#include <TClonesArray.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

#ifdef USEDEBUG
#define showDebug(x) std::cout << x
#else
#define showDebug(x)
#endif

using namespace std;

Tracks * Clusters::findCalorimeterTracks() {
	Tracks *tracks = new Tracks(kEventsPerRun * 5);
	Int_t startOffset = 0;
	Bool_t usedClustersInSeeds = true;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		appendClusterWithoutTrack(At(i));
	}
	
	makeLayerIndex();
	fillMCSRadiusList();
	MCSMultiplicationFactor = 3;

	// first pass, small search cone (3 sigma MCS)
	findTracksFromLayer(tracks, 0, usedClustersInSeeds);

	MCSMultiplicationFactor = 5;
	usedClustersInSeeds = false;
	multiplyRadiusFirstLayers(2);
	
	if (clustersWithoutTrack_.GetEntries() > 0) {
		// second pass, large seach cone (5 sigma MCS)
		findTracksFromLayer(tracks, 0, usedClustersInSeeds);
	
		// third pass, start from layer 1
		findTracksFromLayer(tracks, 1, usedClustersInSeeds);
	}

	Int_t clustersLeft = clustersWithoutTrack_.GetEntries();
	Float_t factor = 100 * (1 - (Float_t) clustersLeft / GetEntriesFast());
	cout << "Found " << tracks->GetEntriesFast() << " tracks. " << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track! (" << factor << " %)\n";

	return tracks;
}

void Clusters::findTracksFromLayer(Tracks * tracks, Int_t layer, Bool_t kUsedClustersInSeeds) {
	Int_t			startOffset = 0;
	Track		 *	bestTrack = nullptr;
	Track		 *	newBestTrack = nullptr;
	Clusters	 *	seeds = nullptr;
	Int_t			minimumLengthOfTrackToPass = 4; // first pass

	if (!kUsedClustersInSeeds) {
		minimumLengthOfTrackToPass = 1; // second / third pass
	}

	seeds = findSeeds(layer, kUsedClustersInSeeds);
	
	for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
		if (!seeds->At(i))
			continue;

		bestTrack = trackPropagation(seeds->At(i));

		if (bestTrack->GetEntriesFast() < minimumLengthOfTrackToPass)
			continue;

		if (!bestTrack->At(0))
			startOffset = 1;

		tracks->appendTrack(bestTrack, startOffset);
		removeTrackFromClustersWithoutTrack(bestTrack);
		markUsedClusters(bestTrack);
	}

	Int_t nTracksWithConflicts = 0;
	for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
		if (!tracks->At(i)) continue;
		nTracksWithConflicts += (int) tracks->isUsedClustersInTrack(i);
	}

	delete bestTrack;
	delete seeds;
}

Clusters * Clusters::findSeeds(Int_t layer, Bool_t kUsedClustersInSeeds) {
	Clusters *seeds = new Clusters(1000);
	
	Int_t layerIdxFrom = getFirstIndexOfLayer(layer);
	Int_t layerIdxTo = getLastIndexOfLayer(layer);

	showDebug("Layer Index From , To: (" << layerIdxFrom << "," << layerIdxTo << ")\n");
	
	if (layerIdxFrom<0)
		return seeds;

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i)) { continue; }
		if (!kUsedClustersInSeeds && isUsed(i)) { continue; }
		
		seeds->appendCluster(At(i));
	}

	return seeds;
}

Track * Clusters::trackPropagation(Cluster *seed) {
	Tracks *	seedTracks = new Tracks(100);
	Track  *	currentTrack = new Track();
	Track  *	longestTrack = new Track();

	Clusters * nextClusters = findNearestClustersInNextLayer(seed);

	for (Int_t i=0; i<nextClusters->GetEntriesFast(); i++) {
		Cluster * nextCluster = nextClusters->At(i);
		Int_t fromLayer = nextCluster->getLayer();

		currentTrack->appendCluster(seed);
		currentTrack->appendCluster(nextCluster);

		growTrackFromLayer(currentTrack, fromLayer);
		
		if (currentTrack->GetEntriesFast())
			seedTracks->appendTrack(currentTrack);

		currentTrack->clearTrack();
	}
	
	if (seedTracks->GetEntriesFast())
		longestTrack = findLongestTrack(seedTracks);

	delete seedTracks;
	delete currentTrack;
	delete nextClusters;

	return longestTrack;
}

Clusters * Clusters::findNearestClustersInNextLayer(Cluster *seed) {
	Clusters *nextClusters = new Clusters(50);
	Clusters *clustersFromThisLayer = 0;

	Int_t layerCounter = 1;
	for (Int_t skipLayers=0; skipLayers<2; skipLayers++) {
		Int_t nextLayer = seed->getLayer() + 1 + skipLayers;
		clustersFromThisLayer = findClustersFromSeedInLayer(seed, nextLayer);

		if (clustersFromThisLayer->GetEntriesFast()) {
			break;
		}
	}

	Int_t nClusters = clustersFromThisLayer->GetEntriesFast();
	for (Int_t i=0; i<nClusters; i++) {
		if (!clustersFromThisLayer->At(i)) { continue; }
		
		nextClusters->appendCluster(clustersFromThisLayer->At(i));
	}

	delete clustersFromThisLayer;

	return nextClusters;
}

Clusters * Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer) {
	Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
	Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
	Clusters *clustersFromThisLayer = new Clusters(50);

	Float_t maxDelta = getSearchRadiusForLayer(nextLayer) * 0.75 * MCSMultiplicationFactor;

	if (layerIdxFrom < 0)
		return clustersFromThisLayer; // empty

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i)) { continue; }

		if (diffmmXY(seed, At(i)) < maxDelta) {
			clustersFromThisLayer->appendCluster(At(i));
		}
	}
	return clustersFromThisLayer;
}

void Clusters::growTrackFromLayer(Track *track, Int_t fromLayer) {
	Cluster   *	projectedPoint = 0;
	Cluster   *	nearestNeighbour = 0;
	Cluster   *	skipNearestNeighbour = 0;
	Float_t		delta, skipDistance;
	Int_t			nSearchLayers = getLastActiveLayer();
	Int_t			lastHitLayer = fromLayer;
	Bool_t		skipCurrentLayer = false;

	for (Int_t layer = lastHitLayer + 1; layer <= nSearchLayers+1; layer++) {
		searchRadius = getSearchRadiusForLayer(layer) * MCSMultiplicationFactor;

		projectedPoint = getTrackExtrapolationToLayer(track, layer);
		if (isPointOutOfBounds(projectedPoint, searchRadius / dx))
			break;

		nearestNeighbour = findNearestNeighbour(projectedPoint);

		delta = diffmmXY(projectedPoint, nearestNeighbour);
		skipCurrentLayer = (delta > searchRadius);

		if (skipCurrentLayer) {
			projectedPoint = getTrackExtrapolationToLayer(track, layer+1);

			if (!isPointOutOfBounds(projectedPoint, searchRadius / dx)) { // repeat process in layer+1
				skipNearestNeighbour = findNearestNeighbour(projectedPoint);
				skipDistance = diffmmXY(projectedPoint, skipNearestNeighbour);

				if (skipDistance * 1.2 < delta  && skipDistance > 0) {
					track->appendCluster(skipNearestNeighbour);
					lastHitLayer = ++layer+1; // don't search next layer...
				}
				
				else if (delta > 0) { // append cluster in layer+1 to track
					// This is a test.. Remove if it doesn't give warnings
					for (Int_t j=0; j<track->GetEntriesFast(); j++) {
						if (!track->At(j)) continue;
						if (track->getLayer(j) == nearestNeighbour->getLayer()) {
							cout << "SOMETHING IS WRONG IN Clusters::growTrackFromLayer()! Track so far = " << *track << " and it wants to add nearestNeighbour " << *nearestNeighbour << ".\n";
						}
					}

					track->appendCluster(nearestNeighbour);
					lastHitLayer = layer+1;
				}
			}
		}

		else if (delta > 0) { // append cluster to track
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

Cluster * Clusters::findNearestNeighbour(Cluster *projectedPoint, Bool_t rejectUsed) {
	Cluster *nearestNeighbour = new Cluster();
	Float_t	delta;
	Bool_t	kFoundNeighbour = false;
	Bool_t	reject = false;
	Int_t		searchLayer = projectedPoint->getLayer();
	Int_t		layerIdxFrom = getFirstIndexOfLayer(searchLayer);
	Int_t		layerIdxTo = getLastIndexOfLayer(searchLayer);
	Float_t	maxDelta = getSearchRadiusForLayer(searchLayer) * MCSMultiplicationFactor;

	if (layerIdxFrom < 0)
		return 0;

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

		delta = diffmmXY(projectedPoint, At(i));
		reject = (At(i)->isUsed() && rejectUsed);

		if (delta < maxDelta && !reject) {
			nearestNeighbour->set(At(i));
			maxDelta = delta;
			kFoundNeighbour = kTRUE;
		}
	}

	if (!kFoundNeighbour)
		return 0;

	return nearestNeighbour;
}

Track * Clusters::findLongestTrack(Tracks *seedTracks) {
	if (!seedTracks->GetEntriesFast())
		return new Track();
	
	Float_t	bestScore = -1;
	Float_t	score = -1;
	Track	 *	longestTrack = new Track();
	Track	 *	track = 0;
	Int_t		startOffset;

	for (Int_t i=0; i<seedTracks->GetEntriesFast(); i++) {
		if (!seedTracks->At(i)) continue;

		score = seedTracks->getTrackScore(i);
		if (score > bestScore) {
			bestScore = score;
			track = seedTracks->At(i);
		}
	}

	startOffset = (track->getLayer(0)) ? 0 : 1;

	longestTrack->setTrack(track, startOffset);

	return longestTrack;
}
