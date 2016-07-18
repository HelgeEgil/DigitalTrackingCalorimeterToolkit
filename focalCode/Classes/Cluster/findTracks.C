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

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		appendClusterWithoutTrack(At(i));
	}
	
	makeLayerIndex();

	Float_t MCSSigma = 3;
	Bool_t usedClustersInSeeds = true;

	fillMCSRadiusList(MCSSigma);
	findTracksFromLayer(tracks, 0, usedClustersInSeeds);

	showDebug( "After first pass, found " << tracks->GetEntriesFast() << " tracks. Size of CWT = " << clustersWithoutTrack_.GetEntries() << endl);

	if (clustersWithoutTrack_.GetEntries() > 0) {

		MCSSigma = 5;
		usedClustersInSeeds = false;
	
		fillMCSRadiusList(MCSSigma);
		multiplyRadiusFirstLayers(2);
		cout << getSearchRadiusForLayer(1) << endl;
		findTracksFromLayer(tracks, 0, usedClustersInSeeds);
	
		showDebug("After second pass, found " << tracks->GetEntriesFast() << " tracks. Size of CWT = " << clustersWithoutTrack_.GetEntries() << endl);
	
		// third pass
		findTracksFromLayer(tracks, 1, usedClustersInSeeds);
	}

	Int_t clustersLeft = clustersWithoutTrack_.GetEntries();
	Float_t factor = 100 * (1 - (Float_t) clustersLeft / GetEntriesFast());
	cout << "Found " << tracks->GetEntriesFast() << " tracks. " << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track! (" << factor << " %)\n";

	return tracks;
}

void Clusters::findTracksFromLayer(Tracks * tracks, Int_t layer, Bool_t kUsedClustersInSeeds) {
	Int_t startOffset = 0;
	Track *bestTrack = nullptr;
	Track *newBestTrack = nullptr;

	Bool_t kIsSecondPass = false;
	if (!kUsedClustersInSeeds) kIsSecondPass = true;

	Int_t minimumLengthOfTrackToPass = 4;

	if (kIsSecondPass) {
		minimumLengthOfTrackToPass = 1;
	}

	Clusters * seeds = findSeeds(layer, kUsedClustersInSeeds);
	Clusters * conflictClusters = nullptr;
	
	for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
		if (!seeds->At(i))
			continue;

		bestTrack = nearestClusterTrackPropagation(seeds->At(i));

		if (bestTrack->GetEntriesFast() >= minimumLengthOfTrackToPass) {
			if (!bestTrack->At(0))
				startOffset = 1;

			tracks->appendTrack(bestTrack, startOffset);
//			removeAllClustersInTrack(bestTrack);
			removeAllClustersInTrackFromClustersWithoutTrack(bestTrack);

			// Get already conflicting clusters with the isUsed tag in bestTrack
			conflictClusters = bestTrack->getConflictClusters();
			appendConflictClusters(conflictClusters); // this TCA may be unneccessary

			// Apply the isUsed tag to the track clusters in clusters_ list
			markUsedClusters(bestTrack);
		}
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

Track * Clusters::nearestClusterTrackPropagation(Cluster *seed) {
	Tracks *seedTracks = new Tracks(100);
	Track *currentTrack = new Track();
	Track *longestTrack = new Track();

	Clusters * nextClusters = findNearestClustersInNextLayer(seed);

	for (Int_t i=0; i<nextClusters->GetEntriesFast(); i++) {
		currentTrack->appendCluster(seed);
		currentTrack->appendCluster(nextClusters->At(i));

		doNearestClusterTrackPropagation(currentTrack, nextClusters->getLayer(i));
		if (currentTrack->GetEntriesFast()) {
			for (Int_t j=0; j<currentTrack->GetEntriesFast(); j++) {
				if (!currentTrack->At(j)) continue;
			}

			seedTracks->appendTrack(currentTrack);
		}

		currentTrack->clearTrack();
	}
	
	if (seedTracks->GetEntriesFast()) {
		longestTrack = findLongestTrack(seedTracks);
		showDebug(Form("Longest track is %d layers deep with an energy of about %.1f MeV.\n", longestTrack->GetEntriesFast(), longestTrack->getEnergy()));
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

	Float_t maxDelta = getSearchRadiusForLayer(nextLayer) * 0.75;

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

void Clusters::doNearestClusterTrackPropagation(Track *track, Int_t lastHitLayer) {
	Int_t nSearchLayers = getLastActiveLayer();
	Cluster * projectedPoint = 0;
	Cluster * nearestNeighbour = 0;
	Cluster * skipNearestNeighbour = 0;
	Float_t delta, skipDistance;

	for (Int_t layer = lastHitLayer + 1; layer <= nSearchLayers+1; layer++) {
		searchRadius = getSearchRadiusForLayer(layer);

		/*		
		if (layer == 2) {
			searchRadius *= 0.5;
		}
		*/

		projectedPoint = getTrackPropagationToLayer(track, layer);

		if (isPointOutOfBounds(projectedPoint, searchRadius / dx)) {
			break;
		}

		nearestNeighbour = findNearestNeighbour(projectedPoint);

		delta = diffmmXY(projectedPoint, nearestNeighbour);

		Bool_t skipCurrentLayer = (delta > searchRadius); // was /2

		if (skipCurrentLayer) {
			projectedPoint = getTrackPropagationToLayer(track, layer+1);

			if (!isPointOutOfBounds(projectedPoint, searchRadius / dx)) { // repeat process in layer+1
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

		else if (delta > 0) {
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
	Float_t delta;
	Bool_t kFoundNeighbour = false;
	Bool_t reject = false;

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
