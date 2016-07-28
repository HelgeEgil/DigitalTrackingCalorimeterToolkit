#include <iostream>
#include <vector>

#include <TClonesArray.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TView.h>
#include <TAttMarker.h>
#include <TAttLine.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>


#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;

Tracks::~Tracks() {
   tracks_.Delete();
   clustersWithoutTrack_.Delete();
}

void Tracks::CompressClusters() {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		At(i)->Compress();
	}
}

void Tracks::Clear(Option_t *) {
	tracks_.Clear("C");
	clustersWithoutTrack_.Clear("C");
}

void Tracks::appendTrack(Track *copyTrack, Int_t startOffset /* default 0 */) {
   Int_t		newIdx = tracks_.GetEntriesFast();
   Track  *	track = (Track*) tracks_.ConstructedAt(newIdx);

   for (Int_t i=0; i<copyTrack->GetEntriesFast(); i++) {
   	if(!copyTrack->At(i))
   		continue;

      track->appendCluster(copyTrack->At(i), startOffset);
   }
}

void Tracks::appendClustersWithoutTrack(TClonesArray *copyCWT) {
	Int_t			idxFrom = clustersWithoutTrack_.GetEntriesFast();
	Cluster   *	newCluster = nullptr;
	Cluster	 *	copyCluster = nullptr;

	for (Int_t i=0; i<copyCWT->GetEntriesFast(); i++) {
		copyCluster = (Cluster*) copyCWT->At(i);
		if (!copyCluster) continue;

		newCluster = (Cluster*) clustersWithoutTrack_.ConstructedAt(idxFrom + i);
		newCluster->set(copyCluster);
	}
}

vector<Int_t> * Tracks::getTracksWithConflictClusters() {
	const Int_t			n = GetEntriesFast();
	Track		       *	thisTrack = nullptr;
	vector<Int_t>   *	trackIdx = new vector<Int_t>;

	trackIdx->reserve(n/10.);

	for (Int_t i=0; i<n; i++) {
		thisTrack = At(i);
		if (!thisTrack) continue;
		if (thisTrack->isUsedClustersInTrack()) { 
			trackIdx->push_back(i);
		}
	}

	return trackIdx;
}

vector<Int_t> * Tracks::getConflictingTracksFromTrack(Int_t trackIdx) {
	Int_t					nClusters, nClustersReal;
	Cluster         *	conflictCluster = nullptr;
	Clusters			 *	conflictClusters = nullptr;
	vector<Int_t>   *	conflictingTracks = new vector<Int_t>;
	vector<Int_t>	 *	possibleTracks = nullptr;
	conflictingTracks->reserve(5);
	
	Track * trackA = At(trackIdx);
	if (!trackA) { return conflictingTracks; }

	conflictClusters = trackA->getConflictClusters();
	nClusters = conflictClusters->GetEntriesFast();
	nClustersReal = conflictClusters->GetEntries();

	for (Int_t i=0; i<nClusters; i++) {
		conflictCluster = conflictClusters->At(i);
		if (!conflictCluster) continue;

		possibleTracks = getTracksFromCluster(conflictCluster);

		for (UInt_t j=0; j<possibleTracks->size(); j++) {
			if (!isItemInVector(possibleTracks->at(j), conflictingTracks)) {
				conflictingTracks->push_back(possibleTracks->at(j));
			}
		}
	}

	return conflictingTracks;
}


void Tracks::extrapolateToLayer0() {
   for (Int_t i=0; i<tracks_.GetEntriesFast(); i++) {
   	if (!At(i)) continue;
      if (!At(i)->At(0)) { // only do if track At(i) misses the first cluster
      	At(i)->extrapolateToLayer0();
		}
   }
}

void Tracks::doFit() {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		At(i)->doFit();
	}
}

Int_t Tracks::getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster) {
	Float_t		minDist = 1e5, dist;
	Int_t			minIdx = -1;
	Cluster	 *	thisCluster = nullptr;

	for (UInt_t i=0; i<clusters.size(); i++) {
		trackCluster thisTrackCluster = clusters.at(i);
		thisCluster = At(thisTrackCluster.track)->At(thisTrackCluster.cluster);
		
		dist = diffmmXY(thisCluster, interpolatedCluster);
		if (dist < minDist) {
			minDist = dist;
			minIdx = i;
		}
	}
	return minIdx;
}

vector<Int_t> * Tracks::getTracksFromCluster(Cluster * cluster) {
	Track			 	 *	thisTrack = nullptr;
	vector<Int_t>	 *	tracksWithCluster = new vector<Int_t>;
	Int_t					idx;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		thisTrack = At(i);
		if (!thisTrack) continue;

		if (thisTrack->isClusterInTrack(cluster)) {
		   idx = thisTrack->getClusterIdx(cluster);
		   thisTrack->At(idx)->markUsed();
			tracksWithCluster->push_back(i);
		}
	}

	return tracksWithCluster;
}	

Int_t Tracks::getTrackIdxFromFirstLayerEID(Int_t eventID) {
	Cluster	* thisCluster = nullptr;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		thisCluster = At(i)->At(0);

		if (!thisCluster) continue;

		if (eventID == thisCluster->getEventID()) return i;
	}

	return -1;
}

Int_t Tracks::getTrackIdxFromCluster(Cluster * cluster) {
	Track * thisTrack = nullptr;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		thisTrack = At(i);
		if (!thisTrack) continue;

		if (thisTrack->isClusterInTrack(cluster)) {
			return i;
		}
	}
}

void Tracks::matchWithEventIDs(Hits * eventIDs) {
	// Use the Monte Carlo truth list eventIDs (x, y, layer, actual event)
	// and find the closes match for each item in list
	// then set truth eventID to the Tracks->Track->Cluster object

	Float_t		minDist = 1e5; // px
	Float_t		thisDist = 0;
	Track		 *	thisTrack = nullptr;
	Cluster	 *	thisCluster = nullptr;
	Hit		 *	thisHit = nullptr;
	Int_t			layer = -1;
	Int_t			minIdx = 0;
	Bool_t		doLoop = true;
	Int_t			nHits = eventIDs->GetEntriesFast();
	Float_t		cX, cY;
	Int_t			nClusters = 0;

	for (Int_t t=0; t<GetEntriesFast(); t++) {
		thisTrack = At(t);
		if (!thisTrack) continue;

		nClusters += thisTrack->GetEntriesFast();

		for (Int_t c=0; c<thisTrack->GetEntriesFast(); c++) {
			thisCluster = thisTrack->At(c);
			layer = thisCluster->getLayer();

			cX = thisCluster->getX();
			cY = thisCluster->getY();

			// Optimization:
			// If thisCluster+1 is also eventIDs+1, don't loop to see if we can find it
			// But instead just use minIdx++ ;-)

			minDist = diffXY(thisCluster, eventIDs->At(minIdx+1));

			if (minDist < 10) {
				minIdx++;
				doLoop = false;
			}

			else {
				doLoop = true;
				minDist = 1e5;
				minIdx = -1;
			}

			if (doLoop) {
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
			}
			
			if (minIdx >= 0 && minDist < 14) {
				thisCluster->setEventID(eventIDs->getEventID(minIdx));
				eventIDs->removeHitAt(minIdx);
			}
		}
	}

	Int_t cWithoutEventID = 0;
	for (Int_t t=0; t<GetEntriesFast(); t++) {
		for (Int_t c=0; c<At(t)->GetEntriesFast(); c++) {
			if (At(t)->getEventID(c) < 0) {
				cWithoutEventID++;
			}
		}
	}

	cout << "Number of clusters without eventID: " << cWithoutEventID << "( " << (float) cWithoutEventID / nClusters * 100 << "%)\n";
}

void Tracks::sortTrackByLayer(Int_t trackIdx) {
	// FIXME Check that this functions does what it should... It looks a bit strange!

	Int_t			lastLayer = 0;
	Bool_t		isSorted = true;
	Cluster	 *	lastCluster = nullptr;
	Track		 *	newTrack = new Track();
	Cluster	 *	nextCluster = nullptr;
	Int_t			firstLayer = 0;
	Int_t			n = GetEntriesFast(trackIdx);

	for (Int_t i=0; i<n; i++) {
		if (!At(i)) continue;
		firstLayer = i;
		break;
	}

	for (Int_t i=0; i<n; i++) {
		if (!At(trackIdx)->At(i)) continue;

		if (lastLayer > At( trackIdx )->getLayer(i))
			isSorted = false;

		lastLayer = At( trackIdx )->getLayer(i);
	}
	
	if (!isSorted) {
		vector<Int_t> sortList;
		sortList.resize(nLayers);
		for (Int_t i=0; i<nLayers; i++) sortList.at(i) = -1;
		for (Int_t i=0; i<n; i++) {
			if (!At(trackIdx)->At(i)) continue;
			sortList.at(At(trackIdx)->getLayer(i)) = i;
		}

		for (UInt_t i=0; i<sortList.size(); i++) {
			if (sortList.at(i) <0) continue;
			nextCluster = At(trackIdx)->At(sortList.at(i));
			newTrack->appendCluster(nextCluster, firstLayer);
		}

		if (newTrack->GetEntriesFast()) {
			appendTrack(newTrack);
			removeTrackAt(trackIdx);
			newTrack->clearTrack();
		}
	}
}

void Tracks::checkLayerOrientation() {
	Float_t		n[nLayers-1], x[nLayers-1], y[nLayers-1];
	Float_t		dx, dy;
	Track		 *	thisTrack = nullptr;
	Cluster	 * lastCluster = nullptr;
	Cluster	 * thisCluster = nullptr;
	Int_t			first = 0;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		thisTrack = At(i);
		if (!At(first)) first = 1;
		if (!At(first)) continue;

		lastCluster = thisTrack->At(first);
		for (Int_t j=first+1; j<thisTrack->GetEntriesFast(); j++) {
			thisCluster = thisTrack->At(j);
			dx = thisCluster->getX() - lastCluster->getX();
			dy = thisCluster->getY() - lastCluster->getY();

			x[j] += dx;
			y[j] += dy;
			n[j]++;

			lastCluster = thisCluster;
		}
	}

	cout << "The average track movement between layers is: ";
	for (Int_t i=0; i<nLayers-1; i++) {
		if (n[i]) {
			cout << Form("Layer %d: (%.1f, %.1f), ", i, x[i]/n[i], y[i]/n[i]);
		}
	}
	cout << ".\n";
}
