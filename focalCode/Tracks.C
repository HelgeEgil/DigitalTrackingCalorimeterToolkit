#include "Tracks.h"
#include "Cluster.h"
#include "Clusters.h"
// #include "TFocal.h"
#include <iostream>
#include <TClonesArray.h>
#include <TObject.h>
#include "Tools.h"
#include <TCanvas.h>
#include <TView.h>
#include <TAttMarker.h>
#include <TAttLine.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>

using namespace std;

Tracks::~Tracks() {
   // Destructor
   clearTracks();
}

void Tracks::appendTrack(Track *copyTrack, Int_t startOffset /* default 0 */) {
   Int_t newIdx = tracks_.GetEntriesFast();
   Track *track = (Track*) tracks_.ConstructedAt(newIdx);

   for (Int_t i=0; i<copyTrack->GetEntriesFast(); i++) {
   	if(!copyTrack->At(i))
   		continue;

      track->appendCluster(copyTrack->At(i), startOffset);
   }
}

void Tracks::appendClustersWithoutTrack(TClonesArray *clustersWithoutTrack) {
	Int_t idxFrom = clustersWithoutTrack_.GetEntriesFast();
	Cluster *newCluster = 0;

	for (Int_t i=0; i<clustersWithoutTrack->GetEntriesFast(); i++) {
		if (!clustersWithoutTrack->At(i))
			continue;

		newCluster = (Cluster*) clustersWithoutTrack_.ConstructedAt(idxFrom + i);
		newCluster->set((Cluster*) clustersWithoutTrack->At(i));
	}
}

void Tracks::extrapolateToLayer0() {
   // do this for all tracks in collection
   for (Int_t i=0; i<tracks_.GetEntriesFast(); i++) {
      if (!At(i)->At(0)) At(i)->extrapolateToLayer0();
   }
}

void Tracks::splitSharedClusters() {
	// to be run as early as possible before joining Tracks* objects

	cout << "Welcome to splitSharedClusters! We'll start out with " << GetEntriesFast() << " track objects.\n\n";

	Cluster *interpolatedCluster;
	Clusters* interpolatedClusters = new Clusters();
	Clusters* interpolatedClosestClusters = new Clusters();

	Float_t dist, minDist = 1e5;
	Int_t minIdx = 0, idx = 0, clusterIdx;
	Bool_t isMissing;

	for (Int_t layer=1; layer<nLayers; layer++) {
		
		vector<trackCluster> clustersThisLayer;
		vector<trackCluster> missingClustersThisLayer;
		
		for (Int_t i=0; i<GetEntriesFast(); i++) {
		Track *thisTrack = At(i);
		
		if (!thisTrack) continue;
		if (thisTrack->getLastLayer() < layer) continue;
		
		trackCluster thisCluster;
		thisCluster.track = i;
		
		if (thisTrack->hasLayer(layer)) {
			thisCluster.cluster = thisTrack->getClusterFromLayer(layer);
			clustersThisLayer.push_back(thisCluster);
		}
		
		else
			missingClustersThisLayer.push_back(thisCluster);
		}
		
		cout << "Layer " << layer << ": Found " << clustersThisLayer.size() << " clusters and " << missingClustersThisLayer.size() << " missing clusters.\n";

		for (UInt_t i=0; i<missingClustersThisLayer.size(); i++) {
			Int_t trackIdx = missingClustersThisLayer.at(i).track;
			Track *missingTrack = At(trackIdx);
			if (!missingTrack ) continue;
		
			interpolatedCluster = missingTrack->getInterpolatedClusterAt(layer);
				if (interpolatedCluster) {
				minIdx = getClosestCluster(clustersThisLayer, interpolatedCluster);
				interpolatedClusters->appendCluster(interpolatedCluster);
				
				trackCluster closestTC= clustersThisLayer.at(minIdx);
				Track *closestTrack = At(closestTC.track);
				Cluster *closestCluster = closestTrack->At(closestTC.cluster);
				
				minDist = diffmm(closestCluster, interpolatedCluster);
				
				Float_t clusterRadius = closestCluster->getRadiusmm();
				if (minDist < clusterRadius * 2) {
					cout << "--> Minimum distance (" << minDist << " mm) is smaller than the cluster radius of " << clusterRadius << " mm (csize = " << closestCluster->getSize() << ")!\n";
					
					// 1) Find interpolated cluster no 2
					Cluster *interpolatedClosestCluster = closestTrack->getInterpolatedClusterAt(layer);
					if (interpolatedClosestCluster) {
						cout << "Distance between interpolated closest cluster " << *interpolatedClosestCluster
							<< " and actual closest cluster " << *closestCluster << " is " << diffmm(interpolatedClosestCluster, closestCluster) << " mm.\n";
						cout << "Distance between interpolated missing cluster " << *interpolatedClosestCluster
							<< " and interpolated closest cluster " << *interpolatedCluster << " is " << diffmm(interpolatedClosestCluster, interpolatedCluster) << " mm.\n";
						
						interpolatedClosestClusters->appendCluster(interpolatedClosestCluster);
						
						Float_t sizeFactor = minDist / closestCluster->getRadiusmm() + 1;
						Float_t x = ( closestCluster->getX() + interpolatedCluster->getX()) / 2;
						Float_t y = ( closestCluster->getY() + interpolatedCluster->getY()) / 2;
						Float_t size = closestCluster->getSize() / sizeFactor;
						missingTrack->appendCluster(new Cluster(x, y, layer, size));
						
						sortTrackByLayer(trackIdx);
						// FIXME: When not tired, do validation of this method
					}
				}
			}
		}
	}
}

Int_t Tracks::getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster) {
	Float_t minDist = 1e5, dist;
	Int_t minIdx = 0;

	for (UInt_t i=0; i<clusters.size(); i++) {
		trackCluster thisTrackCluster = clusters.at(i);
		Cluster *thisCluster = At(thisTrackCluster.track)->At(thisTrackCluster.cluster);
		
		dist = diffmm(thisCluster, interpolatedCluster);
		if (dist < minDist) {
			minDist = dist;
			minIdx = i;
		}
	}
	return minIdx;
}

void Tracks::sortTrackByLayer(Int_t trackIdx) {
	Int_t lastLayer = 0;
	Bool_t isSorted = true;
	Cluster *lastCluster;
	
	Int_t n = GetEntriesFast( trackIdx );
	
	for (Int_t i=0; i<n; i++) {
		if (!At( trackIdx )->At(i)) continue;
		if (lastLayer > At( trackIdx )->getLayer(i)) {
			isSorted = false;
		}
		lastLayer = At( trackIdx )->getLayer(i);
	}
	
	if (!isSorted) {
		vector<Int_t> sortList;
		sortList.resize(nLayers);
		for (Int_t i=0; i<nLayers; i++) sortList.at(i) = -1;
		for (Int_t i=0; i<n; i++) {
			if (!At( trackIdx )->At(i)) continue;
			sortList.at(At( trackIdx)->getLayer(i)) = i;
		}
		
		cout << "A new sorting of the track with layers: ";
		for (Int_t i=0; i<n; i++) {
			if (!At(trackIdx)->At(i)) continue;
			cout << At(trackIdx)->getLayer(i) << ", ";
		}
		cout << " is with indices: ";
		for (UInt_t i=0; i<sortList.size(); i++) cout << sortList.at(i) << ", ";
		cout << ".\n";
		
		Track *newTrack = new Track();
		Cluster *nextCluster = 0;
		for (UInt_t i=0; i<sortList.size(); i++) {
			if (sortList.at(i) <0) continue;
			nextCluster = At(trackIdx)->At(sortList.at(i));
			newTrack->appendCluster(nextCluster);
		}
		if (newTrack->GetEntriesFast()) {
			appendTrack(newTrack);
			removeTrackAt(trackIdx);
		}
		
		cout << "New Track: ";
		for (Int_t i=0; i<newTrack->GetEntriesFast(); i++) { if (!newTrack->At(i)) continue; cout << *newTrack->At(i) << ", "; }
		cout << endl;
	}
}