#include "Tracks.h"
#include "Cluster.h"
#include "Clusters.h"
// #include "TFocal.h"
#include <iostream>
// #include <TClonesArray.h>
// #include <TObject.h>
#include "Tools.h"

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
  Float_t dist, minDist = 1e5;
  Int_t minIdx = 0, idx = 0, clusterIdx;
  Bool_t isMissing;

  for (Int_t layer=1; layer<nLayers; layer++) {
	  
    vector<trackCluster> clustersThisLayer;
    vector<trackCluster> missingClustersThisLayer;
    
    // 1) Find missing clusters in layer
    // 2) Find all other clusters in layer
    
    for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      
      trackCluster thisCluster;
      thisCluster.track = i;

      Track *thisTrack = At(i);

      if (thisTrack->hasLayer(layer)) {
	thisCluster.cluster = thisTrack->getClusterFromLayer(layer);
	clustersThisLayer.push_back(thisCluster);
      }
      
      else {
	missingClustersThisLayer.push_back(thisCluster);
      }
    }
    
    cout << "Layer " << layer << ": Found " << clustersThisLayer.size() << " clusters and " << missingClustersThisLayer.size() << " missing clusters.\n";

    // 3) For all missing clusters, see if any of the clusters are a good match
    // 4) Good match: Small distance and matching enlarged cluster size.
    // 5) If good match: Split them (?) and remove from missingClustersThisLayer.

    for (UInt_t i=0; i<missingClustersThisLayer.size(); i++) {
      Track *thisTrack = At(missingClustersThisLayer.at(i).track);
      if (!thisTrack) continue;
      
      interpolatedCluster = thisTrack->getInterpolatedClusterAt(layer);
      if (interpolatedCluster) {
	cout << "\t Found new interpolated cluster: " << *interpolatedCluster << endl;
	cout << "\t\t From track: ";
	for (Int_t i=0; i<At(idx)->GetEntriesFast(); i++) {
	  if (!At(idx)->At(i)) continue;
	  cout << *At(idx)->At(i) << ", ";
	}
	cout << endl;
	
	// does interpolatedCluster match any in clustersThisLayer?
	
	minDist = 1e5;
	for (UInt_t i=0; i<clustersThisLayer.size(); i++) {
	  trackCluster thisTrackCluster = clustersThisLayer.at(i);
	  Cluster *thisCluster = At(thisTrackCluster.track)->At(thisTrackCluster.cluster);
	  
	  dist = diffmm(thisCluster, interpolatedCluster);
	  if (dist < minDist) {
	    minDist = dist;
	    minIdx = i;
	  }
	}
	
	trackCluster thisTrackCluster = clustersThisLayer.at(minIdx);
	Cluster *thisCluster = At(thisTrackCluster.track)->At(thisTrackCluster.cluster);
	
	cout << "Found minimum distance between " << *interpolatedCluster << " and " << *thisCluster << " of " << minDist << " mm.\n";
      }
    }
  }
}
