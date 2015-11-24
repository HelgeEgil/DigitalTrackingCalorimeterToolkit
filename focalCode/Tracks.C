#include "Tracks.h"
// #include "Cluster.h"
// #include "TFocal.h"
// #include <iostream>
// #include <TClonesArray.h>
// #include <TObject.h>

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

	cout << "Welcome to splitSharedClusters! We'll start out with " << GetEntriesFast() << " track objects.\n";

	Tracks *tracksWithMissingLayers = new Tracks();
	Cluster *interpolatedCluster;
	Float_t minDist = 1e5;
	Int_t idx = 0;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		if (At(i)->getNMissingLayers() > 0) {
			tracksWithMissingLayers->appendTrack(At(i));
		}
	}

	cout << "After looking through the tracks, I've found " << tracksWithMissingLayers->GetEntriesFast() << " tracks with missing layers.\n";

	for (Int_t layer=1; layer<nLayers; layer++) {
		vector<Int_t> * clustersThisLayer = new vector<Int_t>;
		vector<Int_t> * missingClustersThisLayer = new vector<Int_t>;

		clustersThisLayer->reserve(1000);
		missingClustersThisLayer->reserve(200);

		for (Int_t j=0; j<GetEntriesFast(); j++) {
			if (!At(j)) continue;

			if (At(j)->getClusterFromLayer(layer))
				clustersThisLayer->push_back(j);

			else
				missingClustersThisLayer->push_back(j);
		}

		cout << "In layer " << layer << ", there are in total " << clustersThisLayer->size() << "\tclusters, and " << missingClustersThisLayer->size() << "\tmissing clusters.\n";

		for (Int_t j=0; j<missingClustersThisLayer->size(); j++) {
			idx = missingClustersThisLayer->at(j);
			interpolatedCluster = At(idx)->getInterpolatedClusterAt(layer);

			if (interpolatedCluster)
				cout << "\t\tFound new interpolated cluster: " << *interpolatedCluster << endl;
		}
	}
	delete tracksWithMissingLayers;
}
