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

	for (Int_t i=0; i<clustersWithoutTrack->GetEntriesFast(); i++) {
		if (!clustersWithoutTrack->At(i))
			continue;

		Cluster *newCluster = (Cluster*) clustersWithoutTrack_.ConstructedAt(idxFrom + i);
		newCluster->set((Cluster*) clustersWithoutTrack->At(i));
	}
}

void Tracks::extrapolateToLayer0() {
   // do this for all tracks in collection
   for (Int_t i=0; i<tracks_.GetEntriesFast(); i++) {
      if (!At(i)->At(0)) At(i)->extrapolateToLayer0();
   }
}
