#include "TrackerFrame.h"
// #include <iostream>

using namespace std;

TrackerFrame::TrackerFrame() : trackerFrame_("Layer", nTrackers) {
	for (Int_t i=0; i<nTrackers; i++) {
		// make new layer
		new(trackerFrame_[i]) Layer(i, kTracker, kMC);
	}
}

TrackerFrame::~TrackerFrame() {
	trackerFrame_.Delete();
}

Hits * TrackerFrame::findHits() {
	Hits *hits = new Hits();

	for (Int_t layer=0; layer<nTrackers; layer++) {
		At(layer)->findHits(hits);
	}
	return hits;
}

void TrackerFrame::diffuseFrame() {
	for (Int_t layer=0; layer<nTrackers; layer++) {
		At(layer)->diffuseLayer();
	}
}

void TrackerFrame::Reset() {
	for (Int_t layer=0; layer<nTrackers; layer++) {
		At(layer)->Reset();
	}
}
