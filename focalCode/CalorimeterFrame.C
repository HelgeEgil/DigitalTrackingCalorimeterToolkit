#include "CalorimeterFrame.h"
#include "Layer.h"
// #include <iostream>

using namespace std;

Hits * CalorimeterFrame::findHits() {
	Hits *hits = new Hits();

	for (Int_t layer=0; layer<nLayers; layer++) {
		At(layer)->findHits(hits);
	}
	return hits;
}

void CalorimeterFrame::diffuseFrame() {
	for (Int_t layer=0; layer<nLayers; layer++) {
		At(layer)->diffuseLayer();
	}
}

void CalorimeterFrame::Reset() {
	for (Int_t layer=0; layer<nLayers; layer++) {
		At(layer)->Reset();
	}
}
