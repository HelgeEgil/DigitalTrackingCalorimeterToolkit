#include "CalorimeterFrame.h"
#include "Layer.h"
// #include <iostream>

using namespace std;

CalorimeterFrame::CalorimeterFrame() : calorimeterFrame_("Layer", nLayers) {
	for (Int_t i=0; i<nLayers; i++) {
		// make new layer
		new(calorimeterFrame_[i]) Layer(i, kCalorimeter, kMC);
	}
}

CalorimeterFrame::~CalorimeterFrame() {
	calorimeterFrame_.Delete();
}

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
