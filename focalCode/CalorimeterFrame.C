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

void CalorimeterFrame::Clear(Option_t *) {
	calorimeterFrame_.Clear("C");
}

Hits * CalorimeterFrame::findHits() {
	Hits *hits = new Hits();

	Int_t nNoHits = 0;
	Bool_t isHits = false;
	for (Int_t layer=0; layer<nLayers; layer++) {
		isHits = At(layer)->findHits(hits);
		
		if (!isHits) nNoHits++;
		else nNoHits = 0;
		if (nNoHits>2) break;
	}
	
	hits->makeLayerIndex();
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
