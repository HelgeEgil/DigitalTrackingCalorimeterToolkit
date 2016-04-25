#include "CalorimeterFrame.h"
#include "Layer.h"
#include "Constants.h"
#include "MaterialConstants.h"
#include <TH2F.h>
#include <TRandom3.h>
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

Hits * CalorimeterFrame::findHits(Int_t eventID) {
	Hits *hits = new Hits();

	Int_t nNoHits = 0;
	Bool_t isHits = false;
	for (Int_t layer=0; layer<nLayers; layer++) {
		isHits = At(layer)->findHits(hits);

		if (!isHits) nNoHits++;
		else nNoHits = 0;
		if (nNoHits>2) break;
	}
	
	if (eventID > 0) {
		for (Int_t i=0; i<hits->GetEntriesFast(); i++) {
			hits->At(i)->setEventID(eventID);
		}
	}

	hits->makeLayerIndex();
	return hits;
}

void CalorimeterFrame::diffuseFrame(TRandom3 *gRandom) {
	Int_t nHitsInLayer = 1;
	
	for (Int_t layer=0; layer<nLayers; layer++) {
		if (nHitsInLayer)
			nHitsInLayer = At(layer)->diffuseLayer(gRandom);
		
		else break;
	}
}

void CalorimeterFrame::Reset() {
	for (Int_t layer=0; layer<nLayers; layer++) {
		At(layer)->Reset();
	}
}
