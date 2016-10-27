#include <vector>
#include <algorithm>

#include <TStopwatch.h>

#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

Hits::~Hits() {
   hits_.Delete();
}

void Hits::Clear(Option_t *) {
	hits_.Clear("C");
}

void Hits::appendPoint(Int_t x, Int_t y, Int_t layer, Int_t event, Float_t edep) {
   Int_t i = GetEntriesFast();

   Hit *hit = (Hit*) hits_.ConstructedAt(i);
   hit->set(x,y,layer,event, edep);
}

void Hits::appendHits(Hits *hits) {
  Int_t i = GetEntriesFast();

  for (Int_t j=0; j<hits->GetEntriesFast(); j++) {
	Hit *hit = (Hit*) hits_.ConstructedAt(i);
	hit->set(hits->At(i++));
  }
}

Int_t Hits::getI(Int_t x, Int_t y) {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (x == getX(i) && y == getY(i))
			return i;
	}
	return -1;
}

void Hits::makeLayerIndex() {
	Int_t    lastLayer = -1;
	Int_t    startOffset = 0;
	Bool_t   kStarted = true;

	if (layerIndex_.size() == 0) {
      for (Int_t i=0; i<nLayers; i++)
         layerIndex_.push_back(-1);
   }

	else {
		for (Int_t i=0; i<nLayers; i++)
		   layerIndex_.at(i) = -1;
   }

	if (!GetEntriesFast()) return;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i))
			continue;

		if (lastLayer != getLayer(i)) {
			layerIndex_.at(getLayer(i)) = i;
			lastLayer = getLayer(i);
		}
	}

	// set first indices to 0 if getLayer(0)>0
	for (UInt_t i=0; i<layerIndex_.size(); i++) {
		if (layerIndex_.at(i) == -1 && !kStarted) {
			startOffset = i;
		}
		else if (layerIndex_.at(i) != -1) {
			kStarted = true;
		}
	}

	if (startOffset>0) {
		for (Int_t i=0; i<=startOffset; i++)
			layerIndex_.at(i) = 0;
	}
}

Int_t Hits::getFirstIndexOfLayer(UInt_t layer) {
   if (layerIndex_.size() == 0) return -1;
	if (layerIndex_.size() < layer) return -1;

   return layerIndex_.at(layer);
}

Int_t Hits::getLastIndexOfLayer(UInt_t layer) {
   for (Int_t i=layer+1; i<nLayers; i++) {
      if (getFirstIndexOfLayer(i) != -1 ) {
         return getFirstIndexOfLayer(i);
      }
   }
   return GetEntriesFast(); // no hits beyond layer
}

Int_t Hits::getLastActiveLayer() {
	Int_t lastActiveLayer = 0;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i))
			continue;

		if (getLayer(i) > lastActiveLayer)
			lastActiveLayer = getLayer(i);
	}
	return lastActiveLayer;
}

void Hits::makeVerticalIndexOnLayer(Int_t layer) {
	// Run this command when a new layer is to be used
	// It can only hold a single layer

   Int_t    layerIdxFrom, layerIdxTo;
   Int_t    lastY = -1;
   Int_t    startOffset = 0;
   Bool_t   kStarted = false;

	if (!layerIndex_.size())
		makeLayerIndex();

	layerIdxFrom = getFirstIndexOfLayer(layer);
	layerIdxTo = getLastIndexOfLayer(layer);

	if (verticalIndexOfLayer_.size() == 0)
		for (Int_t i=0; i<2*ny+1; i++) verticalIndexOfLayer_.push_back(-1);
	else
		for (Int_t i=0; i<2*ny+1; i++) verticalIndexOfLayer_.at(i) = -1;

	if (!GetEntriesFast()) return;

	for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
		if (!At(i))
			continue;

		if (lastY != getY(i)) {
			verticalIndexOfLayer_.at(getY(i)) = i;
			lastY = getY(i);
		}
	}

	// set first indices to  0 if getLayer(0)>0
	for (UInt_t i=0; i<verticalIndexOfLayer_.size(); i++) {
		if (verticalIndexOfLayer_.at(i) == -1 && !kStarted) {
			startOffset = i;
		}
		else if (verticalIndexOfLayer_.at(i) != -1) {
			kStarted = true;
		}
	}

	if (startOffset>0) {
		for (Int_t i=0; i<=startOffset; i++)
			verticalIndexOfLayer_.at(i) = 0;
	}
}

Int_t Hits::getFirstIndexBeforeY(Int_t y) {
	Int_t idx = 0;

	if (y==0) return idx;

	for (Int_t i=y-1; i>=0; i--) {
		idx = verticalIndexOfLayer_.at(i);
		if (idx>=0) break;
	}

	if (idx == -1) idx = 0;
	return idx;
}

Int_t Hits::getLastIndexAfterY(Int_t y) {
	Int_t idx = GetEntriesFast()-1;

	if (y>2*ny-2) return idx;

	for (Int_t i=y+2; i<2*ny; i++) {
		idx = verticalIndexOfLayer_.at(i);
		if (idx>=0) break;
	}

	if (idx == -1) idx = GetEntriesFast()-1;
	return idx;
}
