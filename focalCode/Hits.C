#include "Hits.h"
#include "Clusters.h"
#include "Tools.h"
#include <vector>
#include <algorithm>

// #include <iostream>
// #include <TClonesArray.h>
// #include <TObject.h>

// ClassImp(HitCollection)

Hits::~Hits() {
   // Destructor
   clearHits();
}

void Hits::appendPoint(Int_t x, Int_t y, Int_t layer, Int_t event) {
   Int_t i = GetEntriesFast();
   Hit *hit = (Hit*) hits_.ConstructedAt(i);
   hit->set(x,y,layer,event);
}

void Hits::appendHits(Hits *hits) {
  Int_t i = GetEntriesFast();
  for (Int_t j=0; j<hits->GetEntriesFast(); j++) {
	Hit *hit = (Hit*) hits_.ConstructedAt(i);
	hit->set(hits->At(i));
	i++;
  }
}

Int_t Hits::getI(Int_t x, Int_t y) {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (x == getX(i) && y == getY(i))
			return i;
	}
	return -1;
}

Clusters * Hits::findClustersFromHits() {
	Clusters *clusters = new Clusters(kEventsPerRun * 20);

	vector<Int_t> *checkedIndices = new vector<Int_t>;
	vector<Int_t> *expandedCluster = 0;
	checkedIndices->reserve(GetEntriesFast());

	for (Int_t i = 0; i < GetEntriesFast(); i++) {
		if (isItemInVector(i, checkedIndices)) continue;

		vector<Int_t> firstHits = findNeighbours(i);
		if (firstHits.size()) {
			expandedCluster = findExpandedCluster(i, checkedIndices);
			appendExpandedClusterToClusters(expandedCluster, clusters);
			delete expandedCluster;
		}
	}
	delete checkedIndices;
// 	delete expandedCluster;

	return clusters;
}

vector<Hits*> * Hits::findClustersHitMap() {
	vector<Hits*> *clusterHitMap = new vector<Hits*>;
	vector<Int_t> *checkedIndices = new vector<Int_t>;
	vector<Int_t> *expandedCluster = 0;
	checkedIndices->reserve(GetEntriesFast());

	for (Int_t i = 0; i < GetEntriesFast(); i++) {
		if (isItemInVector(i, checkedIndices)) continue;

		vector<Int_t> firstHits = findNeighbours(i);
		if (firstHits.size()) {
			expandedCluster = findExpandedCluster(i, checkedIndices);
			appendExpandedClusterToClusterHitMap(expandedCluster, clusterHitMap);
			delete expandedCluster;
		}
	}

	delete checkedIndices;
	return clusterHitMap;
}


void Hits::appendExpandedClusterToClusters(vector<Int_t> *expandedCluster, Clusters *clusters) {
	Int_t cSize = expandedCluster->size();
	Float_t sumX = 0; Float_t sumY = 0;
	Int_t layerNo = getLayer(expandedCluster->at(0));

	for (Int_t j=0; j<cSize; j++) {
		Int_t idx = expandedCluster->at(j);
		sumX += getX(idx) - 0.5; // -0.5  to get
		sumY += getY(idx) - 0.5; // pixel center
	}
	clusters->appendCluster(sumX / cSize, sumY / cSize, layerNo, cSize);
}

void Hits::appendExpandedClusterToClusterHitMap(vector<Int_t> *expandedCluster, vector<Hits*> *clusterHitMap) {
	Float_t sumX = 0; Float_t sumY = 0;
	Float_t x; Float_t y;
	Int_t idx;
	Int_t cSize = expandedCluster->size();
	
	for (Int_t j=0; j<cSize; j++) {
		Int_t idx = expandedCluster->at(j);
		sumX += getX(idx) - 0.5;
		sumY += getY(idx) - 0.5;
	}

	// center cluster hitmap
	Int_t offsetForHistogramX = 5 - sumX/cSize;
	Int_t offsetForHistogramY = 5 - sumY/cSize;

	clusterHitMap->push_back(new Hits(cSize));

	Int_t cidx = clusterHitMap->size() - 1;	
	for (Int_t j=0; j<cSize; j++) {
		idx = expandedCluster->at(j);
		x = offsetForHistogramX + getX(idx);
		y = offsetForHistogramY + getY(idx);

		clusterHitMap->at(cidx)->appendPoint(x, y);
	}
}

vector<Int_t> Hits::findNeighbours(Int_t index) {
	vector<Int_t> neighbours;
	neighbours.reserve(8);
	Int_t xGoal = getX(index);
	Int_t yGoal = getY(index);
	Int_t zGoal = getLayer(index);

	for (Int_t j=0; j < GetEntriesFast(); j++) {
		if (neighbours.size() == 8 || getLayer(j) > zGoal) break;
		if (getLayer(j) < zGoal || index == j) continue;

		if (abs(xGoal - getX(j)) <= 1 && abs(yGoal - getY(j)) <= 1)
			neighbours.push_back(j);
	}
	return neighbours;
}

vector<Int_t> * Hits::findExpandedCluster(Int_t i, vector<Int_t> *checkedIndices) {
	vector<Int_t> *expandedCluster = new vector<Int_t>;
	vector<Int_t> *toCheck = new vector<Int_t>;
	expandedCluster->reserve(80);
	toCheck->reserve(80);

	expandedCluster->push_back(i);
	toCheck->push_back(i);
	
	while (!toCheck->empty()) {
		Int_t currentCandidate = toCheck->back();		
		toCheck->pop_back();
		
		checkedIndices->push_back(currentCandidate);
		
		vector<Int_t> nextCandidates = findNeighbours(currentCandidate);
		checkAndAppendAllNextCandidates(nextCandidates, checkedIndices, toCheck, expandedCluster);
	}

	delete toCheck;

	return expandedCluster;
}

void Hits::checkAndAppendAllNextCandidates(vector<Int_t> nextCandidates, vector<Int_t> *checkedIndices,
			vector<Int_t> *toCheck, vector<Int_t> *expandedCluster) {

	Bool_t inChecked;
	Bool_t inToCheck;
	Int_t nextCandidate;
	
	while (!nextCandidates.empty()) {
		nextCandidate = nextCandidates.back();
		nextCandidates.pop_back();

		inChecked = isItemInVector(nextCandidate, checkedIndices);
		inToCheck = isItemInVector(nextCandidate, toCheck);

		if (!inChecked && !inToCheck) {
			expandedCluster->push_back(nextCandidate);
			toCheck->push_back(nextCandidate);
		}
	}
}

void Hits::makeLayerIndex() {
  
	if (layerIndex_.size() == 0)
		for (Int_t i=0; i<nLayers; i++) layerIndex_.push_back(-1);
	else
		for (Int_t i=0; i<nLayers; i++) layerIndex_.at(i) = -1;

	if (!GetEntriesFast()) return;

	Int_t lastLayer = -1;
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i))
			continue;

		if (lastLayer != getLayer(i)) {
			layerIndex_.at(getLayer(i)) = i;
			lastLayer = getLayer(i);
		}
	}

	// set first indices to  0 if getLayer(0)>0
	Int_t startOffset = 0;
	Bool_t kStarted = kFALSE;
	for (UInt_t i=0; i<layerIndex_.size(); i++) {
		if (layerIndex_.at(i) == -1 && !kStarted) {
			startOffset = i;
		}
		else if (layerIndex_.at(i) != -1) {
			kStarted = kTRUE;
		}
	}

	if (startOffset>0) {
		for (Int_t i=0; i<=startOffset; i++)
			layerIndex_.at(i) = 0;
	}
}

Int_t Hits::getFirstIndexOfLayer(UInt_t layer) {
   // return index in clusters based on fLayerIndex
   // Return value -1 == no cluster in that layer

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
