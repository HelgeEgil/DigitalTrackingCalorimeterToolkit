#include "Hits.h"
#include "Clusters.h"
#include <vector>
// #include <iostream>
// #include <TClonesArray.h>
// #include <TObject.h>

// ClassImp(HitCollection)

Hits::~Hits() {
   // Destructor
   Clear();
}

void Hits::appendPoint(Int_t x, Int_t y, Int_t layer, Int_t event) {
   Int_t i = GetEntriesFast();
   Hit *hit = (Hit*) hits_.ConstructedAt(i);
   hit->set(x,y,layer,event);
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
	vector<Int_t> checkedIndices;
	checkedIndices.reserve(GetEntriesFast());
	Int_t layerNo = -999;

	for (Int_t i = 0; i < GetEntriesFast(); i++) {
		if (layerNo == -999) layerNo = getLayer(i);
		if (isItemInVector(i, checkedIndices)) continue;

		vector<Int_t> firstHits = findNeighbours(i);
		if (firstHits.size()) {
			vector<Int_t> expandedCluster = findExpandedCluster(i, checkedIndices);
			appendExpandedClusterToClusters(expandedCluster, clusters);
		}
	}
	return clusters;
}

void Hits::appendExpandedClusterToClusters(vector<Int_t> expandedCluster, Clusters *clusters) {
	Int_t cSize = expandedCluster.size();
	Float_t sumX = 0; Float_t sumY = 0;
	Int_t layerNo = getLayer(expandedCluster.at(0));

	for (Int_t j=0; j<cSize; j++) {
		Int_t idx = expandedCluster.at(j);
		sumX += getX(idx) - 0.5; // -0.5  to get
		sumY += getY(idx) - 0.5; // pixel center
	}
	clusters->appendCluster(sumX / cSize, sumY / cSize, layerNo, cSize);
}

void Hits::appendExpandedClusterToClusterHitMap(vector<Int_t> expandedCluster, vector<Hits*> *clusterHitMap) {
	Float_t sumX = 0; Float_t sumY = 0;
	Float_t x; Float_t y;
	Int_t idx;
	Int_t cSize = expandedCluster.size();

	for (Int_t j=0; j<cSize; j++) {
		Int_t idx = expandedCluster.at(j);
		sumX += getX(idx) - 0.5; // -0.5  to get
		sumY += getY(idx) - 0.5; // pixel center
	}

	// center cluster hitmap
	Int_t offsetForHistogramX = 5 - sumX/cSize;
	Int_t offsetForHistogramY = 5 - sumY/cSize;

	clusterHitMap->push_back(new Hits(cSize));

	Int_t cidx = clusterHitMap->size() - 1;
	for (Int_t j=0; j<cSize; j++) {
		idx = expandedCluster.at(j);
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
	Int_t jFrom = 0;
	Int_t jTo = GetEntriesFast();

	/* TODO optimization
	 *Find close starting point and end point
	 * Find scan direction and correct algorithm
	 * start on top
	while (true) {
		if (yGoal <= 2*ny-3) {
			// first two columns are defaulted to starting pos. 0
			jFrom = 0;
			break;
		}

		else {
			// Find first index where y <= yGoal + 2
			if (hits->GetY(jFrom) <= yGoal + 2) break;
			else if (jFrom <= 0) break;
			jFrom--;
		}
	}

	while (true) {
		if (yGoal >= 2*ny - 3) {
			jTo = hit_size;
			break;
		}

		else {
			// Find first index where x >= xGoal + 2
			if (hits->GetY(jTo) >= yGoal + 2) break;
			else if (jTo >= hit_size - 1) break;
			jTo++;
		}
	}
	*/

	for (Int_t j=jFrom; j < jTo; j++) {
		if (neighbours.size() == 8) break;
		if (index != j) {
			if (abs(xGoal - getX(j)) <= 1 && abs(yGoal - getY(j)) <= 1)
				neighbours.push_back(j);
		}
	}
	return neighbours;
}

vector<Int_t> Hits::findExpandedCluster(Int_t i, vector<Int_t> checkedIndices) {
	vector<Int_t> expandedCluster;
	vector<Int_t> toCheck;
	expandedCluster.reserve(80);
	toCheck.reserve(80);

	expandedCluster.push_back(i);
	toCheck.push_back(i);
	
	while (!toCheck.empty()) {
		Int_t currentCandidate = toCheck.back();
		toCheck.pop_back();
		checkedIndices.push_back(currentCandidate);

		vector<Int_t> nextCandidates = findNeighbours(currentCandidate);
		checkAndAppendAllNextCandidates(nextCandidates, checkedIndices, toCheck, expandedCluster);
	}

	return expandedCluster;
}

void Hits::checkAndAppendAllNextCandidates(vector<Int_t> nextCandidates, vector<Int_t> checkedIndices,
			vector<Int_t> toCheck, vector<Int_t> expandedCluster) {

	while (!nextCandidates.empty()) {
		Int_t nextCandidate = nextCandidates.back();
		nextCandidates.pop_back();

		Bool_t isValidCandidate;
		if (!isItemInVector(nextCandidate, checkedIndices) && !isItemInVector(nextCandidate, toCheck))
			isValidCandidate = kTRUE;
		else
			isValidCandidate = kFALSE;

		if (isValidCandidate) {
				expandedCluster.push_back(nextCandidate);
				toCheck.push_back(nextCandidate);
		}
	}
}

vector<Hits*> * Hits::findClustersHitMap(Int_t nRuns) {
	Clusters *clusters = new Clusters(nRuns * 20);
	vector<Hits*> *clusterHitMap = new vector<Hits*>;
	vector<Int_t> checkedIndices;
	checkedIndices.reserve(GetEntriesFast());
	Int_t layerNo = -999;

	for (Int_t i = 0; i < GetEntriesFast(); i++) {
		if (layerNo == -999) layerNo = getLayer(i);
		if (isItemInVector(i, checkedIndices)) continue;

		vector<Int_t> firstHits = findNeighbours(i);
		if (firstHits.size()) {
			vector<Int_t> expandedCluster = findExpandedCluster(i, checkedIndices);
			appendExpandedClusterToClusterHitMap(expandedCluster, clusterHitMap);
		}
	}
	return clusterHitMap;
}


Bool_t Hits::isItemInVector(Int_t i, vector<Int_t> vector) {
	return find(vector.begin(), vector.end(), i) != vector.end();
}
