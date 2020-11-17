#ifndef findClusters_cxx
#define findClusters_cxx

#include <vector>
#include <algorithm>

#include <TStopwatch.h>
#include <TH1F.h>

#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

void Hits::findClustersFromHits(Clusters * clusters) {
// Clusters        * clusters = new Clusters(kEventsPerRun * nLayers * 10);
   vector<Int_t>   * expandedCluster = nullptr;
   vector<Int_t>   * checkedIndices = nullptr;
   vector<Int_t>   * firstHits = nullptr;
   Int_t             layerIdxFrom, layerIdxTo;

   for (Int_t layer=0; layer<nLayers; layer++) {
      layerIdxFrom = getFirstIndexOfLayer(layer);
      layerIdxTo = getLastIndexOfLayer(layer);

      if (layerIdxFrom<0) continue;
//      makeVerticalIndexOnLayer(layer); // optimization

      checkedIndices = new vector<Int_t>;
      checkedIndices->reserve(layerIdxTo - layerIdxFrom);

      for (Int_t i = layerIdxFrom; i < layerIdxTo; i++) {
         if (isItemInVector(i, checkedIndices)) continue;

         firstHits = findNeighbours(i);

         if (firstHits->size()) {
            expandedCluster = getAllNeighboursFromCluster(i, checkedIndices); 
            appendNeighboursToClusters(expandedCluster, clusters);
            delete expandedCluster;
         }
         
         delete firstHits;
         
      }
      delete checkedIndices;
   }
}

vector<Int_t> * Hits::findNeighbours(Int_t index) {
   // Find all (of max 8) neighboring hits from the single hit At(index)

   vector<Int_t>   * neighbours = new vector<Int_t>;
   Int_t             yGoal = getY(index);
   Int_t             layer = getLayer(index);
//   Int_t             idxFrom = getFirstIndexBeforeY(yGoal);
//   Int_t             idxTo = getLastIndexAfterY(yGoal);
   Int_t             idxFrom = getFirstIndexOfLayer(layer);
   Int_t             idxTo = getLastIndexOfLayer(layer);
  
   neighbours->reserve(8);

   for (Int_t j=idxFrom; j < idxTo; j++) {
      if (index == j) continue;
      if (neighbours->size() == 8) break;

      Int_t dx = abs(getX(index) - getX(j));
      Int_t dy = abs(yGoal - getY(j));
      if (dx + dy <= 1) { // no diagonal neighbors
         neighbours->push_back(j);
      }
   }
   return neighbours;
}

vector<Int_t> * Hits::getAllNeighboursFromCluster(Int_t i, vector<Int_t> *checkedIndices) {
   // TODO: Optimize this function
   // Use static arrays, and a counter instead of toCheck->empty();

   vector<Int_t>   * expandedCluster = new vector<Int_t>;
   vector<Int_t>   * toCheck = new vector<Int_t>;
   vector<Int_t>   * nextCandidates = nullptr;
   Int_t             currentCandidate = 0;
   
   expandedCluster->reserve(100);
   toCheck->reserve(100);

   expandedCluster->push_back(i);
   toCheck->push_back(i);
   
   while (!toCheck->empty()) {
      currentCandidate = toCheck->back();    
      toCheck->pop_back();
      
      checkedIndices->push_back(currentCandidate);
      
      nextCandidates = findNeighbours(currentCandidate);
      checkAndAppendAllNextCandidates(nextCandidates, checkedIndices, toCheck, expandedCluster);
      delete nextCandidates;
   }

   delete toCheck;

   return expandedCluster;
}

void Hits::appendNeighboursToClusters(vector<Int_t> *expandedCluster, Clusters *clusters) {
   Float_t  sumX = 0;
   Float_t  sumY = 0;
   Int_t    idx = 0;
   Int_t    cSize = expandedCluster->size();
   Hit *    firstHit = At(expandedCluster->at(0));
   Int_t    layerNo = firstHit->getLayer();

   Int_t maxPDG = -1000;
   Int_t maxEventID = -1000; 
   Bool_t maxSecondary = false;
   Int_t otherEventID = -1000;

   for (Int_t j=0; j<cSize; j++) {
      idx = expandedCluster->at(j);
//      if (!isSecondary(idx)) isPrimaryInCluster = true;
      sumX += getX(idx) - 0.5; // -0.5  to get
      sumY += getY(idx) - 0.5; // pixel center
      if (getPDG(idx) > maxPDG) { // use PDG and eventID of "heaviest" particle
         maxPDG = getPDG(idx);
         maxEventID = getEventID(idx);
         maxSecondary = isSecondary(idx);
      }
      if (getPDG(idx) == maxPDG && maxEventID != getEventID(idx)) {
         otherEventID = getEventID(idx);
//         cout << "Merged clusters, eventID " << otherEventID << " suppressed in layer " << getLayer(idx) << endl;
         kSuppressedClustersEventID.push_back(otherEventID);
         kSuppressedClustersLayer.push_back(getLayer(idx));
      }
   }

   clusters->appendCluster(sumX / cSize, sumY / cSize, layerNo, cSize, maxEventID, maxSecondary, maxPDG); //firstHit->getEventID(), firstHit->isSecondary(), firstHit->getPDG());
}

void Hits::checkAndAppendAllNextCandidates(vector<Int_t> * nextCandidates, vector<Int_t> *checkedIndices,
         vector<Int_t> *toCheck, vector<Int_t> *expandedCluster) {

   Bool_t   inChecked;
   Bool_t   inToCheck;
   Int_t    nextCandidate;
   
   while (!nextCandidates->empty()) {
      nextCandidate = nextCandidates->back();
      nextCandidates->pop_back();

      inChecked = isItemInVector(nextCandidate, checkedIndices);
      inToCheck = isItemInVector(nextCandidate, toCheck);

      if (!inChecked && !inToCheck) {
         expandedCluster->push_back(nextCandidate);
         toCheck->push_back(nextCandidate);
      }
   }
}

/////////////////
/////////////////

vector<Hits*> * Hits::findClustersHitMap() {
   // Not part of the standard cluster finder routine
   // return vector<Hits*>, where each Hits contains Hit objects in a single cluster
   vector<Hits*> *clusterHitMap = new vector<Hits*>;
   vector<Int_t> *checkedIndices = new vector<Int_t>;
   vector<Int_t> *expandedCluster = 0;
   checkedIndices->reserve(GetEntriesFast());

   showDebug("Looping through " << GetEntriesFast() << " hits\n");
   for (Int_t i = 0; i < GetEntriesFast(); i++) {
      if (isItemInVector(i, checkedIndices)) continue;

      vector<Int_t> * firstHits = findNeighbours(i);
      if (firstHits->size()) {
         showDebug("Finding neighbours from hit " << i << endl);
         expandedCluster = getAllNeighboursFromCluster(i, checkedIndices);
         appendExpandedClusterToClusterHitMap(expandedCluster, clusterHitMap);
         delete expandedCluster;
         delete firstHits;
      }
   }

   delete checkedIndices;
   return clusterHitMap;
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

#endif
