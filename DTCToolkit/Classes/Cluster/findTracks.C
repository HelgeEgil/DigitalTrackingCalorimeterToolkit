#ifndef findTracks_cxx
#define findTracks_cxx

#include <vector>
#include<fstream>
#include<iostream>

#include <string>
#include <sstream>


//#include <TClonesArray.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Node.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

#include <TH2.h>
#include <TH3.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TPaveStats.h>
#include <TView.h>
#include <TLeaf.h>
#include <TArrow.h>
#include <TF1.h>
#include <Math/ProbFunc.h>


using namespace std;
using namespace DTC;

Tracks * Clusters::findTracksWithRecursiveWeighting() {
   Track     * track = nullptr;
   Tracks    * tracks = new Tracks(50 * 200 * kEventsPerRun);
   Tracks    * trackerTracks = new Tracks(kEventsPerRun);
   Int_t       spotSize = 33;
   Float_t thisMaxTrackScore = kMaxTrackScore;
   makeLayerIndex();
   kMCSFactor = 25;

   thisMaxTrackScore = 0.275;
   if (kHelium) thisMaxTrackScore = 0.175;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      appendClusterWithoutTrack(At(i));
   }

   Int_t FirstLayer = nLayers-1;
   Int_t LastLayer = 0;
      
   if (kDoTrackerPropagation) {
      // Make list of all seed pairs in TRACKER layers where angle a,a,b < 100 mrad
      // Then use the extrapolated vector to match in first interior layer
      
      Cluster * firstTrackerCluster = nullptr;
      Cluster * secondTrackerCluster = nullptr;
      Cluster * firstLayerCluster = nullptr;
      Float_t  angle, bestAngle;
      Int_t    bestIdxI, bestIdxJ;
      Int_t    firstLayerIdxFrom = getFirstIndexOfLayer(2); 
      Int_t    firstLayerIdxTo = getLastIndexOfLayer(2);
      vector<Int_t> * firstTrackerClusters = new vector<Int_t>;
      vector<Int_t> * secondTrackerClusters = new vector<Int_t>;
      
      findSeeds(firstTrackerClusters, 0, false);
      for (Int_t k=firstLayerIdxFrom; k<firstLayerIdxTo; k++) {
         // Find straightest track for each hit in first interior layer

         firstLayerCluster = At(k);
         if (!firstLayerCluster) continue;

         bestIdxI = 0; bestIdxJ = 0;
         bestAngle = 1e5;

         for (UInt_t i=0; i<firstTrackerClusters->size(); i++) {
            // Loop over clusters in first tracker layer

            firstTrackerCluster = At(firstTrackerClusters->at(i));
            if (firstTrackerCluster->isUsed()) continue;
            
            secondTrackerClusters->clear();
            findClustersFromSeedInLayer(firstTrackerCluster, 1, secondTrackerClusters);

            for (UInt_t j=0; j<secondTrackerClusters->size(); j++) {
               secondTrackerCluster = At(secondTrackerClusters->at(j));
               if (secondTrackerCluster->isUsed()) continue;

               angle = getDotProductAngle(firstTrackerCluster, secondTrackerCluster, firstLayerCluster);
               if (angle < bestAngle) {
                  bestAngle = angle; 
                  bestIdxI = firstTrackerClusters->at(i); 
                  bestIdxJ = secondTrackerClusters->at(j);
               }
            }
         }

         if (bestAngle>0.2) continue; // no hit on first 

         firstTrackerCluster = At(bestIdxI);
         secondTrackerCluster = At(bestIdxJ);

         Track * trackerTrack = new Track();
         trackerTrack->appendCluster(firstTrackerCluster);
         trackerTrack->appendCluster(secondTrackerCluster);
         
         markUsedClusters(trackerTrack);

         trackerTrack->appendCluster(firstLayerCluster);
         trackerTracks->appendTrack(trackerTrack);

         delete trackerTrack;
      }
   }

   for (Int_t s=FirstLayer; s>=LastLayer; s--) {
      Cluster   * nextCluster = nullptr;
      Float_t     clusterScore;
      Node      * nextNode = nullptr;
      Node      * seedNode = nullptr;
      Cluster   * seed = nullptr;
      vector<Int_t> * nextClusters = new vector<Int_t>;
      vector<Int_t> * seeds = new vector<Int_t>;
      nextClusters->reserve(200);

      findSeeds(seeds, s, false);
      for (UInt_t i=0; i<seeds->size(); i++) {
         seed = At(seeds->at(i)); 
         findNearestClustersInNextLayer(seed, nextClusters); 
         seedNode = new Node(nullptr, seed, 0); 

         // Explore all identified clusters
         for (UInt_t j=0; j<nextClusters->size(); j++) {
            nextCluster = At(nextClusters->at(j)); 
            clusterScore = seedNode->getNextScore(nextCluster); 
//            clusterScore /= 10;
//            if (std::isnan(clusterScore)) cout << "clusterScore at " << *nextCluster << " isNan!\n";
            
            if (clusterScore < thisMaxTrackScore) {
               nextNode = new Node(seedNode, nextCluster, clusterScore); 
               seedNode->addChild(nextNode);
            }
         }

         seedNode->markExplored();
         vector<Node*> * endNodes = new vector<Node*>;
         endNodes->reserve(200* 50 * kEventsPerRun);
         seedNode->getUnexploredEndNodes(endNodes); 
         doRecursiveWeightedTracking(seedNode, endNodes, thisMaxTrackScore); 
         track = seedNode->getBestTrack();
         
         if (kDoTrackerPropagation) {
            // Match track to trackerTrack
            Cluster * firstClusterInTrack = track->Last(); // Last added, not yet sorted
            for (Int_t j=0; j<trackerTracks->GetEntriesFast(); j++) {
               Cluster *tc = trackerTracks->At(j)->Last();
               Bool_t xOK = firstClusterInTrack->getX() == tc->getX();
               Bool_t yOK = firstClusterInTrack->getY() == tc->getY();
               Bool_t zOK = firstClusterInTrack->getLayer() == tc->getLayer();
               if (xOK && yOK && zOK) {
                  track->appendCluster(trackerTracks->At(j)->At(0));
                  track->appendCluster(trackerTracks->At(j)->At(1));
               }
            }
         }

         track->sortTrack();

         if (track->At(0)->getLayer() < 2) {
            tracks->appendTrack(track);
            removeTrackFromClustersWithoutTrack(track);
            markUsedClusters(track);
         }

         delete track;
         seedNode->deleteNodeTree();
         delete seedNode;
         delete endNodes;
      }
       
      showDebug("Found " << tracks->GetEntries() << " tracks so far\n");
      delete seeds;
      delete nextClusters;
   }
   
   delete trackerTracks;
   return tracks;
}

void Clusters::doRecursiveWeightedTracking(Node * seedNode, vector<Node*> * endNodes, Float_t thisMaxTrackScore) {
   // Input: vector of nodes
   // Find all the next potential nodes segments w/acceptable score
   // Output: Vector of nodes (new end nodes)

   //cout<<"doRecursiveWeightedTracking from layer " << endNodes->at(0)->getCluster()->getLayer() << " with "  << seedNode->getNChildren() << " children" <<endl;
   Node    * thisNode;
   Cluster * nextCluster;
   Float_t   nextScore, nextAngle;
 
   for (UInt_t i=0; i<endNodes->size(); i++) { // All identified endpoints to the tree so far
      thisNode = endNodes->at(i);
      thisNode->markExplored();

      Int_t searchLayer = thisNode->getCluster()->getLayer() - 1;
      Int_t idxFrom = getFirstIndexOfLayer(searchLayer); 
      Int_t idxTo = getLastIndexOfLayer(searchLayer);

      if (kDoTrackerPropagation && searchLayer < 2) break;

      if (idxFrom < 0) continue; 

      for (Int_t j=idxFrom; j<idxTo; j++) { // Loop through possible additions to the track
         nextCluster = At(j);
         if (!nextCluster) continue;
        
         nextScore = thisNode->getNextScore(nextCluster);
         nextAngle = thisNode->getNodeAngle(nextCluster);

         float alwaysAllowAngle = 0.1;
         if (kHelium) alwaysAllowAngle = 0.03;

         if (nextScore < thisMaxTrackScore || nextAngle < alwaysAllowAngle){
            thisNode->addChild(new Node(thisNode, nextCluster, nextScore));
          }
      }
   }
   
   endNodes->clear();
   seedNode->getUnexploredEndNodes(endNodes);

   if (endNodes->size() > 0){
      doRecursiveWeightedTracking(seedNode, endNodes, thisMaxTrackScore);
   }
}

Tracks * Clusters::findCalorimeterTracksWithMCTruth() {
   // Now clusters is sorted by event ID
   // Loop through list, and fill to track when new event ID is encountered

   Cluster   * cluster = nullptr;
   Track     * track = new Track();
   Tracks    * tracks = new Tracks(kEventsPerRun * 5);
   Int_t       lastEventID, eventID;
   Bool_t      hasSecondary = false;
   Int_t       lastLayer;

   lastEventID = At(0)->getEventID();

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      cluster = At(i);
      eventID = cluster->getEventID();

      if (eventID != lastEventID) {
         if (track->GetEntries()) {
            tracks->appendTrack(track);
         }
         track->Clear("C");
         hasSecondary = false;
      }

      // IF SECONDARY LOOP THROUGH ALL, FIND SAME LAYER AND COMPARE PDG, 
      // CHOOSE HIGHEST PDG

      Int_t thisLayer = cluster->getLayer();
      Int_t thatLayer;
      for (Int_t j=0; j<track->GetEntriesFast(); j++) {
         thatLayer = track->getLayer(j);
         if (thisLayer == thatLayer) {
            if (track->At(j)->getPDG() < cluster->getPDG()) {
               // replace
               track->At(j)->set(cluster);
            }
         }
      }

      if (cluster->getPDG() > 1000) {
         track->appendCluster(cluster);
      }
            
      lastLayer = cluster->getLayer();

      if (cluster->isSecondary()) hasSecondary = true;

      lastEventID = eventID;
   }

   if (track->GetEntriesFast()) {
      tracks->appendTrack(track);
   }

   delete track;
   return tracks;
}

void Clusters::findSeeds(vector<Int_t> * seeds, Int_t layer, Bool_t kUsedClustersInSeeds) {
   Int_t layerIdxFrom = getFirstIndexOfLayer(layer);
   Int_t layerIdxTo = getLastIndexOfLayer(layer);

   if (layerIdxFrom>=0) { 
      for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
         if (!At(i)) continue;
         if (!kUsedClustersInSeeds && isUsed(i)) continue;
         seeds->push_back(i);
      }
   }
}

void Clusters::findNearestClustersInNextLayer(Cluster *seed, vector<Int_t> * nextClusters) {
   nextClusters->clear();
   for (Int_t skipLayers=0; skipLayers<2; skipLayers++) { // Loop to enable tracks skipping a layer without data
      Int_t nextLayer = seed->getLayer() - 1 - skipLayers;
      findClustersFromSeedInLayer(seed, nextLayer, nextClusters);
      if (nextClusters->size() > 0) break;
   }
}

void Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer, vector<Int_t> * nextClusters) { 
   Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
   Float_t maxAngle, thisAngle;

   maxAngle = 0.15 * kMCSFactor;
   Cluster * seedCorr = new Cluster();

   if (nextLayer == 1) {
      maxAngle = 0.5; // 500 mrad
      seedCorr->setXmm(seed->getXmm() - tan(getAngleAtSpot(kSpotX)*dz2)); 
      seedCorr->setYmm(seed->getYmm() - tan(getAngleAtSpot(kSpotY)*dz2)); 
      seedCorr->setLayer(-1);
   }
   else {
      seedCorr->setXmm(seed->getXmm()); 
      seedCorr->setYmm(seed->getYmm()); 
      seedCorr->setLayer(seed->getLayer());
   }

   if (layerIdxFrom >= 0) {
      for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
         if (!At(i)) { continue; }
         thisAngle = getDotProductAngle(seedCorr, seed, At(i));
          if ((thisAngle < maxAngle) && !isUsed(i)){
            nextClusters->push_back(i);
         }
      }
   }
}

Cluster * Clusters::findNearestNeighbour(Track *track, Cluster *projectedPoint, Bool_t rejectUsed) {
   // Finds the cluster in layer projectedPoint->getLayer() closest to the projectedPoint,
   // based on the distance from a projected vector from *track. (angle minimization)
   //
   // In some cases we don't know which track to connect the point to - in that case
   // Track * track = nullptr, and the legacy method (diffmmXY) is used

   Cluster *nearestNeighbour = nullptr;
   Float_t  thisAngle, maxAngle;
   Bool_t   kFoundNeighbour = false;
   Bool_t   reject = false;
   Int_t    searchLayer = projectedPoint->getLayer();
   Int_t    layerIdxFrom = getFirstIndexOfLayer(searchLayer);
   Int_t    layerIdxTo = getLastIndexOfLayer(searchLayer);

   maxAngle = 0.05 * kMCSFactor; 

   if (layerIdxFrom < 0) return 0;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
      if (!At(i))
         continue;

      thisAngle = getDotProductAngle(track->At(track->GetEntriesFast()-2), track->Last(), At(i));

      reject = (At(i)->isUsed() && rejectUsed);

      if (thisAngle < maxAngle && !reject) {
         nearestNeighbour = At(i);
         maxAngle = thisAngle;
         kFoundNeighbour = true;
      }
   }
   
   if (kFoundNeighbour)   return nearestNeighbour;
   else                   return nullptr;
}

#endif
