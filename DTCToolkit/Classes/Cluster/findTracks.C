#ifndef findTracks_cxx
#define findTracks_cxx

#include <vector>

#include <TClonesArray.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Node.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;
using namespace DTC;

Tracks * Clusters::findTracksWithRecursiveWeighting() {
   Cluster   * nextCluster = nullptr;
   Track     * track = nullptr;
   Tracks    * tracks = new Tracks(kEventsPerRun * 5);
   Float_t     clusterScore;
   Node      * nextNode = nullptr;
   Node      * seedNode = nullptr;
   vector<Int_t> * nextClusters = new vector<Int_t>;
   vector<Int_t> * seeds = new vector<Int_t>;
  
   Int_t       spotSize = 44;
   // for pencil beams.
   // 22 = 2x2 mm, 33 = 3x3 mm etc.
   // See WoC publication for description of functions below

   if      (spotSize == 22) kMaxTrackScore = 0.433 * pow(kEventsPerRun, -0.172);
   else if (spotSize == 33) kMaxTrackScore = 0.469 * pow(kEventsPerRun, -0.176);
   else if (spotSize == 42) kMaxTrackScore = 0.397 * pow(kEventsPerRun, -0.150);
   else if (spotSize == 55) kMaxTrackScore = 0.455 * pow(kEventsPerRun, -0.162);
   else if (spotSize == 44) kMaxTrackScore = 0.460 * pow(kEventsPerRun, -0.168);

   kMaxTrackScore = 0.445 * pow(kEventsPerRun, -0.176); // With 4 mm geometry, 5% less (!) scattering

   // kMaxTrackScore = 0.746 * pow(kEventsPerRun, -0.242); // No Inelastic scattering
   // kMaxTrackScore = 0.289 * pow(kEventsPerRun, -0.113); // No Elastic Scattering
   // kMaxTrackScore = 0.382 * pow(kEventsPerRun, -0.151); // No nuclear scattering
   // kMaxTrackScore = 0.200 * pow(kEventsPerRun, -0.05); // no scattering
   
   // A constant allowance of 50 mrad per layer is always allowed
   // If that's unwanted uncomment the line below (or see Constants.h)
   // kMaxTrackAngle = 0;
   
   nextClusters->reserve(50);
   
   showDebug("makeLayerIndex...");
   makeLayerIndex();
   showDebug("ok!\nfillMSCradiusList...");
   fillMCSRadiusList();
   showDebug("ok!\n");
   kMCSFactor = 25;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      appendClusterWithoutTrack(At(i));
   }

   findSeeds(seeds, 0, true);

   for (UInt_t i=0; i<seeds->size(); i++) {
      Cluster *seed = At(seeds->at(i));

      nextClusters->clear();
      findNearestClustersInNextLayer(seed, nextClusters);
      seedNode = new Node(nullptr, seed, 0);
      
      for (UInt_t j=0; j<nextClusters->size(); j++) {
         nextCluster = At(nextClusters->at(j));
         
         clusterScore = seedNode->getNextScore(nextCluster);
         clusterScore /= 100;
         if (std::isnan(clusterScore)) cout << "clusterScore at " << *nextCluster << " isNan!\n";

         if (clusterScore < kMaxTrackScore) {
            nextNode = new Node(seedNode, nextCluster, clusterScore); // initial vector is 10 % weighted
            seedNode->addChild(nextNode);
         }
      }

      seedNode->markExplored();

      vector<Node*> * endNodes = new vector<Node*>;
      endNodes->reserve(kEventsPerRun * 5);
      seedNode->getUnexploredEndNodes(endNodes);

      showDebug("doRecursiveWeightedTracking");
      doRecursiveWeightedTracking(seedNode, endNodes);
      showDebug("..ok!\ngetBestTrack...");
      track = seedNode->getBestTrack();
      showDebug("ok!\n");
      
      if (track->GetEntries() >= 3) {
         tracks->appendTrack(track);
         removeTrackFromClustersWithoutTrack(track);
      }

      delete track;
      
      seedNode->deleteNodeTree();
      delete seedNode;
      delete endNodes;
   }

   delete seeds;
   delete nextClusters;

   return tracks;
}

void Clusters::doRecursiveWeightedTracking(Node * seedNode, vector<Node*> * endNodes) {
   // Input: vector of nodes
   // Find all the next potential nodes segments w/acceptable score
   // Output: Vector of nodes (new end nodes)

   Node    * thisNode;
   Cluster * nextCluster;
   Float_t   nextScore, nextAngle;
   
   for (UInt_t i=0; i<endNodes->size(); i++) { // All identified endpoints to the tree so far
      thisNode = endNodes->at(i);
      thisNode->markExplored();

      // Optimization
      Int_t searchLayer = thisNode->getCluster()->getLayer() + 1;
      Int_t idxFrom = getFirstIndexOfLayer(searchLayer);
      Int_t idxTo = getLastIndexOfLayer(searchLayer);
      if (idxFrom < 0) continue; // no more clusters in deeper layers

      for (Int_t j=idxFrom; j<idxTo; j++) { // Loop through possible additions to the track
         nextCluster = At(j);
         if (!nextCluster) continue;
        
         nextScore = thisNode->getNextScore(nextCluster);
         nextAngle = thisNode->getNodeAngle(nextCluster);
         
         if (nextScore < kMaxTrackScore || nextAngle < kMaxTrackAngle) {
            thisNode->addChild(new Node(thisNode, nextCluster, nextScore)); // it is either appended to the tree or deleted if thisNode is full
         }
      }
   }

   endNodes->clear();
   seedNode->getUnexploredEndNodes(endNodes);

   if (endNodes->size() > 0) {
      doRecursiveWeightedTracking(seedNode, endNodes);
   }
}

Tracks * Clusters::findCalorimeterTracksWithMCTruth() {
   // Now clusters is sorted by event ID
   // Loop through list, and fill to track when new event ID is encountered

   Cluster   * cluster = nullptr;
   Track     * track = new Track();
   Tracks    * tracks = new Tracks(kEventsPerRun * 5);
   Int_t       lastEventID, eventID;

   lastEventID = At(0)->getEventID();

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) {
         cout << "Clusters::findCalorimeterTracksWithMCTruth could not find cluster at " << i << endl;
         continue;
      }

      cluster = At(i);
      eventID = cluster->getEventID();

      if (eventID != lastEventID) {
         // Cluster from new event ID encountered, store what we have already
         showDebug("Track at  is full, storing the pointer to tracks having currently " << tracks->GetEntriesFast() << " elements.\n");
         tracks->appendTrack(track);
         track->Clear("C");
      }

      showDebug("Appending cluster at " << *cluster << " with  number " << track->GetEntriesFast() + 1 << " to track.\n");
      track->appendCluster(cluster);

      lastEventID = eventID;
   }

   // Store last track
   tracks->appendTrack(track);

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
   for (Int_t skipLayers=0; skipLayers<2; skipLayers++) { // Loop to enable tracks skipping a layer without data
      Int_t nextLayer = seed->getLayer() + 1 + skipLayers;
      findClustersFromSeedInLayer(seed, nextLayer, nextClusters);
      if (nextClusters->size() > 0) break;
   }
}

void Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer, vector<Int_t> * nextClusters) {
   Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
   Float_t maxAngle, thisAngle;

   if (kUseEmpiricalMCS) {
      maxAngle = 3 * getEmpiricalMCSAngle(nextLayer - 1); // This function is only run once, to find seed candidates in the 2nd layer - here it's okay with a high angle
      if (layerIdxFrom >= 0) {
         for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
            if (!At(i)) { continue; }
            thisAngle = getDotProductAngle(seed, seed, At(i));
            if (thisAngle < maxAngle) {
               nextClusters->push_back(i);
            }
         }
      }
   }

   else { // Use "old" MCS estimation method
      maxAngle = 3 * getSearchRadiusForLayer(nextLayer) * 0.75 * MCSMultiplicationFactor;

      if (layerIdxFrom >= 0) {
         for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
            if (!At(i)) { continue; }
            thisAngle = diffmmXY(seed, At(i));
            if (thisAngle < maxAngle) {
               nextClusters->push_back(i);
            }
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

   if (kUseEmpiricalMCS && track)   maxAngle = getEmpiricalMCSAngle(searchLayer-1);
   else                             maxAngle = getSearchRadiusForLayer(searchLayer) * MCSMultiplicationFactor;


   if (layerIdxFrom < 0) return 0;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
      if (!At(i))
         continue;

      if (kUseEmpiricalMCS && track)   thisAngle = getDotProductAngle(track->At(track->GetEntriesFast()-2), track->Last(), At(i));
      else                             thisAngle = diffmmXY(projectedPoint, At(i));

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
