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
   Cluster   * cluster = nullptr;
   Cluster   * nextCluster = nullptr;
   Int_t       nFirstCandidates = 0;
   Track     * track = nullptr;
   Tracks    * tracks = new Tracks(kEventsPerRun * 5);
   Float_t     clusterScore;
   Node      * nextNode = nullptr;
   Node      * seedNode = nullptr;
  
   ///// THESE VALUES ARE VALID FOR PENCIL BEAMS
   // 2x2 -> 1.49
   // 3x3 -> 2.01
   // 4x4 -> 3.20
   // 5x5 -> 4.22
   // 4x2 -> 1.88

   Int_t spotSize = 33;
   // 22 = 2x2 mm, 33 = 3x3 mm etc.

   if      (spotSize == 22) kMaxTrackScore = 0.433 * pow(kEventsPerRun, -0.172);
   else if (spotSize == 33) kMaxTrackScore = 0.469 * pow(kEventsPerRun, -0.176);
   else if (spotSize == 42) kMaxTrackScore = 0.397 * pow(kEventsPerRun, -0.150);
   else if (spotSize == 55) kMaxTrackScore = 0.455 * pow(kEventsPerRun, -0.162);
   else if (spotSize == 44) kMaxTrackScore = 0.460 * pow(kEventsPerRun, -0.168);

   kMaxTrackScore = 0.445 * pow(kEventsPerRun, -0.176); // With 4 mm geometry, 5% less (!) scattering

//   kMaxTrackScore = 0.3;

   // kMaxTrackScore = 0.746 * pow(kEventsPerRun, -0.242); // No Inelastic scattering
   // kMaxTrackScore = 0.289 * pow(kEventsPerRun, -0.113); // No Elastic Scattering
   // kMaxTrackScore = 0.382 * pow(kEventsPerRun, -0.151); // No nuclear scattering
   // kMaxTrackScore = 0.200 * pow(kEventsPerRun, -0.05); // no scattering

   // UNIFORM FIELD
   // From 0.19 to 0.3 @ 35 per cm2
   /*
   trackDensity = kEventsPerRun / 100;
   if (trackDensity > 140)  kMaxTrackScore = 0.19;
   else                     kMaxTrackScore = 0.30;
   */
   
   showDebug("makeLayerIndex..\n");
   makeLayerIndex();
   showDebug("ok!\nfillMSCradiusList...");
   fillMCSRadiusList();
   showDebug("ok!\n");
   kMCSFactor = 25;
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      appendClusterWithoutTrack(At(i));
   }

   Clusters  * seeds = findSeeds(0, true);
  
   for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
      nFirstCandidates = 0;
      Cluster *seed = seeds->At(i);
      Clusters * nextClusters = findNearestClustersInNextLayer(seed);
      
      seedNode = new Node(nullptr, seed, 0);
      for (Int_t j=0; j<nextClusters->GetEntriesFast(); j++) {
         nextCluster = nextClusters->At(j);

         clusterScore = seedNode->getNextScore(nextCluster);
         clusterScore /= 100;
         if (std::isnan(clusterScore)) {
            cout << "clusterScore for seed at " << *nextCluster << " isNan!\n";
         }

         if (clusterScore < kMaxTrackScore) {
            nextNode = new Node(seedNode, nextCluster, clusterScore); // initial vector is 10 % weighted
            seedNode->addChild(nextNode);
            nFirstCandidates++;
         }
      }
      
      seedNode->markExplored();

      vector<Node*> * endNodes = new vector<Node*>;
      endNodes->reserve(kEventsPerRun * 5);
      seedNode->getUnexploredEndNodes(endNodes);

      showDebug("Performing recursive tracking...");
      doRecursiveWeightedTracking(seedNode, endNodes);
      showDebug("ok\n");
      
      endNodes->clear();
      seedNode->getEndNodes(endNodes);

      showDebug("Finding best track...");
      track = seedNode->getBestTrack();
      showDebug("ok\n");
      if (track->GetEntries() > 3) {
         tracks->appendTrack(track);
         removeTrackFromClustersWithoutTrack(track);
      }

      delete seedNode;
      delete endNodes;
   }

//   kMCSFactor = kMCSFactorLastPass1;
//   findRemainingTracks(tracks);

   return tracks;
}

void Clusters::doRecursiveWeightedTracking(Node * seedNode, vector<Node*> * endNodes) {
   // Input: vector of nodes
   // Find all the next potential nodes segments w/acceptable score
   // Output: Vector of nodes (new end nodes)

   Node    * thisNode;
   Node    * nextNode;
   Cluster * nextCluster;
   Float_t   nextScore, nextAngle;
   Int_t     nPotentials;
   Float_t   bestScore = 1e5;

   for (UInt_t i=0; i<endNodes->size(); i++) {
      thisNode = endNodes->at(i);
      thisNode->markExplored();
      nPotentials = 0;
      bestScore = 1e5;

      // search all clusters for next potential node
      Int_t searchLayer = thisNode->getCluster()->getLayer() + 1;
      Int_t idxFrom = getFirstIndexOfLayer(searchLayer);
      Int_t idxTo = getLastIndexOfLayer(searchLayer);
  
      if (idxFrom < 0) continue; // no more clusters

      for (Int_t j=idxFrom; j<idxTo; j++) {
         // Calculate scores of all potential nodes
         nextCluster = At(j);
         if (!nextCluster) continue;
        
         nextScore = thisNode->getNextScore(nextCluster);
         nextAngle = thisNode->getNodeAngle(nextCluster);
         
         bestScore = fmin(bestScore, nextScore);
         if (nextScore < kMaxTrackScore || nextAngle < kMaxTrackAngle) {
            nPotentials++;
            nextNode = new Node(thisNode, nextCluster, nextScore);
            thisNode->addChild(nextNode);
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
         track->Clear();
      }

      showDebug("Appending cluster at " << *cluster << " with  number " << track->GetEntriesFast() + 1 << " to track.\n");
      track->appendCluster(cluster);

      lastEventID = eventID;
   }

   // Store last track
   tracks->appendTrack(track);

   return tracks;
}

Clusters * Clusters::findSeeds(Int_t layer, Bool_t kUsedClustersInSeeds) {
   Clusters *seeds = new Clusters(1000);
   
   Int_t layerIdxFrom = getFirstIndexOfLayer(layer);
   Int_t layerIdxTo = getLastIndexOfLayer(layer);

   showDebug("Layer Index of layer " << layer << " from " << layerIdxFrom << ", to " << layerIdxTo << ")\n");
   
   if (layerIdxFrom<0)
      return seeds;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
      if (!At(i)) { continue; }
      if (!kUsedClustersInSeeds && isUsed(i)) { continue; }
      
      seeds->appendCluster(At(i));
   }

   return seeds;
}

Clusters * Clusters::findNearestClustersInNextLayer(Cluster *seed) {
   Clusters *nextClusters = new Clusters(50);
   Clusters *clustersFromThisLayer = 0;

   Int_t layerCounter = 1;
   for (Int_t skipLayers=0; skipLayers<2; skipLayers++) {
      Int_t nextLayer = seed->getLayer() + 1 + skipLayers;
      clustersFromThisLayer = findClustersFromSeedInLayer(seed, nextLayer);

      if (clustersFromThisLayer->GetEntriesFast()) {
         break;
      }
   }

   Int_t nClusters = clustersFromThisLayer->GetEntriesFast();
   for (Int_t i=0; i<nClusters; i++) {
      if (!clustersFromThisLayer->At(i)) { continue; }
      
      nextClusters->appendCluster(clustersFromThisLayer->At(i));
   }

   delete clustersFromThisLayer;

   return nextClusters;
}

Clusters * Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer) {
   Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
   Clusters *clustersFromThisLayer = new Clusters(50);

   Float_t maxAngle, thisAngle;

   if (kUseEmpiricalMCS) maxAngle = getEmpiricalMCSAngle(nextLayer - 1);
   else                  maxAngle = getSearchRadiusForLayer(nextLayer) * 0.75 * MCSMultiplicationFactor;

   maxAngle *= 3;

   if (layerIdxFrom < 0)
      return clustersFromThisLayer; // empty

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
      if (!At(i)) { continue; }

      if (kUseEmpiricalMCS)   thisAngle = getDotProductAngle(seed, seed, At(i));
      else                    thisAngle = diffmmXY(seed, At(i));

      if (thisAngle < maxAngle) {
         clustersFromThisLayer->appendCluster(At(i));
      }
   }

   return clustersFromThisLayer;
}

Cluster * Clusters::findNearestNeighbour(Track *track, Cluster *projectedPoint, Bool_t rejectUsed) {
   // Finds the cluster in layer projectedPoint->getLayer() closest to the projectedPoint,
   // based on the distance from a projected vector from *track. (angle minimization)
   //
   // In some cases we don't know which track to connect the point to - in that case
   // Track * track = nullptr, and the legacy method (diffmmXY) is used
   Cluster *nearestNeighbour = new Cluster();
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
         nearestNeighbour->set(At(i));
         maxAngle = thisAngle;
         kFoundNeighbour = true;
      }
   }
   
   if (kFoundNeighbour)   return nearestNeighbour;
   else                   return nullptr;
}

#endif
