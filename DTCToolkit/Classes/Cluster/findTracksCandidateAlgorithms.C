#ifndef findTracksCandidateAlgorithms_cxx
#define findTracksCandidateAlgorithms_cxx

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

   Tracks * Clusters::findTracksWithReverseRecursiveWeighting() {
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

      Float_t trackDensity = kEventsPerRun / 1.88;
      kMaxTrackScore = 0.571 * pow(trackDensity, -0.238); // See "Ulike interaksjoner.xlsx"

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

      Int_t lastLayer = getLastActiveLayer();

      Clusters  * seeds = findSeeds(lastLayer, true);
      printf("Found %d seeds in first pass\n", seeds->GetEntriesFast());
     
      for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
         nFirstCandidates = 0;
         Cluster *seed = seeds->At(i);
         Clusters * nextClusters = findNearestClustersInLastLayer(seed);
         printf("Found %d nextClusters in first pass\n", nextClusters->GetEntriesFast());
         
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
         printf("ncandidates first pass = %d.\n", nFirstCandidates);
         
         seedNode->markExplored();
      
         vector<Node*> * endNodes = new vector<Node*>;
         endNodes->reserve(kEventsPerRun * 5);
         seedNode->getUnexploredEndNodes(endNodes);

         showDebug("Performing recursive tracking...");
         doReverseRecursiveWeightedTracking(seedNode, endNodes);
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

      seeds = findSeeds(lastLayer-1, false);
      printf("Found %d seeds in second pass\n", seeds->GetEntriesFast());
     
      for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
         nFirstCandidates = 0;
         Cluster *seed = seeds->At(i);
         Clusters * nextClusters = findNearestClustersInLastLayer(seed);
         
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
         doReverseRecursiveWeightedTracking(seedNode, endNodes);
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
      
      seeds = findSeeds(lastLayer-2, false);
      printf("Found %d seeds in third pass\n", seeds->GetEntriesFast());
     
      for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
         nFirstCandidates = 0;
         Cluster *seed = seeds->At(i);
         Clusters * nextClusters = findNearestClustersInLastLayer(seed);
         
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
         doReverseRecursiveWeightedTracking(seedNode, endNodes);
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

      seeds = findSeeds(lastLayer-3, false);
      printf("Found %d seeds in fourth pass\n", seeds->GetEntriesFast());
     
      for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
         nFirstCandidates = 0;
         Cluster *seed = seeds->At(i);
         Clusters * nextClusters = findNearestClustersInLastLayer(seed);
         
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
         doReverseRecursiveWeightedTracking(seedNode, endNodes);
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

   void Clusters::doReverseRecursiveWeightedTracking(Node * seedNode, vector<Node*> * endNodes) {
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
         Int_t searchLayer = thisNode->getCluster()->getLayer() - 1;
         if (searchLayer < 0) break;

         Int_t idxFrom = getFirstIndexOfLayer(searchLayer);
         Int_t idxTo = getLastIndexOfLayer(searchLayer);
    
         printf("searchlayer: %d, idxfrom = %d, idxto = %d\n", searchLayer, idxFrom, idxTo);

         if (idxFrom < 0) continue; // no more clusters

         for (Int_t j=idxFrom; j<idxTo; j++) {
            // Calculate scores of all potential nodes
            nextCluster = At(j);
            if (!nextCluster) continue;
           
            nextScore = thisNode->getNextScore(nextCluster);
            nextAngle = thisNode->getNodeAngle(nextCluster);
            
            cout << "nextScore = " << nextScore << ", nextAngle = " << nextAngle << endl;

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
         doReverseRecursiveWeightedTracking(seedNode, endNodes);
      }
   }

Clusters * Clusters::findNearestClustersInLastLayer(Cluster *seed) {
   Clusters *nextClusters = new Clusters(50);
   Clusters *clustersFromThisLayer = 0;

   Int_t layerCounter = 1;
   for (Int_t skipLayers=0; skipLayers<2; skipLayers++) {
      Int_t nextLayer = seed->getLayer() - 1 - skipLayers;
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

void Clusters::growTrackFromLayer(Track *track, Int_t fromLayer) {
   Cluster   * projectedPoint = nullptr;
   Cluster   * nearestNeighbour = nullptr;
   Int_t       nSearchLayers = getLastActiveLayer();
   Int_t       lastHitLayer = fromLayer;

   for (Int_t layer = lastHitLayer + 1; layer <= nSearchLayers+1; layer++) {
      projectedPoint = getTrackExtrapolationToLayer(track, layer);
      if (isPointOutOfBounds(projectedPoint, 10)) break;

      nearestNeighbour = findNearestNeighbour(track, projectedPoint);

      if (nearestNeighbour) {
         track->appendCluster(nearestNeighbour);
         lastHitLayer = layer;
      }

      else { // No matching points, trying next layer
         projectedPoint = getTrackExtrapolationToLayer(track, layer+1);
         if (isPointOutOfBounds(projectedPoint, 10)) break;

         nearestNeighbour = findNearestNeighbour(track, projectedPoint);

         if (nearestNeighbour) { // Found one there! Don't search next layer
            track->appendCluster(nearestNeighbour);
            lastHitLayer = ++layer+1;
         }
      }

      if (layer>lastHitLayer+2) {
         break;
      }
   }

   delete projectedPoint;
   delete nearestNeighbour;
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

Track * Clusters::findLongestTrack(Tracks *seedTracks) {
   if (!seedTracks->GetEntriesFast())
      return new Track();
   
   Float_t  bestScore = -1;
   Float_t  score = -1;
   Track  * longestTrack = new Track();
   Track  * track = nullptr;
   Int_t    startOffset;

   for (Int_t i=0; i<seedTracks->GetEntriesFast(); i++) {
      if (!seedTracks->At(i)) continue;

      score = seedTracks->getTrackScore(i);
      showDebug("Track score is " << score << endl);
      if (score > bestScore) {
         bestScore = score;
         track = seedTracks->At(i);
      }
   }

   if (!track) return longestTrack;

   if (track->At(0)) {
      startOffset = (track->getLayer(0) == 0) ? 0 : 1;
   }
   
   else {
      startOffset = 1;
   }
   
   longestTrack->setTrack(track, startOffset);
   
   return longestTrack;
}

void Clusters::findRemainingTracks(Tracks * tracks) {
   Track *thisTrack = nullptr;

   showDebug("Clusters::findRemainingTracks: There are " << clustersWithoutTrack_.GetEntries() << " unused clusters left.\n");

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      if (!tracks->At(i)) continue;
      thisTrack = tracks->At(i);

      growTrackFromLayer(thisTrack, thisTrack->Last()->getLayer());

      removeTrackFromClustersWithoutTrack(thisTrack);
      markUsedClusters(thisTrack);
   }
   
   showDebug("Clusters::findRemainingTracks: After cleanup, there are " << clustersWithoutTrack_.GetEntries() << " unused clusters left.\n");
}

#endif
