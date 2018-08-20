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
   
//   kMaxTrackScore = 0.4;

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

Tracks * Clusters::findCalorimeterTracksAlpide() {
   Tracks * tracks = new Tracks(kEventsPerRun * 5);
   Int_t    startOffset = 0;
   Bool_t   usedClustersInSeeds = true;
   Int_t    clustersLeft;
   Float_t  factor;
   
   TStopwatch t0, t05, t1, t2, t3, t4, t5;


   t0.Start();
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      appendClusterWithoutTrack(At(i));
   }
   t0.Stop();
   
   t05.Start();
   makeLayerIndex();
   t05.Stop();

   // first pass, small search cone (3 sigma MCS)
   kMCSFactor = kMCSFactorFirstPass;
   t1.Start();
   findTracksFromLayer(tracks, 0, usedClustersInSeeds);
   t1.Stop();
  
   t2.Start();
   if (GetEntriesCWT()) { 
      kMCSFactor = kMCSFactorSecondPass;
      findTracksFromLayer(tracks, 0, usedClustersInSeeds);
   }
   t2.Stop();

   t3.Start();
   if (GetEntriesCWT()) {
      kMCSFactor = kMCSFactorLastPass1;
      findRemainingTracks(tracks);
   }
   t3.Stop();

   t4.Start();
   if (GetEntriesCWT()) {
      kMCSFactor = kMCSFactorLastPass2;
    findRemainingTracks(tracks);
   }
   t4.Stop();
   
   t5.Start();
   if (GetEntriesCWT()) {
      kMCSFactor = kMCSFactorLastPass3;
      findRemainingTracks(tracks);
   }
   t5.Stop();

   printf("Reconstruction TIMING: t0 = %.3f, t05 = %.3f, t1 = %.3fs, t2 = %.3fs, t3 = %.3fs, t4 = %.3fs, t5 = %.3fs.\n", t0.CpuTime(), t05.CpuTime(), t1.CpuTime(), t2.CpuTime(), t3.CpuTime(), t4.CpuTime(), t5.CpuTime());

   clustersLeft = clustersWithoutTrack_.GetEntries();
   factor = 100 * (1 - (Float_t) clustersLeft / GetEntriesFast());
   cout << "Found " << tracks->GetEntriesFast() << " tracks. " << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track! (" << factor << " %)\n";

   return tracks;
}


Tracks * Clusters::findCalorimeterTracks() {
   Tracks * tracks = new Tracks(kEventsPerRun * 5);
   Int_t    startOffset = 0;
   Bool_t   usedClustersInSeeds = true;
   Int_t    clustersLeft;
   Float_t  factor;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      appendClusterWithoutTrack(At(i));
   }
   
   makeLayerIndex();
   fillMCSRadiusList();
   MCSMultiplicationFactor = 4; // 3
   kMCSFactor = kMCSFactorFirstPass * 2.5;

   // first pass, small search cone (3 sigma MCS)
   findTracksFromLayer(tracks, 0, usedClustersInSeeds);

   MCSMultiplicationFactor = 8; // 5
   usedClustersInSeeds = false;
//   multiplyRadiusFirstLayers(2);

   if (clustersWithoutTrack_.GetEntries() > 0) {
      // second pass, large seach cone (5 sigma MCS)
      findTracksFromLayer(tracks, 0, usedClustersInSeeds);
   
      // third pass, start from layer 1
      findTracksFromLayer(tracks, 1, usedClustersInSeeds);
      findTracksFromLayer(tracks, 2, usedClustersInSeeds);
   }

   clustersLeft = clustersWithoutTrack_.GetEntries();
   factor = 100 * (1 - (Float_t) clustersLeft / GetEntriesFast());
   cout << "Found " << tracks->GetEntriesFast() << " tracks. " << clustersLeft << " of total " << GetEntriesFast() << " clusters were not assigned to track! (" << factor << " %)\n";

   return tracks;
}

void Clusters::findTracksFromLayer(Tracks * tracks, Int_t layer, Bool_t kUsedClustersInSeeds) {
   Int_t       startOffset = 0;
   Track     * bestTrack = nullptr;
   Track     * newBestTrack = nullptr;
   Clusters  * seeds = nullptr;
   Int_t       minimumLengthOfTrackToPass = 4; // first pass

   if (!kUsedClustersInSeeds) {
      minimumLengthOfTrackToPass = 1; // second / third pass
   }

   seeds = findSeeds(layer, kUsedClustersInSeeds);
   showDebug("Found " << seeds->GetEntriesFast() << " seeds in layer " << layer << endl);
   
   for (Int_t i=0; i<seeds->GetEntriesFast(); i++) {
      if (!seeds->At(i))
         continue;

      showDebug("Seed(At(i)) = " << *seeds->At(i) << endl);

      bestTrack = trackPropagation(seeds->At(i));
      showDebug("Best track has " << bestTrack->GetEntriesFast() << " entries.\n");

      if (bestTrack->GetEntriesFast() < minimumLengthOfTrackToPass)
         continue;

      if (!bestTrack->At(0))
         startOffset = 1;
      
      tracks->appendTrack(bestTrack, startOffset);
      removeTrackFromClustersWithoutTrack(bestTrack);
      markUsedClusters(bestTrack);
   }

   Int_t nTracksWithConflicts = 0;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      if (!tracks->At(i)) continue;
      nTracksWithConflicts += (int) tracks->isUsedClustersInTrack(i);
   }

   delete bestTrack;
   delete seeds;
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

Track * Clusters::trackPropagation(Cluster *seed) {
   Tracks    * seedTracks = new Tracks(100);
   Track     * currentTrack = new Track();
   Track     * longestTrack = new Track();
   Cluster   * nextCluster = nullptr;
   Int_t       fromLayer;

   Clusters * nextClusters = findNearestClustersInNextLayer(seed);

   for (Int_t i=0; i<nextClusters->GetEntriesFast(); i++) {
      nextCluster = nextClusters->At(i);
      fromLayer = nextCluster->getLayer();

      currentTrack->appendCluster(seed);
      currentTrack->appendCluster(nextCluster);

      growTrackFromLayer(currentTrack, fromLayer);

      if (currentTrack->GetEntriesFast()){
         seedTracks->appendTrack(currentTrack);
      }

      currentTrack->clearTrack();
   }
   
   if (seedTracks->GetEntriesFast()) {
      longestTrack = findLongestTrack(seedTracks);
   }

   delete seedTracks;
   delete currentTrack;
   delete nextClusters;

   return longestTrack;
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
