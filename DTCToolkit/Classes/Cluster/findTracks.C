#include <vector>

#include <TClonesArray.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"


using namespace std;
using namespace DTC;

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
   MCSMultiplicationFactor = 3;

   // first pass, small search cone (3 sigma MCS)
   findTracksFromLayer(tracks, 0, usedClustersInSeeds);

   MCSMultiplicationFactor = 5;
   usedClustersInSeeds = false;
   multiplyRadiusFirstLayers(2);

   if (clustersWithoutTrack_.GetEntries() > 0) {
      // second pass, large seach cone (5 sigma MCS)
      findTracksFromLayer(tracks, 0, usedClustersInSeeds);
   
      // third pass, start from layer 1
      findTracksFromLayer(tracks, 1, usedClustersInSeeds);
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

      showDebug("Growing the track: " << *currentTrack << endl);
      growTrackFromLayer(currentTrack, fromLayer);
      showDebug("Result: " << *currentTrack << endl);

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

Clusters * Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer) {
   Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
   Clusters *clustersFromThisLayer = new Clusters(50);

   Float_t maxAngle, thisAngle;

   if (kUseEmpiricalMCS) maxAngle = getEmpiricalMCSAngle(nextLayer - 1);
   else                  maxAngle = getSearchRadiusForLayer(nextLayer) * 0.75 * MCSMultiplicationFactor;

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
   // Rewrite, we don't have to check if the neighbor is valid both here and in findNearestNeighbour... 
   // So removing a lot of testing.

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

   showDebug("Clusters::findNearestNeighbour using track " << *track << " and with projected point " << *projectedPoint << endl);

   if (kUseEmpiricalMCS && track)   maxAngle = getEmpiricalMCSAngle(searchLayer-1);
   else                             maxAngle = getSearchRadiusForLayer(searchLayer) * MCSMultiplicationFactor;

   showDebug("Clusters::findNearestNeighbour The max angle is " << maxAngle);

   if (layerIdxFrom < 0) return 0;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
      if (!At(i))
         continue;

      if (kUseEmpiricalMCS && track)   thisAngle = getDotProductAngle(track->At(track->GetEntriesFast()-2), track->Last(), At(i));
      else                             thisAngle = diffmmXY(projectedPoint, At(i));

      reject = (At(i)->isUsed() && rejectUsed);

      showDebug(" and found thisAngle " << thisAngle << " from point " << *At(i));

      if (thisAngle < maxAngle && !reject) {
         showDebug(", which is a current best match ... ");
         nearestNeighbour->set(At(i));
         maxAngle = thisAngle;
         kFoundNeighbour = true;
      }
   }
   
   showDebug("OK! Keeping last point.\n");

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
