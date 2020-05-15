#ifndef tracksOptimization_cxx
#define tracksOptimization_cxx

#include <iostream>
#include <vector>

#include <TClonesArray.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TView.h>
#include <TAttMarker.h>
#include <TAttLine.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>


#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

void Tracks::splitSharedClusters() {
   // to be run as early as possible before joining Tracks* objects

   Float_t     dist, minDist = 1e5;
   Float_t     x, y, clusterRadius, size, totalEdep, newSize;
   Int_t       minIdx = 0, idx = 0, clusterIdx;
   Int_t       nSplit = 0, nMissing = 0, nInterpolated = 0;
   Int_t       trackIdx;
   Bool_t      isMissing;
   Cluster   * interpolatedCluster = nullptr;
   Cluster   * interpolatedClosestCluster = nullptr;
   Cluster   * closestCluster = nullptr;
   Track     * missingTrack = nullptr;
   Track     * closestTrack = nullptr;
   Track     * thisTrack = nullptr;
//   Clusters  * interpolatedClusters = new Clusters();
//   Clusters  * interpolatedClosestClusters = new Clusters();

   

   for (Int_t layer=1; layer<nLayers; layer++) {
      
      vector<trackCluster> clustersThisLayer;
      vector<trackCluster> missingClustersThisLayer;
      
      clustersThisLayer.reserve(kEventsPerRun * 20);
      missingClustersThisLayer.reserve(kEventsPerRun * 5);

      for (Int_t i=0; i<GetEntriesFast(); i++) {
         thisTrack = At(i);

         if (!thisTrack) continue;
         if (thisTrack->getLastLayer() < layer) continue;

         trackCluster thisCluster;
         thisCluster.track = i;

         if (thisTrack->hasLayer(layer)) {
            thisCluster.cluster = thisTrack->getClusterFromLayer(layer);
            clustersThisLayer.push_back(thisCluster);
         }

         else {
            missingClustersThisLayer.push_back(thisCluster);
            nMissing++;
         }
      }
      
      for (UInt_t i=0; i<missingClustersThisLayer.size(); i++) {
         trackIdx = missingClustersThisLayer.at(i).track;
         missingTrack = At(trackIdx);
         if (!missingTrack) continue;
      
         interpolatedCluster = missingTrack->getInterpolatedClusterAt(layer);
         if (interpolatedCluster) {
            minIdx = getClosestCluster(clustersThisLayer, interpolatedCluster);
            if (minIdx<0) continue; // no clusters this layer? What else could provoke this error?
            nInterpolated++;

//            interpolatedClusters->appendCluster(interpolatedCluster);
            trackCluster closestTC = clustersThisLayer.at(minIdx); // THE CULPRIT
            closestTrack = At(closestTC.track);
            closestCluster = closestTrack->At(closestTC.cluster);
            
            minDist = diffmmXY(closestCluster, interpolatedCluster);
            clusterRadius = closestCluster->getRadiusmm();
            if (minDist < clusterRadius * 2) {
               
               interpolatedClosestCluster = closestTrack->getInterpolatedClusterAt(layer);
               if (interpolatedClosestCluster) {
//                  interpolatedClosestClusters->appendCluster(interpolatedClosestCluster);
                  
                  // Charge conservation: Each cluster gets half of the 'charge' ~ deposited energy
                  totalEdep = closestCluster->getDepositedEnergy();
                  newSize = getCSFromEdep(totalEdep / 2);
                  x = ( closestCluster->getX() + interpolatedCluster->getX() ) / 2;
                  y = ( closestCluster->getY() + interpolatedCluster->getY() ) / 2;
                  cout << "Moved cluster " << closestCluster->getXmm() - (closestCluster->getXmm() + interpolatedCluster->getXmm())/2. << " mm.\n";
                  missingTrack->appendCluster(new Cluster(x, y, layer, newSize, -1));
                  nSplit++;
                  
                  // set size of closest cluster as well!
                  closestCluster->setSize(newSize);
                  closestCluster->setEventID(-1);

                  sortTrackByLayer(trackIdx);
               }
               delete interpolatedClosestCluster;
            }
         }
         delete interpolatedCluster;
      }
   }

   if (nMissing>0) {
      cout << "From " << nMissing << " missing clusters, " << nInterpolated << " interpolated clusters were tested and " << nSplit << " clusters were split.\n";
   }
}

void Tracks::removeTracksLeavingDetector() {
   Int_t       nTracksRemoved = 0;
   Int_t       lastLayer;
   Cluster   * nextPoint = nullptr;
   Cluster   * lastCluster = nullptr;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      lastCluster = At(i)->Last();

      if (!lastCluster) continue;

      lastLayer = lastCluster->getLayer();
      if (lastLayer <= 1) continue;

      showDebug("Tracks::removeTracksLeavingDetector(): lastLayer = " << lastLayer << endl);
      nextPoint = getTrackExtrapolationToLayer(At(i), lastLayer + 1);
      if (!nextPoint) continue;
      
      if (isPointOutOfBounds(nextPoint, -30)) { // 30 * ~30 um = 1 mm
         removeTrack(At(i));
         nTracksRemoved++;
      }

      delete nextPoint;
   }

   Compress();

//   cout << "Tracks::removeTracksLeavingDetector() has removed " << nTracksRemoved  << " tracks.\n";
}


void Tracks::removeTrackCollisions() {
   // sometimes two tracks end in the same cluster
   // remove the last cluster from the smallest track IF THE SMALLEST TRACK IS <= 2.. delete the track
   // if the final length is 1
   
   vector<Int_t> * conflictPair = nullptr;
   Clusters * conflictClusters[2] = {};
   Track * track[2] = {};
   Int_t len[2] = {0};
   Int_t longTrack, shortTrack;

   Int_t nRemoved = 0;
   Int_t nShortRemoved = 0;

   vector<Int_t> * conflictTracks = getTracksWithConflictClusters();

   for (UInt_t i=0; i<conflictTracks->size(); i++) {
      conflictPair = getConflictingTracksFromTrack(conflictTracks->at(i)); // 2 tracks
      if (conflictPair->size() < 2) continue;

      for (Int_t j=0; j<2; j++) {
         track[j] = At(conflictPair->at(j));
         conflictClusters[j] = track[j]->getConflictClusters(); // 1 cluster
         len[j] = track[j]->GetEntriesFast();
      }
      
      longTrack = (len[0] > len[1]) ? 0 : 1;
      shortTrack = 1 - longTrack;

      Bool_t sameCluster = isSameCluster(conflictClusters[0]->At(0), conflictClusters[1]->At(0));
      if (!sameCluster) { continue; }
      
      Cluster * conflictCluster = conflictClusters[0]->At(0); // only one conflicting cluster!

      if (!conflictCluster) continue;

      Int_t longClusterIdx = track[longTrack]->getClusterIdx(conflictCluster);
      Int_t shortClusterIdx = track[shortTrack]->getClusterIdx(conflictCluster);

      if (shortClusterIdx != track[shortTrack]->GetEntriesFast() - 1) {
         // not a terminating cluster!
         continue;
      }

      track[longTrack]->At(longClusterIdx)->markUnused();
   
      if (len[shortTrack] > 2) {
         track[shortTrack]->removeClusterAt(shortClusterIdx);
      }
      
      else {
         removeTrack(track[shortTrack]);
         nShortRemoved++;
      }
      nRemoved++;
   }

   cout << "Tracks::removeTrackCollisions removed " << nRemoved << " collisions, and removed " << nShortRemoved << " too short tracks.\n";
}

void Tracks::fillOutIncompleteTracks(float angleLimit) {
   Track     * thisTrack = nullptr;
   Int_t       lastLayer, nInTrack;
   Float_t     angle, minAngle, minIdx;
   Cluster   * lastCluster;
   Cluster   * nextToLastCluster;
   Cluster   * clusterWithoutTrack;
   Bool_t      continueSearch;

   for (Int_t i=GetEntriesFast()-1; i>=0; i--) {
      thisTrack = At(i);
      if (!thisTrack) continue;
//      if (thisTrack->Last()->getDepositedEnergy() > 3) continue;
      if (thisTrack->GetEntriesFast() < 4) continue;

//      printf("Working on track %d with length %d\n", i, thisTrack->Last()->getLayer());
      
      nInTrack = thisTrack->GetEntriesFast();
      nextToLastCluster = thisTrack->At(nInTrack-2);
      lastCluster = thisTrack->At(nInTrack-1);
      if (!nextToLastCluster || !lastCluster) continue;
      lastLayer = lastCluster->getLayer();

      // check CWT if a nearby track exists
      continueSearch = true;
      while (continueSearch) {
         continueSearch = false;

         minAngle = 1e5;
         minIdx = -1;
         for (Int_t j=0; j<GetEntriesFastCWT(); j++) {
            clusterWithoutTrack = AtCWT(j);
            if (!clusterWithoutTrack) continue;

            if (clusterWithoutTrack->getLayer() == lastLayer+1) {
               angle = getDotProductAngle(nextToLastCluster, lastCluster, clusterWithoutTrack);
               if (angle < minAngle) {
                  minAngle = angle;
                  minIdx = j;
               }
            }
         }
         
         if (minIdx>=0) {
            if (minAngle <= angleLimit) {
               int eid1 = thisTrack->Last()->getEventID();
               int eid2 = AtCWT(minIdx)->getEventID();

//               printf("Found CWT @ layer %d - angle = %.1f mrad (OK? %d)\n", AtCWT(minIdx)->getLayer(), minAngle*1000, (eid1 == eid2 || eid1*eid2 < 0));
               thisTrack->appendCluster(AtCWT(minIdx));
               removeCWTAt(minIdx);
               continueSearch = true;
               
               nInTrack = thisTrack->GetEntriesFast();
               nextToLastCluster = thisTrack->At(nInTrack-2);
               lastCluster = thisTrack->At(nInTrack-1);
               if (!nextToLastCluster || !lastCluster) continue;
               lastLayer = lastCluster->getLayer();
            }
            else {
               continueSearch = true;
               
               nInTrack = thisTrack->GetEntriesFast();
               nextToLastCluster = thisTrack->At(nInTrack-2);
               lastCluster = thisTrack->At(nInTrack-1);
               if (!nextToLastCluster || !lastCluster) continue;
               lastLayer++;
               if (lastLayer == lastCluster->getLayer() + 1) continueSearch = false;
            }
         }
      }
   }
}

void Tracks::retrogradeTrackImprovement(Clusters * clusters) {
   Track * thisTrack = nullptr;
   Cluster * estimatedRetrogradeCluster = nullptr;
   Cluster * nearestNeighborAllTracks = nullptr;
   Cluster * actualCluster = nullptr;
   Float_t deltaThisTrack = 0;
   Float_t deltaBestTrack = 0;
   Bool_t anotherTrackIsCloserMatch;

   Track * switchTrack = nullptr;
   Int_t switchTrackIdx = -1;
   Int_t switchTrackAtLayer = -1;
   Int_t idxThis, idxSwitch;

   Tracks * newTracks = new Tracks();

   Bool_t ignoreTrack[GetEntriesFast()];
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      ignoreTrack[i] = false;
   }

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (ignoreTrack[i]) continue;

      switchTrack = nullptr;
      switchTrackIdx = -1;
      switchTrackAtLayer = -1;

      thisTrack = At(i);
      if (!thisTrack) continue;

      Int_t lastIdx = thisTrack->GetEntriesFast() - 1;
      Int_t lastLayer = thisTrack->getLayer(lastIdx);
      
      for (Int_t layer = lastLayer - 2; layer>=0; layer--) {
         if (thisTrack->getIdxFromLayer(layer) < 0) continue;

         estimatedRetrogradeCluster = getRetrogradeTrackExtrapolationToLayer(thisTrack, layer);
         if (!estimatedRetrogradeCluster) continue;
         
         actualCluster = thisTrack->At(thisTrack->getClusterFromLayer(layer));
         nearestNeighborAllTracks = clusters->findNearestNeighbour(nullptr, estimatedRetrogradeCluster, false); // BROKEN!!!

         if (!nearestNeighborAllTracks) {
            delete estimatedRetrogradeCluster;
            continue;
         }

         deltaThisTrack = diffXY(estimatedRetrogradeCluster, actualCluster);
         deltaBestTrack = diffXY(estimatedRetrogradeCluster, nearestNeighborAllTracks);
         delete estimatedRetrogradeCluster;
         
         anotherTrackIsCloserMatch = (deltaBestTrack < deltaThisTrack);

         if (anotherTrackIsCloserMatch) {
            switchTrackIdx = getTrackIdxFromCluster(nearestNeighborAllTracks);
            switchTrack = At(switchTrackIdx);
            switchTrackAtLayer = layer;
            break;
         }
      }

      if (switchTrackIdx == i) {
         switchTrack = nullptr;
      }

      if (switchTrack) {
         Track * newTrackA = new Track();
         Track * newTrackB = new Track();
         
         Int_t longestTrack = max(thisTrack->Last()->getLayer(), switchTrack->Last()->getLayer());
         
         for (Int_t layer=0; layer<=longestTrack; layer++) {
            idxThis = thisTrack->getIdxFromLayer(layer);
            idxSwitch = switchTrack->getIdxFromLayer(layer);

            Cluster *clusterThis   = (idxThis>=0) ? thisTrack->At(idxThis) : nullptr;
            Cluster *clusterSwitch = (idxSwitch>=0) ? switchTrack->At(idxSwitch) : nullptr;
            
            if (layer <= switchTrackAtLayer) { // SWITCH TRACKS
               if (clusterThis)   newTrackB->appendCluster(clusterThis);
               if (clusterSwitch) newTrackA->appendCluster(clusterSwitch);
            }

            else { // KEEP ORIGINAL TRACKS
               if (clusterThis)   newTrackA->appendCluster(clusterThis);
               if (clusterSwitch) newTrackB->appendCluster(clusterSwitch);
            }
         }

         Float_t sin1 = thisTrack->getSinuosity();
         Float_t sin2 = switchTrack->getSinuosity();

         Float_t sin3 = newTrackA->getSinuosity();
         Float_t sin4 = newTrackB->getSinuosity();
         
         Float_t metricOld = quadratureAdd(sin1, sin2);
         Float_t metricNew = quadratureAdd(sin3, sin4);

         if (metricNew < metricOld) {
            newTracks->appendTrack(newTrackA);
            newTracks->appendTrack(newTrackB);
            removeTrackAt(switchTrackIdx);
            removeTrackAt(i);
         }

         ignoreTrack[switchTrackIdx] = true;
         delete newTrackA;
         delete newTrackB;
      }
   }

   cout << "Track::retrogradeTrackImprovement found " << newTracks->GetEntries() << " improved tracks:\n";
   for (Int_t i=0; i<newTracks->GetEntriesFast(); i++) {
      cout << newTracks->At(i) << ": " << *newTracks->At(i) << endl;
   }

   if (newTracks->GetEntriesFast()) {
      for (Int_t i=0; i<newTracks->GetEntriesFast(); i++) {
         appendTrack(newTracks->At(i));
      }
   }
   delete newTracks;
}

void Tracks::removeTracksEndingInBadChannels() {
   Track     * thisTrack = nullptr;
   Cluster   * extrapolatedCluster = nullptr;
   Int_t       lastLayer;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;
      lastLayer = thisTrack->Last()->getLayer();
      if (thisTrack->Last()->getDepositedEnergy() > 50) continue; // OK, this may be a bragg peak
      
      extrapolatedCluster = getTrackExtrapolationToLayer(thisTrack, lastLayer + 1);
      if (isBadData(extrapolatedCluster)) {
         cout << "Removing track " << *thisTrack << " due to bad channel at " << *extrapolatedCluster << "\n";
         removeTrackAt(i);
      }
      delete extrapolatedCluster;
   }

   removeEmptyTracks();
   Compress();
}

void Tracks::removeTracksEndingInGapRegion() {
   Track     * thisTrack = nullptr;
   Cluster   * extrapolatedCluster = nullptr;
   Int_t       lastLayer, removedTracks = 0;
   Float_t     y;

   // FIX IN FINAL DETECTOR !!!
   Float_t     gapYfrom = -0.5;
   Float_t     gapYto = 0.5;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;
      lastLayer = thisTrack->Last()->getLayer();

      if (lastLayer <= 1) continue;
      
      extrapolatedCluster = getTrackExtrapolationToLayer(thisTrack, lastLayer + 1);
      y = extrapolatedCluster->getYmm();
      delete extrapolatedCluster;

      if (gapYfrom < y && y < gapYto) {
         removeTrack(thisTrack);
         removedTracks++;
         continue;
      }
   }

   removeEmptyTracks();
   Compress();

   printf("Removed %d tracks ending in gap region\n", removedTracks);
}

void Tracks::removeTracksEndingInHalo(Float_t haloRadius) {
   Track     * thisTrack = nullptr;
   Cluster   * extrapolatedCluster = nullptr;
   Int_t       lastLayer, removedTracks = 0;
   Float_t     x,y;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;
      lastLayer = thisTrack->Last()->getLayer();
      if (thisTrack->Last()->getDepositedEnergy() > 15) continue;

      if (lastLayer <= 1) continue;

      extrapolatedCluster = getTrackExtrapolationToLayer(thisTrack, lastLayer + 1);
      if (!extrapolatedCluster) continue;

      x = extrapolatedCluster->getXmm();
      y = extrapolatedCluster->getYmm();
      delete extrapolatedCluster;

      if (sqrt(x*x+y*y) > haloRadius) {
         removeTrack(thisTrack);
         removedTracks++;
         continue;
      }
   }
   removeEmptyTracks();
   Compress();

   printf("Removed %d tracks ending in halo region.\n", removedTracks);
}

void Tracks::removeHighChiSquare(Float_t chi2limit) {
   Track  * thisTrack = nullptr;
   Int_t    nRemoved = 0;
   Int_t    nRemovedNuclear = 0;
   Float_t  chi2;
   Int_t    nTotal = GetEntriesFast();
   Int_t    nTotalStart = GetEntries();

   for (Int_t i=0; i<nTotal; i++) {
      thisTrack = At(i);
      if (!At(i)) continue;
      chi2 = thisTrack->getFitParameterChiSquare();
      if (chi2>chi2limit) {
         if (thisTrack->Last()->isSecondary()) nRemovedNuclear++;
         removeTrack(thisTrack);
         nRemoved++;
      }
   }

   removeEmptyTracks();
   Compress();

   Float_t fraction = nTotal ? nRemoved/float(nTotalStart) : 0;
   cout << "Tracks::removeHighChiSquare() removed " << nRemoved << " of " << nTotalStart << " (" << fraction * 100 << "%) tracks (" << 100*float(nRemovedNuclear)/nRemoved << "% secondaries).\n";
}

void Tracks::removeNuclearInteractions() {
   Track  * thisTrack = nullptr;
   Int_t    nRemoved = 0;
   Int_t    nRemovedNuclear = 0;
   Int_t    nTotal = GetEntriesFast();
   Int_t    nTotalStart = GetEntries();

   Float_t plateauThreshold = 0;
   if (kHelium) plateauThreshold = 3.5;

   Float_t plateauEdep;

   for (Int_t i=0; i<nTotal; i++) {
      thisTrack = At(i);
      if (!At(i)) continue;
      plateauEdep = thisTrack->getAverageDepositedEnergy(0, thisTrack->GetEntriesFast()-5);
      if (thisTrack->doesTrackEndAbruptly()  || plateauEdep < plateauThreshold) {
         if (thisTrack->Last()->isSecondary()) nRemovedNuclear++;
         removeTrack(thisTrack);
         nRemoved++;
      }
   }

   removeEmptyTracks();
   Compress();

   Float_t fraction = nTotal ? nRemoved/float(nTotalStart) : 0;
//   cout << "Tracks::removeNuclearInteractions() removed " << nRemoved << " of " << nTotalStart << " (" << fraction * 100 << "%) tracks (" << 100*float(nRemovedNuclear)/nRemoved << "% secondaries).\n";
}

void Tracks::removeTracksWithMinWEPL(Float_t minWEPL) {
   Track  * thisTrack = nullptr;
   Int_t    nRemoved = 0;
   Int_t    nRemovedNuclear = 0;
   Int_t    nTotal = GetEntriesFast();
   Int_t    nTotalStart = GetEntries();
   Int_t    fpr = 0;

   for (Int_t i=0; i<nTotal; i++) {
      thisTrack = At(i);
      if (!At(i)) continue;

      fpr = thisTrack->getFitParameterRange();
      if (fpr < minWEPL) {
         if (thisTrack->Last()->isSecondary()) nRemovedNuclear++;
         removeTrack(thisTrack);
         nRemoved++;
      }
   }

   removeEmptyTracks();
   Compress();
   
      Float_t fraction = nTotal ? nRemoved/float(nTotalStart) : 0;
//   cout << "Tracks::removeTracksWithMinWEPL(" << minWEPL << ") removed " << nRemoved << " of " << nTotalStart << " (" << fraction * 100 << "%) tracks (" << 100*float(nRemovedNuclear)/nRemoved << "% secondaries).\n";
}


void Tracks::removeThreeSigmaShortTracks() {
   Track  * thisTrack = nullptr;
   Int_t    nRemoved = 0;
   Int_t    nRemovedNuclear = 0;
   Int_t    nTotal = GetEntriesFast();
   Int_t    nTotalStart = GetEntries();

   Float_t  muRange = 0, sigmaRange = 0, cutRangeLow, cutRangeHigh;
   Int_t    nInSum = 0;

   Float_t fpm;

   TH1I * h = new TH1I("h","h",100,0,330);

   for (Int_t i=0; i<=nTotal; i++) {
      thisTrack = At(i);
      if (!At(i)) continue;
      fpm = getUnitFromTL(thisTrack->getRangemm());
      if (!isnan(fpm)) h->Fill(fpm);
   }

   muRange = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
   sigmaRange = 10;

   delete h;

   cutRangeLow = muRange - 3 * sigmaRange;
   cutRangeHigh = muRange + 3 * sigmaRange;
   printf("muRange = %.3f mm, sigmaRange = %.3f mm, cutRangeLow = %.3f mm, cutRangeHigh = %.3f mm\n", muRange, sigmaRange, cutRangeLow, cutRangeHigh);

   for (Int_t i=0; i<=nTotal; i++) {
      thisTrack = At(i);
      if (!At(i)) continue;
      if (getUnitFromTL(thisTrack->getRangemm()) < cutRangeLow || getUnitFromTL(thisTrack->getRangemm()) > cutRangeHigh || thisTrack->getLayer(0) > 2) {
         if (thisTrack->Last()->isSecondary()) nRemovedNuclear++;
         removeTrack(thisTrack);
         nRemoved++;
      }
   }

   removeEmptyTracks();
   Compress();

   Float_t fraction = nTotal ? nRemoved/float(nTotalStart) : 0;
//   cout << "Tracks::removeThreeSigmaShortTracks() removed " << nRemoved << " of " << nTotalStart << " (" << fraction * 100 << "%) tracks (" << 100*float(nRemovedNuclear)/nRemoved << "% secondaries).\n";
}

void Tracks::removeEmptyTracks() {
   for (int i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (!GetEntriesFast(i)) {
         removeTrackAt(i);
      }
   }
   Compress();
}   

#endif
