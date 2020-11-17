#ifndef Tracks_cxx
#define Tracks_cxx

#include <iostream>
#include <vector>

#include <TClonesArray.h>
#include <TObject.h>
#include <TCollection.h>
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

using namespace std;
using namespace DTC;

Tracks::~Tracks() {
   Clear("C+C");
   tracks_.Delete();
   clustersWithoutTrack_.Delete();
}

void Tracks::CompressClusters() {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      At(i)->Compress();
   }
}

void Tracks::Clear(Option_t *option) {
   tracks_.Clear(option);
   clustersWithoutTrack_.Clear(option);
}

Long64_t Tracks::Merge(TCollection *tlist) {
   if (tlist) {
      Tracks * otherTracks = nullptr;
      TIter nOtherTracks(tlist);
      while ((otherTracks = (Tracks*) nOtherTracks())) {
         for (Int_t i=0; i<otherTracks->GetEntriesFast(); i++) {
            if (otherTracks->At(i)) {
               appendTrack(otherTracks->At(i));
            }
         }
         appendClustersWithoutTrack(otherTracks->getClustersWithoutTrack());
      }
   }
   return 1;
}
      
void Tracks::removeTrack(Track *t) {
   Cluster *thisCluster = nullptr;

   for (Int_t i=0; i<t->GetEntriesFast(); i++) {
      thisCluster = t->At(i);
      if (!thisCluster) continue;
      appendClusterWithoutTrack(thisCluster);
   }

   tracks_.Remove((TObject*) t);
}

TObject*  Tracks::removeTrackAt(Int_t idx) {
   Cluster *thisCluster = nullptr;
   Track *t = At(idx);

   for (Int_t i=0; i<t->GetEntriesFast(); i++) {
      thisCluster = t->At(i);
      if (!thisCluster) continue;
      appendClusterWithoutTrack(thisCluster);
   }

   TObject *to = (TObject*) tracks_.RemoveAt(idx);
   return to;
}

void Tracks::sortTracks() {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      At(i)->sortTrack();
   }
}

void Tracks::appendTrack(Track *copyTrack, Int_t startOffset /* default 0 */) {
   Int_t    newIdx = tracks_.GetEntriesFast();
   Track  * track = (Track*) tracks_.ConstructedAt(newIdx);
   
   showDebug("Copying track ... i = ");
   for (Int_t i=0; i<copyTrack->GetEntriesFast(); i++) {
      showDebug(i << ", ");
      if(!copyTrack->At(i))
         continue;

      track->appendCluster(copyTrack->At(i), startOffset);
   }
   track->setIncomplete(copyTrack->isIncomplete());
   showDebug(endl);
}

void Tracks::appendClusterWithoutTrack(Cluster *cluster) {
   Int_t i = clustersWithoutTrack_.GetEntriesFast();
   Cluster *c = (Cluster*) clustersWithoutTrack_.ConstructedAt(i+1);
   c->set(cluster);
}

void Tracks::appendClustersWithoutTrack(TClonesArray *copyCWT) {
   Int_t       idxFrom = clustersWithoutTrack_.GetEntriesFast();
   Cluster   * newCluster = nullptr;
   Cluster   * copyCluster = nullptr;

   for (Int_t i=0; i<copyCWT->GetEntriesFast(); i++) {
      copyCluster = (Cluster*) copyCWT->At(i);
      if (!copyCluster) continue;

      newCluster = (Cluster*) clustersWithoutTrack_.ConstructedAt(idxFrom + i);
      newCluster->set(copyCluster);
   }
}

vector<Int_t> * Tracks::getTracksWithConflictClusters() {
   const Int_t       n = GetEntriesFast();
   Track           * thisTrack = nullptr;
   vector<Int_t>   * trackIdx = new vector<Int_t>;

   trackIdx->reserve(n/10.);

   for (Int_t i=0; i<n; i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;
      if (thisTrack->isUsedClustersInTrack()) { 
         trackIdx->push_back(i);
      }
   }

   return trackIdx;
}

vector<Int_t> * Tracks::getConflictingTracksFromTrack(Int_t trackIdx) {
   Int_t             nClusters, nClustersReal;
   Cluster         * conflictCluster = nullptr;
   Clusters        * conflictClusters = nullptr;
   vector<Int_t>   * conflictingTracks = new vector<Int_t>;
   vector<Int_t>   * possibleTracks = nullptr;
   conflictingTracks->reserve(5);
   
   Track * trackA = At(trackIdx);
   if (!trackA) { return conflictingTracks; }

   conflictClusters = trackA->getConflictClusters();
   nClusters = conflictClusters->GetEntriesFast();
   nClustersReal = conflictClusters->GetEntries();

   for (Int_t i=0; i<nClusters; i++) {
      conflictCluster = conflictClusters->At(i);
      if (!conflictCluster) continue;

      possibleTracks = getTracksFromCluster(conflictCluster);

      for (UInt_t j=0; j<possibleTracks->size(); j++) {
         if (!isItemInVector(possibleTracks->at(j), conflictingTracks)) {
            conflictingTracks->push_back(possibleTracks->at(j));
         }
      }
   }

   return conflictingTracks;
}


void Tracks::extrapolateToLayer0() {
   for (Int_t i=0; i<tracks_.GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (!At(i)->At(0)) { // only do if track At(i) misses the first cluster
         At(i)->extrapolateToLayer0();
      }
   }
}

void Tracks::doTrackFit(Bool_t freeScale) {
   TGraphErrors * out = nullptr;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      out = At(i)->doTrackFit(freeScale);
      delete out;
   }
}

Int_t Tracks::getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster) {
   Float_t     minDist = 1e5, dist;
   Int_t       minIdx = -1;
   Cluster   * thisCluster = nullptr;

   for (UInt_t i=0; i<clusters.size(); i++) {
      trackCluster thisTrackCluster = clusters.at(i);
      thisCluster = At(thisTrackCluster.track)->At(thisTrackCluster.cluster);
      
      dist = diffmmXY(thisCluster, interpolatedCluster);
      if (dist < minDist) {
         minDist = dist;
         minIdx = i;
      }
   }
   return minIdx;
}

vector<Int_t> * Tracks::getTracksFromCluster(Cluster * cluster) {
   Track           * thisTrack = nullptr;
   vector<Int_t>   * tracksWithCluster = new vector<Int_t>;
   Int_t             idx;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;

      if (thisTrack->isClusterInTrack(cluster)) {
         idx = thisTrack->getClusterIdx(cluster);
         thisTrack->At(idx)->markUsed();
         tracksWithCluster->push_back(i);
      }
   }

   return tracksWithCluster;
}  

Int_t Tracks::getTrackIdxFromFirstLayerEID(Int_t eventID) {
   Cluster  * thisCluster = nullptr;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      thisCluster = At(i)->At(0);

      if (!thisCluster) continue;

      if (eventID == thisCluster->getEventID()) return i;
   }

   return -1;
}

Int_t Tracks::getTrackIdxFromLastLayerEID(Int_t eventID) {
   Cluster  * thisCluster = nullptr;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      thisCluster = At(i)->Last();

      if (!thisCluster) continue;

      if (eventID == thisCluster->getEventID()) return i;
   }

   return -1;
}

Int_t Tracks::getTrackIdxFromCluster(Cluster * cluster) {
   Track * thisTrack = nullptr;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;

      if (thisTrack->isClusterInTrack(cluster)) {
         return i;
      }
   }
   return -1;
}

void Tracks::matchWithEventIDs(Hits * eventIDs) {
   // Use the Monte Carlo truth list eventIDs (x, y, layer, actual event)
   // and find the closes match for each item in list
   // then set truth eventID to the Tracks->Track->Cluster object

   Float_t     minDist = 1e5; // px
   Float_t     thisDist = 0;
   Track     * thisTrack = nullptr;
   Cluster   * thisCluster = nullptr;
   Hit       * thisHit = nullptr;
   Int_t       layer = -1;
   Int_t       minIdx = 0;
   Bool_t      doLoop = true;
   Int_t       nHits = eventIDs->GetEntriesFast();
   Float_t     cX, cY;
   Int_t       nClusters = 0;

   for (Int_t t=0; t<GetEntriesFast(); t++) {
      thisTrack = At(t);
      if (!thisTrack) continue;

      nClusters += thisTrack->GetEntriesFast();

      for (Int_t c=0; c<thisTrack->GetEntriesFast(); c++) {
         thisCluster = thisTrack->At(c);
         layer = thisCluster->getLayer();

         cX = thisCluster->getX();
         cY = thisCluster->getY();

         // Optimization:
         // If thisCluster+1 is also eventIDs+1, don't loop to see if we can find it
         // But instead just use minIdx++ ;-)

         minDist = diffXY(thisCluster, eventIDs->At(minIdx+1));

         if (minDist < 10) {
            minIdx++;
            doLoop = false;
         }

         else {
            doLoop = true;
            minDist = 1e5;
            minIdx = -1;
         }

         if (doLoop) {
            for (Int_t h=0; h<eventIDs->GetEntriesFast(); h++) {
               thisHit = eventIDs->At(h);
               
               if (!thisHit) continue;
               if (thisHit->getLayer() != layer) continue;

               if (fabs(cX - thisHit->getX()) < 10) {
                  if (fabs(cY - thisHit->getY()) < 10) {

                     thisDist = diffXY(thisCluster, thisHit);
                     if (thisDist < minDist) {
                        minDist = thisDist;
                        minIdx = h;
                     }
                  }
               }
            }
         }
         
         if (minIdx >= 0 && minDist < 14) {
            thisCluster->setEventID(eventIDs->getEventID(minIdx));
            eventIDs->removeHitAt(minIdx);
         }
      }
   }

   Int_t cWithoutEventID = 0;
   for (Int_t t=0; t<GetEntriesFast(); t++) {
      for (Int_t c=0; c<At(t)->GetEntriesFast(); c++) {
         if (At(t)->getEventID(c) < 0) {
            cWithoutEventID++;
         }
      }
   }

   cout << "Number of clusters without eventID: " << cWithoutEventID << "( " << (float) cWithoutEventID / nClusters * 100 << "%)\n";
}

void Tracks::checkLayerOrientation() {
   Float_t     n[nLayers-1], x[nLayers-1], y[nLayers-1];
   Float_t     dx, dy;
   Track     * thisTrack = nullptr;
   Cluster   * lastCluster = nullptr;
   Cluster   * thisCluster = nullptr;
   Int_t       first = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!At(first)) first = 1;
      if (!At(first)) continue;

      lastCluster = thisTrack->At(first);
      for (Int_t j=first+1; j<thisTrack->GetEntriesFast(); j++) {
         thisCluster = thisTrack->At(j);
         dx = thisCluster->getX() - lastCluster->getX();
         dy = thisCluster->getY() - lastCluster->getY();

         x[j] += dx;
         y[j] += dy;
         n[j]++;

         lastCluster = thisCluster;
      }
   }

   cout << "The average track movement between layers is: ";
   for (Int_t i=0; i<nLayers-1; i++) {
      if (n[i]) {
         cout << Form("Layer %d: (%.1f, %.1f), ", i, x[i]/n[i], y[i]/n[i]);
      }
   }
   cout << ".\n";
}

Bool_t Tracks::isLastEventIDCloseToFirst(Int_t trackIdx) {
   Track     * track = At(trackIdx);
   Cluster   * comparisonCluster = nullptr;
   Int_t       lastEventID = track->Last()->getEventID();
   Int_t       firstEventID = track->getEventID(0);
   Int_t       comparisonEventID;
   Float_t     deltaXY, deltaPHI, phi1, phi2;
   Int_t       comparisonIdx;
   
   if (!track)                               return false;
   if (!track->At(0))                        return false;
   if (track->isFirstAndLastEventIDEqual())  return true;

   else {
      // check first against last ID
      comparisonIdx = getTrackIdxFromFirstLayerEID(lastEventID);

      if (comparisonIdx < 0) return true; 

      comparisonCluster = At(comparisonIdx)->At(0);
      
      if (!comparisonCluster) return false;

      deltaXY = diffmmXY(track->At(0), comparisonCluster);
      phi1 = track->getSlopeAngleAtLayer(1);
      phi2 = At(comparisonIdx)->getSlopeAngleAtLayer(1);
      deltaPHI = fabs(phi1 - phi2);

      if (deltaXY < 0.5 && deltaPHI < 0.5) {
         return true;
      }

      // check last against first ID
      comparisonIdx = getTrackIdxFromLastLayerEID(firstEventID);
      if (comparisonIdx < 0) return true;
      
      comparisonCluster = At(comparisonIdx)->Last();
      if (!comparisonCluster) return false;

      deltaXY = diffmmXY(track->Last(), comparisonCluster);
      phi1 = track->getSlopeAngleAtLayer(track->GetEntriesFast()-1);
      phi2 = At(comparisonIdx)->getSlopeAngleAtLayer(At(comparisonIdx)->GetEntriesFast()-1);
      deltaPHI = fabs(phi1 - phi2);
      
      if (deltaXY < 0.5 && deltaPHI < 0.5) {
         return true;
      }
   }

   return false;
}

void Tracks::removeShortTracks(Int_t len) {
   Int_t n=0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;

      if (At(i)->GetEntries() <= len) {
         removeTrackAt(i);
         n++;
      }
   }
   Compress();
//   printf("Removed %d tracks with >= %d entries.\n", n, len);
}


Int_t Tracks::getNMissingClustersWithEventID(Int_t eventID, Int_t afterLayer, Int_t beforeLayer) {
   Int_t n = 0;

   for (Int_t i=0; i<GetEntriesFastCWT(); i++) {
      if (!AtCWT(i)) continue;
      if (AtCWT(i)->isSecondary() && AtCWT(i)->getPDG() < 1000) continue;

      Bool_t missingClusterWithEventID = (AtCWT(i)->getEventID() == eventID);
      Bool_t missingInLastLayers = (AtCWT(i)->getLayer() > afterLayer);
      Bool_t missingInFirstLayers = (AtCWT(i)->getLayer() < beforeLayer);

      if (missingClusterWithEventID && (missingInLastLayers || missingInFirstLayers)) {
         n++;
      }
   }
   
   // Check also if kSuppressedClustersEventID contains cluster beyond afterlayer!
   for (UInt_t i=0; i<kSuppressedClustersEventID.size(); i++) {
      if (kSuppressedClustersEventID.at(i) == eventID) {
         if (kSuppressedClustersLayer.at(i) > afterLayer) {
            n++;
         }
      }
   }

   return n;
}

/*
void Tracks::createEIDSortList() {
   Int_t eventID;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      EIDindex_[i] = At(i)->getEventID(0);
      showDebug("Tracks::createEIDSortList: EIDindex_[" << i <<"] = " << At(i)->getEventID(0) << endl);
   }
}

Track * Tracks::getTrackWithEID(Int_t eid) {
   Int_t trackIdx = EIDindex_[eid];
   if (!At(trackIdx)) return nullptr;
   return At(trackIdx);
}
*/

void Tracks::propagateSecondaryStatus() {
   Track   *thisTrack = nullptr;
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!At(i)) continue;
      At(i)->propagateSecondaryStatus();
   }
}


void  Tracks::removeHighAngleTracks(Float_t mradLimit) {
   Cluster *a = nullptr;
   Cluster *b = nullptr;
   Cluster *a2 = nullptr;
   Cluster *b2 = nullptr;
   Track   *thisTrack = nullptr;
   Float_t  incomingAngle, outgoingAngle;
   Int_t    nRemoved = 0;
   Int_t    nRemovedNuclear = 0;
   Int_t    last;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!At(i)) continue;

      last = thisTrack->GetEntriesFast() - 1;

      if (last<1) continue;

      a = thisTrack->At(0);
      b = thisTrack->At(1);

      a2 = thisTrack->At(last);
      b2 = thisTrack->At(last-1);

      if (!a || !b) continue;
      incomingAngle = getDotProductAngle(a, a, b);
      outgoingAngle = getDotProductAngle(b2, b2, a2);
      if (incomingAngle > mradLimit / 1000.) {// || outgoingAngle > mradLimit / 200.) {
         if (thisTrack->Last()->isSecondary()) nRemovedNuclear++;
         removeTrackAt(i);
         nRemoved++;
      }
   }
   printf("Removed %d tracks with higher than %.0f mrad incoming angle (%.2f%% secondaries).\n", nRemoved, mradLimit, 100*float(nRemovedNuclear)/nRemoved);

   removeEmptyTracks();
   Compress();
}

void  Tracks::removeHighAngleTracksRelativeToSpot(Float_t mradLimit, Float_t angleXmrad, Float_t angleYmrad) {
   Cluster *a = nullptr;
   Cluster *b = nullptr;
   Cluster *acorr = nullptr;
   Cluster *a2 = nullptr;
   Cluster *b2 = nullptr;
   Track   *thisTrack = nullptr;
   Float_t  incomingAngle, outgoingAngle;
   Int_t    nRemoved = 0;
   Int_t    nRemovedNuclear = 0;
   Int_t    last;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!At(i)) continue;

      last = thisTrack->GetEntriesFast() - 1;

      if (last<1) continue;

      a = thisTrack->At(0);
      b = thisTrack->At(1);
      if (!a || !b) continue;

      acorr = new Cluster();
      acorr->setXmm(a->getXmm() - tan(angleXmrad/1000) * dz2);
      acorr->setYmm(a->getYmm() - tan(angleYmrad/1000) * dz2);
      acorr->setLayer(-1);
      
      incomingAngle = getDotProductAngle(acorr, a, b);
      delete acorr;

      if (incomingAngle > mradLimit / 1000.) {
         if (thisTrack->Last()->isSecondary()) nRemovedNuclear++;
         removeTrackAt(i);
         nRemoved++;
      }

   }
   printf("Removed %d tracks with higher than %.0f mrad incoming angle (%.2f%% secondaries).\n", nRemoved, mradLimit, 100*float(nRemovedNuclear)/nRemoved);

   removeEmptyTracks();
   Compress();
}

Bool_t Tracks::isTrackIncomplete(Track * originalTrack) {
   Int_t eventID = originalTrack->Last()->getEventID();
   Track * thisTrack = nullptr;

   /*
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;
      if (!thisTrack->At(0)) continue;
      if (thisTrack == originalTrack) continue;
      
      if (thisTrack->At(0)->getEventID() == eventID && !thisTrack->At(0)->isSecondary() && !originalTrack->Last()->isSecondary()) {
         //cout << "Tracks are similar! original = " << *originalTrack << endl << "similar = " << *thisTrack << endl;
         return true;
      }   
   }
   */

   if (getNMissingClustersWithEventID(eventID, originalTrack->Last()->getLayer(), originalTrack->At(0)->getLayer()) > 0) return true;


   return false;
}

void Tracks::removeTracksStartingInDetector() {
   Track * thisTrack = nullptr;
   Int_t n=0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!At(i)) continue;

      if (thisTrack->At(0)->getLayer() > 1) {
         removeTrackAt(i);
         n++;
      }

   }

//   printf("Removed %d tracks starting inside the detector\n.", n);
   Compress();

}


void  Tracks::removeHighAngularChangeTracks(Float_t mradLimit) {
   Track   *thisTrack = nullptr;
   Float_t  maxChange;
   Int_t    nRemoved = 0;
   Int_t    nRemovedNuclear = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!At(i)) continue;

      maxChange = thisTrack->getMaximumSlopeAngleChange() * 3.14159265 / 180 * 1000;

      if (maxChange > mradLimit) {
         if (thisTrack->Last()->isSecondary()) nRemovedNuclear++;
         removeTrackAt(i);
         nRemoved++;
      }
   }
   printf("Removed %d tracks with higher than %.0f mrad max angular change (%.2f%% secondaries).\n", nRemoved, mradLimit, 100*float(nRemovedNuclear)/nRemoved);

   removeEmptyTracks();
   Compress();
}



void Tracks::removeNANs() {
   Track * thisTrack = nullptr;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      thisTrack = At(i);
      if (!thisTrack) continue;

      thisTrack->removeNANs();
   }
   removeEmptyTracks();
}

#endif
