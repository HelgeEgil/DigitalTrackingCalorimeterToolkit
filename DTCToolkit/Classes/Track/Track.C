#include <iostream>
#include <cmath>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TClonesArray.h>
#include <TF1.h>

#include "Classes/Track/Track.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Hit/Hit.h"
#include "HelperFunctions/Tools.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"

using namespace std;

Track::Track() : track_("Cluster", MaxTrackLength) {
   fitEnergy_ = 0;
   fitRange_ = 0;
   fitScale_ = 0;
   fitError_ = 0;
}

Track::Track(Cluster *cluster) : track_("Cluster", MaxTrackLength) {
   appendCluster(cluster);
   fitEnergy_ = 0;
   fitRange_ = 0;
   fitScale_ = 0;
   fitError_ = 0;
}

Track::~Track() {
   track_.Delete();
}

void Track::setTrack(Track *copyTrack, Int_t startOffset /* default 0 */) {
   clearTrack();

   for (Int_t i=0; i<copyTrack->GetEntriesFast(); i++) {
      appendCluster(copyTrack->At(i), startOffset);
   }
}

void Track::appendCluster(Cluster *copyCluster, Int_t startOffset /* default 0 */) {
   // new object, no pointers are copied
   if (!copyCluster) return;
   
   Int_t    i = GetEntriesFast();
   if (i==0) i = startOffset;

   Cluster *c = (Cluster*) track_.ConstructedAt(i);
   c->set(copyCluster); 
}

void Track::appendPoint(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) track_.ConstructedAt(i);
   c->set(x,y,layer,size);
   c->setEventID(eventID);
}

Int_t Track::getModeEventID() {
   // Return the mode (median without averaging similarly frequent points)
   // of the event ID distribuition of a single track
   // Used to determine the 'best' estimate of the identity of the track
   
   const Int_t n = GetEntriesFast();
   Bool_t      metThisID = false;
   Int_t       nUniqueEventIDs = 0;
   Int_t       thisID = -1;
   Int_t       nThisID = 0;
   Int_t       eventIDs[n];
   Int_t       repeatArray[n];
   Int_t       maxIdx = 0;

   for (Int_t i=0; i<n; i++) {
      eventIDs[i] = getEventID(i);
      repeatArray[i] = 0;
   }

   for (Int_t i=0; i<n; i++) {
      thisID = eventIDs[i];
      
      metThisID = false;

      for (Int_t j=0; j<i; j++) {
         if (thisID == eventIDs[j]) {
            metThisID = true;
            break;
         }
      }

      if (!metThisID) {
         nThisID = 0;
         for (Int_t j=0; j<n; j++) {
            if (thisID == eventIDs[j]) nThisID++;
         }
      
         repeatArray[i] = nThisID;
         nUniqueEventIDs++;

      }
   }

   for (Int_t i=0; i<n; i++) {
      if (repeatArray[i] > repeatArray[maxIdx]) {
         maxIdx = i;
      }
   }

   return eventIDs[maxIdx];
}

Bool_t Track::isOneEventID() {
   Int_t eid;

   if (!At(0)) {
      if (!At(1)) return false;
      eid = getEventID(1);
   }

   else { 
      eid = getEventID(0);
      if (eid < 0 && At(1)) { eid = getEventID(1); }
      if (eid < 0) return false;
   }

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (eid != getEventID(i) && getEventID(i) > 0) return false;
   }

   return true;
}

Bool_t Track::isFirstAndLastEventIDEqual() {
   if (!At(0)) return false;

   Int_t first = getEventID(0);
   Int_t last = Last()->getEventID();

   return (first == last || last < 0);
}

ostream& operator<< (ostream &os, Track& t) {
   int n = t.GetEntriesFast();
   os << "[";
   for (int i=0; i<n; i++) {
      if (t.At(i)) {
         os << i << *t.At(i);
         if (i<n-1) os << ", ";
      }
   }
   os << "]";
   return os;
}
