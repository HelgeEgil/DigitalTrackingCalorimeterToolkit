#ifndef Track_cxx
#define Track_cxx

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
using namespace DTC;

Track::Track() : track_("DTC::Cluster", MaxTrackLength) {
   track_.SetOwner(kTRUE);
   track_.SetBit(kCanDelete);
   track_.SetBit(kMustCleanup);

   fitEnergy_ = 0;
   fitRange_ = 0;
   fitScale_ = 0;
   fitError_ = 0;
   isIncomplete_ = false;
}

Track::Track(Cluster *cluster) : track_("DTC::Cluster", MaxTrackLength) {
   track_.SetOwner(kTRUE);
   track_.SetBit(kCanDelete);
   track_.SetBit(kMustCleanup);

   appendCluster(cluster);
   fitEnergy_ = 0;
   fitRange_ = 0;
   fitScale_ = 0;
   fitError_ = 0;
   isIncomplete_ = false;
}

Track::~Track() {
   Clear("C");
   track_.Delete();
}

Int_t Track::Compare(const TObject *obj) const {
   Int_t length = ((Track*) this)->Last()->getLayer();
   if (length == ((Track*) obj)->Last()->getLayer()) return 0;
   else if (length > ((Track*) obj)->Last()->getLayer()) return -1;
   else return 1;
}

void Track::setTrack(Track *copyTrack, Int_t startOffset /* default 0 */) {
   Clear("C");

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

void Track::appendPoint(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID, Bool_t isSecondary, Int_t PDG) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) track_.ConstructedAt(i);
   c->set(x,y,layer,size, eventID, isSecondary, PDG);
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

Bool_t Track::isOneEventID() { // Alba's version
   Int_t eid;
   Int_t last_pos = GetEntriesFast()-1;

    if (!At(last_pos)) {
      if (!At(last_pos-1)) return false;
      eid = getEventID(last_pos-1);
    }

   else { 
      eid = getEventID(last_pos);
      if (eid < 0 && At(last_pos-1)) { eid = getEventID(last_pos-1); }
      if (eid < 0) return false;
   }

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      if (eid != getEventID(i) && getEventID(i) > 0) return false;
   }

   return true;
}

/*
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
*/

Bool_t Track::isFirstAndLastEventIDEqual() {
   if (!At(0)) return false;

   Int_t first = getEventID(0);
   Int_t last = Last()->getEventID();

   return (first == last);
}

void  Track::propagateSecondaryStatus() {
   Bool_t isSecondary = false;

   for (Int_t c=0; c<GetEntriesFast(); c++) {
      if (!At(c)) continue;
      if (At(c)->isSecondary()) {
         isSecondary = true;
         break;
      }
   }
   if (isSecondary) {
      for (Int_t c=0; c<GetEntriesFast(); c++) {
         if (!At(c)) continue;
         At(c)->setSecondary(true);
      }
   }
}

void Track::removeNANs() {
   Cluster *c = nullptr;

   for (Int_t i=0; i<=GetEntriesFast(); i++) {
      c = At(i);
      if (!c) continue;
      if (std::isnan(c->getX()) || std::isnan(c->getY())) {
         removeCluster(c);
      }
   }
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

#endif
