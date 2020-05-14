#ifndef Hits_cxx
#define Hits_cxx

#include <vector>
#include <algorithm>

#include <TStopwatch.h>

#include <TObject.h>
#include <TObjArray.h>
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace DTC;

Hits::~Hits() {
   hits_.Delete();
//   verticalIndexOfLayer_.clear();
   layerIndex_.clear();
}

void Hits::Clear(Option_t *option) {
//   verticalIndexOfLayer_.clear();
   layerIndex_.clear();
   hits_.Clear(option);
}

Long64_t Hits::Merge(TCollection *hlist) {
   if (hlist) {
      Hits * otherHits = nullptr;
      TIter nOtherHits(hlist);
      while ((otherHits = (Hits*) nOtherHits())) {
         appendHits(otherHits);
      }
   }
   return 1;
}

Int_t Hits::findLayerIndex(Int_t findLayer) {
   // Binary search of input ROOT file to file the first occurence of spot X
   // The file is of course sorted with increasing spotX values
   
   Int_t nentries = GetEntriesFast();
   Int_t nextIndex = nentries / 2;
   Int_t maxIndex = nentries-1;
   Int_t firstIndex = 0;

   if (nentries==0) return 0;
   if (findLayer == 0) return 0;

   Int_t nTries = 0;
   Int_t layer;
   while (nTries < 200) {
      nTries++;

      if (maxIndex - firstIndex == 1) break;

      layer = getLayer(nextIndex);

      if (layer >= findLayer) {
         maxIndex = nextIndex;
         nextIndex = (firstIndex + maxIndex) / 2; 
      }

      else if (layer < findLayer) {
         firstIndex = nextIndex;
         nextIndex = (firstIndex + maxIndex) / 2;
      }
   }

   return maxIndex;
}

void Hits::appendPoint(Int_t x, Int_t y, Int_t layer,  Float_t edep, Int_t eventID, Bool_t isSecondary, Int_t PDGEncoding) {
   Int_t i = GetEntriesFast();

   if (kConcatenateHits) {
      Bool_t added = false;
      Bool_t primary;

      // search through hits
      Int_t layerIdxFrom = findLayerIndex(layer);
      for (int j=layerIdxFrom; j<i; j++) {
         if (getLayer(j) > layer) break;

         if (getX(j) == x && getY(j) == y && getLayer(j) == layer) {
            primary = !this->isSecondary(j) || !isSecondary;
            At(j)->setSecondary(!primary);
         
            if (getPDG(j) < PDGEncoding) { // Use highest PDG in position
               At(j)->setEventID(eventID);
            }
            At(j)->setEdep(getEdep(j) + edep);
            added = true;
            break;
         }
      }
   
      if (!added) {
         Hit *hit = (Hit*) hits_.ConstructedAt(i);
         hit->set(x,y,layer,edep,eventID, isSecondary,PDGEncoding);
      }
   }

   else { 
      Hit *hit = (Hit*) hits_.ConstructedAt(i);
      hit->set(x,y,layer,edep,eventID, isSecondary,PDGEncoding);
   }
}

void Hits::appendHits(Hits *hits) {
  Int_t i = GetEntriesFast();

  for (Int_t j=0; j<hits->GetEntriesFast(); j++) {
   Hit *hit = (Hit*) hits_.ConstructedAt(i);
   hit->set(hits->At(i++));
  }
}

void Hits::appendHit(Hit *hit) {
   Int_t i = GetEntriesFast();

   Hit *newHit = (Hit*) hits_.ConstructedAt(i);
   newHit->set(hit);
}


Int_t Hits::getI(Int_t x, Int_t y) {
   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (x == getX(i) && y == getY(i))
         return i;
   }
   return -1;
}

void Hits::makeLayerIndex() {
   Int_t    lastLayer = -1;
   Int_t    startOffset = 0;
   Bool_t   kStarted = true;

   if (layerIndex_.size() == 0) {
      for (Int_t i=0; i<nLayers; i++)
         layerIndex_.push_back(-1);
   }

   else {
      for (Int_t i=0; i<nLayers; i++)
         layerIndex_.at(i) = -1;
   }

   if (!GetEntriesFast()) return;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i))
         continue;

      if (lastLayer != getLayer(i)) {
         layerIndex_.at(getLayer(i)) = i;
         lastLayer = getLayer(i);
      }
   }

   // set first indices to 0 if getLayer(0)>0
   for (UInt_t i=0; i<layerIndex_.size(); i++) {
      if (layerIndex_.at(i) == -1 && !kStarted) {
         startOffset = i;
      }
      else if (layerIndex_.at(i) != -1) {
         kStarted = true;
      }
   }

   if (startOffset>0) {
      for (Int_t i=0; i<=startOffset; i++)
         layerIndex_.at(i) = 0;
   }
}

Int_t Hits::getFirstIndexOfLayer(UInt_t layer) {
   if (layerIndex_.size() == 0) return -1;
   if (layerIndex_.size() < layer) return -1;

   return layerIndex_.at(layer);
}

Int_t Hits::getLastIndexOfLayer(UInt_t layer) {
   for (Int_t i=layer+1; i<nLayers; i++) {
      if (getFirstIndexOfLayer(i) != -1 ) {
         return getFirstIndexOfLayer(i);
      }
   }
   return GetEntriesFast(); // no hits beyond layer
}

Int_t Hits::getLastActiveLayer() {
   Int_t lastActiveLayer = 0;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i))
         continue;

      if (getLayer(i) > lastActiveLayer)
         lastActiveLayer = getLayer(i);
   }
   return lastActiveLayer;
}
/*
void Hits::makeVerticalIndexOnLayer(Int_t layer) {
   // Run this command when a new layer is to be used
   // It can only hold a single layer

   Int_t    layerIdxFrom, layerIdxTo;
   Int_t    lastY = -1;
   Int_t    startOffset = 0;
   Bool_t   kStarted = false;

   if (!layerIndex_.size())
      makeLayerIndex();

   layerIdxFrom = getFirstIndexOfLayer(layer);
   layerIdxTo = getLastIndexOfLayer(layer);

   if (verticalIndexOfLayer_.size() == 0)
      for (Int_t i=0; i<ny+1; i++) verticalIndexOfLayer_.push_back(-1);
   else
      for (Int_t i=0; i<ny+1; i++) verticalIndexOfLayer_.at(i) = -1;

   if (!GetEntriesFast()) return;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
      if (!At(i))
         continue;

      if (lastY != getY(i)) {
         verticalIndexOfLayer_.at(getY(i)) = i;
         lastY = getY(i);
      }
   }

   // set first indices to  0 if getLayer(0)>0
   for (UInt_t i=0; i<verticalIndexOfLayer_.size(); i++) {
      if (verticalIndexOfLayer_.at(i) == -1 && !kStarted) {
         startOffset = i;
      }
      else if (verticalIndexOfLayer_.at(i) != -1) {
         kStarted = true;
      }
   }

   if (startOffset>0) {
      for (Int_t i=0; i<=startOffset; i++)
         verticalIndexOfLayer_.at(i) = 0;
   }
}

Int_t Hits::getFirstIndexBeforeY(Int_t y) {
   Int_t idx = 0;

   if (y==0) return idx;

   for (Int_t i=y-1; i>=0; i--) {
      if (verticalIndexOfLayer_.size() == 0) break;
      idx = verticalIndexOfLayer_.at(i);
      if (idx>=0) break;
   }

   if (idx == -1) idx = 0;
   return idx;
}

Int_t Hits::getLastIndexAfterY(Int_t y) {
   Int_t idx = GetEntriesFast()-1;

   if (y>ny-2) return idx;

   for (Int_t i=y+2; i<ny; i++) {
      if (verticalIndexOfLayer_.size() == 0) break;
      idx = verticalIndexOfLayer_.at(i);
      if (idx>=0) break;
   }

   if (idx == -1) idx = GetEntriesFast()-1;
   return idx;
}
*/
void Hits::propagateSecondaryStatusFromTop(Int_t eventID) {
   if (eventID < 0) eventID = Last()->getEventID();
   Int_t n = 0;

   for (Int_t i=GetEntriesFast()-1; i>0; i--) {
      if (!At(i)) continue;

      if (eventID == At(i)->getEventID()) {
         At(i)->setSecondary(true);
         n++;
      }
      else break;
   }
   printf("Back-converted %d particles with eventID %d as secondaries\n", n, eventID);
}

void Hits::removeHaloAtRadius(Float_t radius) {
   Hit * hit = nullptr;
   float x,y;

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      hit = At(i);
      x = hit->getXmm();
      y = hit->getYmm();

      if (sqrt(x*x+y*y) > radius) {
         removeHitAt(i);
      }
   }
   
   Compress();
}


#endif
