#include <TH2F.h>
#include <TRandom3.h>

#include "Classes/CalorimeterFrame/CalorimeterFrame.h"
#include "Classes/Layer/Layer.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"

using namespace std;
using namespace DTC;

CalorimeterFrame::CalorimeterFrame() : layers_("DTC::Layer", nLayers) {
  layers_.ExpandCreate(nLayers);
   for (Int_t i=0; i<nLayers; i++) {
      At(i)->SetProperties(i, kCalorimeter, kMC);
      // the memory should be contiguous, the object pointers should be
      // in an ascending order with a gap of the size of the Layer class
      // but it is not in the first elements and if placement new is used
      // to create the objects one by one, there is a double delete error
      // when deleting the CalorimeterFrame
      // Needs more investigation
      //std::cout << "CalorimeterFrame::CalorimeterFrame(): Layer " << At(i)->getLayer() << ": " << At(i) << std::endl;
   }
}

CalorimeterFrame::~CalorimeterFrame() {
   layers_.Delete();
}

void CalorimeterFrame::Clear(Option_t *) {
   layers_.Clear("C");
}

void CalorimeterFrame::Reset() {
   for (Int_t layer=0; layer<nLayers; layer++) {
      At(layer)->Reset();
   }
}

Hits * CalorimeterFrame::findHits(Int_t eventID) {
   Hits *hits = new Hits();

   Int_t nNoHits = 0;
   Bool_t isHits = false;
   for (Int_t layer=0; layer<nLayers; layer++) {
      isHits = At(layer)->findHits(hits);

      if (!isHits) nNoHits++;
      else nNoHits = 0;
      if (nNoHits>2) break;
   }
   
   if (eventID > 0) {
      for (Int_t i=0; i<hits->GetEntriesFast(); i++) {
         hits->At(i)->setEventID(eventID);
      }
   }

   hits->makeLayerIndex();
   return hits;
}

void CalorimeterFrame::diffuseFrame(TRandom3 *gRandom) {
   Int_t nHitsInLayer = 1;
   
   if (kUseRefinedClustering) {
      for (Int_t layer=0; layer<nLayers; layer++) {
         if (nHitsInLayer) {
            nHitsInLayer = At(layer)->refinedDiffuseLayer(gRandom);
         }
         else break;
      }
   }

   else {
      for (Int_t layer=0; layer<nLayers; layer++) {
         if (nHitsInLayer) {
            nHitsInLayer = At(layer)->diffuseLayer(gRandom);
         }
         else break;
      }
   }
}

Float_t CalorimeterFrame::getOccupancyLastLayer() {
   Float_t occupancy = 0;
   Float_t lastOccupancy = 0;
   for (Int_t layer=0; layer<nLayers; layer++) {
      occupancy = At(layer)->getOccupancy();

      if (occupancy>lastOccupancy*0.9) { lastOccupancy = occupancy; }
      else { break; }
   }

   return lastOccupancy;
}
