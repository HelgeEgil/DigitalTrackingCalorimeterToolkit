#ifndef Layer_cxx
#define Layer_cxx

#include <math.h>

#include <TH2F.h>
#include <TRandom3.h>

#include "GlobalConstants/MaterialConstants.h"
#include "Classes/Layer/Layer.h"
#include "Classes/Hit/Hits.h"

DTC::Layer::Layer(Int_t layerNo, Bool_t frameType, Bool_t dataType) : frame2D_(Form("%d_frame2D_%i",frameType, layerNo), Form("frame2D_layer_%i", layerNo), nx, 0, nx, ny, 0, ny) {
   dataType_ = dataType;
   frameType_ = frameType;
   layerNo_ = layerNo;
}

DTC::Layer::~Layer() {
}

Int_t DTC::Layer::diffuseLayer(TRandom3 *gRandom) {
   if (kUseRefinedClustering) {
      return refinedDiffuseLayer(gRandom);
   }
   else {
      return oldDiffuseLayer(gRandom);
   }
}


Int_t DTC::Layer::oldDiffuseLayer(TRandom3 *gRandom) {
   Int_t repeatFactor, x, y, z, randX, randY;
   Float_t EnergyFactor = 2.67;
   Float_t newSigma, eDep;

   Hits *hits = new Hits();
   Bool_t isHits = findHits(hits);
   Int_t nHits = hits->GetEntriesFast();
   
   frame2D_.Reset();

   showDebug("Diffusing layer  " << layerNo_ << ". Number of hits = " << nHits << endl);

   for (Int_t h=0; h<nHits; h++) {
      x = hits->getX(h);
      y = hits->getY(h);
      eDep = hits->getEdep(h);
      repeatFactor = eDep * EnergyFactor;
      newSigma = pow(repeatFactor, 0.35) / 6;

      for (Int_t j=0; j<repeatFactor; j++) {
         randX = gRandom->Gaus(x, newSigma);
         randY = gRandom->Gaus(y, newSigma);
         frame2D_.Fill(randX, randY, EnergyFactor);
      }
   }

   delete hits;
   return nHits;
}


Int_t DTC::Layer::refinedDiffuseLayer(TRandom3 *gRandom) {
   Int_t repeatFactor, x, y, z, randX, randY;
   Float_t EnergyFactor = 2.5;
   Float_t newSigma, eDep;

   Hits *hits = new Hits();
   Bool_t isHits = findHits(hits);
   Int_t nHits = hits->GetEntriesFast();
   frame2D_.Reset();

   showDebug("Diffusing layer " << layerNo_ << " (refined). Number of hits = " << nHits << endl);

   for (Int_t h=0; h<nHits; h++) {
      x = hits->getX(h);
      y = hits->getY(h);
      eDep = hits->getEdep(h); // per um
//      showDebug("Hit at " << x << ", " << y << ": " << eDep << "keV/um\n");
      repeatFactor = eDep * EnergyFactor;
      newSigma = pow(repeatFactor, 0.41) / 6;

      for (Int_t j=0; j<repeatFactor; j++) {
         randX = gRandom->Gaus(x, newSigma);
         randY = gRandom->Gaus(y, newSigma);
//         frame2D_.Fill(randX, randY, EnergyFactor);
         frame2D_.SetBinContent(randX, randY, 1);
      }
   }

   delete hits;
   return nHits;
}



Bool_t DTC::Layer::findHits(Hits* hits) {
   Int_t    x, y, z, nBins;
   Bool_t   isHits = false;
   Float_t  edep;
   Int_t    dummyEventID = -1;

   nBins = frame2D_.GetBin(frame2D_.GetNbinsX(), frame2D_.GetNbinsY());

   for (Int_t i=1; i<nBins+1; i++) {
      edep = frame2D_.GetBinContent(i);

      if (edep) {
         frame2D_.GetBinXYZ(i,x,y,z);
         hits->appendPoint(x,y,layerNo_, edep, dummyEventID);
         isHits = true;
      }
   }

   return isHits;
}

Float_t DTC::Layer::getOccupancy() {
   Int_t    numberOfActivatedPixels = 0;
   Int_t    nBins = frame2D_.GetBin(frame2D_.GetNbinsX(), frame2D_.GetNbinsY());
   Float_t  occupancy;

   for (int i=1; i<nBins+1; i++) {
      if (frame2D_.GetBinContent(i)) numberOfActivatedPixels++;
   }

   occupancy = (Float_t) numberOfActivatedPixels / nBins;

   return occupancy;
}

#endif
