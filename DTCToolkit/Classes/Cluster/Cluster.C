#ifndef Cluster_c
#define Cluster_c

#include <iostream>

#include <math.h>

#include "Classes/Cluster/Cluster.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"

using namespace DTC;

Cluster::Cluster() {
   x_ = -1;
   y_ = -1;
   layerNo_ = -1;
   clusterSize_ = -1;
   eventID_ = -1;
   isUsed_ = false;
   isSecondary_ = false;
   PDG_ = 0;
}

Cluster::Cluster(Cluster* cluster) {
   x_ = cluster->getX();
   y_ = cluster->getY();
   layerNo_ = cluster->getLayer();
   clusterSize_ = cluster->getSize();
   eventID_ = cluster->getEventID();
   isUsed_ = cluster->isUsed();
   isSecondary_ = cluster->isSecondary();
   PDG_ = cluster->getPDG();
}

Cluster::Cluster(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID, Bool_t isSecondary, Int_t PDG) {
   x_ = x;
   y_ = y;
   layerNo_ = layer;
   clusterSize_ = size;
   eventID_ = eventID;
   isUsed_ = false;
   isSecondary_ = isSecondary;
   PDG_ = PDG;
}

Cluster::~Cluster() {
   // Destructor
}

Int_t Cluster::Compare(const TObject *obj) const {
//   if (layerNo_ == ((Cluster*) obj)->getLayer()) return 0;
//   else if (layerNo_ < ((Cluster*) obj)->getLayer()) return -1;
   if (eventID_ == ((Cluster*) obj)->getEventID() && layerNo_ == ((Cluster*) obj)->getLayer()) return 0;
   else if (eventID_ < ((Cluster*) obj)->getEventID()) return -1;
   else if (eventID_ == ((Cluster*) obj)->getEventID() && layerNo_ < ((Cluster*) obj)->getLayer()) return -1;
   else return 1;
}

Double_t Cluster::getLayermm() {
   return getLayerPositionmm(layerNo_);
}

Int_t Cluster::getChip() {
   Int_t    layer = getLayer();
   Int_t    chipIdx = 0;
   Float_t  y = getYmm();
   Float_t  x = getXmm();
   
    if   ( x>0 &&  y>0) chipIdx = 0;
   else if  (x<=0 &&  y>0) chipIdx = 3;
   else if  (x<=0 && y<=0) chipIdx = 2;
   else if  ( x>0 && y<=0) chipIdx = 1;

   chipIdx += layer*4;

    return chipIdx;
}

Float_t Cluster::getRadiusmm() {
  Float_t A = getSize();
  Float_t r = sqrt( A / 3.141592653);
  Float_t r_mm = r * dx;
  
  return r_mm;
}  

Float_t Cluster::getDepositedEnergy(Bool_t correctSensitivity) {
   // An inverse fit to the charge diffusion algorithm defined in Layer.C
   // Returns the local dEdx value in keV.

   Float_t  edep = 0;
   Int_t    n = clusterSize_;

   edep = getEdepFromCS(n);

   if (correctSensitivity) {
      edep *= getChipCalibrationFactor(getChip());
   }

   return edep; // layer is 14 modelled as um thick
}

Float_t Cluster::getCalibratedSize() {
   Int_t n = clusterSize_;

   Float_t calibration = getChipCalibrationFactor(getChip());
   Float_t edep = getEdepFromCS(n);

   if (calibration != 0) { 
      edep *= calibration;
   }

   Float_t n_corrected = getCSFromEdep(edep);

   return n_corrected;
}

Float_t Cluster::getDepositedEnergyError(Bool_t correctSensitivity) {
   // sigma_E = sigma_N (dE / dN)
   // TODO: FIX THIS
   
   Float_t  sn = sqrt(clusterSize_);
   Int_t    n = clusterSize_;
   Float_t  diff_n;
  
   diff_n = 3.883 - 2 * 0.0124 * n + 3 * 0.00114 * pow(n,2) - 4 * 1.420e-5 * pow(n,3);

   if (correctSensitivity) {
      diff_n *= getChipCalibrationFactor(getChip());
   }

   return sn * diff_n / 14.; // layer is modelled as 14 um thick
}

void Cluster::set(Cluster* copyCluster) {
   x_ = copyCluster->getX();
   y_ = copyCluster->getY();
   layerNo_ = copyCluster->getLayer();
   clusterSize_ = copyCluster->getSize();
   eventID_ = copyCluster->getEventID();
   isUsed_ = copyCluster->isUsed();
   isSecondary_ = copyCluster->isSecondary(); 
   PDG_ = copyCluster->getPDG();
}

void Cluster::set(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID, Bool_t isSecondary, Int_t PDG) {
   x_ = x;
   y_ = y;
   if (layer>=0) layerNo_ = layer;
   if(size>=0) clusterSize_ = size;
   if (eventID>=0) eventID_ = eventID;
   isSecondary_ = isSecondary;
   PDG_ = PDG;

}

ostream& operator<< (ostream &os, Cluster& c) {
//   os << "(" << c.getXmm() << ", " << c.getYmm() << ", " << c.getLayermm() << ", EID " << c.getEventID() << ", CS " << c.getSize() << ")";
   os << "(" << c.getXmm() << ", " << c.getYmm() << ", " << c.getLayer() << ", EID " << c.getEventID() << ", CS " << c.getSize() << ", is2nd " << c.isSecondary() << ", PDG " << c.getPDG() << ")";
   return os;
}

#endif
