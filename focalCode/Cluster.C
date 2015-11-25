#include <iostream>
#include "Cluster.h"
#include <math.h>

// ClassImp(Cluster)

Cluster::~Cluster() {
   // Destructor
   Clear();
}

Cluster::Cluster(Cluster* cluster) {
   x_ = cluster->getX();
   y_ = cluster->getY();
   layerNo_ = cluster->getLayer();
   clusterSize_ = cluster->getSize();
}

Cluster::Cluster(Float_t x, Float_t y, Int_t layer, Int_t size) {
   x_ = x;
   y_ = y;
   layerNo_ = layer;
   clusterSize_ = size;
}

Double_t Cluster::getLayermm() {
	Double_t z = 0;

	if (layerNo_==-4) z = -67.68;
	if (layerNo_==-3) z = -63.84;
	if (layerNo_==-2) z = -13.84;
	if (layerNo_==-1) z = -10;

	if (layerNo_>0) z = layerNo_*dz;

	return z;
}

Float_t Cluster::getRadiusmm() {
  Float_t A = getSize();
  Float_t r = sqrt( A / 3.141592653);
  Float_t r_mm = r * dx;
  
  return r_mm;
}  

Float_t Cluster::getDepositedEnergy() {
	// An inverse fit to the charge diffusion algorithm defined in Layer.C
	// Returns the local dEdx value in keV.
	return -13.6 + 7.524*clusterSize_ + 0.248 * pow(clusterSize_,2);
}


void Cluster::set(Cluster* copyCluster) {
   x_ = copyCluster->getX();
   y_ = copyCluster->getY();
   layerNo_ = copyCluster->getLayer();
   clusterSize_ = copyCluster->getSize();
}

void Cluster::set(Float_t x, Float_t y, Int_t layer, Int_t size) {
   	x_ = x;
      y_ = y;
      if (layer>=0) layerNo_ = layer;
      if(size>=0) clusterSize_ = size;
   }

   ostream& operator<< (ostream &os, Cluster& c) {
      os << "(" << c.x_ << "," << c.y_ << "," << c.layerNo_ << ")";
      return os;
}

