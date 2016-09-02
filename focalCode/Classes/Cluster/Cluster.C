#include <iostream>

#include <math.h>

#include "Classes/Cluster/Cluster.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"

Cluster::Cluster() {
	x_ = -1;
	y_ = -1;
	layerNo_ = -1;
	clusterSize_ -1;
	eventID_ = -1;
	isUsed_ = false;
}

Cluster::Cluster(Cluster* cluster) {
   x_ = cluster->getX();
   y_ = cluster->getY();
   layerNo_ = cluster->getLayer();
   clusterSize_ = cluster->getSize();
   eventID_ = -1;
	isUsed_ = false;
}

Cluster::Cluster(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID) {
   x_ = x;
   y_ = y;
   layerNo_ = layer;
   clusterSize_ = size;
   eventID_ = eventID;
	isUsed_ = false;
}

Cluster::~Cluster() {
   // Destructor
}

Double_t Cluster::getLayermm() {
	Double_t z = 0;
	Float_t discriminator = 0;

	if (kDataType == kMC) {
		discriminator = getYmm();
	}
	else {
		discriminator = getXmm();
	}

	if (layerNo_ > 0) {
		
		if (discriminator > 0)	{ z = firstUpperLayerZ + layerNo_*dz; }
		else							{ z = firstLowerLayerZ + layerNo_*dz; }
	}

	return z;
}

Int_t Cluster::getChip() {
    Int_t   layer = getLayer();
    Int_t   chipIdx = 0;
	Float_t	y = getYmm();
	Float_t	x = getXmm();
	
    if  	( x>0 &&  y>0) chipIdx = 0;
	else if	(x<=0 &&  y>0) chipIdx = 3;
	else if	(x<=0 && y<=0) chipIdx = 2;
	else if	( x>0 && y<=0) chipIdx = 1;

	chipIdx += layer*4;

    return chipIdx;
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
	
	return 2.854 * clusterSize_ + 0.06641 * pow(clusterSize_, 2) - 0.001156 * pow(clusterSize_, 3) + 8.246e-6 * pow(clusterSize_, 4);
}

Float_t Cluster::getDepositedEnergyError() {
	// sigma_E = sigma_N (dE / dN)
	
	Float_t sigma_n = sqrt(clusterSize_);
	Float_t diff_n = 2.854 + 2 * 0.06641 * clusterSize_ - 3 * 0.001156 * pow(clusterSize_, 2) + 4 * 8.246e-6 * pow(clusterSize_, 3);

	return sigma_n * diff_n;
}

void Cluster::set(Cluster* copyCluster) {
   x_ = copyCluster->getX();
   y_ = copyCluster->getY();
   layerNo_ = copyCluster->getLayer();
   clusterSize_ = copyCluster->getSize();
   eventID_ = copyCluster->getEventID();
	isUsed_ = copyCluster->isUsed();	
}

void Cluster::set(Float_t x, Float_t y, Int_t layer, Int_t size) {
   x_ = x;
   y_ = y;
   if (layer>=0) layerNo_ = layer;
   if(size>=0) clusterSize_ = size;
}

ostream& operator<< (ostream &os, Cluster& c) {
   os << "(" << c.x_ << ", " << c.y_ << ", " << c.layerNo_ << ", EID " << c.eventID_ << ", CS " << c.clusterSize_ << ")";
   return os;
}

