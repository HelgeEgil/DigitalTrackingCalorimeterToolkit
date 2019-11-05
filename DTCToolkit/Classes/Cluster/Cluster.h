#ifndef Cluster_h
#define Cluster_h

#include <iostream>

#include <TObject.h>
#include <math.h>

#include "GlobalConstants/Constants.h"

namespace DTC {
class Cluster : public TObject {
private:
   Float_t  x_, y_;
   Int_t    layerNo_;
   Int_t    clusterSize_;
   Int_t    eventID_;
   Bool_t   isUsed_;
   Bool_t   isSecondary_;
   Int_t    PDG_;

public:
   Cluster(); 
   Cluster(Cluster *cluster);
   Cluster(Float_t x, Float_t y, Int_t layer = -1, Int_t size = -1, Int_t eventID = -1, Bool_t isSecondary = false, Int_t PDG = 0);
   virtual ~Cluster();
   
   // Sorting
   Int_t    Compare(const TObject *obj) const;
   Bool_t   IsSortable() const { return kTRUE; }

   // inline getters and setters
   Float_t  getX()            { return x_; }
   Float_t  getY()            { return y_; }
   Int_t    getLayer()        { return layerNo_; }
   Int_t    getSize()         { return clusterSize_; }
   Int_t    getEventID()      { return eventID_; }
   Int_t    getError()        { return sqrt(clusterSize_); }
   Int_t    getPDG()          { return PDG_; }
   Double_t getXmm()          { return (x_ - nx/2) * dx; } // from -nx*dx/2 to nx*dx/2
   Double_t getYmm()          { return (y_ - ny/2) * dy; } // from -ny*dy to ny*dy
   Bool_t   isSecondary()     { return isSecondary_; }
   Bool_t   isUsed()          { return isUsed_; }

   void     setX(Float_t x)      { x_ = x; }
   void     setY(Float_t y)      { y_ = y; }
   void     setLayer(Int_t l)    { layerNo_ = l; }
   void     setSize(Int_t size)  { clusterSize_ = size; }
   void     setEventID(Int_t id) { eventID_ = id; }
   void     setSecondary(Bool_t s){ isSecondary_ = s; }
   void     setPDG(Int_t PDG)    { PDG_ = PDG; }
   void     setXmm(Float_t x)    { x_ = x / dx + nx/2; }
   void     setYmm(Float_t y)    { y_ = y / dy + ny/2; }
   void     markUsed()           { isUsed_ = true; }
   void     markUnused()         { isUsed_ = false; }
   
   // In .C
   Double_t getLayermm(); // take into account placement of sensor layer, depends on y
   Float_t  getCalibratedSize();
   Int_t    getChip();
   Float_t  getRadiusmm(); // cluster size as circle radius [mm]
   Float_t  getDepositedEnergy(Bool_t checkResistivity = false);
   Float_t  getDepositedEnergyError(Bool_t checkResistivity = false);
   void     set(Float_t x, Float_t y, Int_t layer = -1, Int_t size = -1, Int_t eventID = -1, Bool_t isSecondary = false, Int_t PDG = 0);
   void     set(Cluster *copyCluster); // copy properties, not pointer

   ClassDef(Cluster,5)
};
}

ostream& operator<<(ostream &os, DTC::Cluster &c);

#endif
