#ifndef Hit_cxx
#define Hit_cxx

#include <iostream>

#include "Classes/Hit/Hit.h"

using namespace DTC;

Hit::Hit() {
   x_ = -1;
   y_ = -1;
   layerNo_ = -1;
   eventID_ = -1;
   edep_ = 0;
}


Hit::Hit(Int_t x, Int_t y, Int_t layerNo, Int_t eventID, Float_t edep) {
   x_ = x;
   y_ = y;
   layerNo_ = layerNo;
   eventID_ = eventID;
   edep_ = edep;
}

Hit::Hit(Hit* hit) {
   x_ = hit->getX();
   y_ = hit->getY();
   layerNo_ = hit->getLayer();
   eventID_ = hit->getEventID();
   edep_ = hit->getEdep();
}

Hit::~Hit() {
}

Int_t Hit::getChip() {
   Int_t x = getX();
   Int_t y = getY();
   Int_t layer = getLayer();
   Int_t chipIdx = 0;

   if       (x >  nx/2 && y >  ny/2) chipIdx = 0;
   else if  (x <= nx/2 && y >  ny/2) chipIdx = 3;
   else if  (x <= nx/2 && y <= ny/2)   chipIdx = 2;
   else if  (x >  nx/2 && y <= ny/2)   chipIdx = 1;

   chipIdx += layer*4;

   return chipIdx;
}

void Hit::set(Int_t x, Int_t y, Int_t layerNo, Int_t eventID, Float_t edep) {
   x_ = x;
   y_ = y;
   layerNo_ = layerNo;
   eventID_ = eventID;
   edep_ = edep;
}

void Hit::set(Hit* hit) {
   x_ = hit->getX();
   y_ = hit->getY();
   layerNo_ = hit->getLayer();
   eventID_ = hit->getEventID();
   edep_ = hit->getEdep();
}

ostream &operator<< (ostream &os, Hit &h) {
   os << "(" << h.getX() << "," << h.getY() << "," << h.getLayer() << ")";
   return os;
}

#endif
