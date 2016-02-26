#include <iostream>
#include "Hit.h"

// ClassImp(Hit)

Hit::Hit() {
	x_ = -1;
	y_ = -1;
	layerNo_ = -1;
	eventNo_ = -1;
	edep_ = 0;
}

Hit::~Hit() {
  // Destructor
}

Hit::Hit(Int_t x, Int_t y, Int_t layerNo, Int_t eventNo, Float_t edep) {
  x_ = x;
  y_ = y;
  layerNo_ = layerNo;
  eventNo_ = eventNo;
  edep_ = edep;
}

Hit::Hit(Hit* hit) {
  x_ = hit->getX();
  y_ = hit->getY();
  layerNo_ = hit->getLayer();
  eventNo_ = hit->getEvent();
  edep_ = hit->getEdep();
}

void Hit::set(Int_t x, Int_t y, Int_t layerNo, Int_t eventNo, Float_t edep) {
  x_ = x;
  y_ = y;
  layerNo_ = layerNo;
  eventNo_ = eventNo;
  edep_ = edep;
}

void Hit::set(Hit* hit) {
  x_ = hit->getX();
  y_ = hit->getY();
  layerNo_ = hit->getLayer();
  eventNo_ = hit->getEvent();
  edep_ = hit->getEdep();
}

ostream &operator<< (ostream &os, Hit &h) {
	os << "(" << h.x_ << "," << h.y_ << "," << h.layerNo_ << ")";
	return os;
}

