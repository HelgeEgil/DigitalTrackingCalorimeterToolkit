#include <iostream>
#include "Hit.h"

// ClassImp(Hit)

Hit::Hit() {
	x_ = -1;
	y_ = -1;
	layerNo_ = -1;
	eventNo_ = -1;
}

Hit::~Hit() {
  // Destructor
}

Hit::Hit(Int_t x, Int_t y, Int_t layerNo, Int_t eventNo) {
  x_ = x;
  y_ = y;
  layerNo_ = layerNo;
  eventNo_ = eventNo;
}

Hit::Hit(Hit* hit) {
  x_ = hit->getX();
  y_ = hit->getY();
  layerNo_ = hit->getLayer();
  eventNo_ = hit->getEvent();
}

void Hit::set(Int_t x, Int_t y, Int_t layerNo, Int_t eventNo) {
  x_ = x;
  y_ = y;
  layerNo_ = layerNo;
  eventNo_ = eventNo;
}

void Hit::set(Hit* hit) {
  x_ = hit->getX();
  y_ = hit->getY();
  layerNo_ = hit->getLayer();
  eventNo_ = hit->getEvent();
}

ostream &operator<< (ostream &os, Hit &h) {
	os << "(" << h.x_ << "," << h.y_ << "," << h.layerNo_ << ")";
	return os;
}

