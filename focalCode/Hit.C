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
   Clear();
}

Hit::Hit(Int_t x, Int_t y, Int_t layerNo, Int_t eventNo) {
   x_ = x;
   y_ = y;
   layerNo_ = layerNo;
	eventNo_ = eventNo;
}

void Hit::set(Int_t x, Int_t y, Int_t layerNo, Int_t eventNo) {
   x_ = x;
   y_ = y;
   layerNo_ = layerNo;
   eventNo_ = eventNo;
}

ostream &operator<< (ostream &os, Hit &h) {
	os << "(" << h.x_ << "," << h.y_ << "," << h.layerNo_ << ")";
	return os;
}

