#ifndef Hit_h
#define Hit_h

#include "TObject.h"
#include "Constants.h"

class Hit : public TObject {
  private:
      Int_t x_, y_;
	  Float_t edep_;
      Int_t layerNo_;
      Int_t eventID_;

  public:
      Hit();
      Hit(Int_t x, Int_t y, Int_t layer = -1, Int_t event = -1, Float_t edep_ = 0);
      Hit(Hit* hit);
      virtual ~Hit(); 

      Int_t getX() { return x_; }
      Int_t getY() { return y_; }
      Int_t getLayer() { return layerNo_; }
      Int_t getEventID() { return eventID_; }
      Float_t getEdep() { return edep_; }

      Double_t getXmm() { return x_ * dx; }
      Double_t getYmm() { return y_ * dy; }

      void set(Int_t x, Int_t y, Int_t layerNo = -1, Int_t eventNo = -1, Float_t edep = 0);
      void set(Hit* hit);
      void setEventID(Int_t event) { eventID_ = event; }
      void setEdep(Float_t edep) { edep_ = edep; }

      friend ostream &operator<< (ostream &os, Hit &h);

      ClassDef(Hit,1);
};
#endif
