#ifndef Hit_h
#define Hit_h

#include <TObject.h>

#include "GlobalConstants/Constants.h"

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

      Int_t		getX()			{ return x_; }
      Int_t		getY()			{ return y_; }
      Double_t	getXmm()			{ return (x_ - nx) * dx; }
      Double_t	getYmm()			{ return (y_ - ny) * dy; }
      Int_t		getLayer()		{ return layerNo_; }
      Int_t		getEventID()	{ return eventID_; }
      Float_t	getEdep()		{ return edep_; }

      void		setEventID(Int_t event)	{ eventID_ = event; }
      void		setEdep(Float_t edep)		{ edep_ = edep; }
      void		set(Int_t x, Int_t y, Int_t layerNo = -1, Int_t eventNo = -1, Float_t edep = 0);
      void		set(Hit* hit);

      friend ostream &operator<< (ostream &os, Hit &h);

      ClassDef(Hit,2)
};
#endif
