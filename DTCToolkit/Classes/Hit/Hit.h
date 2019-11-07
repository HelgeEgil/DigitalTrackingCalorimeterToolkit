#ifndef Hit_h
#define Hit_h

#include <TObject.h>

#include "GlobalConstants/Constants.h"

namespace DTC {
class Hit : public TObject {
  private:
      Int_t    x_, y_;
      Float_t  edep_;
      Int_t    layerNo_;
      Int_t    eventID_;
      Bool_t   isSecondary_;
      Int_t    pdg_;

  public:
      Hit();
      Hit(Int_t x, Int_t y, Int_t layer = -1, Float_t edep_ = 0, Int_t eventID = 0, Bool_t isSecondary_ = false, Int_t PDG = 0);
      Hit(Hit* hit);
      virtual ~Hit(); 

      Int_t    getX()         { return x_; }
      Int_t    getY()         { return y_; }
      Double_t getXmm()       { return (x_ - nx/2) * dx; }
      Double_t getYmm()       { return (y_ - ny/2) * dy; }
      Int_t    getLayer()     { return layerNo_; }
      Int_t    getEventID()   { return eventID_; }
      Float_t  getEdep()      { return edep_; }
      Int_t    getPDG()       { return pdg_; }
      Bool_t   isSecondary()  { return isSecondary_; }
      Int_t    getChip();

      void     setEventID(Int_t event)          { eventID_ = event; }
      void     setEdep(Float_t edep)            { edep_ = edep; }
      void     setSecondary(Bool_t isSecondary) { isSecondary_ = isSecondary; }
      void     setPDG(Int_t PDG)                { pdg_ = PDG; }
      void     set(Int_t x, Int_t y, Int_t layerNo = -1, Float_t edep = 0, Int_t eventID = -1, Bool_t isSecondary = false, Int_t PDG = 0);
      void     set(Hit* hit);
      void     Clear(Option_t *);
      Int_t    Compare(const TObject *obj) const;
      Bool_t   IsSortable() const { return kTRUE; }

      ClassDef(Hit,3)
};
}

ostream &operator<< (ostream &os, DTC::Hit &h);

#endif
