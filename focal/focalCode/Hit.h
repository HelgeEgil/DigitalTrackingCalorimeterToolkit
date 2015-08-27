#ifndef Hit_h
#define Hit_h

#include "TObject.h"
#include "Constants.h"

class Hit : public TObject {
   private:
      Int_t x_, y_;
      Int_t layerNo_;
		Int_t eventNo_;

   public:
      Hit();
      Hit(Int_t x, Int_t y, Int_t layer = -1, Int_t event = -1);
      virtual ~Hit(); 

      Int_t getX() { return x_; }
      Int_t getY() { return y_; }
      Int_t getLayer() { return layerNo_; }
		Int_t getEvent() { return eventNo_; }

		Int_t getXmm() { return x_ * dx; }
		Int_t getYmm() { return y_ * dy; }

      void set(Int_t x, Int_t y, Int_t layerNo = -1, Int_t eventNo = -1);
		void setEvent(Int_t event) { eventNo_ = event; }

		friend ostream &operator<< (ostream &os, Hit &h);

   ClassDef(Hit,1);
};
#endif
