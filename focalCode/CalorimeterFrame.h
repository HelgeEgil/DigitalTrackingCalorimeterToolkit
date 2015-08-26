#ifndef FrameCalorimeter_h
#define FrameCalorimeter_h
#include "TClonesArray.h"

#include <vector>
using namespace std;

// A Frame contains many layers

class CalorimeterFrame : public TClonesArray {

   private:
      TClonesArray frameCalorimeter3D_;

   public:
      CalorimeterFrame() : frameCalorimeter3D_("Layer", nLayers) { }

      Layer * At(Int_t i) { return (Layer*) frameCalorimeter3D_.At(i); }
      TH2F * getTH2F(Int_t i) { return (TH2F*) At(i)->getTH2F(); }
      virtual void fillAt(Int_t i, Float_t x, Float_t y, Float_t val) { getTH2F(i)->Fill(x, y, val); }

      virtual Hits *findHits();
      virtual void diffuseFrame();
      virtual void Reset();

   ClassDef(CalorimeterFrame,1);
};
#endif
