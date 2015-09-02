#ifndef FrameCalorimeter_h
#define FrameCalorimeter_h

#include "TClonesArray.h"
#include "Constants.h"
#include "Hits.h"
#include "Layer.h"

#include <vector>

using namespace std;

// A Frame contains many layers

class CalorimeterFrame : public TObject {

   private:
      TClonesArray calorimeterFrame_;

   public:
      CalorimeterFrame();
      virtual ~CalorimeterFrame();

      virtual Layer * At(Int_t i) { return (Layer*) calorimeterFrame_.At(i); }
      virtual TH2F * getTH2F(Int_t i) { return (TH2F*) At(i)->getTH2F(); }
      virtual void fillAt(Int_t i, Float_t x, Float_t y, Float_t val = 1) { getTH2F(i)->Fill(x, y, val); }

      virtual Hits *findHits();
      virtual void diffuseFrame();
      virtual void Reset();

   ClassDef(CalorimeterFrame,1);
};
#endif
