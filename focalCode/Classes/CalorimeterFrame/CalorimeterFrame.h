#ifndef FrameCalorimeter_h
#define FrameCalorimeter_h

#include <vector>
#include <TClonesArray.h>
#include <TH2F.h>

#include "GlobalConstants/Constants.h"
#include "Classes/Layer/Layer.h"

using namespace std;

class Hits;

// A Frame contains many layers

class CalorimeterFrame : public TObject {

   private:
      TClonesArray calorimeterFrame_;

   public:
      CalorimeterFrame();
      virtual ~CalorimeterFrame();

      virtual void Clear(Option_t * = "");

      virtual Layer * At(Int_t i) { return (Layer*) calorimeterFrame_.At(i); }
      virtual TH2F * getTH2F(Int_t i) { return (TH2F*) At(i)->getTH2F(); }
      virtual void fillAt(Int_t i, Float_t x, Float_t y, Float_t val = 1) { getTH2F(i)->Fill(x, y, val); }

      virtual Hits *findHits(Int_t eventID = -1);
      virtual void diffuseFrame(TRandom3 *gRandom);
      virtual void Reset();

   ClassDef(CalorimeterFrame,2);
};
#endif
