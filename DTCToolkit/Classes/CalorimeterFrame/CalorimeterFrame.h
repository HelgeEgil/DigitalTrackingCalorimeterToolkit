#ifndef FrameCalorimeter_h
#define FrameCalorimeter_h

#include <vector>
#include <TClonesArray.h>
#include <TH2F.h>

#include "GlobalConstants/Constants.h"
#include "Classes/Layer/Layer.h"

using namespace std;

namespace DTC {
class Hits;

// A Frame contains many layers

class CalorimeterFrame : public TObject {

   private:
      TClonesArray calorimeterFrame_;

   public:
      CalorimeterFrame();
      virtual ~CalorimeterFrame();

      // ROOT functions
      virtual void   Clear(Option_t * = "");
      virtual void   Reset();
      virtual Layer *At(Int_t i) { return (DTC::Layer*) calorimeterFrame_.At(i); }

      // inline functions
      TH2F         * getTH2F(Int_t i) { return (TH2F*) At(i)->getTH2F(); }
      void           fillAt(Int_t i, Float_t x, Float_t y, Float_t val = 1) { getTH2F(i)->Fill(x, y, val); }

      // bigger functions in .C
      virtual Hits * findHits(Int_t eventID = -1);
      virtual void   diffuseFrame(TRandom3 *gRandom);
      Float_t        getOccupancyLastLayer();

   ClassDef(CalorimeterFrame,2)
};
}
#endif
