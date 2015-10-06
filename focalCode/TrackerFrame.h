#ifndef TrackerFrame_h
#define TrackerFrame_h

#include "TClonesArray.h"
#include "Constants.h"
#include "Layer.h"
#include <vector>

using namespace std;

// A Frame contains many layers

class TrackerFrame : public TObject {

   private:
      TClonesArray trackerFrame_;

   public:
      TrackerFrame();
      virtual ~TrackerFrame();

      virtual Layer * At(Int_t i) { return (Layer*) trackerFrame_.At(i); }
      virtual TH2F * getTH2F(Int_t i) { return (TH2F*) At(i)->getTH2F(); }
      virtual void fillAt(Int_t i, Float_t x, Float_t y, Float_t val = 1) { getTH2F(i)->Fill(x, y, val); }

      virtual Hits *findHits();
      virtual void diffuseFrame();
      virtual void Reset();

   ClassDef(TrackerFrame,1);
};
#endif
