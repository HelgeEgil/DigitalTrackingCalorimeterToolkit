#ifndef Layer_h
#define Layer_h

#include "Constants.h"
#include "Clusters.h"
#include "Hits.h"
#include <TH2F.h>

class Layer : public TObject {
	private:
		TH2F frame2D_;
		Int_t layerNo_;
		Bool_t frameType_;
		Bool_t dataType_;

	public:
		// frameType: kCalorimeter og kTracker
		// dataType: kMC or kData

		Layer(Int_t layerNo, Bool_t frameType, Bool_t dataType);
		virtual ~Layer();

		void diffuseLayer();
		void oldDiffuseLayer();
		virtual Bool_t findHits(Hits *hits);

		virtual void Reset() { frame2D_.Reset(); }

		virtual void Fill(Float_t x, Float_t y, Float_t val = 1) { frame2D_.Fill(x,y,val); }
		virtual TH2F * getTH2F() { return (TH2F*) &frame2D_; }

		ClassDef(Layer,1);
};

#endif /* Layer_h */
