/*
 * Layer.h
 *
 *  Created on: Aug 11, 2015
 *      Author: local1
 */

#ifndef Layer_h
#define Layer_h

#include "Constants.h"

class Layer : public TH2F {
	private:
		TH2F frame2D_;
		const Int_t layerNo_;
		const Bool_t frameType_;
		const Bool_t dataType_;

	public:
		// frameType: kCalorimeter og kTracker
		// dataType: kMC or kData

		Layer(Int_t layerNo, Bool_t frameType, Bool_t dataType);
		virtual ~Layer();

		void diffuseLayer();
		virtual void findHits(Hits *hits);
		Clusters *findClustersFromHits();

		virtual void Reset() { frame2D_.Reset(); }

		virtual void Fill(Float_t x, Float_t y, Float_t val = 1) { frame2D_.Fill(x,y,val); }
		TH2F * getTH2F() { return (TH2F*) &frame2D_; }

};

#endif /* Layer_h */
