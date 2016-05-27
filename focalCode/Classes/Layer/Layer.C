#include <math.h>

#include <TH2F.h>
#include <TRandom3.h>

#include "Classes/Layer/Layer.h"
#include "Classes/Hit/Hits.h"

Layer::Layer(Int_t layerNo, Bool_t frameType, Bool_t dataType) : frame2D_(Form("%d_frame2D_%i",frameType, layerNo), Form("frame2D_layer_%i", layerNo), 2*nx, 0, 2*nx, 2*ny, 0, 2*ny) {
	dataType_ = dataType;
	frameType_ = frameType;
	layerNo_ = layerNo;
}

Layer::~Layer() {
	// TODO Auto-generated destructor stub
}

Int_t Layer::diffuseLayer(TRandom3 *gRandom) {
	Int_t repeatFactor, x, y, z, randX, randY;
	Float_t EnergyFactor = 1000 / 30. * SpreadNumber;
	Float_t newSigma, eDep;

	Hits *hits = new Hits();
	Bool_t isHits = findHits(hits);
	Int_t nHits = hits->GetEntriesFast();
	
	frame2D_.Reset();

	if (kDebug) {
		cout << "Diffusing layer  " << layerNo_ << ". Number of hits = " << nHits << endl;
	}

	for (Int_t h=0; h<nHits; h++) {
		x = hits->getX(h);
		y = hits->getY(h);
		eDep = hits->getEdep(h);
		
		repeatFactor = eDep * EnergyFactor;
		newSigma = pow(repeatFactor, 0.25) / 6;

		for (Int_t j=0; j<repeatFactor; j++) {
			randX = gRandom->Gaus(x, newSigma);
			randY = gRandom->Gaus(y, newSigma);
			frame2D_.Fill(randX, randY, EnergyFactor);
		}
	}

	delete hits;
	return nHits;
}

Bool_t Layer::findHits(Hits* hits) {
	Int_t x, y, z;
	Bool_t isHits = false;
	Float_t edep;

	Int_t nBins = frame2D_.GetBin(frame2D_.GetNbinsX(), frame2D_.GetNbinsY());
	for (int i=1; i<nBins+1; i++) {
		edep = frame2D_.GetBinContent(i);
		if (edep) {
			frame2D_.GetBinXYZ(i,x,y,z);			
			hits->appendPoint(x,y,layerNo_, -1, edep);
			isHits = true;
		}
	}
	return isHits;
}
