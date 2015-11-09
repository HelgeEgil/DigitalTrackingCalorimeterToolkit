#include "Layer.h"
#include <TH2F.h>
#include <TRandom3.h>

Layer::Layer(Int_t layerNo, Bool_t frameType, Bool_t dataType) : frame2D_(Form("frame2D_%i",layerNo), Form("frame2D_layer_%i", layerNo), 2*nx, 0, 2*nx, 2*ny, 0, 2*ny) {
	dataType_ = dataType;
	frameType_ = frameType;
	layerNo_ = layerNo;
}

Layer::~Layer() {
	// TODO Auto-generated destructor stub
}

void Layer::diffuseLayer() {
	Int_t repeatFactor, x, y, z, randX, randY;
   Float_t EnergyFactor = 1000 / 30. * SpreadNumber;
	Float_t newSigma, eDep;
   gRandom = new TRandom3(0);

	TH2F *frame2DCopy = (TH2F*) frame2D_.Clone();
	frame2DCopy->SetName("frame2DCopy");
	frame2D_.Reset();

	Int_t nBins = frame2D_.GetBin(frame2D_.GetNbinsX(), frame2D_.GetNbinsY());
	for (Int_t b=1; b<nBins+1; b++) {
		eDep = frame2DCopy->GetBinContent(b);
		if (eDep == 0.0) continue;

		repeatFactor = eDep * EnergyFactor;
		newSigma = pow(repeatFactor, 0.25) / 6;

		frame2DCopy->GetBinXYZ(b, x, y, z);

		for (Int_t j = 0; j < repeatFactor; j++) {
			randX = gRandom->Gaus(x, newSigma);
			randY = gRandom->Gaus(y, newSigma);
			frame2D_.Fill(randX, randY, EnergyFactor);
		}
	}
	delete frame2DCopy;
}

void Layer::findHits(Hits* hits) {
	Int_t x; Int_t y; Int_t z;

	Int_t nBins = frame2D_.GetBin(frame2D_.GetNbinsX(), frame2D_.GetNbinsY());
	for (int i=1; i<nBins+1; i++) {
		if (frame2D_.GetBinContent(i)) {
			frame2D_.GetBinXYZ(i,x,y,z);
			
			// TODO: Check if this point is already in hits... This is where the crux might be
			
			hits->appendPoint(x,y,layerNo_);
		}
	}
}
