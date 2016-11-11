#include <TROOT.h>

using namespace std;

void makeBraggPeakPDF(Int_t energy) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("drawBraggPeakGraphFit(1, 0, 0, %d)", energy)); 
	gROOT->ProcessLine(".q");
	// 	Tracks * tracks = loadOrCreateTracks(1, Runs, dataType, energy);
}
