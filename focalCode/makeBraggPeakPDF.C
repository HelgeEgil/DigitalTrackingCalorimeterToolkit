#include <TROOT.h>

using namespace std;

void makeBraggPeakPDF(Float_t energy) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("drawBraggPeakFit(1, 0, 0, %.2f)", energy)); 
	gROOT->ProcessLine(".q");
	// 	Tracks * tracks = loadOrCreateTracks(1, Runs, dataType, energy);
}