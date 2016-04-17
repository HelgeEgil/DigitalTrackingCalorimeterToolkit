#include <TROOT.h>

using namespace std;

void makeTracksFile(Float_t energy) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("drawBraggPeakGraphFit(1, 0, 0, %.2f)", energy));
	gROOT->ProcessLine(".q");
}
