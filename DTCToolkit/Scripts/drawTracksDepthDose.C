#include <TROOT.h>

using namespace std;

void drawTracksDepthDose(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness, Int_t eventsPerRun) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("drawTracksDepthDose(%d, %d, %d, %.0f, %.0f, %d)", Runs, dataType, recreate, energy, degraderThickness, eventsPerRun));
//	gROOT->ProcessLine(".q");
}
