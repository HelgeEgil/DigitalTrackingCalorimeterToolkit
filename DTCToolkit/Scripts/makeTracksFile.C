#include <TROOT.h>

using namespace std;

void makeTracksFile(Int_t Runs, Float_t energy) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("loadOrCreateTracks(1, %d, 0, %.2f)", Runs, energy));
	gROOT->ProcessLine(".q");
}
