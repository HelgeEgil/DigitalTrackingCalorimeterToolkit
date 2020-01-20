#include <TROOT.h>

using namespace std;

void makeInputFileForImageReconstruction(Int_t Runs, Int_t eventsPerRun, Int_t spotPosX) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeOutputFileForImageReconstruction(%d, %d, %d)", Runs, eventsPerRun, spotPosX));
//	gROOT->ProcessLine(".q");
}
