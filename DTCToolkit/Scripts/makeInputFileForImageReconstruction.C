#include <TROOT.h>

using namespace std;

void makeInputFileForImageReconstruction(Int_t Runs, Int_t eventsPerRun, Int_t spotPosY) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeOutputFileForImageReconstruction(%d, %d, %d)", Runs, eventsPerRun, spotPosY));
//	gROOT->ProcessLine(".q");
}
