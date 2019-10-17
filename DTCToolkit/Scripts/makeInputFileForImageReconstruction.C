#include <TROOT.h>

using namespace std;

void makeInputFileForImageReconstruction(Int_t Runs, Int_t eventsPerRun, Int_t spotPosX, Int_t spotPosY) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeInputFileForImageReconstruction(%d, %d, %d, %d)", Runs, eventsPerRun, spotPosX, spotPosY));
//	gROOT->ProcessLine(".q");
}
