#include <TROOT.h>

using namespace std;

void makeInputFileForImageReconstructionRad(Int_t Runs, Int_t eventsPerRun, Int_t rotation, Int_t spotPosX) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeOutputFileForImageReconstructionRad(%d, %d, %d, %d)", Runs, eventsPerRun, rotation, spotPosX));
//	gROOT->ProcessLine(".q");
}
