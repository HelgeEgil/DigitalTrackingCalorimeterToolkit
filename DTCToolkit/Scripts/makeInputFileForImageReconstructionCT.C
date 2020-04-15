#include <TROOT.h>

using namespace std;

void makeInputFileForImageReconstructionCT(Int_t Runs, Int_t eventsPerRun, Int_t rotation) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeOutputFileForImageReconstructionCT(%d, %d, %d)", Runs, eventsPerRun, rotation));
//	gROOT->ProcessLine(".q");
}
