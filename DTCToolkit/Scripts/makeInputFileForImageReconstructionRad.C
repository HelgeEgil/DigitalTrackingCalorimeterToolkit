#include <TROOT.h>

using namespace std;

void makeInputFileForImageReconstructionRad(Int_t Runs, Int_t eventsPerRun, Int_t rotation, Int_t spotPosX, TString phantomName) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeOutputFileForImageReconstructionRad(%d, %d, %d, %d, \"%s\")", Runs, eventsPerRun, rotation, spotPosX, phantomName.Data()));
//	gROOT->ProcessLine(".q");
}
