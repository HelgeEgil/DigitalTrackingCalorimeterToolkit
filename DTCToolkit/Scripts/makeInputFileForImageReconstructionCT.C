#include <TROOT.h>

using namespace std;

void makeInputFileForImageReconstructionCT(Int_t Runs, Int_t eventsPerRun, Int_t rotation, TString phantomName) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeOutputFileForImageReconstructionCT(%d, %d, %d, \"%s\")", Runs, eventsPerRun, rotation, phantomName.Data()));
//	gROOT->ProcessLine(".q");
}
