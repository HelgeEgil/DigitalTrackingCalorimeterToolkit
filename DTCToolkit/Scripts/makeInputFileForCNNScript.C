#include <TROOT.h>

using namespace std;

void makeInputFileForCNNScript(Int_t spotPosX) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("makeInputFileForCNN(%d)", spotPosX));
}
