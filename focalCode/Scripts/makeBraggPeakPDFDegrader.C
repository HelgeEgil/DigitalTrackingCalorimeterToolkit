#include <TROOT.h>

using namespace std;

void makeBraggPeakPDFDegrader(Int_t thickness) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("drawBraggPeakGraphFit(1, 0, 1, 250, %d)", thickness)); 
	gROOT->ProcessLine(".q");
}
