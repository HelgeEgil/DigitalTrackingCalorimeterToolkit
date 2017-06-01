#include <TROOT.h>

using namespace std;

void makeBraggPeakPDFDegrader(Int_t thickness, Int_t idx = -1);
void makeBraggPeakPDFDegrader(Int_t thickness, Int_t idx) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("drawBraggPeakGraphFit(1, 0, 1, 250, %d, %d)", thickness, idx)); 
	gROOT->ProcessLine(".q");
}
