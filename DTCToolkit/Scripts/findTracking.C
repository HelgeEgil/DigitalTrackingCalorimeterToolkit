#include <TROOT.h>

using namespace std;

void findTracking() {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine("drawTrackAngleAtVaryingRunNumbers(0,250,50)"); 
	gROOT->ProcessLine(".q");
}
