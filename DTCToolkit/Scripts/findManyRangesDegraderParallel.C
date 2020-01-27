#include <TROOT.h>
//#include <TSpline.h>
//#include <iostream>
//#include <fstream>

//#include <vector>
//#include "findRange.h"
//#include "findRange.C"

using namespace std;

void findManyRangesDegraderParallel(Int_t mm, Int_t degrader, Int_t fileIdx);
void findManyRangesDegraderParallel(Int_t mm, Int_t degrader, Int_t fileIdx) {
   gROOT->ProcessLine(".L Scripts/findRange.C+");
   gROOT->ProcessLine(Form("findRange f(917, %d, %d);", mm, degrader));
   gROOT->ProcessLine(Form("f.Run(917, 0, %d, %d, %d);", mm, degrader, fileIdx));//, expectedEnergy, expectedEnergySpread));
   
   gROOT->ProcessLine(".q");
}
