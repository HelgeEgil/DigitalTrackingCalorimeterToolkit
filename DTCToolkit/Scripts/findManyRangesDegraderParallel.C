#include <TROOT.h>
//#include <TSpline.h>
//#include <iostream>
//#include <fstream>

//#include <vector>
//#include "findRange.h"
//#include "findRange.C"

using namespace std;

void findManyRangesDegraderParallel(Int_t degrader, Int_t fileIdx);
void findManyRangesDegraderParallel(Int_t degrader, Int_t fileIdx) {
   gROOT->ProcessLine(".L Scripts/findRange.C+");
   gROOT->ProcessLine(Form("findRange f(230, 3, %d);", degrader));
   gROOT->ProcessLine(Form("f.Run(230, 0, 3, %d, %d);", degrader, fileIdx));//, expectedEnergy, expectedEnergySpread));
   
   gROOT->ProcessLine(".q");
}
