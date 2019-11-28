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

/*
   // Load phase space spline
   Double_t phaseSpaceDegraderthickness[500];
   Double_t phaseSpaceEnergy[500];
   Double_t dt, e, es;
   Int_t idx = 0;
   ifstream in;

   in.open("../Data/Ranges/EnergyAfterDegraderPSTAR.csv");

   while (1) {
      in >> dt >> e >> es;
      if (!in.good()) break;
      phaseSpaceDegraderthickness[idx] = dt;
      phaseSpaceEnergy[idx++] = e;
   }
   in.close();

   TSpline3 *phaseSpaceSpline = new TSpline3("phaseSpaceSpline", phaseSpaceDegraderthickness, phaseSpaceEnergy, idx);
   vector<Float_t> resultVector;
   
   Float_t expectedEnergy = phaseSpaceSpline->Eval(degrader);
   Float_t expectedEnergySpread = 6.11e-14*pow(mm,6) - 5.59e-11*pow(mm,5) + 1.90e-8*pow(mm,4) - 2.84e-6*pow(mm,3) + 1.57e-4*pow(mm,2) + 6.88e-3*mm + 2.23e-1;

   */

   gROOT->ProcessLine(".L Scripts/findRange.C+");
   gROOT->ProcessLine(Form("findRange f(230, %d, %d);", mm, degrader));
   gROOT->ProcessLine(Form("f.Run(100, 0, %d, %d, %d);", mm, degrader, fileIdx));//, expectedEnergy, expectedEnergySpread));
   
   /*
   Float_t expectedRange = resultVector.at(0);
   Float_t expectedSigma = resultVector.at(1);
   Float_t attenuation   = resultVector.at(2);
   
   std::ofstream filename(Form("../OutputFiles/findManyRangesDegrader_idx%d.csv", fileIdx));// , std::ofstream::out | std::ofstream::app);
   filename << degrader << " " << mm << " " << expectedRange << " " << expectedSigma << " " << attenuation << " " <<  expectedEnergy << " " << expectedEnergySpread  << endl;
   
   delete phaseSpaceSpline;
*/
   gROOT->ProcessLine(".q");
}
