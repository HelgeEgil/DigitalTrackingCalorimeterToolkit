#include <TROOT.h>
#include <TSpline.h>
#include <iostream>
#include <fstream>

#include <vector>
#include "findRange.h"
#include "findRange.C"

using namespace std;

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement = 1, Int_t mmTo = -1);

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement, Int_t mmTo) {
   if (mmTo < 0) mmTo = mmFrom;

   // Load phase space spline
   Double_t phaseSpaceDegraderthickness[500];
   Double_t phaseSpaceEnergy[500];
   Double_t dt, e, es;
   Int_t idx = 0;
   ifstream in;
   printf("a\n");
   in.open("../Data/Ranges/EnergyAfterDegraderPSTAR.csv");

   while (1) {
      in >> dt >> e >> es;
      if (!in.good()) break;
      phaseSpaceDegraderthickness[idx] = dt;
      phaseSpaceEnergy[idx++] = e;
   }
   in.close();

   printf("b\n");
   TSpline3 *phaseSpaceSpline = new TSpline3("phaseSpaceSpline", phaseSpaceDegraderthickness, phaseSpaceEnergy, idx);

   vector<Float_t> resultVector;
   std::ofstream filename("../OutputFiles/findManyRangesDegrader.csv");// , std::ofstream::out | std::ofstream::app);

   printf("c\n");
//   gROOT->ProcessLine(".L findRange.C+");
   for (Int_t mm=mmFrom; mm<=mmTo; mm += mmIncrement) {
      for (Int_t degrader=degraderFrom; degrader<=degraderTo; degrader += degraderIncrement) {
         findRange f(250, mm, degrader);
         printf("e\n");
         resultVector = f.Run();
         printf("f\n");
         if (resultVector.size() > 1) {
            Float_t expectedRange = resultVector.at(0);
            Float_t expectedSigma = resultVector.at(1);
            Float_t attenuation   = resultVector.at(2);
            Float_t expectedEnergy = phaseSpaceSpline->Eval(degrader);
            Float_t expectedEnergySpread = 6.11e-14*pow(mm,6) - 5.59e-11*pow(mm,5) + 1.90e-8*pow(mm,4) - 2.84e-6*pow(mm,3) + 1.57e-4*pow(mm,2) + 6.88e-3*mm + 2.23e-1; // parameterized from energy spectrum EnergyAfterDegrader.csv, might be a bit wrong but not much
            filename << degrader << " " << mm << " " << expectedRange << " " << expectedSigma << " " << attenuation << " " <<  expectedEnergy << " " << expectedEnergySpread  << endl;
         }
      }
   }

   delete phaseSpaceSpline;

   gROOT->ProcessLine(".q");
}
