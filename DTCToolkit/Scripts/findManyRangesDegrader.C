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

   vector<Float_t> resultVector;
   std::ofstream filename("../OutputFiles/findManyRangesDegraderHelium_final.csv");// , std::ofstream::out | std::ofstream::app);

   for (Int_t mm=mmFrom; mm<=mmTo; mm += mmIncrement) {
      for (Int_t degrader=degraderFrom; degrader<=degraderTo; degrader += degraderIncrement) {
         findRange f(917, mm, degrader);
         resultVector = f.Run();
         if (resultVector.size() > 1) {
            Float_t expectedRange = resultVector.at(0);
            Float_t expectedSigma = resultVector.at(1);
            Float_t expectedEnergy = resultVector.at(2);
            Float_t expectedEnergySpread = resultVector.at(3);
            filename << degrader << " " << mm << " " << expectedRange << " " << expectedSigma << " " <<  expectedEnergy << " " << expectedEnergySpread  << endl;
         }
      }
   }

   gROOT->ProcessLine(".q");
}
