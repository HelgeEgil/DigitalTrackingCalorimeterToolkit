#include <TROOT.h>
#include <iostream>

#include <vector>

using namespace std;

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement = 1, Int_t mmTo = -1);

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement, Int_t mmTo) {
   if (mmTo < 0) mmTo = mmFrom;

   // Load phase space spline
   Double_t phaseSpaceDegraderthickness[300];
   Double_t phaseSpaceEnergy[300];
   Double_t phaseSpaceEnergySpread[300];
   Double_t dt, e, es;
   Int_t idx = 0;
   ifstream in;
   in.open("../Data/Ranges/EnergyAfterDegrader.csv");

   while (1) {
      in >> dt >> e >> es;
      if (!in.good()) break;
      phaseSpaceDegraderthickness[idx] = dt;
      phaseSpaceEnergy[idx] = e;
      phaseSpaceEnergySpread[idx++] = es;
   }
   in.close();

   TSpline3 *phaseSpaceSpline = new TSpline3("phaseSpaceSpline", phaseSpaceDegraderthickness, phaseSpaceEnergy, idx);
   TSpline3 *phaseSpaceSpreadSpline = new TSpline3("phaseSpaceSpreadSpline", phaseSpaceDegraderthickness, phaseSpaceEnergySpread, idx);

   vector<Float_t> resultVector;
   ofstream file("../OutputFiles/findManyRangesDegrader.csv", ofstream::out || ofstream:app);

   gROOT->ProcessLine(".L findRange.C+");
   for (Int_t degrader=degraderFrom; degrader<=degraderTo; degrader += degraderIncrement) {
      for (Int_t mm=mmFrom; mm<=mmTo; mm += mmIncrement) {
      	findRange f(250, mm, degrader);
         resultVector = f.Run();
         if (resultVector.size() > 1) {
            Float_t expectedRange = resultVector.at(0);
            Float_t expectedSigma = resultVector.at(1);
            Float_t attenuation   = resultVector.at(2);
            Float_t expectedEnergy = phaseSpaceSpline->Eval(degrader);
            Float_t expectedEnergySpread = phaseSpaceSpreadSpline->Eval(degrader);

            // Previously, the expected energy was from resultVector.at(3) and 4 (spread)
            // However, this calculation was from an edep sum, which gave some weird results
            // Now the results are precalculated using a scoring plane after a water phantom
            // of thickness 'degrader', and interpolated to this value using splines.

            file << degrader << " " << mm << " " << expectedRange << " " << expectedSigma << " " << attenuation << " " <<  expectedEnergy << " " << expectedEnergySpread  << endl;
         }
      }
   }

   gROOT->ProcessLine(".q");
}
