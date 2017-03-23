#include <TROOT.h>
#include <vector>

using namespace std;

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement = 1, Int_t mmTo = -1);

void findManyRanges(Int_t degraderFrom, Int_t degraderIncrement, Int_t degraderTo, Int_t mmFrom, Int_t mmIncrement, Int_t mmTo) {
   if (mmTo < 0) mmTo = mmFrom;

   vector<Float_t> resultVector;
   ofstream file("../OutputFiles/findManyRangesDegrader.csv", ofstream::out || ofstream:app);

   gROOT->ProcessLine(".L findRange.C");
   for (Int_t degrader=degraderFrom; degrader<=degraderTo; degrader += degraderIncrement) {
      for (Int_t mm=mmFrom; mm<=mmTo; mm += mmIncrement) {
      	findRange f(250, mm, degrader);
         resultVector = f.Run();
         if (resultVector.size() > 1) {
            file << degrader << " " << mm << " " << resultVector.at(0) << " " << resultVector.at(1) << " " << resultVector.at(2) << " " <<  resultVector.at(3) << " " << resultVector.at(4) << endl;
         }
      }
   }

   gROOT->ProcessLine(".q");
}
