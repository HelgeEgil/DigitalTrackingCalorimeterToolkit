#include <TROOT.h>
#include <vector>

using namespace std;

void findManyRanges(Int_t energyFrom, Int_t energyIncrement, Int_t energyTo, Int_t mmFrom, Int_t mmIncrement, Int_t mmTo) {
   vector<Float_t> resultVector;
   ofstream file("../OutputFiles/findManyRangesHelium.csv", ofstream::out || ofstream:app);

   gROOT->ProcessLine(".L findRange.C");
   for (Int_t energy=energyFrom; energy<=energyTo; energy += energyIncrement) {
      for (Int_t mm=mmFrom; mm<=mmTo; mm += mmIncrement) {
      	findRange f(energy, mm);
         resultVector = f.Run();
         file << energy << " " << mm << " " << resultVector.at(0) << " " << resultVector.at(1) << " " << resultVector.at(2) << endl;
      }
   }

   gROOT->ProcessLine(".q");
}
