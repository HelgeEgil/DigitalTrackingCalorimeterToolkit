#include <TTree.h>
#include <TFile.h>

using namespace std;

void Run() {
   TFile *f = new TFile("gate_simulation.root");
   TTree *tree = (TTree*) f->Get("Hits");

   Float_t x,y,z,edep, time;
   Int_t eventID, parentID;

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("time", &time);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if (eventID > 5) break;
      if (parentID != 0) continue;
      printf("Primary particle with eventID %d has an interaction with %.2f MeV energy loss at (x,y,z) = (%.2f, %.2f, %.2f) with time stamp %.2e s.\n", eventID, edep, x, y, z, time);
   }

   delete f;
}
