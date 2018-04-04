#include <TTree.h>
#include <TFile.h>

using namespace std;

void openROOTFile() {

   // Open file
   TFile *f = new TFile("Data/DTC_Aluminium_Absorber3mm_Degrader10mm_250MeV.root");

   // Open Tree
   TTree *tree = (TTree*) f->Get("Hits");

   // Define variables to read out
   Float_t  x,y,z,edep;
   Int_t    eventID, parentID;

   // Connect variables to Branches of ROOT tree
   tree->SetBranchAddress("posX", &x); // X position in mm
   tree->SetBranchAddress("posY", &y); // Y position in mm
   tree->SetBranchAddress("posZ", &z); // Z position in mm
   tree->SetBranchAddress("edep", &edep); // Energy deposited in interaction, in MeV
   tree->SetBranchAddress("eventID", &eventID); // Primary particle number
   tree->SetBranchAddress("parentID", &parentID); // Primary (0) or secondary (>0) particle?

   // Loop through Tree
   for (Int_t i=0; i<tree->GetEntries(); ++i) {

      // Retrieve one entry (one interaction)
      tree->GetEntry(i);

      if (eventID > 2) break; // Limit the amount of data
      if (parentID == 0) {
         printf("Primary particle with eventID %d has an interaction with %.2f keV energy loss at (x,y,z) = (%.2f, %.2f, %.2f).\n", eventID, 1e3*edep, x, y, z);
      }

      else {
         printf("Secondary particle with eventID %d has an interaction with %.2f keV energy loss at (x,y,z) = (%.2f, %.2f, %.2f).\n", eventID, 1e3*edep, x, y, z);
      }
   }

   delete f;
}
