#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>

using namespace std;

void Run() {
   TFile  * f = new TFile("gate_simulation.root");
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  z,lastZ = -1;
   Int_t    eventID, parentID, lastEventID = -1;
   TH1F   * rangeHistogram = new TH1F("rangeHistogram", "Stopping position for protons;Range [mm];Entries", 800, 0, 400);

   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if (parentID != 0) continue;

      if (eventID != lastEventID && lastEventID >= 0) {
         rangeHistogram->Fill(lastZ);
      }

      lastZ = z;
      lastEventID = eventID;
   }

   rangeHistogram->Draw();

   TF1 *fit = new TF1("fit", "gaus");
   rangeHistogram->Fit(fit, "", "", 350, 400);

   printf("The range of the proton beam is %.3f mm +- %.3f mm.\n", fit->GetParameter(1), fit->GetParameter(2));

}
