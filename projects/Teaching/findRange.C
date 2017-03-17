#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>

using namespace std;

void Run(int energy) {
   TFile  * f = new TFile(Form("output/gate_simulation_%dMeV.root", energy));
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  z,lastZ = -1;
   Int_t    eventID, parentID, lastEventID = -1;
   Float_t expectedRange = 0.0262 * pow(energy, 1.736);

   TH1F   * rangeHistogram = new TH1F("rangeHistogram", "Stopping position for protons;Range [mm];Entries", 800, fmax(expectedRange - 20, 0), expectedRange + 20);

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
   rangeHistogram->Fit(fit, "Q", "", fmax(expectedRange - 20, 0), expectedRange + 20);

   //printf("The range of the %d MeV proton beam is %.3f mm +- %.3f mm.\n", energy, fit->GetParameter(1), fit->GetParameter(2));
   printf("%.1f %.3f\n", float(energy), fit->GetParameter(1));

}
