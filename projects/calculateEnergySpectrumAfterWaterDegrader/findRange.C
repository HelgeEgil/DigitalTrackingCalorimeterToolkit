#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>

using namespace std;

void Run(int energy) {
//   TFile  * f = new TFile(Form("output/gate_simulation_%dMeV.root", energy));
   TFile  * f = new TFile("rootoutput.root");
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, lastEventID = -1;
   Float_t expectedRange = 0.0262 * pow(energy, 1.736);
   printf("Expected range from 230 MeV beam in 160 cm water is %.2f\n", expectedRange);

   TH1F   * rangeHistogram = new TH1F("rangeHistogram", "Stopping position for protons;Range [mm];Entries", 800, fmax(expectedRange - 20, 0), expectedRange + 20);
   TH1F   * energyHistogram = new TH1F("energyHistogram", "Remaining energy protons;Energy [MeV];Entries", 200, 100, 150);

   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if (parentID != 0) continue;
      
      if (eventID != lastEventID && lastEventID >= 0) {
         rangeHistogram->Fill(lastZ);
         energyHistogram->Fill(250 - dE);
         dE = 0;
      }

      dE += edep;
      lastZ = z;
      lastEventID = eventID;
   }


   TCanvas *c = new TCanvas("c", "Ranges and energies", 1200, 600);
   c->Divide(2,1,0.001,0.001);

   c->cd(1);
   rangeHistogram->Draw();
   TF1 *fit = new TF1("fit", "gaus");
   rangeHistogram->Fit(fit, "Q", "", fmax(expectedRange - 20, 0), expectedRange + 20);

   c->cd(2);
   energyHistogram->Draw();

   //printf("The range of the %d MeV proton beam is %.3f mm +- %.3f mm.\n", energy, fit->GetParameter(1), fit->GetParameter(2));
   printf("%.1f %.3f\n", float(energy), fit->GetParameter(1));

}
