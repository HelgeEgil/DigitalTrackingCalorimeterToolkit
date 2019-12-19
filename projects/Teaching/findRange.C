#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>

using namespace std;

void Run() {
   TFile  * f = new TFile("output/waterphantom_78eV.root");
   Int_t energy = 200;
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, lastEventID = -1;
   Float_t  expectedRange = 0.022 * pow(energy, 1.77);

   TCanvas *c = new TCanvas("c", "c", 1000, 500);
   c->Divide(2, 1,1e-5,1e-5);

   TH1F   * doseHistogram = new TH1F("doseHistogram", "Energy deposition;Range [mm];MeV/proton", 800, 0, 400);
   TH1F   * rangeHistogram = new TH1F("rangeHistogram", "Stopping position;Range [mm];Entries", 800, 0, 400);

   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);

   printf("Found %d entries.\n", tree->GetEntries());

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if (parentID != 0) continue;
      
      doseHistogram->Fill(z, edep);

      if (eventID != lastEventID && lastEventID >= 0) {
         rangeHistogram->Fill(lastZ);
         printf("range = %.2f\n", lastZ);
         dE = 0;

      }

      dE += edep;
      lastZ = z;
      lastEventID = eventID;
   }

   doseHistogram->SetFillColor(kGray);
   doseHistogram->SetLineColor(kGray);
   doseHistogram->SetLineWidth(2);
   rangeHistogram->SetFillColor(kAzure+1);
   rangeHistogram->SetLineColor(kBlack);
   rangeHistogram->SetLineWidth(2);


   //   rangeHistogram->Scale(doseHistogram->GetMaximum() / rangeHistogram->GetMaximum() * 0.7);

   c->cd(1);
   doseHistogram->Draw();
   c->cd(2);
   rangeHistogram->Draw();


//   TF1 *fit = new TF1("fit", "gaus");
//   rangeHistogram->Fit(fit, "Q", "", fmax(expectedRange - 20, 0), expectedRange + 20);

   //printf("The range of the %d MeV proton beam is %.3f mm +- %.3f mm.\n", energy, fit->GetParameter(1), fit->GetParameter(2));
  // printf("%.1f %.3f\n", float(energy), fit->GetParameter(1));

}
