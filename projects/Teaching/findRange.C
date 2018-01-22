#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>

using namespace std;

void Run(int energy) {
   TFile  * f = new TFile(Form("output/gate_simulation_%dMeV.root", energy));
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, lastEventID = -1;
   Float_t expectedRange = 0.0262 * pow(energy, 1.736);

   TCanvas *c = new TCanvas("c", "c", 1000, 500);
   c->Divide(2, 1,1e-5,1e-5);

   gStyle->SetTitleFont(22);
   gStyle->SetLabelFont(22);
   gStyle->SetTextFont(22);
   gStyle->SetLabelFont(22, "Y");
   gStyle->SetTitleFont(22, "Y");
   gStyle->SetTitleSize(0.06);
   gStyle->SetLabelSize(0.06);
   gStyle->SetTextSize(0.06);
   gStyle->SetLabelSize(0.06, "Y");
   gStyle->SetTitleSize(0.06, "Y");

   TH1F   * doseHistogram = new TH1F("doseHistogram", "Energy deposition;Range [mm];MeV/proton", 800, fmax(expectedRange * 0.85*0, 0), expectedRange*1.07);
   TH1F   * rangeHistogram = new TH1F("rangeHistogram", "Stopping position;Range [mm];Entries", 800, fmax(expectedRange * 0.85*0, 0), expectedRange*1.07);

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
         dE = 0;

      }

      dE += edep;
      lastZ = z;
      lastEventID = eventID;
   }

   doseHistogram->SetFillColor(kAzure+1);
   doseHistogram->SetLineColor(kBlack);
   doseHistogram->SetLineWidth(2);
   rangeHistogram->SetFillColor(kAzure+1);
   rangeHistogram->SetLineColor(kBlack);
   rangeHistogram->SetLineWidth(2);


   c->cd(1);
   doseHistogram->Draw();
   c->cd(2);
   rangeHistogram->Draw();


//   TF1 *fit = new TF1("fit", "gaus");
//   rangeHistogram->Fit(fit, "Q", "", fmax(expectedRange - 20, 0), expectedRange + 20);

   //printf("The range of the %d MeV proton beam is %.3f mm +- %.3f mm.\n", energy, fit->GetParameter(1), fit->GetParameter(2));
  // printf("%.1f %.3f\n", float(energy), fit->GetParameter(1));

}
