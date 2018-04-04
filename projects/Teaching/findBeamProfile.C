#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>

using namespace std;

void Run() {
   TFile  * f = new TFile("../../DTCToolkit/Data/MonteCarlo/DTC_Aluminium_Absorber35mm_Degrader160mm_250MeV_4x2.root");
   Int_t energy = 200;
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  x,y,z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, baseID, level1ID, lastEventID = -1;

   TCanvas *c = new TCanvas("c", "c", 1000, 500);
   c->Divide(2, 1, 1e-6, 1e-6);

   gStyle->SetOptStat(0);

   gStyle->SetTitleFont(22);
   gStyle->SetTitleSize(0.06);
   gStyle->SetLabelFont(22, "XYZ");
   gStyle->SetTitleFont(22, "XYZ");
   gStyle->SetLabelSize(0.06, "XYZ");
   gStyle->SetTitleSize(0.06, "XYZ");

   Int_t nbins = 175;
   Int_t xrange = 30;
   Int_t yrange = 30;

   TH2F   * beamProfile1 = new TH2F("beamProfile1", "1st layer;X position [mm];Y position [mm]", nbins, -xrange, xrange, nbins, -yrange, yrange);
   TH2F   * beamProfile2 = new TH2F("beamProfile2", "31st layer;X position [mm];Y position [mm]", nbins, -xrange, xrange, nbins, -yrange, yrange);

   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("level1ID", &level1ID);
   tree->SetBranchAddress("baseID", &baseID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
   
      if (baseID == 1)     beamProfile1->Fill(x,y);
      if (level1ID == 29)  beamProfile2->Fill(x,y);

   }
   c->cd(1);
   beamProfile1->Draw("colz");
   TEllipse *e1 = new TEllipse(0, 0, 8, 4);
   e1->SetLineWidth(3);
   e1->SetLineColor(kRed-4);
   e1->SetFillStyle(0);
   e1->Draw();
   c->cd(2);
   beamProfile2->Draw("colz");
   TEllipse *e2 = new TEllipse(0, 0, 8, 4);
   e2->SetLineWidth(3);
   e2->SetLineColor(kRed-4);
   e2->SetFillStyle(0);
   e2->Draw();
}
