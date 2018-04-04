#include <TTree.h>
#include <TFile.h>
#include <TH2F.h>

using namespace std;

void makeDoseProfile2D() {
   TFile *f = new TFile("Data/DTC_Full_Aluminium_Absorber3mm_Degrader250mm_250MeV.root");
   TTree *tree = (TTree*) f->Get("Hits");
   Float_t x,y,z,edep;
   
   TCanvas *canvas = new TCanvas("canvas", "Dose profile", 800, 800);

   // TH2F(NAME, TITLE (MAIN;X;Y), NBINSX, XFROM, XTO, NBINSY, YFROM, YTO)
   TH2F *hDoseProfile = new TH2F("hDoseProfile", "Lateral dose profile in calorimeter -- All interactions;Depth [mm];Deposited energy [MeV]", 400, -10, 10, 400, -10, 10);

   hDoseProfile->SetFillColor(kGreen-2);
   hDoseProfile->SetLineColor(kBlack);

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("edep", &edep);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      hDoseProfile->Fill(x,y, edep);
   }

   hDoseProfile->Draw();

//   delete f;
}
