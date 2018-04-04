#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

using namespace std;

void makeDoseProfileChip() {
   TFile *f = new TFile("Data/DTC_Aluminium_Absorber3mm_Degrader250mm_250MeV.root");
   TTree *tree = (TTree*) f->Get("Hits");
   Float_t x,y,z,edep;
   
   TCanvas *canvas = new TCanvas("canvas", "Dose profile", 1200, 600);

   // TH1F(NAME, TITLE (MAIN;X;Y), NBINS, FROM, TO)
   TH1F *hDoseProfile = new TH1F("hDoseProfile", "Dose profile in calorimeter -- chip only;Depth [mm];Deposited energy [MeV]", 200, 0, 80);

   hDoseProfile->SetFillColor(kGreen-2);
   hDoseProfile->SetLineColor(kBlack);

   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      hDoseProfile->Fill(z, edep);
   }

   hDoseProfile->Draw();

//   delete f;
}
