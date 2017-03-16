#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TText.h>
#include <iostream>
#include <fstream>
#include <TLatex.h>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

void findAPAndStraggling(Int_t absorberthickness) {
   TCanvas *c1 = new TCanvas("c1", "Ranges", 1200, 1000);
   c1->Divide(2,2, 0.0001, 0.0001);

   gStyle->SetOptStat(0);

   Float_t  arrayE[500] = {0}; 
   Float_t  arrayRange[500] = {0};
   Float_t  arrayEMaterial[500] = {0};
   Float_t  arrayStraggling[500] = {0};
   Float_t  arrayWEPLStraggling[500] = {0};
   Float_t  arrayRangeWater[500] = {0};
   Float_t  arrayEnergyStraggling[500] = {0};
   Float_t  range, straggling, inelasticfraction;
   Float_t  weplfactor, energyWater;
   Float_t  a, p, aw, pw;
   Int_t    energy, thickness;
   Float_t  energyFloat;
   Int_t    nlinesWater = 0, nlinesMaterial = 0;
   Bool_t   useDegrader = true;
   Float_t  energyStraggling;
   Int_t    degraderThickness;
   ifstream inWater, inMaterial;

   a = 0.0098; aw = 0.0239;
   p = 1.7806; pw = 1.7548;

   inWater.open("../Data/Ranges/Water.csv");
   while (1) {
      inWater >> energyWater >> range;
      
      if (!inWater.good()) {
         break;
      }

      arrayE[nlinesWater] = energyWater;
      arrayRangeWater[nlinesWater++] = range * 10; // file in CM
   }

   inWater.close();

   if (useDegrader) {
      inMaterial.open("../OutputFiles/findManyRangesDegrader.csv");
   }
   else {
      inMaterial.open("../OutputFiles/findManyRanges.csv");
   }

   while (1) {
      if (!useDegrader) {
         inMaterial >> energy >> thickness >> range >> straggling >> inelasticfraction;
      }
      else {
         inMaterial >> degraderThickness >> thickness >> range >> straggling >> inelasticfraction >> energyFloat >> energyStraggling;
      }

      if (!inMaterial.good()) {
         printf("!inMaterial.good() in line %d.\n", nlinesMaterial);
         printf("degraderThickness %d, thickness %d, range %.2f, straggling %.2f, inelasticfraction %.2f, energy %.2f, energyStraggling %.2f.\n", degraderThickness, thickness, range, straggling, inelasticfraction, energyFloat, energyStraggling);
         break;
      }

      if (thickness != absorberthickness) {
         continue;
      }

      arrayEMaterial[nlinesMaterial] = energyFloat;
      arrayRange[nlinesMaterial] = range;
      arrayStraggling[nlinesMaterial] = straggling;
      arrayEnergyStraggling[nlinesMaterial] = energyStraggling;

//      weplfactor = arrayRangeWater[nlinesMaterial] / range;
      weplfactor = aw / a * pow(range / aw, 1-p/pw);
      arrayRangeWater[nlinesMaterial] = range * weplfactor; 
      arrayWEPLStraggling[nlinesMaterial++] = straggling * weplfactor;
   }
   inMaterial.close();


   Int_t nLinesChange = 381;
   TGraph *gRange = new TGraph(nlinesMaterial, arrayEMaterial, arrayRange);
   TGraph *gStraggling = new TGraph(nLinesChange, arrayRange, arrayStraggling);
   TGraph *gStraggling2 = new TGraph(nlinesMaterial-nLinesChange, arrayRange+nLinesChange, arrayStraggling+nLinesChange);
   TGraph *gWEStraggling = new TGraph(nlinesMaterial, arrayRangeWater, arrayWEPLStraggling);
   TGraph *gEnergyStraggling = new TGraph(nlinesMaterial, arrayEMaterial, arrayEnergyStraggling);

   gRange->SetTitle(Form("Ranges for material in %d mm Al;Energy [MeV];Range [mm]", absorberthickness));
   gStraggling->SetTitle(Form("Straggling in Al for %d mm absorber;Range [mm];Straggling [mm]", absorberthickness));
   gWEStraggling->SetTitle(Form("WE straggling in Al for %d mm absorber;Water Equivalent Range [mm];Water Equivalent Straggling [mm]", absorberthickness));
   gEnergyStraggling->SetTitle(Form("Energy straggling in DTC for %d mm absorber;Energy [MeV];Energy straggling [MeV]", absorberthickness));

   c1->cd(1);
   gRange->SetMarkerStyle(7);
   gRange->SetMarkerColor(kBlue);
   gRange->Draw("AP");
   TF1 *BK = new TF1("BK", "[0] * pow(x, [1])");
   BK->SetParameters(0.01, 1.78);
   gRange->Fit("BK", "B,Q,M");


   c1->cd(3);
   gStraggling->SetMarkerStyle(7);
   gStraggling->SetMarkerColor(kBlue);
   gStraggling->Draw("AP");
//   TF1 *Straggling = new TF1("Straggling", "[0]*x + [1]*pow(x,2)");
   TF1 *Straggling = new TF1("Straggling", "[0] + [1]*x");
   gStraggling->GetYaxis()->SetRangeUser(0,6);
   gStraggling->Fit("Straggling", "Q,M");
   
   gStraggling2->SetMarkerStyle(7);
   gStraggling2->SetMarkerColor(kRed);
   gStraggling2->Draw("P");

   c1->cd(4);
   gWEStraggling->SetMarkerStyle(7);
   gWEStraggling->SetMarkerColor(kBlue);
   gWEStraggling->Draw("AP");
//   TF1 *WEStraggling = new TF1("WEStraggling", "[0]*x + [1]*pow(x,2)");
   TF1 *WEStraggling = new TF1("WEStraggling", "[0] + [1]*x", 0, 290);
   gWEStraggling->GetYaxis()->SetRangeUser(0, 6);
   gWEStraggling->Fit("WEStraggling", "Q,M,B", "", 0, 290);

   c1->cd(2);
   gEnergyStraggling->SetMarkerStyle(7);
   gEnergyStraggling->SetMarkerColor(kBlue);
   gEnergyStraggling->Draw("AP");
   TF1 *fLandau = new TF1("fLandau", "landau");
   gEnergyStraggling->Fit("fLandau", "Q,M");


   printf("-----------------------------------\n");
   printf("Bragg-Kleeman parameters: R = %.6f E ^ %.6f\n", BK->GetParameter(0), BK->GetParameter(1));
   printf("Straggling = %.5f + %.6fR \n", Straggling->GetParameter(0), Straggling->GetParameter(1));
   printf("WE Straggling = %.5f + %.6fR \n", WEStraggling->GetParameter(0));
   printf("Energy straggling  = Landau w/ constants %.1f,  MPV %.6f and sigma %.6f.\n", fLandau->GetParameter(0), fLandau->GetParameter(1), fLandau->GetParameter(2)); 

}
