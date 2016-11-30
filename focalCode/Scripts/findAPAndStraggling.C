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
   TCanvas *c1 = new TCanvas("c1", "Ranges", 1200, 800);
   TCanvas *c2 = new TCanvas("c2", "Straggling", 1200, 800);
   TCanvas *c3 = new TCanvas("c3", "WE straggling", 1200, 800);

   gStyle->SetOptStat(0);

   Float_t  arrayE[200] = {0}; 
   Float_t  arrayRange[200] = {0}; 
   Float_t  arrayStraggling[200] = {0};
   Float_t  arrayWEPLStraggling[200] = {0};
   Float_t  arrayRangeWater[200] = {0};
   Float_t  range, straggling, inelasticfraction;
   Float_t  weplfactor, energyWater;
   Int_t    energy, thickness;
   Int_t    nlinesWater = 0, nlinesMaterial = 0;
   ifstream inWater, inMaterial;


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

   inMaterial.open("../OutputFiles/findManyRanges.csv");
   while (1) {
      inMaterial >> energy >> thickness >> range >> straggling >> inelasticfraction;

      if (!inMaterial.good()) {
         break;
      }

      if (thickness != absorberthickness) {
         continue;
      }

      arrayRange[nlinesMaterial] = range;
      arrayStraggling[nlinesMaterial] = straggling;

      weplfactor = arrayRangeWater[nlinesMaterial] / range;
      arrayWEPLStraggling[nlinesMaterial++] = straggling * weplfactor;
   }
   inMaterial.close();

   TGraph *gRange = new TGraph(nlinesMaterial, arrayE, arrayRange);
   TGraph *gStraggling = new TGraph(nlinesMaterial, arrayRange, arrayStraggling);
   TGraph *gWEStraggling = new TGraph(nlinesMaterial, arrayRangeWater, arrayWEPLStraggling);

   gRange->SetTitle(Form("Ranges for material in %d mm Al;Energy [MeV];Range [mm]", absorberthickness));
   gStraggling->SetTitle(Form("Straggling in Al for %d mm absorber;Range [mm];Straggling [mm]", absorberthickness));
   gWEStraggling->SetTitle(Form("WE straggling in Al for %d mm absorber;Water Equivalent Range [mm];Water Equivalent Straggling [mm]", absorberthickness));

   c1->cd();
   gRange->SetMarkerStyle(7);
   gRange->SetMarkerColor(kBlue);
   gRange->Draw("AP");
   TF1 *BK = new TF1("BK", "[0] * pow(x, [1])");
   BK->SetParameters(0.01, 1.78);
   gRange->Fit("BK", "B,Q,M");

   c2->cd();
   gStraggling->SetMarkerStyle(7);
   gStraggling->SetMarkerColor(kBlue);
   gStraggling->Draw("AP");
   TF1 *Straggling = new TF1("Straggling", "[0]*x + [1]*pow(x,2)");
   gStraggling->Fit("Straggling", "Q,M");

   c3->cd();
   gWEStraggling->SetMarkerStyle(7);
   gWEStraggling->SetMarkerColor(kBlue);
   gWEStraggling->Draw("AP");
   TF1 *WEStraggling = new TF1("WEStraggling", "[0]*x + [1]*pow(x,2)");
   gWEStraggling->Fit("WEStraggling", "Q,M");

   printf("-----------------------------------\n");
   printf("Bragg-Kleeman parameters: R = %.4f E ^ %.4f\n", BK->GetParameter(0), BK->GetParameter(1));
   printf("Straggling = %.2e * R + %.2e * R^2\n", Straggling->GetParameter(0), Straggling->GetParameter(1));
   printf("WE Straggling = %.2e * R + %.2e * R^2\n", WEStraggling->GetParameter(0), WEStraggling->GetParameter(1));

}
