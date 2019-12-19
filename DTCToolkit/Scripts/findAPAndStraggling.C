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
#include <TSpline.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

void findAPAndStraggling(Int_t absorberthickness) {
   TCanvas *c1 = new TCanvas("c1", "Ranges", 1200, 1000);
   c1->Divide(2,2, 0.0001, 0.0001);

   gStyle->SetOptStat(0);

   Float_t  arrayE[1500] = {0}; 
   Float_t  arrayRange[1500] = {0};
   Float_t  arrayWET[1500] = {0};
   Float_t  arrayDegrader[1500] = {0};
   Float_t  arrayEMaterial[1500] = {0};
   Float_t  arrayStraggling[1500] = {0};
   Float_t  arrayWEPLStraggling[1500] = {0};
   Float_t  arrayRangeWater[1500] = {0};
   Float_t  arrayEnergyStraggling[1500] = {0};
   Float_t  range, straggling, inelasticfraction;
   Float_t  weplfactor, energyWater;
   Float_t  a, p, aw, pw;
   Int_t    energy, thickness;
   Float_t  energyFloat;
   Int_t    nlinesWater = 0, nlinesMaterial = 0;
   Bool_t   useDegrader = true;
   Float_t  energyStraggling;
   Int_t    degraderThickness;
   ifstream inWater, inMaterial, inEnergy;


   // 2 mm: 0.0096, 1.784
   // 3 mm: 0.0097, 1.7825
   // 4 mm: 0.0098, 1.7806
   // Focal: 0.0004461, 1.6677
   // H20:  0.0239, 1.7548

   a = 0.004461, p = 1.6677;
   aw = 0.0239, pw = 1.7548;

   inWater.open("../Data/Ranges/Water.csv");
   while (1) {
      inWater >> energyWater >> range;
      
      if (!inWater.good()) {
         break;
      }

      arrayE[nlinesWater] = energyWater;
      arrayRangeWater[nlinesWater++] = range; // file in CM
   }

   inWater.close();

   inEnergy.open("../Data/Ranges/EnergyAfterDegrader230MeV.csv");

   Double_t ead_d[500] = {};
   Double_t ead_e[500] = {};
   Double_t ead_es[500] = {};
   Int_t ead_idx = 0;

   while (1) {
      inEnergy >> degraderThickness >> energyFloat; //  >> energyStraggling;
      if (!inEnergy.good()) break;

      ead_d[ead_idx] = degraderThickness;
      ead_es[ead_idx] = energyStraggling;
      ead_e[ead_idx++] = energyFloat;
   }

   TSpline3 *spline_e = new TSpline3("spline_e", ead_d, ead_e, ead_idx);
   TSpline3 *spline_es = new TSpline3("spline_es", ead_d, ead_es, ead_idx);
   
   inEnergy.close();

   cout << "Degrader 100 -> energy = " << spline_e->Eval(100) << endl;

   if (useDegrader) {
      inMaterial.open("../OutputFiles/findManyRangesDegrader_final.csv");
   }
   else {
      inMaterial.open("../OutputFiles/findManyRanges.csv");
   }

   Float_t WET;
   while (1) {
      if (!useDegrader) {
         inMaterial >> energy >> thickness >> range >> straggling >> inelasticfraction;
      }
      else {
         inMaterial >> degraderThickness >> WET >> range >> straggling >> energyFloat >> energyStraggling;
      }

      if (!inMaterial.good()) {
         break;
      }

      if (thickness != absorberthickness) {
         continue;
      }

//      arrayEMaterial[nlinesMaterial] = energyFloat;
      arrayRange[nlinesMaterial] = range;
      arrayWET[nlinesMaterial] = WET;
      arrayDegrader[nlinesMaterial] = degraderThickness;
      arrayStraggling[nlinesMaterial] = straggling;
//      arrayEnergyStraggling[nlinesMaterial] = energyStraggling;
      arrayEMaterial[nlinesMaterial] = spline_e->Eval(degraderThickness);
      arrayEnergyStraggling[nlinesMaterial] = spline_es->Eval(degraderThickness);

//      weplfactor = arrayRangeWater[nlinesMaterial] / range;
//      weplfactor = aw / a * pow(range / aw, 1-p/pw);
      weplfactor = 2.10;
      arrayRangeWater[nlinesMaterial] = range * weplfactor; 
      arrayWEPLStraggling[nlinesMaterial++] = straggling * weplfactor;
   }
   inMaterial.close();

   std::ofstream outMaterials("../OutputFiles/final_Proton.csv");
   for (Int_t ii=0; ii<nlinesMaterial; ii++) {
      outMaterials << arrayEMaterial[ii] << " " << arrayRange[ii] << endl;
   }
   outMaterials.close();


   Int_t nLinesChange = 381;
   TGraph *gRange = new TGraph(nlinesMaterial, arrayEMaterial, arrayRange);
   TGraph *gWET = new TGraph(nlinesMaterial, arrayEMaterial, arrayWET);
   TGraph *gRangeWET = new TGraph(nlinesMaterial, arrayWET, arrayRange);
   TGraph *gEnergy = new TGraph(nlinesMaterial, arrayDegrader, arrayEMaterial);
   TGraph *gStraggling = new TGraph(nLinesChange, arrayRange, arrayStraggling);
   TGraph *gStraggling2 = new TGraph(nlinesMaterial-nLinesChange, arrayRange+nLinesChange, arrayStraggling+nLinesChange);
   TGraph *gWEStraggling = new TGraph(nlinesMaterial, arrayRangeWater, arrayWEPLStraggling);
   TGraph *gEnergyStraggling = new TGraph(nlinesMaterial, arrayEMaterial, arrayEnergyStraggling);

   if (absorberthickness <= 10) {
      gRange->SetTitle(Form("Ranges for material in %d mm Al;Energy [MeV];Range [mm]", absorberthickness));
      gStraggling->SetTitle(Form("Straggling in Al for %d mm absorber;Range [mm];Straggling [mm]", absorberthickness));
      gWEStraggling->SetTitle(Form("WE straggling in Al for %d mm absorber;Water Equivalent Range [mm];Water Equivalent Straggling [mm]", absorberthickness));
      gEnergyStraggling->SetTitle(Form("Energy straggling in DTC for %d mm absorber;Energy [MeV];Energy straggling [MeV]", absorberthickness));
   }
   else {
      gRange->SetTitle(Form("Ranges for material in %.1f mm Al;Energy [MeV];Range [mm]", float(absorberthickness)/10));
      gStraggling->SetTitle(Form("Straggling in Al for %.1f mm absorber;Range [mm];Straggling [mm]", float(absorberthickness)/10));
      gWEStraggling->SetTitle(Form("WE straggling in Al for %.1f mm absorber;Water Equivalent Range [mm];Water Equivalent Straggling [mm]", float(absorberthickness)/10));
      gEnergyStraggling->SetTitle(Form("Energy straggling in DTC for %.1f mm absorber;Energy [MeV];Energy straggling [MeV]", float(absorberthickness)/10));
   }

   gWET->SetTitle("WET range in detector; Initial energy [MeV];Water Equivalent Thickness [mm]");
   gRangeWET->SetTitle("Range in detector;Water Equivalent Thickness [mm];Physical range [mm]");

   TCanvas *c2 = new TCanvas();
   c2->Divide(2,1,1e-5,1e-5);
   c2->cd(1);
   gWET->SetMarkerStyle(7);
   gWET->SetMarkerColor(kBlue);
   gWET->Draw("AP");
   c2->cd(1);
   gRangeWET->SetMarkerStyle(7);
   gRangeWET->SetMarkerColor(kBlue);
   gRangeWET->Draw("AP");

   c1->cd(1);
   gRange->SetMarkerStyle(7);
   gRange->SetMarkerColor(kBlue);
   gRange->Draw("AP");
   TF1 *BK = new TF1("BK", "[0] * pow(x, [1])");
   BK->SetParameters(0.01, 1.77);
   gRange->Fit("BK", "B,Q,M", "", 15, 920);

   c1->cd(3);
   gStraggling->SetMarkerStyle(7);
   gStraggling->SetMarkerColor(kBlue);
   gStraggling->Draw("AP");
//   TF1 *Straggling = new TF1("Straggling", "[0]*x + [1]*pow(x,2)");
   TF1 *Straggling = new TF1("Straggling", "[0] + [1]*x");
   gStraggling->GetYaxis()->SetRangeUser(0,6);
   gStraggling->Fit("Straggling", "Q,M", "", 10, 200);
   
   gStraggling2->SetMarkerStyle(7);
   gStraggling2->SetMarkerColor(kRed);
   gStraggling2->Draw("P");

   c1->cd(4);
//   gWEStraggling->SetMarkerStyle(7);
//   gWEStraggling->SetMarkerColor(kBlue);
//   gWEStraggling->Draw("AP");
////   TF1 *WEStraggling = new TF1("WEStraggling", "[0]*x + [1]*pow(x,2)");
//   TF1 *WEStraggling = new TF1("WEStraggling", "[0] + [1]*x", 0, 290);
//   gWEStraggling->GetYaxis()->SetRangeUser(0, 6);
//   gWEStraggling->Fit("WEStraggling", "Q,M,B", "", 0, 330);

   c1->cd(4);
   gEnergyStraggling->SetMarkerStyle(7);
   gEnergyStraggling->SetMarkerColor(kBlue);
   gEnergyStraggling->Draw("AP");
   TF1 *fLandau = new TF1("fLandau", "landau");
   gEnergyStraggling->Fit("fLandau", "Q,M");
   
   c1->cd(2);
   gEnergy->SetMarkerStyle(7);
   gEnergy->SetMarkerColor(kBlue);
   gEnergy->Draw("AP");

   c1->SaveAs("../OutputFiles/Proton_final.png");


   printf("-----------------------------------\n");
   printf("Bragg-Kleeman parameters: R = %.6f E ^ %.6f\n", BK->GetParameter(0), BK->GetParameter(1));
   printf("Straggling = %.5f + %.6fR \n", Straggling->GetParameter(0), Straggling->GetParameter(1));
//   printf("Energy straggling  = Landau w/ constants %.1f,  MPV %.6f and sigma %.6f.\n", fLandau->GetParameter(0), fLandau->GetParameter(1), fLandau->GetParameter(2)); 

}
