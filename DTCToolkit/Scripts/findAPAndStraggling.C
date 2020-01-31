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

   Double_t  arrayE[1500] = {0}; 
   Double_t  arrayRange[1500] = {0};
   Double_t  arrayWET[1500] = {0};
   Double_t  arrayDegrader[1500] = {0};
   Double_t  arrayEMaterial[1500] = {0};
   Double_t  arrayStraggling[1500] = {0};
   Double_t  arrayWEPLStraggling[1500] = {0};
   Double_t  arrayRangeWater[1500] = {0};
   Double_t  arrayEnergyStraggling[1500] = {0};
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
   
   inMaterial.open("../OutputFiles/findManyRangesDegrader_final_Helium.csv");

   Float_t WET;
   while (1) {
      inMaterial >> degraderThickness >> WET >> range >> straggling >> energyFloat >> energyStraggling;

      if (!inMaterial.good()) {
         break;
      }

      arrayEMaterial[nlinesMaterial] = energyFloat;
      arrayRange[nlinesMaterial] = range;
      arrayWET[nlinesMaterial] = WET + 0.6;
      arrayDegrader[nlinesMaterial] = degraderThickness;
      arrayStraggling[nlinesMaterial] = straggling;
      arrayEnergyStraggling[nlinesMaterial] = energyStraggling;
      weplfactor = 2.10;
      arrayRangeWater[nlinesMaterial] = range * weplfactor; 
      arrayWEPLStraggling[nlinesMaterial++] = straggling * weplfactor;
   }
   inMaterial.close();

   std::ofstream outMaterials("../OutputFiles/final_Helium.csv");
   for (Int_t ii=0; ii<nlinesMaterial; ii++) {
      outMaterials << arrayEMaterial[ii] << " " << arrayRange[ii] << endl;
   }
   outMaterials.close();

   std::ofstream outEnergyAfterDegrader("../OutputFiles/energyAfterDegrader_Helium.csv");
   for (Int_t ii=0; ii<nlinesMaterial; ii++) {
      outEnergyAfterDegrader << arrayDegrader[ii] << " " << arrayEMaterial[ii] << " " << arrayEnergyStraggling[ii] << endl;
   }
   outEnergyAfterDegrader.close();

   Int_t nLinesChange = 381;
   TGraph *gRange = new TGraph(nlinesMaterial, arrayEMaterial, arrayRange);
   TGraph *gWET = new TGraph(nlinesMaterial, arrayEMaterial, arrayWET);
   TGraph *gRangeWET = new TGraph(nlinesMaterial, arrayRange, arrayWET);
   TGraph *gEnergy = new TGraph(nlinesMaterial, arrayDegrader, arrayEMaterial);
   TGraph *gStraggling = new TGraph(nLinesChange, arrayRange, arrayStraggling);
   TGraph *gStraggling2 = new TGraph(nlinesMaterial-nLinesChange, arrayRange+nLinesChange, arrayStraggling+nLinesChange);
   TGraph *gWEStraggling = new TGraph(nlinesMaterial, arrayRangeWater, arrayWEPLStraggling);
   TGraph *gEnergyStraggling = new TGraph(nlinesMaterial, arrayEMaterial, arrayEnergyStraggling);

   Double_t arrayRangeInv[1500], arrayWETInv[1500];
   Int_t idx = 0;
   for (Int_t i=nlinesMaterial-1; i>=0; i--) {
      arrayRangeInv[i] = arrayRange[idx];
      arrayWETInv[i] = arrayWET[idx++];
   }

   // Make RangeWET spline
   TSpline3 *sRangeWET = new TSpline3("sRangeWET", arrayRangeInv, arrayWETInv, nlinesMaterial);
   Double_t arrayLayerNumber[45];
   Double_t arrayLayerWET[45];
   for (Int_t l=0; l<45; l++) {
      Double_t lDepth = 0;
      if (l==0) arrayLayerWET[l] = 6; // minus 0.5 MeV
      else if (l==1) arrayLayerWET[l] = 7.8; // minus 1 MeV
      else {
         lDepth = 52.4 * 2 + 3.749 + (l-2) * 5.5;
         arrayLayerWET[l] = sRangeWET->Eval(lDepth);
      }
      arrayLayerNumber[l] = l;
      printf("At layer %d, physical depth is %.2f mm and WET is %.2f mm\n", l, lDepth, sRangeWET->Eval(lDepth));
   }
   TSpline3 *sLayerWET = new TSpline3("sLayerWET", arrayLayerNumber, arrayLayerWET, 45);
   TGraph *gLayerWET = new TGraph(45, arrayLayerNumber, arrayLayerWET);


   std::ofstream outLW("../OutputFiles/layer_wet_wepl_Helium.csv");
   for (Int_t ii=0; ii<45; ii++) {
      outLW << arrayLayerNumber[ii] << " " << arrayLayerWET[ii] << " " << 330.9 - arrayLayerWET[ii]  << endl;
   }
   outLW.close();

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
   gRangeWET->SetTitle("Range in detector;Physical range [mm];Water Equivalent Thickness [mm]");
   gLayerWET->SetTitle("Layer vs WET;Layer number; WET [mm]");
   gEnergy->SetTitle("Residual energy;Degrader thickness;Residual energy");

   TCanvas *c2 = new TCanvas();
   c2->Divide(3,1,1e-5,1e-5);
   c2->cd(1);
   gWET->SetMarkerStyle(7);
   gWET->SetMarkerColor(kBlue);
   gWET->Draw("AP");
   c2->cd(2);
   gRangeWET->SetMarkerStyle(7);
   gRangeWET->SetMarkerColor(kBlue);
   gRangeWET->Draw("AP");
   sRangeWET->Draw("same");
   c2->cd(3);
   gLayerWET->SetMarkerStyle(7);
   gLayerWET->SetMarkerColor(kBlue);
   gLayerWET->Draw("AP");
   sLayerWET->Draw("same");

   c1->cd(1);
   gRange->SetMarkerStyle(7);
   gRange->SetMarkerColor(kBlue);
   gRange->Draw("AP");
   TF1 *BK = new TF1("BK", "[0] * pow(x, [1])");
   BK->SetParameters(0.01, 1.77);
   gRange->Fit("BK", "B,Q,M", "", 0, 230);

   c1->cd(3);
   gStraggling->SetMarkerStyle(7);
   gStraggling->SetMarkerColor(kBlue);
   gStraggling->Draw("AP");
//   TF1 *Straggling = new TF1("Straggling", "[0]*x + [1]*pow(x,2)");
   TF1 *Straggling = new TF1("Straggling", "[0] + [1]*x");
   gStraggling->GetYaxis()->SetRangeUser(0,6);
   gStraggling->Fit("Straggling", "Q,M", "", 100, 350);
   
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

   printf("-----------------------------------\n");
   printf("Bragg-Kleeman parameters: R = %.6f E ^ %.6f\n", BK->GetParameter(0), BK->GetParameter(1));
   printf("Straggling = %.5f + %.6fR \n", Straggling->GetParameter(0), Straggling->GetParameter(1));
//   printf("Energy straggling  = Landau w/ constants %.1f,  MPV %.6f and sigma %.6f.\n", fLandau->GetParameter(0), fLandau->GetParameter(1), fLandau->GetParameter(2)); 

}
