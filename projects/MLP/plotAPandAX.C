#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TText.h>
#include <iostream>
#include <fstream>
#include <TVirtualFFT.h>
#include <TLatex.h>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSpline.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

const Int_t arraySize = 3500;

void plotAPandAX() {
  
   Int_t    energy = 230;

   Float_t  arrayPhantomSize200[arraySize] = {0}; // energy MC
   Float_t  arrayPhantomSizeB100200[arraySize] = {0}; // energy MC
   Float_t  arrayResidualEnergy200[arraySize] = {0};
   Float_t  arrayResidualEnergyB100200[arraySize] = {0};
   Float_t  arrayWetWepl200[arraySize] = {0};
   Float_t  arrayWetWeplB100200[arraySize] = {0};
   Float_t  arrayError200[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAX200[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAP200[arraySize] = {0};
   Float_t  arrayAXB100200[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAPB100200[arraySize] = {0};
   Float_t  arrayPhantomSize230[arraySize] = {0}; // energy MC
   Float_t  arrayResidualEnergy230[arraySize] = {0};
   Float_t  arrayWetWepl230[arraySize] = {0};
   Float_t  arrayError230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAX230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAP230[arraySize] = {0};
   Float_t  arrayAXspot[arraySize] = {0};
   Float_t  arrayAPspot[arraySize] = {0};
   Float_t  arraySpotsize[arraySize] = {0};
   
   Double_t rangesWater[arraySize] = {};
   Double_t energiesWater[arraySize] = {};
   Double_t rangesB100[arraySize] = {};
   Double_t energiesB100[arraySize] = {};
   Int_t    idxB100 = 0;
   Int_t    idxWater = 0;
   Int_t    arrayIdx200 = 0;
   Int_t    arrayIdxB100200 = 0;
   Int_t    arrayIdx230 = 0;
   Int_t    arrayIdxSpot = 0;

   
   ifstream in;
   in.open("Data/WaterPSTAR.csv");
   Float_t energy_, range_;
   while (1) {
      in >> energy_ >> range_;
      if (!in.good()) break;
      rangesWater[idxWater] = range_*10; // [mm]
      energiesWater[idxWater++] = energy_;
   }
   in.close();

   in.open("Data/B100.csv");
   while (1) {
      in >> energy_ >> range_;
      if (!in.good()) break;
      rangesB100[idxB100] = range_*10;
      energiesB100[idxB100++] = energy_;
   }
   in.close();

   TSpline3 *splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);
   TSpline3 *splineWaterInv = new TSpline3("splineWaterInv", rangesWater, energiesWater, idxWater);
   
   TSpline3 *splineB100 = new TSpline3("splineB100", energiesB100, rangesB100, idxB100);
   TSpline3 *splineB100Inv = new TSpline3("splineB100Inv", rangesB100, energiesB100, idxB100);
 
   gStyle->SetOptStat(0);
   gStyle->SetTitleFont(22);
   gStyle->SetLabelFont(22);
   gStyle->SetTextFont(22);
   gStyle->SetLabelFont(22, "Y");
   gStyle->SetTitleFont(22, "Y");
   gStyle->SetTitleYOffset(1);
   gStyle->SetLabelSize(0.045);
   gStyle->SetLabelSize(0.045, "Y");
   gStyle->SetTitleSize(0.045);
   gStyle->SetTitleSize(0.045, "Y");
   gStyle->SetTextSize(0.045);

   Float_t  phantomSize_, error_, AX_, AP_; 

   Float_t wepl = splineWater->Eval(200);
   Float_t residualEnergy, wet;

   in.open("Output/accuracy_energy200MeV_Water_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_;
      if (!in.good()) break;

      residualEnergy = 200 - splineWaterInv->Eval(phantomSize_);
      wet = wepl - splineWater->Eval(residualEnergy);

      arrayPhantomSize200[arrayIdx200] = phantomSize_;
      arrayWetWepl200[arrayIdx200] = pow(wet/wepl, 2);
      arrayResidualEnergy200[arrayIdx200] = wet;
      arrayError200[arrayIdx200] = error_;
      arrayAX200[arrayIdx200] = AX_;
      arrayAP200[arrayIdx200++] = AP_;
   }
   in.close();
   
   Float_t spotsize_;
   in.open("Output/accuracy_energy200MeV_Water_phantom_spot.csv");
   while (1) {
      if (!in.good()) break;
      in >> phantomSize_ >> spotsize_ >> error_ >> AX_ >> AP_;

      arrayAXspot[arrayIdxSpot] = AX_;
      arrayAPspot[arrayIdxSpot] = AP_;
      arraySpotsize[arrayIdxSpot++] = spotsize_;
   }
   in.close();

   in.open("Output/accuracy_energy200MeV_B100_phantom.csv");
   wepl = splineB100->Eval(200);
   while (1) {
      if (!in.good()) break;
      in >> phantomSize_ >> error_ >> AX_ >> AP_ ;

      residualEnergy = 200 - splineB100Inv->Eval(phantomSize_);
      wet = wepl - splineWater->Eval(residualEnergy);

      arrayPhantomSizeB100200[arrayIdxB100200] = phantomSize_;
      arrayWetWeplB100200[arrayIdxB100200] = pow(wet/wepl, 2);
      arrayResidualEnergyB100200[arrayIdxB100200] = wet;
      arrayAXB100200[arrayIdxB100200] = AX_;
      arrayAPB100200[arrayIdxB100200++] = AP_;
   }
   in.close();

   wepl = splineWater->Eval(230);

   in.open("Output/accuracy_energy230MeV_Water_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_;
      if (!in.good()) break;
      
      residualEnergy = 230 - splineWaterInv->Eval(phantomSize_);
      wet = wepl - splineWater->Eval(residualEnergy);

      arrayPhantomSize230[arrayIdx230] = phantomSize_;
      arrayWetWepl230[arrayIdx230] = pow(wet/wepl, 2);
      arrayResidualEnergy230[arrayIdx230] = wet;
      arrayError230[arrayIdx230] = error_;
      arrayAX230[arrayIdx230] = AX_;
      arrayAP230[arrayIdx230++] = AP_;
   }
   in.close();


   TCanvas *c1 = new TCanvas("c1", "Fit results", 1200, 800);
   c1->Divide(2,2,1e-4,1e-4);

   TGraph *gAX200 = new TGraph(arrayIdx200, arrayWetWepl200, arrayAX200);
   TGraph *gAP200 = new TGraph(arrayIdx200, arrayWetWepl200, arrayAP200);
   TGraph *gAXB100200 = new TGraph(arrayIdxB100200, arrayWetWeplB100200, arrayAXB100200);
   TGraph *gAPB100200 = new TGraph(arrayIdxB100200, arrayWetWeplB100200, arrayAPB100200);

   gAXB100200->SetMarkerColor(kRed);
   gAPB100200->SetMarkerColor(kRed);
   gAX200->SetMarkerColor(kBlue);
   gAP200->SetMarkerColor(kBlue);
   gAX200->SetMarkerStyle(21);
   gAP200->SetMarkerStyle(21);
   gAPB100200->SetMarkerStyle(21);
   gAXB100200->SetMarkerStyle(21);
   gAX200->SetMarkerSize(0.8);
   gAP200->SetMarkerSize(0.8);
   gAXB100200->SetMarkerSize(0.8);
   gAPB100200->SetMarkerSize(0.8);

   gAX200->SetTitle("Optimal A_{X} parameter for a 200 MeV beam;(WET/WEPL)^2;A_{X}");
   gAP200->SetTitle("Optimal A_{P} parameter for a 200 MeV beam;(WET/WEPL)^2;A_{P}");

   c1->cd(1);
   gAX200->Draw("AP");
   gAXB100200->Draw("P");

   c1->cd(2);
   gAP200->Draw("AP");
   gAPB100200->Draw("P");

   TLegend *leg = new TLegend(.18, .7, .68, .85);
   leg->AddEntry(gAP200, "Water", "P");
   leg->AddEntry(gAPB100200, "B100 Bone", "P");
   leg->Draw();
   
   TGraph *gAX230 = new TGraph(arrayIdx230, arrayWetWepl230, arrayAX230);
   TGraph *gAP230 = new TGraph(arrayIdx230, arrayWetWepl230, arrayAP230);

   gAX230->SetMarkerColor(kBlue);
   gAP230->SetMarkerColor(kRed);
   gAX230->SetMarkerStyle(21);
   gAP230->SetMarkerStyle(21);
   gAX230->SetMarkerSize(0.8);
   gAP230->SetMarkerSize(0.8);

   gAX230->SetTitle("Optimal A_{X} parameter for a 230 MeV beam;(WET/WEPL)^2;A_{X}");
   gAP230->SetTitle("Optimal A_{P} parameter for a 230 MeV beam;(WET/WEPL)^2;A_{P}");

   c1->cd(3);
   gAX230->Draw("AP");

   c1->cd(4);
   gAP230->Draw("AP");

   TCanvas *c2 = new TCanvas("c2", "params vs spot sizes", 1200, 800);
   c2->Divide(1,2,1e-5,1e-5);

   TGraph *gAXspot = new TGraph(arrayIdxSpot, arraySpotsize, arrayAXspot);
   TGraph *gAPspot = new TGraph(arrayIdxSpot, arraySpotsize, arrayAPspot);

   gAXspot->SetMarkerColor(kBlue);
   gAPspot->SetMarkerColor(kRed);
   gAXspot->SetMarkerStyle(21);
   gAPspot->SetMarkerStyle(21);
   gAXspot->SetMarkerSize(0.8);
   gAPspot->SetMarkerSize(0.8);

   gAXspot->SetTitle("Optimal A_{X} parameter for a 200 MeV beam;Spot size sigma [mm];A_{X}");
   gAPspot->SetTitle("Optimal A_{P} parameter for a 200 MeV beam;Spot size sigma [mm];A_{P}");
   c2->cd(1);
   gAXspot->Draw("AP");
   c2->cd(2);
   gAPspot->Draw("AP");

}
