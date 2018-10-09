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
   Float_t  arrayResidualEnergy200[arraySize] = {0};
   Float_t  arrayWetWepl200[arraySize] = {0};
   Float_t  arrayAX200[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAP200[arraySize] = {0};
   Float_t  arrayError200[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayPhantomSizeB100230[arraySize] = {0}; // energy MC
   Float_t  arrayResidualEnergyB100230[arraySize] = {0};
   Float_t  arrayWetWeplB100230[arraySize] = {0};
   Float_t  arrayAXB100230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAPB100230[arraySize] = {0};
   Float_t  arrayPhantomSizeAdipose230[arraySize] = {0}; // energy MC
   Float_t  arrayResidualEnergyAdipose230[arraySize] = {0};
   Float_t  arrayWetWeplAdipose230[arraySize] = {0};
   Float_t  arrayAXAdipose230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAPAdipose230[arraySize] = {0};
   Float_t  arrayPhantomSizeCorticalBone230[arraySize] = {0}; // energy MC
   Float_t  arrayResidualEnergyCorticalBone230[arraySize] = {0};
   Float_t  arrayWetWeplCorticalBone230[arraySize] = {0};
   Float_t  arrayAXCorticalBone230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAPCorticalBone230[arraySize] = {0};
   Float_t  arrayPhantomSizeA150230[arraySize] = {0}; // energy MC
   Float_t  arrayResidualEnergyA150230[arraySize] = {0};
   Float_t  arrayWetWeplA150230[arraySize] = {0};
   Float_t  arrayAXA150230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAPA150230[arraySize] = {0};
   Float_t  arrayPhantomSize230[arraySize] = {0}; // energy MC
   Float_t  arrayResidualEnergy230[arraySize] = {0};
   Float_t  arrayWetWepl230[arraySize] = {0};
   Float_t  arrayError230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAX230[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayAP230[arraySize] = {0};
   Float_t  arrayAXspot[arraySize] = {0};
   Float_t  arrayAPspot[arraySize] = {0};
   Float_t  arraySpotsize[arraySize] = {0};

   Float_t  arrayDivergence[arraySize] = {0};
   Float_t  arrayBeamSpotErrorNoTrk[arraySize] = {0};
   Float_t  arrayBeamSpotErrorEst[arraySize] = {0};
   Float_t  arrayMLPmidErrorNoTrk[arraySize] = {0};
   Float_t  arrayMLPmidErrorEst[arraySize] = {0};
   Float_t  arrayMLPstartErrorNoTrk[arraySize] = {0};
   Float_t  arrayMLPstartErrorEst[arraySize] = {0};
   
   Float_t  arrayA150Divergence[arraySize] = {0};
   Float_t  arrayA150BeamSpotErrorNoTrk[arraySize] = {0};
   Float_t  arrayA150BeamSpotErrorEst[arraySize] = {0};
   Float_t  arrayA150MLPmidErrorNoTrk[arraySize] = {0};
   Float_t  arrayA150MLPmidErrorEst[arraySize] = {0};
   Float_t  arrayA150MLPstartErrorNoTrk[arraySize] = {0};
   Float_t  arrayA150MLPstartErrorEst[arraySize] = {0};
   
   Float_t  arrayB100Divergence[arraySize] = {0};
   Float_t  arrayB100BeamSpotErrorNoTrk[arraySize] = {0};
   Float_t  arrayB100BeamSpotErrorEst[arraySize] = {0};
   Float_t  arrayB100MLPmidErrorNoTrk[arraySize] = {0};
   Float_t  arrayB100MLPmidErrorEst[arraySize] = {0};
   Float_t  arrayB100MLPstartErrorNoTrk[arraySize] = {0};
   Float_t  arrayB100MLPstartErrorEst[arraySize] = {0};
   
   Float_t  arrayAdiposeDivergence[arraySize] = {0};
   Float_t  arrayAdiposeBeamSpotErrorNoTrk[arraySize] = {0};
   Float_t  arrayAdiposeBeamSpotErrorEst[arraySize] = {0};
   Float_t  arrayAdiposeMLPmidErrorNoTrk[arraySize] = {0};
   Float_t  arrayAdiposeMLPmidErrorEst[arraySize] = {0};
   Float_t  arrayAdiposeMLPstartErrorNoTrk[arraySize] = {0};
   Float_t  arrayAdiposeMLPstartErrorEst[arraySize] = {0};
   
   Float_t  arrayCorticalBoneDivergence[arraySize] = {0};
   Float_t  arrayCorticalBoneBeamSpotErrorNoTrk[arraySize] = {0};
   Float_t  arrayCorticalBoneBeamSpotErrorEst[arraySize] = {0};
   Float_t  arrayCorticalBoneMLPmidErrorNoTrk[arraySize] = {0};
   Float_t  arrayCorticalBoneMLPmidErrorEst[arraySize] = {0};
   Float_t  arrayCorticalBoneMLPstartErrorNoTrk[arraySize] = {0};
   Float_t  arrayCorticalBoneMLPstartErrorEst[arraySize] = {0};

   Double_t rangesWater[arraySize] = {};
   Double_t energiesWater[arraySize] = {};
   Double_t rangesB100[arraySize] = {};
   Double_t energiesB100[arraySize] = {};
   Double_t rangesAdipose[arraySize] = {};
   Double_t energiesAdipose[arraySize] = {};
   Double_t rangesCorticalBone[arraySize] = {};
   Double_t energiesCorticalBone[arraySize] = {};
   Double_t rangesA150[arraySize] = {};
   Double_t energiesA150[arraySize] = {};
   Int_t    idxB100 = 0;
   Int_t    idxWater = 0;
   Int_t    idxAdipose = 0;
   Int_t    idxCorticalBone = 0;
   Int_t    idxA150 = 0;
   Int_t    arrayIdx200 = 0;
   Int_t    arrayIdxB100230 = 0;
   Int_t    arrayIdxAdipose230 = 0;
   Int_t    arrayIdxCorticalBone230 = 0;
   Int_t    arrayIdxA150230 = 0;
   Int_t    arrayIdx230 = 0;
   Int_t    arrayIdxSpot = 0;
   Int_t    arrayIdxDivergence = 0;
   Int_t    arrayA150IdxDivergence = 0;
   Int_t    arrayB100IdxDivergence = 0;
   Int_t    arrayAdiposeIdxDivergence = 0;
   Int_t    arrayCorticalBoneIdxDivergence = 0;

   
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

   in.open("Data/myAdipose.csv");
   while (1) {
      in >> energy_ >> range_;
      if (!in.good()) break;
      rangesAdipose[idxAdipose] = range_*10;
      energiesAdipose[idxAdipose++] = energy_;
   }
   in.close();
   
   in.open("Data/CorticalBone.csv");
   while (1) {
      in >> energy_ >> range_;
      if (!in.good()) break;
      rangesCorticalBone[idxCorticalBone] = range_*10;
      energiesCorticalBone[idxCorticalBone++] = energy_;
   }
   in.close();
   
   in.open("Data/A150.csv");
   while (1) {
      in >> energy_ >> range_;
      if (!in.good()) break;
      rangesA150[idxA150] = range_*10;
      energiesA150[idxA150++] = energy_;
   }
   in.close();

   TSpline3 *splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);
   TSpline3 *splineWaterInv = new TSpline3("splineWaterInv", rangesWater, energiesWater, idxWater);
   
   TSpline3 *splineB100 = new TSpline3("splineB100", energiesB100, rangesB100, idxB100);
   TSpline3 *splineB100Inv = new TSpline3("splineB100Inv", rangesB100, energiesB100, idxB100);
 
   TSpline3 *splineAdipose = new TSpline3("splineAdipose", energiesAdipose, rangesAdipose, idxAdipose);
   TSpline3 *splineAdiposeInv = new TSpline3("splineAdiposeInv", rangesAdipose, energiesAdipose, idxAdipose);
   
   TSpline3 *splineCorticalBone = new TSpline3("splineCorticalBone", energiesCorticalBone, rangesCorticalBone, idxCorticalBone);
   TSpline3 *splineCorticalBoneInv = new TSpline3("splineCorticalBoneInv", rangesCorticalBone, energiesCorticalBone, idxCorticalBone);
   
   TSpline3 *splineA150 = new TSpline3("splineA150", energiesA150, rangesA150, idxA150);
   TSpline3 *splineA150Inv = new TSpline3("splineA150Inv", rangesA150, energiesA150, idxA150);

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

   Float_t  phantomSize_, error_, AX_, AP_, residualEnergy_; 

   Float_t  wepl = splineWater->Eval(200);
   Float_t  residualEnergy, wet;
   Float_t  weplPower = 0.1;

   Float_t weplFactor = 1;
   in.open("Output/accuracy_energy200MeV_Water_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_;
      if (!in.good()) break;

      residualEnergy = 200 - splineWaterInv->Eval(phantomSize_);
      wet = wepl - splineWater->Eval(residualEnergy);

      arrayPhantomSize200[arrayIdx200] = phantomSize_;
      arrayWetWepl200[arrayIdx200] = pow(wet/wepl, 2);
      arrayResidualEnergy200[arrayIdx200] = pow(wet/wepl,10) * weplFactor;
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

   wepl = splineWater->Eval(230);
   in.open("Output/accuracy_energy230MeV_B100_phantom.csv");
   while (1) {
      if (!in.good()) break;
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;

      wet = wepl - splineWater->Eval(residualEnergy_);

      arrayPhantomSizeB100230[arrayIdxB100230] = phantomSize_;
      arrayWetWeplB100230[arrayIdxB100230] = pow(wet/wepl, 2); // * pow(weplFactor, weplPower);
      arrayResidualEnergyB100230[arrayIdxB100230] = pow(wet/wepl, 4); // * weplFactor;
      arrayAXB100230[arrayIdxB100230] = AX_;
      arrayAPB100230[arrayIdxB100230++] = AP_;
   }
   in.close();
   
   in.open("Output/accuracy_energy230MeV_Adipose_phantom.csv");
   while (1) {
      if (!in.good()) break;
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;

      wet = wepl - splineWater->Eval(residualEnergy_);

      arrayPhantomSizeAdipose230[arrayIdxAdipose230] = phantomSize_;
      arrayWetWeplAdipose230[arrayIdxAdipose230] = pow(wet/wepl, 2); // * pow(weplFactor, weplPower);
      arrayResidualEnergyAdipose230[arrayIdxAdipose230] = pow(wet/wepl, 4); // * weplFactor;
      arrayAXAdipose230[arrayIdxAdipose230] = AX_;
      arrayAPAdipose230[arrayIdxAdipose230++] = AP_;
   }
   in.close();
   
   in.open("Output/accuracy_energy230MeV_CorticalBone_phantom.csv");
   while (1) {
      if (!in.good()) break;
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;

      wet = wepl - splineWater->Eval(residualEnergy_);
      weplFactor = phantomSize_ / wet;

      arrayPhantomSizeCorticalBone230[arrayIdxCorticalBone230] = phantomSize_;
      arrayWetWeplCorticalBone230[arrayIdxCorticalBone230] = pow(wet/wepl, 2); // * pow(weplFactor, weplPower);
      arrayResidualEnergyCorticalBone230[arrayIdxCorticalBone230] = pow(wet/wepl, 4); // * weplFactor;
      arrayAXCorticalBone230[arrayIdxCorticalBone230] = AX_;
      arrayAPCorticalBone230[arrayIdxCorticalBone230++] = AP_;
   }
   in.close();
   
   in.open("Output/accuracy_energy230MeV_A150_phantom.csv");
   while (1) {
      if (!in.good()) break;
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;

      wet = wepl - splineWater->Eval(residualEnergy_);
      weplFactor = phantomSize_ / wet;

      arrayPhantomSizeA150230[arrayIdxA150230] = phantomSize_;
      arrayWetWeplA150230[arrayIdxA150230] = pow(wet/wepl, 2); // * pow(weplFactor, weplPower);
      arrayResidualEnergyA150230[arrayIdxA150230] = pow(wet/wepl, 4); // * weplFactor;
      arrayAXA150230[arrayIdxA150230] = AX_;
      arrayAPA150230[arrayIdxA150230++] = AP_;
   }
   in.close();

   in.open("Output/accuracy_energy230MeV_Water_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;
      if (!in.good()) break;
      
      wet = wepl - splineWater->Eval(residualEnergy_);
      printf("Residual energy for a %.2f mm phantom is %.2f MeV, wet/wepl = %.3f .\n", phantomSize_, residualEnergy_, wet/wepl);

      arrayPhantomSize230[arrayIdx230] = phantomSize_;
      arrayWetWepl230[arrayIdx230] = pow(wet/wepl, 2);
      arrayResidualEnergy230[arrayIdx230] = pow(wet/wepl,4);
      arrayError230[arrayIdx230] = error_;
      arrayAX230[arrayIdx230] = AX_;
      arrayAP230[arrayIdx230++] = AP_;
   }
   in.close();

   Float_t divergence_, mlpMidNoTrk_, mlpMidEst_, mlpStartNoTrk_, mlpStartEst_, bsNoTrk_, bsEst_;
   in.open("Output/MLPerror_energy230MeV_Water_degrader.csv");
   while (1) {
      in >> phantomSize_ >> divergence_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> bsNoTrk_ >> bsEst_;
      if (!in.good()) break;

      arrayDivergence[arrayIdxDivergence] = phantomSize_; //divergence_;
      arrayMLPmidErrorNoTrk[arrayIdxDivergence] = mlpMidNoTrk_;
      arrayMLPmidErrorEst[arrayIdxDivergence] = mlpMidEst_;
      arrayMLPstartErrorNoTrk[arrayIdxDivergence] = mlpStartNoTrk_;
      arrayMLPstartErrorEst[arrayIdxDivergence] = mlpStartEst_;
      arrayBeamSpotErrorNoTrk[arrayIdxDivergence] = bsNoTrk_;
      arrayBeamSpotErrorEst[arrayIdxDivergence++] = bsEst_;
   }
   in.close();
  
   Int_t dummy;
   in.open("Output/MLPerror_energy230MeV_A150_degrader.csv");
   while (1) {
      in >> phantomSize_ >> dummy >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> bsNoTrk_ >> bsEst_;
      if (!in.good()) break;

      arrayA150Divergence[arrayA150IdxDivergence] = phantomSize_; //divergence_;
      arrayA150MLPmidErrorNoTrk[arrayA150IdxDivergence] = mlpMidNoTrk_;
      arrayA150MLPmidErrorEst[arrayA150IdxDivergence] = mlpMidEst_;
      arrayA150MLPstartErrorNoTrk[arrayA150IdxDivergence] = mlpStartNoTrk_;
      arrayA150MLPstartErrorEst[arrayA150IdxDivergence] = mlpStartEst_;
      arrayA150BeamSpotErrorNoTrk[arrayA150IdxDivergence] = bsNoTrk_;
      arrayA150BeamSpotErrorEst[arrayA150IdxDivergence++] = bsEst_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_B100_degrader.csv");
   while (1) {
      in >> phantomSize_ >> dummy >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> bsNoTrk_ >> bsEst_;
      if (!in.good()) break;

      arrayB100Divergence[arrayB100IdxDivergence] = phantomSize_; //divergence_;
      arrayB100MLPmidErrorNoTrk[arrayB100IdxDivergence] = mlpMidNoTrk_;
      arrayB100MLPmidErrorEst[arrayB100IdxDivergence] = mlpMidEst_;
      arrayB100MLPstartErrorNoTrk[arrayB100IdxDivergence] = mlpStartNoTrk_;
      arrayB100MLPstartErrorEst[arrayB100IdxDivergence] = mlpStartEst_;
      arrayB100BeamSpotErrorNoTrk[arrayB100IdxDivergence] = bsNoTrk_;
      arrayB100BeamSpotErrorEst[arrayB100IdxDivergence++] = bsEst_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_Adipose_degrader.csv");
   while (1) {
      in >> phantomSize_ >> dummy >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> bsNoTrk_ >> bsEst_;
      if (!in.good()) break;

      arrayAdiposeDivergence[arrayAdiposeIdxDivergence] = phantomSize_; //divergence_;
      arrayAdiposeMLPmidErrorNoTrk[arrayAdiposeIdxDivergence] = mlpMidNoTrk_;
      arrayAdiposeMLPmidErrorEst[arrayAdiposeIdxDivergence] = mlpMidEst_;
      arrayAdiposeMLPstartErrorNoTrk[arrayAdiposeIdxDivergence] = mlpStartNoTrk_;
      arrayAdiposeMLPstartErrorEst[arrayAdiposeIdxDivergence] = mlpStartEst_;
      arrayAdiposeBeamSpotErrorNoTrk[arrayAdiposeIdxDivergence] = bsNoTrk_;
      arrayAdiposeBeamSpotErrorEst[arrayAdiposeIdxDivergence++] = bsEst_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_CorticalBone_degrader.csv");
   while (1) {
      in >> phantomSize_ >> dummy >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> bsNoTrk_ >> bsEst_;
      if (!in.good()) break;

      arrayCorticalBoneDivergence[arrayCorticalBoneIdxDivergence] = phantomSize_; //divergence_;
      arrayCorticalBoneMLPmidErrorNoTrk[arrayCorticalBoneIdxDivergence] = mlpMidNoTrk_;
      arrayCorticalBoneMLPmidErrorEst[arrayCorticalBoneIdxDivergence] = mlpMidEst_;
      arrayCorticalBoneMLPstartErrorNoTrk[arrayCorticalBoneIdxDivergence] = mlpStartNoTrk_;
      arrayCorticalBoneMLPstartErrorEst[arrayCorticalBoneIdxDivergence] = mlpStartEst_;
      arrayCorticalBoneBeamSpotErrorNoTrk[arrayCorticalBoneIdxDivergence] = bsNoTrk_;
      arrayCorticalBoneBeamSpotErrorEst[arrayCorticalBoneIdxDivergence++] = bsEst_;
   }
   in.close();

   TCanvas *c1 = new TCanvas("c1", "Fit results", 1200, 800);
   c1->Divide(2,1,1e-4,1e-4);

   TGraph *gAX200 = new TGraph(arrayIdx200, arrayResidualEnergy200, arrayAX200);
   TGraph *gAP200 = new TGraph(arrayIdx200, arrayWetWepl200, arrayAP200);
   TGraph *gAXB100230 = new TGraph(arrayIdxB100230, arrayResidualEnergyB100230, arrayAXB100230);
   TGraph *gAPB100230 = new TGraph(arrayIdxB100230, arrayWetWeplB100230, arrayAPB100230);
   TGraph *gAXAdipose230 = new TGraph(arrayIdxAdipose230, arrayResidualEnergyAdipose230, arrayAXAdipose230);
   TGraph *gAPAdipose230 = new TGraph(arrayIdxAdipose230, arrayWetWeplAdipose230, arrayAPAdipose230);
   TGraph *gAXCorticalBone230 = new TGraph(arrayIdxCorticalBone230, arrayResidualEnergyCorticalBone230, arrayAXCorticalBone230);
   TGraph *gAPCorticalBone230 = new TGraph(arrayIdxCorticalBone230, arrayWetWeplCorticalBone230, arrayAPCorticalBone230);
   TGraph *gAXA150230 = new TGraph(arrayIdxA150230, arrayResidualEnergyA150230, arrayAXA150230);
   TGraph *gAPA150230 = new TGraph(arrayIdxA150230, arrayWetWeplA150230, arrayAPA150230);

   gAX200->SetMarkerColor(kBlue);
   gAP200->SetMarkerColor(kBlue);
   gAX200->SetMarkerStyle(21);
   gAP200->SetMarkerStyle(21);
   gAX200->SetMarkerSize(0.8);
   gAP200->SetMarkerSize(0.8);
   gAXB100230->SetMarkerColor(kRed);
   gAPB100230->SetMarkerColor(kRed);
   gAPB100230->SetMarkerStyle(21);
   gAXB100230->SetMarkerStyle(21);
   gAXB100230->SetMarkerSize(0.8);
   gAPB100230->SetMarkerSize(0.8);
   gAXAdipose230->SetMarkerColor(kGreen);
   gAPAdipose230->SetMarkerColor(kGreen);
   gAPAdipose230->SetMarkerStyle(21);
   gAXAdipose230->SetMarkerStyle(21);
   gAXAdipose230->SetMarkerSize(0.8);
   gAPAdipose230->SetMarkerSize(0.8);
   gAXCorticalBone230->SetMarkerColor(kBlack);
   gAPCorticalBone230->SetMarkerColor(kBlack);
   gAPCorticalBone230->SetMarkerStyle(21);
   gAXCorticalBone230->SetMarkerStyle(21);
   gAXCorticalBone230->SetMarkerSize(0.8);
   gAPCorticalBone230->SetMarkerSize(0.8);
   gAXA150230->SetMarkerColor(kOrange+2);
   gAPA150230->SetMarkerColor(kOrange+2);
   gAPA150230->SetMarkerStyle(21);
   gAXA150230->SetMarkerStyle(21);
   gAXA150230->SetMarkerSize(0.8);
   gAPA150230->SetMarkerSize(0.8);

   gAX200->SetTitle("Optimal A_{X} parameter for a 200 MeV beam;(WET/WEPL)^{4};A_{X}");
   gAP200->SetTitle("Optimal A_{P} parameter for a 200 MeV beam;(WET/WEPL)^{2};A_{P}");

   /*
   c1->cd(1);
   gAX200->Draw("AP");

   c1->cd(2);
   gAP200->Draw("AP");
*/
   
   TGraph *gAX230 = new TGraph(arrayIdx230, arrayResidualEnergy230, arrayAX230);
   TGraph *gAP230 = new TGraph(arrayIdx230, arrayWetWepl230, arrayAP230);

   gAX230->SetMarkerColor(kBlue);
   gAP230->SetMarkerColor(kBlue);
   gAX230->SetMarkerStyle(21);
   gAP230->SetMarkerStyle(21);
   gAX230->SetMarkerSize(0.8);
   gAP230->SetMarkerSize(0.8);

   gAX230->SetTitle("Optimal A_{X} parameter for a 230 MeV beam;(WET/WEPL)^{4};A_{X}");
   gAP230->SetTitle("Optimal A_{P} parameter for a 230 MeV beam;(WET/WEPL)^{2};A_{P}");

   c1->cd(1);
   gAX230->Draw("AP");
   gAXA150230->Draw("P");
   gAXB100230->Draw("P");
   gAXCorticalBone230->Draw("P");
   gAXAdipose230->Draw("P");

   c1->cd(2);
   gAP230->Draw("AP");
   gAPA150230->Draw("P");
   gAPB100230->Draw("P");
   gAPCorticalBone230->Draw("P");
   gAPAdipose230->Draw("P");
   
   TLegend *leg = new TLegend(.3, .66, .64, .8655);
   leg->AddEntry(gAP200, "Water", "P");
   leg->AddEntry(gAPB100230, "ICRU B100 Bone", "P");
   leg->AddEntry(gAPAdipose230, "ICRU Adipose", "P");
   leg->AddEntry(gAPCorticalBone230, "ICRU Cortical Bone", "P");
   leg->AddEntry(gAPA150230, "ICRU A150 T.E.P.", "P");
   leg->Draw();

/*
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
*/
   
   TCanvas *c3 = new TCanvas("c3", "Error vs divergence", 1200, 400);
   c3->Divide(3,1,1e-3,1e-3);

   c3->cd(1);
   TGraph *gMLPstartErrorNoTrk = new TGraph(arrayIdxDivergence, arrayDivergence, arrayMLPstartErrorNoTrk);
   TGraph *gMLPstartErrorEst = new TGraph(arrayIdxDivergence, arrayDivergence, arrayMLPstartErrorEst);
   TGraph *gA150MLPstartErrorEst = new TGraph(arrayA150IdxDivergence, arrayA150Divergence, arrayA150MLPstartErrorEst);
   TGraph *gCorticalBoneMLPstartErrorEst = new TGraph(arrayCorticalBoneIdxDivergence, arrayCorticalBoneDivergence, arrayCorticalBoneMLPstartErrorEst);
   gMLPstartErrorNoTrk->SetMarkerColor(kRed);
   gMLPstartErrorNoTrk->SetMarkerStyle(21);
   gMLPstartErrorNoTrk->SetMarkerSize(0.8);
   gMLPstartErrorEst->SetMarkerColor(kBlue);
   gMLPstartErrorEst->SetMarkerStyle(21);
   gMLPstartErrorEst->SetMarkerSize(0.8);
   gA150MLPstartErrorEst->SetMarkerColor(kGreen);
   gA150MLPstartErrorEst->SetMarkerStyle(21);
   gA150MLPstartErrorEst->SetMarkerSize(0.8);
   gCorticalBoneMLPstartErrorEst->SetMarkerColor(kBlack);
   gCorticalBoneMLPstartErrorEst->SetMarkerStyle(21);
   gCorticalBoneMLPstartErrorEst->SetMarkerSize(0.8);
   gMLPstartErrorNoTrk->SetTitle("MLP - MC error at z=0;Phantom size [mm];MLP - MC [mm]");
   gMLPstartErrorNoTrk->Draw("AP");
   gMLPstartErrorEst->Draw("P");
//   gA150MLPstartErrorEst->Draw("P");
//   gCorticalBoneMLPstartErrorEst->Draw("P");
   gMLPstartErrorNoTrk->GetYaxis()->SetRangeUser(0,5);
   
  
   c3->cd(2);
   TGraph *gMLPmidErrorEst = new TGraph(arrayIdxDivergence, arrayDivergence, arrayMLPmidErrorEst);
   TGraph *gA150MLPmidErrorEst = new TGraph(arrayA150IdxDivergence, arrayA150Divergence, arrayA150MLPmidErrorEst);
   TGraph *gCorticalBoneMLPmidErrorEst = new TGraph(arrayCorticalBoneIdxDivergence, arrayCorticalBoneDivergence, arrayCorticalBoneMLPmidErrorEst);
   TGraph *gMLPmidErrorNoTrk = new TGraph(arrayIdxDivergence, arrayDivergence, arrayMLPmidErrorNoTrk);
   gMLPmidErrorNoTrk->SetMarkerColor(kRed);
   gMLPmidErrorNoTrk->SetMarkerStyle(21);
   gMLPmidErrorNoTrk->SetMarkerSize(0.8);
   gMLPmidErrorEst->SetMarkerColor(kBlue);
   gMLPmidErrorEst->SetMarkerStyle(21);
   gMLPmidErrorEst->SetMarkerSize(0.8);
   gA150MLPmidErrorEst->SetMarkerColor(kGreen);
   gA150MLPmidErrorEst->SetMarkerStyle(21);
   gA150MLPmidErrorEst->SetMarkerSize(0.8);
   gCorticalBoneMLPmidErrorEst->SetMarkerColor(kBlack);
   gCorticalBoneMLPmidErrorEst->SetMarkerStyle(21);
   gCorticalBoneMLPmidErrorEst->SetMarkerSize(0.8);
   gMLPmidErrorNoTrk->SetTitle("MLP - MC error at z=middle of phantom;Phantom size [mm];MLP - MC [mm]");
   gMLPmidErrorNoTrk->Draw("AP");
   gMLPmidErrorEst->Draw("P");
//   gA150MLPmidErrorEst->Draw("P");
//   gCorticalBoneMLPmidErrorEst->Draw("P");
   gMLPmidErrorNoTrk->GetYaxis()->SetRangeUser(0,5);

   c3->cd(3);
   TGraph *gBeamspotNoTrk = new TGraph(arrayIdxDivergence, arrayDivergence, arrayBeamSpotErrorNoTrk);
   TGraph *gBeamspotEst = new TGraph(arrayIdxDivergence, arrayDivergence, arrayBeamSpotErrorEst);
   TGraph *gA150BeamspotEst = new TGraph(arrayA150IdxDivergence, arrayA150Divergence, arrayA150BeamSpotErrorEst);
   TGraph *gCorticalBoneBeamspotEst = new TGraph(arrayCorticalBoneIdxDivergence, arrayCorticalBoneDivergence, arrayCorticalBoneBeamSpotErrorEst);
   gBeamspotNoTrk->SetMarkerColor(kRed);
   gBeamspotNoTrk->SetMarkerStyle(21);
   gBeamspotNoTrk->SetMarkerSize(0.8);
   gBeamspotEst->SetMarkerColor(kBlue);
   gBeamspotEst->SetMarkerStyle(21);
   gBeamspotEst->SetMarkerSize(0.8);
   gA150BeamspotEst->SetMarkerColor(kGreen);
   gA150BeamspotEst->SetMarkerStyle(21);
   gA150BeamspotEst->SetMarkerSize(0.8);
   gCorticalBoneBeamspotEst->SetMarkerColor(kBlack);
   gCorticalBoneBeamspotEst->SetMarkerStyle(21);
   gCorticalBoneBeamspotEst->SetMarkerSize(0.8);
   gBeamspotNoTrk->SetTitle("Beamspot sigma;Phantom size [mm];Beamspot sigma [mm]");
   gBeamspotNoTrk->Draw("AP");
   gBeamspotEst->Draw("P");
//   gA150BeamspotEst->Draw("P");
//   gCorticalBoneBeamspotEst->Draw("P");
   gBeamspotNoTrk->GetYaxis()->SetRangeUser(0,5);
   

   TLegend *leg2 = new TLegend(.52, .73, .87, .87);
   leg2->AddEntry(gBeamspotNoTrk, "Use X0 = (0,0)", "P");
   leg2->AddEntry(gBeamspotEst, "Estimate X0", "P");
//   leg2->AddEntry(gA150BeamspotEst, "Estimate X0 A150", "P");
//   leg2->AddEntry(gCorticalBoneBeamspotEst, "Estimate X0 CorticalBone", "P");
   leg2->Draw();


   TCanvas *c4 = new TCanvas("c4", "Error different materials vs w", 600,600);

   TGraph *gMLPErrorVsWWater = new TGraph(arrayIdxDivergence, arrayWetWepl230, arrayMLPmidErrorEst);
   TGraph *gMLPErrorVsWA150 = new TGraph(arrayIdxA150230, arrayWetWeplA150230, arrayA150MLPmidErrorEst);
   TGraph *gMLPErrorVsWB100 = new TGraph(arrayIdxB100230, arrayWetWeplB100230, arrayB100MLPmidErrorEst);
   TGraph *gMLPErrorVsWCorticalBone = new TGraph(arrayIdxCorticalBone230, arrayWetWeplCorticalBone230, arrayCorticalBoneMLPmidErrorEst);
   TGraph *gMLPErrorVsWAdipose = new TGraph(arrayIdxAdipose230, arrayWetWeplAdipose230, arrayAdiposeMLPmidErrorEst);

   gMLPErrorVsWWater->SetMarkerColor(kBlue);
   gMLPErrorVsWA150->SetMarkerColor(kOrange+2);
   gMLPErrorVsWB100->SetMarkerColor(kRed);
   gMLPErrorVsWCorticalBone->SetMarkerColor(kBlack);
   gMLPErrorVsWAdipose->SetMarkerColor(kGreen);

   gMLPErrorVsWWater->SetMarkerStyle(21);
   gMLPErrorVsWA150->SetMarkerStyle(21);
   gMLPErrorVsWB100->SetMarkerStyle(21);
   gMLPErrorVsWAdipose->SetMarkerStyle(21);
   gMLPErrorVsWCorticalBone->SetMarkerStyle(21);
   gMLPErrorVsWWater->SetMarkerSize(0.8);
   gMLPErrorVsWA150->SetMarkerSize(0.8);
   gMLPErrorVsWB100->SetMarkerSize(0.8);
   gMLPErrorVsWAdipose->SetMarkerSize(0.8);
   gMLPErrorVsWCorticalBone->SetMarkerSize(0.8);

   gMLPErrorVsWWater->SetTitle("MLP estimation error with #hat{X_{0}}, 230 MeV beam;(WET/WEPL)^{2};Error MLP-MC at middle of phantom[mm]");

   gMLPErrorVsWWater->Draw("AP");
   gMLPErrorVsWA150->Draw("P");
   gMLPErrorVsWB100->Draw("P");
   gMLPErrorVsWCorticalBone->Draw("P");
   gMLPErrorVsWAdipose->Draw("P");
   
   TLegend *leg3 = new TLegend(.3, .66, .64, .8655);
   leg3->AddEntry(gMLPErrorVsWWater, "Water", "P");
   leg3->AddEntry(gMLPErrorVsWB100, "ICRU B100 Bone", "P");
   leg3->AddEntry(gMLPErrorVsWAdipose, "ICRU Adipose", "P");
   leg3->AddEntry(gMLPErrorVsWCorticalBone, "ICRU Cortical Bone", "P");
   leg3->AddEntry(gMLPErrorVsWA150, "ICRU A150 T.E.P.", "P");
   leg3->Draw();
}
