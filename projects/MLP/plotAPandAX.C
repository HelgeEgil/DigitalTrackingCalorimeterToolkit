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
   
   Float_t  arrayA150krahwetwepl[arraySize];
   Float_t  arrayA150krahDivergence[arraySize];
   Float_t  arrayA150krahMLPstartErrorEst[arraySize];
   Float_t  arrayA150krahMLPstartErrorKrah[arraySize];
   Float_t  arrayCorticalBonekrahwetwepl[arraySize];
   Float_t  arrayCorticalBonekrahDivergence[arraySize];
   Float_t  arrayCorticalBonekrahMLPstartErrorEst[arraySize];
   Float_t  arrayCorticalBonekrahMLPstartErrorKrah[arraySize];
   
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
   
   float arrayP2Water200w[arraySize];
   float arrayP2Water200Emp[arraySize];
   float arrayP2Water200Fit[arraySize];
   float arrayP2Water230w[arraySize];
   float arrayP2Water230Emp[arraySize];
   float arrayP2Water230Fit[arraySize];
   float arrayP2CorticalBone230w[arraySize];
   float arrayP2CorticalBone230Emp[arraySize];
   float arrayP2CorticalBone230Fit[arraySize];
   float arrayP2A150230w[arraySize];
   float arrayP2A150230Emp[arraySize];
   float arrayP2A150230Fit[arraySize];
   float arrayP2B100230w[arraySize];
   float arrayP2B100230Emp[arraySize];
   float arrayP2B100230Fit[arraySize];
   float arrayP2Adipose230w[arraySize];
   float arrayP2Adipose230Emp[arraySize];
   float arrayP2Adipose230Fit[arraySize];

   float array230UncBergenEst[arraySize];
   float array230UncBergenKrahNoModel[arraySize];
   float array230UncBergenKrahModel[arraySize];
   float array230UncLLEst[arraySize];
   float array230UncLLKrahNoModel[arraySize];
   float array230UncLLKrahModel[arraySize];

   Float_t  arrayRotation[arraySize] = {0};
   Float_t  arrayRotationError[arraySize] = {0};
   
   Float_t  arraySpotsize[arraySize] = {0};
   Float_t  arraySpotsizeError[arraySize] = {0};
   
   Float_t  arraySpotsizeKrah[arraySize] = {0};
   Float_t  arraySpotsizeErrorKrah[arraySize] = {0};

   Float_t  array200wetwepl[arraySize];
   Float_t  array200Phantomsize[arraySize];
   Float_t  array200NoTrk[arraySize];
   Float_t  array200Est[arraySize];
   Float_t  array200Krah[arraySize];
   
   Float_t  array230wetwepl10mrad[arraySize];
   Float_t  array230Phantomsize10mrad[arraySize];
   Float_t  array230NoTrk10mrad[arraySize];
   Float_t  array230Est10mrad[arraySize];
   Float_t  array230Krah10mrad[arraySize];
   
   Float_t  array230wetwepl10mrad66um[arraySize];
   Float_t  array230Phantomsize10mrad66um[arraySize];
   Float_t  array230NoTrk10mrad66um[arraySize];
   Float_t  array230Est10mrad66um[arraySize];
   Float_t  array230Krah10mrad66um[arraySize];
   
   Float_t  array230wetwepl[arraySize];
   Float_t  array230Phantomsize[arraySize];
   Float_t  array230NoTrk[arraySize];
   Float_t  array230Est[arraySize];
   Float_t  array230Krah[arraySize];

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
   Int_t    idxKrahA150 = 0;
   Int_t    idxKrahCorticalBone = 0;
   Int_t    idxRotation = 0;
   Int_t    idxSpotsize = 0;
   Int_t    idxSpotsizeKrah = 0;
   Int_t    idx200 = 0;
   Int_t    idx23010mrad = 0;
   Int_t    idx23010mrad66um = 0;
   Int_t    idx230 = 0;
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
   Int_t    idx230UncLLNoModel = 0;
   Int_t    idx230UncLLModel = 0;
   Int_t idxP2CorticalBone230 = 0;
   Int_t idxP2Water200 = 0;
   Int_t idxP2Water230 = 0;
   Int_t idxP2A150230 = 0;
   Int_t idxP2B100230 = 0;
   Int_t idxP2Adipose230 = 0;

   Float_t  theoryW[arraySize];
   Float_t  theoryAX[arraySize];
   Float_t  theoryAP[arraySize];
   Int_t    theoryIdx = 0;

   float fitMu_, fitSigma_, empMu_, empSigma_;


   ifstream inTheory;
   inTheory.open("Data/theoryParams230MeV.txt");
   Float_t w_, a_, p_;
   while (1) {
      inTheory >> w_ >> a_ >> p_;
      if (!inTheory.good()) break;
      theoryW[theoryIdx] = w_;
      theoryAX[theoryIdx] = a_;
      theoryAP[theoryIdx++] = p_;
   }
   
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
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;
      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);

      arrayPhantomSize200[arrayIdx200] = phantomSize_;
      arrayWetWepl200[arrayIdx200] = pow(wet/wepl, 1);
      arrayResidualEnergy200[arrayIdx200] = pow(wet/wepl, 1) * weplFactor;
      arrayError200[arrayIdx200] = error_;
      arrayAX200[arrayIdx200] = AX_;
      arrayAP200[arrayIdx200++] = AP_;
   }
   in.close();
   
   wepl = splineWater->Eval(230);
   in.open("Output/accuracy_energy230MeV_B100_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;
      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);

      arrayPhantomSizeB100230[arrayIdxB100230] = phantomSize_;
      arrayWetWeplB100230[arrayIdxB100230] = pow(wet/wepl, 1); // * pow(weplFactor, weplPower);
      arrayResidualEnergyB100230[arrayIdxB100230] = pow(wet/wepl, 1); // * weplFactor;
      arrayAXB100230[arrayIdxB100230] = AX_;
      arrayAPB100230[arrayIdxB100230++] = AP_;

      cout << "Adding wet/wepl = " << wet/wepl << endl;
   }
   in.close();
   
   in.open("Output/accuracy_energy230MeV_Adipose_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;
      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);

      arrayPhantomSizeAdipose230[arrayIdxAdipose230] = phantomSize_;
      arrayWetWeplAdipose230[arrayIdxAdipose230] = pow(wet/wepl, 1); // * pow(weplFactor, weplPower);
      arrayResidualEnergyAdipose230[arrayIdxAdipose230] = pow(wet/wepl, 1); // * weplFactor;
      arrayAXAdipose230[arrayIdxAdipose230] = AX_;
      arrayAPAdipose230[arrayIdxAdipose230++] = AP_;
   }
   in.close();

   in.open("Output/accuracy_energy230MeV_CorticalBone_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;
      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      weplFactor = phantomSize_ / wet;

      arrayPhantomSizeCorticalBone230[arrayIdxCorticalBone230] = phantomSize_;
      arrayWetWeplCorticalBone230[arrayIdxCorticalBone230] = pow(wet/wepl, 1); // * pow(weplFactor, weplPower);
      arrayResidualEnergyCorticalBone230[arrayIdxCorticalBone230] = pow(wet/wepl, 1); // * weplFactor;
      arrayAXCorticalBone230[arrayIdxCorticalBone230] = AX_;
      arrayAPCorticalBone230[arrayIdxCorticalBone230++] = AP_;
   }
   in.close();
   
   in.open("Output/accuracy_energy230MeV_A150_phantom.csv");
   while (1) {
      in >> phantomSize_ >> error_ >> AX_ >> AP_ >> residualEnergy_;
      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      weplFactor = phantomSize_ / wet;

      arrayPhantomSizeA150230[arrayIdxA150230] = phantomSize_;
      arrayWetWeplA150230[arrayIdxA150230] = pow(wet/wepl, 1); // * pow(weplFactor, weplPower);
      arrayResidualEnergyA150230[arrayIdxA150230] = pow(wet/wepl, 1); // * weplFactor;
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
      arrayWetWepl230[arrayIdx230] = pow(wet/wepl, 1);
      arrayResidualEnergy230[arrayIdx230] = pow(wet/wepl, 1);
      arrayError230[arrayIdx230] = error_;
      arrayAX230[arrayIdx230] = AX_;
      arrayAP230[arrayIdx230++] = AP_;
   }
   in.close();

   Float_t divergence_, mlpMidNoTrk_, mlpMidEst_, mlpStartNoTrk_, mlpStartEst_, bsNoTrk_, bsEst_;
   Int_t dummy;
   Float_t fdummy;
   Float_t mlpStartKrah_, mlpMidKrah_, rotation_;
   in.open("Output/MLPerror_energy230MeV_Water_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;
      if (!in.good()) break;

      arrayDivergence[arrayIdxDivergence] = phantomSize_; //divergence_;
      arrayMLPmidErrorNoTrk[arrayIdxDivergence] = mlpMidNoTrk_;
      arrayMLPmidErrorEst[arrayIdxDivergence] = mlpMidEst_;
      arrayMLPstartErrorNoTrk[arrayIdxDivergence] = mlpStartKrah_;
      arrayMLPstartErrorEst[arrayIdxDivergence++] = mlpStartEst_;
   }
   in.close();
  
   in.open("Output/MLPerror_energy230MeV_A150_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;
      if (!in.good()) break;

      arrayA150Divergence[arrayA150IdxDivergence] = phantomSize_; //divergence_;
      arrayA150MLPmidErrorNoTrk[arrayA150IdxDivergence] = mlpMidNoTrk_;
      arrayA150MLPmidErrorEst[arrayA150IdxDivergence] = mlpMidEst_;
      arrayA150MLPstartErrorNoTrk[arrayA150IdxDivergence] = mlpStartKrah_;
      arrayA150MLPstartErrorEst[arrayA150IdxDivergence++] = mlpStartEst_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_B100_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;
      if (!in.good()) break;

      arrayB100Divergence[arrayB100IdxDivergence] = phantomSize_; //divergence_;
      arrayB100MLPmidErrorNoTrk[arrayB100IdxDivergence] = mlpMidNoTrk_;
      arrayB100MLPmidErrorEst[arrayB100IdxDivergence] = mlpMidEst_;
      arrayB100MLPstartErrorNoTrk[arrayB100IdxDivergence] = mlpStartKrah_;
      arrayB100MLPstartErrorEst[arrayB100IdxDivergence++] = mlpStartEst_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_Adipose_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;
      if (!in.good()) break;

      arrayAdiposeDivergence[arrayAdiposeIdxDivergence] = phantomSize_; //divergence_;
      arrayAdiposeMLPmidErrorNoTrk[arrayAdiposeIdxDivergence] = mlpMidNoTrk_;
      arrayAdiposeMLPmidErrorEst[arrayAdiposeIdxDivergence] = mlpMidEst_;
      arrayAdiposeMLPstartErrorNoTrk[arrayAdiposeIdxDivergence] = mlpStartKrah_;
      arrayAdiposeMLPstartErrorEst[arrayAdiposeIdxDivergence++] = mlpStartEst_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_CorticalBone_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;
      if (!in.good()) break;

      arrayCorticalBoneDivergence[arrayCorticalBoneIdxDivergence] = phantomSize_; //divergence_;
      arrayCorticalBoneMLPmidErrorNoTrk[arrayCorticalBoneIdxDivergence] = mlpMidNoTrk_;
      arrayCorticalBoneMLPmidErrorEst[arrayCorticalBoneIdxDivergence] = mlpMidEst_;
      arrayCorticalBoneMLPstartErrorNoTrk[arrayCorticalBoneIdxDivergence] = mlpStartKrah_;
      arrayCorticalBoneMLPstartErrorEst[arrayCorticalBoneIdxDivergence++] = mlpStartEst_;
   }
   in.close();

   in.open("Output/MLPerror_energy230MeV_Water_rotation.csv");
   while (1) {
      in >> dummy >> rotation_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> bsNoTrk_ >> bsEst_;
      if (!in.good()) break;

      arrayRotation[idxRotation] = rotation_;
      arrayRotationError[idxRotation++] = mlpStartEst_;
   }
   in.close();

   /*
   in.open("Output/MLPerror_energy230MeV_Water_spotsize.csv");
   while (1) {
      in >> dummy >> rotation_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> bsNoTrk_ >> bsEst_;
      if (!in.good()) break;

   }
   in.close();
*/

   in.open("Output/MLPerror_energy230MeV_Water_spotsize_krah.csv");
   while (1) {
      in >> dummy >> rotation_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_;
      if (!in.good()) break;

      arraySpotsizeKrah[idxSpotsizeKrah] = rotation_;
      arraySpotsizeErrorKrah[idxSpotsizeKrah++] = mlpStartKrah_;
      arraySpotsize[idxSpotsize] = rotation_;
      arraySpotsizeError[idxSpotsize++] = mlpStartEst_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_Water_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;

      if (!in.good()) break;
         
      wet = wepl - splineWater->Eval(residualEnergy_);
      array230wetwepl[idx230] = wet/wepl;
      array230Phantomsize[idx230] = phantomSize_;
      array230NoTrk[idx230] = mlpMidNoTrk_;
      array230Est[idx230] = mlpStartEst_;
      array230Krah[idx230++] = mlpStartKrah_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy230MeV_Water_krah_uib_nomodel.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;

      if (!in.good()) break;
         
      wet = wepl - splineWater->Eval(residualEnergy_);
      array230wetwepl10mrad[idx23010mrad] = wet/wepl;
      array230Phantomsize10mrad[idx23010mrad] = phantomSize_;
      array230NoTrk10mrad[idx23010mrad] = mlpStartNoTrk_;
      array230UncBergenEst[idx23010mrad] = mlpStartEst_;
      array230UncBergenKrahNoModel[idx23010mrad++] = mlpStartKrah_;
   }
   in.close();

   in.open("Output/MLPerror_energy230MeV_Water_krah_uib_model.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;

      if (!in.good()) break;
         
      wet = wepl - splineWater->Eval(residualEnergy_);
      array230wetwepl10mrad66um[idx23010mrad66um] = wet/wepl;
      array230Phantomsize10mrad66um[idx23010mrad66um] = phantomSize_;
      array230NoTrk10mrad66um[idx23010mrad66um] = mlpStartNoTrk_;
//      array230UncBergenEst[idx23010mrad66um] = mlpStartEst_;
      array230UncBergenKrahModel[idx23010mrad66um++] = mlpStartKrah_;
   }
   in.close();

   in.open("Output/MLPerror_energy230MeV_Water_krah_ll_nomodel.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;

      if (!in.good()) break;
         
      wet = wepl - splineWater->Eval(residualEnergy_);
      array230wetwepl10mrad[idx230UncLLNoModel] = wet/wepl;
      array230Phantomsize10mrad[idx230UncLLNoModel] = phantomSize_;
      array230UncLLEst[idx230UncLLNoModel] = mlpStartEst_;
      array230UncLLKrahNoModel[idx230UncLLNoModel++] = mlpStartKrah_;
   }
   in.close();

   in.open("Output/MLPerror_energy230MeV_Water_krah_ll_model.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;

      if (!in.good()) break;
         
      wet = wepl - splineWater->Eval(residualEnergy_);
      array230wetwepl10mrad66um[idx230UncLLModel] = wet/wepl;
      array230Phantomsize10mrad66um[idx230UncLLModel] = phantomSize_;
      array230UncLLKrahModel[idx230UncLLModel++] = mlpStartKrah_;
   }
   in.close();


   wepl = splineWater->Eval(200);
   in.open("Output/MLPerror_energy200MeV_Water_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_ >> fdummy;

      if (!in.good()) break;
         
      wet = wepl - splineWater->Eval(residualEnergy_);
      array200wetwepl[idx200] = wet/wepl;
      array200Phantomsize[idx200] = phantomSize_;
      array200NoTrk[idx200] = mlpMidNoTrk_;
      array200Est[idx200] = mlpStartEst_;
      array200Krah[idx200++] = mlpStartKrah_;
   }
   in.close();
   
   
   in.open("Output/MLPerror_energy200MeV_A150_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayA150krahwetwepl[idxKrahA150] = wet/wepl;
      arrayA150krahDivergence[idxKrahA150] = phantomSize_;
      arrayA150krahMLPstartErrorEst[idxKrahA150] = mlpMidEst_;
      arrayA150krahMLPstartErrorKrah[idxKrahA150++] = mlpMidKrah_;
   }
   in.close();
   
   in.open("Output/MLPerror_energy200MeV_CorticalBone_krah.csv");
   while (1) {
      in >> phantomSize_ >> mlpStartNoTrk_ >> mlpMidNoTrk_ >> mlpStartEst_ >> mlpMidEst_ >> mlpStartKrah_ >> mlpMidKrah_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayCorticalBonekrahwetwepl[idxKrahCorticalBone] = wet/wepl;
      arrayCorticalBonekrahDivergence[idxKrahCorticalBone] = phantomSize_;
      arrayCorticalBonekrahMLPstartErrorEst[idxKrahCorticalBone] = mlpMidEst_;
      arrayCorticalBonekrahMLPstartErrorKrah[idxKrahCorticalBone++] = mlpMidKrah_;
   }
   in.close();
   
   wepl = splineWater->Eval(200);
   in.open("Output/P2_energy200MeV_Water_krah.csv");
   while (1) {
      in >> phantomSize_ >> empMu_ >> empSigma_ >> fitMu_ >> fitSigma_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayP2Water200w[idxP2Water200] = pow(wet/wepl, 2);
      arrayP2Water200Emp[idxP2Water200] = empMu_ + 2.5* empSigma_;
      arrayP2Water200Fit[idxP2Water200++] = fitMu_ + 2.5*fitSigma_;
   }
   in.close();
   
   wepl = splineWater->Eval(230);
   in.open("Output/P2_energy230MeV_Water_krah.csv");
   while (1) {
      in >> phantomSize_ >> empMu_ >> empSigma_ >> fitMu_ >> fitSigma_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayP2Water230w[idxP2Water230] = pow(wet/wepl, 2);
      arrayP2Water230Emp[idxP2Water230] = empMu_ + 2.5* empSigma_;
      arrayP2Water230Fit[idxP2Water230++] = fitMu_ + 2.5*fitSigma_;
   }
   in.close();

   in.open("Output/P2_energy230MeV_Adipose_krah.csv");
   while (1) {
      in >> phantomSize_ >> empMu_ >> empSigma_ >> fitMu_ >> fitSigma_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayP2Adipose230w[idxP2Adipose230] = pow(wet/wepl, 2);
      arrayP2Adipose230Emp[idxP2Adipose230] = empMu_ + 2.5* empSigma_;
      arrayP2Adipose230Fit[idxP2Adipose230++] = fitMu_ + 2.5*fitSigma_;
   }
   in.close();

   in.open("Output/P2_energy230MeV_B100_krah.csv");
   while (1) {
      in >> phantomSize_ >> empMu_ >> empSigma_ >> fitMu_ >> fitSigma_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayP2B100230w[idxP2B100230] = pow(wet/wepl, 2);
      arrayP2B100230Emp[idxP2B100230] = empMu_ + 2.5* empSigma_;
      arrayP2B100230Fit[idxP2B100230++] = fitMu_ + 2.5*fitSigma_;
   }
   in.close();

   in.open("Output/P2_energy230MeV_A150_krah.csv");
   while (1) {
      in >> phantomSize_ >> empMu_ >> empSigma_ >> fitMu_ >> fitSigma_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayP2A150230w[idxP2A150230] = pow(wet/wepl, 2);
      arrayP2A150230Emp[idxP2A150230] = empMu_ + 2.5* empSigma_;
      arrayP2A150230Fit[idxP2A150230++] = fitMu_ + 2.5*fitSigma_;
   }
  in.close();
 
   in.open("Output/P2_energy230MeV_CorticalBone_krah.csv");
   while (1) {
      in >> phantomSize_ >> empMu_ >> empSigma_ >> fitMu_ >> fitSigma_ >> residualEnergy_;

      if (!in.good()) break;

      wet = wepl - splineWater->Eval(residualEnergy_);
      arrayP2CorticalBone230w[idxP2CorticalBone230] = pow(wet/wepl, 2);
      arrayP2CorticalBone230Emp[idxP2CorticalBone230] = empMu_ + 2.5* empSigma_;
      arrayP2CorticalBone230Fit[idxP2CorticalBone230++] = fitMu_ + 2.5*fitSigma_;
   }
   in.close();

   
   TCanvas *c1 = new TCanvas("c1", "Fit results", 1500, 1200);
//   c1->Divide(2,1,1e-4,1e-4);

   Float_t  arrayAPall[arraySize];
   Float_t  arrayAXall[arraySize];
   Float_t  arrayWetWeplAll[arraySize];
   Float_t  arrayResidualEnergyAll[arraySize];
   Int_t fullIdx = 0;
   for (Int_t i=0; i<arrayIdxAdipose230; i++) {
      if (i<arrayIdx230) { // water
         arrayAPall[fullIdx] = arrayAP230[i];
         arrayWetWeplAll[fullIdx] = arrayWetWepl230[i];
         arrayAXall[fullIdx] = arrayAX230[i];
         arrayResidualEnergyAll[fullIdx++] = arrayResidualEnergy230[i];
      }

      if (i<arrayIdxB100230) { // B100
         arrayAPall[fullIdx] = arrayAPB100230[i];
         arrayWetWeplAll[fullIdx] = arrayWetWeplB100230[i];
         arrayAXall[fullIdx] = arrayAXB100230[i];
         arrayResidualEnergyAll[fullIdx++] = arrayResidualEnergyB100230[i];
      }
      
      if (i<arrayIdxA150230) { // A150
         arrayAPall[fullIdx] = arrayAPA150230[i];
         arrayWetWeplAll[fullIdx] = arrayWetWeplA150230[i];
         arrayAXall[fullIdx] = arrayAXA150230[i];
         arrayResidualEnergyAll[fullIdx++] = arrayResidualEnergyA150230[i];
      }
      
      if (i<arrayIdxAdipose230) { // Adipose
         arrayAPall[fullIdx] = arrayAPAdipose230[i];
         arrayWetWeplAll[fullIdx] = arrayWetWeplAdipose230[i];
         arrayAXall[fullIdx] = arrayAXAdipose230[i];
         arrayResidualEnergyAll[fullIdx++] = arrayResidualEnergyAdipose230[i];
      }
      
      if (i<arrayIdxCorticalBone230) { // A150
         arrayAPall[fullIdx] = arrayAPCorticalBone230[i];
         arrayWetWeplAll[fullIdx] = arrayWetWeplCorticalBone230[i];
         arrayAXall[fullIdx] = arrayAXCorticalBone230[i];
         arrayResidualEnergyAll[fullIdx++] = arrayResidualEnergyCorticalBone230[i];
      }
   }

   printf("In total %d indices in arrayA[X,P]all.\n", fullIdx);

   gStyle->SetMarkerSize(1.2);

   TGraph *gAX230 = new TGraph(arrayIdx230, arrayResidualEnergy230, arrayAX230);
   TGraph *gAP230 = new TGraph(arrayIdx230, arrayWetWepl230, arrayAP230);
   TGraph *gAXB100230 = new TGraph(arrayIdxB100230, arrayResidualEnergyB100230, arrayAXB100230);
   TGraph *gAPB100230 = new TGraph(arrayIdxB100230, arrayWetWeplB100230, arrayAPB100230);
   TGraph *gAXAdipose230 = new TGraph(arrayIdxAdipose230, arrayResidualEnergyAdipose230, arrayAXAdipose230);
   TGraph *gAPAdipose230 = new TGraph(arrayIdxAdipose230, arrayWetWeplAdipose230, arrayAPAdipose230);
   TGraph *gAXCorticalBone230 = new TGraph(arrayIdxCorticalBone230, arrayResidualEnergyCorticalBone230, arrayAXCorticalBone230);
   TGraph *gAPCorticalBone230 = new TGraph(arrayIdxCorticalBone230, arrayWetWeplCorticalBone230, arrayAPCorticalBone230);
   TGraph *gAXA150230 = new TGraph(arrayIdxA150230, arrayResidualEnergyA150230, arrayAXA150230);
   TGraph *gAPA150230 = new TGraph(arrayIdxA150230, arrayWetWeplA150230, arrayAPA150230);

   TGraph *gAX200 = new TGraph(arrayIdx200, arrayResidualEnergy200, arrayAX200);
   TGraph *gAP200 = new TGraph(arrayIdx200, arrayResidualEnergy200, arrayAP200);

   TGraph *gAXall = new TGraph(fullIdx, arrayResidualEnergyAll, arrayAXall);
   TGraph *gAPall = new TGraph(fullIdx, arrayWetWeplAll, arrayAPall);

   TGraph *gAXtheory = new TGraph(theoryIdx, theoryW, theoryAX);
   TGraph *gAPtheory = new TGraph(theoryIdx, theoryW, theoryAP);
      
   gAXall->SetTitle("All materials;wet/wepl;A_{X}");
   gAPall->SetTitle("All materials;wet/wepl;A_{P}");

   gAXall->SetMarkerStyle(21); gAPall->SetMarkerStyle(21);
   gAX200->SetMarkerStyle(21); gAP200->SetMarkerStyle(21);

   gAX200->SetMarkerColor(kRed);
   gAP200->SetMarkerColor(kRed);
/*   
   c1->cd(3);
//   gAXall->Draw("AP");
   gAX200->Draw("AP");
   c1->cd(4);  
//   gAPall->Draw("AP");
   gAP200->Draw("AP");
*/
   gAXB100230->SetMarkerColor(kRed);
   gAPB100230->SetMarkerColor(kRed);
   gAPB100230->SetMarkerStyle(21);
   gAXB100230->SetMarkerStyle(21);
   gAXAdipose230->SetMarkerColor(kGreen);
   gAPAdipose230->SetMarkerColor(kGreen);
   gAPAdipose230->SetMarkerStyle(21);
   gAXAdipose230->SetMarkerStyle(21);
   gAXCorticalBone230->SetMarkerColor(kBlack);
   gAPCorticalBone230->SetMarkerColor(kBlack);
   gAPCorticalBone230->SetMarkerStyle(21);
   gAXCorticalBone230->SetMarkerStyle(21);
   gAXA150230->SetMarkerColor(kOrange+2);
   gAPA150230->SetMarkerColor(kOrange+2);
   gAPA150230->SetMarkerStyle(21);
   gAXA150230->SetMarkerStyle(21);


   gAX230->SetMarkerColor(kBlue);
   gAP230->SetMarkerColor(kBlue);
   gAX230->SetMarkerStyle(21);
   gAP230->SetMarkerStyle(21);

   gAX230->SetTitle(";WET/WEPL;A_{X} parameter");
   gAP230->SetTitle(";WET/WEPL;A_{P} parameter");
//   c1->cd(1);
   gAX230->Draw("AP");
   gAXA150230->Draw("P");
   gAXB100230->Draw("P");
   gAXCorticalBone230->Draw("P");
   gAXAdipose230->Draw("P");
  
   TF1 *fitX = new TF1("fitX", "pol4");
   fitX->SetParameters(1, -0.30, 1.53, -4.25, 2.30);
   fitX->SetLineColor(kBlack);
   fitX->SetLineStyle(9);
   fitX->SetLineWidth(4);
   fitX->Draw("same");

   gAXtheory->SetLineColor(kRed);
   gAXtheory->SetLineStyle(9);
   gAXtheory->SetLineWidth(4);
   gAXtheory->Draw("same");

   TLegend *legX = new TLegend(.3, .66, .64, .8655);
   legX->SetTextFont(22);
   legX->AddEntry(gAX230, "Water", "P");
   legX->AddEntry(gAXB100230, "ICRU B100 Bone", "P");
   legX->AddEntry(gAXAdipose230, "ICRU Adipose", "P");
   legX->AddEntry(gAXCorticalBone230, "ICRU Cortical Bone", "P");
   legX->AddEntry(gAXA150230, "ICRU A150 T.E.P.", "P");
   legX->AddEntry(fitX, "Polynomial fit", "L");
   legX->AddEntry(gAXtheory, "Theoretical prediction", "L");
   legX->Draw();
   
   /*
   //c1->cd(2);
   gAP230->Draw("AP");
   gAPA150230->Draw("P");
   gAPB100230->Draw("P");
   gAPCorticalBone230->Draw("P");
   gAPAdipose230->Draw("P");

   TF1 *fitP = new TF1("fitP", "pol5");
   fitP->SetParameters(-0.659, 2.072, -8.66, 18.14, -16.27, 5.34);
   fitP->SetLineColor(kBlack);
   fitP->SetLineWidth(4);
   fitP->SetLineStyle(9);
   fitP->Draw("same");

   gAPtheory->SetLineColor(kRed);
   gAPtheory->SetLineStyle(9);
   gAPtheory->SetLineWidth(4);
   gAPtheory->Draw("same");

   TLegend *leg = new TLegend(.3, .66, .64, .8655);
   leg->SetTextFont(22);
   leg->AddEntry(gAP230, "Water", "P");
   leg->AddEntry(gAPB100230, "ICRU B100 Bone", "P");
   leg->AddEntry(gAPAdipose230, "ICRU Adipose", "P");
   leg->AddEntry(gAPCorticalBone230, "ICRU Cortical Bone", "P");
   leg->AddEntry(gAPA150230, "ICRU A150 T.E.P.", "P");
   leg->AddEntry(fitP, "Polynomial fit", "L");
   leg->AddEntry(gAPtheory, "Theoretical prediction", "L");
   leg->SetTextFont(22);
   leg->Draw(); 
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

   gStyle->SetMarkerSize(1.2);

   TCanvas *c4 = new TCanvas("c4", "Error different materials vs w", 1000,800);

   TGraph *gMLPErrorVsWWater = new TGraph(arrayIdxDivergence, arrayWetWepl230, arrayMLPstartErrorEst);
   TGraph *gMLPErrorVsWA150 = new TGraph(arrayIdxA150230, arrayWetWeplA150230, arrayA150MLPstartErrorEst);
   TGraph *gMLPErrorVsWB100 = new TGraph(arrayIdxB100230, arrayWetWeplB100230, arrayB100MLPstartErrorEst);
   TGraph *gMLPErrorVsWCorticalBone = new TGraph(arrayIdxCorticalBone230, arrayWetWeplCorticalBone230, arrayCorticalBoneMLPstartErrorEst);
   TGraph *gMLPErrorVsWAdipose = new TGraph(arrayIdxAdipose230, arrayWetWeplAdipose230, arrayAdiposeMLPstartErrorEst);
   
   TGraph *gMLPErrorVsWWaterK = new TGraph(arrayIdxDivergence, arrayWetWepl230, arrayMLPstartErrorNoTrk);
   TGraph *gMLPErrorVsWA150K = new TGraph(arrayIdxA150230, arrayWetWeplA150230, arrayA150MLPstartErrorNoTrk);
   TGraph *gMLPErrorVsWB100K = new TGraph(arrayIdxB100230, arrayWetWeplB100230, arrayB100MLPstartErrorNoTrk);
   TGraph *gMLPErrorVsWCorticalBoneK = new TGraph(arrayIdxCorticalBone230, arrayWetWeplCorticalBone230, arrayCorticalBoneMLPstartErrorNoTrk);
   TGraph *gMLPErrorVsWAdiposeK = new TGraph(arrayIdxAdipose230, arrayWetWeplAdipose230, arrayAdiposeMLPstartErrorNoTrk);

   gMLPErrorVsWWater->SetMarkerColor(kBlue);
   gMLPErrorVsWA150->SetMarkerColor(kOrange+2);
   gMLPErrorVsWB100->SetMarkerColor(kRed);
   gMLPErrorVsWCorticalBone->SetMarkerColor(kBlack);
   gMLPErrorVsWAdipose->SetMarkerColor(kGreen);
   
   gMLPErrorVsWWaterK->SetMarkerColor(kBlue);
   gMLPErrorVsWA150K->SetMarkerColor(kOrange+2);
   gMLPErrorVsWB100K->SetMarkerColor(kRed);
   gMLPErrorVsWCorticalBoneK->SetMarkerColor(kBlack);
   gMLPErrorVsWAdiposeK->SetMarkerColor(kGreen);

   gMLPErrorVsWWater->SetMarkerStyle(21);
   gMLPErrorVsWA150->SetMarkerStyle(21);
   gMLPErrorVsWB100->SetMarkerStyle(21);
   gMLPErrorVsWAdipose->SetMarkerStyle(21);
   gMLPErrorVsWCorticalBone->SetMarkerStyle(21);
   
   gMLPErrorVsWWaterK->SetMarkerStyle(22);
   gMLPErrorVsWA150K->SetMarkerStyle(22);
   gMLPErrorVsWB100K->SetMarkerStyle(22);
   gMLPErrorVsWAdiposeK->SetMarkerStyle(22);
   gMLPErrorVsWCorticalBoneK->SetMarkerStyle(22);

   gMLPErrorVsWWater->SetTitle(";WET/WEPL;Error |X_{0} - X_{0}^{opt}| [mm]");

   gMLPErrorVsWWater->Draw("AP");
   gMLPErrorVsWA150->Draw("P");
   gMLPErrorVsWB100->Draw("P");
   gMLPErrorVsWCorticalBone->Draw("P");
   gMLPErrorVsWAdipose->Draw("P");
  
   /*
   gMLPErrorVsWWaterK->Draw("P");
   gMLPErrorVsWA150K->Draw("P");
   gMLPErrorVsWB100K->Draw("P");
   gMLPErrorVsWCorticalBoneK->Draw("P");
   gMLPErrorVsWAdiposeK->Draw("P");
*/
   
   TLegend *leg3 = new TLegend(.3, .66, .64, .8655);
   leg3->AddEntry(gMLPErrorVsWWater, "Water", "P");
//   leg3->AddEntry(gMLPErrorVsWCorticalBoneK, "MLP", "P");
   leg3->AddEntry(gMLPErrorVsWAdipose, "Adipose", "P");
   leg3->AddEntry(gMLPErrorVsWB100, "B100 Bone", "P");
   leg3->AddEntry(gMLPErrorVsWA150, "A150 T.E.P.", "P");
   leg3->AddEntry(gMLPErrorVsWCorticalBone, "Cortical Bone", "P");
   leg3->SetTextFont(22);
   leg3->Draw();

   /*
   TLatex *lat = new TLatex();
   lat->SetTextSize(0.04);
   lat->SetTextColor(kRed);
   lat->DrawLatex(0.6, 0.9, "B100 Bone");
   lat->SetTextColor(kOrange+2);
   lat->DrawLatex(0.6, 0.7, "A150 TEP");
   lat->SetTextColor(kBlack);
   lat->DrawLatex(0.6, 0.5, "Cortical Bone");
   lat->SetTextColor(kGreen);
   lat->DrawLatex(0.6, 0.3, "Adipose");
   lat->SetTextColor(kBlue);
   lat->DrawLatex(0.6, 0.1, "Water");
   */

   TCanvas *c5 = new TCanvas("c5", "Errors vs rotation", 600, 600);
   TGraph *gRotationError = new TGraph(idxRotation, arrayRotation, arrayRotationError);

   gRotationError->SetTitle("MLP estimation error with #hat{X_{0}}, 230 MeV beam/160 mm phantom;Beam rotation angle [mrad];Error |X_{0}^{MC} - X_{0}^{est}| [mm]");

   gRotationError->SetMarkerColor(kBlue);
   gRotationError->SetMarkerStyle(21);
   gRotationError->SetMarkerSize(0.8);

   gRotationError->Draw("PA");
   gRotationError->GetYaxis()->SetTitleOffset(1.2);
   
   
   TCanvas *c6 = new TCanvas("c6", "Errors vs rotation", 1000, 800);
   TGraph *gSpotsizeError = new TGraph(idxSpotsize, arraySpotsize, arraySpotsizeError);
   TGraph *gSpotsizeErrorKrah = new TGraph(idxSpotsizeKrah, arraySpotsizeKrah, arraySpotsizeErrorKrah);

   gSpotsizeError->SetTitle("MLP estimation error with #hat{X_{0}}, 200 MeV beam/160 mm phantom;Beam spotsize [mm];Error |X_{0}^{MC} - X_{0}^{est}| [mm]");

   gSpotsizeError->SetMarkerColor(kBlue-7);
   gSpotsizeError->SetMarkerStyle(21);
   gSpotsizeError->SetMarkerSize(1.2);
   
   gSpotsizeErrorKrah->SetMarkerColor(kRed);
   gSpotsizeErrorKrah->SetMarkerStyle(21);
   gSpotsizeErrorKrah->SetMarkerSize(1.2);

   gSpotsizeError->Draw("PA");
   gSpotsizeErrorKrah->Draw("P");
//   gSpotsizeError->GetYaxis()->SetTitleOffset(1.2);

   TLegend *legss = new TLegend(.3, .66, .64, .8655);
   legss->AddEntry(gSpotsizeError, "Linear Projection Model", "P");
   legss->AddEntry(gSpotsizeErrorKrah, "Bayesian MLP", "P");
   legss->Draw();
   legss->SetTextFont(22);

   TCanvas *c7 = new TCanvas("c7", "Different estimation models", 1000, 800);

   TGraph *g200NoTrk = new TGraph(idx230, array230wetwepl, array230NoTrk);
   TGraph *g200Est = new TGraph(idx200, array200wetwepl, array200Est);
   TGraph *g230Est = new TGraph(idx230, array230wetwepl, array230Est);
   TGraph *g200Krah = new TGraph(idx200, array200wetwepl, array200Krah);
   TGraph *g230Krah = new TGraph(idx230, array230wetwepl, array230Krah);

   g200Est->SetMarkerColor(kBlue);
   g200Est->SetMarkerStyle(21);
   g200Est->SetMarkerSize(1.5);
   g200Krah->SetMarkerColor(kBlue-7);
   g200Krah->SetMarkerStyle(21);
   g200Krah->SetMarkerSize(1);
   g230Est->SetMarkerColor(kRed);
   g230Est->SetMarkerStyle(21);
   g230Est->SetMarkerSize(1.5);
   g230Krah->SetMarkerColor(kRed-9);
   g230Krah->SetMarkerStyle(21);
   g230Krah->SetMarkerSize(1);

   g200Est->SetTitle(";WET/WEPL;X_{0} estimation error [mm]");

   g200Est->Draw("AP");
   g200Krah->Draw("P");
   g230Est->Draw("P");
   g230Krah->Draw("P");

//   g200NoTrk->GetYaxis()->SetRangeUser(0, 6);

   TLegend *leg4 = new TLegend(.3, .66, .64, .8655);
   leg4->AddEntry(g230Krah, "MLP 230 MeV", "P");
   leg4->AddEntry(g230Est, "LPM 230 MeV", "P");
   leg4->AddEntry(g200Krah, "MLP 200 MeV", "P");
   leg4->AddEntry(g200Est, "LPM 200 MeV", "P");
   leg4->Draw();
   leg4->SetTextFont(22);

   TCanvas *c8 = new TCanvas("c8", "Error in different materials", 600, 600);

   TGraph *g200WaterEst = new TGraph(idx200, array200Phantomsize, array200Est);
   TGraph *g200WaterKrah = new TGraph(idx200, array200Phantomsize, array200Krah);
   TGraph *g200CorticalBoneEst = new TGraph(idxKrahCorticalBone, arrayCorticalBonekrahDivergence, arrayCorticalBonekrahMLPstartErrorEst);
   TGraph *g200CorticalBoneKrah = new TGraph(idxKrahCorticalBone, arrayCorticalBonekrahDivergence, arrayCorticalBonekrahMLPstartErrorKrah);
   TGraph *g200A150Est = new TGraph(idxKrahA150, arrayA150krahDivergence, arrayA150krahMLPstartErrorEst);
   TGraph *g200A150Krah = new TGraph(idxKrahA150, arrayA150krahDivergence, arrayA150krahMLPstartErrorKrah);

   g200WaterEst->SetMarkerColor(kBlue);
   g200WaterEst->SetMarkerStyle(21);
   g200WaterEst->SetMarkerSize(0.8);
   g200WaterKrah->SetMarkerColor(kBlue);
   g200WaterKrah->SetMarkerStyle(22);
   g200WaterKrah->SetMarkerSize(1.4);
   
   g200CorticalBoneEst->SetMarkerColor(kRed);
   g200CorticalBoneEst->SetMarkerStyle(21);
   g200CorticalBoneEst->SetMarkerSize(0.8);
   g200CorticalBoneKrah->SetMarkerColor(kRed);
   g200CorticalBoneKrah->SetMarkerStyle(22);
   g200CorticalBoneKrah->SetMarkerSize(1.4);
   
   g200A150Est->SetMarkerColor(kBlack);
   g200A150Est->SetMarkerStyle(21);
   g200A150Est->SetMarkerSize(0.8);
   g200A150Krah->SetMarkerColor(kBlack);
   g200A150Krah->SetMarkerStyle(22);
   g200A150Krah->SetMarkerSize(1.4);

   g200WaterEst->SetTitle(";WET/WEPL;Error | X_{0}^{MC} - X_{0}^{est} | [mm]");

   g200WaterEst->Draw("AP");
   g200WaterKrah->Draw("P");
   g200CorticalBoneEst->Draw("P");
   g200CorticalBoneKrah->Draw("P");
   g200A150Est->Draw("P");
   g200A150Krah->Draw("P");
   
   TLegend *leg5 = new TLegend(.3, .66, .64, .8655);
   leg5->AddEntry(g200WaterEst, "Projection Model (Water)", "P");
   leg5->AddEntry(g200CorticalBoneEst, "Projection Model (Cortical Bone)", "P");
   leg5->AddEntry(g200A150Est, "Projection Model (A150)", "P");
   leg5->AddEntry(g200WaterKrah, "Bayesian Model (Water)", "P");
   leg5->AddEntry(g200CorticalBoneKrah, "Bayesian Model (Cortical Bone)", "P");
   leg5->AddEntry(g200A150Krah, "Bayesian Model (A150)", "P");
   leg5->Draw();
   leg5->SetTextFont(22);

   TCanvas *c9 = new TCanvas("c9", "X0 accuracy vs different tracker uncertainties", 1000, 800);

   for (int i=0; i<idx230; i++) {
      array230UncBergenKrahNoModel[i] -= array230Krah[i];
      array230UncBergenKrahModel[i] -= array230Krah[i];
      array230UncBergenEst[i] -= array230Krah[i];
      array230UncLLKrahNoModel[i] -= array230Krah[i];
      array230UncLLKrahModel[i] -= array230Krah[i];
      array230UncLLEst[i] -= array230Krah[i];
   }

   TGraph *g230WaterEst0mrad = new TGraph(idx230, array230Phantomsize, array230Krah);
   TGraph *g230WaterUncLLKrahNoModel = new TGraph(idx230UncLLNoModel, array230Phantomsize10mrad, array230UncLLKrahNoModel);
   TGraph *g230WaterUncLLEst = new TGraph(idx230UncLLModel, array230Phantomsize10mrad, array230UncLLEst);
   TGraph *g230WaterUncLLKrahModel = new TGraph(idx230UncLLModel, array230Phantomsize10mrad66um, array230UncLLKrahModel);
   TGraph *g230WaterUncBergenKrahNoModel = new TGraph(idx23010mrad, array230Phantomsize10mrad, array230UncBergenKrahNoModel);
   TGraph *g230WaterUncBergenEst = new TGraph(idx23010mrad, array230Phantomsize10mrad, array230UncBergenEst);
   TGraph *g230WaterUncBergenKrahModel = new TGraph(idx23010mrad66um, array230Phantomsize10mrad66um, array230UncBergenKrahModel);

   g230WaterUncLLKrahNoModel->SetLineColor(kRed);
   g230WaterUncLLKrahNoModel->SetLineWidth(3);
   g230WaterUncLLEst->SetLineColor(kOrange+2);
   g230WaterUncLLEst->SetLineWidth(3);
   g230WaterUncLLKrahModel->SetLineColor(kBlack);
   g230WaterUncLLKrahModel->SetLineWidth(3);
   
   g230WaterUncBergenKrahNoModel->SetLineColor(kRed);
   g230WaterUncBergenKrahNoModel->SetLineWidth(3);
   g230WaterUncBergenKrahNoModel->SetLineStyle(9);
   g230WaterUncBergenEst->SetLineColor(kOrange+2);
   g230WaterUncBergenEst->SetLineWidth(3);
   g230WaterUncBergenEst->SetLineStyle(9);
   g230WaterUncBergenKrahModel->SetLineColor(kBlack);
   g230WaterUncBergenKrahModel->SetLineWidth(3);
   g230WaterUncBergenKrahModel->SetLineStyle(9);

   g230WaterUncLLEst->SetTitle(";Phantom size [mm]; Difference from ideal detector (with MLP+CSP) [mm]");

//   g230WaterEst0mrad->Draw("AL");
   g230WaterUncLLEst->Draw("AL");
   g230WaterUncLLKrahModel->Draw("L");
   g230WaterUncLLKrahNoModel->Draw("L");
   g230WaterUncBergenEst->Draw("L");
   g230WaterUncBergenKrahModel->Draw("L");
   g230WaterUncBergenKrahNoModel->Draw("L");

   TLegend *leg6 = new TLegend(.3, .66, .64, .8655);
//   leg6->AddEntry(g230WaterEst0mrad, "Ideal trackers MLP", "L")
   leg6->AddEntry(g230WaterUncLLKrahModel, "Modeled LL uncertainties MLP", "L");
   leg6->AddEntry(g230WaterUncLLKrahNoModel, "Unmodeled LL uncertainties MLP", "L");
   leg6->AddEntry(g230WaterUncLLEst, "LL uncertainties LPM", "L");
   leg6->AddEntry(g230WaterUncBergenKrahModel, "Modeled Bergen uncertainties MLP", "L");
   leg6->AddEntry(g230WaterUncBergenKrahNoModel, "Unmodeled Bergen uncertainties MLP", "L");
   leg6->AddEntry(g230WaterUncBergenEst, "Bergen uncertainties LPM", "L");
   leg6->Draw();
   leg6->SetTextFont(22);

   gStyle->SetMarkerStyle(21);
   gStyle->SetMarkerSize(1.2);
   TCanvas *c10 = new TCanvas("c10", "3 sigma estimation vs w", 1500, 1200);
   TGraph *gP2Water200 = new TGraph(idxP2Water200, arrayP2Water200w, arrayP2Water200Fit);
   TGraph *gP2Water230 = new TGraph(idxP2Water230, arrayP2Water230w, arrayP2Water230Fit);
   TGraph *gP2A150230 = new TGraph(idxP2A150230, arrayP2A150230w, arrayP2A150230Fit);
   TGraph *gP2B100230 = new TGraph(idxP2B100230, arrayP2B100230w, arrayP2B100230Fit);
   TGraph *gP2Adipose230 = new TGraph(idxP2Adipose230, arrayP2Adipose230w, arrayP2Adipose230Fit);
   TGraph *gP2CorticalBone230 = new TGraph(idxP2CorticalBone230, arrayP2CorticalBone230w, arrayP2CorticalBone230Fit);

   gP2Water200->SetMarkerColor(kBlue-7);
   gP2Water230->SetMarkerColor(kBlue);
   gP2A150230->SetMarkerColor(kRed);
   gP2CorticalBone230->SetMarkerColor(kBlack);
   gP2Adipose230->SetMarkerColor(kGreen);
   gP2B100230->SetMarkerColor(kOrange+2);

   gP2Water230->SetTitle(";(WET/WEPL)^{2};#mu + 2.5#sigma from Gaussian fit [mrad]");
   gP2Water230->Draw("AP");
    gP2Water200->Draw("P");
   gP2A150230->Draw("P");
   gP2B100230->Draw("P");
   gP2CorticalBone230->Draw("P");
   gP2Adipose230->Draw("P");

   TF1 * gP2fit = new TF1("gP2fit", "pol3");
   gP2Water200->Fit("gP2fit");
   gP2Water200->SetLineWidth(4);
   
   TLegend *leg7 = new TLegend(.3, .66, .64, .8655);
   leg7->AddEntry(gP2Water200, "Water 200 MeV", "P");
   leg7->AddEntry(gP2Water230, "Water 230 MeV", "P");
   leg7->AddEntry(gP2CorticalBone230, "C. B. 230 MeV", "P");
   leg7->AddEntry(gP2A150230, "A150 230 MeV", "P");
   leg7->AddEntry(gP2B100230, "B100 230 MeV", "P");
   leg7->AddEntry(gP2Adipose230, "Adipose 230 MeV", "P");
   leg7->AddEntry(gP2fit, "3^{rd} order pol. fit", "L");
   leg7->SetTextFont(22);
   leg7->Draw();
   
}
