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
#include <TLatex.h>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>
#include <TColor.h>

using namespace std;

Int_t absorberThickness = 2;

Color_t kColorBK = kRed-4;
Color_t kColorUlmer = kCyan+1;
Color_t kColorSpline = kGreen+1;
Color_t kColorLinear = kBlack;

Int_t kLineBK = 1;
Int_t kLineUlmer = 9;
Int_t kLineLinear = 7;
Int_t kLineSpline = 10;

void makePlots() {
   TCanvas *c1 = new TCanvas("c1", "Fit results", 1200, 800);
   TPad *pad1 = new TPad("pad1", "70", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad2 = new TPad("pad2", "30", 0.0, 0.0, 1.0, 0.3, 0);
   pad1->Draw();
   pad2->Draw();
   
   TCanvas *c2 = new TCanvas("c2", "Correct Tracks fraction focal", 1200, 800);
   TCanvas *c22 = new TCanvas("c22", "Correct Tracks fraction PB", 1200, 900);
   TCanvas *c222 = new TCanvas("c222", "Correct Tracks fraction uniform", 1200, 900);
   TCanvas *c3 = new TCanvas("c3", "Reconstruction efficiency", 1200, 800);
   TCanvas *c4 = new TCanvas("c4", "Chip alignment", 1200, 800);
   TCanvas *c5 = new TCanvas("c5", "Chip sensitivity calibration", 1200, 800);
   TCanvas *c6 = new TCanvas("c6", "Resolution", 1200, 800);
   c6->Divide(3,1,0.0001,0.0001);
   TCanvas *c7 = new TCanvas("c7", "Parameterization accuracy", 1200, 900);

   Float_t  arrayE[500] = {0}; // energy MC
   Float_t  arrayMCActualSigma[500] = {0}; // Measured range straeggling from full MC
   Float_t  arrayMCActualSigmaRatio[500] = {0}; // Measured range straeggling from full MC
   Float_t  arrayMCActualResidualRange[500] = {0}; // Measured range straeggling from full MC
   Float_t  arrayEE[500] = {0}; // error on energy MC
   Float_t  arrayMC[500] = {0}; // range MC
   Float_t  arrayEMC[500] = {0}; // error on range MC
   Float_t  arrayEMCRatio[500] = {0}; // error on range MC
   Float_t  arrayEMCSub[500] = {0}; // error on range MC
   Float_t  arrayEMCSubRatio[500] = {0}; // error on range MC
   Float_t  arrayMCDelta[500] = {0}; // range MC
   Float_t  arrayPSTAR[500] = {0};
   Float_t  arrayPSTARDelta[500] = {0};
   Float_t  arrayPSTARshade[400] = {0};
   Float_t  arrayPSTARmin[500] = {0};
   Float_t  arrayPSTARmax[500] = {0};
   Float_t  arrayPSTARminDelta[500] = {0};
   Float_t  arrayPSTARmaxDelta[500] = {0};
   Float_t  arrayPSTARmaxDeltaRatio[500] = {0};
   Float_t  arrayEPstar[500] = {0};
   Float_t  arrayE2[500] = {0}; // energy data 
   Float_t  arrayEE2[500] = {0}; // error on energy data
   Float_t  arrayEData[500] = {0};  // range data
   Float_t  arrayData[500] = {0}; // error on range data
   Float_t  arrayEfficiencyEnergy[300];
   Float_t  arrayEfficiencyFinal[300];
   Float_t  arrayEfficiencyFinal2[300];
   Float_t  alignmentChipXOrig[96] = {0};
   Float_t  alignmentChipYOrig[96] = {0};
   Float_t  alignmentChipXMine[96] = {0};
   Float_t  alignmentChipYMine[96] = {0};
   Float_t  alignmentChipsMine[96] = {0};
   Float_t  alignmentChipsOrig[96] = {0};
   Float_t  chipCalibrationChip150[27] = {0};
   Float_t  chipCalibrationChip160[27] = {0};
   Float_t  chipCalibrationChip170[27] = {0};
   Float_t  chipCalibrationChip180[27] = {0};
   Float_t  chipCalibrationChip188[27] = {0};
   Float_t  chipCalibrationFactor150[27] = {0};
   Float_t  chipCalibrationFactor160[27] = {0};
   Float_t  chipCalibrationFactor170[27] = {0};
   Float_t  chipCalibrationFactor180[27] = {0};
   Float_t  chipCalibrationFactor188[27] = {0};
   Float_t  chipCalibrationNumber150[27] = {0};
   Float_t  chipCalibrationNumber160[27] = {0};
   Float_t  chipCalibrationNumber170[27] = {0};
   Float_t  chipCalibrationNumber180[27] = {0};
   Float_t  chipCalibrationNumber188[27] = {0};
   
   Float_t  chipCalibrationError150[27] = {0};
   Float_t  chipCalibrationError160[27] = {0};
   Float_t  chipCalibrationError170[27] = {0};
   Float_t  chipCalibrationError180[27] = {0};
   Float_t  chipCalibrationError188[27] = {0};

   // Parameterization (l = 1st quartile, m = median, h = 3rd quartile)
   Double_t paraNPoints[100] = {};
   Double_t paraBKl[100] = {};
   Double_t paraBKm[100] = {};
   Double_t paraBKh[100] = {};
   Double_t paraBKInvl[100] = {};
   Double_t paraBKInvm[100] = {};
   Double_t paraBKInvh[100] = {};
   Double_t paraUlmerl[100] = {};
   Double_t paraUlmerm[100] = {};
   Double_t paraUlmerh[100] = {};
   Double_t paraUlmerInvl[100] = {};
   Double_t paraUlmerInvm[100] = {};
   Double_t paraUlmerInvh[100] = {};
   Double_t paraSplinel[100] = {};
   Double_t paraSplinem[100] = {};
   Double_t paraSplineh[100] = {};
   Double_t paraSplineInvl[100] = {};
   Double_t paraSplineInvm[100] = {};
   Double_t paraSplineInvh[100] = {};
   Double_t paraLinearl[100] = {};
   Double_t paraLinearm[100] = {};
   Double_t paraLinearh[100] = {};
   Double_t paraLinearInvl[100] = {};
   Double_t paraLinearInvm[100] = {};
   Double_t paraLinearInvh[100] = {};

   Int_t nThisEnergy = 0, lastEnergy = 0;
   Float_t  mmAbsorbator;

   gStyle->SetOptStat(0);

   Float_t nomrange_, estrange_, sigmaRange_, lastRange_, nomsigma_, waterphantomthickness_, dummy0;
   Int_t energy_, thickness_;
   Float_t estimatedStraggling;

   ifstream in0;
   in0.open("OutputFiles/findManyRangesDegrader.csv");
   Int_t nlines0 = 0;
   Float_t a_dtc = 0, p_dtc = 0;
   Float_t a_wtr = 0.02387;
   Float_t p_wtr = 1.7547;

   if       (absorberThickness == 2) {
      a_dtc = 0.012511;
      p_dtc = 1.730529;
   }
   else if  (absorberThickness == 3) {
      a_dtc = 0.010746;
      p_dtc = 1.758228;
   }

   Float_t wepl_ratio0 = a_wtr / a_dtc * pow(250 / a_wtr, 1 - p_dtc / p_wtr);
//   Float_t wtr_range = a_wtr * pow(250, p_wtr);
   Float_t wtr_range = 379.4; // PSTAR

   while (1) {
      in0 >>  waterphantomthickness_ >> thickness_ >> nomrange_ >> nomsigma_ >> dummy0 >> dummy0 >> dummy0;

      if (!in0.good()) {
         break;
      }

      if (thickness_ != absorberThickness) {
         continue;
      }

      arrayMCActualSigma[nlines0] = nomsigma_ * wepl_ratio0;
      arrayMCActualSigmaRatio[nlines0] = nomsigma_ * 100 * wepl_ratio0 / wtr_range;
      arrayMCActualResidualRange[nlines0++] = nomrange_ * wepl_ratio0;
   }
   in0.close();

   ifstream in;
   in.open("OutputFiles/result_makebraggpeakfit.csv");

   cout << "Opened file.\n";

   Int_t nlines = 0;
   TNtuple *ntuple = new TNtuple("ntuple", "data from file", "energy_:nomrange_:estrange_:sigmaRange:lastRange_");
   
   Int_t MC2Data = -1;

   Float_t meanError = 0;
   Float_t meanAbsError = 0;
   Float_t meanSigma = 0;
   Float_t aprime_dtc = 0;
   Float_t aprime_wtr = 0.0087; // MeV^2 / mm

   while (1) {
      in >> thickness_ >> energy_ >> nomrange_ >> estrange_ >> nomsigma_ >> sigmaRange_;

      if (!in.good()) {
         break;
      }

      if (thickness_ != absorberThickness) continue;

      mmAbsorbator = thickness_;
      meanError += ( estrange_ - nomrange_ ) / nomrange_;
      meanAbsError += fabs(( estrange_ - nomrange_ ) / nomrange_);
      meanSigma += sigmaRange_;

      
      if (mmAbsorbator == 2) {
         a_dtc = 0.0096;
         p_dtc = 1.784;
         aprime_dtc = 0.01938;
      }
      
      else if (mmAbsorbator == 3) {
         a_dtc = 0.010746;
         p_dtc = 1.758228;
         aprime_dtc = 0.01971;
      }

      else if (mmAbsorbator == 4) {
         a_dtc = 0.0117;
         p_dtc = 1.7450;
         aprime_dtc = 0.01988;
      }

      Float_t zprime = energy_;
      
      Float_t wtr_term = aprime_wtr * (pow(p_wtr, 2) * pow(a_wtr, 2/p_wtr)) / (3 - 2/p_wtr);
      Float_t wtr_strag = wtr_term * (pow(wtr_range, 3-2/p_wtr) - pow(wtr_range - zprime, 3-2/p_wtr));
     
      Float_t wepl_ratio = a_wtr / a_dtc * pow(nomrange_ / a_wtr, 1 - p_dtc / p_wtr);
      Float_t dtc_term = aprime_dtc * (pow(p_dtc, 2) * pow(a_dtc, 2/p_dtc)) / (3 - 2/p_dtc);
      Float_t dtc_strag = pow(wepl_ratio, 2) * dtc_term * pow(nomrange_/wepl_ratio, 3-2/p_dtc);

      /*
      printf("nomrange_ = %.2f. Calculated WEPL range = %.2f.\n", nomrange_, wtr_range - zprime);
      printf("Bad method of calculating DTC straggling = %.2f mm.\n", sqrt(wtr_term * pow(nomrange_, 3-2/p_wtr)));
      printf("Using a %.2f mm water phantom, the estimated straggling from water is %.2f mm and from DTC is %.2f mm. WEPL ratio is %.2f.\n", zprime, sqrt(wtr_strag), sqrt(dtc_strag), wepl_ratio);
      */

      estimatedStraggling = sqrt(dtc_strag + wtr_strag);   
      estimatedStraggling = nomsigma_; 

      if (nlines < MC2Data || MC2Data<0) {
         arrayE[nlines] = energy_;
         arrayEE[nlines] = 0;
         arrayMC[nlines] = estrange_;
         arrayMCDelta[nlines] = estrange_ - nomrange_;
         arrayEMC[nlines] = sigmaRange_;
         arrayEMCRatio[nlines] = sigmaRange_ / wtr_range * 100; // used wtr_range
         arrayEMCSub[nlines] = sqrt(abs(pow(estimatedStraggling, 2) - pow(sigmaRange_, 2)));
         arrayEMCSubRatio[nlines] = sqrt(abs(pow(estimatedStraggling, 2) - pow(sigmaRange_, 2))) / nomrange_ * 100;
         arrayPSTAR[nlines] = nomrange_;
         arrayPSTARmin[nlines] = nomrange_ - estimatedStraggling;
         arrayPSTARmax[nlines] = nomrange_ + estimatedStraggling;
         arrayPSTARDelta[nlines] = 0;
         arrayPSTARminDelta[nlines] = -estimatedStraggling;
         arrayPSTARmaxDelta[nlines] = estimatedStraggling;
         arrayPSTARmaxDeltaRatio[nlines] = estimatedStraggling / wtr_range * 100; // used wtr_range
         arrayEPstar[nlines] = 0;
      }

      else {
         arrayE2[nlines-MC2Data] = energy_;
         arrayEE2[nlines-MC2Data] = 0;
         arrayData[nlines-MC2Data] = estrange_;
         arrayEData[nlines-MC2Data] = sigmaRange_;
      }

      nlines++;
   }

   if (MC2Data<0) {
      MC2Data = nlines;
   }

   meanError /= nlines;
   meanSigma /= nlines;

   cout << "Mean error on fit range is " << meanError << " mm.\n";
   cout << "Mean | error | on fit range is " << meanAbsError << " mm.\n";
   cout << "Mean SIGMA on fit range is " << meanSigma << " mm.\n";

   in.close();
   
   ifstream in2;
   ifstream in2uni;
//   in2.open("OutputFiles/lastLayerCorrect_different_nRuns_elastic_noinelastic.csv");
   in2.open("OutputFiles/lastLayerCorrect_different_nRuns_ProtonHelium_5cm_2.csv");
   in2uni.open("OutputFiles/lastLayerCorrect_different_nRuns_uniform.csv");
   Float_t np, correctLast, correctWhole, lastIsFirst, lastIsAlmostFirst;
   Int_t   mmAbsorber_;
   Float_t arrayFractionX[200] = {0};
   Float_t arrayFractionY[200] = {0};
   Float_t arrayFractionY2[200] = {0};
   Float_t arrayFractionY3[200] = {0};
   
   Float_t arrayFraction2mmX[200] = {0};
   Float_t arrayFraction2mmY[200] = {0};
   Float_t arrayFraction2mmY2[200] = {0};
   Float_t arrayFraction3mmX[200] = {0};
   Float_t arrayFraction3mmY[200] = {0};
   Float_t arrayFraction3mmY2[200] = {0};
   Float_t arrayFraction35mmX[200] = {0};
   Float_t arrayFraction35mmY[200] = {0};
   Float_t arrayFraction4mmX[200] = {0};
   Float_t arrayFraction4mmY[200] = {0};
   Float_t arrayFraction5mmX[200] = {0};
   Float_t arrayFraction5mmY[200] = {0};
   Float_t arrayFraction6mmX[200] = {0};
   Float_t arrayFraction6mmY[200] = {0};
   Float_t arrayFraction2mmNoDiffusionX[200] = {0};
   Float_t arrayFraction2mmNoDiffusionY[200] = {0};
   Float_t arrayFraction3mmNoDiffusionX[200] = {0};
   Float_t arrayFraction3mmNoDiffusionY[200] = {0};
   Float_t arrayFraction35mmNoDiffusionX[200] = {0};
   Float_t arrayFraction35mmNoDiffusionY[200] = {0};
   Float_t arrayFraction4mmNoDiffusionX[200] = {0};
   Float_t arrayFraction4mmNoDiffusionY[200] = {0};
   Float_t arrayFraction5mmNoDiffusionX[200] = {0};
   Float_t arrayFraction5mmNoDiffusionY[200] = {0};
   Float_t arrayFraction6mmNoDiffusionX[200] = {0};
   Float_t arrayFraction6mmNoDiffusionY[200] = {0};
   Float_t arrayFractionUniformX[200] = {0};
   Float_t arrayFractionUniformY[200] = {0};

   Int_t nlinesF = 0, nlinesF2 = 0, nlinesF3 = 0, nlinesF4 = 0, nlinesF5 = 0, nlinesF6 = 0, nlinesF35 = 0;
   Int_t nlinesNoDiffusionF = 0, nlinesNoDiffusionF2 = 0, nlinesNoDiffusionF3 = 0, nlinesNoDiffusionF4 = 0, nlinesNoDiffusionF5 = 0, nlinesNoDiffusionF6 = 0, nlinesNoDiffusionF35 = 0;
   Int_t nlinesUniform = 0;
   
   Float_t lastIsFirstAllTracks, lastIsFirstAllTracksAfterFiltering, lastIsFirst2ndOK;

   Int_t dropIdx = 0;

   while (1) {
      in2 >> mmAbsorber_ >> np >> lastIsFirstAllTracks >> lastIsFirstAllTracksAfterFiltering;

      /*
      if (np < 50) {
         if (dropIdx++ % 5 != 0) continue;
      }
      */

      if (!in2.good()) break;
     
      if (mmAbsorber_ == 0) {  // p 5
         arrayFraction2mmX[nlinesF2] = np;
         arrayFraction2mmY[nlinesF2] = lastIsFirstAllTracks * 100;
         arrayFraction2mmY2[nlinesF2++] = lastIsFirstAllTracksAfterFiltering * 100;
      }

      if (mmAbsorber_ == 1) { // he 5
         arrayFraction3mmX[nlinesF3] = np;
         arrayFraction3mmY[nlinesF3] = lastIsFirstAllTracks * 100;
         arrayFraction3mmY2[nlinesF3++] = lastIsFirstAllTracksAfterFiltering * 100;
      }
      if (mmAbsorber_ == 2) { // p 16
         arrayFraction4mmX[nlinesF4] = np;
         arrayFraction4mmY[nlinesF4++] = lastIsFirstAllTracksAfterFiltering * 100;
      }

      if (mmAbsorber_ == 3) { // he 16
         arrayFraction5mmX[nlinesF5] = np;
         arrayFraction5mmY[nlinesF5++] = lastIsFirstAllTracksAfterFiltering * 100;
      }
/*

      if (mmAbsorber_ == 6) {
         arrayFraction6mmX[nlinesF6] = np;
         arrayFraction6mmY[nlinesF6++] = lastIsFirstAllTracks * 100;
      }
      
      if (mmAbsorber_ == 35) {
         arrayFraction35mmX[nlinesF35] = np;
         arrayFraction35mmY[nlinesF35++] = lastIsFirstAllTracks * 100;
      }
      */
   }

   in2.close();
   
   in2.open("OutputFiles/lastLayerCorrect_different_nRuns_nodiffusion.csv");
   while (1) {
      in2 >> mmAbsorber_ >> np >> lastIsFirstAllTracks;

      /*
      if (np < 50) {
         if (dropIdx++ % 5 != 0) continue;
      }
      */

      if (!in2.good()) break;
     
      if (mmAbsorber_ == 2) { 
         arrayFraction2mmNoDiffusionX[nlinesNoDiffusionF2] = np;
         arrayFraction2mmNoDiffusionY[nlinesNoDiffusionF2++] = lastIsFirstAllTracks * 100;
      }

      if (mmAbsorber_ == 3) {
         arrayFraction3mmNoDiffusionX[nlinesNoDiffusionF3] = np;
         arrayFraction3mmNoDiffusionY[nlinesNoDiffusionF3++] = lastIsFirstAllTracks * 100;
      }

      if (mmAbsorber_ == 4) {
         arrayFraction4mmNoDiffusionX[nlinesNoDiffusionF4] = np;
         arrayFraction4mmNoDiffusionY[nlinesNoDiffusionF4++] = lastIsFirstAllTracks * 100;
      }

      if (mmAbsorber_ == 5) {
         arrayFraction5mmNoDiffusionX[nlinesNoDiffusionF5] = np;
         arrayFraction5mmNoDiffusionY[nlinesNoDiffusionF5++] = lastIsFirstAllTracks * 100;
      }

      if (mmAbsorber_ == 6) {
         arrayFraction6mmNoDiffusionX[nlinesNoDiffusionF6] = np;
         arrayFraction6mmNoDiffusionY[nlinesNoDiffusionF6++] = lastIsFirstAllTracks * 100;
      }
      
      if (mmAbsorber_ == 35) {
         arrayFraction35mmX[nlinesNoDiffusionF35] = np;
         arrayFraction35mmY[nlinesNoDiffusionF35++] = lastIsFirstAllTracks * 100;
      }
   }

   in2.close();

   
   while (1) {
      in2uni >> mmAbsorber_ >> np >> correctWhole >> lastIsFirst >> lastIsAlmostFirst >> lastIsFirstAllTracks;

      if (!in2uni.good()) break;
      arrayFractionUniformX[nlinesUniform] = np;
      arrayFractionUniformY[nlinesUniform++] = lastIsFirstAllTracks * 100;
   }
   in2uni.close();

   ifstream in3;
   in3.open("OutputFiles/efficiency_500.csv");
   Int_t energy, nRecon, nNotLeaving, nFinal;
   Int_t nEnergies = 0;

   while (1) {
      in3 >> energy >> np >> nRecon >> nNotLeaving >> nFinal;

      if (!in3.good()) break;
      if (!lastEnergy) {
         arrayEfficiencyEnergy[0] = energy;
         arrayEfficiencyFinal[0] = nFinal;
         nThisEnergy = np;
         lastEnergy = energy;
      }

      else if (lastEnergy == energy) {
         arrayEfficiencyFinal[nEnergies] += nFinal;
         nThisEnergy += np;
         lastEnergy = energy;
      }

      else if (lastEnergy != energy) {
         arrayEfficiencyFinal[nEnergies] /= nThisEnergy;
         nEnergies++;
         arrayEfficiencyEnergy[nEnergies] = energy;
         arrayEfficiencyFinal[nEnergies] = nFinal;
         nThisEnergy = np;
         lastEnergy = energy;
      }
   }

   in3.close();

   ifstream in4;
   in4.open("OutputFiles/efficiency_200.csv");
   nEnergies = 0;
   lastEnergy = 0;
   while (1) {
      in4 >> energy >> np >> nRecon >> nNotLeaving >> nFinal;

      if (!in4.good()) break;
      if (!lastEnergy) {
         arrayEfficiencyFinal2[0] = nFinal;
         nThisEnergy = np;
         lastEnergy = energy;
      }

      else if (lastEnergy == energy) {
         arrayEfficiencyFinal2[nEnergies] += nFinal;
         nThisEnergy += np;
         lastEnergy = energy;
      }

      else if (lastEnergy != energy) {
         arrayEfficiencyFinal2[nEnergies] /= nThisEnergy;
         nEnergies++;
         arrayEfficiencyFinal2[nEnergies] = nFinal;
         nThisEnergy = np;
         lastEnergy = energy;
      }
   }

   in4.close();

   ifstream in5;
   in5.open("Data/ExperimentalData/Alignment.txt");
   Int_t chip, nMine;
   Float_t deltaX, deltaY, theta;
   while (1) {
      in5 >> chip >> deltaX >> deltaY >> theta;
      if (!in5.good()) break;

      alignmentChipXOrig[chip] = deltaX * 10000;
      alignmentChipYOrig[chip] = deltaY * 10000;
      alignmentChipsMine[chip] = chip;// - 0.1;
      alignmentChipsOrig[chip] = chip;// + 0.1;
   }

   in5.close();

   ifstream in6;
   in6.open("Data/ExperimentalData/Alignment_mine.txt");
   while (1) {
      in6 >> chip >> deltaX >> deltaY >> theta;
      if (!in6.good()) break;

      alignmentChipXMine[chip] = deltaX * 10000;
      alignmentChipYMine[chip] = deltaY * 10000;
      nMine++;
   }
   in6.close();
      
   // chip 11 is dead
   alignmentChipXOrig[11] = 1e5;
   alignmentChipYOrig[11] = 1e5;
   alignmentChipXMine[11] = 1e5;
   alignmentChipYMine[11] = 1e5;

   ifstream in7;
   Int_t    chipNumber;
   Float_t  calibration150;
   Float_t  calibration160;
   Float_t  calibration170;
   Float_t  calibration180;
   Float_t  calibration188;
   Int_t    numberOfClusters150;
   Int_t    numberOfClusters160;
   Int_t    numberOfClusters170;
   Int_t    numberOfClusters180;
   Int_t    numberOfClusters188;
   Float_t  error150;
   Float_t  error160;
   Float_t  error170;
   Float_t  error180;
   Float_t  error188;
   Float_t  calibrationErrorX[27] = {0};

   in7.open("Data/ExperimentalData/calibration_factors2.txt");
   while (1) {
      in7 >> chipNumber >> calibration150 >> numberOfClusters150 >> error150 >> calibration160 >> numberOfClusters160 >> error160 >> calibration170 >> numberOfClusters170 >> error170 >> calibration180 >> numberOfClusters180 >> error180 >> calibration188 >> numberOfClusters188 >> error188;

      if (!in7.good()) break;

      cout << Form("Found: Chipnumber %d, calibration 150 %.2f, calibration160 %.2f, calibration170 %.2f, calibration180 %.2f, calibration188 %.2f.\n", chipNumber, calibration150, calibration160, calibration170, calibration180, calibration188);

      chipCalibrationNumber150[chipNumber] = numberOfClusters150;
      chipCalibrationNumber160[chipNumber] = numberOfClusters160;
      chipCalibrationNumber170[chipNumber] = numberOfClusters170;
      chipCalibrationNumber180[chipNumber] = numberOfClusters180;
      chipCalibrationNumber188[chipNumber] = numberOfClusters188;

      chipCalibrationChip150[chipNumber] = chipNumber - 0.4;
      chipCalibrationChip160[chipNumber] = chipNumber - 0.2;
      chipCalibrationChip170[chipNumber] = chipNumber;
      chipCalibrationChip180[chipNumber] = chipNumber + 0.2;
      chipCalibrationChip188[chipNumber] = chipNumber + 0.4;

      if (calibration150 == 0) calibration150 = -5;
      if (calibration160 == 0) calibration160 = -5;
      if (calibration170 == 0) calibration170 = -5;
      if (calibration180 == 0) calibration180 = -5;
      if (calibration188 == 0) calibration188 = -5;

      chipCalibrationFactor150[chipNumber] = calibration150;
      chipCalibrationFactor160[chipNumber] = calibration160;
      chipCalibrationFactor170[chipNumber] = calibration170;
      chipCalibrationFactor180[chipNumber] = calibration180;
      chipCalibrationFactor188[chipNumber] = calibration188;

      chipCalibrationError150[chipNumber] = error150;
      chipCalibrationError160[chipNumber] = error160;
      chipCalibrationError170[chipNumber] = error170;
      chipCalibrationError180[chipNumber] = error180;
      chipCalibrationError188[chipNumber] = error188;
   }

   ifstream in8;
   Double_t    pN;
   Int_t    pi = 0;
   Double_t pBKl, pBKm, pBKh, pBKIl, pBKIm, pBKIh;
   Double_t pUlmerl, pUlmerm, pUlmerh, pUlmerIl, pUlmerIm, pUlmerIh;
   Double_t pLinearl, pLinearm, pLinearh, pLinearIl, pLinearIm, pLinearIh;
   Double_t pSplinel, pSplinem, pSplineh, pSplineIl, pSplineIm, pSplineIh;
   in8.open("OutputFiles/MedianValuesForParameterization.csv");
   while (1) {
      in8 >> pN >> pBKl >> pBKm >> pBKh >> pBKIl >> pBKIm >> pBKIh >> pUlmerl >> pUlmerm >> pUlmerh >> pUlmerIl >> pUlmerIm >> pUlmerIh >> pSplinel >> pSplinem >> pSplineh >> pSplineIl >> pSplineIm >> pSplineIh >> pLinearl >> pLinearm >> pLinearh >> pLinearIl >> pLinearIm >> pLinearIh;

      if (!in8.good()) break;

      paraNPoints[pi] = pN;
      paraBKm[pi] = pBKm;
      paraBKl[pi] = pBKm - pBKl;
      paraBKh[pi] = pBKh - pBKm;
      paraBKInvm[pi] = pBKIm;
      paraBKInvl[pi] = pBKIm - pBKIl;
      paraBKInvh[pi] = pBKIh - pBKIm;
      paraUlmerm[pi] = pUlmerm;
      paraUlmerl[pi] = pUlmerm - pUlmerl;
      paraUlmerh[pi] = pUlmerh - pUlmerm;
      paraUlmerInvm[pi] = pUlmerIm;
      paraUlmerInvl[pi] = pUlmerIm - pUlmerIl;
      paraUlmerInvh[pi] = pUlmerIh - pUlmerIm;
      paraSplinem[pi] = pSplinem;
      paraSplinel[pi] = pSplinem - pSplinel;
      paraSplineh[pi] = pSplineh - pSplinem;
      paraSplineInvm[pi] = pSplineIm;
      paraSplineInvl[pi] = pSplineIm - pSplineIl;
      paraSplineInvh[pi] = pSplineIh - pSplineIm;
      paraLinearm[pi] = pLinearm;
      paraLinearl[pi] = pLinearm - pLinearl;
      paraLinearh[pi] = pLinearh - pLinearm;
      paraLinearInvm[pi] = pLinearIm;
      paraLinearInvl[pi] = pLinearIm - pLinearIl;
      paraLinearInvh[pi] = pLinearIh - pLinearIm;

      pi++;
   }
   
   printf("pi = %d\n", pi);


   TGraphErrors *cal150 = new TGraphErrors(27, chipCalibrationChip150, chipCalibrationFactor150, calibrationErrorX, chipCalibrationError150);
   TGraphErrors *cal160 = new TGraphErrors(27, chipCalibrationChip160, chipCalibrationFactor160, calibrationErrorX, chipCalibrationError160);
   TGraphErrors *cal170 = new TGraphErrors(27, chipCalibrationChip170, chipCalibrationFactor170, calibrationErrorX, chipCalibrationError170);
   TGraphErrors *cal180 = new TGraphErrors(27, chipCalibrationChip180, chipCalibrationFactor180, calibrationErrorX, chipCalibrationError180);
   TGraphErrors *cal188 = new TGraphErrors(27, chipCalibrationChip188, chipCalibrationFactor188, calibrationErrorX, chipCalibrationError188);

   pad1->cd();
   
   TGraphErrors *hMC = new TGraphErrors(MC2Data, arrayE, arrayMC, arrayEE, arrayEMC);
   TGraphErrors *hData = new TGraphErrors(nlines-MC2Data, arrayE2, arrayData, arrayEE2, arrayEData);
//   TGraphErrors *pstar = new TGraphErrors(MC2Data, arrayE, arrayPSTAR, arrayEE, arrayEPstar);
   TGraph *pstar = new TGraph(MC2Data, arrayE, arrayPSTAR);
   TGraph *pstarmin = new TGraph(MC2Data, arrayE, arrayPSTARmin);
   TGraph *pstarmax = new TGraph(MC2Data, arrayE, arrayPSTARmax);
   TGraph *pstarshade = new TGraph(MC2Data*2);
   
   for (Int_t i=0; i<MC2Data; i++) {
      pstarshade->SetPoint(i, arrayE[i], arrayPSTARmax[i]);
      pstarshade->SetPoint(MC2Data+i, arrayE[MC2Data-i-1], arrayPSTARmin[MC2Data-i-1]);
   }

   pstarshade->GetXaxis()->SetTitleFont(22);
   pstarshade->GetYaxis()->SetTitleFont(22);
   pstarshade->GetXaxis()->SetTitleOffset(0.9);
   pstarshade->GetYaxis()->SetTitleOffset(0.9);
   pstarshade->GetXaxis()->SetLabelFont(22);
   pstarshade->GetXaxis()->SetTitleSize(0.05);
   pstarshade->GetXaxis()->SetLabelSize(0.05);
   pstarshade->GetYaxis()->SetLabelFont(22);
   pstarshade->GetYaxis()->SetTitleSize(0.05);
   pstarshade->GetYaxis()->SetLabelSize(0.05);

   hMC->SetMarkerColor(kBlue);
   hMC->SetMarkerStyle(21);
   hMC->SetMarkerSize(0.8);

   hData->SetMarkerColor(kRed);
   hData->SetMarkerStyle(22);
   hData->SetMarkerSize(2);
   hData->SetLineColor(kBlack);
   
   pstar->SetLineWidth(3);
   pstar->SetLineColor(kMagenta-10);
   pstarshade->SetTitle(Form("Reconstructed ranges of proton beams with %d mm Al absorbator;Energy [MeV];Reconstructed WET range [mm]", mmAbsorbator));
   
   pstarshade->SetFillColor(kMagenta-10);
   pstarshade->Draw("FA");
   pstarmin->Draw("L");
   pstarmax->Draw("L");

   hMC->SetTitle("Reconstructed ranges #LT#hat{R_{0}}#GT of proton tracks;Energy [MeV];Projected range [mm]");
   hMC->Draw("P");
   hData->Draw("P");


   TLegend *leg = new TLegend(0.15, 0.68, 0.40, 0.85);
   leg->SetTextSize(0.03);
   leg->SetTextFont(22);
   leg->AddEntry(pstar, "Range straggling from MC", "L");
   leg->AddEntry(hMC, "Monte Carlo", "PE");
//   leg->AddEntry(hData, "Experimental data", "PE");
   leg->Draw();

   pad1->Update();
   
   pad2->cd();

   TGraphErrors *hMCD = new TGraphErrors(MC2Data, arrayE, arrayMCDelta, arrayEE, arrayEMC);
   TGraph *pstarD = new TGraph(MC2Data, arrayE, arrayPSTARDelta);
   TGraph *pstarminD = new TGraph(MC2Data, arrayE, arrayPSTARminDelta);
   TGraph *pstarmaxD = new TGraph(MC2Data, arrayE, arrayPSTARmaxDelta);
   TGraph *pstarshadeD = new TGraph(MC2Data*2);

   for (Int_t i=0; i<MC2Data; i++) {
      pstarshadeD->SetPoint(i, arrayE[i], arrayPSTARmaxDelta[i]);
      pstarshadeD->SetPoint(MC2Data+i, arrayE[MC2Data-i-1], arrayPSTARminDelta[MC2Data-i-1]);
   }

   pstarshadeD->GetXaxis()->SetTitleFont(22);
   pstarshadeD->GetYaxis()->SetTitleFont(22);
   pstarshadeD->GetXaxis()->SetTitleOffset(0.9);
   pstarshadeD->GetYaxis()->SetTitleOffset(0.9);
   pstarshadeD->GetXaxis()->SetLabelFont(22);
   pstarshadeD->GetXaxis()->SetTitleSize(0.05);
   pstarshadeD->GetXaxis()->SetLabelSize(0.05);
   pstarshadeD->GetYaxis()->SetLabelFont(22);
   pstarshadeD->GetYaxis()->SetTitleSize(0.05);
   pstarshadeD->GetYaxis()->SetLabelSize(0.05);

   pstarshadeD->GetYaxis()->SetRangeUser(-10, 10); // 145 - 270

   hMCD->SetMarkerColor(kBlue);
   hMCD->SetMarkerStyle(21);
   hMCD->SetMarkerSize(1);

   pstarD->SetLineWidth(3);
   pstarD->SetLineColor(kMagenta-10);
   pstarshadeD->SetTitle(Form("Reconstructed ranges of proton beams with %d mm Al absorbator;Energy [MeV];WET error [mm]", mmAbsorbator));
   
   pstarshadeD->SetFillColor(kMagenta-10);
   pstarshadeD->Draw("FA");
   pstarminD->Draw("L");
   pstarmaxD->Draw("L");
   
   pad2->Update();
   TLine *zeroLine = new TLine(pad2->GetUxmin(), 0, pad2->GetUxmax(), 0);
   zeroLine->Draw();

   hMCD->SetTitle("Reconstructed ranges #LT#hat{R_{0}}#GT of proton tracks;WET error [mm];Projected range [mm]");
   hMCD->Draw("P");

   c6->cd(1);
   TGraph *gResolution = new TGraph(MC2Data, arrayMC, arrayEMC);
   gResolution->SetTitle(Form("WEPL resolution using %d mm Al absorber;Depth in detector [WEPL mm];#Delta WEPL [mm]", mmAbsorbator));
   gResolution->GetYaxis()->SetRangeUser(0,7);
   gResolution->SetMarkerColor(kBlue);
   gResolution->SetMarkerStyle(21);
   gResolution->SetMarkerSize(1);
   gResolution->Draw("PA");
   
   TGraph *gResolutionSub = new TGraph(MC2Data, arrayMC, arrayEMCSub);
   gResolutionSub->SetMarkerColor(kRed);
   gResolutionSub->SetMarkerStyle(21);
   gResolutionSub->SetMarkerSize(1);
//   gResolutionSub->Draw("P");
   
//   TGraph *gResolutionStraggling = new TGraph(MC2Data, arrayMC, arrayPSTARmaxDelta);
   TGraph *gResolutionStraggling = new TGraph(nlines0, arrayMCActualResidualRange, arrayMCActualSigma);
   gResolutionStraggling->SetLineWidth(2);
   gResolutionStraggling->SetLineColor(kRed);
   gResolutionStraggling->Draw("L");

   TLegend *legRes = new TLegend(0.17, 0.77, 0.69, 0.88);
   legRes->SetTextSize(0.03);
   legRes->SetTextFont(22);
   legRes->AddEntry(gResolution, "Measured range spread", "P");
//   legRes->AddEntry(gResolutionSub, "#Delta WEPL - Range straggling", "P");
   legRes->AddEntry(gResolutionStraggling, "MC truth range straggling", "L");
   legRes->Draw();

   c6->cd(2);
   TGraph *gResolutionRatio = new TGraph(MC2Data, arrayMC, arrayEMCRatio);
   gResolutionRatio->SetTitle(Form("WEPL resolution using %d mm Al absorber;Depth in detector [WEPL mm];#Delta WEPL / WEPL [%]", mmAbsorbator));
   gResolutionRatio->GetYaxis()->SetTitleOffset(1.5);
   gResolutionRatio->SetMarkerColor(kBlue);
   gResolutionRatio->SetMarkerStyle(21);
   gResolutionRatio->SetMarkerSize(1);
   gResolutionRatio->Draw("PA");
   
   TGraph *gResolutionSubRatio = new TGraph(MC2Data, arrayMC, arrayEMCSubRatio);
   gResolutionSubRatio->SetMarkerColor(kRed);
   gResolutionSubRatio->SetMarkerStyle(21);
   gResolutionSubRatio->SetMarkerSize(1);
//   gResolutionSubRatio->Draw("P");
   
// TGraph *gResolutionStragglingRatio = new TGraph(MC2Data, arrayMC, arrayPSTARmaxDeltaRatio);
   TGraph *gResolutionStragglingRatio = new TGraph(nlines0, arrayMCActualResidualRange, arrayMCActualSigmaRatio);
   gResolutionStragglingRatio->SetLineWidth(2);
   gResolutionStragglingRatio->SetLineColor(kRed);
   gResolutionStragglingRatio->Draw("L");


   c6->cd(3);
   Float_t arrayStragglingsRatio[500] = {0};
   for (Int_t i=0; i<nlines0; i++) {
      arrayStragglingsRatio[i] = arrayEMCRatio[i] / arrayMCActualSigmaRatio[i];
   }

   TGraph *gStragglingsRatio = new TGraph(nlines0, arrayMCActualResidualRange, arrayStragglingsRatio);
   gStragglingsRatio->SetTitle(Form("Ratio resolution / straggling  using %d mm Al absorber;Depth in detector [WEPL mm]; Resolution / straggling", mmAbsorbator));
   gStragglingsRatio->SetLineWidth(2);
   gStragglingsRatio->SetLineColor(kRed);
   gStragglingsRatio->Draw("AL");

   gPad->Update();
   TLine *stragglingline = new TLine(0, 1, gPad->GetUxmax(), 1);
   stragglingline->Draw("same");


   TLegend *legResRatio = new TLegend(0.36, 0.76, 0.88, 0.88);
   legResRatio->SetTextSize(0.03);
   legResRatio->SetTextFont(22);
   legResRatio->AddEntry(gResolutionRatio, "Measured range spread", "P");
//   legResRatio->AddEntry(gResolutionSubRatio, "#Delta WEPL - Range straggling", "P");
   legResRatio->AddEntry(gResolutionStragglingRatio, "MC truth range straggling", "L");
   legResRatio->Draw();

   c2->cd();
   TGraph *gFraction = new TGraph(nlinesF, arrayFractionX, arrayFractionY);
   TGraph *gFraction2 = new TGraph(nlinesF, arrayFractionX, arrayFractionY2);
   TGraph *gFraction3 = new TGraph(nlinesF, arrayFractionX, arrayFractionY3);
   gFraction->GetXaxis()->SetRangeUser(5, 1200);
   gFraction->SetMaximum(100);
   gFraction->SetMinimum(0);
   gFraction->GetXaxis()->SetTitleSize(0.045);
   gFraction->GetYaxis()->SetTitleSize(0.045);
   gFraction->GetXaxis()->SetLabelSize(0.045);
   gFraction->GetYaxis()->SetLabelSize(0.045);
   gFraction->GetXaxis()->SetTitleFont(22);
   gFraction->GetYaxis()->SetTitleFont(22);
   gFraction->GetXaxis()->SetTitleOffset(0.9);
   gFraction->GetYaxis()->SetTitleOffset(1.2);
   gFraction->GetXaxis()->SetLabelFont(22);
   gFraction->GetYaxis()->SetLabelFont(22);
   gFraction->SetTitle("f;Number of protons in frame;Fraction of correctly reconstructed tracks");
   gFraction->SetLineColor(kGreen+2);
   gFraction->SetLineWidth(3);
   gFraction2->SetLineColor(kAzure+4);
   gFraction2->SetLineWidth(3);
   gFraction3->SetLineColor(kPink+4);
   gFraction3->SetLineWidth(3);

   gFraction->Draw("AL");
   gFraction2->Draw("L");
   gFraction3->Draw("L");
   
   TText *t = new TText();
   gFraction->GetYaxis()->SetLabelOffset(5);
   t->SetTextAlign(32);
   t->SetTextSize(0.04);
   t->SetTextFont(22);
   for (Int_t i=0; i<6;i++) {
      cout << "Drawing text at " << -0.42 << ", " << i*20 << endl;
      t->DrawText(-0.42, i*20, Form("%d%%", i*20));
   }


   gPad->SetLogx();
   gFraction->GetXaxis()->SetNoExponent();
   
   TLegend * leg2 = new TLegend(0.67, 0.76, 0.97, 0.93);
   leg2->SetTextSize(0.04);
   leg2->SetTextFont(22);
   leg2->AddEntry(gFraction, "Whole track correct", "L");
   leg2->AddEntry(gFraction2, "Correct endpoints (ID_{first} = ID_{last})", "L");
   leg2->AddEntry(gFraction3, "Close endpoints (#pm 0.5 mm, #pm 0.5#circ)", "L");
   leg2->Draw();

   c2->Update();

   printf("BEFORE\n");

   //c22->Divide(2,1,1e-5,1e-5);

   c22->cd();

//   c22->cd(2);
   TGraph *gFraction2mm = new TGraph(nlinesF2, arrayFraction2mmX, arrayFraction2mmY);
   TGraph *gFraction3mm = new TGraph(nlinesF3, arrayFraction3mmX, arrayFraction3mmY);
   TGraph *gFraction2mmFilter = new TGraph(nlinesF2, arrayFraction2mmX, arrayFraction2mmY2);
   TGraph *gFraction3mmFilter = new TGraph(nlinesF3, arrayFraction3mmX, arrayFraction3mmY2);
//   TGraph *gFraction35mm = new TGraph(nlinesF35, arrayFraction35mmX, arrayFraction35mmY);
   TGraph *gFraction4mm = new TGraph(nlinesF4, arrayFraction4mmX, arrayFraction4mmY);
   TGraph *gFraction5mm = new TGraph(nlinesF5, arrayFraction5mmX, arrayFraction5mmY);
   TGraph *gFraction6mm = new TGraph(nlinesF6, arrayFraction6mmX, arrayFraction6mmY);


   printf("Some values in arrayFraction2mmXY: (%.2f, %.2f), (%.2f, %.2f)\n", arrayFraction2mmX[0], arrayFraction2mmY[0], arrayFraction3mmY[1], arrayFraction4mmY[1]);
   printf("nlinesF2,F3,F4 = %d, %d, %d.\n", nlinesF2, nlinesF3, nlinesF4);
//   c22->SetLogx();
   gPad->SetGridy();
   gPad->SetGridx();
   gFraction2mmFilter->SetTitle(";Protons / readout frame;Fraction of correctly reconstructed tracks");
   gFraction2mmFilter->GetXaxis()->SetRangeUser(1, 520);
   gFraction2mmFilter->SetMaximum(100);
   gFraction2mmFilter->SetMinimum(40);
   gFraction2mmFilter->GetXaxis()->SetTitleSize(0.04);
   gFraction2mmFilter->GetYaxis()->SetTitleSize(0.04);
   gFraction2mmFilter->GetXaxis()->SetLabelSize(0.04);
   gFraction2mmFilter->GetXaxis()->SetNdivisions(9);
   gFraction2mmFilter->GetYaxis()->SetLabelSize(0.04);
   gFraction2mmFilter->GetXaxis()->SetLabelSize(0.04);
   gFraction2mmFilter->GetXaxis()->SetTitleFont(22);
   gFraction2mmFilter->GetYaxis()->SetTitleFont(22);
   gFraction2mmFilter->GetXaxis()->SetTitleOffset(1);
   gFraction2mmFilter->GetYaxis()->SetTitleOffset(1.3);
   gFraction2mmFilter->GetYaxis()->SetLabelOffset(5);
   gFraction2mmFilter->GetXaxis()->SetLabelFont(22);
   gFraction2mmFilter->GetXaxis()->SetNoExponent();
   gFraction2mmFilter->GetXaxis()->SetMoreLogLabels();
   gFraction2mmFilter->GetYaxis()->SetLabelFont(22);

   gFraction2mm->SetLineWidth(3);
   gFraction3mm->SetLineWidth(3);
   gFraction2mmFilter->SetLineWidth(3);
   gFraction3mmFilter->SetLineWidth(3);
   gFraction2mm->SetMarkerStyle(21);
   gFraction3mm->SetMarkerStyle(21);
   gFraction2mmFilter->SetMarkerStyle(21);
   gFraction3mmFilter->SetMarkerStyle(21);

   gFraction2mm->SetLineColor(kRed+1);
   gFraction3mm->SetLineColor(kRed-4);
   gFraction2mmFilter->SetLineColor(kRed+2);
   gFraction3mmFilter->SetLineColor(kRed);

//   gFraction2mm->Draw("LPA");
//   gFraction3mm->Draw("LP");
   gFraction2mmFilter->Draw("LPA");
   gFraction3mmFilter->Draw("LP");
//   gFraction35mm->Draw("L");
   gFraction4mm->Draw("LP");
   gFraction5mm->Draw("LP");
//   gFraction6mm->Draw("LP");
   
   TText *tt = new TText();   
   tt->SetTextAlign(32);
   tt->SetTextSize(0.04);
   tt->SetTextFont(22);
   for (Int_t i=4; i<11;i++) {
      cout << "Drawing text at " << -0.42 << ", " << i*10 << endl;
      tt->DrawText(-0.0052, i*10, Form("%d%%", i*10));
   }
   
   gPad->Update();
   TGaxis *axis = new TGaxis(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax(), 0.1, 52, 510, "-L");
   axis->Draw("same");
   axis->SetTitleFont(22);
   axis->SetLabelFont(22);
   axis->SetTitleSize(0.04);
   axis->SetLabelSize(0.04);
//   axis->SetTitleOffset(0.95);
   axis->SetTitle("Million primaries/s");
   
   TText *curveLabel = new TText();
   curveLabel->SetTextFont(22);
   curveLabel->SetTextSize(0.04);
//   curveLabel->DrawText(535, gFraction2mm->Eval(535)*0.98, "Ideal");
   curveLabel->DrawText(220, gFraction5mm->Eval(220)*0.98, "Helium 16 cm");
//   curveLabel->DrawText(220, gFraction35mm->Eval(220)*0.98, "5x5 mm^{2}");
   curveLabel->DrawText(220, gFraction4mm->Eval(220)*0.98, "Proton 16 cm");
   curveLabel->DrawText(220, gFraction2mmFilter->Eval(220)*0.98, "Helium 5 cm");
   curveLabel->DrawText(220, gFraction3mmFilter->Eval(220)*0.98, "Proton 5 cm");
//   curveLabel->DrawText(220, 91, "Pixel Diffusion Model");

   /*
   TLegend *legEfficiency = new TLegend(0.7, 0.6, 0.9, 0.9);
   legEfficiency->AddEntry(gFraction2mm, "Ideal", "L");
   legEfficiency->AddEntry(gFraction3mm, "Charge Diffusion", "L");
   legEfficiency->AddEntry(gFraction2mmFilter, "Ideal + 2nd filter", "L");
   legEfficiency->AddEntry(gFraction3mmFilter, "Charge Diffusion + 2nd filters", "L");
   legEfficiency->SetTextFont(22);
   legEfficiency->Draw();
*/
//   c22->cd(1);
   TGraph *gFraction2mmNoDiffusion = new TGraph(nlinesNoDiffusionF2, arrayFraction2mmNoDiffusionX, arrayFraction2mmNoDiffusionY);
   TGraph *gFraction3mmNoDiffusion = new TGraph(nlinesNoDiffusionF3, arrayFraction3mmNoDiffusionX, arrayFraction3mmNoDiffusionY);
   TGraph *gFraction4mmNoDiffusion = new TGraph(nlinesNoDiffusionF4, arrayFraction4mmNoDiffusionX, arrayFraction4mmNoDiffusionY);
   TGraph *gFraction5mmNoDiffusion = new TGraph(nlinesNoDiffusionF5, arrayFraction5mmNoDiffusionX, arrayFraction5mmNoDiffusionY);
   TGraph *gFraction6mmNoDiffusion = new TGraph(nlinesNoDiffusionF6, arrayFraction6mmNoDiffusionX, arrayFraction6mmNoDiffusionY);

   gPad->SetGridy();
   gPad->SetGridx();
   gFraction2mmNoDiffusion->SetTitle(";Protons / readout frame;Fraction of correctly reconstructed tracks");
   gFraction2mmNoDiffusion->GetXaxis()->SetRangeUser(1, 520);
   gFraction2mmNoDiffusion->SetMaximum(100);
   gFraction2mmNoDiffusion->SetMinimum(20);
   gFraction2mmNoDiffusion->GetXaxis()->SetTitleSize(0.05);
   gFraction2mmNoDiffusion->GetYaxis()->SetTitleSize(0.05);
   gFraction2mmNoDiffusion->GetXaxis()->SetLabelSize(0.05);
   gFraction2mmNoDiffusion->GetXaxis()->SetNdivisions(9);
   gFraction2mmNoDiffusion->GetYaxis()->SetLabelSize(0.05);
   gFraction2mmNoDiffusion->GetXaxis()->SetLabelSize(0.05);
   gFraction2mmNoDiffusion->GetXaxis()->SetTitleFont(22);
   gFraction2mmNoDiffusion->GetYaxis()->SetTitleFont(22);
   gFraction2mmNoDiffusion->GetXaxis()->SetTitleOffset(1);
   gFraction2mmNoDiffusion->GetYaxis()->SetTitleOffset(1.6);
   gFraction2mmNoDiffusion->GetYaxis()->SetLabelOffset(5);
   gFraction2mmNoDiffusion->GetXaxis()->SetLabelFont(22);
   gFraction2mmNoDiffusion->GetXaxis()->SetNoExponent();
   gFraction2mmNoDiffusion->GetXaxis()->SetMoreLogLabels();
   gFraction2mmNoDiffusion->GetYaxis()->SetLabelFont(22);
   gFraction2mmNoDiffusion->SetLineWidth(3);
   gFraction3mmNoDiffusion->SetLineWidth(3);
   gFraction4mmNoDiffusion->SetLineWidth(3);
   gFraction5mmNoDiffusion->SetLineWidth(3);
   gFraction6mmNoDiffusion->SetLineWidth(3);
   gFraction2mmNoDiffusion->SetMarkerStyle(21);
   gFraction3mmNoDiffusion->SetMarkerStyle(21);
   gFraction4mmNoDiffusion->SetMarkerStyle(21);
   gFraction5mmNoDiffusion->SetMarkerStyle(21);
   gFraction6mmNoDiffusion->SetMarkerStyle(21);

   gFraction2mmNoDiffusion->SetLineColor(kRed-4);
   gFraction3mmNoDiffusion->SetLineColor(kRed);
   gFraction4mmNoDiffusion->SetLineColor(kRed+1);
   gFraction5mmNoDiffusion->SetLineColor(kRed+2);
   gFraction6mmNoDiffusion->SetLineColor(kRed+3);

   /*
   gFraction2mmNoDiffusion->Draw("LPA");
   gFraction3mmNoDiffusion->Draw("LP");
   gFraction4mmNoDiffusion->Draw("LP");
   gFraction5mmNoDiffusion->Draw("LP");
   gFraction6mmNoDiffusion->Draw("LP");
   */
   
  /* 
   TText *tt2 = new TText();   
   tt2->SetTextAlign(32);
   tt2->SetTextSize(0.045);
   tt2->SetTextFont(22);
   for (Int_t i=2; i<11;i++) {
      cout << "Drawing text at " << -0.42 << ", " << i*10 << endl;
      tt2->DrawText(-0.0052, i*10, Form("%d%%", i*10));
   }
   gPad->Update();
   TGaxis *axis2 = new TGaxis(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax(), 0.1, 52, 510, "-L");
   axis2->Draw("same");
   axis2->SetTitleFont(22);
   axis2->SetLabelFont(22);
   axis2->SetTitleSize(0.05);
   axis2->SetLabelSize(0.05);
//   axis->SetTitleOffset(0.95);
   axis2->SetTitle("Million protons/s");
   
   TText *curveLabel2 = new TText();
   curveLabel2->SetTextFont(22);
   curveLabel2->SetTextSize(0.045);
   curveLabel2->DrawText(535, gFraction2mmNoDiffusion->Eval(535)*0.98, "2 mm");
   curveLabel2->DrawText(535, gFraction3mmNoDiffusion->Eval(535)*0.98, "3 mm");
   curveLabel2->DrawText(535, gFraction4mmNoDiffusion->Eval(535)*0.98, "4 mm");
   curveLabel2->DrawText(535, gFraction5mmNoDiffusion->Eval(535)*0.98, "5 mm");
   curveLabel2->DrawText(535, gFraction6mmNoDiffusion->Eval(535)*0.98, "6 mm");
   curveLabel2->DrawText(260, 91, "No Diffusion Model");

   TLegend *leg22 = new TLegend(0.21, 0.20, 0.50, 0.51);
   leg22->SetTextSize(0.04);
   leg22->SetTextFont(22);
   leg22->AddEntry(gFraction2mm, "2 mm Al absorber", "L");
   leg22->AddEntry(gFraction3mm, "3 mm Al absorber", "L");
   leg22->AddEntry(gFraction4mm, "4 mm Al absorber", "L");
   leg22->AddEntry(gFraction5mm, "5 mm Al absorber", "L");
   leg22->AddEntry(gFraction6mm, "6 mm Al absorber", "L");
   leg22->Draw();
   */

   c222->cd();
   TGraph *gFractionUniform = new TGraph(nlinesUniform, arrayFractionUniformX, arrayFractionUniformY);
   printf("Some values in arrayFractionUniformXY: (%.2f, %.2f), (%.2f, %.2f)\n", arrayFractionUniformX[0], arrayFractionUniformY[0], arrayFraction3mmY[1], arrayFraction4mmY[1]);
   printf("nlinesUniform = %d.\n", nlinesUniform);
   c222->SetGridy();
   c222->SetGridx();
   gPad->SetLogx();
   gFractionUniform->SetTitle(";Protons/100 cm^{2}/frame;Fraction of correctly reconstructed tracks");
   gFractionUniform->GetXaxis()->SetRangeUser(0, 15000);
   gFractionUniform->SetMaximum(100);
   gFractionUniform->SetMinimum(20);
   gFractionUniform->GetXaxis()->SetTitleSize(0.045);
   gFractionUniform->GetYaxis()->SetTitleSize(0.045);
   gFractionUniform->GetXaxis()->SetLabelSize(0.045);
   gFractionUniform->GetXaxis()->SetNdivisions(16);
   gFractionUniform->GetYaxis()->SetLabelSize(0.045);
   gFractionUniform->GetXaxis()->SetLabelSize(0.045);
   gFractionUniform->GetXaxis()->SetTitleFont(22);
   gFractionUniform->GetYaxis()->SetTitleFont(22);
   gFractionUniform->GetXaxis()->SetTitleOffset(0.9);
   gFractionUniform->GetYaxis()->SetTitleOffset(1.2);
   gFractionUniform->GetYaxis()->SetLabelOffset(5);
   gFractionUniform->GetXaxis()->SetLabelFont(22);
   gFractionUniform->GetXaxis()->SetNoExponent();
   gFractionUniform->GetXaxis()->SetMoreLogLabels();
   gFractionUniform->GetYaxis()->SetLabelFont(22);
   gFractionUniform->SetLineWidth(3);
   gFractionUniform->SetLineColor(3);
   gFractionUniform->Draw("LA");
   
   TText *ttt = new TText();   
   ttt->SetTextAlign(32);
   ttt->SetTextSize(0.045);
   ttt->SetTextFont(22);
   for (Int_t i=0; i<11;i++) {
      cout << "Drawing text at " << -0.42 << ", " << i*10 << endl;
      ttt->DrawText(-0.0052, i*10, Form("%d%%", i*10));
   }


   c3->cd();
   TGraph *gEfficiency = new TGraph(nEnergies, arrayEfficiencyEnergy, arrayEfficiencyFinal);
   TGraph *gEfficiency2 = new TGraph(nEnergies, arrayEfficiencyEnergy, arrayEfficiencyFinal2);
   gEfficiency->SetTitle("Efficiency of tracking algorithm at n_{p} = 500;Energy [MeV];Tracks reconstructed / n_{p}");
   gEfficiency->SetMinimum(0.7);
   gEfficiency->SetMaximum(1);
   gEfficiency->GetXaxis()->SetTitleFont(22);
   gEfficiency->GetYaxis()->SetTitleFont(22);
   gEfficiency->GetXaxis()->SetTitleOffset(1.2);
   gEfficiency->GetYaxis()->SetTitleOffset(1.2);
   gEfficiency->GetXaxis()->SetLabelFont(22);
   gEfficiency->GetYaxis()->SetLabelFont(22);
   gEfficiency->SetLineColor(kGreen+3);
   gEfficiency->SetLineWidth(3);
   gEfficiency2->SetLineColor(kGreen-3);
   gEfficiency2->SetLineWidth(3);
   gEfficiency->Draw("AL");
   gEfficiency2->Draw("L");

   c4->Divide(1, 2, 0.01, 0.001);
   c4->cd(1);

   TGraph *gAlignmentXMine = new TGraph(nMine, alignmentChipsMine, alignmentChipXMine);
   TGraph *gAlignmentYMine = new TGraph(nMine, alignmentChipsMine, alignmentChipYMine);
   TGraph *gAlignmentXOrig = new TGraph(96, alignmentChipsOrig, alignmentChipXOrig);
   TGraph *gAlignmentYOrig = new TGraph(96, alignmentChipsOrig, alignmentChipYOrig);

   gAlignmentXMine->SetMarkerStyle(21);
   gAlignmentYMine->SetMarkerStyle(21);
   gAlignmentXOrig->SetMarkerStyle(22);
   gAlignmentYOrig->SetMarkerStyle(22);
   gAlignmentXMine->SetMarkerColor(kRed);
   gAlignmentYMine->SetMarkerColor(kRed);
   gAlignmentXOrig->SetMarkerColor(kBlue);
   gAlignmentYOrig->SetMarkerColor(kBlue);
   
   gAlignmentXMine->SetTitle("Alignment correction for all chips");
   gAlignmentXMine->GetXaxis()->SetTitle("Chip number");
   gAlignmentXMine->GetXaxis()->SetLabelFont(22);
   gAlignmentXMine->GetXaxis()->SetTitleFont(22);
   gAlignmentXMine->GetYaxis()->SetTitleFont(22);
   gAlignmentXMine->GetYaxis()->SetLabelFont(22);
   gAlignmentXMine->GetYaxis()->SetTitle("Correction value in X direction [#mum]");
   gAlignmentXMine->GetXaxis()->SetNdivisions(54);
   gAlignmentXMine->Draw("AP");
   gAlignmentXOrig->Draw("P");
   
   gAlignmentXMine->GetYaxis()->SetRangeUser(-1000, 1000);
   gAlignmentXMine->GetXaxis()->SetRangeUser(0, 27.5);

   Float_t x_value;
   for (Int_t i=1; i<7; i++) {
      x_value = i*4 - 0.5;
      TLine *l = new TLine(x_value, -1000, x_value, 1000);
      l->Draw();
   }

   TLine *vl = new TLine(0, 0, 27.5, 0);
   vl->SetLineStyle(7);
   vl->SetLineWidth(2);
   vl->Draw();
   
   TLegend * leg3 = new TLegend(0.77, 0.71, 0.985, 0.94);
   leg3->SetTextSize(0.035);
   leg3->SetTextFont(22);
   leg3->AddEntry(gAlignmentXOrig, "Original correction values", "P");
   leg3->AddEntry(gAlignmentXMine, "My correction values", "P");
   leg3->Draw();

   c4->Update();

   c4->cd(2);
   gAlignmentYMine->SetTitle();
   gAlignmentYMine->GetXaxis()->SetTitle("Chip number");
   gAlignmentYMine->GetXaxis()->SetTitleFont(22);
   gAlignmentYMine->GetXaxis()->SetLabelFont(22);
   gAlignmentYMine->GetYaxis()->SetTitle("Correction value in Y direction [#mum]");
   gAlignmentYMine->GetYaxis()->SetTitleFont(22);
   gAlignmentYMine->GetYaxis()->SetLabelFont(22);
   gAlignmentYMine->GetXaxis()->SetNdivisions(54);
   gAlignmentYMine->Draw("AP");
   gAlignmentYOrig->Draw("P");
   gAlignmentYMine->GetYaxis()->SetRangeUser(-1000, 1000);
   gAlignmentYMine->GetXaxis()->SetRangeUser(0, 27.5);
   
   for (Int_t i=1; i<7; i++) {
      x_value = i*4 - 0.5;
      TLine *l = new TLine(x_value, -1000, x_value, 1000);
      l->Draw();
   }
   vl = new TLine(0, 0, 27.5, 0);
   vl->SetLineStyle(7);
   vl->SetLineWidth(2);
   vl->Draw();

   c5->cd();
   cal150->SetTitle("Sensitivity calibration factors for the different datasets;Chip number;Calibration factor");
   cal150->SetMinimum(0);
   cal150->SetMaximum(3);
   cal150->GetXaxis()->SetRangeUser(-1, 27);
   cal150->GetYaxis()->SetRangeUser(0.5, 2.7);
   cal150->GetXaxis()->SetNdivisions(54);
   cal150->GetXaxis()->SetTitleFont(22);   
   cal150->GetXaxis()->SetLabelFont(22);
   cal150->GetYaxis()->SetTitleFont(22);
   cal150->GetYaxis()->SetLabelFont(22);
   cal150->GetXaxis()->SetNdivisions(54);
   cal150->SetMarkerStyle(21);
   cal150->SetMarkerSize(0.7);
   cal150->SetMarkerColor(kBlue);
   cal160->SetMarkerStyle(21);
   cal160->SetMarkerSize(0.7);
   cal160->SetMarkerColor(kRed);
   cal170->SetMarkerStyle(21);
   cal170->SetMarkerSize(0.7);
   cal170->SetMarkerColor(kGreen);
   cal180->SetMarkerStyle(21);
   cal180->SetMarkerSize(0.7);
   cal180->SetMarkerColor(kOrange);
   cal188->SetMarkerStyle(21);
   cal188->SetMarkerSize(0.7);
   cal188->SetMarkerColor(kYellow+3);
   cal150->Draw("AP");
   cal160->Draw("P");
   cal170->Draw("P");
   cal180->Draw("P");
   cal188->Draw("P");

   for (Int_t i=-1; i<26; i++) {
      TLine *l = new TLine(i+0.5, 0.5, i+0.5, 2.7);
      if ((i+1)%4 != 0) l->SetLineColor(kGray);
      l->Draw();
   }

   for (Int_t i=0; i<28; i++) {
      Int_t    thisNumber = chipCalibrationNumber150[i] + chipCalibrationNumber160[i] + 
                            chipCalibrationNumber170[i] + chipCalibrationNumber180[i] + chipCalibrationNumber188[i];
      Float_t  thisCalibration = chipCalibrationFactor150[i] * chipCalibrationNumber150[i] + 
                                 chipCalibrationFactor160[i] * chipCalibrationNumber160[i] + 
                                 chipCalibrationFactor170[i] * chipCalibrationNumber170[i] + 
                                 chipCalibrationFactor180[i] * chipCalibrationNumber180[i] + 
                                 chipCalibrationFactor188[i] * chipCalibrationNumber188[i];
      if (thisNumber)
         thisCalibration /= thisNumber;
      cout << Form("Calibration factor for chip %d is %.3f, based on %d hits.\n", i, thisCalibration, thisNumber);
   }
   
   TLegend * leg4 = new TLegend(0.15, 0.65, 0.28, 0.88);
   leg4->SetTextSize(0.035);
   leg4->SetTextFont(22);
   leg4->AddEntry(cal150, "150 MeV", "P");
   leg4->AddEntry(cal160, "160 MeV", "P");
   leg4->AddEntry(cal170, "170 MeV", "P");
   leg4->AddEntry(cal180, "180 MeV", "P");
   leg4->AddEntry(cal188, "188 MeV", "P");
   leg4->Draw();

   c5->Update();

   c7->cd();
   Double_t dummy[21] = {};
   
   gPad->SetLogy();
   gPad->SetLogx();

   /*
   TGraphAsymmErrors *tgaBK     = new TGraphAsymmErrors(21, paraNPoints, paraBKm, dummy, dummy, paraBKl, paraBKh);
   TGraphAsymmErrors *tgaUlmer  = new TGraphAsymmErrors(21, paraNPoints, paraUlmerm, dummy, dummy, paraUlmerl, paraUlmerh);
   TGraphAsymmErrors *tgaSpline = new TGraphAsymmErrors(21, paraNPoints, paraSplinem, dummy, dummy, paraSplinel, paraSplineh);
   TGraphAsymmErrors *tgaLinear = new TGraphAsymmErrors(21, paraNPoints, paraLinearm, dummy, dummy, paraLinearl, paraLinearh);
   */

   TGraph *tgaBK     = new TGraph(pi, paraNPoints, paraBKh);
   TGraph *tgaUlmer  = new TGraph(pi, paraNPoints, paraUlmerh);
   TGraph *tgaSpline = new TGraph(pi, paraNPoints, paraSplineh);
   TGraph *tgaLinear = new TGraph(pi, paraNPoints, paraLinearh);

   tgaBK->SetTitle(";Number of data points in training group; 75% percentile of errors in range calculation [%]");

   tgaBK->SetLineColor(kColorBK);
   tgaBK->SetLineWidth(3);
   tgaBK->SetLineStyle(kLineBK);
   tgaUlmer->SetLineColor(kColorUlmer);
   tgaUlmer->SetLineStyle(kLineUlmer);
   tgaUlmer->SetLineWidth(3);
   tgaSpline->SetLineColor(kColorSpline);
   tgaSpline->SetLineStyle(kLineSpline);
   tgaSpline->SetLineWidth(3);
   tgaLinear->SetLineColor(kColorLinear);
   tgaLinear->SetLineStyle(kLineLinear);
   tgaLinear->SetLineWidth(3);

   tgaBK->GetXaxis()->SetTitleFont(22);
   tgaBK->GetXaxis()->SetLabelFont(22);
   tgaBK->GetYaxis()->SetTitleFont(22);
   tgaBK->GetYaxis()->SetLabelFont(22);

   tgaBK->Draw("AL");
   tgaUlmer->Draw("same, L");
   tgaSpline->Draw("same, L");
   tgaLinear->Draw("same, L");
   
   
   tgaBK->GetYaxis()->SetRangeUser(0.0005, 200);
   tgaBK->GetYaxis()->SetNoExponent();
   tgaBK->GetXaxis()->SetNoExponent();
   tgaBK->GetYaxis()->SetTitleOffset(1.2);

   TLegend *leg5 = new TLegend(0.72, 0.7, 0.94, 0.93);
   leg5->SetTextSize(0.035);
   leg5->SetTextFont(22);
   leg5->AddEntry(tgaBK, "Bragg-Kleeman", "L");
   leg5->AddEntry(tgaUlmer, "Sum of exponentials", "L");
   leg5->AddEntry(tgaLinear, "Linear interpolation", "L");
   leg5->AddEntry(tgaSpline, "Spline interpolation", "L");
   leg5->Draw();
   
   c1->SaveAs("OutputFiles/figures/finalPlotsForArticle/estimated_ranges_all_energies.eps");
   c1->SaveAs("OutputFiles/figures/finalPlotsForArticle/estimated_ranges_all_energies.root");
}
