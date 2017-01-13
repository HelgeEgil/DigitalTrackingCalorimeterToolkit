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

using namespace std;

void makePlots() {
   TCanvas *c1 = new TCanvas("c1", "Fit results", 1200, 800);
   TPad *pad1 = new TPad("pad1", "70", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad2 = new TPad("pad2", "30", 0.0, 0.0, 1.0, 0.3, 0);
   pad1->Draw();
   pad2->Draw();
   
   TCanvas *c2 = new TCanvas("c2", "Correct Tracks fraction", 1200, 800);
   TCanvas *c3 = new TCanvas("c3", "Reconstruction efficiency", 1200, 800);
   TCanvas *c4 = new TCanvas("c4", "Chip alignment", 1200, 800);
   TCanvas *c5 = new TCanvas("c5", "Chip sensitivity calibration", 1200, 800);
   TCanvas *c6 = new TCanvas("c6", "Resolution", 1200, 800);
   c6->Divide(2,1,0.0001,0.0001);
   TCanvas *c7 = new TCanvas("c7", "Parameterization accuracy", 1200, 800);

   Float_t  arrayE[200] = {0}; // energy MC
   Float_t  arrayEE[200] = {0}; // error on energy MC
   Float_t  arrayMC[200] = {0}; // range MC
   Float_t  arrayEMC[200] = {0}; // error on range MC
   Float_t  arrayEMCRatio[200] = {0}; // error on range MC
   Float_t  arrayEMCSub[200] = {0}; // error on range MC
   Float_t  arrayEMCSubRatio[200] = {0}; // error on range MC
   Float_t  arrayMCDelta[200] = {0}; // range MC
   Float_t  arrayPSTAR[200] = {0};
   Float_t  arrayPSTARDelta[200] = {0};
   Float_t  arrayPSTARshade[400] = {0};
   Float_t  arrayPSTARmin[200] = {0};
   Float_t  arrayPSTARmax[200] = {0};
   Float_t  arrayPSTARminDelta[200] = {0};
   Float_t  arrayPSTARmaxDelta[200] = {0};
   Float_t  arrayPSTARmaxDeltaRatio[200] = {0};
   Float_t  arrayEPstar[200] = {0};
   Float_t  arrayE2[200] = {0}; // energy data 
   Float_t  arrayEE2[200] = {0}; // error on energy data
   Float_t  arrayEData[200] = {0};  // range data
   Float_t  arrayData[200] = {0}; // error on range data
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
   Double_t paraNPoints[21] = {};
   Double_t paraBKl[21] = {};
   Double_t paraBKm[21] = {};
   Double_t paraBKh[21] = {};
   Double_t paraBKInvl[21] = {};
   Double_t paraBKInvm[21] = {};
   Double_t paraBKInvh[21] = {};
   Double_t paraUlmerl[21] = {};
   Double_t paraUlmerm[21] = {};
   Double_t paraUlmerh[21] = {};
   Double_t paraUlmerInvl[21] = {};
   Double_t paraUlmerInvm[21] = {};
   Double_t paraUlmerInvh[21] = {};
   Double_t paraSplinel[21] = {};
   Double_t paraSplinem[21] = {};
   Double_t paraSplineh[21] = {};
   Double_t paraSplineInvl[21] = {};
   Double_t paraSplineInvm[21] = {};
   Double_t paraSplineInvh[21] = {};
   Double_t paraLinearl[21] = {};
   Double_t paraLinearm[21] = {};
   Double_t paraLinearh[21] = {};
   Double_t paraLinearInvl[21] = {};
   Double_t paraLinearInvm[21] = {};
   Double_t paraLinearInvh[21] = {};

   Int_t nThisEnergy = 0, lastEnergy = 0, mmAbsorbator;

   gStyle->SetOptStat(0);

   ifstream in;
   in.open("OutputFiles/result_makebraggpeakfit.csv");

   cout << "Opened file.\n";

   Float_t nomrange_, estrange_, sigmaRange_, lastRange_;
   Int_t energy_, thickness_;
   Float_t estimatedStraggling;

   Int_t nlines = 0;
   TNtuple *ntuple = new TNtuple("ntuple", "data from file", "energy_:nomrange_:estrange_:sigmaRange:lastRange_");
   
   Int_t MC2Data = -1;

   Float_t meanError = 0;
   Float_t meanAbsError = 0;
   Float_t meanSigma = 0;
   Float_t a_dtc = 0;
   Float_t p_dtc = 0;
   Float_t aprime_dtc = 0;
   Float_t a_wtr = 0.02387;
   Float_t p_wtr = 1.7547;
   Float_t aprime_wtr = 0.0087; // MeV^2 / mm

   while (1) {
      in >> thickness_ >> energy_ >> nomrange_ >> estrange_ >> sigmaRange_;

      if (!in.good()) {
         break;
      }

      mmAbsorbator = thickness_;
      meanError += ( estrange_ - nomrange_ ) / nomrange_;
      meanAbsError += fabs(( estrange_ - nomrange_ ) / nomrange_);
      meanSigma += sigmaRange_;
      if       (mmAbsorbator == 5) estimatedStraggling = 1.57e-2 * nomrange_ + 7.54e-6 * pow(nomrange_,2);
      else if  (mmAbsorbator == 4) estimatedStraggling = 1.73e-2 * nomrange_ + 2.28e-6 * pow(nomrange_,2);
      else if  (mmAbsorbator == 3) estimatedStraggling = 1.52e-2 * nomrange_ + 1.02e-5 * pow(nomrange_,2);
      else if  (mmAbsorbator == 2) estimatedStraggling = 1.53e-2 * nomrange_ + 9.59e-6 * pow(nomrange_,2);

//      estimatedStraggling = 0.012 * pow(nomrange_, 0.935); // using water values for minimized straggling limit

      if (mmAbsorbator == 2) {
         a_dtc = 0.0096;
         p_dtc = 1.784;
         aprime_dtc = 0.01938;
      }
      
      else if (mmAbsorbator == 3) {
         a_dtc = 0.0173389;
         p_dtc = 1.7751;
         aprime_dtc = 0.01971;
      }

      else if (mmAbsorbator == 4) {
         a_dtc = 0.0117;
         p_dtc = 1.7450;
         aprime_dtc = 0.01988;
      }

      Float_t zprime = energy_;
      Float_t wtr_range = a_wtr * pow(250, p_wtr);
      Float_t dtc_term = aprime_dtc * (p_dtc * pow(a_dtc, 2/p_dtc)) / (3 - 2/p_dtc);
      Float_t wtr_term = aprime_wtr * (p_wtr * pow(a_wtr, 2/p_wtr)) / (3 - 2/p_wtr);
      Float_t wepl_ratio = pow(a_wtr / a_dtc * pow(nomrange_ / a_wtr, 1 - p_dtc / p_wtr), 2);
      Float_t dtc_strag = wepl_ratio * dtc_term * pow(nomrange_/wepl_ratio, 3-2/p_dtc);
      Float_t wtr_strag = wtr_term * (pow(wtr_range, 3-2/p_wtr) - pow(wtr_range - zprime, 3-2/p_wtr));

      printf("Using a %.2f mm water phantom, the estimated straggling from water is %.2f mm and from DTC is %.2f mm. WEPL ratio is %.2f.\n", zprime, sqrt(wtr_strag), sqrt(dtc_strag), sqrt(wepl_ratio));

      estimatedStraggling = sqrt(dtc_strag + wtr_strag);   

      if (nlines < MC2Data || MC2Data<0) {
         cout << "Line " << nlines << ", energy " << energy_ << ",  MC" << endl;
         arrayE[nlines] = energy_;
         arrayEE[nlines] = 0;
         arrayMC[nlines] = estrange_;
         arrayMCDelta[nlines] = estrange_ - nomrange_;
         arrayEMC[nlines] = sigmaRange_;
         arrayEMCRatio[nlines] = sigmaRange_ / nomrange_ * 100;
         arrayEMCSub[nlines] = sqrt(abs(pow(estimatedStraggling, 2) - pow(sigmaRange_, 2)));
         arrayEMCSubRatio[nlines] = sqrt(abs(pow(estimatedStraggling, 2) - pow(sigmaRange_, 2))) / nomrange_ * 100;
         arrayPSTAR[nlines] = nomrange_;
         arrayPSTARmin[nlines] = nomrange_ - estimatedStraggling;
         arrayPSTARmax[nlines] = nomrange_ + estimatedStraggling;
         arrayPSTARDelta[nlines] = 0;
         arrayPSTARminDelta[nlines] = -estimatedStraggling;
         arrayPSTARmaxDelta[nlines] = estimatedStraggling;
         arrayPSTARmaxDeltaRatio[nlines] = estimatedStraggling / nomrange_ * 100;
         arrayEPstar[nlines] = 0;
      }

      else {
         cout << "Line " << nlines << ", (or " << nlines-MC2Data << "), energy " << energy_ << ", data" << endl;
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
   in2.open("OutputFiles/lastLayerCorrect_different_nRuns.csv");
   Float_t factor, np, correctLast, correctWhole, lastIsFirst, lastIsAlmostFirst;
   Float_t arrayFractionX[200] = {0};
   Float_t arrayFractionY[200] = {0};
   Float_t arrayFractionY2[200] = {0};
   Float_t arrayFractionY3[200] = {0};
   Int_t nlines2 = 0;
   while (1) {
      in2 >> factor >> np >> correctWhole >> lastIsFirst >> lastIsAlmostFirst;

      if (!in2.good()) break;
      
      arrayFractionX[nlines2] = np;
      arrayFractionY[nlines2] = correctWhole * 100;
      arrayFractionY2[nlines2] = lastIsFirst * 100;
      arrayFractionY3[nlines2] = lastIsAlmostFirst * 100;

      nlines2++;
   }
   cout << "Found " << nlines2 << " lines in lastLayerCorrect.\n";
   
   in2.close();
   
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
   Int_t    pN;
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

   gPad->Update();
   TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();

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

   gPad->Update();
   TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   title->SetTextSize(0.08);
   gPad->Modified();
   pad2->Update();

   c6->cd(1);
   TGraph *gResolution = new TGraph(MC2Data, arrayMC, arrayEMC);
   gResolution->SetTitle(Form("WEPL resolution using %d mm Al absorber;WEPL [mm];#Delta WEPL [mm]", mmAbsorbator));
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
   
   TGraph *gResolutionStraggling = new TGraph(MC2Data, arrayMC, arrayPSTARmaxDelta);
   gResolutionStraggling->SetLineWidth(2);
   gResolutionStraggling->SetLineColor(kRed);
   gResolutionStraggling->Draw("L");

   TLegend *legRes = new TLegend(0.17, 0.77, 0.69, 0.88);
   legRes->SetTextSize(0.03);
   legRes->SetTextFont(22);
   legRes->AddEntry(gResolution, "#Delta WEPL", "P");
//   legRes->AddEntry(gResolutionSub, "#Delta WEPL - Range straggling", "P");
   legRes->AddEntry(gResolutionStraggling, "Range straggling", "L");
   legRes->Draw();

   c6->cd(2);
   TGraph *gResolutionRatio = new TGraph(MC2Data, arrayMC, arrayEMCRatio);
   gResolutionRatio->SetTitle(Form("WEPL resolution using %d mm Al absorber;WEPL [mm];#Delta WEPL / WEPL [%]", mmAbsorbator));
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
   
   TGraph *gResolutionStragglingRatio = new TGraph(MC2Data, arrayMC, arrayPSTARmaxDeltaRatio);
   gResolutionStragglingRatio->SetLineWidth(2);
   gResolutionStragglingRatio->SetLineColor(kRed);
   gResolutionStragglingRatio->Draw("L");

   TLegend *legResRatio = new TLegend(0.36, 0.76, 0.88, 0.88);
   legResRatio->SetTextSize(0.03);
   legResRatio->SetTextFont(22);
   legResRatio->AddEntry(gResolutionRatio, "#Delta WEPL", "P");
//   legResRatio->AddEntry(gResolutionSubRatio, "#Delta WEPL - Range straggling", "P");
   legResRatio->AddEntry(gResolutionStragglingRatio, "Range straggling", "L");
   legResRatio->Draw();

   c2->cd();
   TGraph *gFraction = new TGraph(nlines2, arrayFractionX, arrayFractionY);
   TGraph *gFraction2 = new TGraph(nlines2, arrayFractionX, arrayFractionY2);
   TGraph *gFraction3 = new TGraph(nlines2, arrayFractionX, arrayFractionY3);
   gFraction->GetXaxis()->SetRangeUser(10, 6000);
   gFraction->SetMaximum(100);
   gFraction->SetMinimum(0);
   gFraction->GetXaxis()->SetTitleSize(0.045);
   gFraction->GetYaxis()->SetTitleSize(0.045);
   gFraction->GetXaxis()->SetLabelSize(0.045);
   gFraction->GetYaxis()->SetLabelSize(0.045);
   gFraction->GetXaxis()->SetTitleFont(22);
   gFraction->GetYaxis()->SetTitleFont(22);
   gFraction->GetXaxis()->SetTitleOffset(0.9);
   gFraction->GetYaxis()->SetTitleOffset(1.1);
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

   gPad->Update();
   title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();
   
   TLegend * leg2 = new TLegend(0.67, 0.76, 0.97, 0.93);
   leg2->SetTextSize(0.04);
   leg2->SetTextFont(22);
   leg2->AddEntry(gFraction, "Whole track correct", "L");
   leg2->AddEntry(gFraction2, "Correct endpoints (ID_{first} = ID_{last})", "L");
   leg2->AddEntry(gFraction3, "Close endpoints (#pm 0.5 mm, #pm 0.5#circ)", "L");
   leg2->Draw();

   c2->Update();

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
   
   gPad->Update();
   title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();

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
   TLine *vl = new TLine(0, 0, 27.5, 0);
   vl->SetLineStyle(7);
   vl->SetLineWidth(2);
   vl->Draw();

   c5->cd();
   cal150->SetTitle("Sensitivity calibration factors for the different datasets;Chip number;Calibration factor");
   cal150->SetMinimum(0);
   cal150->SetMaximum(3);
   cal150->GetXaxis()->SetRangeUser(-1, 27);
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
      TLine *l = new TLine(i+0.5, 0, i+0.5, 3);
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
   
   gPad->Update();
   title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();
   
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
   TGraphAsymmErrors *tgaBK     = new TGraphAsymmErrors(21, paraNPoints, paraBKm, dummy, dummy, paraBKl, paraBKh);
   TGraphAsymmErrors *tgaUlmer  = new TGraphAsymmErrors(21, paraNPoints, paraUlmerm, dummy, dummy, paraUlmerl, paraUlmerh);
   TGraphAsymmErrors *tgaSpline = new TGraphAsymmErrors(21, paraNPoints, paraSplinem, dummy, dummy, paraSplinel, paraSplineh);
   TGraphAsymmErrors *tgaLinear = new TGraphAsymmErrors(21, paraNPoints, paraLinearm, dummy, dummy, paraLinearl, paraLinearh);

   tgaBK->SetTitle("Median error (+- 1st & 4rd quartile) over all energies of Bragg Curve parameterizations;Number of data points for model fit; Median error [%]");

   tgaBK->SetLineColor(kRed);
   tgaBK->SetMarkerColor(kRed);
   tgaUlmer->SetLineColor(kBlue);
   tgaUlmer->SetMarkerColor(kBlue);
   tgaSpline->SetLineColor(kBlack);
   tgaSpline->SetMarkerColor(kBlack);
   tgaLinear->SetLineColor(kGreen);
   tgaLinear->SetMarkerColor(kGreen);

   tgaBK->Draw("ALP");
   tgaUlmer->Draw("same, LP");
   tgaSpline->Draw("same, LP");
   tgaLinear->Draw("same, LP");
   
   c1->SaveAs("OutputFiles/figures/finalPlotsForArticle/estimated_ranges_all_energies.eps");
   c1->SaveAs("OutputFiles/figures/finalPlotsForArticle/estimated_ranges_all_energies.root");
   c2->SaveAs("OutputFiles/figures/finalPlotsForArticle/fraction_of_correct_tracks.eps");
   c2->SaveAs("OutputFiles/figures/finalPlotsForArticle/fraction_of_correct_tracks.root");
}
