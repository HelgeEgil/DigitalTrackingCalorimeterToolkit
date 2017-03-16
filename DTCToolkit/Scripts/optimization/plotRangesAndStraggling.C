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

Int_t absorberThickness = 2;
Bool_t kFilterData = true;
Int_t filterSize = 15;
const Int_t arraySize = 500;

void filterArray(Float_t *array, Int_t filterSize) {
   Float_t tempArray[arraySize];
   Float_t value;
   Int_t n;

   for (Int_t i=0; i<arraySize; i++) {
      value = 0;
      n = 0;

      for (Int_t j=i-filterSize/2; j<=i+filterSize/2; j++) {
         if (j<0 || j>=arraySize) {
            continue;
         }
         value += array[j];
         n++;
      }

      value /= n;
      tempArray[i] = value;
   }

   for (Int_t i=0; i<arraySize; i++) {
      array[i] = tempArray[i];
   }
}

void plotRangesAndStraggling() {
   TCanvas *c1 = new TCanvas("c1", "Fit results", 1200, 800);
   TPad *pad1 = new TPad("pad1", "70", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad2 = new TPad("pad2", "30", 0.0, 0.0, 1.0, 0.3, 0);
   pad1->Draw();
   pad2->Draw();
   
   TCanvas *c6 = new TCanvas("c6", "Resolution", 1800, 600);
   c6->Divide(4,1,0.0001,0.0001);

   Float_t  arrayE[arraySize] = {0}; // energy MC
   Float_t  arrayMCActualSigma[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayMCActualSigmaRatio[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayMCActualResidualRange[arraySize] = {0}; // Measured range straeggling from full MC
   Float_t  arrayStragglingsRatio[arraySize] = {0};
   Float_t  arrayEE[arraySize] = {0}; // error on energy MC
   Float_t  arrayMC[arraySize] = {0}; // range MC
   Float_t  arrayEMC[arraySize] = {0}; // error on range MC
   Float_t  arrayEMCRatio[arraySize] = {0}; // error on range MC
   Float_t  arrayEMCSub[arraySize] = {0}; // error on range MC
   Float_t  arrayEMCSubRatio[arraySize] = {0}; // error on range MC
   Float_t  arrayMCDelta[arraySize] = {0}; // range MC
   Float_t  arrayPSTAR[arraySize] = {0};
   Float_t  arrayPSTARDelta[arraySize] = {0};
   Float_t  arrayPSTARshade[400] = {0};
   Float_t  arrayPSTARmin[arraySize] = {0};
   Float_t  arrayPSTARmax[arraySize] = {0};
   Float_t  arrayPSTARminDelta[arraySize] = {0};
   Float_t  arrayPSTARmaxDelta[arraySize] = {0};
   Float_t  arrayPSTARmaxDeltaRatio[arraySize] = {0};
   Float_t  arrayEPstar[arraySize] = {0};
   Float_t  arrayE2[arraySize] = {0}; // energy data 
   Float_t  arrayEE2[arraySize] = {0}; // error on energy data
   Float_t  arrayEData[arraySize] = {0};  // range data
   Float_t  arrayData[arraySize] = {0}; // error on range data

   Int_t nThisEnergy = 0, lastEnergy = 0, mmAbsorbator;

   gStyle->SetOptStat(0);

   Float_t nomrange_, estrange_, sigmaRange_, lastRange_, nomsigma_, waterphantomthickness_, dummy0;
   Int_t energy_, thickness_;
   Float_t estimatedStraggling;

   ifstream in0;
   in0.open("../../OutputFiles/findManyRangesDegrader.csv");
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
   else if (absorberThickness == 4) {
      a_dtc = 0.018745;
      p_dtc = 1.660262;
   }

   Float_t wepl_ratio0 = a_wtr / a_dtc * pow(250 / a_wtr, 1 - p_dtc / p_wtr);
//   Float_t wtr_range = a_wtr * pow(250, p_wtr);
   Float_t wtr_range = 378.225; // GATE

   while (1) {
      in0 >>  waterphantomthickness_ >> thickness_ >> nomrange_ >> nomsigma_ >> dummy0 >> dummy0 >> dummy0;

      if (!in0.good()) {
         break;
      }

      if (thickness_ != absorberThickness) {
         continue;
      }

      arrayMCActualSigma[nlines0] = nomsigma_; //  * wepl_ratio0;
      arrayMCActualSigmaRatio[nlines0] = nomsigma_ * 100 * wepl_ratio0 / wtr_range;
      arrayMCActualResidualRange[nlines0++] = nomrange_; //  * wepl_ratio0;
   }
   in0.close();

   ifstream in;
   in.open("../../OutputFiles/result_makebraggpeakfit.csv");

   cout << "Opened file.\n";

   Int_t nlines = 0;
   
   Float_t meanError = 0;
   Float_t meanAbsError = 0;
   Float_t meanSigma = 0;

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
      
      estimatedStraggling = nomsigma_; 

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
      arrayEPstar[nlines++] = 0;

   }

   meanError /= nlines;
   meanSigma /= nlines;

   cout << "Mean error on fit range is " << meanError << " mm.\n";
   cout << "Mean | error | on fit range is " << meanAbsError << " mm.\n";
   cout << "Mean SIGMA on fit range is " << meanSigma << " mm.\n";

   in.close();

   pad1->cd();
   
   TGraphErrors *hMC = new TGraphErrors(nlines, arrayE, arrayMC, arrayEE, arrayEMC);
   TGraphErrors *hData = new TGraphErrors(nlines-nlines, arrayE2, arrayData, arrayEE2, arrayEData);
//   TGraphErrors *pstar = new TGraphErrors(nlines, arrayE, arrayPSTAR, arrayEE, arrayEPstar);
   TGraph *pstar = new TGraph(nlines, arrayE, arrayPSTAR);
   TGraph *pstarmin = new TGraph(nlines, arrayE, arrayPSTARmin);
   TGraph *pstarmax = new TGraph(nlines, arrayE, arrayPSTARmax);
   TGraph *pstarshade = new TGraph(nlines*2);
   
   for (Int_t i=0; i<nlines; i++) {
      pstarshade->SetPoint(i, arrayE[i], arrayPSTARmax[i]);
      pstarshade->SetPoint(nlines+i, arrayE[nlines-i-1], arrayPSTARmin[nlines-i-1]);
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
   pstarshade->SetTitle(Form("Reconstructed ranges of proton beams with %d mm Al absorbator;Degrader thickness [mm];Reconstructed WET range [mm]", mmAbsorbator));
   
   pstarshade->SetFillColor(kMagenta-10);
   pstarshade->Draw("FA");
   pstarmin->Draw("L");
   pstarmax->Draw("L");

   hMC->SetTitle("Reconstructed ranges #LT#hat{R_{0}}#GT of proton tracks;Degrader thickness [mm];Projected range [mm]");
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

   TGraphErrors *hMCD = new TGraphErrors(nlines, arrayE, arrayMCDelta, arrayEE, arrayEMC);
   TGraph *pstarD = new TGraph(nlines, arrayE, arrayPSTARDelta);
   TGraph *pstarminD = new TGraph(nlines, arrayE, arrayPSTARminDelta);
   TGraph *pstarmaxD = new TGraph(nlines, arrayE, arrayPSTARmaxDelta);
   TGraph *pstarshadeD = new TGraph(nlines*2);

   for (Int_t i=0; i<nlines; i++) {
      pstarshadeD->SetPoint(i, arrayE[i], arrayPSTARmaxDelta[i]);
      pstarshadeD->SetPoint(nlines+i, arrayE[nlines-i-1], arrayPSTARminDelta[nlines-i-1]);
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

   if (kFilterData) {
      filterArray(arrayEMC, filterSize);
      filterArray(arrayMCActualSigma, filterSize);
      filterArray(arrayEMCRatio, filterSize);
      filterArray(arrayEMCSubRatio, filterSize);
      filterArray(arrayMCActualSigmaRatio, filterSize);
      filterArray(arrayStragglingsRatio, filterSize);
   }

   c6->cd(1);
   TGraph *gResolution = new TGraph(nlines, arrayMC, arrayEMC);
   gResolution->SetTitle(Form("WEPL resolution using %d mm Al absorber;Depth in detector [WEPL mm];#Delta WEPL [mm]", mmAbsorbator));
   gResolution->GetYaxis()->SetRangeUser(2,5.5);
   gResolution->SetMarkerColor(kBlue);
   gResolution->SetMarkerStyle(21);
   gResolution->SetMarkerSize(1);
   gResolution->SetLineWidth(3);
   gResolution->SetLineColor(kBlue);
   gResolution->Draw("LA");
   
   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->InsertText(Form("Data filtered, size %d.", filterSize));
      title->SetTextSize(0.04);
      gPad->Modified();
   }
   
   TGraph *gResolutionStraggling = new TGraph(nlines0, arrayMCActualResidualRange, arrayMCActualSigma);
   gResolutionStraggling->SetLineWidth(3);
   gResolutionStraggling->SetLineColor(kRed);
   gResolutionStraggling->Draw("L");

   gPad->Update();
   // Janni: Straggling in water is 1.063 percent
   // Range in water is 379.4 mm (PSTAR) or 382.57 (Janni)
   // Straggling in water is 3.791 mm as measured in GATE
   Float_t waterStraggling = 3.791; // GATE
   TLine *lWaterStraggling = new TLine(0, waterStraggling, gPad->GetUxmax(), waterStraggling);
   lWaterStraggling->SetLineWidth(2);
   lWaterStraggling->SetLineColor(kGreen);
   lWaterStraggling->Draw();

   TLegend *legRes = new TLegend(0.17, 0.77, 0.69, 0.88);
   legRes->SetTextSize(0.03);
   legRes->SetTextFont(22);
   legRes->AddEntry(gResolution, "Measured range spread", "L");
   legRes->AddEntry(gResolutionStraggling, "MC truth range straggling", "L");
   legRes->AddEntry(lWaterStraggling, "Straggling in water", "L");
   legRes->Draw();

   c6->cd(2);
   TGraph *gResolutionRatio = new TGraph(nlines, arrayMC, arrayEMCRatio);
   gResolutionRatio->SetTitle(Form("WEPL resolution using %d mm Al absorber;Depth in detector [WEPL mm];#Delta WEPL / WEPL [%]", mmAbsorbator));
   gResolutionRatio->GetYaxis()->SetTitleOffset(1.5);
   gResolutionRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
   gResolutionRatio->SetMarkerColor(kBlue);
   gResolutionRatio->SetLineColor(kBlue);
   gResolutionRatio->SetLineWidth(3);
   gResolutionRatio->SetMarkerStyle(21);
   gResolutionRatio->SetMarkerSize(1);
   gResolutionRatio->Draw("LA");

   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->SetTextSize(0.04);
      title->InsertText(Form("Data filtered, size %d.", filterSize));
      gPad->Modified();
   }
   
   TGraph *gResolutionSubRatio = new TGraph(nlines, arrayMC, arrayEMCSubRatio);
   gResolutionSubRatio->SetMarkerColor(kRed);
   gResolutionSubRatio->SetMarkerStyle(21);
   gResolutionSubRatio->SetMarkerSize(1);
//   gResolutionSubRatio->Draw("P");
   
   gPad->Update();
   // Janni: Straggling in water is 1.063 percent
   // Range in water is 379.4 mm (PSTAR) or 382.57 (Janni)
   // Straggling in water is 3.791 mm as measured in GATE
   waterStraggling = 3.791 * 100 / 378.225; // GATE
   TLine *lWaterStraggling = new TLine(0, waterStraggling, gPad->GetUxmax(), waterStraggling);
   lWaterStraggling->SetLineWidth(2);
   lWaterStraggling->SetLineColor(kGreen);
   lWaterStraggling->Draw();
   
   TGraph *gResolutionStragglingRatio = new TGraph(nlines0, arrayMCActualResidualRange, arrayMCActualSigmaRatio);
   gResolutionStragglingRatio->SetLineWidth(3);
   gResolutionStragglingRatio->SetLineColor(kRed);
   gResolutionStragglingRatio->Draw("L");
   
   
   TLegend *legResRatio = new TLegend(0.36, 0.76, 0.88, 0.88);
   legResRatio->SetTextSize(0.03);
   legResRatio->SetTextFont(22);
   legResRatio->AddEntry(gResolutionRatio, "Measured range spread", "L");
   legResRatio->AddEntry(gResolutionStragglingRatio, "MC truth range straggling", "L");
   legResRatio->AddEntry(lWaterStraggling, "Straggling in water", "L");
   legResRatio->Draw();

   c6->cd(3);
   for (Int_t i=0; i<nlines0; i++) {
      arrayStragglingsRatio[i] = arrayEMCRatio[i] / arrayMCActualSigmaRatio[i];
   }

   TGraph *gStragglingsRatio = new TGraph(nlines0, arrayMCActualResidualRange, arrayStragglingsRatio);
   gStragglingsRatio->SetTitle(Form("Resolution/straggling, %d mm Al absorber;Depth in detector [WEPL mm]; Resolution / straggling", mmAbsorbator));
   gStragglingsRatio->GetYaxis()->SetRangeUser(0.3, 1.4);
   gStragglingsRatio->SetLineWidth(2);
   gStragglingsRatio->SetLineColor(kRed);
   gStragglingsRatio->Draw("AL");
   
   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->SetTextSize(0.04);
      title->InsertText(Form("Data filtered, size %d.", filterSize));
      gPad->Modified();
   }


   gPad->Update();
   TLine *stragglingline = new TLine(0, 1, gPad->GetUxmax(), 1);
   stragglingline->Draw("same");

   Float_t subtractMCStraggling[arraySize];
   Float_t subtractWaterStraggling[arraySize];
   Float_t waterStraggling = 379.4 * 0.01063;

   for (Int_t i=0; i<arraySize; i++) {
      if (arrayEMC[i] > arrayMCActualSigma[i]) {
         subtractMCStraggling[i] = sqrt(pow(arrayEMC[i], 2) - pow(arrayMCActualSigma[i], 2));
      }
      else {
         subtractMCStraggling[i] = 0;
      }

      if (arrayEMC[i] > waterStraggling) {
         subtractWaterStraggling[i] = sqrt(pow(arrayEMC[i], 2) - pow(waterStraggling, 2));
      }
      else {
         subtractWaterStraggling[i] = 0;
      }
   }

   c6->cd(4);
   TGraph *gSubtractFullResolution = new TGraph(nlines, arrayMC, arrayEMC);
   TGraph *gSubtractMCResolution = new TGraph(nlines, arrayMC, subtractMCStraggling);
   TGraph *gSubtractWaterResolution = new TGraph(nlines, arrayMC, subtractWaterStraggling);

   gSubtractFullResolution->SetTitle(Form("Resolution deconstruction %d mm Al absorber; Depth in detector [WEPL mm]; Deconstructed resolution [WEPL mm]", mmAbsorbator));

   gSubtractFullResolution->GetYaxis()->SetRangeUser(0, 6);
   gSubtractFullResolution->SetLineColor(kRed);
   gSubtractFullResolution->SetLineWidth(3);
   gSubtractMCResolution->SetLineColor(kBlue);
   gSubtractMCResolution->SetLineWidth(3);
   gSubtractWaterResolution->SetLineColor(kGreen);
   gSubtractWaterResolution->SetLineWidth(3);

   gSubtractFullResolution->Draw("LA");
   gSubtractMCResolution->Draw("same");
   gSubtractWaterResolution->Draw("same");

   TLegend *legSubtractResolution = new TLegend(0.50, 0.78, 0.96, 0.93);
   legSubtractResolution->SetTextSize(0.03);
   legSubtractResolution->SetTextFont(22);
   legSubtractResolution->AddEntry(gSubtractFullResolution, "Range resolution", "L");
   legSubtractResolution->AddEntry(gSubtractMCResolution, "#sqrt{- MC straggling^{2}}", "L");
   legSubtractResolution->AddEntry(gSubtractWaterResolution, "#sqrt{- water straggling^{2}}", "L");
   legSubtractResolution->Draw();
   
   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->SetTextSize(0.04);
      title->InsertText(Form("Data filtered, size %d.", filterSize));
      gPad->Modified();
   }

}
