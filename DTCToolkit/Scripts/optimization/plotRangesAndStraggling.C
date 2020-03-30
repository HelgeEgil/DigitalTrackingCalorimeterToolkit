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

// IF FLOAT VALUES OF ABSORBER THICKNESS: USE *10 VALUES (3.5 -> 35)
Int_t absorberThickness = 1;
Bool_t kFilterData = false;
Bool_t kUseCarbon = false;
Int_t filterSize = 3;
const Int_t arraySize = 3500;
const Int_t xFrom = 10;
const Int_t xTo = 350;

void filterArray(Float_t *array, Int_t filterSize) {
   Float_t tempArray[arraySize];
   Float_t value;
   Int_t n;


   for (Int_t i=0; i<arraySize; i++) {
      value = 0;
      n = 0;

      for (Int_t j=i-filterSize/2; j<=i+filterSize/2; j++) {
         if (j<0 || j>=arraySize || array[j] == 0) {
            continue;
         }
         value += array[j];
         n++;
      }

      if (n>0) value /= n;
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
  
   gStyle->SetTitleH(0.1);

   TCanvas *c6 = new TCanvas("c6", "Resolution", 1800, 600);
   c6->Divide(4,1,0.0001,0.0001);

   TCanvas *c6red = new TCanvas("c6red", "Resolution", 500, 500);

   TCanvas *cFourier = new TCanvas("cFourier", "1D fourier of accuracy", 1000, 500);

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
   

   Float_t a_dtc = 0, p_dtc = 0;
   ifstream in;
   in.open(Form("../../Data/Ranges/%dmm_Al_Helium.csv", absorberThickness));

   Double_t dtcRanges[500];
   Double_t dtcEnergies[500];
   Double_t rangesWater[500];
   Double_t energiesWater[500];
   Double_t dtcRange = 0;
   Double_t dtcEnergy = 0;
   Int_t    dtcIdx = 0;

   while (1) {
      in >> dtcEnergy >> dtcRange;
      if (!in.good()) break;

      dtcEnergies[dtcIdx] = dtcEnergy;
      dtcRanges[dtcIdx++] = dtcRange;
   }
   in.close();

   TGraph * range_energy = new TGraph(dtcIdx, dtcEnergies, dtcRanges);
   TF1    * range_energy_fit = new TF1("range_energy_fit", "[0] * pow(x, [1])");
   range_energy_fit->SetParameters(0.02, 1.7);
   range_energy->Fit("range_energy_fit", "M,Q");
   a_dtc = range_energy_fit->GetParameter(0);
   p_dtc = range_energy_fit->GetParameter(1);


   Int_t idxWater = 0;
   Float_t energy, range;
   in.open("../../Data/Ranges/Water.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;

      rangesWater[idxWater] = range;
      energiesWater[idxWater++] = energy * (917/230.);
   }
   in.close();

   TSpline3 * splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);
   TSpline3 * splineWaterInv = new TSpline3("splineWaterInv", rangesWater, energiesWater, idxWater);
   TSpline3 * splineDTC = new TSpline3("splineDTC", dtcEnergies, dtcRanges, dtcIdx);
   TSpline3 * splineDTCInv = new TSpline3("splineDTCInv", dtcRanges, dtcEnergies, dtcIdx);

   if (!kUseCarbon) {
      in.open("../../OutputFiles/findManyRangesDegrader.csv");
   }
   else {
      in.open("../../OutputFiles/findManyRangesDegraderCarbon.csv");
   }

   Int_t nlines0 = 0;
   Float_t a_wtr = 0.02387;
   Float_t p_wtr = 1.7547;

   Float_t wepl_ratio = splineWater->Eval(900) / splineDTC->Eval(900);
   Float_t wepl_ratio0 = a_wtr / a_dtc * pow(225 / a_wtr, 1 - p_dtc / p_wtr);

   cout << "WEPL ratio (spline) = " << wepl_ratio << endl;
   cout << "WEPL ratio (BK) = " << wepl_ratio0 << endl;
   
   Float_t wtr_range = splineWater->Eval(225); // GATE

   while (1) {
      in >>  waterphantomthickness_ >> thickness_ >> nomrange_ >> nomsigma_ >> energy >> dummy0;

      if (!in.good()) {
         break;
      }

      if (thickness_ != absorberThickness) {
         continue;
      }

      //wepl_ratio = splineWater->Eval(energy) / splineDTC->Eval(energy);
//      wepl_ratio = 2.09;

//      arrayMCActualSigma[nlines0] = nomsigma_ * wepl_ratio0;
//      arrayMCActualSigmaRatio[nlines0] = nomsigma_ * 100 * wepl_ratio0 / wtr_range;
//      arrayMCActualResidualRange[nlines0++] = nomrange_ * wepl_ratio0;
      arrayMCActualSigma[nlines0] = nomsigma_ * wepl_ratio;
      arrayMCActualSigmaRatio[nlines0] = nomsigma_ * 100 * wepl_ratio / wtr_range;
      arrayMCActualResidualRange[nlines0++] = nomrange_ * wepl_ratio;
      
   }
   in.close();

   if (!kUseCarbon) {
      in.open("../../OutputFiles/result_makebraggpeakfit_ProtonHelium100k.csv");
   }
   else {
      in.open("../../OutputFiles/result_makebraggpeakfitCarbon.csv");
   }


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

   Float_t DCcorrection = 0;
   for (int i=0; i<nlines; i++) DCcorrection += arrayMCDelta[i];
   DCcorrection /= nlines;

   TGraphErrors *hMCD = new TGraphErrors(nlines, arrayMC, arrayMCDelta, arrayEE, arrayEMC);
   TGraph *pstarD = new TGraph(nlines, arrayMC, arrayPSTARDelta);
   TGraph *pstarminD = new TGraph(nlines, arrayMC, arrayPSTARminDelta);
   TGraph *pstarmaxD = new TGraph(nlines, arrayMC, arrayPSTARmaxDelta);
   TGraph *pstarshadeD = new TGraph(nlines*2);

   for (Int_t i=0; i<nlines; i++) {
      pstarshadeD->SetPoint(i, arrayMC[i], arrayPSTARmaxDelta[i]);
      pstarshadeD->SetPoint(nlines+i, arrayMC[nlines-i-1], arrayPSTARminDelta[nlines-i-1]);
   }

   pstarshadeD->GetXaxis()->SetTitleFont(22);
   pstarshadeD->GetYaxis()->SetTitleFont(22);
   pstarshadeD->GetXaxis()->SetTitleOffset(0.9);
   pstarshadeD->GetYaxis()->SetTitleOffset(0.2);
   pstarshadeD->GetXaxis()->SetLabelFont(22);
   pstarshadeD->GetXaxis()->SetTitleSize(0.05);
   pstarshadeD->GetXaxis()->SetLabelSize(0.05);
   pstarshadeD->GetYaxis()->SetLabelFont(22);
   pstarshadeD->GetYaxis()->SetTitleSize(0.08);
   pstarshadeD->GetYaxis()->SetLabelSize(0.05);

   pstarshadeD->GetYaxis()->SetRangeUser(-10, 10); // 145 - 270

   hMCD->SetMarkerColor(kBlue);
   hMCD->SetMarkerStyle(21);
   hMCD->SetMarkerSize(1);

   pstarD->SetLineWidth(3);
   pstarD->SetLineColor(kMagenta-10);
   pstarshadeD->SetTitle(Form("Reconstructed ranges of proton beams with %d mm Al absorbator;Detector depth [WEPL mm];WET error [mm]", mmAbsorbator));
   
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
   TPaveText *title2 = (TPaveText*) gPad->GetPrimitive("title");
   title2->SetTextFont(22);
   title2->SetTextSize(0.08);
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
//   gResolution->SetTitle(Form("WEPL resolution using %d mm Al absorber;Depth in detector [WEPL mm];#Delta WEPL [mm]", mmAbsorbator));
   gResolution->SetTitle(";Proton range [mm WET];Range straggling [mm WET]");
   gResolution->GetYaxis()->SetRangeUser(3.5,6);
   gResolution->GetYaxis()->SetTitleOffset(1.2);
   gResolution->GetYaxis()->SetDecimals();
   gResolution->GetXaxis()->SetRangeUser(xFrom,xTo);
   gResolution->GetXaxis()->SetTitleFont(22);
   gResolution->GetYaxis()->SetTitleFont(22);
   gResolution->GetXaxis()->SetLabelFont(22);
   gResolution->GetYaxis()->SetLabelFont(22);
   gResolution->GetXaxis()->SetTitleSize(0.05);
   gResolution->GetYaxis()->SetTitleSize(0.05);
   gResolution->GetXaxis()->SetLabelSize(0.05);
   gResolution->GetYaxis()->SetLabelSize(0.05);
   gResolution->SetMarkerColor(kBlue);
   gResolution->SetMarkerStyle(21);
   gResolution->SetMarkerSize(1);
   gResolution->SetLineWidth(3);
   gResolution->SetLineColor(kBlue);
   gResolution->Draw("LA");


   c6->cd(1);
  
   TF1 * fRes = new TF1("fRes", "pol0");
   fRes->SetLineColor(kBlack);
   gResolution->Fit(fRes, "B,Q,N", "", xFrom, xTo);
/*
   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->InsertText(Form("Data filtered, size %d. Mean value = %.2f mm.", filterSize, fRes->GetParameter(0)));
      title->SetTextSize(0.04);
      gPad->Modified();
   }
   */
   
   TGraph *gResolutionStraggling = new TGraph(nlines0, arrayMCActualResidualRange, arrayMCActualSigma);
   gResolutionStraggling->SetLineWidth(3);
   gResolutionStraggling->SetLineColor(kRed);
   gResolutionStraggling->Draw("L");

   gPad->Update();
   // Janni: Straggling in water is 1.063 percent
   // Range in water is 379.4 mm (PSTAR) or 382.57 (Janni)
   // Straggling in water is 3.791 mm as measured in GATE
   // Need new measurement for Helium
   Float_t waterStraggling = 3.791; // GATE
   waterStraggling = 1.84;
   TLine *lWaterStraggling = new TLine(gPad->GetUxmin(), waterStraggling, gPad->GetUxmax(), waterStraggling);
   lWaterStraggling->SetLineWidth(2);
   lWaterStraggling->SetLineColor(kGreen);
   lWaterStraggling->Draw();

   TLegend *legRes = new TLegend(0.17, 0.77, 0.69, 0.88);
   legRes->SetTextSize(0.03);
   legRes->SetTextFont(22);
   legRes->AddEntry(gResolution, "Measured straggling", "L");
//   legRes->AddEntry(fRes, "   + Mean value", "L");
   legRes->AddEntry(gResolutionStraggling, "MC truth straggling", "L");
   legRes->AddEntry(lWaterStraggling, "Straggling in water", "L");
   legRes->Draw();
   
   c6red->cd();
   gResolution->Draw("LA");
   gResolutionStraggling->Draw("L");
   lWaterStraggling->Draw();

   TLatex *t = new TLatex();
   t->SetTextSize(0.045);
   t->SetTextFont(22);
   t->SetTextAngle(16.35);
   t->DrawLatex(82.5, 5.36, "Measured straggling #LT#hat#sigma_{R}#GT");
   t->SetTextAngle(13);
   t->DrawLatex(200.8, 4, "MC truth straggling");
   t->SetTextAngle(0);
   t->DrawLatex(232, 3.679, "Straggling in water");

//   legRes->Draw();

   c6->cd(2);
   TGraph *gResolutionRatio = new TGraph(nlines, arrayMC, arrayEMCRatio);
   gResolutionRatio->SetTitle(Form("WEPL resolution using %d mm Al absorber;Depth in detector [WEPL mm];#Delta WEPL / WEPL [%]", mmAbsorbator));
   gResolutionRatio->GetYaxis()->SetTitleOffset(1.35);
   gResolutionRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
   gResolutionRatio->GetXaxis()->SetRangeUser(xFrom, xTo);
   gResolutionRatio->SetMarkerColor(kBlue);
   gResolutionRatio->SetLineColor(kBlue);
   gResolutionRatio->SetLineWidth(3);
   gResolutionRatio->SetMarkerStyle(21);
   gResolutionRatio->SetMarkerSize(1);
   gResolutionRatio->Draw("LA");
   
   TF1 * fResRatio = new TF1("fResRatio", "pol0");
   fResRatio->SetLineColor(kBlack);
   gResolutionRatio->Fit(fResRatio, "B,Q", "", xFrom, xTo);

   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->SetTextSize(0.04);
      title->InsertText(Form("Data filtered, size %d. Mean value = %.2f.", filterSize, fResRatio->GetParameter(0)));
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
   //waterStraggling = 3.791 * 100 / 378.225; // GATE
   TLine *lWaterStraggling2 = new TLine(gPad->GetUxmin(), waterStraggling, gPad->GetUxmax(), waterStraggling);
   lWaterStraggling2->SetLineWidth(2);
   lWaterStraggling2->SetLineColor(kGreen);
   lWaterStraggling2->Draw();
   
   TGraph *gResolutionStragglingRatio = new TGraph(nlines0, arrayMCActualResidualRange, arrayMCActualSigmaRatio);
   gResolutionStragglingRatio->SetLineWidth(3);
   gResolutionStragglingRatio->SetLineColor(kRed);
   gResolutionStragglingRatio->Draw("L");
   
   TLegend *legResRatio = new TLegend(0.36, 0.76, 0.88, 0.88);
   legResRatio->SetTextSize(0.03);
   legResRatio->SetTextFont(22);
   legResRatio->AddEntry(gResolutionRatio, "Measured range spread", "L");
   legResRatio->AddEntry(fResRatio, "   + Mean value", "L");
   legResRatio->AddEntry(gResolutionStragglingRatio, "MC truth range straggling", "L");
   legResRatio->AddEntry(lWaterStraggling2, "Straggling in water", "L");
   legResRatio->Draw();

   c6->cd(3);
   Int_t idx=0;
   for (Int_t i=0; i<nlines0; i++) {
      if (arrayEMCRatio[i]>0) arrayStragglingsRatio[idx++] = arrayEMCRatio[i] / arrayMCActualSigmaRatio[i];
   }

   TGraph *gStragglingsRatio = new TGraph(idx, arrayMCActualResidualRange, arrayStragglingsRatio);
   gStragglingsRatio->SetTitle(Form("Resolution/straggling, %d mm Al absorber;Depth in detector [WEPL mm]; Resolution / straggling", mmAbsorbator));
   gStragglingsRatio->GetYaxis()->SetRangeUser(0.3, 1.4);
   gStragglingsRatio->GetYaxis()->SetTitleOffset(1.2);
   gStragglingsRatio->GetXaxis()->SetRangeUser(xFrom, xTo);
   gStragglingsRatio->SetLineWidth(3);
   gStragglingsRatio->SetLineColor(kBlue);
   gStragglingsRatio->Draw("AL");
   
   TF1 * fStragRatio = new TF1("fStragRatio", "pol0");
   fStragRatio->SetLineColor(kBlack);
   gStragglingsRatio->Fit(fStragRatio, "B,Q", "", xFrom, xTo);
   
   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->SetTextSize(0.04);
      title->InsertText(Form("Data filtered, size %d. Mean excess = %.2f %%", filterSize, 100*(fStragRatio->GetParameter(0)-1)));
      gPad->Modified();
   }


   gPad->Update();
   TLine *stragglingline = new TLine(gPad->GetUxmin(), 1, gPad->GetUxmax(), 1);
   stragglingline->SetLineStyle(7);
   stragglingline->Draw("same");

   Float_t subtractMCStraggling[arraySize];
   Float_t subtractWaterStraggling[arraySize];
   // waterStraggling = 3.791;

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
   gSubtractFullResolution->GetYaxis()->SetTitleOffset(1.2);
   gSubtractFullResolution->GetXaxis()->SetRangeUser(xFrom, xTo);
   gSubtractFullResolution->SetLineColor(kBlue);
   gSubtractFullResolution->SetLineWidth(3);
   gSubtractMCResolution->SetLineColor(kRed);
   gSubtractMCResolution->SetLineWidth(3);
   gSubtractWaterResolution->SetLineColor(kGreen);
   gSubtractWaterResolution->SetLineWidth(3);

   
   TF1 * fResMC = new TF1("fResMC", "pol0");
   fResMC->SetLineColor(kBlack);
   gSubtractMCResolution->Fit(fResMC, "B,Q", "", xFrom, xTo);

   gSubtractFullResolution->Draw("LA");
   gSubtractMCResolution->Draw("same");
   gSubtractWaterResolution->Draw("same");

   TLegend *legSubtractResolution = new TLegend(0.41, 0.73, 0.86, 0.88);
   legSubtractResolution->SetTextSize(0.03);
   legSubtractResolution->SetTextFont(22);
   legSubtractResolution->AddEntry(gSubtractFullResolution, "Range resolution", "L");
   legSubtractResolution->AddEntry(gSubtractMCResolution, "#sqrt{- MC straggling^{2}}", "L");
   legSubtractResolution->AddEntry(fResMC, "   + Mean value", "L");
   legSubtractResolution->AddEntry(gSubtractWaterResolution, "#sqrt{- water straggling^{2}}", "L");
   legSubtractResolution->Draw();
   if (kFilterData) {
      gPad->Update();
      TPaveText * title = (TPaveText *)gPad->FindObject("title");
      title->SetTextSize(0.04);
      title->InsertText(Form("Data filtered, size %d. Mean value = %.2f mm.", filterSize, fResMC->GetParameter(0)));
      gPad->Modified();
   }
   
   cFourier->Divide(2,1,0.0000001,0.000001);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.05);
   gPad->SetTopMargin(0.05);
   gPad->SetBottomMargin(0.15);
   gPad->Update();

   cFourier->cd(1);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.05);
   gPad->SetTopMargin(0.05);
   gPad->SetBottomMargin(0.12);
   Float_t from = arrayMC[0];
   Float_t to = arrayMC[nlines-1];
   // Want 4pi/sqrt(nlines)
   // from = 0
   // to = 4pi/sqrt(nlines)

   DCcorrection = 0;

   Float_t  RANGE = 50;
   Int_t    PAD = 500;
   TH1D *hRangeError = new TH1D("hRangeError", Form("Range Error Histogram (%d mm Al)",absorberThickness), 250+2*PAD, 50-PAD, 300+PAD);
   if (absorberThickness>10) hRangeError->SetTitle(Form("Range Error Histogram (%.1f mm Al)", float(absorberThickness)/10));

   printf("Making 251 bins from %.2f to %.2f\n", arrayMC[308], arrayMC[58]);
   for (int i=308; i>57; i--) {
      hRangeError->Fill(RANGE, arrayMCDelta[i] - DCcorrection);// - ( 0.333 - 0.00454799 + 0.0567 - (0.003867 +0.000814089)*RANGE + (4.012e-6 + 2.32894e-06) * pow(RANGE,2) ));
      RANGE++;
   }
   hRangeError->GetXaxis()->SetRangeUser(50,300);
   hRangeError->SetLineColor(kBlack);
//   hRangeError->SetFillColor(kOrange-3);
   hRangeError->SetTitle(";Range [mm WEPL];Range deviation [mm WEPL]");
   hRangeError->Draw();
   hRangeError->GetXaxis()->SetLabelFont(22);
   hRangeError->GetXaxis()->SetTitleFont(22);
   hRangeError->GetYaxis()->SetLabelFont(22);
   hRangeError->GetYaxis()->SetTitleFont(22);
   hRangeError->GetXaxis()->SetLabelSize(0.05);
   hRangeError->GetXaxis()->SetTitleSize(0.05);
   hRangeError->GetYaxis()->SetLabelSize(0.05);
   hRangeError->GetYaxis()->SetTitleSize(0.05);
   hRangeError->GetYaxis()->SetTitleOffset(1.4);

   cFourier->cd(2);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.05);
   gPad->SetTopMargin(0.05);
   gPad->SetBottomMargin(0.12);
   TH1 *hm = 0;
   TVirtualFFT::SetTransform(0);
   hm = hRangeError->FFT(hm, "MAG");
   hm->Scale(1/sqrt(250+2*PAD));
   hm->GetXaxis()->SetRangeUser(10,250);
   hm->SetTitle(";Frequency [A.U.];Fourier magnitude [A.U.]");
   hm->SetFillColor(kRed-4);
   hm->SetLineColor(kBlack);
   hm->GetXaxis()->SetLabelFont(22);
   hm->GetXaxis()->SetTitleFont(22);
   hm->GetYaxis()->SetLabelFont(22);
   hm->GetYaxis()->SetTitleFont(22);
   hm->GetXaxis()->SetLabelSize(0.05);
   hm->GetXaxis()->SetTitleSize(0.05);
   hm->GetYaxis()->SetLabelSize(0.05);
   hm->GetYaxis()->SetTitleSize(0.05);
   hm->GetYaxis()->SetTitleOffset(1);
//   hm->GetYaxis()->SetRangeUser(0, 5);
   hm->Draw();
   gPad->Update();

   TLegend *fLeg = new TLegend(0.41, 0.83, 0.94, 0.92);
   if (absorberThickness < 10) {
      fLeg->AddEntry(hm, Form("%d mm absorber thickness", absorberThickness), "F");
   }
   else {
      fLeg->AddEntry(hm, Form("%.1f mm absorber thickness", float(absorberThickness)/10), "F");
   }
   fLeg->SetTextFont(22);
   fLeg->SetTextSize(0.04);
   fLeg->Draw();

}
