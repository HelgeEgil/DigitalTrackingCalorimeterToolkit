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
#include <TStyle.h>
#include <TSpline.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

using namespace std;

Bool_t kUseCarbon = false;
Bool_t kUseDCCorrection = true;
Bool_t kUseNominalValuesFromRMBPF = true;
Bool_t kUseWEPL = false;

const Int_t arraySize = 1500;
const Int_t xFrom = 40;

void plotEnergyVsRange() {
   TCanvas *c1 = new TCanvas("c1", "Range accuracy", 1100, 800);
   
   
   TPaveLabel *Ytitle = new TPaveLabel(0.01, 0.05, 0.03, 0.9, "Range deviation [mm WEPL]");
   Ytitle->SetBorderSize(0);
   Ytitle->SetFillColor(kWhite);
   Ytitle->SetTextAngle(90);
   Ytitle->SetTextFont(22);
   Ytitle->SetTextSize(0.04);
   Ytitle->Draw();
   
   TPaveLabel *Xtitle = new TPaveLabel(0.5, 0.05, 0.95, 0.1, "Proton range [mm WET]");
   Xtitle->SetBorderSize(0);
   Xtitle->SetFillColor(kWhite);
   Xtitle->SetTextFont(22);
   Xtitle->SetTextSize(0.6);
   Xtitle->Draw();
   
   gStyle->SetPadLeftMargin(0.035);
   gStyle->SetPadRightMargin(0.025);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadTopMargin(0.23);

   TPad *graphPad = new TPad("Graphs", "Graphs", 0.05, 0.1, 0.95, 0.95);
   graphPad->Draw();
   graphPad->cd();
   graphPad->Divide(1,5,0.00001,0.00001);

   Float_t  arrayE2[arraySize] = {0}; // energy MC
   Float_t  arrayE3[arraySize] = {0}; // energy MC
   Float_t  arrayE35[arraySize] = {0}; // energy MC
   Float_t  arrayE4[arraySize] = {0}; // energy MC
   Float_t  arrayE5[arraySize] = {0}; // energy MC
   Float_t  arrayE6[arraySize] = {0}; // energy MC
   Float_t  arrayMC2[arraySize] = {0}; // range MC
   Float_t  arrayMC3[arraySize] = {0}; // range MC
   Float_t  arrayMC35[arraySize] = {0}; // range MC
   Float_t  arrayMC4[arraySize] = {0}; // range MC
   Float_t  arrayMC5[arraySize] = {0}; // range MC
   Float_t  arrayMC6[arraySize] = {0}; // range MC
   Double_t  arrayR2[arraySize] = {0}; // truth range MC
   Double_t  arrayR3[arraySize] = {0}; // truth range MC
   Double_t  arrayR4[arraySize] = {0}; // truth range MC
   Double_t  arrayR5[arraySize] = {0}; // truth range MC
   Double_t  arrayR6[arraySize] = {0}; // truth range MC
   Double_t  arrayRD2[arraySize] = {0}; // truth range MC
   Double_t  arrayRD3[arraySize] = {0}; // truth range MC
   Double_t  arrayRD4[arraySize] = {0}; // truth range MC
   Double_t  arrayRD5[arraySize] = {0}; // truth range MC
   Double_t  arrayRD6[arraySize] = {0}; // truth range MC
   Double_t  arrayRE2[arraySize] = {0}; // truth range MC
   Double_t  arrayRE3[arraySize] = {0}; // truth range MC
   Double_t  arrayRE4[arraySize] = {0}; // truth range MC
   Double_t  arrayRE5[arraySize] = {0}; // truth range MC
   Double_t  arrayRE6[arraySize] = {0}; // truth range MC

   Double_t energies[arraySize] = {0};
   Double_t thicknesses[arraySize] = {0};
   Double_t energiesWater[arraySize] = {0};
   Double_t rangesWater[arraySize] = {0};

   Float_t correction_2 = 0.02, correction_3 = 0.415, correction_35 = 0.74, correction_4 = 0.74, correction_5 = 1.39 , correction_6 = 2.04;

   gStyle->SetOptStat(0);

   Float_t nomrange_, estrange_, sigmaRange_, lastRange_, nomsigma_, waterphantomthickness_, dummy0;
   Int_t energy_, thickness_, degrader_;
   Float_t estimatedStraggling;

   ifstream in0;
   Double_t energy, range;
   Int_t idxWater = 0;
   in0.open("../../Data/Ranges/Water.csv");
   while (1) {
      in0 >> energy >> range;
      if (!in0.good()) break;

      energiesWater[idxWater] = energy;
      rangesWater[idxWater++] = range;
   }

   TSpline3 * waterSpline = new TSpline3("waterSpline", energiesWater, rangesWater, idxWater); // get WEPL from ENERGY


   ifstream in1;
   in1.open("../../Data/Ranges/EnergyAfterDegraderPSTAR.csv");
   Int_t thick, n=0;
   while (1) {
      in1 >> thick >> energy;
      if (!in1.good()) break;   
      thicknesses[n] = thick;
      energies[n++] = energy;
   }
   in1.close();
   printf("Found %d lines in EnergyAfterDegraderPSTAR.csv\n", n);

   TSpline3 *energySpline = new TSpline3("energySpline", thicknesses, energies, n);

   Int_t nlinesR6 = 0, nlinesR2 = 0, nlinesR3 = 0, nlinesR4 = 0, nlinesR5 = 0;
   Float_t dummy, floatenergy_;
   ifstream in2;
   in2.open("../../OutputFiles/findManyRangesDegrader.csv");
   while (1) {
      in2 >> degrader_ >> thickness_ >> nomrange_ >> dummy >> dummy >> floatenergy_ >> dummy;

      if (!in2.good()) break;

      if (thickness_ == 2) {
         arrayR2[nlinesR2] = nomrange_;
         arrayRE2[nlinesR2] = floatenergy_;
         arrayRD2[nlinesR2++] = degrader_;
      }
      else if (thickness_ == 3) {
         arrayR3[nlinesR3] = nomrange_;
         arrayRE3[nlinesR3] = floatenergy_;
         arrayRD3[nlinesR3++] = degrader_;
      }
      else if (thickness_ == 4) {
         arrayR4[nlinesR4] = nomrange_;
         arrayRE4[nlinesR4] = floatenergy_;
         arrayRD4[nlinesR4++] = degrader_;
      }
      else if (thickness_ == 5) {
         arrayR5[nlinesR5] = nomrange_;
         arrayRE5[nlinesR5] = floatenergy_;
         arrayRD5[nlinesR5++] = degrader_;
      }
      else if (thickness_ == 6) {
         arrayR6[nlinesR6] = nomrange_;
         arrayRE6[nlinesR6] = floatenergy_;
         arrayRD6[nlinesR6++] = degrader_;
      }
   }
   in2.close();
/*
   Double_t arrayR2reverse[arraySize] = {};
   Double_t arrayR3reverse[arraySize] = {};
   Double_t arrayR4reverse[arraySize] = {};
   Double_t arrayR5reverse[arraySize] = {};
   Double_t arrayR6reverse[arraySize] = {};
   Double_t arrayRE2reverse[arraySize] = {};
   Double_t arrayRE3reverse[arraySize] = {};
   Double_t arrayRE4reverse[arraySize] = {};
   Double_t arrayRE5reverse[arraySize] = {};
   Double_t arrayRE6reverse[arraySize] = {};
   for (int i=0; i<nlinesR2; i++) arrayR2reverse[i] = arrayR2[nlinesR2-1-i];
   for (int i=0; i<nlinesR2; i++) arrayR3reverse[i] = arrayR2[nlinesR3-1-i];
   for (int i=0; i<nlinesR2; i++) arrayR4reverse[i] = arrayR2[nlinesR4-1-i];
   for (int i=0; i<nlinesR2; i++) arrayR5reverse[i] = arrayR2[nlinesR5-1-i];
   for (int i=0; i<nlinesR2; i++) arrayR6reverse[i] = arrayR2[nlinesR6-1-i];
   for (int i=0; i<nlinesR2; i++) arrayRE2reverse[i] = arrayRE2[nlinesR2-1-i];
   for (int i=0; i<nlinesR2; i++) arrayRE3reverse[i] = arrayRE2[nlinesR3-1-i];
   for (int i=0; i<nlinesR2; i++) arrayRE4reverse[i] = arrayRE2[nlinesR4-1-i];
   for (int i=0; i<nlinesR2; i++) arrayRE5reverse[i] = arrayRE2[nlinesR5-1-i];
   for (int i=0; i<nlinesR2; i++) arrayRE6reverse[i] = arrayRE2[nlinesR6-1-i];

   TSpline3 * energySpline2 = new TSpline3("energySpline2", arrayR2reverse, arrayRE2reverse, nlinesR2); // get ENERGY from RANGE
   TSpline3 * energySpline3 = new TSpline3("energySpline3", arrayR3reverse, arrayRE3reverse, nlinesR3); // get ENERGY from RANGE
   TSpline3 * energySpline4 = new TSpline3("energySpline4", arrayR4reverse, arrayRE4reverse, nlinesR4); // get ENERGY from RANGE
   TSpline3 * energySpline5 = new TSpline3("energySpline5", arrayR5reverse, arrayRE5reverse, nlinesR5); // get ENERGY from RANGE
   TSpline3 * energySpline6 = new TSpline3("energySpline6", arrayR6reverse, arrayRE6reverse, nlinesR6); // get ENERGY from RANGE

   printf("linesR2 = %d, R3 = %d, R4 = %d, R5 = %d, R6 = %d.\n", nlinesR2, nlinesR3, nlinesR4, nlinesR5, nlinesR6);

   std::unordered_map<int, int> lookup2;
   std::unordered_map<int, int> lookup3;
   std::unordered_map<int, int> lookup4;
   std::unordered_map<int, int> lookup5;
   std::unordered_map<int, int> lookup6;
   for (int index = 0; index < nlinesR2; index++) {
      lookup2[arrayRD2[index]] = index;
   }
   for (int index = 0; index < nlinesR3; index++) lookup3[arrayRD3[index]] = index;
   for (int index = 0; index < nlinesR4; index++) lookup4[arrayRD4[index]] = index;
   for (int index = 0; index < nlinesR5; index++) lookup5[arrayRD5[index]] = index;
   for (int index = 0; index < nlinesR6; index++) lookup6[arrayRD6[index]] = index;
  */

   ifstream in;
   if (!kUseCarbon) {
      in.open("../../OutputFiles/result_makebraggpeakfit_proj.csv");
   }
   else {
      in.open("../../OutputFiles/result_makebraggpeakfitCarbon.csv");
   }

   Int_t nlines6 = 0, nlines2 = 0, nlines3 = 0, nlines35 = 0, nlines4 = 0, nlines5 = 0;
   

   while (1) {
      in >> thickness_ >> energy_ >> nomrange_ >> estrange_ >> nomsigma_ >> sigmaRange_;
      // energy_ here is actually degraderthickness

      if (!in.good()) {
         break;
      }

      energy = energySpline->Eval(energy_);

      if (kUseNominalValuesFromRMBPF) {
         if (thickness_ == 2) arrayE2[nlines2] = nomrange_;
         if (thickness_ == 3) arrayE3[nlines3] = nomrange_;
         if (thickness_ == 35) arrayE35[nlines35] = nomrange_;
         if (thickness_ == 4) arrayE4[nlines4] = nomrange_;
         if (thickness_ == 5) arrayE5[nlines5] = nomrange_;
         if (thickness_ == 6) arrayE6[nlines6] = nomrange_;

         if (thickness_ == 2) arrayMC2[nlines2++] = -nomrange_ + estrange_ + correction_2 * kUseDCCorrection;
         if (thickness_ == 3) arrayMC3[nlines3++] = estrange_ - nomrange_ + correction_3 * kUseDCCorrection;
         if (thickness_ == 35) arrayMC35[nlines35++] = estrange_ - nomrange_ + correction_35 * kUseDCCorrection;
         if (thickness_ == 4) arrayMC4[nlines4++] = -nomrange_ + estrange_ + correction_4 * kUseDCCorrection;
         if (thickness_ == 5) arrayMC5[nlines5++] = estrange_ - nomrange_ + correction_5 * kUseDCCorrection;
         if (thickness_ == 6) arrayMC6[nlines6++] = -nomrange_ + estrange_ + correction_6 * kUseDCCorrection;
      }
      else {
         pass;
         /*
         if (thickness_ == 2) {
            int index = lookup2[energy_];
            float tl = arrayR2[index];
            double e = energySpline2->Eval(tl);
            double wepl = (kUseWEPL) ? waterSpline->Eval(e) : tl;
            printf("tl = %.2f mm; energy = %.2f MeV; wepl = %.2f mm.\n", tl, e, wepl);
            arrayE2[nlines2] = wepl;
            arrayMC2[nlines2++] = estrange_ - wepl + correction_2 * kUseDCCorrection;
         }
         else if (thickness_ == 3) {
            int index = lookup3[energy_];
            float tl = arrayR3[index];
            double e = energySpline3->Eval(tl);
            double wepl = (kUseWEPL) ? waterSpline->Eval(e) : tl;
            arrayE3[nlines3] = wepl;
            arrayMC3[nlines3++] = estrange_ - wepl + correction_3 * kUseDCCorrection;
         }
         else if (thickness_ == 4) {
            int index = lookup4[energy_];
            float tl = arrayR4[index];
            double e = energySpline4->Eval(tl);
            double wepl = (kUseWEPL) ? waterSpline->Eval(e) : tl;
            arrayE4[nlines4] = wepl;
            arrayMC4[nlines4++] = estrange_ - wepl + correction_4 * kUseDCCorrection;
         }
         else if (thickness_ == 5) {
            int index = lookup5[energy_];
            float tl = arrayR5[index];
            double e = energySpline5->Eval(tl);
            double wepl = (kUseWEPL) ? waterSpline->Eval(e) : tl;
            arrayE5[nlines5] = wepl;
            arrayMC5[nlines5++] = estrange_ - wepl + correction_5 * kUseDCCorrection;
         }
         else if (thickness_ == 6) {
            int index = lookup6[energy_];
            float tl = arrayR6[index];
            double e = energySpline6->Eval(tl);
            double wepl = (kUseWEPL) ? waterSpline->Eval(e) : tl;
            arrayE6[nlines6] = wepl;
            arrayMC6[nlines6++] = estrange_ - wepl + correction_6 * kUseDCCorrection;
         }
      */
      }
   }
   
   in.close();

   printf("Found the following number of lines for the different geometries:\n2 mm: %d lines\n3 mm: %d lines\n3.5 mm: %d lines\n4 mm: %d lines\n5 mm: %d lines\n6 mm: %d lines\n", nlines2, nlines3, nlines35, nlines4, nlines5, nlines6);

   TGraph *hMC2 = new TGraph(nlines2, arrayE2, arrayMC2);
   TGraph *hMC3 = new TGraph(nlines3, arrayE3, arrayMC3);
   TGraph *hMC35 = new TGraph(nlines35, arrayE35, arrayMC35);
   TGraph *hMC4 = new TGraph(nlines4, arrayE4, arrayMC4);
   TGraph *hMC5 = new TGraph(nlines5, arrayE5, arrayMC5);
   TGraph *hMC6 = new TGraph(nlines6, arrayE6, arrayMC6);

   hMC2->SetTitle(";Range [mm WEPL];");
   hMC3->SetTitle(";Range [mm WEPL];");
   hMC35->SetTitle(";Range [mm WEPL];");
   hMC4->SetTitle(";Range [mm WEPL];");
   hMC5->SetTitle(";Range [mm WEPL];");
   hMC6->SetTitle(";Range [mm WEPL];");

   Float_t ysize = 0.18;

   hMC6->GetXaxis()->SetLabelSize(ysize);
   hMC6->GetXaxis()->SetLabelFont(22);
   hMC4->GetYaxis()->SetLabelSize(ysize);
   hMC4->GetYaxis()->SetLabelFont(22);
   hMC2->GetYaxis()->SetLabelSize(ysize);
   hMC2->GetYaxis()->SetLabelFont(22);
   hMC6->GetYaxis()->SetLabelSize(ysize);
   hMC6->GetYaxis()->SetLabelFont(22);
   hMC3->GetYaxis()->SetLabelSize(ysize);
   hMC3->GetYaxis()->SetLabelFont(22);
   hMC35->GetYaxis()->SetLabelSize(ysize);
   hMC35->GetYaxis()->SetLabelFont(22);
   hMC5->GetYaxis()->SetLabelSize(ysize);
   hMC5->GetYaxis()->SetLabelFont(22);
   hMC4->GetXaxis()->SetLabelSize(ysize);
   hMC4->GetXaxis()->SetLabelFont(22);
   hMC2->GetXaxis()->SetLabelSize(ysize);
   hMC2->GetXaxis()->SetLabelFont(22);
   hMC3->GetXaxis()->SetLabelSize(ysize);
   hMC3->GetXaxis()->SetLabelFont(22);
   hMC35->GetXaxis()->SetLabelSize(ysize);
   hMC35->GetXaxis()->SetLabelFont(22);
   hMC5->GetXaxis()->SetLabelSize(ysize);
   hMC5->GetXaxis()->SetLabelFont(22);
   hMC2->GetXaxis()->SetTitleOffset(3);
   hMC3->GetXaxis()->SetTitleOffset(3);
   hMC4->GetXaxis()->SetTitleOffset(3);
   hMC5->GetXaxis()->SetTitleOffset(3);
   hMC6->GetXaxis()->SetTitleOffset(3);
   hMC35->GetXaxis()->SetTitleOffset(3);

   hMC2->SetLineColor(kRed+4);
   hMC3->SetLineColor(kRed+3);
   hMC35->SetLineColor(kRed+2);
   hMC4->SetLineColor(kRed+1);
   hMC5->SetLineColor(kRed);
   hMC6->SetLineColor(kRed-1);
   hMC2->SetLineWidth(3);
   hMC3->SetLineWidth(3);
   hMC35->SetLineWidth(3);
   hMC4->SetLineWidth(3);
   hMC5->SetLineWidth(3);
   hMC6->SetLineWidth(3);

   Float_t yfrom = -2.5;
   Float_t yto = 2.5;

   Float_t xfrom = 0;
   Float_t xto = 380;

   hMC2->GetXaxis()->SetRangeUser(xfrom, xto);
   hMC3->GetXaxis()->SetRangeUser(xfrom, xto);
   hMC35->GetXaxis()->SetRangeUser(xfrom, xto);
   hMC4->GetXaxis()->SetRangeUser(xfrom, xto);
   hMC5->GetXaxis()->SetRangeUser(xfrom, xto);
   hMC6->GetXaxis()->SetRangeUser(xfrom, xto);
   hMC2->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC3->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC35->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC4->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC5->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC6->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC2->GetYaxis()->SetNdivisions(404);
   hMC3->GetYaxis()->SetNdivisions(404);
   hMC35->GetYaxis()->SetNdivisions(404);
   hMC4->GetYaxis()->SetNdivisions(404);
   hMC5->GetYaxis()->SetNdivisions(404);
   hMC6->GetYaxis()->SetNdivisions(404);

   if (!kUseDCCorrection) {
      hMC2->GetYaxis()->SetRangeUser(yfrom - correction_2, yto - correction_2);
      hMC3->GetYaxis()->SetRangeUser(yfrom - correction_3, yto - correction_3);
      hMC35->GetYaxis()->SetRangeUser(yfrom - correction_35, yto - correction_35);
      hMC4->GetYaxis()->SetRangeUser(yfrom - correction_4, yto - correction_4);
      hMC5->GetYaxis()->SetRangeUser(yfrom - correction_5, yto - correction_5);
      hMC6->GetYaxis()->SetRangeUser(yfrom - correction_6, yto - correction_6);
   }

   Float_t textX = 7.22;
   Float_t textY = 3;

   graphPad->cd(1);
   gPad->SetGridy();
   hMC2->Draw("LA");
   TText *t2 = new TText();
   t2->SetTextSize(0.15);
   t2->SetTextFont(22);
   t2->DrawText(textX, textY - correction_2 * (!kUseDCCorrection), Form("2 mm Al absorber (calibration constant = +%.1f mm)", correction_2)); 

   graphPad->cd(2);
   gPad->SetGridy();
   hMC3->Draw("LA");
   t2->DrawText(textX, textY - correction_3 * (!kUseDCCorrection), Form("3 mm Al absorber (calibration constant = +%.1f mm)", correction_3));
   /*
   graphPad->cd(3);
   gPad->SetGridy();
   hMC35->Draw("LA");
   t2->DrawText(textX, textY, "3.5 mm Al absorber");
*/
   graphPad->cd(3);
   gPad->SetGridy();
   hMC4->Draw("LA");
   t2->DrawText(textX, textY - correction_4 * (!kUseDCCorrection), Form("4 mm Al absorber (calibration constant = +%.1f mm)", correction_4));
   
   graphPad->cd(4);
   gPad->SetGridy();
   hMC5->Draw("LA");
   t2->DrawText(textX, textY - correction_5 * (!kUseDCCorrection), Form("5 mm Al absorber (calibration constant = +%.1f mm)", correction_5));
   
   graphPad->cd(5);
   gPad->SetGridy();
   hMC6->Draw("LA");
   t2->DrawText(textX, textY - correction_6 * (!kUseDCCorrection), Form("6 mm Al absorber (calibration constant = +%.1f mm)", correction_6));


}
