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
   
   TPaveLabel *Xtitle = new TPaveLabel(0.5, 0.05, 0.95, 0.1, "Phantom thickness [mm WEPL]");
   Xtitle->SetBorderSize(0);
   Xtitle->SetFillColor(kWhite);
   Xtitle->SetTextFont(22);
   Xtitle->SetTextSize(0.6);
   Xtitle->Draw();
   
   gStyle->SetPadLeftMargin(0.035);
   gStyle->SetPadRightMargin(0.025);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadTopMargin(0.23);
   

   //   TCanvas *cRangeAccuracy = new TCanvas("cRangeAccuracy", "Range accuracy", 1200, 300);

   TPad *graphPad = new TPad("Graphs", "Graphs", 0.05, 0.1, 0.95, 0.95);
   graphPad->Draw();
   graphPad->cd();
   graphPad->Divide(1,2,0.00001,0.00001);
   
   TCanvas *c2 = new TCanvas("c2", "Range straggling", 1100, 800);
   TPaveLabel *Ytitle2 = new TPaveLabel(0.01, 0.05, 0.03, 0.9, "Meas. range straggling (#pm2#sigma) [mm WEPL]");
   Ytitle2->SetBorderSize(0);
   Ytitle2->SetFillColor(kWhite);
   Ytitle2->SetTextAngle(90);
   Ytitle2->SetTextFont(22);
   Ytitle2->SetTextSize(0.04);
   Ytitle2->Draw();
   
   TPaveLabel *Xtitle2 = new TPaveLabel(0.5, 0.05, 0.95, 0.1, "Phantom thickness [mm WEPL]");
   Xtitle2->SetBorderSize(0);
   Xtitle2->SetFillColor(kWhite);
   Xtitle2->SetTextFont(22);
   Xtitle2->SetTextSize(0.6);
   Xtitle2->Draw();

   //   TCanvas *cRangeAccuracy = new TCanvas("cRangeAccuracy", "Range accuracy", 1200, 300);

   TPad *graphPad2 = new TPad("Graphs", "Graphs", 0.05, 0.1, 0.95, 0.95);
   graphPad2->Draw();
   graphPad2->cd();
   graphPad2->Divide(1,2,0.00001,0.00001);

   Double_t  arrayE2[arraySize] = {0}; // energy MC
   Double_t  arrayE3[arraySize] = {0}; // energy MC
   Double_t  arrayE35[arraySize] = {0}; // energy MC
   Double_t  arrayE4[arraySize] = {0}; // energy MC
   Double_t  arrayE5[arraySize] = {0}; // energy MC
   Double_t  arrayE6[arraySize] = {0}; // energy MC
   Double_t  arrayMC2[arraySize] = {0}; // range MC
   Double_t  arrayMC3[arraySize] = {0}; // range MC
   Double_t  arrayMC35[arraySize] = {0}; // range MC
   Double_t  arrayMC4[arraySize] = {0}; // range MC
   Double_t  arrayMC5[arraySize] = {0}; // range MC
   Double_t  arrayMC6[arraySize] = {0}; // range MC
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
   Double_t  arrayRS35[arraySize] = {0}; // Range Straggling
   Double_t  arrayRS4[arraySize] = {0}; // Range Straggling

   Double_t energies[arraySize] = {0};
   Double_t thicknesses[arraySize] = {0};
   Double_t energiesWater[arraySize] = {0};
   Double_t rangesWater[arraySize] = {0};

   Float_t correction_2 = 0.02, correction_3 = 0.415, correction_35 = 0.82, correction_4 = 1.45, correction_5 = 1.39 , correction_6 = 2.04;

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
   in1.open("../../Data/Ranges/EnergyAfterDegraderProton.csv");
   Int_t thick, n=0;
   while (1) {
      in1 >> thick >> energy;
      if (!in1.good()) break;   
      thicknesses[n] = thick;
      energies[n++] = energy;
   }
   in1.close();
   printf("Found %d lines in EnergyAfterDegrader.csv\n", n);


   Int_t nlinesR6 = 0, nlinesR2 = 0, nlinesR3 = 0, nlinesR4 = 0, nlinesR5 = 0;
   Float_t dummy, floatenergy_, sigmarange_;
   ifstream in2;
   in2.open("../../OutputFiles/findManyRangesDegraderHelium.csv");
   while (1) {
      in2 >> degrader_ >> thickness_ >> nomrange_ >> dummy >> dummy >> floatenergy_ >> dummy;

      if (!in2.good()) break;

      if (thickness_ == 2) {
         arrayR2[nlinesR2] = nomrange_;
         arrayRE2[nlinesR2] = floatenergy_;
         arrayRD2[nlinesR2++] = degrader_;
      }
      else if (thickness_ == 1) {
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
   
   ifstream in;
   if (!kUseCarbon) {
      in.open("../../OutputFiles/result_makebraggpeakfit.csv");
   }
   else {
      in.open("../../OutputFiles/result_makebraggpeakfitCarbon.csv");
   }

   Int_t nlines6 = 0, nlines2 = 0, nlines3 = 0, nlines35 = 0, nlines4 = 0, nlines5 = 0;

   while (1) {
      in >> thickness_ >> energy_ >> nomrange_ >> estrange_ >> nomsigma_ >> sigmaRange_;
      
      if (!in.good()) {
         break;
      }

      if (thickness_ == 0) arrayE35[nlines35] = 333.7 - nomrange_;  // Proton
      if (thickness_ == 0) arrayRS35[nlines35] = sigmaRange_ * 4;  // Proton
      if (thickness_ == 1) arrayE4[nlines4] = 332.3 - nomrange_; // Helium
      if (thickness_ == 1) arrayRS4[nlines4] = sigmaRange_ * 4; // Helium
      if (thickness_ == 0) arrayMC35[nlines35++] = estrange_ - nomrange_ + (-0.41 + 0.0024 * estrange_) * kUseDCCorrection; // Proton 
      if (thickness_ == 1) arrayMC4[nlines4++] = -nomrange_ + estrange_ + correction_4 * kUseDCCorrection; // Helium
   }
   
   in.close();

   printf("Found the following number of lines for the different geometries:\n2 mm: %d lines\n3 mm: %d lines\n3.5 mm: %d lines\n4 mm: %d lines\n5 mm: %d lines\n6 mm: %d lines\n", nlines2, nlines3, nlines35, nlines4, nlines5, nlines6);

  /*
   gStyle->SetTitleFont(22, "xy");
   gStyle->SetTitleFont(22, "t");
   gStyle->SetLabelFont(22, "xy");
   gStyle->SetTitleSize(0.06, "xy");
   gStyle->SetLabelSize(0.1, "xy");
   gStyle->SetTitleOffset(0.5, "y");
   */

   TGraph *hMC35 = new TGraph(nlines35, arrayE35, arrayMC35);
   TGraph *hMC4 = new TGraph(nlines4, arrayE4, arrayMC4);

   hMC35->SetLineColor(kRed+2);
   hMC4->SetLineColor(kRed+1);
   hMC35->SetLineWidth(3);
   hMC4->SetLineWidth(3);

   Float_t yfrom = 0;
   Float_t yto = 7;

   Float_t xfrom = 0;
   Float_t xto = 380;

/*
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
*/

   if (!kUseDCCorrection) {
//      hMC2->GetYaxis()->SetRangeUser(yfrom - correction_2, yto - correction_2);
//      hMC3->GetYaxis()->SetRangeUser(yfrom - correction_3, yto - correction_3);
      hMC35->GetYaxis()->SetRangeUser(yfrom - correction_35, yto - correction_35);
      hMC4->GetYaxis()->SetRangeUser(yfrom - correction_4, yto - correction_4);
//      hMC5->GetYaxis()->SetRangeUser(yfrom - correction_5, yto - correction_5);
//      hMC6->GetYaxis()->SetRangeUser(yfrom - correction_6, yto - correction_6);
   }

   Float_t textX = 7.22;
   Float_t textY = 1.7;
   
   hMC35->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC4->GetYaxis()->SetRangeUser(yfrom, yto);
   hMC35->GetXaxis()->SetRangeUser(xfrom, xto);
   hMC4->GetXaxis()->SetRangeUser(xfrom, xto);

   graphPad->cd(1);
   gPad->SetGridy();
   hMC35->SetTitle("");
   hMC35->Draw("LA");
   TText *t2 = new TText();
   t2->SetTextSize(0.15);
   t2->SetTextFont(22);
   t2->DrawText(textX, textY - correction_3 * (!kUseDCCorrection), "230 MeV proton beam");

   graphPad->cd(2);
   gPad->SetGridy();
   hMC4->SetTitle("");
   hMC4->Draw("LA");
   t2->DrawText(textX, textY - correction_2 * (!kUseDCCorrection), "229.25 MeV/u Helium beam");

   TGraph *hMCRS35 = new TGraph(nlines35, arrayE35, arrayRS35);
   TGraph *hMCRS4 = new TGraph(nlines4, arrayE4, arrayRS4);

   hMCRS35->SetLineColor(kRed+2);
   hMCRS4->SetLineColor(kRed+1);
   hMCRS35->SetLineWidth(3);
   hMCRS4->SetLineWidth(3);

   yfrom = 0;
   yto = 20;

   xfrom = 10;
   xto = 310;

   textX = 7.22;
   textY = 8;

   Float_t layerWEPL = 8.12;

   hMCRS35->GetYaxis()->SetRangeUser(yfrom, yto);
   hMCRS4->GetYaxis()->SetRangeUser(yfrom, yto);
   hMCRS35->GetXaxis()->SetRangeUser(xfrom, xto);
   hMCRS4->GetXaxis()->SetRangeUser(xfrom, xto);

   graphPad2->cd(1);
   gPad->SetGridy();
   hMCRS35->SetTitle("");
   hMCRS35->Draw("LA");
   TText *t3 = new TText();
   t3->SetTextSize(0.15);
   t3->SetTextFont(22);
   t3->DrawText(textX, textY - correction_3 * (!kUseDCCorrection), "230 MeV proton beam");
   TLine *l1 = new TLine(xfrom, layerWEPL, xto, layerWEPL);
   l1->SetLineWidth(2);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   TLine *l2 = new TLine(xfrom, layerWEPL*2, xto, layerWEPL*2);
   l2->SetLineWidth(2);
   l2->SetLineStyle(7);
   l2->SetLineColor(kBlack);
   l2->Draw();


   graphPad2->cd(2);
   gPad->SetGridy();
   hMCRS4->SetTitle("");
   hMCRS4->Draw("LA");
   t3->DrawText(textX, textY - correction_2 * (!kUseDCCorrection), "229.25 MeV/u Helium beam");
   l1->Draw();
   l2->Draw();
}
