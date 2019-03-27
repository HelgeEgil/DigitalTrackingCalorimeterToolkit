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
Bool_t kUseDCCorrection = false;
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

   Double_t energies[arraySize] = {0};
   Double_t thicknesses[arraySize] = {0};

//   Float_t correction_2 = -0.132, correction_3 = 0.271, correction_35 = 0.813, correction_4 = 0.734, correction_5 = 1.292, correction_6 = 1.879; // 250 MeV
   Float_t correction_2 = 0.4, correction_3 = 0.98, correction_35 = 0.813, correction_4 = 1.37, correction_5 = 1.97 , correction_6 = 2.54;
   if (!kUseDCCorrection) {
      correction_2 = 0;
      correction_3 = 0;
      correction_35 = 0;
      correction_4 = 0;
      correction_5 = 0;
      correction_6 = 0;
   }

   gStyle->SetOptStat(0);

   Float_t nomrange_, estrange_, sigmaRange_, lastRange_, nomsigma_, waterphantomthickness_, dummy0;
   Int_t energy_, thickness_;
   Float_t estimatedStraggling;

   ifstream in1;
   in1.open("../../Data/Ranges/EnergyAfterDegrader230MeV.csv");
   Int_t thick, n=0;
   Double_t energy;
   while (1) {
      in1 >> thick >> energy;
      if (!in1.good()) break;   
      thicknesses[n] = thick;
      energies[n++] = energy;
   }
   in1.close();
   printf("Found %d lines in EnergyAfterDegrader230MeV.csv\n", n);

   TSpline3 *energySpline = new TSpline3("energySpline", thicknesses, energies, n);

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

      energy = energySpline->Eval(energy_);

      if (thickness_ == 2) arrayE2[nlines2] = nomrange_;
      if (thickness_ == 3) arrayE3[nlines3] = nomrange_;
      if (thickness_ == 35) arrayE35[nlines35] = nomrange_;
      if (thickness_ == 4) arrayE4[nlines4] = nomrange_;
      if (thickness_ == 5) arrayE5[nlines5] = nomrange_;
      if (thickness_ == 6) arrayE6[nlines6] = nomrange_;

      if (thickness_ == 2) arrayMC2[nlines2++] = -nomrange_ + estrange_ + correction_2;
      if (thickness_ == 3) arrayMC3[nlines3++] = estrange_ - nomrange_ + correction_3;
      if (thickness_ == 35) arrayMC35[nlines35++] = estrange_ - nomrange_ + correction_35;
      if (thickness_ == 4) arrayMC4[nlines4++] = -nomrange_ + estrange_ + correction_4;
      if (thickness_ == 5) arrayMC5[nlines5++] = estrange_ - nomrange_ + correction_5;
      if (thickness_ == 6) arrayMC6[nlines6++] = -nomrange_ + estrange_ + correction_6;
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

   Float_t xfrom = 5;
   Float_t xto = 370;

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

   Float_t textX = 7.22;
   Float_t textY = 3.11;

   graphPad->cd(1);
   gPad->SetGridy();
   hMC2->Draw("LA");
   TText *t2 = new TText();
   t2->SetTextSize(0.2);
   t2->SetTextFont(22);
   t2->DrawText(textX, textY, "2 mm Al absorber"); 

   graphPad->cd(2);
   gPad->SetGridy();
   hMC3->Draw("LA");
   t2->DrawText(textX, textY, "3 mm Al absorber");
   /*
   graphPad->cd(3);
   gPad->SetGridy();
   hMC35->Draw("LA");
   t2->DrawText(textX, textY, "3.5 mm Al absorber");
*/
   graphPad->cd(3);
   gPad->SetGridy();
   hMC4->Draw("LA");
   t2->DrawText(textX, textY, "4 mm Al absorber");
   
   graphPad->cd(4);
   gPad->SetGridy();
   hMC5->Draw("LA");
   t2->DrawText(textX, textY, "5 mm Al absorber");
   
   graphPad->cd(5);
   gPad->SetGridy();
   hMC6->Draw("LA");
   t2->DrawText(textX, textY, "6 mm Al absorber");


}
