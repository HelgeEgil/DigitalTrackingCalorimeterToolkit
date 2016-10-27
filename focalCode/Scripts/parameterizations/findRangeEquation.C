#include <vector>
#include <algorithm>

#include <TObject.h>
#include <TSpline.h>
#include <TAxis.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TLine.h>
#include <TSpline.h>

using namespace std;

void Run();
Double_t fitfunc_DBP(Double_t *v, Double_t *par);
Double_t fitfunc_Ulmer(Double_t *v, Double_t *par);

Double_t fitfunc_DBP(Double_t *v, Double_t *par) {

   Float_t depth = v[0];
   Float_t alpha = par[0];
   Float_t p = par[1];
   Float_t range = par[2];

   Double_t fitval = 0;
   fitval = 1 / (p * pow(alpha, 1/p)) * pow(range - depth, 1/p - 1);
   
   if (isnan(fitval)) fitval = 0;

   return fitval;
}

Double_t fitfunc_Ulmer(Double_t *v, Double_t *par) {
   Float_t depth = v[0];
   Float_t range = par[10];
   
   Float_t c1 = par[0];
   Float_t c2 = par[2];
   Float_t c3 = par[4];
   Float_t c4 = par[6];
   Float_t c5 = par[8];
   Float_t l1  = par[1];
   Float_t l2  = par[3];
   Float_t l3  = par[5];
   Float_t l4  = par[7];
   Float_t l5  = par[9];

   Float_t fitval = 0;

   fitval = (1 - l1 * (range - depth)) * c1 * exp(-l1 * (range - depth)) +
            (1 - l2 * (range - depth)) * c2 * exp(-l2 * (range - depth)) +
            (1 - l3 * (range - depth)) * c3 * exp(-l3 * (range - depth)) +
            (1 - l4 * (range - depth)) * c4 * exp(-l4 * (range - depth)) +
            (1 - l5 * (range - depth)) * c5 * exp(-l5 * (range - depth));

    if (depth > range) fitval = 0;

   return fitval;
}

void Run() {
   TCanvas          *c1 = new TCanvas("c1", "Bragg-Kleeman", 1200, 900);
   TCanvas          *c2 = new TCanvas("c2", "Ulmer 1", 1200, 900);
//   TCanvas          *c3 = new TCanvas("c3", "Ulmer 2", 1200, 900);
   TCanvas          *c4 = new TCanvas("c4", "Linear interpolation", 1200, 900);
   TCanvas          *c5 = new TCanvas("c5", "Spline interpolation", 1200, 900);
   TCanvas          *c1Inv = new TCanvas("c1Inv", "Inverse Bragg-Kleeman", 1200, 900);
   TCanvas          *c2Inv = new TCanvas("c2Inv", "Inverse Ulmer 1", 1200, 900);
//   TCanvas          *c3Inv = new TCanvas("c3Inv", "Inverse Ulmer 2", 1200, 900);
   TCanvas          *c4Inv = new TCanvas("c4Inv", "Inverse Linear interpolation", 1200, 900);
   TCanvas          *c5Inv = new TCanvas("c5Inv", "Inverse Spline interpolation", 1200, 900);


   c1->cd();
   TPad *pad11 = new TPad("pad11", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad12 = new TPad("pad12", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad11->Draw();
   pad12->Draw();

   c1Inv->cd();
   TPad *pad11Inv = new TPad("pad11Inv", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad12Inv = new TPad("pad12Inv", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad11Inv->Draw();
   pad12Inv->Draw();

   c2->cd();
   TPad *pad21 = new TPad("pad21", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad22 = new TPad("pad22", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad21->Draw();
   pad22->Draw();

   c2Inv->cd();
   TPad *pad21Inv = new TPad("pad21Inv", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad22Inv = new TPad("pad22Inv", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad21Inv->Draw();
   pad22Inv->Draw();

   c4->cd();
   TPad *pad41 = new TPad("pad41", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad42 = new TPad("pad42", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad41->Draw();
   pad42->Draw();

   c4Inv->cd();
   TPad *pad41Inv = new TPad("pad41Inv", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad42Inv = new TPad("pad42Inv", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad41Inv->Draw();
   pad42Inv->Draw();

   c5->cd();
   TPad *pad51 = new TPad("pad51", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad52 = new TPad("pad52", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad51->Draw();
   pad52->Draw();

   c5Inv->cd();
   TPad *pad51Inv = new TPad("pad51Inv", "The pad 80% of the height", 0.0, 0.3, 1.0, 1.0, 0);
   TPad *pad52Inv = new TPad("pad52Inv", "The pad 20% of the height", 0.0, 0.00, 1.0, 0.3, 0);
   pad51Inv->Draw();
   pad52Inv->Draw();

   Float_t           energy, range, sigma, nuclearfraction, range_csda;
   Float_t           chi2_BK, chi2_Ulmer, chi2_UlmerInv;
   Int_t             absorber;
   Int_t             idx;
   ifstream          in;
   TF1             * braggKleeman = nullptr;
   TF1             * braggKleemanInv = nullptr;
   TF1             * fitFunction = nullptr;
   TF1             * fitInvFunction = nullptr;
   const Int_t       numberOfEnergies = 132;
   Double_t          ranges[numberOfEnergies] = {};
   Double_t          sigmas[numberOfEnergies] = {};
   Double_t          energies[numberOfEnergies] = {}; // {0, 10, 30, 50, 70, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
   Double_t          ranges_control[numberOfEnergies] = {};
   Double_t          energies_control[numberOfEnergies] = {};
   Double_t          al_alpha, al_p, al_a1, al_b1, al_g1, al_b2, al_g2;
   Double_t          al_c1, al_c2, al_c3, al_c4, al_c5;
   Double_t          al_l1, al_l2, al_l3, al_l4, al_l5;

   Double_t           deltaBK[numberOfEnergies] = {};
   Double_t           deltaBKInv[numberOfEnergies] = {};
   Double_t           deltaUlmer[numberOfEnergies] = {};
   Double_t           deltaUlmerInv[numberOfEnergies] = {};
   Double_t           deltaSpline[numberOfEnergies] = {};
   Double_t           deltaSplineInv[numberOfEnergies] = {};
   Double_t           deltaLinear[numberOfEnergies] = {};
   Double_t           deltaLinearInv[numberOfEnergies] = {};

   cout << "Reading file\n";
//   in.open("ranges_water_pstar.csv");
   in.open("Ranges_2mm_Al.csv");
   idx = 0;
   Bool_t useCSDA = true;
   Bool_t useCompleteDataset = true;
   Int_t idxPara = 0;
   Int_t idxCtrl = 0;
   
   Int_t counter = 0;
   while (1) {

      in >> energy >> range_csda >> range;
//      if (counter++%10 != 0) continue; // Reduce N of data points
      if (!in.good()) break;
      
      if (energy > 1200) {
         break;
      }

      if (useCSDA) range = range_csda;

      if (idx%2 == 0 || useCompleteDataset) { // even
         ranges[idxPara] = range;
         energies[idxPara] = energy;
         idxPara++;
      }

      if (idx%2 != 0 || useCompleteDataset) { // odd
         ranges_control[idxCtrl] = range;
         energies_control[idxCtrl] = energy;
         idxCtrl++;
      }
      idx++;
   }

   TGraph *gBK = new TGraph(idxPara, energies, ranges);
   TGraph *gBKInv = new TGraph(idxPara, ranges, energies);
   TGraph *gUlmer = new TGraph(idxPara, energies, ranges);
   TGraph *gUlmerInv = new TGraph(idxPara, ranges, energies);
   TGraph *gSpline = new TGraph(idxPara, energies, ranges);
   TGraph *gSplineInv = new TGraph(idxPara, ranges, energies);
   TGraph *gLinear = new TGraph(idxPara, energies, ranges);
   TGraph *gLinearInv = new TGraph(idxPara, ranges, energies);
   TSpline3 *spline = new TSpline3("spline", energies, ranges, idxPara);
   TSpline3 *splineInv = new TSpline3("splineInv", ranges, energies, idxPara);

   pad11->cd();
   gBK->SetTitle("Bragg-Kleeman fit;Energy [MeV];Range [cm]");
   gPad->SetLogx(); gPad->SetLogy();
   gBK->Draw("A*");

   pad11Inv->cd();
   gBKInv->SetTitle("Inverse Bragg-Kleeman fit;Range [cm];Energy [MeV]");
   gPad->SetLogx(); gPad->SetLogy();
   gBKInv->Draw("A*");

   pad21->cd(); 
   gUlmer->SetTitle("Ulmer fit;Energy [MeV];Range [cm]");
   gPad->SetLogx(); gPad->SetLogy();
   gUlmer->Draw("A*");
   
   pad21Inv->cd();
   gUlmerInv->SetTitle("Inverse Ulmer fit;Range [cm];Energy [cm]");
   gPad->SetLogx(); gPad->SetLogy();
   gUlmerInv->Draw("A*");

   pad41->cd();
   gLinear->SetTitle("Linear interpolation;Energy [MeV];Range [cm]");
   gLinear->SetLineColor(kRed);
   gLinear->SetLineWidth(2);
   gPad->SetLogx(); gPad->SetLogy();
   gLinear->Draw("LA*");

   pad41Inv->cd();
   gLinearInv->SetTitle("Inverse linear intepolation;Range [cm];Energy [MeV]");
   gLinearInv->SetLineColor(kRed);
   gLinearInv->SetLineWidth(2);
   gPad->SetLogx(); gPad->SetLogy();
   gLinearInv->Draw("LA*");

   pad51->cd();
   gSpline->SetTitle("Spline fit;Energy [MeV];Range [cm]");
   gSpline->Draw("*A");
   spline->SetLineColor(kRed);
   spline->SetLineWidth(2);
   gPad->SetLogx(); gPad->SetLogy();
   spline->Draw("same");
   
   pad51Inv->cd();
   gSplineInv->SetTitle("Inverse Spline fit;Range [cm];Energy [MeV]");
   gSplineInv->Draw("*A");
   splineInv->SetLineColor(kRed);
   splineInv->SetLineWidth(2);
   gPad->SetLogx(); gPad->SetLogy();
   splineInv->Draw("same");

   cout << "Fitting functions\n";
   braggKleeman = new TF1("braggKleeman", "[0] * pow(x, [1])", 0, 500);
   braggKleemanInv = new TF1("braggKleemanInv", "pow(x/[0], 1/[1])", 0, 500);
   fitFunction = new TF1("fitFunction", "[0] * x * (1 + ([1] - [1]* exp(-[2] * x)) + ([3] - [3] * exp(-[4]*x)))", 0, 500);
   fitInvFunction = new TF1("fitInvFunction", "x * ([0] * exp ( - [1] * x) + [2] * exp( - [3] * x) + [4] * exp( - [5] * x) + [6] * exp(-[7] * x) + [8] * exp(-[9] * x))", 0, 500);
 

   braggKleeman->SetParameters(0.45, 1.7);
   braggKleemanInv->SetParameters(0.45, 1.7);
   fitFunction->SetParameters(6.94656e-2/2, 15.14450027, 0.001260021, 29.84400076, 0.003260031);
   fitInvFunction->SetParameters(9.663872*2, 1/0.975*2, 2.50472*2, 1/12.4999*2, 0.880745*2, 1/57.001*2, 0.419001*2, 1/106.501*2, 0.92732*2, 1/1067.2784*2);

   fitFunction->SetParLimits(0, 0.005, 0.5);
   fitFunction->SetParLimits(1, 1, 50);
   fitFunction->SetParLimits(2, 0.0001, 0.0045);
   fitFunction->SetParLimits(3, 10, 90);
   fitFunction->SetParLimits(4, 0.0006, 0.01);
   

   fitInvFunction->SetParLimits(0, 0.1, 100);
   fitInvFunction->SetParLimits(1, 0.0001, 1);
   fitInvFunction->SetParLimits(2, 0.1, 100);
   fitInvFunction->SetParLimits(3, 0.0001, 1);
   fitInvFunction->SetParLimits(4, 0.1, 100);
   fitInvFunction->SetParLimits(5, 0.0001, 1);
   fitInvFunction->SetParLimits(6, 0.1, 100);
   fitInvFunction->SetParLimits(7, 0.0001, 1);
   fitInvFunction->SetParLimits(8, 0.1, 100);
   fitInvFunction->SetParLimits(9, 0.0001, 1);

   braggKleeman->SetNpx(1000);
   fitFunction->SetNpx(1000);
   fitInvFunction->SetNpx(1000);


   c1->cd();
   gBK->Fit("braggKleeman", "M, B, Q");
   c1Inv->cd();
   gBKInv->Fit("braggKleemanInv", "M, B, Q");
   c2->cd(); 
   gUlmer->Fit("fitFunction", "M, Q, B");
   c2Inv->cd();
   gUlmerInv->Fit("fitInvFunction", "M, B, Q");


   // calculate the difference between the models and the CONTROL points
   Float_t controlValue;
   for (Int_t i=0; i<idxCtrl; i++) {
      energy = energies_control[i];
      controlValue = ranges_control[i];

      deltaBK[i] = 100*fabs(braggKleeman->Eval(energy) - controlValue) / controlValue;
      deltaBKInv[i] = 100*fabs(braggKleemanInv->Eval(controlValue) - energy) / energy;
      deltaUlmer[i] = 100*fabs(fitFunction->Eval(energy) - controlValue) / controlValue;
      deltaUlmerInv[i] = 100*fabs(fitInvFunction->Eval(controlValue) - energy) / energy;
      deltaSpline[i] = 100*fabs(spline->Eval(energy) - controlValue) / controlValue ;
      deltaSplineInv[i] = 100*fabs(splineInv->Eval(controlValue) - energy) / energy;
      deltaLinear[i] = 100*fabs(gLinear->Eval(energy) - controlValue) / controlValue;
      deltaLinearInv[i] = 100*fabs(gLinearInv->Eval(controlValue) - energy) / energy;
   }

   cout << "THE NUMBER OF DATA POINTS TO OBTAIN THE MODEL IS THUS: " << idxPara << endl;

   TGraph *gBKControl = new TGraph(idxCtrl, energies_control, deltaBK);
   TGraph *gBKInvControl = new TGraph(idxCtrl, ranges_control, deltaBKInv);
   TGraph *gUlmerControl = new TGraph(idxCtrl, energies_control, deltaUlmer);
   TGraph *gUlmerInvControl = new TGraph(idxCtrl, ranges_control, deltaUlmerInv);
   TGraph *gSplineControl = new TGraph(idxCtrl, energies_control, deltaSpline);
   TGraph *gSplineInvControl = new TGraph(idxCtrl, ranges_control, deltaSplineInv);
   TGraph *gLinearControl = new TGraph(idxCtrl, energies_control, deltaLinear);
   TGraph *gLinearInvControl = new TGraph(idxCtrl, ranges_control, deltaLinearInv);

   pad12->cd();
   gBKControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gBKControl->SetLineWidth(2);
   gBKControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gBKControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gBKControl->GetYaxis()->SetNdivisions(5);
   gBKControl->GetYaxis()->SetNoExponent();
   gBKControl->Draw("AL");
//   TLine *line = new TLine(0, 0, 1045, 0); line->Draw();

   pad12Inv->cd(); 
   gBKInvControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gBKInvControl->SetLineWidth(2);
   gBKInvControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gBKInvControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gBKInvControl->GetYaxis()->SetNdivisions(5);
   gBKInvControl->GetYaxis()->SetNoExponent();
   gBKInvControl->Draw("AL");
//   line = new TLine(0, 0, 333, 0); line->Draw();

   pad22->cd(); 
   gUlmerControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gUlmerControl->SetLineWidth(2);
   gUlmerControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gUlmerControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gUlmerControl->GetYaxis()->SetNdivisions(5);
   gUlmerControl->GetYaxis()->SetNoExponent();
   gUlmerControl->Draw("AL");
//   line = new TLine(0, 0, 1045, 0); line->Draw();

   pad22Inv->cd();
   gUlmerInvControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gUlmerInvControl->SetLineWidth(2);
   gUlmerInvControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gUlmerInvControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gUlmerInvControl->GetYaxis()->SetNdivisions(5);
   gUlmerInvControl->GetYaxis()->SetNoExponent();
   gUlmerInvControl->Draw("AL");
//   line = new TLine(0, 0, 333, 0); line->Draw();

   pad52->cd(); 
   gSplineControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gSplineControl->SetLineWidth(2);
   gSplineControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gSplineControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gSplineControl->GetYaxis()->SetNdivisions(5);
   gSplineControl->GetYaxis()->SetNoExponent();
   gSplineControl->Draw("AL");
//   line = new TLine(0, 0, 1045, 0); line->Draw();

   pad52Inv->cd(); 
   gSplineInvControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gSplineInvControl->SetLineWidth(2);
   gSplineInvControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gSplineInvControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gSplineInvControl->GetYaxis()->SetNdivisions(5);
   gSplineInvControl->GetYaxis()->SetNoExponent();
   gSplineInvControl->Draw("AL");
//   line = new TLine(0, 0, 333, 0); line->Draw();

   pad42->cd();
   gLinearControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gLinearControl->SetLineWidth(2);
   gLinearControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gLinearControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gLinearControl->GetYaxis()->SetNdivisions(5);
   gLinearControl->GetYaxis()->SetNoExponent();
   gLinearControl->Draw("AL");
//   line = new TLine(0, 0, 1045, 0); line->Draw();

   pad42Inv->cd();  
   gLinearInvControl->SetTitle("Deviation from PSTAR;Energy [MeV];Range deviation [%]"); 
   gLinearInvControl->SetLineWidth(2);
   gLinearInvControl->SetLineColor(kRed);
   gPad->SetLogx(); gPad->SetLogy();
   gLinearInvControl->GetYaxis()->SetRangeUser(1e-5, 500);
   gLinearInvControl->GetYaxis()->SetNdivisions(5);
   gLinearInvControl->GetYaxis()->SetNoExponent();
   gLinearInvControl->Draw("AL");
//   line = new TLine(0, 0, 333, 0); line->Draw();


   // CALCULATE CHI SQUARE FOR THE DIFFERENT MODELS
   cout << "Calculating accuracies\n";

   chi2_BK = gBK->Chisquare(braggKleeman);
   Float_t chi2_BKInv = gBKInv->Chisquare(braggKleemanInv);
   chi2_Ulmer = gUlmer->Chisquare(fitFunction);
   chi2_UlmerInv = gUlmerInv->Chisquare(fitInvFunction);

   Float_t rms_BK = 0, rms_BKInv = 0, rms_Ulmer = 0, rms_UlmerInv = 0;
   Float_t rms_Spline = 0, rms_SplineInv = 0, rms_Linear = 0, rms_LinearInv = 0;

   for (Int_t i=0; i<idxCtrl; i++) {
      energy = energies_control[i];
      controlValue = ranges_control[i];

      rms_BK += pow(deltaBK[i]/100*controlValue, 2);
      rms_BKInv += pow(deltaBKInv[i]/100*energy, 2);
      rms_Ulmer += pow(deltaUlmer[i]/100*controlValue, 2);
      rms_UlmerInv += pow(deltaUlmerInv[i]/100*energy, 2);
      rms_Spline += pow(deltaSpline[i]/100*controlValue, 2);
      rms_SplineInv += pow(deltaSplineInv[i]/100*energy, 2);
      rms_Linear += pow(deltaLinear[i]/100*controlValue, 2);
      rms_LinearInv += pow(deltaLinearInv[i]/100*energy, 2);
   }

   rms_BK /= idxCtrl;
   rms_BKInv /= idxCtrl;
   rms_Ulmer /= idxCtrl;
   rms_UlmerInv /= idxCtrl;
   rms_Spline /= idxCtrl;
   rms_SplineInv /= idxCtrl;
   rms_Linear /= idxCtrl;
   rms_LinearInv /= idxCtrl;

   rms_BK = sqrt(rms_BK);
   rms_BKInv = sqrt(rms_BKInv);
   rms_Ulmer = sqrt(rms_Ulmer);
   rms_UlmerInv = sqrt(rms_UlmerInv);
   rms_Spline = sqrt(rms_Spline);
   rms_SplineInv = sqrt(rms_SplineInv);
   rms_Linear = sqrt(rms_Linear);
   rms_LinearInv = sqrt(rms_LinearInv);

   cout << "ACCURACY OF  DIFFERENT MODELS:::::::\n";
   cout << "BRAGG-KLEEMAN MODEL: " << endl;
   cout << "CHI SQUARE = " << chi2_BK << endl;
   cout << "RMS = " << rms_BK << endl;
   cout << "INVERSE BRAGG-KLEEMAN MODEL: " << endl;
   cout << "CHI SQUARE = " << chi2_BKInv << endl;
   cout << "RMS = " << rms_BKInv << endl;
   
   cout << endl;
   cout << "ULMER MODEL: " << endl;
   cout << "CHI SQUARE = " << chi2_Ulmer << endl;
   cout << "RMS = " << rms_Ulmer << endl;
   cout << "INVERSE ULMER MODEL: " << endl;
   cout << "CHI SQUARE = " << chi2_UlmerInv << endl;
   cout << "RMS = " << rms_UlmerInv << endl;

   cout << "SPLINE MODEL: " << endl;
   cout << "RMS = " << rms_Spline << endl;
   cout << endl;
   cout << "INVERSE SPLINE MODEL: " << endl;
   cout << "RMS = " << rms_SplineInv << endl;
   cout << endl;
   cout << "LINEAR MODEL: " << endl;
   cout << "RMS = " << rms_Linear << endl;
   cout << endl;
   cout << "INVERSE LINEAR MODEL : " << endl;
   cout << "RMS = " << rms_LinearInv << endl;
   cout << endl;

   TCanvas *cCompare = new TCanvas("cCompare", "Error comparison", 800, 800);
   cCompare->Divide(1,2,0.01,0.01);

   TGraph *gBKCompare = new TGraph(idxCtrl, energies_control, deltaBK);
   TGraph *gBKInvCompare = new TGraph(idxCtrl, ranges_control, deltaBKInv);
   TGraph *gUlmerCompare = new TGraph(idxCtrl, energies_control, deltaUlmer);
   TGraph *gUlmerInvCompare = new TGraph(idxCtrl, ranges_control, deltaUlmerInv);
   TGraph *gSplineCompare = new TGraph(idxCtrl, energies_control, deltaSpline);
   TGraph *gSplineInvCompare = new TGraph(idxCtrl, ranges_control, deltaSplineInv);
   TGraph *gLinearCompare = new TGraph(idxCtrl, energies_control, deltaLinear);
   TGraph *gLinearInvCompare = new TGraph(idxCtrl, ranges_control, deltaLinearInv);

   gBKCompare->SetLineColor(kRed);
   gBKInvCompare->SetLineColor(kRed);
   gUlmerCompare->SetLineColor(kBlue);
   gUlmerInvCompare->SetLineColor(kBlue);
   gSplineCompare->SetLineColor(kGreen);
   gSplineInvCompare->SetLineColor(kGreen);
   gLinearCompare->SetLineColor(kBlack);
   gLinearInvCompare->SetLineColor(kBlack);
   gBKCompare->SetLineWidth(3);
   gBKInvCompare->SetLineWidth(3);
   gUlmerCompare->SetLineWidth(3);
   gUlmerInvCompare->SetLineWidth(3);
   gSplineCompare->SetLineWidth(3);
   gSplineInvCompare->SetLineWidth(3);
   gLinearCompare->SetLineWidth(3);
   gLinearInvCompare->SetLineWidth(3);

   gBKCompare->SetTitle("Range-from-energy accuracy comparison;Energy [MeV];Range error [%]");
   gBKInvCompare->SetTitle("Energy-from-range accuracy comparison;Range [MeV];Energy error [%]");

   cCompare->cd(1);
   gPad->SetLogy(); gPad->SetLogx();
   gBKCompare->Draw("LA");
   gUlmerCompare->Draw("L");
   gSplineCompare->Draw("L");
   gLinearCompare->Draw("L)");
   TLegend *leg1 = new TLegend(0.15, 0.6, 0.35, 0.88);
   leg1->AddEntry(gBKCompare, "Bragg-Kleeman", "L");
   leg1->AddEntry(gUlmerCompare, "Ulmer", "L");
   leg1->AddEntry(gLinearCompare, "Linear interpolation", "L");
   leg1->AddEntry(gSplineCompare, "Spline interpolation", "L");
   leg1->Draw();
   gBKCompare->GetYaxis()->SetRangeUser(1e-5, 1000);

   cCompare->cd(2);
   gPad->SetLogy(); gPad->SetLogx();
   gBKInvCompare->Draw("LA");
   gUlmerInvCompare->Draw("L");
   gSplineInvCompare->Draw("L");
   gLinearInvCompare->Draw("L)");
   TLegend *leg2 = new TLegend(0.15, 0.6, 0.35, 0.88);
   leg2->AddEntry(gBKInvCompare, "Bragg-Kleeman", "L");
   leg2->AddEntry(gUlmerInvCompare, "Ulmer", "L");
   leg2->AddEntry(gLinearInvCompare, "Linear interpolation", "L");
   leg2->AddEntry(gSplineInvCompare, "Spline interpolation", "L");
   leg2->Draw();
   gBKInvCompare->GetYaxis()->SetRangeUser(1e-5, 1000);

   cout << "Finding depth dose curves\n";
   // ASSIGN THE PARAMETERS
   Double_t alpha, p, a1, b1, g1, b2, g2;
   Double_t cc1, cc2, cc3, cc4, cc5;
   Double_t l1, l2, l3, l4, l5;
   Double_t E = 190;

   range = gLinear->Eval(E);

   alpha = braggKleeman->GetParameter(0);
   p     = braggKleeman->GetParameter(1);
   a1    = fitFunction->GetParameter(0);
   b1    = fitFunction->GetParameter(1);
   g1    = fitFunction->GetParameter(2);
   b2    = fitFunction->GetParameter(3);
   g2    = fitFunction->GetParameter(4);
   cc1    = fitInvFunction->GetParameter(0);
   cc2    = fitInvFunction->GetParameter(2);
   cc3    = fitInvFunction->GetParameter(4);
   cc4    = fitInvFunction->GetParameter(6);
   cc5    = fitInvFunction->GetParameter(8);
   l1    = fitInvFunction->GetParameter(1);
   l2    = fitInvFunction->GetParameter(3);
   l3    = fitInvFunction->GetParameter(5);
   l4    = fitInvFunction->GetParameter(7);
   l5    = fitInvFunction->GetParameter(9);

   cout << "ALPHA = " << alpha << ", p = " << p << endl;

   TCanvas *cDepth = new TCanvas("cDepth", "depth-dose curves", 800, 800);
   cDepth->Divide(1,5,0.001,0.001);
   
   // the depth of the semi-empirical models is easy to find, since we have formulas for it
   // make depth dose functions
   
   TF1 *DDBK = new TF1("DDBK", fitfunc_DBP, 0, range+3, 3);
   DDBK->SetParameters(alpha, p, range);
   DDBK->SetLineColor(kBlue);
  
   TF1 *DDUlmer = new TF1("DDUlmer", fitfunc_Ulmer, 0, range+3, 11); 
   DDUlmer->SetParameters(cc1, l1, cc2, l2, cc3, l3, cc4, l4, cc5, l5, range);

   // The depth dose of the interpolation models is a bit harder to find...
   // For a given R0, what is the R0-z at a given z?
   // Find the dE/dx at this point
   // dE = change in energy (delta y at dx)
   // dx = change in range (constant)
   
   Float_t z = 0;
   idx = 0;
   Float_t  thisZ = 0;
   Float_t  dx = 0.001; // 0.01 mm
   Float_t  E1, E2, dE, dEdx;
   const Int_t size = 190 / 0.001 + 100;
   Float_t  dEdx_Linear[size] = {};
   Float_t  dEdx_Spline[size] = {};
   Float_t  x_Linear[size] = {};
   Float_t  x_Spline[size] = {};

   while (1) {
      thisZ = range - z;
      E1 = gLinearInv->Eval(thisZ + dx/2);
      E2 = gLinearInv->Eval(thisZ - dx/2);
      dE = E1 - E2;
      dEdx = dE / dx;

      dEdx_Linear[idx] = dEdx;
      x_Linear[idx] = z;

      // STEP
      z += dx;
      idx++;

      if (z > range) {
         dEdx_Linear[idx] = 0;
         x_Linear[idx] = z;
         break;
      }
   }

   TGraph *DDLinear = new TGraph(idx+1, x_Linear, dEdx_Linear);

   z = 0;
   idx = 0;
   while (1) {
      thisZ = range - z;
      E1 = gSplineInv->Eval(thisZ + dx/2);
      E2 = gSplineInv->Eval(thisZ - dx/2);
      dE = E1 - E2;
      dEdx = dE / dx;

      dEdx_Spline[idx] = dEdx;
      x_Spline[idx] = z;

      // STEP
      z += dx;
      idx++;

      if (z > range) {
         dEdx_Spline[idx] = 0;
         x_Spline[idx] = z;
         break;
      }
   }
   TGraph *DDSpline = new TGraph(idx+1, x_Spline, dEdx_Spline);


   cDepth->cd(1);
   DDBK->Draw();
   TLine *r1 = new TLine(range, 0, range, 60); r1->Draw();
   DDBK->GetHistogram()->GetYaxis()->SetRangeUser(0, 60);
   DDBK->GetHistogram()->SetTitle("Depth dose curve from Bragg-Kleeman;Depth [cm];dE/dx [A.U.]");

   cDepth->cd(2);
   DDUlmer->Draw();
   TLine *r2 = new TLine(range, 0, range, 60); r2->Draw();
   DDUlmer->GetHistogram()->GetYaxis()->SetRangeUser(0, 60);
   DDUlmer->GetHistogram()->GetXaxis()->SetRangeUser(0, range+3);
   DDUlmer->GetHistogram()->SetTitle("Depth dose curve from Ulmer;Depth [cm];dE/dx [A.U.]");

   cDepth->cd(3);
   DDLinear->Draw("LA");
   DDLinear->GetYaxis()->SetRangeUser(0, 60);
   DDLinear->GetXaxis()->SetLimits(0, range+3);
   DDLinear->SetTitle("Depth dose curve from linear interpolation;Depth [cm];dE/dx [A.U.]");

   cDepth->cd(4);
   DDSpline->Draw("LA");
   DDSpline->GetYaxis()->SetRangeUser(0, 60);
   DDSpline->GetXaxis()->SetLimits(0, range+3);
   DDSpline->SetTitle("Depth dose curve from Spline interpolation;Depth [cm];dE/dx [A.U.]");
   
   cDepth->Update();

   cDepth->cd(5);
   DDBK->Draw();
   DDUlmer->Draw("same");
   DDLinear->Draw("L");
   DDSpline->Draw("L");


   /*
   TLegend *leg3 = new TLegend(0.18, 0.75, 0.46, 0.87);
   leg3->AddEntry(DDBK, "Bragg-Kleeman", "L");
   leg3->AddEntry(DDUlmer, "Ulmer", "L");
   leg3->Draw();
   */
}

