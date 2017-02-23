#include <vector>
#include <algorithm>

#include <TObject.h>
#include <TSpline.h>
#include <TAxis.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>

using namespace std;

void getPValues() {
   TCanvas *c1 = new TCanvas("c1", "Comparison of depth dose functions", 1200, 900);
   vector<TCanvas*>  canvasVector;
   vector<TCanvas*>  canvasInvVector;
   vector<TGraph*>   graphVector;
   vector<TGraph*>   graphInvVector;
   Float_t           energy, range, sigma, nuclearfraction;
   Float_t           chi2_BK, chi2_Ulmer, chi2_UlmerInv;
   Int_t             absorber;
   Int_t             idx;
   ifstream          in;
   TF1             * braggKleeman = nullptr;
   TF1             * fitFunction = nullptr;
   TF1             * fitInvFunction = nullptr;
   Double_t           range_1mm[17] = {};
   Double_t           range_2mm[17] = {};
   Double_t           range_3mm[17] = {};
   Double_t           range_4mm[17] = {};
   Double_t           range_5mm[17] = {};
   Double_t           range_6mm[17] = {};
   Double_t           range_7mm[17] = {};
   Double_t           range_8mm[17] = {};
   Double_t           range_9mm[17] = {};
   Double_t           range_10mm[17] = {};
   
   Double_t           sigma_1mm[17] = {};
   Double_t           sigma_2mm[17] = {};
   Double_t           sigma_3mm[17] = {};
   Double_t           sigma_4mm[17] = {};
   Double_t           sigma_5mm[17] = {};
   Double_t           sigma_6mm[17] = {};
   Double_t           sigma_7mm[17] = {};
   Double_t           sigma_8mm[17] = {};
   Double_t           sigma_9mm[17] = {};
   Double_t           sigma_10mm[17] = {};
   Double_t           energies[17] = {0, 10, 30, 50, 70, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
   Double_t           energyError[17] = {};
      
   Double_t al_alpha, al_p, al_a1, al_b1, al_g1, al_b2, al_g2;
   Double_t al_c1, al_c2, al_c3, al_c4, al_c5;
   Double_t al_l1, al_l2, al_l3, al_l4, al_l5;

   for (Int_t i=1; i<=10; i++) {
      canvasVector.push_back(new TCanvas(Form("canvas_%d_mm_absorber", i), Form("%d mm absorber", i), 1200, 900));
      canvasInvVector.push_back(new TCanvas(Form("invCanvas_%d_mm_absorber", i), Form("Inverse %d mm absorber", i), 1200, 900));
   }
   
   in.open("../OutputFiles/findManyRanges.csv");
   while (1) {
      in >> energy >> absorber >> range >> sigma >> nuclearfraction;
      if (!in.good()) break;

      if (energy<90) idx = energy/20+1;
      else           idx = energy/10-4;

      switch (absorber) {
         case 1:
            range_1mm[idx] = range;
            sigma_1mm[idx] = sigma;
            break;
         case 2:
            range_2mm[idx] = range;
            sigma_2mm[idx] = sigma;
            break;
         case 3:
            range_3mm[idx] = range;
            sigma_3mm[idx] = sigma;
            break;
         case 4:
            range_4mm[idx] = range;
            sigma_4mm[idx] = sigma;
            break;
         case 5:
            range_5mm[idx] = range;
            sigma_5mm[idx] = sigma;
            break;
         case 6:
            range_6mm[idx] = range;
            sigma_6mm[idx] = sigma;
            break;
         case 7:
            range_7mm[idx] = range;
            sigma_7mm[idx] = sigma;
            break;
         case 8:
            range_8mm[idx] = range;
            sigma_8mm[idx] = sigma;
            break;
         case 9:
            range_9mm[idx] = range;
            sigma_9mm[idx] = sigma;
            break;
         case 10:
            range_10mm[idx]= range;
            sigma_10mm[idx]= sigma;
            break;
      }
   }

   /*
   graphVector.push_back(new TGraphErrors(16, energies, range_1mm, energyError, sigma_1mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_2mm, energyError, sigma_2mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_3mm, energyError, sigma_3mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_4mm, energyError, sigma_4mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_5mm, energyError, sigma_5mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_6mm, energyError, sigma_6mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_7mm, energyError, sigma_7mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_8mm, energyError, sigma_8mm));
   graphVector.push_back(new TGraphErrors(17, energies, range_9mm, energyError, sigma_9mm));
   graphVector.push_back(new TGraphErrors(13, energies+4, range_10mm+4, energyError+4, sigma_10mm+4));

   graphInvVector.push_back(new TGraphErrors(16, range_1mm, energies, sigma_1mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_2mm, energies, sigma_2mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_3mm, energies, sigma_3mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_4mm, energies, sigma_4mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_5mm, energies, sigma_5mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_6mm, energies, sigma_6mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_7mm, energies, sigma_7mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_8mm, energies, sigma_8mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_9mm, energies, sigma_9mm, energyError));
   graphInvVector.push_back(new TGraphErrors(17, range_10mm, energies, sigma_10mm, energyError));

   */

   graphVector.push_back(new TGraph(16, energies, range_1mm));
   graphVector.push_back(new TGraph(17, energies, range_2mm));
   graphVector.push_back(new TGraph(17, energies, range_3mm));
   graphVector.push_back(new TGraph(17, energies, range_4mm));
   graphVector.push_back(new TGraph(17, energies, range_5mm));
   graphVector.push_back(new TGraph(17, energies, range_6mm));
   graphVector.push_back(new TGraph(17, energies, range_7mm));
   graphVector.push_back(new TGraph(17, energies, range_8mm));
   graphVector.push_back(new TGraph(17, energies, range_9mm));
   graphVector.push_back(new TGraph(17, energies, range_10mm));

   graphInvVector.push_back(new TGraph(16, range_1mm, energies));
   graphInvVector.push_back(new TGraph(17, range_2mm, energies));
   graphInvVector.push_back(new TGraph(17, range_3mm, energies));
   graphInvVector.push_back(new TGraph(17, range_4mm, energies));
   graphInvVector.push_back(new TGraph(17, range_5mm, energies));
   graphInvVector.push_back(new TGraph(17, range_6mm, energies));
   graphInvVector.push_back(new TGraph(17, range_7mm, energies));
   graphInvVector.push_back(new TGraph(17, range_8mm, energies));
   graphInvVector.push_back(new TGraph(17, range_9mm, energies));
   graphInvVector.push_back(new TGraph(17, range_10mm, energies));


   for (Int_t i=0; i<10; i++) {
      canvasVector.at(i)->cd();
      graphVector.at(i)->SetTitle(Form("Range fit for DTC using %d mm absorbator", i+1));
      graphVector.at(i)->GetXaxis()->SetRangeUser(0, 220);
      graphVector.at(i)->Draw("A*");

      canvasInvVector.at(i)->cd();
      graphInvVector.at(i)->SetTitle(Form("Energy fit for DTC using %d mm absorbator", i+1));
      graphInvVector.at(i)->GetXaxis()->SetRangeUser(0, 220);
      graphInvVector.at(i)->Draw("A*");

      braggKleeman = new TF1("braggKleeman", "[0] * pow(x, [1])", 10, 250);
      
      fitFunction = new TF1("fitFunction", "[0] * x * (1 + ([1] - [1]* exp(-[2] * x)) + ([3] - [3] * exp(-[4]*x)))", 10, 250);
      fitInvFunction = new TF1("fitInvFunction", "x * ([0] * exp ( - [1] * x) + [2] * exp( - [3] * x) + [4] * exp( - [5] * x) + [6] * exp(-[7] * x) + [8] * exp(-[9] * x))", 1, 400);
    
      braggKleeman->SetParameters(0.45, 1.7);
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
   
      canvasVector.at(i)->cd();
      graphVector.at(i)->Fit("fitFunction", "M, Q, B");
      graphVector.at(i)->Fit("braggKleeman", "M, B");
      canvasInvVector.at(i)->cd();
      graphInvVector.at(i)->Fit("fitInvFunction", "M, B");

      Double_t alpha, p, a1, b1, g1, b2, g2;
      Double_t c1, c2, c3, c4, c5;
      Double_t l1, l2, l3, l4, l5;

      alpha = braggKleeman->GetParameter(0);
      p     = braggKleeman->GetParameter(1);
      a1    = fitFunction->GetParameter(0);
      b1    = fitFunction->GetParameter(1);
      g1    = fitFunction->GetParameter(2);
      b2    = fitFunction->GetParameter(3);
      g2    = fitFunction->GetParameter(4);
      c1    = fitInvFunction->GetParameter(0);
      c2    = fitInvFunction->GetParameter(2);
      c3    = fitInvFunction->GetParameter(4);
      c4    = fitInvFunction->GetParameter(6);
      c5    = fitInvFunction->GetParameter(8);
      l1    = fitInvFunction->GetParameter(1);
      l2    = fitInvFunction->GetParameter(3);
      l3    = fitInvFunction->GetParameter(5);
      l4    = fitInvFunction->GetParameter(7);
      l5    = fitInvFunction->GetParameter(9);

      if (i==9) {
         al_alpha = alpha; al_p = p;
         al_a1 = a1; al_b1 = b1; al_g1 = g1; al_b2 = b2; al_g2 = g2;
         al_c1 = c1; al_c2 = c2; al_c3 = c3; al_c4 = c4; al_c5 = c5;
         al_l1 = l1; al_l2 = l2; al_l3 = l3; al_l4 = l4; al_l5 = l5;
      }

      ofstream file("../OutputFiles/RangeEnergyParameters.csv", ofstream::out | ofstream::app);
      file << i+1 << " " << alpha << " " << p << " ";
      file << a1 << " " << b1 << " " << g1 << " " << b2 << " " << g2 << " ";
      file << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c5 << " ";
      file << l1 << " " << l2 << " " << l3 << " " << l4 << " " << l5 << endl;
      file.close();

      // CALCULATE CHI SQUARE FOR THE DIFFERENT MODELS

      chi2_BK = graphVector.at(i)->Chisquare(braggKleeman);
      chi2_Ulmer = graphVector.at(i)->Chisquare(fitFunction);
      chi2_UlmerInv = graphInvVector.at(i)->Chisquare(fitInvFunction);

      cout << " --- \033[1m VALUES FOR " << i+1 << " mm ABSORBER \033[0m --- " << endl;
      cout << "BRAGG-KLEEMAN MODEL: " << endl;
      cout << "CHI SQUARE = " << chi2_BK << endl;
      cout << "alpha   = " << braggKleeman->GetParameter(0) << endl;
      cout << "p       = " << braggKleeman->GetParameter(1) << endl;
      cout << endl;
      cout << "ULMER MODEL: " << endl;
      cout << "CHI SQUARE = " << chi2_Ulmer << endl;
      cout << "a1      = " << fitFunction->GetParameter(0) << endl;
      cout << "b1      = " << fitFunction->GetParameter(1) << endl;
      cout << "g1      = " << fitFunction->GetParameter(2) << endl;
      cout << "b2      = " << fitFunction->GetParameter(3) << endl;
      cout << "g2      = " << fitFunction->GetParameter(4) << endl;
      cout << endl;
      cout << "ULMER INVESE MODEL: " << endl;
      cout << "CHI SQUARE = " << chi2_UlmerInv << endl;
      cout << "c1      = " << fitInvFunction->GetParameter(0) << endl;
      cout << "lambda1 = " << fitInvFunction->GetParameter(1) << endl;
      cout << "c2      = " << fitInvFunction->GetParameter(2) << endl;
      cout << "lambda2 = " << fitInvFunction->GetParameter(3) << endl;
      cout << "c3      = " << fitInvFunction->GetParameter(4) << endl;
      cout << "lambda3 = " << fitInvFunction->GetParameter(5) << endl;
      cout << "c4      = " << fitInvFunction->GetParameter(6) << endl;
      cout << "lambda4 = " << fitInvFunction->GetParameter(7) << endl;
      cout << "c5      = " << fitInvFunction->GetParameter(8) << endl;
      cout << "lambda5 = " << fitInvFunction->GetParameter(9) << endl;

      cout << endl << endl;      
      delete fitFunction;
      delete fitInvFunction;
   }

   // find SPLINE from the pure AL case
   canvasVector.at(9)->cd();
   TSpline3 *sp = new TSpline3("Cubic spline", energies, range_10mm, 17);
   sp->SetLineColor(kBlue);
   sp->Draw("lsame");

   c1->cd();
   Double_t E = 0, lastE = 0;
   Double_t energiesSpline[2000];
   Double_t rangeSpline[2000] = {};
   Double_t spline_derivative[2000];
   Double_t ulmer[2000];
   Double_t bortfeld[2000];
   Double_t dE = 0, dx = 0;
   Double_t dxUlmer;
   Double_t ELUlmer;
   Double_t EL = 0;

   for (Int_t j=1; j<2000; j++) {
      E = 200 - j/10.;
      lastE = 200 - (j-1)/10.;

      // MeV / cm
      dx = sp->Eval(lastE) - sp->Eval(E);
      dxUlmer = fitFunction->Eval(lastE) - fitFunction->Eval(E);
      dE = lastE-E;

      EL = dE / dx;
      ELUlmer = dE / dxUlmer;

      spline_derivative[j-1] = EL;
      ulmer[j-1] = ELUlmer;
      rangeSpline[j] = rangeSpline[j-1] + dx;
   }

   TGraph *splineGraph = new TGraph(1999, rangeSpline, spline_derivative);
   splineGraph->SetLineColor(kRed);
   splineGraph->Draw("LA");
   TGraph *bortfeldGraph = new TGraph(1999, rangeSpline, ulmer);
   bortfeldGraph->SetLineColor(kBlue);
   bortfeldGraph->Draw("L");

         
}
