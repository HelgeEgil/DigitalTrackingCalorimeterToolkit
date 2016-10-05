#include <vector>
#include <algorithm>

#include <TObject.h>
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
   vector<TCanvas*>  canvasVector;
   vector<TCanvas*>  canvasInvVector;
   vector<TGraphErrors*>   graphVector;
   vector<TGraphErrors*>   graphInvVector;
   Float_t           energy, range, sigma, nuclearfraction;
   Float_t           chi2_BK, chi2_Ulmer, chi2_UlmerInv;
   Int_t             absorber;
   Int_t             idx;
   ifstream          in;
   TF1             * braggKleeman = nullptr;
   TF1             * fitFunction = nullptr;
   TF1             * fitInvFunction = nullptr;
   Float_t           range_1mm[13] = {};
   Float_t           range_2mm[13] = {};
   Float_t           range_3mm[13] = {};
   Float_t           range_4mm[13] = {};
   Float_t           range_5mm[13] = {};
   Float_t           range_6mm[13] = {};
   Float_t           range_7mm[13] = {};
   Float_t           range_8mm[13] = {};
   Float_t           range_9mm[13] = {};
   Float_t           range_10mm[13] = {};
   
   Float_t           sigma_1mm[13] = {};
   Float_t           sigma_2mm[13] = {};
   Float_t           sigma_3mm[13] = {};
   Float_t           sigma_4mm[13] = {};
   Float_t           sigma_5mm[13] = {};
   Float_t           sigma_6mm[13] = {};
   Float_t           sigma_7mm[13] = {};
   Float_t           sigma_8mm[13] = {};
   Float_t           sigma_9mm[13] = {};
   Float_t           sigma_10mm[13] = {};
   Float_t           energies[13] = {0, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
   Float_t           energyError[13] = {};

   for (Int_t i=1; i<=10; i++) {
      canvasVector.push_back(new TCanvas(Form("canvas_%d_mm_absorber", i), Form("%d mm absorber", i), 1200, 900));
      canvasInvVector.push_back(new TCanvas(Form("invCanvas_%d_mm_absorber", i), Form("Inverse %d mm absorber", i), 1200, 900));
   }
   
   in.open("../OutputFiles/findManyRanges.csv");
   while (1) {
      in >> energy >> absorber >> range >> sigma >> nuclearfraction;
      if (!in.good()) break;

      idx = energy/10-9 + 1;

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

   graphVector.push_back(new TGraphErrors(12, energies, range_1mm, energyError, sigma_1mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_2mm, energyError, sigma_2mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_3mm, energyError, sigma_3mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_4mm, energyError, sigma_4mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_5mm, energyError, sigma_5mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_6mm, energyError, sigma_6mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_7mm, energyError, sigma_7mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_8mm, energyError, sigma_8mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_9mm, energyError, sigma_9mm));
   graphVector.push_back(new TGraphErrors(13, energies, range_10mm, energyError, sigma_10mm));
   
   graphInvVector.push_back(new TGraphErrors(12, range_1mm, energies, sigma_1mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_2mm, energies, sigma_2mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_3mm, energies, sigma_3mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_4mm, energies, sigma_4mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_5mm, energies, sigma_5mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_6mm, energies, sigma_6mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_7mm, energies, sigma_7mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_8mm, energies, sigma_8mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_9mm, energies, sigma_9mm, energyError));
   graphInvVector.push_back(new TGraphErrors(13, range_10mm, energies, sigma_10mm, energyError));

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

      fitFunction->SetParLimits(0, 0.01, 0.1);
      fitFunction->SetParLimits(1, 5, 45);
      fitFunction->SetParLimits(2, 0.0005, 0.0045);
      fitFunction->SetParLimits(3, 10, 90);
      fitFunction->SetParLimits(4, 0.001, 0.01);
      
      fitInvFunction->SetParLimits(0, 0.1, 100);
      fitInvFunction->SetParLimits(1, 0.001, 5);
      fitInvFunction->SetParLimits(2, 0.1, 100);
      fitInvFunction->SetParLimits(3, 0.001, 5);
      fitInvFunction->SetParLimits(4, 0.1, 100);
      fitInvFunction->SetParLimits(5, 0.001, 5);
      fitInvFunction->SetParLimits(6, 0.1, 100);
      fitInvFunction->SetParLimits(7, 0.001, 5);
      fitInvFunction->SetParLimits(8, 0.1, 100);
      fitInvFunction->SetParLimits(9, 0.001, 5);

      braggKleeman->SetNpx(500);
      fitFunction->SetNpx(500);
      fitInvFunction->SetNpx(500);
   
      canvasVector.at(i)->cd();
      graphVector.at(i)->Fit("fitFunction", "M, Q, B");
      graphVector.at(i)->Fit("braggKleeman", "M, Q, B");
      canvasInvVector.at(i)->cd();
      graphInvVector.at(i)->Fit("fitInvFunction", "M, Q, B");

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
      cout << "lambda1 = " << 1/fitInvFunction->GetParameter(1) << endl;
      cout << "c2      = " << fitInvFunction->GetParameter(2) << endl;
      cout << "lambda2 = " << 1/fitInvFunction->GetParameter(3) << endl;
      cout << "c3      = " << fitInvFunction->GetParameter(4) << endl;
      cout << "lambda3 = " << 1/fitInvFunction->GetParameter(5) << endl;
      cout << "c4      = " << fitInvFunction->GetParameter(6) << endl;
      cout << "lambda4 = " << 1/fitInvFunction->GetParameter(7) << endl;
      cout << "c5      = " << fitInvFunction->GetParameter(8) << endl;
      cout << "lambda5 = " << 1/fitInvFunction->GetParameter(9) << endl;

      cout << endl << endl;      
      delete fitFunction;
      delete fitInvFunction;
   }
}
