#define findLET_cxx

#include <iostream>
#include <string.h>

#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>

#include "findLET.h"
#include "../GlobalConstants/Constants.h"

using namespace std;


void findLET::BinLogX(TH1 *h) {
   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;

   Axis_t *new_bins = new Axis_t[bins+1];

   for (int i=0; i <= bins; i++) {
      new_bins[i] = TMath::Power(10, from + i * width);
   }

   axis->Set(bins, new_bins);
   delete new_bins;
}


void findLET::Loop(Double_t energy)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   Int_t nLayers = 10;

   vector<TCanvas*> *cLET = new vector<TCanvas*>;
   cLET->reserve(nLayers);

   vector<TH1F*> *hLET = new vector<TH1F*>;
   hLET->reserve(nLayers);

   TCanvas *cLETAllLayers = new TCanvas("cLetAllLayers", "LET Landau fit for all layers", 1000, 800);
   TH1F *hLETAllLayers = new TH1F("hLETAllLayers", "LET Landau Fit for all layers", nLayers, 0, nLayers-1);

   for (Int_t layer=0; layer<nLayers; layer++) {
      cLET->push_back(new TCanvas(Form("c%d", layer), Form("LET for layer %d", layer), 1000, 800));
      hLET->push_back(new TH1F(Form("h%d", layer), Form("LET distribution for layer %d", layer), 250, 0, 10));

      hLET->at(layer)->SetXTitle("Linear Energy Transfer [kev/#mum]");
      hLET->at(layer)->SetYTitle("Number of protons");
      hLET->at(layer)->SetFillColor(kGreen-4);
      hLET->at(layer)->SetLineColor(kBlack);
   //   BinLogX(hLET->at(layer));
   //   cLET->at(layer)->SetLogx();
   }

   cout << "Made " << nLayers << " layers\n";

   gStyle->SetOptStat(0);

   Int_t lastLayer = 0;
   Float_t sumEdep = 0;

   TRandom3 *gRandom = new TRandom3();

   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }
      
      // if (parentID != 0) continue;
      // if (volumeID[4] != 4) continue;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (level1ID>=nLayers) {
         cout << "UH OH! level1ID = " << level1ID << endl;
         continue;
      }

      if (level1ID == lastLayer) { // same layer
         sumEdep += edep;
      }

      else { // new layer
         if (!sumEdep) {
            sumEdep = edep; // first layer
         }
         else {
            hLET->at(level1ID)->Fill(sumEdep*1000 / 30.);
            sumEdep = edep;
         }
      }

      lastLayer = level1ID;
   }

   for (Int_t layer=0; layer<nLayers; layer++) {
      cLET->at(layer)->cd();
      hLET->at(layer)->SetFillColor(kBlue-7);
      hLET->at(layer)->SetLineColor(kBlack);
      hLET->at(layer)->Draw();

      TF1 *landau = new TF1("landau", "landau(0)", 0, 20);
      hLET->at(layer)->Fit(landau, "Q, M, WW", "");
      hLET->at(layer)->SaveAs(Form("figures/let_degrader_layer%dMeV.root", layer));

      hLETAllLayers->Fill(layer, landau->GetParameter(1));

      cout << layer << ", " << landau->GetParameter(1) << ", " <<  landau->GetParameter(2) << endl;
   }
   cLETAllLayers->cd();
   hLETAllLayers->Draw();
}
