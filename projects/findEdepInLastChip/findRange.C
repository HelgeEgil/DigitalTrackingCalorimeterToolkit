#define findRange_cxx
#include "findRange.h"
#include <TPaveStats.h>
#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <vector>
#include <iostream>
#include <math.h>
#include <cmath>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <iostream>
#include <fstream>
#include <TPad.h>
#include <TMath.h>

using namespace std;

vector<Float_t> findRange::Run(Double_t energy, Double_t sigma_mev)
{
   vector<Float_t> returnValues;
   if (fChain == 0) return returnValues;
   if (fChain->GetEntries() == 0) return returnValues;

   Float_t degraderThickness = run_degraderThickness;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   Bool_t useDegrader = (degraderThickness > 0) ? true : false;

   TCanvas *c = new TCanvas("c", "edep in last chip", 1200, 900);
//   c->Divide(2,1,0.0001,0.0001);

   TCanvas *c2 = new TCanvas("c2", "Range of protons", 1200, 900);
   c2->Divide(2,1,0.0001,0.0001);

   TH1F *hInelastic = new TH1F("hInelastic", "Protons ending in inelastic collisions;Edep in last traversed chip[keV/#mum];Entries", 400, 0, 20);
      TH1F *hStopping = new TH1F("hStopping", Form("Energy deposition in last traversed chip for %.0f mm degrader;Edep in last traversed chip [keV/#mum];Entries", degraderThickness), 400, 0, 20);

   TH1F *hRangeAll = new TH1F("hRangeAll", "Range of all protons;Range [mm];Entries", 200, 0, 200);
   TH1F *hRangeStopping = new TH1F("hRangeStopping", "Range of stopping protons;Range [mm];Entries", 200, 0, 200);

   Int_t lastEvent = -1;
   Float_t lastEdep = -1;
   Float_t lastZ = 0;
   Float_t lastRange = -1;

   Int_t lastID = -1;
   Float_t dE = 0;
   Char_t lastProcessName[17];
   
   Long64_t ientry = LoadTree(0);
   fChain->GetEntry(0);
   
   if (nentries == 0) return returnValues;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (posZ < 0) continue;

      if (lastID < 0) lastID = eventID;

      if (parentID == 0) {
         if (eventID != lastID) {
            hRangeAll->Fill(lastRange);
            if (lastProcessName[0] == 'P')   {
               hInelastic->Fill(lastEdep);
            }

            else {
               hStopping->Fill(lastEdep);
               hRangeStopping->Fill(lastRange);
            }
         }

         lastID = eventID;
         if (volumeID[4] == 0) { // chip
            lastZ = posZ;
            lastEdep = edep * 1000 / 14;
         }
         
         for (Int_t j=0; j<17; j++) {
            lastProcessName[j] = processName[j];
         }
         lastRange = posZ;
      }
   }

   cout << lastID << endl;
   c->cd();

   gStyle->SetOptStat(0);
   
   hStopping->GetXaxis()->SetTitleFont(22);
   hStopping->GetXaxis()->SetLabelFont(22);
   hStopping->GetYaxis()->SetTitleFont(22);
   hStopping->GetYaxis()->SetLabelFont(22);

   hStopping->SetFillColor(kRed+2);
   hStopping->Draw();

   hInelastic->SetFillColor(kRed-7);
   hInelastic->Draw("same");

   hStopping->GetYaxis()->SetRangeUser(0, fmax(hStopping->GetMaximum(), hInelastic->GetMaximum()) * 1.2);

   Int_t ixbin = 0;
   Int_t fromBin = hStopping->GetXaxis()->FindBin(1.2);
   for (Int_t i=fromBin; i<hStopping->GetNbinsX(); i++) {
      if (hStopping->GetBinContent(i) > hInelastic->GetBinContent(i)) {
         ixbin = i;
         break;
      }
   }
   gPad->Update();
   Float_t ixbinpos = hStopping->GetXaxis()->GetBinCenter(ixbin) - 0.025;
   printf("ixbin = %d, pos = %.1f.\n", ixbin, ixbinpos);
   TLine *line = new TLine(ixbinpos, 0, ixbinpos, gPad->GetUymax());
   line->Draw();

   TF1 *fStop = new TF1("fStop", "landau");
   TF1 *fInel = new TF1("fInel", "landau");
   hStopping->Fit(fStop);
   hInelastic->Fit(fInel);

/*
   c->Update();
   TPaveStats *ps = (TPaveStats*) c->GetPrimitive("stats");
   hStopping->SetBit(TH1::kNoStats);
   ps->SetY1NDC(0.53); ps->SetX1NDC(0.72);
//   ps->SetTextFont(22);
   ps->AddText(Form("Stopping MPV = %.2f keV/#mum", fStop->GetParameter(1)));
   ps->AddText(Form("Inelastic MPV = %.2f keV/#mum", fInel->GetParameter(1)));
   ps->AddText(Form("Separation point = %.2f keV/#mum", ixbinpos));
*/
   TLegend *l = new TLegend(0.5, 0.5, 0.7, 0.7);
   l->AddEntry(hStopping, "Stopping protons", "F");
   l->AddEntry(hInelastic, "Nuclear interactions", "F");
   l->SetTextFont(22);
   l->Draw();


   c2->cd(1);
   hRangeAll->SetFillColor(kRed+2);
   hRangeAll->Draw();
   c2->cd(2);
   hRangeStopping->SetFillColor(kRed-7);
   hRangeStopping->Draw();
}
