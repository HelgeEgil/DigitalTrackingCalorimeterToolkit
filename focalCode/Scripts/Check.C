#define Check_cxx

#include <iostream>
#include <string.h>

#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGraph.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TMath.h>

#include "Check.h"

using namespace std;

Double_t alpha, p, alpha_prime, beta, gamma_p;

void Check::BinLogY(TH2 *h) {
   TAxis *axis = h->GetYaxis();
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

void Check::Loop(Double_t energy, Double_t sigma_mev)
{

   alpha = 0.0014467;
   alpha_prime = 0.203815;
   beta = 0.08; // validate!
   p = 1.707283;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
   TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
   TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
   c3->Divide(2,1,0.01,0.01);

   TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
   TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
   TCanvas *c6 = new TCanvas("c6", "c6", 800, 600);

   TCanvas *c7 = new TCanvas("c7", "c7", 800, 600);
   TCanvas *c8 = new TCanvas("c8", "c8", 800, 600);
   TCanvas *c9 = new TCanvas("c9", "c9", 800, 600);
   TCanvas *c11 = new TCanvas("c11", "c11", 800, 600);

   Int_t nbinsx = 1000;
   Int_t xfrom = -50;
   Int_t xto = 500;
   Int_t nbinsy = 400;
   Int_t yfrom = 1e-3;
   Int_t yto = 3;

   Int_t yfromlog = -4;
   Int_t ytolog = 1;

   Float_t x_compensate = 0; // was 45.7

   const Int_t indN = 5000;
   Float_t indX[indN];
   Float_t indY[indN];
   Float_t indZ[indN];
   Float_t indEdep[indN];
   Float_t indRange[indN];

   for (Int_t i=0; i<indN; i++) {
      indX[i] = 0; indY[i] = 0; indZ[i] = 0;
      indEdep[i] = 0; indRange[i] = 0;
   }

   Int_t indIdx = 0;

   TH2I *h2DAll = new TH2I("h2DAll", "E_{dep} vs z for all processes", nbinsx, xfrom, xto, nbinsy, yfromlog, ytolog);
   TH2I *h2DMCS = new TH2I("h2DMCS", "Multiple Coulomb Scattering", nbinsx, xfrom, xto, nbinsy, yfromlog, ytolog);
   TH2I *h2DhIoni = new TH2I("h2DhIonI", "Ionisation by hadrons", nbinsx, xfrom, xto, nbinsy, yfromlog, ytolog);
   TH2I *h2DionIoni = new TH2I("h2DionIoni", "Ionisation by ions", nbinsx, xfrom, xto, nbinsy, yfromlog, ytolog);
   TH2I *h2DPI = new TH2I("h2DPI", "Inelastic scattering of protons", nbinsx, xfrom, xto, nbinsy, yfromlog, ytolog);
   TH2I *h2DTransportation = new TH2I("h2DTransportation", "Transportation", nbinsx, xfrom, xto, nbinsy, yfromlog, ytolog);
   TH2I *h2DSL = new TH2I("h2DSL", "Step limiter", nbinsx, xfrom, xto, nbinsy, yfromlog, ytolog);

   TH1F *hZ = new TH1F("hZ", "Z profile", nbinsx, xfrom + x_compensate, xto + x_compensate);
   TH1F *hRange = new TH1F("hRange", "Primary ranges", nbinsx*3, xfrom + x_compensate, xto + x_compensate);
   TH1F *hRange2 = new TH1F("hRange2", "Primary ranges", nbinsx*3, xfrom + x_compensate, xto + x_compensate);

   BinLogY(h2DAll); BinLogY(h2DMCS); BinLogY(h2DhIoni); BinLogY(h2DionIoni); BinLogY(h2DPI); BinLogY(h2DTransportation); BinLogY(h2DSL);

   TH1I *hAllProcesses = new TH1I("hAllProcesses", "Energy deposition for all processes", 1000, 0, 15);
   TH1I *hhIoni = new TH1I("hhIoni", "Ionization and energy loss by hadrons", 1000, 0, 15);
   TH1I *hionIoni = new TH1I("hionIoni", "Ionization and energy loss by ions", 1000, 0, 15);
   TH1I *hProtonInelastic = new TH1I("hProtonInelastic", "Inelastic scattering of protons", 1000, 0, 15);
   TH1I *hMCS = new TH1I("hMCS", "Multiple Coulomb Scattering energy loss", 1000, 0, 15);

   TH1I *hFractionhIoni = new TH1I("hFractionhIoni", "Fraction of edep by energy loss", 24, 0, 24);
   TH1I *hFractionProtonInelastic = new TH1I("hFractionProtonInelastic", "Fraction of edep by inelastic scattering", 24, 0, 24);
   TH1I *hFractionMCS = new TH1I("hFractionMCS", "Fraction of edep by MCS", 24, 0, 24);
   
   TH2F *hYZ = new TH2F("hYZ", "E_{dep} for primary particles", nbinsx*10, xfrom, xto, nbinsx, -40, 40);
   TH2F *hYZ2 = new TH2F("hYZ2", "E_{dep} for secondary particles", nbinsx*10, xfrom, xto, nbinsx, -40, 40);

   TH1F *hBPpp = new TH1F("hBPpp", "Depth dose distribution by various particle types", nbinsx, xfrom, xto); // primary proton
   TH1F *hBPsp = new TH1F("hBPsp", "E_{dep} vs Z", nbinsx, xfrom, xto); // secondary proton
   TH1F *hBPsn = new TH1F("hBPsn", "E_{dep} vs Z", nbinsx, xfrom, xto); // secondary neutron
   TH1F *hBPsa = new TH1F("hBPsa", "E_{dep} vs Z", nbinsx, xfrom, xto); // secondary alpha
   TH1F *hBPs3he = new TH1F("hBPs3he", "E_{dep} vs Z", nbinsx, xfrom, xto); // secondary 3He
   TH1F *hBPsd = new TH1F("hBPsd", "E_{dep} vs Z", nbinsx, xfrom, xto); // secondary deuteron
   TH1F *hBPst = new TH1F("hBPst", "E_{dep} vs Z", nbinsx, xfrom, xto); // secondary triton
   TH1F *hBPsw = new TH1F("hBPsw", "E_{dep} vs Z", nbinsx, xfrom, xto); // secondary triton

   TH1F *hEdepSample = new TH1F("hEdepSample", "E_{dep} in sampling layers", nbinsx/2, xfrom, xto);
   TH1F *hEdepAbsorb = new TH1F("hEdepAbsorb", "E_{dep} in absorbtion layers", nbinsx/2, xfrom, xto);

   TH2F *h2DLP = new TH2F("h2DLP", "Local position for pixels", 1000, -20, 20, 1000, -20, 20);

   hFractionhIoni->SetLineColor(kBlack);
   hFractionProtonInelastic->SetLineColor(kBlack);
   hFractionMCS->SetLineColor(kBlack);

   Bool_t absorber;
   Bool_t mcs;
   Bool_t pixel;

   Bool_t khIoni;
   Bool_t kionIoni;
   Bool_t kProtonInelastic;
   Bool_t kMCS;
   Bool_t kTransportation;
   Bool_t kStepLimiter;
   
   gStyle->SetOptStat(0);

   Float_t minX = 100000;
   Float_t maxX = -1;
   Float_t minY = 100000;
   Float_t maxY = -1;

   Int_t lastEvent = -1;

   Float_t lastRange = 0;
   Int_t lastID = -1;

   cout << nentries << " entries.\n";

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) {
         cout << "Aborting run at jentry = " << jentry << endl;
         break;
      }

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      Float_t z = posZ;
      Float_t y = posY;
      Float_t x = posX;

      cout << "Hit at z = " << z << " with volumeIds " << volumeID[1] << ", " << volumeID[2] << ", " << volumeID[3] << ", " << volumeID[4] << "and level1ID " << level1ID << endl;

      if (lastID < 0) {
         lastID = eventID;
      }

      if (volumeID[4] == 0) absorber = kTRUE;
      else absorber = kFALSE;

      if (level1ID > 0) pixel = kTRUE;
      else pixel = kFALSE;

      if (!pixel)     hEdepAbsorb->Fill(z, edep);
      else           hEdepSample->Fill(z, edep);

      if (volumeID[2] == 1) {
         h2DLP->Fill(posX, posY, edep);
         hZ->Fill(z + x_compensate, edep);
      }

      if (parentID == 0) {
         if (eventID != lastID) {
            hRange->Fill(lastRange + x_compensate);
         }

         hRange2->Fill(posZ + x_compensate);
         
         lastRange = posZ;
         lastID = eventID; 
         
         if (eventID == lastID) {
            if (edep > 0.2) { // find good number
               indX[indIdx] = posX;
               indY[indIdx] = posY;
               indZ[indIdx] = posZ;
               indEdep[indIdx] = edep;
               indIdx++;
            }
         }
            // same particle
         else {
            if (lastID != -1) {
               // first event of new track, fit
               indRange[0] = 0;
               for (Int_t i=1; i<indIdx; i++) {
                  Float_t posDiff = sqrt(pow(indX[i-1] - indX[i], 2) +
                                         pow(indY[i-1] - indY[i], 2) +
                                         pow(indZ[i-1] - indZ[i], 2));
                  indRange[i] = indRange[i-1] + posDiff;
               }

               static Bool_t first = true;
               TGraph *indGraph = new TGraph(indIdx, indRange, indEdep);
               
               // cleanup
               if (first) {
                  c11->cd();
                  indGraph->Draw("A*");
                  first = false;

                  for (Int_t i=0; i<indIdx; i++) {
                     cout << "[" << indRange[i] << ", " << indEdep[i] << "], ";
                  }
                  cout << endl;

               }

               else {
                  delete indGraph;
               }

               for (Int_t i=0; i<indIdx; i++) {
                  indX[i] = 0; indY[i] = 0; indZ[i] = 0; indEdep[i] = 0; indRange[i] = 0;
               }
               
               indIdx = 0;    
            } 
         }
      }

      // 2212 protons
      // 2112 neutrons
      // 11 electrons (-11 positrons)
      // 22 photons
      // nuclei: 10LZZZAAAI
      //
      Int_t part = 0;
      Int_t partZ = 0;
      Int_t partA = 0;
      if (PDGEncoding != 22 && PDGEncoding != 2212 && PDGEncoding != 11 && PDGEncoding != 2112 && PDGEncoding != -11) {
        part = PDGEncoding - 1000000000;
        if (part > 10000000) {
           part -= part - (part%10000000);
        }

         partZ = part/10000;
         part -= partZ*10000;

         partA = part/10;
         part -= partA*10;
      }

      // then onto the processes

      TString pn = TString(processName);

      if (pn.Contains("hIoni")) khIoni = kTRUE;
      else khIoni = kFALSE;

      if (pn.Contains("ionIoni")) kionIoni = kTRUE;
      else kionIoni = kFALSE;

      if (pn.Contains("ProtonInelastic")) kProtonInelastic = kTRUE;
      else kProtonInelastic = kFALSE;

      if (pn.Contains("msc")) kMCS = kTRUE; // note spelling error in ROOT file...
      else kMCS = kFALSE;

      if (pn.Contains("Transportation")) kTransportation = kTRUE;
      else kTransportation = kFALSE;

      if (pn.Contains("StepLimiter")) kStepLimiter = kTRUE;
      else kStepLimiter = kFALSE;

      if (kMCS) h2DMCS->Fill(z, edep);
      if (kionIoni) h2DionIoni->Fill(z,edep);
      if (khIoni) h2DhIoni->Fill(z,edep);
      if (kProtonInelastic) h2DPI->Fill(z,edep);
      if (kTransportation) h2DTransportation->Fill(z,edep);
      if (kStepLimiter) h2DSL->Fill(z,edep);

      h2DAll->Fill(z,edep);

      if (parentID > 0) {
         hYZ2->Fill(z,y,edep); // was z,y
      }
      else {
         hYZ->Fill(z,y,edep);
      }

      // differentiate on particle parents here
      // Need to tag:
      // pp     Primary protons
      // sp     Secondary protons
      // sa     Secondary alphas
      // s3he   Secondary 3He
      // sd     Secondary d (Z=1 A=2)
      // st     Secondary t (Z=1 A=3)
      //
      // WHAT about electrons? Maybe I should create a new simulation
      // and dump all energy locally (OK)
      // ... or dump electron energy as their parent species?
      // Must create a LUT in that case

      if (PDGEncoding == 2212 && parentID == 0) hBPpp->Fill(z, edep); // primary proton
      if (PDGEncoding == 2212 && parentID != 0) hBPsp->Fill(z, edep); // secondary proton
      if (PDGEncoding == 2112)                  hBPsn->Fill(z, edep); // secondary neutron
      if (partZ == 2 && partA == 4)             hBPsa->Fill(z, edep); // secondary alpha
      if (partZ == 2 && partA == 3)             hBPs3he->Fill(z, edep); // secondary 3He
      if (partZ == 1 && partA == 2)             hBPsd->Fill(z, edep); // secondary deuteron
      if (partZ == 1 && partA == 3)             hBPst->Fill(z, edep); // secondary triton
      if (partZ == 74)                          hBPsw->Fill(z, edep); // scattered tungsten

      if (posX > maxX) maxX = posX;
      if (posX < minX) minX = posX;
      if (posY > maxY) maxY = posY;
      if (posY < minY) minY = posY;

      if (absorber) {
         hAllProcesses->Fill(edep);
         if (khIoni) hhIoni->Fill(edep);
         if (kionIoni) hionIoni->Fill(edep);
         if (kProtonInelastic) hProtonInelastic->Fill(edep);
         if (kMCS) hMCS->Fill(edep);
      }
   }

   cout << "X = [" << minX << "," << maxX << "].\n";
   cout << "Y = [" << minY << "," << maxY << "].\n";

   // range 1: 80 % of maximum for bragg peak on distal edge
   Double_t range_1 = hZ->GetXaxis()->GetBinCenter(hZ->FindLastBinAbove(hZ->GetMaximum() * 0.8));

   // range 2: 50 % of maximum for remainding protons plots
   Double_t range_2 = hRange2->GetXaxis()->GetBinCenter(hRange2->FindLastBinAbove(hRange2->GetMaximum() * 0.5));

   Double_t range_3 = hRange->GetXaxis()->GetBinCenter(hRange->GetMaximumBin());
   
   Double_t energy_1 = pow(range_1 / alpha, 1/p);
   Double_t energy_2 = pow(range_2 / alpha, 1/p);
   Double_t energy_3 = pow(range_3 / alpha, 1/p);

   cout << "Maximum from hRange plot: Max bin " << range_3 << " mm. (" << energy_3 << " MeV).\n";

   Double_t sigma = 0.012 * pow(range_3/10, 0.935) * 10;

   cout << "Estimated straggling from 0.012*pow(range plot range, 0.935): " << sigma << " mm.\n";

   // fit on hRange to get both mean and sigma
   TF1 *fRange = new TF1("fit_range", "gaus", 0, 900);

   hRange->Fit("fit_range", "Q, W", "", 0, 900);

   cout << "Mean range from individual range distribution: " << fRange->GetParameter(1) << " mm.\n";
   cout << "Sigma from individual range distribution: " << fRange->GetParameter(2) << " mm.\n";

   TPad *pad1  = new TPad("pad1",  "Pad 1",       0.0, 0.0, 0.4, 1.0);
   TPad *pad2d = new TPad("pad2d", "Lower pad 2", 0.4, 0.0, 0.6, 0.5);
   TPad *pad2u = new TPad("pad2u", "Upper pad 2", 0.4, 0.5, 0.6, 1.0);
   TPad *pad3d = new TPad("pad3d", "Lower pad 3", 0.6, 0.0, 0.8, 0.5);
   TPad *pad3u = new TPad("pad3u", "Upper pad 3", 0.6, 0.5, 0.8, 1.0);
   TPad *pad4u = new TPad("pad4u", "Upper pad 4", 0.8, 0.5, 1.0, 1.0);
   TPad *pad4d = new TPad("pad4d", "Lower pad 4", 0.8, 0.0, 1.0, 0.5);

   c2->cd();
   pad1->Draw(); pad2d->Draw(); pad2u->Draw(); pad3d->Draw(); 
   pad3u->Draw(); pad4u->Draw(); pad4d->Draw();
   
   h2DAll->GetZaxis()->SetLabelFont(22);
   h2DMCS->GetZaxis()->SetLabelFont(22);
   h2DionIoni->GetZaxis()->SetLabelFont(22);
   h2DhIoni->GetZaxis()->SetLabelFont(22);
   h2DPI->GetZaxis()->SetLabelFont(22);
   h2DTransportation->GetZaxis()->SetLabelFont(22);
   h2DSL->GetZaxis()->SetLabelFont(22);

   pad1->cd();
      pad1->SetLogz();
      pad1->SetLogy();
      h2DAll->SetXTitle("Local Z position [mm]");
      h2DAll->SetYTitle("Energy deposition [AU]");
      h2DAll->Draw("COLZ");

   pad2u->cd();
      pad2u->SetLogz();
      pad2u->SetLogy();
      h2DMCS->SetXTitle("Local Z position [mm]");
      h2DMCS->SetYTitle("Energy deposition [AU]");
      h2DMCS->Draw("COLZ");
      
   pad2d->cd();
      pad2d->SetLogz();
      pad2d->SetLogy();
      h2DionIoni->SetXTitle("Local Z position [mm]");
      h2DionIoni->SetYTitle("Energy deposition [AU]");
      h2DionIoni->Draw("COLZ");

   pad3u->cd();
      pad3u->SetLogz();
      pad3u->SetLogy();
      h2DhIoni->SetXTitle("Local Z position [mm]");
      h2DhIoni->SetYTitle("Energy deposition [AU]");
      h2DhIoni->Draw("COLZ");

   pad3d->cd();
      pad3d->SetLogz();
      pad3d->SetLogy();
      h2DPI->SetXTitle("Local Z position [mm]");
      h2DPI->SetYTitle("Energy deposition [AU]");
      h2DPI->Draw("COLZ");

   pad4u->cd();
      pad4u->SetLogz();
      pad4u->SetLogy();
      h2DTransportation->SetXTitle("Local Z position [mm]");
      h2DTransportation->SetYTitle("Energy deposition [AU]");
      h2DTransportation->Draw("COLZ");

   pad4d->cd();
      pad4d->SetLogz();
      pad4d->SetLogy();
      h2DSL->SetXTitle("Local Z position [mm]");
      h2DSL->SetYTitle("Energy deposition [AU]");
      h2DSL->Draw("COLZ");

   c1->cd();
   hAllProcesses->SetXTitle("Energy deposition");
   hAllProcesses->SetYTitle("Number of interactions");
   hhIoni->SetLineColor(kRed);
   hionIoni->SetLineColor(kRed-2);
   hProtonInelastic->SetLineColor(kGreen);
   hMCS->SetLineColor(kYellow+2);
   c1->SetLogy();

   TLegend *leg = new TLegend(0.7, 0.8, 0.95, 0.95);
   leg->AddEntry(hAllProcesses, "All processes", "f");
   leg->AddEntry(hhIoni, "Energy loss and ionization of hadrons", "f");
   leg->AddEntry(hionIoni, "Energy loss and ionization of ions", "f");
   leg->AddEntry(hProtonInelastic, "Inelastic scattering of protons", "f");
   leg->AddEntry(hMCS, "Energy loss from multiple coulomb scattering", "f");


//   THStack *hs = new THStack("hs", "Stacked 1D histograms");
//   hs->Add(hFractionhIoni);
//   hs->Add(hFractionProtonInelastic);
//   hs->Add(hFractionMCS);
   
   
   hAllProcesses->Draw();
   hhIoni->Draw("same");
   hionIoni->Draw("same");
   hProtonInelastic->Draw("same");
   hMCS->Draw("same");
   
   leg->Draw();

   c3->cd(1);
      c3->SetLogz();
      hYZ->SetXTitle("Z position [mm]");
      hYZ->SetYTitle("Y position [mm]");
      hYZ->Draw("COLZ");
   c3->cd(2);
      c3->SetLogz();
      hYZ2->SetXTitle("Z position [mm]");
      hYZ2->SetYTitle("Y position [mm]");
      hYZ2->Draw("COLZ");

   c4->cd();
      c4->SetLogy();

      hBPpp->SetXTitle("Z position [mm]");
      hBPpp->SetYTitle("Relative maximum E_{dep} [%]");

      // scale everything so that hBPpp at max is 100 %

      Float_t max_edep = hBPpp->GetMaximum();
   
      hBPpp->Scale(100./max_edep);
      hBPsp->Scale(100./max_edep);
      hBPsn->Scale(100./max_edep);
      hBPsa->Scale(100./max_edep);
      hBPs3he->Scale(100./max_edep);
      hBPsd->Scale(100./max_edep);
      hBPst->Scale(100./max_edep);
      hBPsw->Scale(100./max_edep);
      
      hBPpp->GetYaxis()->SetRangeUser(1e-5, 2e2);
      hBPpp->GetYaxis()->SetNoExponent();

      hBPpp->SetFillColor(kOrange+10);
      hBPsp->SetFillColor(kOrange+8);
      hBPsn->SetFillColor(kOrange+6);
      hBPsa->SetFillColor(kOrange+4);
      hBPsd->SetFillColor(kOrange+2);
      hBPs3he->SetFillColor(kPink+10);
      hBPst->SetFillColor(kPink+8);
      hBPsw->SetFillColor(kPink+6);
      hBPpp->SetLineColor(kBlack);
      hBPsp->SetLineColor(kBlack);
      hBPsn->SetLineColor(kBlack);
      hBPsa->SetLineColor(kBlack);
      hBPs3he->SetLineColor(kBlack);
      hBPsd->SetLineColor(kBlack);
      hBPst->SetLineColor(kBlack);
      hBPsw->SetLineColor(kBlack);
      
      TLegend *leg3 = new TLegend(0.7, 0.8, 0.95, 0.95);
      leg3->AddEntry(hBPpp, "Primary protons", "f");
      leg3->AddEntry(hBPsp, "Secondary protons", "f");
      leg3->AddEntry(hBPsa, "Secondary alphas", "f");
      leg3->AddEntry(hBPs3he, "Secondary ^{3}He", "f");
      leg3->AddEntry(hBPsd, "Secondary deuterons", "f");
      leg3->AddEntry(hBPst, "Secondary tritons", "f");
      leg3->AddEntry(hBPsw, "Scattered tungsten", "f");
      leg3->AddEntry(hBPsn, "Secondary neutrons", "f");

      hBPpp->Draw();
      hBPsp->Draw("same");
      hBPsa->Draw("same");
      hBPs3he->Draw("same");
      hBPsd->Draw("same");
      hBPst->Draw("same");
      hBPsw->Draw("same");
      hBPsn->Draw("same");
      leg3->Draw();

   TLegend *leg2 = new TLegend(0.7, 0.8, 0.95, 0.95);
   leg2->AddEntry(hEdepAbsorb, "E_{dep} outside sampling layers", "f");
   leg2->AddEntry(hEdepSample, "E_{dep} inside sampling layess", "f");

   c5->cd();
      c5->SetLogy();
      hEdepAbsorb->SetXTitle("Z position [mm]");
      hEdepAbsorb->SetYTitle("Total E_{dep} [MeV]");
      hEdepAbsorb->SetFillColor(kRed-2);
      hEdepAbsorb->SetLineColor(kBlack);
      hEdepSample->SetFillColor(kYellow-2);
      hEdepSample->SetLineColor(kBlack);
      hEdepAbsorb->Draw();
      hEdepSample->Draw("same");
      leg2->Draw();

   c6->cd();
      h2DLP->SetXTitle("X [mm]");
      h2DLP->SetYTitle("Y [mm]");
      h2DLP->Draw("COLZ");



   c7->cd();
      hZ->SetXTitle("Z [mm]");
      hZ->SetYTitle("Edep [MeV]");
      hZ->SetFillColor(kBlue-7);
      hZ->SetLineColor(kBlack);
      hZ->Draw();

   c8->cd();
      hRange->SetXTitle("Range [mm]");
      hRange->SetYTitle("Number of primaries");
      hRange->SetFillColor(kBlue-7);
      hRange->SetLineColor(kBlack);
      hRange->Draw();

   c9->cd();
      hRange2->SetXTitle("Range [mm]");
      hRange2->SetYTitle("Number of primaries");
      hRange2->SetFillColor(kBlue-7);
      hRange2->SetLineColor(kBlack);
      hRange2->Draw();

//   c2->cd();
//   hs->Draw();

   /*
   delete hAllProcesses;
   delete hhIoni;
   delete hionIoni;
   delete hProtonInelastic;
   delete hMCS;
   delete hYZ;
   delete hYZ2;
   delete hBP;
   
   delete h2DAll;
   delete h2DMCS;
   delete h2DionIoni;
   delete h2DhIoni;
   delete h2DPI;
   delete h2DTransportation;
   delete h2DSL;
   
   delete hFractionhIoni;
   delete hFractionProtonInelastic;
   delete hFractionMCS;
   */
}
