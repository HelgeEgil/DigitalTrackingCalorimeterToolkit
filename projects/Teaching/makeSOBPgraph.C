#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>

using namespace std;

float getw(float k, float p) {
   float w = 0;
   float n = 10;
   if (k == 0) {
      w = 1 - pow(1 - 1/(2*n), 1-1/p);
   }
   else if (k == 10) {
      w = pow(1/(2*n), 1 - 1/p);
   }
   else {
      w = pow(1 - (1/n) * (k - 0.5), 1-1/p) - pow(1 - (1/n) * (k + 0.5), 1-1/p);
      
//      if (k == 10) w *= 0.9;
      if (k == 9) w *= 0.95;

//      if (k<5) w*= 0.8;
//      if (k<2) w*= 0.8;
/*
      if (k < 10) w *= 0.99;
      if (k < 9) w *= 0.99;
      if (k < 8) w *= 0.98;
      if (k < 7) w *= 0.95;
*/
      if (k < 6) w *= 0.98;
      if (k < 5) w *= 0.97;
      if (k < 4) w *= 0.96;
      if (k < 3) w *= 0.93;
      if (k < 2) w *= 0.92;
      if (k < 1) w *= 0.9;

      
   }

   return w;
}

void makeSOBPgraph() {
   TTree * tree;

   gStyle->SetTextFont(22);
   gStyle->SetOptStat(0);

   TFile  * f0 = new TFile(Form("output/sobp_10um_200MeV.root"));
   TFile  * f1 = new TFile(Form("output/sobp_200MeV.root"));
   TFile  * f2 = new TFile(Form("output/sobp_195.44011MeV.root"));
   TFile  * f3 = new TFile(Form("output/sobp_190.7968MeV.root"));
   TFile  * f4 = new TFile(Form("output/sobp_186.0648MeV.root"));
   TFile  * f5 = new TFile(Form("output/sobp_181.2383MeV.root"));
   TFile  * f6 = new TFile(Form("output/sobp_176.3107MeV.root"));
   TFile  * f7 = new TFile(Form("output/sobp_171.2746MeV.root"));
   TFile  * f8 = new TFile(Form("output/sobp_166.1219MeV.root"));
   TFile  * f9 = new TFile(Form("output/sobp_160.8430MeV.root"));
   TFile  * f10 = new TFile(Form("output/sobp_155.4272MeV.root"));
   TFile  * f11 = new TFile(Form("output/sobp_149.8621MeV.root"));

   Int_t    maxEntries = -5000; //300000;

   Float_t  z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, lastEventID = -1;
   Float_t  expectedRange = 0.022 * pow(200, 1.77);
   Float_t  actualRange;
   Float_t  rangeFactor;

   Float_t  p = 1.55;
   Float_t  w = 0;
   Float_t  k = 0, n = 10;

//   Float_t w1 = 0.50091301, w2 = 0.14454384, w3 = 0.08075458, w4 = 0.05863545, w5 = 0.04686350, w6 = 0.03942108, w7 = 0.03423867, w8 = 0.03039743, w9 = 0.02742265, w10 = 0.02504266, w11 = 0.011761713;
//   Float_t w1 = 0.272, w2 = 0.166, w3 = 0.109, w4 = 0.0862, w5 = 0.0731, w6 = 0.0645, w7 = 0.0581, w8 = 0.0532, w9 = 0.0494, w10 = 0.0462, w11 = 0.0221;

   printf("expected range*1.07 = %.2f mm\n", expectedRange*1.07);
   Int_t nbins = 278;


   TH1F   * individualDose = new TH1F("individualDose", "Single proton;Range [mm];Relative dose [%]", nbins, 0, expectedRange*1.07);
   TH1F   * beamDose = new TH1F("beamDose", "Proton beam;Range [mm];Relative dose [%]", nbins, 0, expectedRange*1.07);
   TH1F   * spreadOutDose = new TH1F("spreadOutDose", "Spread Out Bragg Peak;Range [mm];Relative dose [%]", nbins, 0, expectedRange*1.07);

   TH1F   * sobp0 = new TH1F("sobp0", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp1 = new TH1F("sobp1", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp2 = new TH1F("sobp2", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp3 = new TH1F("sobp3", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp4 = new TH1F("sobp4", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp5 = new TH1F("sobp5", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp6 = new TH1F("sobp6", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp7 = new TH1F("sobp7", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp8 = new TH1F("sobp8", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp9 = new TH1F("sobp9", ";;;", nbins, 0, expectedRange*1.07);
   TH1F   * sobp10 = new TH1F("sobp10", ";;;", nbins, 0, expectedRange*1.07);

   tree = (TTree*) f0->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      individualDose->Fill(z, edep);
   }
   printf("Finished tree 0\n");


   tree = (TTree*) f1->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);

   printf("Starting tree 1 with %d entries...\n", tree->GetEntries());
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if (maxEntries > 0 && eventID > maxEntries) break;

      beamDose->Fill(z, edep);
      spreadOutDose->Fill(z, edep * getw(10, p));
      sobp10->Fill(z, edep * getw(10, p));
   }
   printf("Finished tree 1\n");

   tree = (TTree*) f2->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(9, p));
      sobp9->Fill(z, edep * getw(9, p));
   }
   printf("Finished tree 2\n");

   tree = (TTree*) f3->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(8, p));
      sobp8->Fill(z, edep * getw(8, p));
   }

   tree = (TTree*) f4->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(7, p));
      sobp7->Fill(z, edep * getw(7, p));
   }
   printf("Finished tree 4\n");

   tree = (TTree*) f5->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(6, p));
      sobp6->Fill(z, edep * getw(6, p));
   }

   tree = (TTree*) f6->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(5, p));
      sobp5->Fill(z, edep * getw(5, p));
   }

   tree = (TTree*) f7->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(4, p));
      sobp4->Fill(z, edep * getw(4, p));
   }
   printf("Finished tree 7\n");

   tree = (TTree*) f8->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(3, p));
      sobp3->Fill(z, edep * getw(3, p));
   }

   tree = (TTree*) f9->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(2, p));
      sobp2->Fill(z, edep * getw(2, p));
   }

   tree = (TTree*) f10->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(1, p));
      sobp1->Fill(z, edep * getw(1, p));
   }

   tree = (TTree*) f11->Get("Hits");
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i); 
      if (maxEntries > 0 && eventID > maxEntries) break;
      spreadOutDose->Fill(z, edep * getw(0, p));
      sobp0->Fill(z, edep * getw(0, p));
   }

   TCanvas *c = new TCanvas("c", "Different proton dose distributions", 1300, 500);
   c->Divide(3, 0, 0.001, 0.001);
   c->cd(1);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.03);
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.03);
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   gPad->SetRightMargin(0.03);
   
   sobp0->SetLineColor(kGray+3);
   sobp1->SetLineColor(kGray+3);
   sobp2->SetLineColor(kGray+3);
   sobp3->SetLineColor(kGray+3);
   sobp4->SetLineColor(kGray+3);
   sobp5->SetLineColor(kGray+3);
   sobp6->SetLineColor(kGray+3);
   sobp7->SetLineColor(kGray+3);
   sobp8->SetLineColor(kGray+3);
   sobp9->SetLineColor(kGray+3);
   sobp10->SetLineColor(kGray+3);

   individualDose->GetXaxis()->SetTitleFont(22);
   individualDose->GetYaxis()->SetTitleFont(22);
   individualDose->GetXaxis()->SetLabelFont(22);
   individualDose->GetYaxis()->SetLabelFont(22);
   individualDose->GetXaxis()->SetTitleSize(0.05);
   individualDose->GetYaxis()->SetTitleSize(0.05);
   individualDose->GetXaxis()->SetLabelSize(0.05);
   individualDose->GetYaxis()->SetLabelSize(0.05);
   individualDose->GetYaxis()->SetTitleOffset(1.3);

   beamDose->GetXaxis()->SetTitleFont(22);
   beamDose->GetYaxis()->SetTitleFont(22);
   beamDose->GetXaxis()->SetLabelFont(22);
   beamDose->GetYaxis()->SetLabelFont(22);
   beamDose->GetXaxis()->SetTitleSize(0.05);
   beamDose->GetYaxis()->SetTitleSize(0.05);
   beamDose->GetXaxis()->SetLabelSize(0.05);
   beamDose->GetYaxis()->SetLabelSize(0.05);
   beamDose->GetYaxis()->SetTitleOffset(1.3);

   spreadOutDose->GetXaxis()->SetTitleFont(22);
   spreadOutDose->GetYaxis()->SetTitleFont(22);
   spreadOutDose->GetXaxis()->SetLabelFont(22);
   spreadOutDose->GetYaxis()->SetLabelFont(22);
   spreadOutDose->GetXaxis()->SetTitleSize(0.05);
   spreadOutDose->GetYaxis()->SetTitleSize(0.05);
   spreadOutDose->GetXaxis()->SetLabelSize(0.05);
   spreadOutDose->GetYaxis()->SetLabelSize(0.05);
   spreadOutDose->GetYaxis()->SetTitleOffset(1.3);

   individualDose->Scale(100/individualDose->GetMaximum());
   beamDose->Scale(100/beamDose->GetMaximum());
   sobp0->Scale(100/spreadOutDose->GetMaximum());
   sobp1->Scale(100/spreadOutDose->GetMaximum());
   sobp2->Scale(100/spreadOutDose->GetMaximum());
   sobp3->Scale(100/spreadOutDose->GetMaximum());
   sobp4->Scale(100/spreadOutDose->GetMaximum());
   sobp5->Scale(100/spreadOutDose->GetMaximum());
   sobp6->Scale(100/spreadOutDose->GetMaximum());
   sobp7->Scale(100/spreadOutDose->GetMaximum());
   sobp8->Scale(100/spreadOutDose->GetMaximum());
   sobp9->Scale(100/spreadOutDose->GetMaximum());
   sobp10->Scale(100/spreadOutDose->GetMaximum());
   spreadOutDose->Scale(100/spreadOutDose->GetMaximum());

   individualDose->SetFillColor(kRed-4);
   beamDose->SetFillColor(kRed-4);
   spreadOutDose->SetFillColor(kRed-4);
   
   individualDose->SetLineColor(kBlack);
   beamDose->SetLineColor(kBlack);
   spreadOutDose->SetLineColor(kBlack);

   c->cd(1);
   individualDose->Draw();

   c->cd(2);
   beamDose->Draw();

   c->cd(3);
   spreadOutDose->Draw();
   sobp0->Draw("same");
   sobp1->Draw("same");
   sobp2->Draw("same");
   sobp3->Draw("same");
   sobp4->Draw("same");
   sobp5->Draw("same");
   sobp6->Draw("same");
   sobp7->Draw("same");
   sobp8->Draw("same");
   sobp9->Draw("same");
   sobp10->Draw("same");

}
