#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>

using namespace std;

void plotSecondaryProfileFromGATE() {
   TFile  * f = new TFile("../../DTCToolkit/Data/MonteCarlo/DTC_Final_Helium_defThres_Degrader160mm_917MeV.root");
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  x,y,z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, baseID, level1ID, PDG, lastEventID = -1;

   gStyle->SetOptStat(0);
   
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("level1ID", &level1ID);
   tree->SetBranchAddress("baseID", &baseID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("PDGEncoding", &PDG);

   TCanvas *cPDG = new TCanvas("cPDG", "PDG(z)");
   cPDG->Divide(2,1,1e-5,1e-5);
   cPDG->cd(1);
   float xfrom = -0.5;
   float xto = 40.5;
   int xlen = 41;
   
   gStyle->SetOptStat(0);

   TH1F *hPDGBeforee = new TH1F("hPDGBeforee", ";Depth in detector [Layer number];Fraction of secondary particles", xlen, xfrom, xto);
   TH1F *hPDGBeforepos = new TH1F("hPDGBeforepos", "Positrons", xlen, xfrom, xto);
   TH1F *hPDGBeforepro = new TH1F("hPDGBeforepro", ";Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGBeforeneu = new TH1F("hPDGBeforeneu", ";Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGBeforehe = new TH1F("hPDGBeforehe", "Helium", xlen, xfrom, xto);
   TH1F *hPDGBeforehe3 = new TH1F("hPDGBeforehe3", ";Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGBeforegam = new TH1F("hPDGBeforegam", "Gamma", xlen, xfrom, xto);
   TH1F *hPDGBeforeli = new TH1F("hPDGBeforeli", "Litium", xlen, xfrom, xto);
   TH1F *hPDGBeforedeu = new TH1F("hPDGBeforedeu", "Deuterium", xlen, xfrom, xto);
   TH1F *hPDGBeforetri = new TH1F("hPDGBeforetri", "Tritium", xlen, xfrom, xto);
   TH1F *hPDGBeforesi = new TH1F("hPDGBeforesi", "Silicon", xlen, xfrom, xto);
   TH1F *hPDGBeforemg = new TH1F("hPDGBeforemg", "Magnesium", xlen, xfrom, xto);
   TH1F *hPDGBeforeal = new TH1F("hPDGBeforeal", "Aluminum", xlen, xfrom, xto);
   
   TH1F *hPDGAftere = new TH1F("hPDGAftere", ";Depth in detector [Layer number];Fraction of secondary particles", xlen, xfrom, xto);
   TH1F *hPDGAfterpos = new TH1F("hPDGAfterpos", "Positrons", xlen, xfrom, xto);
   TH1F *hPDGAfterpro = new TH1F("hPDGAfterpro", ";Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGAfterneu = new TH1F("hPDGAfterneu", ";Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGAfterhe = new TH1F("hPDGAfterhe", "Helium", xlen, xfrom, xto);
   TH1F *hPDGAfterhe3 = new TH1F("hPDGAfterhe3", ";Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGAftergam = new TH1F("hPDGAftergam", "Gamma", xlen, xfrom, xto);
   TH1F *hPDGAfterli = new TH1F("hPDGAfterli", "Litium", xlen, xfrom, xto);
   TH1F *hPDGAfterdeu = new TH1F("hPDGAfterdeu", "Deuterium", xlen, xfrom, xto);
   TH1F *hPDGAftertri = new TH1F("hPDGAftertri", "Tritium", xlen, xfrom, xto);
   TH1F *hPDGAftersi = new TH1F("hPDGAftersi", "Silicon", xlen, xfrom, xto);
   TH1F *hPDGAftermg = new TH1F("hPDGAftermg", "Magnesium", xlen, xfrom, xto);
   TH1F *hPDGAfteral = new TH1F("hPDGAfteral", "Aluminum", xlen, xfrom, xto);


   Int_t l, p;

   Int_t eidHere;

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);

//      if (i>100000) break;
      l = baseID*2 + level1ID - 2;

      if (parentID > 0 && PDG == 1000020040) {
         if (edep>8e-3) {
            hPDGAfterhe->Fill(l);
            eidHere = eventID;
         }
         hPDGBeforehe->Fill(l);
         if (eidHere != eventID) {
            cout << "Missing He cluster ... l = " << l << " and eventID = " << eventID << endl;
         }

         continue;
      }

      p = PDG;

      if (p == 11) hPDGBeforee->Fill(l);
      else if (p == -11) hPDGBeforepos->Fill(l);
      else if (p == 2212) hPDGBeforepro->Fill(l);
      else if (p == 2112) hPDGBeforeneu->Fill(l);
      else if (p == 22) hPDGBeforegam->Fill(l);
      else if (p == 1000030040) hPDGBeforeli->Fill(l);
      else if (p == 1000030060) hPDGBeforeli->Fill(l);
      else if (p == 1000020030) hPDGBeforehe3->Fill(l);
      else if (p == 1000010030) hPDGBeforetri->Fill(l);
      else if (p == 1000010020) hPDGBeforedeu->Fill(l);
      else if (p == 1000130270) hPDGBeforeal->Fill(l);
      else if (p == 1000140280) hPDGBeforesi->Fill(l);
      else if (p == 1000140260) hPDGBeforesi->Fill(l);
      else if (p == 1000140299) hPDGBeforesi->Fill(l);
      else if (p == 1000120240) hPDGBeforemg->Fill(l);
      
      if (edep > 8e-3) {
         if (p == 11) hPDGAftere->Fill(l);
         else if (p == -11) hPDGAfterpos->Fill(l);
         else if (p == 2212) hPDGAfterpro->Fill(l);
         else if (p == 2112) hPDGAfterneu->Fill(l);
         else if (p == 22) hPDGAftergam->Fill(l);
         else if (p == 1000030040) hPDGAfterli->Fill(l);
         else if (p == 1000030060) hPDGAfterli->Fill(l);
         else if (p == 1000020030) hPDGAfterhe3->Fill(l);
         else if (p == 1000010030) hPDGAftertri->Fill(l);
         else if (p == 1000010020) hPDGAfterdeu->Fill(l);
         else if (p == 1000130270) hPDGAfteral->Fill(l);
         else if (p == 1000140280) hPDGAftersi->Fill(l);
         else if (p == 1000140260) hPDGAftersi->Fill(l);
         else if (p == 1000140299) hPDGAftersi->Fill(l);
         else if (p == 1000120240) hPDGAftermg->Fill(l);
      }
   }

   // NORMALIZE
   Float_t totalParticles = 100000;
   
   hPDGBeforehe->SetLineColor(kOrange); hPDGBeforepro->SetLineColor(1); hPDGBeforedeu->SetLineColor(2); hPDGBeforehe3->SetLineColor(3);
   hPDGBeforetri->SetLineColor(4); hPDGBeforee->SetLineColor(kMagenta); hPDGBeforegam->SetLineColor(6); hPDGBeforepos->SetLineColor(7);
   hPDGBeforeal->SetLineColor(14); hPDGBeforesi->SetLineColor(16); hPDGBeforemg->SetLineColor(18); hPDGBeforeli->SetLineColor(12);
   hPDGBeforeneu->SetLineColor(28); 
   hPDGBeforehe->SetLineWidth(3); hPDGBeforepro->SetLineWidth(3); hPDGBeforedeu->SetLineWidth(3); hPDGBeforehe3->SetLineWidth(3);
   hPDGBeforetri->SetLineWidth(3); hPDGBeforee->SetLineWidth(3); hPDGBeforeli->SetLineWidth(3); hPDGBeforegam->SetLineWidth(3);
   hPDGBeforeal->SetLineWidth(3); hPDGBeforesi->SetLineWidth(3); hPDGBeforemg->SetLineWidth(3); hPDGBeforepos->SetLineWidth(3);
   hPDGBeforeneu->SetLineWidth(3);

   hPDGBeforehe->Scale(1 / totalParticles);
   hPDGBeforepro->Scale(1 / totalParticles);
   hPDGBeforedeu->Scale(1 / totalParticles);
   hPDGBeforehe3->Scale(1 / totalParticles);
   hPDGBeforetri->Scale(1 / totalParticles);
   hPDGBeforee->Scale(1 / totalParticles);
   hPDGBeforegam->Scale(1 / totalParticles);
   hPDGBeforepos->Scale(1 / totalParticles);
   hPDGBeforeal->Scale(1 / totalParticles);
   hPDGBeforesi->Scale(1 / totalParticles);
   hPDGBeforemg->Scale(1 / totalParticles);
   hPDGBeforeli->Scale(1 / totalParticles);
   hPDGBeforeneu->Scale(1 / totalParticles);

   hPDGBeforepro->GetYaxis()->SetTitleOffset(1.5);
   hPDGBeforehe3->GetYaxis()->SetTitleOffset(1.5);
   
   hPDGAfterhe->SetLineColor(kOrange); hPDGAfterpro->SetLineColor(1); hPDGAfterdeu->SetLineColor(2); hPDGAfterhe3->SetLineColor(3);
   hPDGAftertri->SetLineColor(4); hPDGAftere->SetLineColor(kMagenta); hPDGAftergam->SetLineColor(6); hPDGAfterpos->SetLineColor(7);
   hPDGAfteral->SetLineColor(14); hPDGAftersi->SetLineColor(16); hPDGAftermg->SetLineColor(18); hPDGAfterli->SetLineColor(12);
   hPDGAfterneu->SetLineColor(28); 
   hPDGAfterhe->SetLineWidth(3); hPDGAfterpro->SetLineWidth(3); hPDGAfterdeu->SetLineWidth(3); hPDGAfterhe3->SetLineWidth(3);
   hPDGAftertri->SetLineWidth(3); hPDGAftere->SetLineWidth(3); hPDGAfterli->SetLineWidth(3); hPDGAftergam->SetLineWidth(3);
   hPDGAfteral->SetLineWidth(3); hPDGAftersi->SetLineWidth(3); hPDGAftermg->SetLineWidth(3); hPDGAfterpos->SetLineWidth(3);
   hPDGAfterneu->SetLineWidth(3);

   hPDGAfterhe->Scale(1 / totalParticles);
   hPDGAfterpro->Scale(1 / totalParticles);
   hPDGAfterdeu->Scale(1 / totalParticles);
   hPDGAfterhe3->Scale(1 / totalParticles);
   hPDGAftertri->Scale(1 / totalParticles);
   hPDGAftere->Scale(1 / totalParticles);
   hPDGAftergam->Scale(1 / totalParticles);
   hPDGAfterpos->Scale(1 / totalParticles);
   hPDGAfteral->Scale(1 / totalParticles);
   hPDGAftersi->Scale(1 / totalParticles);
   hPDGAftermg->Scale(1 / totalParticles);
   hPDGAfterli->Scale(1 / totalParticles);
   hPDGAfterneu->Scale(1 / totalParticles);

   hPDGAfterpro->GetYaxis()->SetTitleOffset(1.5);
   hPDGAfterhe3->GetYaxis()->SetTitleOffset(1.5);

   cPDG->cd(1);
   TLegend *legBefore = new TLegend(0.68,0.64,0.98,0.99);

   if (hPDGBeforee->Integral()) {
      hPDGBeforee->Draw("hist"); 
      legBefore->AddEntry(hPDGBeforee, "Electrons", "L");
   }
   
   hPDGBeforepro->Draw("same hist");
   legBefore->AddEntry(hPDGBeforepro, "Protons", "L");
   
   if (hPDGBeforedeu->Integral()) { 
      hPDGBeforedeu->Draw("same hist"); 
      legBefore->AddEntry(hPDGBeforedeu, "Deuterium", "L");
   }
   
   hPDGBeforehe3->Draw("same hist"); 
   legBefore->AddEntry(hPDGBeforehe3, "Helium3", "L");
   
   if (hPDGBeforetri->Integral()) { 
      hPDGBeforetri->Draw("same hist"); 
      legBefore->AddEntry(hPDGBeforetri, "Tritium", "L");
   }
   
   hPDGBeforehe->Draw("same hist"); 
   legBefore->AddEntry(hPDGBeforehe, "Helium", "L"); 
   
   if (hPDGBeforeneu->Integral()) {
      hPDGBeforeneu->Draw("same hist");
      legBefore->AddEntry(hPDGBeforeneu, "Neutrons", "L");
   }

   
   legBefore->SetTextFont(22);
   legBefore->SetTextSize(0.045);
   legBefore->Draw();
   
   /*
   TText *tBefore = new TText();
   hPDGBeforepro->GetYaxis()->SetLabelOffset(5);
   hPDGBeforehe3->GetYaxis()->SetLabelOffset(5);
   tBefore->SetTextAlign(32);
   tBefore->SetTextSize(0.05);
   tBefore->SetTextFont(22);
   Float_t at;
   for (Int_t i=0; i<=6;i++) {
      // max > i*6
      // i < max/6
      
      at = i/6. * 20; // 10*round(max(hPDGhe->GetMaximum(), max(hPDGpro->GetMaximum(), hPDGhe3->GetMaximum())/10);
      tBefore->DrawText(-0.42, at, Form("%.2f%%", at));
   }

   */
  
   cPDG->cd(2);
   TLegend *legAfter = new TLegend(0.68,0.64,0.98,0.99);
   
   if (hPDGAftere->Integral()) {
      hPDGAftere->Draw("hist"); 
      legAfter->AddEntry(hPDGAftere, "Electrons", "L");
   }
   
   hPDGAfterpro->Draw("same hist");
   legAfter->AddEntry(hPDGAfterpro, "Protons", "L");
   
   if (hPDGAfterdeu->Integral()) { 
      hPDGAfterdeu->Draw("same hist"); 
      legAfter->AddEntry(hPDGAfterdeu, "Deuterium", "L");
   }
   
   hPDGAfterhe3->Draw("same hist"); 
   legAfter->AddEntry(hPDGAfterhe3, "Helium3", "L");
   
   if (hPDGAftertri->Integral()) { 
      hPDGAftertri->Draw("same hist"); 
      legAfter->AddEntry(hPDGAftertri, "Tritium", "L");
   }
   
   hPDGAfterhe->Draw("same hist"); 
   legAfter->AddEntry(hPDGAfterhe, "Helium", "L"); 
   
   if (hPDGAfterneu->Integral()) {
      hPDGAfterneu->Draw("same hist");
      legAfter->AddEntry(hPDGAfterneu, "Neutrons", "L");
   }

   hPDGBeforee->GetYaxis()->SetRangeUser(1e-5, 300); 
   hPDGAftere->GetYaxis()->SetRangeUser(1e-5, 300); 
   
   legAfter->SetTextFont(22);
   legAfter->SetTextSize(0.045);
   legAfter->Draw();



   /*
   if (hPDGli->Integral()) { 
      hPDGli->Draw("same hist"); 
      leg->AddEntry(hPDGli, "Litium", "L");
   }
*/
/*
   if (hPDGsi->Integral()) { 
      hPDGsi->Draw("same hist"); 
      leg->AddEntry(hPDGsi, "Silicon", "L");
   }

   if (hPDGmg->Integral()) { 
      hPDGmg->Draw("same hist"); 
      leg->AddEntry(hPDGmg, "Magnesium", "L");
   }

   if (hPDGal->Integral()) { 
      hPDGal->Draw("same hist"); 
      leg->AddEntry(hPDGal, "Aluminum", "L");
   }
*/
   /*
   // ADD % TO AXIS
   TText *tAfter = new TText();
   hPDGAfterpro->GetYaxis()->SetLabelOffset(5);
   hPDGAfterhe3->GetYaxis()->SetLabelOffset(5);
   tAfter->SetTextAlign(32);
   tAfter->SetTextSize(0.05);
   tAfter->SetTextFont(22);
   for (Int_t i=0; i<=6;i++) {
      // max > i*6
      // i < max/6
      
      at = i/6. * 20; // 10*round(max(hPDGhe->GetMaximum(), max(hPDGpro->GetMaximum(), hPDGhe3->GetMaximum())/10);
      tAfter->DrawText(-0.42, at, Form("%.2f%%", at));
   }
   */

}
