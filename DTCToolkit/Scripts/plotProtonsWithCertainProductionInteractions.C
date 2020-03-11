#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TString.h>
#include <TEllipse.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSpline.h>
#include <TStyle.h>

using namespace std;

void plotProtonsWithCertainProductionInteractions() {
   TFile  * f = new TFile("../../DTCToolkit/Data/MonteCarlo/DTC_Full_Final_Helium_Degrader160mm_917MeV.root");
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  x,y,z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, baseID, level1ID, PDG, trackID, lastEventID = -1;
   Char_t   processName[17];

   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("level1ID", &level1ID);
   tree->SetBranchAddress("trackID", &trackID);
   tree->SetBranchAddress("baseID", &baseID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("PDGEncoding", &PDG);
   tree->SetBranchAddress("processName", processName);

   TCanvas *cPDG = new TCanvas("cPDG", "PDG(z)", 1200, 700);
   cPDG->Divide(2,1,1e-5,1e-5);
   float xfrom = -0.5;
   float xto = 40.5;
   int xlen = 41;
   
   gStyle->SetOptStat(0);

   TH1F *hPDGAllPro = new TH1F("hPDGAllPro", "Protons;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGNeutronPro = new TH1F("hPDNeutronPro", "Protons from nuclear recoil;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGHePro = new TH1F("hPDGHePro", "Protons from projectile fragmentation;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGNeutronInelPro = new TH1F("hPDNeutronInelPro", "Protons from nuclear recoil;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   
   TH1F *hPDGAllE = new TH1F("hPDGAllE", "Electrons;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGGammaE = new TH1F("hPDGGammaE", "Electrons from gamma conversion;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGGammaPhotE = new TH1F("hPDGGammaPhotE", "Electrons from gamma conversion;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGGammaRaylE = new TH1F("hPDGGammaRaylE", "Electrons from gamma conversion;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGGammaComptE = new TH1F("hPDGGammaComptE", "Electrons from gamma conversion;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);
   TH1F *hPDGGammaConvE = new TH1F("hPDGGammaConvE", "Electrons from gamma conversion;Depth in detector [layer number];Percentage of secondary species", xlen, xfrom, xto);

   Int_t l, p;
   Int_t gammaEID, neutronEID = -1;

   Int_t gammaParent, neutronParent, heliumParent;
   Int_t goodGammaEvent, goodNeutronEvent, goodGammaPhotEvent, goodGammaComptEvent, goodGammaRaylEvent, goodGammaConvEvent, goodHeliumEvent, goodNeutronInelasticEvent;

   Double_t xp[1000];
   Double_t yp[1000];
   Double_t zp[1000];
   Int_t idxp = 0;
   TSpline3 *splinexp = nullptr;
   TSpline3 *splineyp = nullptr;
   Double_t diffxy;

   for (Int_t i=0; i<tree->GetEntries(); i++) {
      // Must find correct trackid also
      
      if (i>10000) break;
      
      tree->GetEntry(i);
   

      if (parentID == 0) {
         xp[idxp] = x;
         yp[idxp] = y;
         zp[idxp++] = z;
      }
      else if (parentID>0 && idxp>0) {
         
         if (splinexp) {
            delete splinexp;
            delete splineyp;
         }

         splinexp = new TSpline3("splinexp", zp, xp, idxp);
         splineyp = new TSpline3("splineyp", zp, yp, idxp);
         idxp = 0;
      }

//      cout << "z " << z << ", parentID " << parentID << ", trackID " << trackID << ", processname " << processName << ", PDG " << PDG << endl;

      l = baseID*2 + level1ID - 2; // control this

      if (parentID > 0) {
         diffxy = sqrt(pow(splinexp->Eval(z) - x, 2) + pow(splineyp->Eval(z) - y, 2));
         if (diffxy < 0.025) continue;
      }

      if (PDG == 1000020040) {
         heliumParent = trackID;
         if (TString(processName) == "alphaInelastic") goodHeliumEvent = eventID;
      }

      if (PDG == 2112) { 
         neutronParent = trackID;
         goodNeutronEvent = eventID;
         if (TString(processName) == "neutronInelastic") goodNeutronInelasticEvent = eventID;
      }
      if (PDG == 22) {
         gammaParent = trackID;
         goodGammaEvent = eventID;
         if (TString(processName) == "phot") goodGammaPhotEvent = eventID;
         if (TString(processName) == "comp") goodGammaComptEvent = eventID;
         if (TString(processName) == "Rayl") goodGammaRaylEvent = eventID;
         if (TString(processName) == "conv") goodGammaConvEvent = eventID;
      }

      if (heliumParent == parentID && eventID == goodHeliumEvent) {
         if (PDG == 2212) {
            hPDGHePro->Fill(l);
         }
      }

      if (gammaParent == parentID && eventID == goodGammaEvent) {
         if (PDG == 11) {
            hPDGGammaE->Fill(l);
            if (eventID == goodGammaPhotEvent) {
               hPDGGammaPhotE->Fill(l);
            }
            if (eventID == goodGammaComptEvent) {
               hPDGGammaComptE->Fill(l);
            }
            if (eventID == goodGammaRaylEvent) {
               hPDGGammaRaylE->Fill(l);
            }
            if (eventID == goodGammaConvEvent) {
               hPDGGammaConvE->Fill(l);
            }
         }
      }
      if (neutronParent == parentID && eventID == goodNeutronEvent) {
         if (PDG == 2212) {
            if (eventID != goodNeutronInelasticEvent) {
               hPDGNeutronPro->Fill(l);
            }
            else {
               hPDGNeutronInelPro->Fill(l);
            }

         }
      }

      if (PDG == 2212) {
         hPDGAllPro->Fill(l);
      }
      if (PDG == 11) {
         hPDGAllE->Fill(l);
      }
   }

   // NORMALIZE
   Float_t totalParticles = eventID;
   
   hPDGAllPro->Scale(1 / totalParticles);
   hPDGNeutronPro->Scale(1 / totalParticles);
   hPDGHePro->Scale(1 / totalParticles);
   hPDGNeutronInelPro->Scale(1 / totalParticles);
   hPDGAllE->Scale(1 / totalParticles);
   hPDGGammaE->Scale(1 / totalParticles);
   hPDGGammaPhotE->Scale(1 / totalParticles);
   hPDGGammaConvE->Scale(1 / totalParticles);
   hPDGGammaComptE->Scale(1 / totalParticles);
   hPDGGammaRaylE->Scale(1 / totalParticles);
   
   hPDGAllPro->SetLineColor(kBlue);
   hPDGHePro->SetLineColor(kGreen);
   hPDGNeutronPro->SetLineColor(28);
   hPDGNeutronInelPro->SetLineColor(kRed);
   hPDGAllPro->SetLineWidth(3);
   hPDGHePro->SetLineWidth(3);
   hPDGNeutronPro->SetLineWidth(3);
   hPDGNeutronInelPro->SetLineWidth(3);
   
   hPDGAllE->SetLineColor(kBlue);
   hPDGGammaE->SetLineColor(kRed);
   hPDGGammaPhotE->SetLineColor(kGreen);
   hPDGGammaConvE->SetLineColor(kOrange);
   hPDGGammaRaylE->SetLineColor(kYellow);
   hPDGGammaComptE->SetLineColor(28);
   hPDGAllE->SetLineWidth(3);
   hPDGGammaE->SetLineWidth(3);
   hPDGGammaPhotE->SetLineWidth(3);
   hPDGGammaConvE->SetLineWidth(3);
   hPDGGammaRaylE->SetLineWidth(3);
   hPDGGammaComptE->SetLineWidth(3);

   cPDG->cd(1);
   TLegend *legPro = new TLegend(0.39,0.77,0.95,0.9);

   hPDGAllPro->Draw("hist"); 
   legPro->AddEntry(hPDGAllPro, "All protons", "L");

   hPDGHePro->Draw("same hist");
   legPro->AddEntry(hPDGHePro, "Helium fragmentation", "L");

   hPDGAllPro->GetYaxis()->SetRangeUser(1e-5, 300);
   hPDGAllE->GetYaxis()->SetRangeUser(1e-5, 300);

   hPDGNeutronPro->Draw("same hist");
   legPro->AddEntry(hPDGNeutronPro, "Neutron recoil", "L");
   
   hPDGNeutronInelPro->Draw("same hist");
   legPro->AddEntry(hPDGNeutronInelPro, "Neutron inelastic", "L");
   
   legPro->SetTextFont(22);
   legPro->SetTextSize(0.045);
   legPro->Draw();

   cPDG->cd(2);
   TLegend *legE = new TLegend(0.39,0.67,0.95,0.9);

   hPDGAllE->Draw("hist"); 
   legE->AddEntry(hPDGAllE, "All electrons", "L");

   hPDGGammaE->Draw("same hist");
   legE->AddEntry(hPDGGammaE, "Gamma parent", "L");
  
   hPDGGammaRaylE->Draw("same hist");
   legE->AddEntry(hPDGGammaRaylE, "-> Rayleigh", "L");
   hPDGGammaPhotE->Draw("same hist");
   legE->AddEntry(hPDGGammaPhotE, "-> Photoabsorption", "L");
   hPDGGammaConvE->Draw("same hist");
   legE->AddEntry(hPDGGammaConvE, "-> Photoconversion", "L");
   hPDGGammaComptE->Draw("same hist");
   legE->AddEntry(hPDGGammaComptE, "-> Compton", "L");

   legE->SetTextFont(22);
   legE->SetTextSize(0.045);
   legE->Draw();
}
