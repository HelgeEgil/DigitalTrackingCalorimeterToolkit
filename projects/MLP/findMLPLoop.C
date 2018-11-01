#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <fstream>
#include <TBranch.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TStyle.h>
#include <Math/Vector3D.h>

void findMLPLoop(Float_t phantomSize = 200, Int_t eventsToUse = -1, Float_t spotSize = -1);

using namespace std;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > XYZVector;

XYZVector SplineMLP(Double_t t, XYZVector X0, XYZVector X2, XYZVector P0, XYZVector P2, Double_t Lambda0, Double_t Lambda1) {
   XYZVector P0Lambda, P2Lambda, X2mX0, S;
   Float_t tt = pow(t, 2), ttt = pow(t, 3);

   X2mX0 = X2 - X0;
   P0Lambda = P0 * Lambda0 * sqrt(X2mX0.Mag2());
   P2Lambda = P2 * Lambda1 * sqrt(X2mX0.Mag2());

   S = (2*ttt - 3*tt + 1) * X0 + (ttt - 2*tt + t) * P0Lambda + (-2*ttt + 3*tt) * X2 + (ttt - tt) * P2Lambda;

   return S;
}

void findMLPLoop(Float_t phantomSize, Int_t eventsToUse, Float_t spotSize) {
   Float_t     initialEnergy = 230;
   Float_t     differenceArrayDZ = 3;
   if (eventsToUse < 0) eventsToUse = 100;
   Float_t     x, y, z, edep, sum_edep = 0, residualEnergy = 0;
   Int_t       eventID, parentID, lastEID = -1;
   XYZVector   Xp0, Xp1, Xp2, Xp3, X0, X2, X0est, X0err, X0NoTrk, P0, P0NoTrk, P2, S; // Xp are the plane coordinates, X are the tracker coordinates (X2 = (Xp1 + Xp2) / 2)
   XYZVector   projectToHullX0, projectToHullX2;
   Float_t     wepl, wet, Lambda0, Lambda1;
   ifstream    in;
   TSpline3  * splineMCx, * splineMCy, * splineMLPx, * splineMLPy;
   TSpline3  * splineMLPNoTrkx, * splineMLPNoTrky, * splineMLPestx, * splineMLPesty;
   Double_t    arSplineMCx[1000], arSplineMCy[1000], arSplineMCz[1000], arSplineMLPx[1000], arSplineMLPy[1000], arSplineMLPz[1000];
   Double_t    arSplineMLPNoTrkx[1000], arSplineMLPNoTrky[1000], arSplineMLPestx[1000], arSplineMLPesty[1000];
   Int_t       idxSplineMC = 0, idxSplineMLP = 0;
   Double_t    differenceArrayDiff[1000] = {};
   Double_t    differenceArrayDiffNoTrk[1000] = {};
   Double_t    differenceArrayDiffest[1000] = {};
   Int_t       idxWater = 0;
   Int_t       d_entry = 15;
   Int_t       d_exit = 15;
   Double_t    energiesWater[500], rangesWater[500];
   Float_t     energy, range;
   Float_t     sigmaFilter = 0.0;
   Double_t    aPosMCx[10000], aPosMCz[10000];
   Double_t    aPosMLPx[10000], aPosMLPz[10000];
   Double_t    aPosMLPNoTrkx[10000];
   Double_t    aPosMLPestx[10000];
   Int_t       volumeID[10];
   TBranch   * b_volumeID;
   Int_t       aIdxMC = 0;
   Int_t       aIdxMLP = 0;
   Bool_t      stop = false, stopacc = false;
   Double_t    diff_x, diff_y;

   gStyle->SetTitleFont(22);
   gStyle->SetLabelFont(22);
   gStyle->SetTextFont(22);
   gStyle->SetLabelFont(22, "Y");
   gStyle->SetTitleFont(22, "Y");
   gStyle->SetTitleYOffset(1);
   gStyle->SetLabelSize(0.045);
   gStyle->SetLabelSize(0.045, "Y");
   gStyle->SetTitleSize(0.045);
   gStyle->SetTitleSize(0.045, "Y");
   gStyle->SetTextSize(0.045);

   TFile *f;

   if (spotSize <0) {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%03.0fmm.root", initialEnergy, phantomSize)); 
   }
   else {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%03.0fmm_spotsize%04.1fmm.root", initialEnergy, phantomSize, spotSize)); 
   }

   TTree *tree = (TTree*) f->Get("Hits");

   if (!tree) exit(0);

   Float_t  AXlow = 0.2;
   Float_t  AXhigh = 1.05;
   Float_t  APlow = -0.7;
   Float_t  APhigh = 0;
   
   Float_t  AXdelta = 0.005;
   Float_t  APdelta = 0.005;

   Int_t    AXbins = (AXhigh - AXlow) / AXdelta;
   Int_t    APbins = (APhigh - APlow) / APdelta;

   printf("Using %d AXbins and %d APbins.\n", AXbins, APbins);

   TH2F * hErrorMatrix = new TH2F("hErrorMatrix", Form("AUC of Errors between MC and MLP, %.0f mm B100 phantom;A_{X} parameter;A_{P} parameter", phantomSize), AXbins, AXlow, AXhigh, APbins, APlow, APhigh);
   TH2I * hIdxMatrix = new TH2I("hIdxMatrix", "Normalization matrix;A_{X} parameter;A_{P} parameter", AXbins, AXlow, AXhigh, APbins, APlow, APhigh);
   TH1I * hResidualEnergy = new TH1I("residualEnergy", "Residual Energy", 300, 0, 240);

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("volumeID", volumeID, &b_volumeID);
   
   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);

      if (eventID > eventsToUse) stop = true;

      if (lastEID < 0) {
         lastEID = eventID;
         sum_edep = edep;
         residualEnergy = 0;
      }

      if (lastEID != eventID) {
         if (lastEID % 1000 == 0) printf("%.0f mm phantom: Particle %d/%d\n", phantomSize, lastEID, eventsToUse);
         
         if (hResidualEnergy->GetEntries() > 25)   sigmaFilter = hResidualEnergy->GetMean() * 0.9;
         else                                      sigmaFilter = 0;

         if (residualEnergy > sigmaFilter) {
            hResidualEnergy->Fill(residualEnergy);
            // Find vectors as defined in paper
            X0 = (Xp0 + Xp1) / 2;
            X2 = (Xp2 + Xp3) / 2;
            
            P0 = Xp1 - Xp0;
            P2 = Xp3 - Xp2;

            P0 /= (Xp1.Z() - Xp0.Z());
            P2 /= (Xp3.Z() - Xp2.Z());

            if (std::isnan(P2.Z())) continue; // Don't ask why this can happen ...

            projectToHullX0.SetCoordinates(0, 0, d_entry);
            projectToHullX2.SetCoordinates(d_exit  * P2.X(), d_exit  * P2.Y(), d_exit );

            X0 += projectToHullX0;
            X2 -= projectToHullX2;
            
            // INNER MINIMIZATION LOOP
            for (Float_t AX = AXlow; AX <= AXhigh; AX += AXdelta) {
               for (Float_t AP = APlow; AP <= APhigh; AP += APdelta) {
                  X0est = X2 * AX + phantomSize * AP * P2;
                  diff_x = fabs(X0.X() - X0est.X());
                  diff_y = fabs(X0.Y() - X0est.Y());
                  hErrorMatrix->Fill(AX, AP, sqrt(pow(diff_x, 2) + pow(diff_y, 2)));
                  hIdxMatrix->Fill(AX, AP);
               }
            }
         }
         
         if (stop) break;

         // Reset counters for next primary
         sum_edep = edep;
         residualEnergy = 0;
         lastEID = eventID;
         idxSplineMC = 0;
      }

      else { // Still following the same particle
         if (parentID == 0) { sum_edep += edep; }
         lastEID = eventID;
      }

      if (parentID == 0) {
         if       (volumeID[2] == 0) Xp0.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 1) Xp1.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 2) Xp2.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 3) Xp3.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 5) residualEnergy += edep;
      }
   }
   cout << endl;

   gStyle->SetOptStat(0);
   gStyle->SetNumberContours(99);

   hErrorMatrix->Divide(hIdxMatrix);
  
//   Int_t nbins = hErrorMatrix->GetNcells();
   Int_t nbins = hErrorMatrix->GetNbinsX() * hErrorMatrix->GetNbinsY();
   Float_t mincont = 1e5;
   Int_t mincell = 0;
   printf("There are %d cells in the matrix\n", nbins);
   Float_t cont;
   for (Int_t i=0; i<nbins; i++) {
      cont = hErrorMatrix->GetBinContent(i);
      if (cont > 0) {
         if (cont < mincont) {
            mincont = cont;
            mincell = i;
         }
      }
   }
   
   Int_t binx, biny, binz;
   hErrorMatrix->GetBinXYZ(mincell, binx, biny, binz);
   Float_t minXvalue = hErrorMatrix->GetXaxis()->GetBinCenter(binx);
   Float_t minYvalue = hErrorMatrix->GetYaxis()->GetBinCenter(biny);
   Float_t minZvalue = hErrorMatrix->GetZaxis()->GetBinCenter(binz);
   printf("The bin with the minimum nonzero content is %.2f. AX = %.3f and AP = %.3f.\n", mincont, minXvalue, minYvalue);

   TCanvas *c2 = new TCanvas();
   hErrorMatrix->Draw("COLZ");

//   c2->SaveAs(Form("Output/accuracy_energy%.0fMeV_%.0fmm_A150.pdf", initialEnergy, phantomSize));

   // Phantom size, error , AX , AP
   ofstream file(Form("Output/accuracy_energy%.0fMeV_Water_phantom.csv", initialEnergy), ofstream::out | ofstream::app); 
   file << phantomSize << " " << mincont << " " << minXvalue << " " << minYvalue << " " <<  hResidualEnergy->GetMean() << endl;
   file.close();
}
