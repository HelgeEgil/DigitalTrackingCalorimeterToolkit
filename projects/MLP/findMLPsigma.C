#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <fstream>
#include <TBranch.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TStyle.h>
#include <Math/Vector3D.h>
#include <TRandom3.h>

#define X_0 36.1

enum eMat {kWater, kA150, kB100, kCorticalBone, kAdipose};

const int kAdiposeMax = 350;
const int kA150Max = 290;
const int kWaterMax = 330;
const int kWater200Max = 270;
const int kB100Max = 250;
const int kCorticalBoneMax = 200;

using namespace std;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > XYZVector;

XYZVector SplineMLP(Double_t t, XYZVector X0, XYZVector X2, XYZVector P0, XYZVector P2, Double_t Lambda0, Double_t Lambda2) {
   XYZVector P0Lambda, P2Lambda, X2mX0, S;
   Float_t tt = pow(t, 2), ttt = pow(t, 3);

   X2mX0 = X2 - X0;
   P0Lambda = P0 * Lambda0 * sqrt(X2mX0.Mag2());
   P2Lambda = P2 * Lambda2 * sqrt(X2mX0.Mag2());

   S = (2*ttt - 3*tt + 1) * X0 + (ttt - 2*tt + t) * P0Lambda + (-2*ttt + 3*tt) * X2 + (ttt - tt) * P2Lambda;

   return S;
}


void findMLP(Float_t phantomSize = 200, Float_t rotation = -1, Float_t spotsize = -1, Float_t initialEnergy = 230, Int_t material = kWater) {
   Int_t      printed = 0;
   const Int_t eventsToUse = 100000;
   Float_t     Xp3sigma = 0; // scattering between last two tracker layers, dz * X mrad
   Float_t     x, y, z, edep, sum_edep = 0, residualEnergy = 0, AP, AX;
   Int_t       eventID, parentID, lastEID = -1;
   XYZVector   Xp0, Xp1, Xp2, Xp3, X0, X2, X0est, X0err, X0tps, P0, P0tps, P2, S, P2prime, X2prime, scat; // Xp are the plane coordinates, X are the tracker coordinates (X2 = (Xp1 + Xp2) / 2)
   XYZVector   X0errNaive, X0errKrah, X2cm, X0cm;
   Float_t     wepl, wet, Lambda0, Lambda2;
   Float_t     theta_x_0, theta_y_0, theta_x_2, theta_y_2, tps_theta_x_0, tps_theta_y_0;
   ifstream    in;
   TSpline3  * splineMCx, * splineMCy, * splineMLPx, * splineMLPy;
   TSpline3  * splineMLPKrahx, * splineMLPKrahy, * splineMLPestx, * splineMLPesty;
   TSpline3  * splineMLPNoTrkx, * splineMLPNoTrky; 
   Double_t    arSplineMCx[1000], arSplineMCy[1000], arSplineMCz[1000], arSplineMLPx[1000], arSplineMLPy[1000], arSplineMLPz[1000];
   Double_t    arSplineMLPNoTrkx[1000], arSplineMLPNoTrky[1000], arSplineMLPestx[1000], arSplineMLPesty[1000];
   Double_t    arSplineMLPKrahx[1000], arSplineMLPKrahy[1000]; 
   Int_t       idxSplineMC = 0, idxSplineMLP = 0;
   Double_t    differenceArrayZ[1000];
   Double_t    differenceArrayDiff[1000] = {};
   Double_t    differenceArrayDiffNoTrk[1000] = {};
   Double_t    differenceArrayDiffKrah[1000] = {};
   Double_t    differenceArrayDiffest[1000] = {};
   Float_t     sourceToX0dist = 100;
   Float_t     d_entry = 15;
   Float_t     d_exit = 15;
   Int_t       nInDifferenceArray = 0;
   Int_t       idxDifferenceArray = 0;
   Int_t       idxWater = 0;
   Double_t    energiesWater[500], rangesWater[500];
   Float_t     energy, range;
   Float_t     energyFilter = 0; // 128.47;
   Float_t     sigmaFilter = 1e5;
   Double_t    aPosMCx[10000], aPosMCy[10000], aPosMCz[10000];
   Double_t    aPosMLPx[10000], aPosMLPy[10000], aPosMLPz[10000];
   Double_t    aPosMLPNoTrkx[10000], aPosMLPNoTrky[10000];
   Double_t    aPosMLPKrahx[10000], aPosMLPKrahy[10000];
   Double_t    aPosMLPestx[10000], aPosMLPesty[10000];
   Int_t       volumeID[10];
   TBranch   * b_volumeID;
   Int_t       aIdxMC = 0;
   Int_t       aIdxMLP = 0;
   Bool_t      stop = false; 
   TRandom3  * gRandom = new TRandom3(0);
   Float_t     f10xAvg = 0;
   Float_t     f10yAvg = 0;
   Int_t       f10xN = 0;
   Int_t       f10yN = 0;
   Int_t       maxAcc = 7;
   char      * sMaterial;

   if      (material == kWater)        sMaterial = (char*) "Water";
   else if (material == kB100)         sMaterial = (char*) "B100";
   else if (material == kA150)         sMaterial = (char*) "A150";
   else if (material == kAdipose)      sMaterial = (char*) "myAdipose";
   else if (material == kCorticalBone) sMaterial = (char*) "CorticalBone";
   else {
      cout << "Material " << material << " not defined!\n";
      exit(0);
   }

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

   TFile *f = nullptr;

   if (rotation >= 0) {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_%s_phantom%03.0fmm_rotation%02.0fmrad.root", initialEnergy, sMaterial, phantomSize, rotation));
   }

   else if (spotsize >= 0)  {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_%s_phantom%03.0fmm_spotsize%04.1fmm.root", initialEnergy, sMaterial, phantomSize, spotsize));
   }
   
   else {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_%s_phantom%03.0fmm.root", initialEnergy, sMaterial, phantomSize));
   }
   
   if (material == kAdipose) sMaterial = (char*) "Adipose";

   TTree *tree = (TTree*) f->Get("Hits");

   if (!tree) exit(0);
   // Load Energy <-> Range spline
   in.open("Data/WaterPSTAR.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;
      rangesWater[idxWater] = range*10; // [mm]
      energiesWater[idxWater++] = energy;
   }
   in.close();
   TSpline3 *splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);

   TH1F *hResEnergy = new TH1F("hResEnergy", "Residual energy in calorimeter;Energy [MeV];Entries", 300, 0, 250);
   TH1F *hAngle = new TH1F("hAngle", "Outgoing angle;Energy [MeV];Entries", 500, 0, 1000);

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("volumeID", volumeID, &b_volumeID);
  
   if (rotation >= 0) { // defined, use actual number (mrad rotated about X)
      tps_theta_x_0 = 0;
      tps_theta_y_0 = -rotation / 1000;
   }
   else { // rotation = -1 means not defined, not interesting, etc ... it's zero
      tps_theta_x_0 = 0;
      tps_theta_y_0 = 0;
   }
   
   P0tps.SetCoordinates(atan(tps_theta_x_0), atan(tps_theta_y_0), 1);
   X0tps.SetCoordinates(P0tps.X() * (sourceToX0dist + d_entry), P0tps.Y() * (sourceToX0dist + d_entry), 0);


   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);

      if (eventID > eventsToUse) stop = true;

      if (lastEID < 0) {
         lastEID = eventID;
         sum_edep = edep;
         residualEnergy = 0;
      }

      if (lastEID != eventID) {
         // Add scattering
         scat.SetCoordinates(gRandom->Gaus(0, Xp3sigma), gRandom->Gaus(0, Xp3sigma), 0);
         Xp3 += scat;

         X0 = (Xp0 + Xp1) / 2;
         X2 = (Xp2 + Xp3) / 2;
         
         P0 = (Xp1 - Xp0) / (Xp1.Z() - Xp0.Z());
         P2 = (Xp3 - Xp2) / (Xp3.Z() - Xp2.Z());

         P2prime = P2 - P0tps;
         P2prime.SetZ(1);

         wepl = splineWater->Eval(initialEnergy);
         wet = wepl - splineWater->Eval(residualEnergy);
         Float_t w = wet / wepl;

         Float_t dxy = sqrt(pow(P2prime.X(), 2) + pow(P2prime.Y(), 2));
         Float_t angle = fabs(atan2(dxy,1)) * 1000;

//         sigmaFilter = 139.15 * w + 13.991;
         hResEnergy->Fill(residualEnergy);
         hAngle->Fill(angle);
         
         if (stop) break;

         // Reset counters for next primary
         sum_edep = edep;
         residualEnergy = 0;
         lastEID = eventID;
         idxSplineMLP = 0;
         idxSplineMC = 0;
      }

      else { // Still following the same particle
         if (parentID == 0) { sum_edep += edep; }
         lastEID = eventID;
      }

      if (parentID == 0) {
         if       (volumeID[2] == 0) Xp0.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 1) Xp1.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 2) Xp2.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 3) Xp3.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 5) residualEnergy += edep;


         if  (volumeID[2] < 5) {
            arSplineMCx[idxSplineMC] = x;
            arSplineMCy[idxSplineMC] = y;
            arSplineMCz[idxSplineMC++] = z+phantomSize/2;
            
            if (eventID < maxAcc) {
               aPosMCx[aIdxMC] = x;
               aPosMCy[aIdxMC] = y;
               aPosMCz[aIdxMC++] = z+phantomSize/2;
            }
         }
      }
   }

   float empMu = hAngle->GetMean();
   float empSigma = hAngle->GetStdDev();

   TF1 * fit = new TF1("fit", "gaus");
   hAngle->Fit("fit");
   float fitMu = fit->GetParameter(1);
   float fitSigma = fit->GetParameter(2);

   printf("Empirical values: %.3f +- %.3f. Fit values: %.3f +- %.3f.\n", empMu, empSigma, fitMu, fitSigma);

   ofstream file(Form("Output/P2_energy%.0fMeV_%s_krah.csv", initialEnergy, sMaterial), ofstream::out | ofstream::app);
   file << phantomSize << " " << empMu << " " << empSigma << " " << fitMu << " " << fitSigma  << " " << hResEnergy->GetMean() << endl;
   file.close();
}

