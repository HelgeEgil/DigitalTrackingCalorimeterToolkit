#include <TTree.h>
#include <TFile.h>
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

// #define USEHIGHENERGY

#ifdef  USEHIGHENERGY
#define azero  5.77619e-6
#define aone   2.19784e-7
#define atwo   (-1.23920e-8)
#define athree 3.41725e-9
#define afour  (-2.20283e-10)
#define afive  5.68267e-12
#else
#define azero   7.556e-6
#define aone    4.548e-7
#define atwo    (-5.777e-8)
#define athree  1.301e-8
#define afour   (-9.228e-10)
#define afive   2.687e-11
#endif

#define X_0 36.1

enum eMat {kWater, kA150, kB100, kCorticalBone, kAdipose};
enum eMod {kUiB, kLL};


const int kAdiposeMax = 350;
const int kA150Max = 290;
const int kWaterMax = 330;
const int kWater200Max = 270;
const int kB100Max = 250;
const int kCorticalBoneMax = 200;

using namespace std;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > XYZVector;

float Sigmat1(float);
float Sigmaz1(float);
float Sigmatz1(float);
float Sigmat2(float, float);
float Sigmaz2(float, float);
float Sigmatz2(float, float);

XYZVector SplineMLP(Double_t t, XYZVector X0, XYZVector X2, XYZVector P0, XYZVector P2, Double_t Lambda0, Double_t Lambda2) {
   XYZVector P0Lambda, P2Lambda, X2mX0, S;
   Float_t tt = pow(t, 2), ttt = pow(t, 3);

   X2mX0 = X2 - X0;
   P0Lambda = P0 * Lambda0 * sqrt(X2mX0.Mag2());
   P2Lambda = P2 * Lambda2 * sqrt(X2mX0.Mag2());

   S = (2*ttt - 3*tt + 1) * X0 + (ttt - 2*tt + t) * P0Lambda + (-2*ttt + 3*tt) * X2 + (ttt - tt) * P2Lambda;

   return S;
}

void findMLP(Float_t phantomSize = 200, Float_t rotation = -1, Float_t spotsize = -1, Float_t initialEnergy = 230, Int_t material = kWater, Int_t kModelType = kUiB, Bool_t kModelUncertainties = false) {
   printf("---------\nRunning with phantomSize = %.0f, rotation = %.0f, spotsize = %.1f, energy = %.0f and material =%d\n-----------\n", phantomSize, rotation, spotsize, initialEnergy, material);
   Int_t      printed = 0;
   const Int_t eventsToUse = 100000;
   Float_t     sigmaTheta = 0; // scattering between last two tracker layers, rad
   Float_t     sigmaPos = 0; // Position unc. in tracking layers, mm
   Bool_t      kDeleteGraphics = false; // batch mode
   Float_t     x, y, z, edep, sum_edep = 0, residualEnergy = 0, AP, AX;
   Int_t       eventID, parentID, lastEID = -1;
   XYZVector   Xp0, Xp1, Xp2, Xp2prime, Xp3, Xp3prime, X0, X2, X0est, X0err, X0tps, P0, P0tps, P2, S, P2gaus, P2prime, X2prime, scat, scatXp2, scatXp3; // Xp are the plane coordinates, X are the tracker coordinates (X2 = (Xp1 + Xp2) / 2)
   XYZVector   X0errNaive, X0errKrah, X2cm, X0cm;
   Float_t     theta_x_gaus, theta_y_gaus;
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
   Float_t     d_entry = 10;
   Float_t     d_exit = 10;
   Float_t     d_T = 10; // mm between trackers
   Float_t     d_Tcm = d_T/10; // mm between trackers
   Float_t     d_entry_cm = 1;
   Float_t     d_exit_cm = 1;
   Int_t       nInDifferenceArray = 0;
   Int_t       idxDifferenceArray = 0;
   Int_t       idxWater = 0;
   Double_t    energiesWater[500], rangesWater[500];
   Float_t     energy, range;
   Float_t     energyFilter = 0; // 128.47;
   Float_t     sigmaFilter = 1e5;
   Double_t    aPosMCx[1000], aPosMCy[1000], aPosMCz[1000];
   Double_t    aPosMLPx[1000], aPosMLPy[1000], aPosMLPz[1000];
   Double_t    aPosMLPNoTrkx[1000], aPosMLPNoTrky[1000];
   Double_t    aPosMLPKrahx[1000], aPosMLPKrahy[1000];
   Double_t    aPosMLPestx[1000], aPosMLPesty[1000];
   Int_t       volumeID[10];
   TBranch   * b_volumeID;
   Int_t       aIdxMC = 0;
   Int_t       aIdxMLP = 0;
   Bool_t      stop = false; 
   TRandom3  * gRandom = new TRandom3(0);
   Float_t     f10xAvg = 0;
   Float_t     f10yAvg = 0;
   Float_t     spotSizeAtX0 = 0;
   Int_t       f10xN = 0;
   Int_t       f10yN = 0;
   Int_t       maxAcc = 5;
   Int_t       firstMCidx = 0;
   char      * sMaterial;
   
   float st1, sz1, stz1, st2, sz2, stz2;
   float determinant_1, determinant_2, determinant_C12;
   float X_mlp, Y_mlp;

   // MLP matrices
   TVector3 p;
   int a; // iterator for matrix operations
   float scatter_1[4] = {0};
   float scatter_2[4] = {0};
   float sigma_beam[4] = {};
   double R_0[4] = {0};
   double R_0_transpose[4] = {0};
   double R_1_inverse[4] = {0};
   double R_1_inverse_transpose[4] = {0};
   double y_0[2] = {0};
   double y_2[2] = {0};
   double C1_1[4] = {0};
   double C1_2[4] = {0};
   double C1[4] = {0};
   double C2_1[4] = {0};
   double C2_2[4] = {0};
   double C2[4] = {0};
   double C12[4];
   double C12_inverse[4] = {0};
   double first_first[4] = {0};
   double second_first[4] = {0};
   double first_second[2] = {0};
   double second_second[2] = {0};
   double first[2] = {0};
   double second[2] = {0};
   double T_out[4];
   double T_out_transpose[4];
   double sigma_out[4];
   double S_out_inverse[4];
   double S_out_inverse_transpose[4];
   double track_uncert_1[4];
   double track_uncert[4];

   if (kModelType == kLL) {
      sigmaPos = 0.066;
      sigmaTheta = 0.01;
   }

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
   gStyle->SetPadGridX(kTRUE);
   gStyle->SetPadGridY(kTRUE);

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

   if (spotsize < 0) spotSizeAtX0 = 3.14; // undefined, 3 mm at source
   else              spotSizeAtX0 = 0.9869 * spotsize + 0.1985; // correct for scattering in air + divergence from source to X0 plane

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

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("volumeID", volumeID, &b_volumeID);
  
   P0tps.SetCoordinates(0, 0, 1);
   X0tps.SetCoordinates(0, 0, 0);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);

      if (eventID > eventsToUse) stop = true;

      if (lastEID < 0) {
         lastEID = eventID;
         sum_edep = edep;
         residualEnergy = 0;
      }

      if (lastEID != eventID) {
         X0 = Xp1;
         X2 = Xp2;
         P0 = (Xp1 - Xp0) / (Xp1.Z() - Xp0.Z());
         P2 = (Xp3 - Xp2) / (Xp3.Z() - Xp2.Z());

         P2prime = P2 - P0tps;
         P2prime.SetZ(1);
         
         XYZVector projectToHullX0(d_entry * tan(P0.X()), d_entry * tan(P0.Y()), d_entry);
         XYZVector projectToHullX2(d_exit * tan(P2.X()), d_exit * tan(P2.Y()), d_exit);
         X0 += projectToHullX0;
         X2 -= projectToHullX2;
         
         wepl = splineWater->Eval(initialEnergy);
         wet = wepl - splineWater->Eval(residualEnergy);
         Float_t w = wet / wepl;
         Float_t w2 = pow(wet / wepl, 2);
  
         // NEW VERSION OF LPM !!!!!
         AX = exp(2.9143 - 5.5692 * w - 1.4734 * w2 + 1.2822 * w*w2);
         AP = exp(2.4946 - 8.0676 * w + 4.2200 * w2 - 3.2835 * w*w2);

         Float_t dxy = sqrt(pow(P2prime.X(), 2) + pow(P2prime.Y(), 2));
         Float_t angle = fabs(atan2(dxy,1)) * 1000;

         sigmaFilter = 21 + 280 * w2 - 515 * pow(w2,2) + 410 * pow(w2,3);
         
         if (hResEnergy->GetEntries() > 5) {
            energyFilter = hResEnergy->GetMean() * 0.9;
         }
         else energyFilter = 0;

         if (residualEnergy > energyFilter && angle < sigmaFilter) {
            hResEnergy->Fill(residualEnergy);
            
            X2prime = X2 - X0tps - phantomSize * P0tps;
            X0est = X2prime * AX/(pow(spotSizeAtX0,-2)+AX) - P2prime * AP/(pow(spotSizeAtX0,-2)+AX) * phantomSize; // LPM

            X0est += X0tps;
            X0est.SetZ(0);

            // Find lambda values using polynomial in Fig. 4 in Fekete et al. 2015
            Lambda0 = 1.01 + 0.43 * w2;
            Lambda2 = 0.99 - 0.46 * w2;

            for (Float_t t=0; t<1; t += 0.01) {
               S = SplineMLP(t, X0est, X2, P0tps, P2, Lambda0, Lambda2); // mine
               arSplineMLPestx[idxSplineMLP] = S.X();
               arSplineMLPesty[idxSplineMLP] = S.Y();
            }

            // Compare MC and MLP here
            // 1) Make splines to compare at same Z
            splineMCx  = new TSpline3("splineMCx", arSplineMCz, arSplineMCx, idxSplineMC);
            splineMCy  = new TSpline3("splineMCy", arSplineMCz, arSplineMCy, idxSplineMC);
            splineMLPestx = new TSpline3("splineMLPestx", arSplineMLPz, arSplineMLPestx, idxSplineMLP);
            splineMLPesty = new TSpline3("splineMLPesty", arSplineMLPz, arSplineMLPesty, idxSplineMLP);

            // 2) Sweep Z values and add the absolute difference to an array
            float rvalue = splineMCx->Eval(0);
            if (!isnan(rvalue)) {
               Double_t diff_x, diff_y;
               idxDifferenceArray = 0; // Sweep each time
               nInDifferenceArray++; // To average at end
               for (Double_t zSweep = 0; zSweep <= phantomSize; zSweep += 0.5) { // Keep Sweep increment high to speed up!
                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPestx->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPesty->Eval(zSweep));
                  differenceArrayZ[idxDifferenceArray] = zSweep;
                  differenceArrayDiffest[idxDifferenceArray++] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));
               }
            }

            delete splineMCx;
            delete splineMCy;
            delete splineMLPestx;
            delete splineMLPesty;
         }
         
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

         if       (volumeID[2] == 0) Xp0.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 1) Xp1.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 2) Xp2.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 3) Xp3.SetCoordinates(x,y,z+phantomSize/2);
         else if  (volumeID[2] == 5) residualEnergy += edep;
      }
   }

   TCanvas *c2 = new TCanvas("c2", "MC vs MLP difference", 600, 600);
   TGraph *gDifferenceest = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffest);
   TGraph *gDifferenceestHe = new TGraph(idxDifferenceArrayHe, differenceArrayZHe, differenceArrayDiffestHe);

   gDifferenceest->SetTitle(";Depth in phantom [mm];Error | X_{1}^{opt} - X_{1}^{MC} | [mm]");
   gDifferenceest->SetLineColor(kBlue);
   gDifferenceest->SetLineWidth(3);
   gDifferenceestHe->SetLineColor(kRed-7);
   gDifferenceestHe->SetLineWidth(3);
   gDifferenceest->Draw("LA");
   gDifferenceestHe->Draw("A");

   TLegend *l = new TLegend(.6, .64, .85, .84);
   l->AddEntry(gDifferenceest, "LPM protons", "L");
   l->AddEntry(gDifferenceestHe, "LPM Helium", "L");
   l->Draw();

   printf("The proton error is %.2f mm. The Helium error is %.2f mm.\n", gDifferenceest->Eval(0), gDifferenceestHe->Eval(0));

}
