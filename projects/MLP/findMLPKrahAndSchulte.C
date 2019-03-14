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

#define USEHIGHENERGY
// #define USEVERYHIGHENERGY

#ifdef  USEHIGHENERGY
#define azero  5.77619e-6
#define aone   2.19784e-7
#define atwo   (-1.23920e-8)
#define athree 3.41725e-9
#define afour  (-2.20283e-10)
#define afive  5.68267e-12
#elif defined USEVERYHIGHENERGY
#define azero 3.1102e-6
#define aone 1.22439e-7
#define atwo (-9.9263e-9)
#define athree 8.25145e-10
#define afour (-2.61843e-11)
#define afive 3.29311e-13
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

double Sigmat1(double);
double Sigmaz1(double);
double Sigmatz1(double);
double Sigmat2(double, double);
double Sigmaz2(double, double);
double Sigmatz2(double, double);

XYZVector SplineMLP(Double_t t, XYZVector X0, XYZVector X2, XYZVector P0, XYZVector P2, Double_t Lambda0, Double_t Lambda2) {
   XYZVector P0Lambda, P2Lambda, X2mX0, S;
   Double_t tt = pow(t, 2), ttt = pow(t, 3);

   X2mX0 = X2 - X0;
   P0Lambda = P0 * Lambda0 * sqrt(X2mX0.Mag2());
   P2Lambda = P2 * Lambda2 * sqrt(X2mX0.Mag2());

   S = (2*ttt - 3*tt + 1) * X0 + (ttt - 2*tt + t) * P0Lambda + (-2*ttt + 3*tt) * X2 + (ttt - tt) * P2Lambda;

   return S;
}

void findMLP(Double_t phantomSize = 200, Double_t rotation = -1, Double_t spotsize = -1, Double_t initialEnergy = 230, Int_t material = kWater, Int_t kModelType = kUiB, Bool_t kModelUncertainties = false) {
   printf("---------\nRunning with phantomSize = %.0f, rotation = %.0f, spotsize = %.1f, energy = %.0f and material =%d\n-----------\n", phantomSize, rotation, spotsize, initialEnergy, material);
   Int_t      printed = 0;
   const Int_t eventsToUse = 100;
   Double_t     sigmaTheta = 0; // scattering between last two tracker layers, rad
   Double_t     sigmaPos = 0; // Position unc. in tracking layers, mm
   Bool_t      kDeleteGraphics = false, firstPass; // batch mode
   Float_t     x, y, z, edep, sum_edep = 0, residualEnergy = 0, AP, AX;
   Int_t       eventID, parentID, lastEID = -1;
   XYZVector   Xp0, Xp1, Xp2, Xp2prime, Xp3, Xp3prime, X0, X2, X0tps, X0tpsOrig, P0, P0tps, P2, S, P2gaus, P2prime, scat, scatXp2, scatXp3; // Xp are the plane coordinates, X are the tracker coordinates (X2 = (Xp1 + Xp2) / 2)
   XYZVector   X0errScuhlte, X0errKrah, X2cm, X0cm;
   Double_t     theta_x_gaus, theta_y_gaus;
   Double_t     wepl, wet, Lambda0, Lambda2;
   Double_t     theta_x_0, theta_y_0, theta_x_2, theta_y_2, tps_theta_x_0, tps_theta_y_0;
   Double_t     X_mlp_start, Y_mlp_start, theta_X_mlp_start, theta_Y_mlp_start;
   Double_t     theta_X_mlp, theta_Y_mlp, X_mlp, Y_mlp, X_mlp_sigma, Y_mlp_sigma;
   ifstream    in;
   TSpline3  * splineMCx, * splineMCy, * splineMLPx, * splineMLPy;
   TSpline3  * splineMLPKrahx, * splineMLPKrahy;
   TSpline3  * splineMLPSchultex, * splineMLPSchultey; 
   Double_t    arSplineMCx[1000], arSplineMCy[1000], arSplineMCz[1000], arSplineMLPx[1000], arSplineMLPy[1000], arSplineMLPz[1000];
   Double_t    arSplineMLPSchultex[1000], arSplineMLPSchultey[1000], arSplineMLPScuhltez[1000];
   Double_t    arSplineMLPKrahx[1000], arSplineMLPKrahy[1000], arSplineMLPKrahz[1000];
   Double_t    arMLPSchultex[1000], arMLPSchultey[1000], arMLPSchultez[1000]; 
   Double_t    arMLPKrahx[1000], arMLPKrahy[1000], arMLPKrahz[1000];
   Int_t       idxSplineMC = 0, idxSplineMLP = 0, idxMLPKrah, idxMLPSchulte;
   Int_t       aIdxMLPKrah = 0, aIdxMLPSchulte = 0;
   Double_t    differenceArrayZ[10000];
   Double_t    differenceArrayDiff[10000] = {};
   Double_t    differenceArrayDiffSchulte[10000] = {};
   Double_t    differenceArrayDiffKrah[10000] = {};
   Double_t    differenceArrayKrahSchulte[10000] = {};
   Double_t     sourceToX0dist = 100;
   Double_t     d_entry = 10;
   Double_t     d_exit = 10;
   Double_t     d_T = 10; // mm between trackers
   Double_t     d_Tcm = d_T/10; // mm between trackers
   Double_t     d_entry_cm = 1;
   Double_t     d_exit_cm = 1;
   Int_t       nInDifferenceArray = 0;
   Int_t       idxDifferenceArray = 0;
   Int_t       idxWater = 0;
   Double_t    energiesWater[500], rangesWater[500];
   Double_t     energy, range;
   Double_t     energyFilter = 0; // 128.47;
   Double_t     sigmaFilter = 1e5;
   Double_t    aPosMCx[10000], aPosMCy[10000], aPosMCz[10000];
   Double_t    aPosMLPx[10000], aPosMLPy[10000], aPosMLPz[10000];
   Double_t    aPosMLPSchultex[10000], aPosMLPSchultey[10000], aPosMLPSchultez[10000];
   Double_t    aPosMLPKrahx[10000], aPosMLPKrahy[10000], aPosMLPKrahz[10000];
   Int_t       volumeID[10];
   TBranch   * b_volumeID;
   Int_t       aIdxMC = 0;
   Int_t       aIdxMLP = 0;
   Bool_t      stop = false; 
   TRandom3  * gRandom = new TRandom3(0);
   Double_t     f10xAvg = 0;
   Double_t     f10yAvg = 0;
   Double_t     spotSizeAtX0 = 0;
   Int_t       f10xN = 0;
   Int_t       f10yN = 0;
   Int_t       maxAcc = 10;
   Int_t       firstMCidx = 0;
   char      * sMaterial;
   
   double st1, sz1, stz1, st2, sz2, stz2;
   double determinant_1, determinant_2, determinant_C12;

   // MLP matrices
   TVector3 p;
   int a; // iterator for matrix operations
   double scatter_1[4] = {0};
   double scatter_2[4] = {0};
   double sigma_beam[4] = {};
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
   f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_%s_phantom%03.0fmm.root", initialEnergy, sMaterial, phantomSize));
   
   if (spotsize < 0) spotSizeAtX0 = 3.14; // undefined, 3 mm at source
   else              spotSizeAtX0 = 0.9869 * spotsize + 0.1985; // correct for scattering in air + divergence from source to X0 plane

   Double_t spos = pow(spotSizeAtX0, -2);
   Double_t sposCM = pow(spotSizeAtX0/10, -2);

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
   TH2F *hErrorSchulte = new TH2F("hErrorSchulte", "Beamspot uncertainty (X_{0}^{Scuhlte});X position [mm];Y position [mm]", 100, -25, 25, 100, -25, 25);
   TH2F *hErrorKrah = new TH2F("hErrorKrah", "Beamspot uncertainty (X_{0}^{MLP});X position [mm];Y position [mm]", 100, -25, 25, 100, -25, 25);

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
         scatXp2.SetCoordinates(gRandom->Gaus(0, sigmaPos), gRandom->Gaus(0, sigmaPos), 0);
         scatXp3.SetCoordinates(gRandom->Gaus(0, sigmaPos), gRandom->Gaus(0, sigmaPos), 0);
         P2 = (Xp3 - Xp2) / (Xp3.Z() - Xp2.Z());
         theta_x_gaus = gRandom->Gaus(tan(P2.X()), sigmaTheta);
         theta_y_gaus = gRandom->Gaus(tan(P2.Y()), sigmaTheta);
         P2gaus.SetCoordinates(atan(theta_x_gaus), atan(theta_y_gaus), 1);
         Xp2prime = Xp2 + scatXp2;
         Xp3prime = Xp2 + P2gaus * d_T + scatXp3;

         // Redefine X2, P2 to account for different scattering uncertainties
         X0 = Xp1;
         X2 = Xp2prime;
         P0 = (Xp1 - Xp0) / (Xp1.Z() - Xp0.Z());
         P2 = (Xp3prime - Xp2prime) / (Xp3prime.Z() - Xp2prime.Z());
         
         P2prime = P2 - P0tps;
         P2prime.SetZ(1);
         
         XYZVector projectToHullX0(d_entry * tan(P0.X()), d_entry * tan(P0.Y()), d_entry);
         XYZVector projectToHullX2(d_exit * tan(P2.X()), d_exit * tan(P2.Y()), d_exit);
         X0 += projectToHullX0;
         X2 -= projectToHullX2;
         
         X0cm = X0tps / 10;
         X2cm = X2 / 10;

         wepl = splineWater->Eval(initialEnergy);
         wet = wepl - splineWater->Eval(residualEnergy);
         Double_t w = wet / wepl;
         Double_t w2 = pow(wet / wepl, 2);

         Double_t dxy = sqrt(pow(P2prime.X(), 2) + pow(P2prime.Y(), 2));
         Double_t angle = fabs(atan2(dxy,1)) * 1000;

         sigmaFilter = 21 + 280 * w2 - 515 * pow(w2,2) + 410 * pow(w2,3);
         
         if (hResEnergy->GetEntries() > 5) {
            energyFilter = hResEnergy->GetMean() * 0.9;
         }
         else energyFilter = 0;

         if (residualEnergy > energyFilter && angle < sigmaFilter) {
            hResEnergy->Fill(residualEnergy);
            double d_source = 200; // assume to first tracker layer -- is this valid for all phantom sizes / spot sizes
            double s_pos = pow(spotSizeAtX0/10, 2); // cm
            double s_angle = pow(0.002, 2); // div. so 2 mrad
            
            sigma_beam[0] = s_pos;
            sigma_beam[1] = s_pos/d_source;
            sigma_beam[2] = s_pos/d_source;
            sigma_beam[3] = s_pos/(pow(d_source,2)) + s_angle;

            double s_pos_out = sigmaPos / 10; // in cm LLU 66 um, Bergen 5 um
            double s_scat_out = sigmaTheta; // rad LLU 10 mrad, Bergen 7 mrad

            if (!kModelUncertainties) {
               s_pos_out = 0;
               s_scat_out = 0;   
            }

            T_out[0] = 0;
            T_out[1] = 0;
            T_out[2] = -1/d_Tcm;
            T_out[3] = 1/d_Tcm;

            T_out_transpose[0] = T_out[0];
            T_out_transpose[1] = T_out[2];
            T_out_transpose[2] = T_out[1];
            T_out_transpose[3] = T_out[3];

            sigma_out[0] = pow(s_pos_out, 2) * ((T_out_transpose[0] * T_out[0]) + (T_out_transpose[1]*T_out[2]));
            sigma_out[1] = pow(s_pos_out, 2) * ((T_out_transpose[0] * T_out[1]) + (T_out_transpose[1]*T_out[3]));
            sigma_out[2] = pow(s_pos_out, 2) * ((T_out_transpose[2] * T_out[0]) + (T_out_transpose[3]*T_out[2]));
            sigma_out[3] = pow(s_pos_out, 2) * ((T_out_transpose[2] * T_out[1]) + (T_out_transpose[3]*T_out[3])) + pow(s_scat_out, 2);

            S_out_inverse[0] = 1;
            S_out_inverse[1] = -d_exit_cm;
            S_out_inverse[2] = 0;
            S_out_inverse[3] = 1;

            S_out_inverse_transpose[0] = 1;
            S_out_inverse_transpose[1] = 0;
            S_out_inverse_transpose[2] = -d_exit_cm;
            S_out_inverse_transpose[3] = 1;

            track_uncert_1[0] = (sigma_out[0] * S_out_inverse_transpose[0]) + (sigma_out[1] * S_out_inverse_transpose[2]);
            track_uncert_1[1] = (sigma_out[0] * S_out_inverse_transpose[1]) + (sigma_out[1] * S_out_inverse_transpose[3]);
            track_uncert_1[2] = (sigma_out[2] * S_out_inverse_transpose[0]) + (sigma_out[3] * S_out_inverse_transpose[2]);
            track_uncert_1[3] = (sigma_out[2] * S_out_inverse_transpose[1]) + (sigma_out[3] * S_out_inverse_transpose[3]);

            track_uncert[0] = (S_out_inverse[0] * track_uncert_1[0]) + (S_out_inverse[1] * track_uncert_1[2]);
            track_uncert[1] = (S_out_inverse[0] * track_uncert_1[1]) + (S_out_inverse[1] * track_uncert_1[3]);
            track_uncert[2] = (S_out_inverse[2] * track_uncert_1[0]) + (S_out_inverse[3] * track_uncert_1[2]);
            track_uncert[3] = (S_out_inverse[2] * track_uncert_1[1]) + (S_out_inverse[3] * track_uncert_1[3]);
         
            double step_length = (X2cm.Z() - X0cm.Z()) / 512;
            firstPass = true;
            for (double posz = X0cm.Z() + step_length; posz < X2cm.Z(); posz += step_length) { 
               sz1 = Sigmaz1(posz - X0cm.Z());
               stz1 = Sigmatz1(posz - X0cm.Z());
               st1 = Sigmat1(posz - X0cm.Z());

               sz2 = Sigmaz2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());
               stz2 = Sigmatz2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());
               st2 = Sigmat2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());

               R_0[0] = 1;
               R_0[1] = posz - X0cm.Z();
               R_0[2] = 0;
               R_0[3] = 1;

               R_0_transpose[0] = 1;
               R_0_transpose[1] = 0;
               R_0_transpose[2] = posz - X0cm.Z();
               R_0_transpose[3] = 1;

               R_1_inverse[0] = 1;
               R_1_inverse[1] = -(X2cm.Z() - posz); // ok
               R_1_inverse[2] = 0;
               R_1_inverse[3] = 1;

               R_1_inverse_transpose[0] = 1;
               R_1_inverse_transpose[1] = 0;
               R_1_inverse_transpose[2] = -(X2cm.Z() - posz);
               R_1_inverse_transpose[3] = 1;
               
               scatter_1[0] = sz1;
               scatter_1[1] = stz1;
               scatter_1[2] = stz1;
               scatter_1[3] = st1;
         
               scatter_2[0] = sz2;
               scatter_2[1] = stz2;
               scatter_2[2] = stz2;
               scatter_2[3] = st2;

               // pre-factors C1 + C2 as in Krah et al. (2018)
               C1_1[0] = (sigma_beam[0] * R_0_transpose[0]) + (sigma_beam[1] * R_0_transpose[2]);
               C1_1[1] = (sigma_beam[0] * R_0_transpose[1]) + (sigma_beam[1] * R_0_transpose[3]);
               C1_1[2] = (sigma_beam[2] * R_0_transpose[0]) + (sigma_beam[3] * R_0_transpose[2]);
               C1_1[3] = (sigma_beam[2] * R_0_transpose[1]) + (sigma_beam[3] * R_0_transpose[3]);

               C1_2[0] = (R_0[0] * C1_1[0]) + (R_0[1] * C1_1[2]);
               C1_2[1] = (R_0[0] * C1_1[1]) + (R_0[1] * C1_1[3]);
               C1_2[2] = (R_0[2] * C1_1[0]) + (R_0[3] * C1_1[2]);
               C1_2[3] = (R_0[2] * C1_1[1]) + (R_0[3] * C1_1[3]);

               for (a=0; a<4; a++) C1[a] = C1_2[a] + scatter_1[a];

               C2_1[0] = (track_uncert[0] * R_1_inverse_transpose[0]) + (track_uncert[1] * R_1_inverse_transpose[2]);
               C2_1[1] = (track_uncert[0] * R_1_inverse_transpose[1]) + (track_uncert[1] * R_1_inverse_transpose[3]);
               C2_1[2] = (track_uncert[2] * R_1_inverse_transpose[0]) + (track_uncert[3] * R_1_inverse_transpose[2]);
               C2_1[3] = (track_uncert[2] * R_1_inverse_transpose[1]) + (track_uncert[3] * R_1_inverse_transpose[3]);

               C2_2[0] = (scatter_2[0] * R_1_inverse_transpose[0]) + (scatter_2[1] * R_1_inverse_transpose[2]);
               C2_2[1] = (scatter_2[0] * R_1_inverse_transpose[1]) + (scatter_2[1] * R_1_inverse_transpose[3]);
               C2_2[2] = (scatter_2[2] * R_1_inverse_transpose[0]) + (scatter_2[3] * R_1_inverse_transpose[2]);
               C2_2[3] = (scatter_2[2] * R_1_inverse_transpose[1]) + (scatter_2[3] * R_1_inverse_transpose[3]);
           
               C2[0] = (R_1_inverse[0] * C2_1[0]) + (R_1_inverse[1] * C2_1[2]) + (R_1_inverse[0] * C2_2[0]) + (R_1_inverse[1] * C2_2[2]);
               C2[1] = (R_1_inverse[0] * C2_1[1]) + (R_1_inverse[1] * C2_1[3]) + (R_1_inverse[0] * C2_2[1]) + (R_1_inverse[1] * C2_2[3]);
               C2[2] = (R_1_inverse[2] * C2_1[0]) + (R_1_inverse[3] * C2_1[2]) + (R_1_inverse[2] * C2_2[0]) + (R_1_inverse[3] * C2_2[2]);
               C2[3] = (R_1_inverse[2] * C2_1[1]) + (R_1_inverse[3] * C2_1[3]) + (R_1_inverse[2] * C2_2[1]) + (R_1_inverse[3] * C2_2[3]);

               for (a=0; a<4; a++) C12[a] = C1[a] + C2[a];

               // invert to get the complete prefactor (C1 + C2)^-1
               determinant_C12 = (C12[0] * C12[3]) - (C12[1] * C12[2]);
               C12_inverse[0] =  C12[3] / determinant_C12;
               C12_inverse[1] = -C12[1] / determinant_C12;
               C12_inverse[2] = -C12[2] / determinant_C12;
               C12_inverse[3] =  C12[0] / determinant_C12;

               first_first[0] = (C2[0] * C12_inverse[0]) + (C2[1] * C12_inverse[2]);
               first_first[1] = (C2[0] * C12_inverse[1]) + (C2[1] * C12_inverse[3]);
               first_first[2] = (C2[2] * C12_inverse[0]) + (C2[3] * C12_inverse[2]);
               first_first[3] = (C2[2] * C12_inverse[1]) + (C2[3] * C12_inverse[3]);

               second_first[0] = (C1[0] * C12_inverse[0]) + (C1[1] * C12_inverse[2]);
               second_first[1] = (C1[0] * C12_inverse[1]) + (C1[1] * C12_inverse[3]);
               second_first[2] = (C1[2] * C12_inverse[0]) + (C1[3] * C12_inverse[2]);
               second_first[3] = (C1[2] * C12_inverse[1]) + (C1[3] * C12_inverse[3]);

               y_0[0] = X0cm.X();
               y_0[1] = tan(P0tps.X());
               y_2[0] = X2cm.X();
               y_2[1] = tan(P2.X());
      
               first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
               first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);
               second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
               second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

               first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
               first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);
               second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
               second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

               X_mlp = first[0] + second[0];
               theta_X_mlp = (first[1] + second[1]);
               X_mlp_sigma = C12_inverse[2] * C2[1] + C12_inverse[3] * C2[3];

               // now do the y value
               y_0[0] = X0cm.Y();
               y_0[1] = tan(P0tps.Y());
               y_2[0] = X2cm.Y();
               y_2[1] = tan(P2.Y());
               
               first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
               first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);

               second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
               second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

               first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
               first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);

               second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
               second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

               Y_mlp = first[0] + second[0];
               theta_Y_mlp = first[1] + second[1];
               Y_mlp_sigma = C12_inverse[2] * C2[1] + C12_inverse[3] * C2[3];
               
               arMLPKrahz[idxMLPKrah] = posz*10;
               arMLPKrahx[idxMLPKrah] = X_mlp*10;
               arMLPKrahy[idxMLPKrah++] = Y_mlp*10;

               if (lastEID < maxAcc) {
                  aPosMLPKrahz[aIdxMLPKrah] = posz*10;
                  aPosMLPKrahx[aIdxMLPKrah] = X_mlp*10;
                  aPosMLPKrahy[aIdxMLPKrah++] = Y_mlp*10;
               }
            } /// End loop over KRAH phantom

            // CALCULATED OPTIMIZED VALUES
            double posz = X0cm.Z();

            sz2 = Sigmaz2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());
            stz2 = Sigmatz2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());
            st2 = Sigmat2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());

            R_0[0] = 1;
            R_0[1] = posz - X0cm.Z();
            R_0[2] = 0;
            R_0[3] = 1;

            R_0_transpose[0] = 1;
            R_0_transpose[1] = 0;
            R_0_transpose[2] = posz - X0cm.Z();
            R_0_transpose[3] = 1;

            R_1_inverse[0] = 1;
            R_1_inverse[1] = -(X2cm.Z() - posz); // ok
            R_1_inverse[2] = 0;
            R_1_inverse[3] = 1;

            R_1_inverse_transpose[0] = 1;
            R_1_inverse_transpose[1] = 0;
            R_1_inverse_transpose[2] = -(X2cm.Z() - posz);
            R_1_inverse_transpose[3] = 1;
      
            scatter_2[0] = sz2;
            scatter_2[1] = stz2;
            scatter_2[2] = stz2;
            scatter_2[3] = st2;

            // pre-factors C1 + C2 as in Krah et al. (2018)
            C1_1[0] = (sigma_beam[0] * R_0_transpose[0]) + (sigma_beam[1] * R_0_transpose[2]);
            C1_1[1] = (sigma_beam[0] * R_0_transpose[1]) + (sigma_beam[1] * R_0_transpose[3]);
            C1_1[2] = (sigma_beam[2] * R_0_transpose[0]) + (sigma_beam[3] * R_0_transpose[2]);
            C1_1[3] = (sigma_beam[2] * R_0_transpose[1]) + (sigma_beam[3] * R_0_transpose[3]);

            C1_2[0] = (R_0[0] * C1_1[0]) + (R_0[1] * C1_1[2]);
            C1_2[1] = (R_0[0] * C1_1[1]) + (R_0[1] * C1_1[3]);
            C1_2[2] = (R_0[2] * C1_1[0]) + (R_0[3] * C1_1[2]);
            C1_2[3] = (R_0[2] * C1_1[1]) + (R_0[3] * C1_1[3]);

            for (a=0; a<4; a++) C1[a] = C1_2[a];

            C2_1[0] = (track_uncert[0] * R_1_inverse_transpose[0]) + (track_uncert[1] * R_1_inverse_transpose[2]);
            C2_1[1] = (track_uncert[0] * R_1_inverse_transpose[1]) + (track_uncert[1] * R_1_inverse_transpose[3]);
            C2_1[2] = (track_uncert[2] * R_1_inverse_transpose[0]) + (track_uncert[3] * R_1_inverse_transpose[2]);
            C2_1[3] = (track_uncert[2] * R_1_inverse_transpose[1]) + (track_uncert[3] * R_1_inverse_transpose[3]);

            C2_2[0] = (scatter_2[0] * R_1_inverse_transpose[0]) + (scatter_2[1] * R_1_inverse_transpose[2]);
            C2_2[1] = (scatter_2[0] * R_1_inverse_transpose[1]) + (scatter_2[1] * R_1_inverse_transpose[3]);
            C2_2[2] = (scatter_2[2] * R_1_inverse_transpose[0]) + (scatter_2[3] * R_1_inverse_transpose[2]);
            C2_2[3] = (scatter_2[2] * R_1_inverse_transpose[1]) + (scatter_2[3] * R_1_inverse_transpose[3]);
        
            C2[0] = (R_1_inverse[0] * C2_1[0]) + (R_1_inverse[1] * C2_1[2]) + (R_1_inverse[0] * C2_2[0]) + (R_1_inverse[1] * C2_2[2]);
            C2[1] = (R_1_inverse[0] * C2_1[1]) + (R_1_inverse[1] * C2_1[3]) + (R_1_inverse[0] * C2_2[1]) + (R_1_inverse[1] * C2_2[3]);
            C2[2] = (R_1_inverse[2] * C2_1[0]) + (R_1_inverse[3] * C2_1[2]) + (R_1_inverse[2] * C2_2[0]) + (R_1_inverse[3] * C2_2[2]);
            C2[3] = (R_1_inverse[2] * C2_1[1]) + (R_1_inverse[3] * C2_1[3]) + (R_1_inverse[2] * C2_2[1]) + (R_1_inverse[3] * C2_2[3]);

            for (a=0; a<4; a++) C12[a] = C1[a] + C2[a];

            // invert to get the complete prefactor (C1 + C2)^-1
            determinant_C12 = (C12[0] * C12[3]) - (C12[1] * C12[2]);
            C12_inverse[0] =  C12[3] / determinant_C12;
            C12_inverse[1] = -C12[1] / determinant_C12;
            C12_inverse[2] = -C12[2] / determinant_C12;
            C12_inverse[3] =  C12[0] / determinant_C12;

            first_first[0] = (C2[0] * C12_inverse[0]) + (C2[1] * C12_inverse[2]);
            first_first[1] = (C2[0] * C12_inverse[1]) + (C2[1] * C12_inverse[3]);
            first_first[2] = (C2[2] * C12_inverse[0]) + (C2[3] * C12_inverse[2]);
            first_first[3] = (C2[2] * C12_inverse[1]) + (C2[3] * C12_inverse[3]);

            second_first[0] = (C1[0] * C12_inverse[0]) + (C1[1] * C12_inverse[2]);
            second_first[1] = (C1[0] * C12_inverse[1]) + (C1[1] * C12_inverse[3]);
            second_first[2] = (C1[2] * C12_inverse[0]) + (C1[3] * C12_inverse[2]);
            second_first[3] = (C1[2] * C12_inverse[1]) + (C1[3] * C12_inverse[3]);

            y_0[0] = X0cm.X();
            y_0[1] = tan(P0tps.X());
            y_2[0] = X2cm.X();
            y_2[1] = tan(P2.X());
   
            first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
            first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);
            second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
            second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

            first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
            first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);
            second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
            second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

            X_mlp = first[0] + second[0];
            theta_X_mlp = (first[1] + second[1]);

            // now do the y value
            y_0[0] = X0cm.Y();
            y_0[1] = tan(P0tps.Y());
            y_2[0] = X2cm.Y();
            y_2[1] = tan(P2.Y());
            
            first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
            first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);

            second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
            second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

            first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
            first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);

            second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
            second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

            Y_mlp = first[0] + second[0];
            theta_Y_mlp = first[1] + second[1];
            Y_mlp_sigma = C12_inverse[2] * C2[1] + C12_inverse[3] * C2[3];
            
            Y_mlp_start = Y_mlp;
            X_mlp_start = X_mlp;
            theta_Y_mlp_start = theta_Y_mlp;
            theta_X_mlp_start = theta_X_mlp;

            // Now loop over Schulte's version
            for (double posz = X0cm.Z() + step_length; posz < X2cm.Z(); posz += step_length) { 
               sz1 = Sigmaz1(posz - X0cm.Z());
               stz1 = Sigmatz1(posz - X0cm.Z());
               st1 = Sigmat1(posz - X0cm.Z());

               sz2 = Sigmaz2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());
               stz2 = Sigmatz2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());
               st2 = Sigmat2(X2cm.Z() - X0cm.Z(), posz - X0cm.Z());

               R_0[0] = 1;
               R_0[1] = posz - X0cm.Z();
               R_0[2] = 0;
               R_0[3] = 1;

               R_0_transpose[0] = 1;
               R_0_transpose[1] = 0;
               R_0_transpose[2] = posz - X0cm.Z();
               R_0_transpose[3] = 1;

               R_1_inverse[0] = 1;
               R_1_inverse[1] = -(X2cm.Z() - posz); // ok
               R_1_inverse[2] = 0;
               R_1_inverse[3] = 1;

               R_1_inverse_transpose[0] = 1;
               R_1_inverse_transpose[1] = 0;
               R_1_inverse_transpose[2] = -(X2cm.Z() - posz);
               R_1_inverse_transpose[3] = 1;
               
               scatter_1[0] = sz1;
               scatter_1[1] = stz1;
               scatter_1[2] = stz1;
               scatter_1[3] = st1;
         
               scatter_2[0] = sz2;
               scatter_2[1] = stz2;
               scatter_2[2] = stz2;
               scatter_2[3] = st2;

               // pre-factors C1 + C2 as in Krah et al. (2018)
               // NOW we use Schulte 2008, so set C1 = scatter_1 and C2 = R_1^-1 scatter_2 R_1^-1^T
               for (a=0; a<4; a++) C1[a] = scatter_1[a];

               C2_2[0] = (scatter_2[0] * R_1_inverse_transpose[0]) + (scatter_2[1] * R_1_inverse_transpose[2]);
               C2_2[1] = (scatter_2[0] * R_1_inverse_transpose[1]) + (scatter_2[1] * R_1_inverse_transpose[3]);
               C2_2[2] = (scatter_2[2] * R_1_inverse_transpose[0]) + (scatter_2[3] * R_1_inverse_transpose[2]);
               C2_2[3] = (scatter_2[2] * R_1_inverse_transpose[1]) + (scatter_2[3] * R_1_inverse_transpose[3]);
           
               C2[0] = (R_1_inverse[0] * C2_2[0]) + (R_1_inverse[1] * C2_2[2]);
               C2[1] = (R_1_inverse[0] * C2_2[1]) + (R_1_inverse[1] * C2_2[3]);
               C2[2] = (R_1_inverse[2] * C2_2[0]) + (R_1_inverse[3] * C2_2[2]);
               C2[3] = (R_1_inverse[2] * C2_2[1]) + (R_1_inverse[3] * C2_2[3]);

               for (a=0; a<4; a++) C12[a] = C1[a] + C2[a];

               // invert to get the complete prefactor (C1 + C2)^-1
               determinant_C12 = (C12[0] * C12[3]) - (C12[1] * C12[2]);
               C12_inverse[0] =  C12[3] / determinant_C12;
               C12_inverse[1] = -C12[1] / determinant_C12;
               C12_inverse[2] = -C12[2] / determinant_C12;
               C12_inverse[3] =  C12[0] / determinant_C12;

               first_first[0] = (C2[0] * C12_inverse[0]) + (C2[1] * C12_inverse[2]);
               first_first[1] = (C2[0] * C12_inverse[1]) + (C2[1] * C12_inverse[3]);
               first_first[2] = (C2[2] * C12_inverse[0]) + (C2[3] * C12_inverse[2]);
               first_first[3] = (C2[2] * C12_inverse[1]) + (C2[3] * C12_inverse[3]);

               second_first[0] = (C1[0] * C12_inverse[0]) + (C1[1] * C12_inverse[2]);
               second_first[1] = (C1[0] * C12_inverse[1]) + (C1[1] * C12_inverse[3]);
               second_first[2] = (C1[2] * C12_inverse[0]) + (C1[3] * C12_inverse[2]);
               second_first[3] = (C1[2] * C12_inverse[1]) + (C1[3] * C12_inverse[3]);

               // SCHULTE - base this on the KRAH optimised value
               y_0[0] = X_mlp_start; // X0cm.X();
               y_0[1] = theta_X_mlp_start; // tan(P0tps.X());
               y_2[0] = X2cm.X();
               y_2[1] = tan(P2.X());
      
               first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
               first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);
               second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
               second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

               first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
               first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);
               second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
               second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

               X_mlp = first[0] + second[0];
               theta_X_mlp = (first[1] + second[1]);
               X_mlp_sigma = C12_inverse[2] * C2[1] + C12_inverse[3] * C2[3];

               // now do the y value
               // SCHULTE - base this on the KRAH optimised value
               y_0[0] = Y_mlp_start; // X0cm.Y();
               y_0[1] = theta_Y_mlp_start; // tan(P0tps.Y());
               y_2[0] = X2cm.Y();
               y_2[1] = tan(P2.Y());
               
               first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
               first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);

               second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
               second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

               first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
               first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);

               second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
               second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

               Y_mlp = first[0] + second[0];
               theta_Y_mlp = first[1] + second[1];
               Y_mlp_sigma = C12_inverse[2] * C2[1] + C12_inverse[3] * C2[3];
               
               arMLPSchultez[idxMLPSchulte] = posz*10;
               arMLPSchultex[idxMLPSchulte] = X_mlp*10;
               arMLPSchultey[idxMLPSchulte++] = Y_mlp*10;

               if (lastEID < maxAcc) {
                  aPosMLPSchultez[aIdxMLPSchulte] = posz*10;
                  aPosMLPSchultex[aIdxMLPSchulte] = X_mlp*10;
                  aPosMLPSchultey[aIdxMLPSchulte++] = Y_mlp*10;
               }
            } /// End loop over phantom (SCHULTE)

            XYZVector X0Krah(X_mlp*10, Y_mlp*10, 0);
            XYZVector P0Krah(atan(theta_X_mlp), atan(theta_Y_mlp), 1);

            // Find lambda values using polynomial in Fig. 4 in Fekete et al. 2015
            Lambda0 = 1.01 + 0.43 * w2;
            Lambda2 = 0.99 - 0.46 * w2;

            for (Double_t t=0; t<1; t += 0.01) {
               S = SplineMLP(t, X0, X2, P0, P2, Lambda0, Lambda2); // perfect info
               if (lastEID < maxAcc) {
                  aPosMLPx[aIdxMLP] = S.X();
                  aPosMLPy[aIdxMLP] = S.Y();
                  aPosMLPz[aIdxMLP++] = S.Z();
               }
               
               arSplineMLPx[idxSplineMLP] = S.X(); 
               arSplineMLPy[idxSplineMLP] = S.Y();
               arSplineMLPz[idxSplineMLP] = S.Z();

               S = SplineMLP(t, X0Krah, X2, P0tps, P2, Lambda0, Lambda2); // (0,0)
               arSplineMLPSchultex[idxSplineMLP] = S.X();
               arSplineMLPSchultey[idxSplineMLP] = S.Y();

               S = SplineMLP(t, X0Krah, X2, P0Krah, P2, Lambda0, Lambda2); // niels'
               arSplineMLPKrahx[idxSplineMLP] = S.X();
               arSplineMLPKrahy[idxSplineMLP++] = S.Y();
            }

            // Compare MC and MLP here
            // 1) Make splines to compare at same Z
            splineMCx  = new TSpline3("splineMCx", arSplineMCz, arSplineMCx, idxSplineMC);
            splineMCy  = new TSpline3("splineMCy", arSplineMCz, arSplineMCy, idxSplineMC);
            splineMLPx = new TSpline3("splineMLPx", arSplineMLPz, arSplineMLPx, idxSplineMLP);
            splineMLPy = new TSpline3("splineMLPy", arSplineMLPz, arSplineMLPy, idxSplineMLP);
            splineMLPSchultex = new TSpline3("splineMLPSchultex", arMLPSchultez, arMLPSchultex, idxMLPSchulte);
            splineMLPSchultey = new TSpline3("splineMLPSchultey", arMLPSchultez, arMLPSchultey, idxMLPSchulte);
            splineMLPKrahx = new TSpline3("splineMLPKrahx", arMLPKrahz, arMLPKrahx, idxMLPKrah);
            splineMLPKrahy = new TSpline3("splineMLPKrahy", arMLPKrahz, arMLPKrahy, idxMLPKrah);

            // 2) Sweep Z values and add the absolute difference to an array
            double rvalue = splineMCx->Eval(0);
            if (!isnan(rvalue)) {
               Double_t diff_x, diff_y;
               idxDifferenceArray = 0; // Sweep each time
               nInDifferenceArray++; // To average at end
               for (Double_t zSweep = 0; zSweep < phantomSize; zSweep += 0.5) { // Keep Sweep increment high to speed up!
                  differenceArrayZ[idxDifferenceArray] = zSweep;
                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPx->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPy->Eval(zSweep));
                  differenceArrayDiff[idxDifferenceArray] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));

                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPKrahx->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPKrahy->Eval(zSweep));
                  differenceArrayDiffKrah[idxDifferenceArray] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));
                  
                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPSchultex->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPSchultey->Eval(zSweep));
                  differenceArrayDiffSchulte[idxDifferenceArray++] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));
               }
            }

            delete splineMCx;
            delete splineMCy;
            delete splineMLPx;
            delete splineMLPy;
            delete splineMLPSchultex;
            delete splineMLPSchultey;
            delete splineMLPKrahx;
            delete splineMLPKrahy;
         }
         
         if (stop) break;

         // Reset counters for next primary
         sum_edep = edep;
         residualEnergy = 0;
         lastEID = eventID;
         idxSplineMLP = 0;
         idxSplineMC = 0;
         idxMLPKrah = 0;
         idxMLPSchulte = 0;
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

   TCanvas *c1 = new TCanvas("c1", "MLP estimation", 1500, 1200);
   TPad *pad2 = new TPad("pad2", "pad2", .005, 0, 0.5235, .995);
   TPad *pad3 = new TPad("pad3", "pad3", .5245, 0, 0.995, .995);


//   pad1->Divide(1, 2, 1e-5, 1e-5);
   pad2->Divide(1, 2, 1e-5, 1e-5);
   pad3->Divide(1, 2, 1e-5, 1e-5);
//   pad4->Divide(1, 2, 1e-5, 1e-5);

  // pad1->Draw(); pad2->Draw(); pad3->Draw(); pad4->Draw();
   pad2->Draw(); pad3->Draw();

   /*
   pad1->cd(1);
   TGraph *gMC = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMCy = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLP = new TGraph(aIdxMLP, aPosMLPz, aPosMLPx);
   TGraph *gMLPy = new TGraph(aIdxMLP, aPosMLPz, aPosMLPy);
   gMC->SetMarkerStyle(21);
   gMLP->SetMarkerStyle(21);
   gMC->SetMarkerSize(0.3);
   gMLP->SetMarkerSize(0.3);
   gMC->SetMarkerColor(kRed);
   gMLP->SetMarkerColor(kBlue);
   gMC->Draw("AP");
   gMLP->Draw("P");
   gMC->SetTitle("Perfect knowledge;Depth [mm];X position [mm]");

   pad1->cd(2);
   gMCy->SetMarkerStyle(21);
   gMLPy->SetMarkerStyle(21);
   gMCy->SetMarkerSize(0.3);
   gMLPy->SetMarkerSize(0.3);
   gMCy->SetMarkerColor(kRed);
   gMLPy->SetMarkerColor(kBlue);
   gMCy->Draw("AP");
   gMLPy->Draw("P");
   gMCy->SetTitle("Perfect knowledge;Depth [mm];Y position [mm]");
*/

   // 2 and 6: Schulte model
   pad2->cd(1);
   TGraph *gMCSchulte = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPSchulte = new TGraph(aIdxMLPSchulte, aPosMLPSchultez, aPosMLPSchultex);
   gMLPSchulte->SetMarkerStyle(21);
   gMLPSchulte->SetMarkerColor(kBlue);
   gMCSchulte->SetMarkerStyle(21);
   gMCSchulte->SetMarkerColor(kRed);
   gMCSchulte->SetMarkerSize(0.3);
   gMLPSchulte->SetMarkerSize(0.3);
   gMCSchulte->Draw("AP");
   gMLPSchulte->Draw("P");
   gMCSchulte->SetTitle("No Front Tracker - Opt X_{0} + Schulte;Depth [mm];X position [mm]");
   
   pad2->cd(2);
   TGraph *gMCSchultey = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPSchultey = new TGraph(aIdxMLPSchulte, aPosMLPSchultez, aPosMLPSchultey);
   gMLPSchultey->SetMarkerStyle(21);
   gMLPSchultey->SetMarkerColor(kBlue);
   gMCSchultey->SetMarkerStyle(21);
   gMCSchultey->SetMarkerColor(kRed);
   gMCSchultey->SetMarkerSize(0.3);
   gMLPSchultey->SetMarkerSize(0.3);
   gMCSchultey->Draw("AP");
   gMLPSchultey->Draw("P");
   gMCSchultey->SetTitle("No Front Tracker - Opt X_{0} + Schulte;Depth [mm];X position [mm]");

   pad3->cd(1);
   TGraph *gMCKrah = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPKrah = new TGraph(aIdxMLPKrah, aPosMLPKrahz, aPosMLPKrahx);
   gMLPKrah->SetMarkerStyle(21);
   gMLPKrah->SetMarkerColor(kBlue);
   gMCKrah->SetMarkerStyle(21);
   gMCKrah->SetMarkerColor(kRed);
   gMCKrah->SetMarkerSize(0.3);
   gMLPKrah->SetMarkerSize(0.3);
   gMCKrah->Draw("AP");
   gMLPKrah->Draw("P");
   gMCKrah->SetTitle("No Front Tracker - Full extended MLP;Depth [mm];X position [mm]");
   
   pad3->cd(2);
   TGraph *gMCKrahy = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPKrahy = new TGraph(aIdxMLPKrah, aPosMLPKrahz, aPosMLPKrahy);
   gMLPKrahy->SetMarkerStyle(21);
   gMLPKrahy->SetMarkerColor(kBlue);
   gMCKrahy->SetMarkerStyle(21);
   gMCKrahy->SetMarkerColor(kRed);
   gMCKrahy->SetMarkerSize(0.3);
   gMLPKrahy->SetMarkerSize(0.3);
   gMCKrahy->Draw("AP");
   gMLPKrahy->Draw("P");
   gMCKrahy->SetTitle("No Front Tracker - Full extended MLP;Depth [mm];Y position [mm]");

   for (Int_t i=0; i<idxDifferenceArray; i++) {
      differenceArrayDiff[i] /= nInDifferenceArray;
      differenceArrayDiffSchulte[i]  /= nInDifferenceArray;
      differenceArrayDiffKrah[i]  /= nInDifferenceArray;
      differenceArrayKrahSchulte[i] = 1000*fabs(differenceArrayDiffSchulte[i] - differenceArrayDiffKrah[i]);
   }

   TCanvas *c2 = new TCanvas("c2", "MC vs MLP difference", 600, 600);
//   TGraph *gDifferencePerfect = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiff);
   TGraph *gDifferenceSchulte = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffSchulte);
   TGraph *gDifferenceKrah = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffKrah);
   TGraph *gDifferenceKrahSchulte = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayKrahSchulte);

   gDifferenceSchulte->SetTitle(";Depth in phantom [mm];Error | X_{1}^{opt} - X_{1}^{MC} | [mm]");
//   gDifferencePerfect->SetLineColor(kRed);
//   gDifferencePerfect->SetLineWidth(2);
   gDifferenceSchulte->SetLineColor(kBlack);
   gDifferenceSchulte->SetLineWidth(2);
   gDifferenceKrah->SetLineColor(kOrange+2);
   gDifferenceKrah->SetLineWidth(2);
//   gDifferencePerfect->Draw("L");
   gDifferenceSchulte->Draw("AL");
   gDifferenceKrah->Draw("L");
   gDifferenceSchulte->GetYaxis()->SetRangeUser(0,4.5);

   TCanvas *c3 = new TCanvas();
   gDifferenceKrahSchulte->SetTitle(";Depth in phantom [mm];Difference krah vs krah+schulte [um]");
   gDifferenceKrahSchulte->Draw("AL");

   TLegend *l = new TLegend(.6, .64, .85, .84);
//   l->AddEntry(gDifferencePerfect, "Measure X_{0}", "L");
   l->AddEntry(gDifferenceSchulte, "X_{0} Krah + Schulte MLP", "L");
   l->AddEntry(gDifferenceKrah, "Krah MLP", "L");
   l->Draw();

   char *sIsModel = (char*) "";
   if (!kModelUncertainties) sIsModel = (char*) "no";

   char *sModelType = (char*) "ll";
   if (kModelType == kUiB) sModelType = (char*) "uib";

/*
   ofstream file(Form("Output/MLPerror_energy%.0fMeV_%s_Krah.csv", initialEnergy, sMaterial), ofstream::out | ofstream::app);
   file << phantomSize << " " <<  gDifferenceSchulte->Eval(0) << " " << gDifferenceSchulte->Eval(phantomSize/2) << " " << Differenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << gDifferenceKrah->Eval(0) << " " << gDifferenceKrah->Eval(phantomSize/2) << " " << hResEnergy->GetMean() << endl;
   file.close();
 */

   if (kDeleteGraphics) {
      delete c1;
      delete c2;
      delete gMLPSchulte;
      delete gMLPSchultey;
      delete gMLPKrah;
      delete gMLPKrahy;
      delete gMCSchulte;
      delete gMCSchultey;
      delete gMCKrah;
      delete gMCKrahy;
      delete l;
      delete gDifferenceSchulte;
      delete gDifferenceKrah;
   }
}


double Sigmat1(double position)
{
	double p = position;
	double sigt1 = (azero*p)+(aone*p*p/2)+(atwo*p*p*p/3)+(athree*p*p*p*p/4)+(afour*p*p*p*p*p/5)+(afive*p*p*p*p*p*p/6);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigt1/X_0);
}

double Sigmaz1(double position)
{
	double p = position;
	double sigz1 = (azero*p*p*p/3)+(aone*p*p*p*p/12)+(atwo*p*p*p*p*p/30)+(athree*p*p*p*p*p*p/60)+(afour*p*p*p*p*p*p*p/105)+(afive*p*p*p*p*p*p*p*p/168);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigz1/X_0);
}

double Sigmatz1(double position)
{
	double p = position;
	double sigtz1 = (azero*p*p/2)+(aone*p*p*p/6)+(atwo*p*p*p*p/12)+(athree*p*p*p*p*p/20)+(afour*p*p*p*p*p*p/30)+(afive*p*p*p*p*p*p*p/42);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigtz1/X_0);
}

double Sigmat2(double sep, double position)
{
	double p = position;
	double s = sep;
	double sigt2 = ((azero*s)+(aone*s*s/2)+(atwo*s*s*s/3)+(athree*s*s*s*s/4)+(afour*s*s*s*s*s/5)+(afive*s*s*s*s*s*s/6))-((azero*p)+(aone*p*p/2)+(atwo*p*p*p/3)+(athree*p*p*p*p/4)+(afour*p*p*p*p*p/5)+(afive*p*p*p*p*p*p/6));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigt2/X_0);
}

double Sigmaz2(double sep, double position)
{
	double p = position;
	double s = sep;
	double sigz2 = (azero*s*s*s/3)+(aone*s*s*s*s/12)+(atwo*s*s*s*s*s/30)+(athree*s*s*s*s*s*s/60)+(afour*s*s*s*s*s*s*s/105)+(afive*s*s*s*s*s*s*s*s/168)-((azero*s*s*p)+(((aone*s*s/2)-(azero*s))*p*p)+(((atwo*s*s/3)-(2*aone*s/3)+(azero/3))*p*p*p)+(((athree*s*s/4)-(atwo*s/2)+(aone/4))*p*p*p*p)+(((afour*s*s/5)-(2*athree*s/5)+(atwo/5))*p*p*p*p*p)+(((afive*s*s/6)-(afour*s/3)+(athree/6))*p*p*p*p*p*p)+(((afour/7)-(2*afive*s/7))*p*p*p*p*p*p*p)+(afive*p*p*p*p*p*p*p*p/8));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigz2/X_0);
}

double Sigmatz2(double sep, double position)
{
	double p = position;
	double s = sep;
	double sigtz2 = ((azero*s*s/2)+(aone*s*s*s/6)+(atwo*s*s*s*s/12)+(athree*s*s*s*s*s/20)+(afour*s*s*s*s*s*s/30)+(afive*s*s*s*s*s*s*s/42))-((azero*s*p)+(((aone*s)-azero)*p*p/2)+(((atwo*s)-aone)*p*p*p/3)+(((athree*s)-atwo)*p*p*p*p/4)+(((afour*s)-athree)*p*p*p*p*p/5)+(((afive*s)-afour)*p*p*p*p*p*p/6)-(afive*p*p*p*p*p*p*p/7));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigtz2/X_0);
}
