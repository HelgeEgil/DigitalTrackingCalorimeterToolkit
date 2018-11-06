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

#ifndef  USEHIGHENERGY
#define azero   7.457e-6
#define aone    4.548e-7
#define atwo   (-5.777e-8)
#define athree  1.301e-8
#define afour  (-9.228e-10)
#define afive   2.687e-11
#else
#define azero   5.77619e-6
#define aone    2.19784e-7
#define atwo    (-1.23920e-8)
#define athree  3.41725e-9
#define afour   (-2.20283e-10)
#define afive   5.68267e-12
#endif

#define X_0 36.1

enum eMat {kWater, kA150, kB100, kCorticalBone, kAdipose};

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


void findMLP(Float_t phantomSize = 200, Float_t rotation = -1, Float_t spotsize = -1, Float_t initialEnergy = 230, Int_t material = kWater) {
   Int_t      printed = 0;
   const Int_t eventsToUse = 25000;
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
   TH2F *hErrorNaive = new TH2F("hErrorNaive", "Beamspot uncertainty (assume point beam);X position [mm];Y position [mm]", 100, -25, 25, 100, -25, 25);
   TH2F *hErrorKrah = new TH2F("hErrorKrah", "Beamspot uncertainty (assume point beam);X position [mm];Y position [mm]", 100, -25, 25, 100, -25, 25);
   TH1F *hErrorKrah1D = new TH1F("hErrorKrah1D", "Beamspot uncertainty (assume point beam);X position [mm];frequency", 100, -2, 2);
   TH2F *hErrorSigmaScale = new TH2F("hErrorSigmaScale", Form("Beamspot uncertainty (est);X position [mm];Y position [mm]"), 100, -25, 25, 100, -25, 25);

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
   
   P0tps.SetCoordinates(tps_theta_x_0, tps_theta_y_0, sqrt(1-pow(tps_theta_x_0, 2) - pow(tps_theta_y_0, 2)));
   X0tps.SetCoordinates(tan(tps_theta_x_0) * (sourceToX0dist + d_entry), tan(tps_theta_y_0) * (sourceToX0dist + d_entry), 0);


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
         
         XYZVector projectToHullX0(d_entry * tan(P0tps.X()), d_entry * tan(P0tps.Y()), d_entry);
         XYZVector projectToHullX2(d_exit * tan(P2.X()), d_exit * tan(P2.Y()), d_exit);
         X0 += projectToHullX0;
         X2 -= projectToHullX2;
         
         X0cm = X0tps / 10;
         X2cm = X2 / 10;

         wepl = splineWater->Eval(initialEnergy);
         wet = wepl - splineWater->Eval(residualEnergy);
         Float_t w = wet / wepl;
         Float_t w2 = pow(wet / wepl, 2);
            
         AX = 1 - 0.30 * w + 1.53 * w2 - 4.25 * pow(w,3) + 2.30 * pow(w,4);
         AP = -0.659 + 2.072* w - 8.66 * w2 + 18.14 * pow(w,3) - 16.27 * pow(w,4) + 5.34 * pow(w,5);
         
         if (spotsize >= 0) {
            AX =  -3.984e-2 + 5.928e-1 * spotsize - 1.436e-1 * pow(spotsize,2) + 1.699e-2 * pow(spotsize,3) - 9.634e-4 * pow(spotsize,4) + 2.093e-5 * pow(spotsize,5);
            AP = 1.789e-2 - 2.604e-1 * spotsize + 6.261e-2  * pow(spotsize,2) - 7.354e-3 * pow(spotsize,3) + 4.145e-4 * pow(spotsize,4) - 8.953e-6 * pow(spotsize,5);
         }

         if (printed++ < 5) {
            printf("w = %2f -> AX = %.2f, AP = %.2f\n", w, AX, AP);
         }

         Float_t dxy = sqrt(pow(P2prime.X(), 2) + pow(P2prime.Y(), 2));
         Float_t angle = fabs(atan2(dxy,1)) * 1000;

         sigmaFilter = 139.15 * w + 13.991;
         
         if (hResEnergy->GetEntries() > 5) {
            energyFilter = hResEnergy->GetMean() * 0.9;
         }
         else energyFilter = 0;

         if (residualEnergy > energyFilter && angle < sigmaFilter) {
            hResEnergy->Fill(residualEnergy);
         
            double step_length = (X2cm.Z() - X0cm.Z()) / 512;
            double posz = X0cm.Z() + step_length;

            float st2, sz2, stz2;
            float determinant_1, determinant_2, determinant_C12;
            float d_source = 10000; // assume to first tracker layer -- is this valid for all phantom sizes / spot sizes
            float s_pos = pow(0.32, 2); // cm
            float s_angle = pow(1e-6, 2); // div. so 2 mrad

            if (spotsize >= 0) s_pos = pow(spotsize / 10 + 0.02, 2);
      
            // init MLP
            float X_mlp, Y_mlp;
            TVector3 p;

            int a; // iterator for matrix operations
            float scatter_2[4] = {0};

            float sigma_beam[4] = {};
            sigma_beam[0] = s_pos;
            sigma_beam[1] = s_pos/d_source;
            sigma_beam[2] = s_pos/d_source;
            sigma_beam[3] = s_pos/(pow(d_source,2)) + s_angle;

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
            double C2[4] = {0};
            double C12[4];
            double C12_inverse[4] = {0};
            
            double first_first[4] = {0};
            double second_first[4] = {0};

            double first_second[2] = {0};
            double second_second[2] = {0};

            double first[2] = {0};
            double second[2] = {0};

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

            C1[0] = (R_0[0] * C1_1[0]) + (R_0[1] * C1_1[2]);
            C1[1] = (R_0[0] * C1_1[1]) + (R_0[1] * C1_1[3]);
            C1[2] = (R_0[2] * C1_1[0]) + (R_0[3] * C1_1[2]);
            C1[3] = (R_0[2] * C1_1[1]) + (R_0[3] * C1_1[3]);

            C2_1[0] = (scatter_2[0] * R_1_inverse_transpose[0]) + (scatter_2[1] * R_1_inverse_transpose[2]);
            C2_1[1] = (scatter_2[0] * R_1_inverse_transpose[1]) + (scatter_2[1] * R_1_inverse_transpose[3]);
            C2_1[2] = (scatter_2[2] * R_1_inverse_transpose[0]) + (scatter_2[3] * R_1_inverse_transpose[2]);
            C2_1[3] = (scatter_2[2] * R_1_inverse_transpose[1]) + (scatter_2[3] * R_1_inverse_transpose[3]);
         
            C2[0] = (R_1_inverse[0] * C2_1[0]) + (R_1_inverse[1] * C2_1[2]);
            C2[1] = (R_1_inverse[0] * C2_1[1]) + (R_1_inverse[1] * C2_1[3]);
            C2[2] = (R_1_inverse[2] * C2_1[0]) + (R_1_inverse[3] * C2_1[2]);
            C2[3] = (R_1_inverse[2] * C2_1[1]) + (R_1_inverse[3] * C2_1[3]);

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
            y_0[1] = atan(P0.X());
            y_2[0] = X2cm.X();
            y_2[1] = atan(P2.X());
   
            first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
            first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);
            second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
            second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

            first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
            first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);
            second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
            second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

            X_mlp = first[0] + second[0];
            double theta_X_mlp = (first[1] + second[1]);
            double X_mlp_sigma = C12_inverse[2] * C2[1] + C12_inverse[3] * C2[3];

            // now do the y value
            y_0[0] = X0cm.Y();
            y_0[1] = atan(P0.Y());
            y_2[0] = X2cm.Y();
            y_2[1] = atan(P2.Y());
            
            first_second[0] = (R_0[0] * y_0[0]) + (R_0[1] * y_0[1]);
            first_second[1] = (R_0[2] * y_0[0]) + (R_0[3] * y_0[1]);

            second_second[0] = (R_1_inverse[0] * y_2[0]) + (R_1_inverse[1] * y_2[1]);
            second_second[1] = (R_1_inverse[2] * y_2[0]) + (R_1_inverse[3] * y_2[1]);

            first[0] = (first_first[0] * first_second[0]) + (first_first[1] * first_second[1]);
            first[1] = (first_first[2] * first_second[0]) + (first_first[3] * first_second[1]);

            second[0] = (second_first[0] * second_second[0]) + (second_first[1] * second_second[1]);
            second[1] = (second_first[2] * second_second[0]) + (second_first[3] * second_second[1]);

            Y_mlp = first[0] + second[0];
            double theta_Y_mlp = first[1] + second[1];
            double Y_mlp_sigma = C12_inverse[2] * C2[1] + C12_inverse[3] * C2[3];

            XYZVector X0krah(X_mlp*10, Y_mlp*10, 0);
            XYZVector P0krah(theta_X_mlp, theta_Y_mlp, sqrt(1 - pow(theta_X_mlp, 2) - pow(theta_Y_mlp, 2)));

            X2prime = X2 - X0tps - phantomSize * P0tps;

            X0est = X2prime * AX + P2prime * AP * phantomSize;
            X0est += X0tps + phantomSize * P0tps;
            X0est.SetZ(0);

            X0err = X0est - X0;
            X0errNaive = X0tps - X0;
            X0errKrah = X0krah - X0;

            hErrorSigmaScale->Fill(X0err.X(), X0err.Y());
            hErrorNaive->Fill(X0errNaive.X(), X0errNaive.Y());
            hErrorKrah->Fill(X0errKrah.X(), X0errKrah.Y());
            hErrorKrah1D->Fill(X0errKrah.X() / 10 ) ;

            float f10x = sqrt(2 * log(10)) / (2 * 3.141592 * X_mlp_sigma);
            float f10y = sqrt(2 * log(10)) / (2 * 3.141592 * Y_mlp_sigma);

            f10xAvg += f10x;
            f10yAvg += f10y;
            f10xN++;
            f10yN++;

            // Find lambda values using polynomial in Fig. 4 in Fekete et al. 2015
            Lambda0 = 1.01 + 0.43 * w2;
            Lambda2 = 0.99 - 0.46 * w2;

            for (Float_t t=0; t<1; t += 0.01) {
               S = SplineMLP(t, X0, X2, P0, P2, Lambda0, Lambda2); // perfect info
               if (lastEID < maxAcc) {
                  aPosMLPx[aIdxMLP] = S.X();
                  aPosMLPy[aIdxMLP] = S.Y();
                  aPosMLPz[aIdxMLP] = S.Z();
               }
               
               arSplineMLPx[idxSplineMLP] = S.X(); 
               arSplineMLPy[idxSplineMLP] = S.Y();
               arSplineMLPz[idxSplineMLP] = S.Z();

               S = SplineMLP(t, X0tps, X2, P0tps, P2, Lambda0, Lambda2); // (0,0)
               arSplineMLPNoTrkx[idxSplineMLP] = S.X();
               arSplineMLPNoTrky[idxSplineMLP] = S.Y();

               if (lastEID < maxAcc) {
                  aPosMLPNoTrkx[aIdxMLP] = S.X();
                  aPosMLPNoTrky[aIdxMLP] = S.Y();
               }
               
               S = SplineMLP(t, X0est, X2, P0tps, P2, Lambda0, Lambda2); // mine
               arSplineMLPestx[idxSplineMLP] = S.X();
               arSplineMLPesty[idxSplineMLP] = S.Y();
               
               if (lastEID < maxAcc) {
                  aPosMLPesty[aIdxMLP] = S.Y();
                  aPosMLPestx[aIdxMLP] = S.X();
               }

               S = SplineMLP(t, X0krah, X2, P0krah, P2, Lambda0, Lambda2); // niels'
               arSplineMLPKrahx[idxSplineMLP] = S.X();
               arSplineMLPKrahy[idxSplineMLP++] = S.Y();

               if (lastEID < maxAcc) {
                  aPosMLPKrahx[aIdxMLP] = S.X();
                  aPosMLPKrahy[aIdxMLP++] = S.Y();
               }
            }

            // Compare MC and MLP here
            // 1) Make splines to compare at same Z
            splineMCx  = new TSpline3("splineMCx", arSplineMCz, arSplineMCx, idxSplineMC);
            splineMCy  = new TSpline3("splineMCy", arSplineMCz, arSplineMCy, idxSplineMC);
            splineMLPx = new TSpline3("splineMLPx", arSplineMLPz, arSplineMLPx, idxSplineMLP);
            splineMLPy = new TSpline3("splineMLPy", arSplineMLPz, arSplineMLPy, idxSplineMLP);
            splineMLPestx = new TSpline3("splineMLPestx", arSplineMLPz, arSplineMLPestx, idxSplineMLP);
            splineMLPesty = new TSpline3("splineMLPesty", arSplineMLPz, arSplineMLPesty, idxSplineMLP);
            splineMLPNoTrkx = new TSpline3("splineMLPNoTrkx", arSplineMLPz, arSplineMLPNoTrkx, idxSplineMLP);
            splineMLPNoTrky = new TSpline3("splineMLPNoTrky", arSplineMLPz, arSplineMLPNoTrky, idxSplineMLP);
            splineMLPKrahx = new TSpline3("splineMLPKrahx", arSplineMLPz, arSplineMLPKrahx, idxSplineMLP);
            splineMLPKrahy = new TSpline3("splineMLPKrahy", arSplineMLPz, arSplineMLPKrahy, idxSplineMLP);

            // 2) Sweep Z values and add the absolute difference to an array
            float rvalue = splineMCx->Eval(0);
            if (!isnan(rvalue)) {
               Double_t diff_x, diff_y;
               idxDifferenceArray = 0; // Sweep each time
               nInDifferenceArray++; // To average at end
               for (Double_t zSweep = 0; zSweep <= phantomSize; zSweep += 2) { // Keep Sweep increment high to speed up!
                  differenceArrayZ[idxDifferenceArray] = zSweep;
                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPx->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPy->Eval(zSweep));
                  differenceArrayDiff[idxDifferenceArray] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));

                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPKrahx->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPKrahy->Eval(zSweep));
                  differenceArrayDiffKrah[idxDifferenceArray] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));
                  
                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPNoTrkx->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPNoTrky->Eval(zSweep));
                  differenceArrayDiffNoTrk[idxDifferenceArray] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));

                  diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPestx->Eval(zSweep));
                  diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPesty->Eval(zSweep));
                  differenceArrayDiffest[idxDifferenceArray++] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));
               }
            }


            delete splineMCx;
            delete splineMCy;
            delete splineMLPx;
            delete splineMLPy;
            delete splineMLPKrahx;
            delete splineMLPKrahy;
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

   f10xAvg /= f10xN;
   f10yAvg /= f10yN;

   printf("The 10%% MTF in X direction is %.3f lp/mm.\n", f10xAvg);
   printf("The 10%% MTF in Y direction is %.3f lp/mm.\n", f10yAvg);


   TCanvas *c0 = new TCanvas("c0", "Residual Energy", 1500, 600);
   hResEnergy->Draw();

   TCanvas *c1 = new TCanvas("c1", "MLP estimation", 1500, 1200);
   TPad *pad1 = new TPad("pad1", "pad1", 0.005, 0, 0.2855, .995);
   TPad *pad2 = new TPad("pad2", "pad2", .2865, 0, 0.5235, .995);
   TPad *pad3 = new TPad("pad3", "pad3", .5245, 0, 0.7615, .995);
   TPad *pad4 = new TPad("pad4", "pad4", 0.7625, 0, 0.995, .995);

   pad1->Divide(1, 2, 1e-5, 1e-5);
   pad2->Divide(1, 2, 1e-5, 1e-5);
   pad3->Divide(1, 2, 1e-5, 1e-5);
   pad4->Divide(1, 2, 1e-5, 1e-5);

   pad1->Draw(); pad2->Draw(); pad3->Draw(); pad4->Draw();

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

   // 2 and 6: Naive model
   pad2->cd(1);
   TGraph *gMCNoTrk = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPNoTrk = new TGraph(aIdxMLP, aPosMLPz, aPosMLPNoTrkx);
   gMLPNoTrk->SetMarkerStyle(21);
   gMLPNoTrk->SetMarkerColor(kBlue);
   gMCNoTrk->SetMarkerStyle(21);
   gMCNoTrk->SetMarkerColor(kRed);
   gMCNoTrk->SetMarkerSize(0.3);
   gMLPNoTrk->SetMarkerSize(0.3);
   gMCNoTrk->Draw("AP");
   gMLPNoTrk->Draw("P");
   gMCNoTrk->SetTitle("No Front Tracker - #hat{X}_{0} = (0,0,z_{0});Depth [mm];X position [mm]");
   
   pad2->cd(2);
   TGraph *gMCNoTrky = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPNoTrky = new TGraph(aIdxMLP, aPosMLPz, aPosMLPNoTrky);
   gMLPNoTrky->SetMarkerStyle(21);
   gMLPNoTrky->SetMarkerColor(kBlue);
   gMCNoTrky->SetMarkerStyle(21);
   gMCNoTrky->SetMarkerColor(kRed);
   gMCNoTrky->SetMarkerSize(0.3);
   gMLPNoTrky->SetMarkerSize(0.3);
   gMCNoTrky->Draw("AP");
   gMLPNoTrky->Draw("P");
   gMCNoTrky->SetTitle("No Front Tracker - #hat{X}_{0} = (0,0,z_{0});Depth [mm];X position [mm]");

   pad3->cd(1);
   TGraph *gMCKrah = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPKrah = new TGraph(aIdxMLP, aPosMLPz, aPosMLPKrahx);
   gMLPKrah->SetMarkerStyle(21);
   gMLPKrah->SetMarkerColor(kBlue);
   gMCKrah->SetMarkerStyle(21);
   gMCKrah->SetMarkerColor(kRed);
   gMCKrah->SetMarkerSize(0.3);
   gMLPKrah->SetMarkerSize(0.3);
   gMCKrah->Draw("AP");
   gMLPKrah->Draw("P");
   gMCKrah->SetTitle("No Front Tracker - Bayesian MLP;Depth [mm];X position [mm]");
   
   pad3->cd(2);
   TGraph *gMCKrahy = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPKrahy = new TGraph(aIdxMLP, aPosMLPz, aPosMLPKrahy);
   gMLPKrahy->SetMarkerStyle(21);
   gMLPKrahy->SetMarkerColor(kBlue);
   gMCKrahy->SetMarkerStyle(21);
   gMCKrahy->SetMarkerColor(kRed);
   gMCKrahy->SetMarkerSize(0.3);
   gMLPKrahy->SetMarkerSize(0.3);
   gMCKrahy->Draw("AP");
   gMLPKrahy->Draw("P");
   gMCKrahy->SetTitle("No Front Tracker - Bayesian MLP;Depth [mm];Y position [mm]");

   pad4->cd(1);
   TGraph *gMCest = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPest = new TGraph(aIdxMLP, aPosMLPz, aPosMLPestx);
   gMLPest->SetMarkerStyle(21);
   gMLPest->SetMarkerColor(kBlue);
   gMCest->SetMarkerStyle(21);
   gMCest->SetMarkerColor(kRed);
   gMCest->SetMarkerSize(0.3);
   gMLPest->SetMarkerSize(0.3);
   gMCest->Draw("AP");
   gMLPest->Draw("P");
   gMCest->SetTitle("No Front Tracker - Projection Model;Depth [mm];X position [mm]");
   
   pad4->cd(2);
   TGraph *gMCesty = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPesty = new TGraph(aIdxMLP, aPosMLPz, aPosMLPesty);
   gMLPesty->SetMarkerStyle(21);
   gMLPesty->SetMarkerColor(kBlue);
   gMCesty->SetMarkerStyle(21);
   gMCesty->SetMarkerColor(kRed);
   gMCesty->SetMarkerSize(0.3);
   gMLPesty->SetMarkerSize(0.3);
   gMCesty->Draw("AP");
   gMLPesty->Draw("P");
   gMCesty->SetTitle("No Front Tracker - Projection Model;Depth [mm];Y position [mm]");

   TCanvas *c = new TCanvas("c", "Beam spot estimation", 1570, 580);
   c->Divide(3, 1, 0.001, 0.001);

   c->cd(1);
   hErrorNaive->Draw("COLZ");
   c->cd(2);
   hErrorSigmaScale->Draw("COLZ");
   c->cd(3);
   hErrorKrah->Draw("COLZ");

   TCanvas *c1D = new TCanvas("c1D", "beam spot estimation", 1000, 600);
   hErrorKrah1D->Draw();

   for (Int_t i=0; i<idxDifferenceArray; i++) {
      differenceArrayDiff[i] /= nInDifferenceArray;
      differenceArrayDiffNoTrk[i]  /= nInDifferenceArray;
      differenceArrayDiffKrah[i]  /= nInDifferenceArray;
      differenceArrayDiffest[i] /= nInDifferenceArray;
//      differenceArrayZ[i] += phantomSize/2;
   }

   TCanvas *c2 = new TCanvas("c2", "MC vs MLP difference", 600, 600);
   TGraph *gDifferencePerfect = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiff);
   TGraph *gDifferenceNoTrk = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffNoTrk);
   TGraph *gDifferenceKrah = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffKrah);
   TGraph *gDifferenceest = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffest);

   gDifferencePerfect->SetTitle(";Depth in phantom [mm];Error | X_{1}^{opt} - X_{1}^{MC} | [mm]");
   gDifferencePerfect->SetLineColor(kRed);
   gDifferencePerfect->SetLineWidth(2);
   gDifferenceNoTrk->SetLineColor(kBlack);
   gDifferenceNoTrk->SetLineWidth(2);
   gDifferenceKrah->SetLineColor(kOrange+2);
   gDifferenceKrah->SetLineWidth(2);
   gDifferenceest->SetLineColor(kBlue);
   gDifferenceest->SetLineWidth(2);
   gDifferencePerfect->Draw("AL");
   gDifferenceNoTrk->Draw("L");
   gDifferenceKrah->Draw("L");
   gDifferenceest->Draw("L");
   gDifferencePerfect->GetYaxis()->SetRangeUser(0,4.5);

   TLegend *l = new TLegend(.6, .64, .85, .84);
   l->AddEntry(gDifferencePerfect, "Measure X_{0}", "L");
   l->AddEntry(gDifferenceNoTrk, "X_{0}^{opt} = X_{0}^{TPS}", "L");
   l->AddEntry(gDifferenceKrah, "Krah + MLP", "L");
   l->AddEntry(gDifferenceest, "LPM + CSP", "L");
   l->Draw();

   Float_t sigmaNoTrk = (hErrorNaive->GetStdDev(1) + hErrorNaive->GetStdDev(2)) / 2;
   Float_t sigmaEst = (hErrorSigmaScale->GetStdDev(1) + hErrorSigmaScale->GetStdDev(2)) / 2;

   if (spotsize >= 0) {
      ofstream file(Form("Output/MLPerror_energy%.0fMeV_%s_spotsize_krah.csv", initialEnergy, sMaterial), ofstream::out | ofstream::app);
//      file << phantomSize << " " << spotsize << " " <<  gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << sigmaNoTrk << " " << sigmaEst << endl;
      file << phantomSize << " " <<  spotsize << " " << gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << gDifferenceKrah->Eval(0) << " " << gDifferenceKrah->Eval(phantomSize/2) << " " << hResEnergy->GetMean() << endl;
      file.close();
   }
   else if (rotation >= 0) {
      ofstream file(Form("Output/MLPerror_energy%.0fMeV_%s_rotation.csv", initialEnergy, sMaterial), ofstream::out | ofstream::app);
      file << phantomSize << " " << rotation << " " <<  gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << sigmaNoTrk << " " << sigmaEst << endl;
      file.close();
   }
   else {
      ofstream file(Form("Output/MLPerror_energy%.0fMeV_%s_krah.csv", initialEnergy, sMaterial), ofstream::out | ofstream::app);
      file << phantomSize << " " <<  gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << gDifferenceKrah->Eval(0) << " " << gDifferenceKrah->Eval(phantomSize/2) << " " << hResEnergy->GetMean() << " " << f10xAvg << endl;
      file.close();
   }

   printf("The Krah error is %.2f mm. The LPM error is %.2f mm.\n", gDifferenceKrah->Eval(0), gDifferenceest->Eval(0));

}


float Sigmat1(float position)
{
	float p = position;
	float sigt1 = (azero*p)+(aone*p*p/2)+(atwo*p*p*p/3)+(athree*p*p*p*p/4)+(afour*p*p*p*p*p/5)+(afive*p*p*p*p*p*p/6);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigt1/X_0);
}

float Sigmaz1(float position)
{
	float p = position;
	float sigz1 = (azero*p*p*p/3)+(aone*p*p*p*p/12)+(atwo*p*p*p*p*p/30)+(athree*p*p*p*p*p*p/60)+(afour*p*p*p*p*p*p*p/105)+(afive*p*p*p*p*p*p*p*p/168);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigz1/X_0);
}

float Sigmatz1(float position)
{
	float p = position;
	float sigtz1 = (azero*p*p/2)+(aone*p*p*p/6)+(atwo*p*p*p*p/12)+(athree*p*p*p*p*p/20)+(afour*p*p*p*p*p*p/30)+(afive*p*p*p*p*p*p*p/42);
	
	return (13.6*13.6*pow((1+0.038*log(position/X_0)),2)*sigtz1/X_0);
}

float Sigmat2(float sep, float position)
{
	float p = position;
	float s = sep;
	float sigt2 = ((azero*s)+(aone*s*s/2)+(atwo*s*s*s/3)+(athree*s*s*s*s/4)+(afour*s*s*s*s*s/5)+(afive*s*s*s*s*s*s/6))-((azero*p)+(aone*p*p/2)+(atwo*p*p*p/3)+(athree*p*p*p*p/4)+(afour*p*p*p*p*p/5)+(afive*p*p*p*p*p*p/6));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigt2/X_0);
}

float Sigmaz2(float sep, float position)
{
	float p = position;
	float s = sep;
	float sigz2 = (azero*s*s*s/3)+(aone*s*s*s*s/12)+(atwo*s*s*s*s*s/30)+(athree*s*s*s*s*s*s/60)+(afour*s*s*s*s*s*s*s/105)+(afive*s*s*s*s*s*s*s*s/168)-((azero*s*s*p)+(((aone*s*s/2)-(azero*s))*p*p)+(((atwo*s*s/3)-(2*aone*s/3)+(azero/3))*p*p*p)+(((athree*s*s/4)-(atwo*s/2)+(aone/4))*p*p*p*p)+(((afour*s*s/5)-(2*athree*s/5)+(atwo/5))*p*p*p*p*p)+(((afive*s*s/6)-(afour*s/3)+(athree/6))*p*p*p*p*p*p)+(((afour/7)-(2*afive*s/7))*p*p*p*p*p*p*p)+(afive*p*p*p*p*p*p*p*p/8));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigz2/X_0);
}

float Sigmatz2(float sep, float position)
{
	float p = position;
	float s = sep;
	float sigtz2 = ((azero*s*s/2)+(aone*s*s*s/6)+(atwo*s*s*s*s/12)+(athree*s*s*s*s*s/20)+(afour*s*s*s*s*s*s/30)+(afive*s*s*s*s*s*s*s/42))-((azero*s*p)+(((aone*s)-azero)*p*p/2)+(((atwo*s)-aone)*p*p*p/3)+(((athree*s)-atwo)*p*p*p*p/4)+(((afour*s)-athree)*p*p*p*p*p/5)+(((afive*s)-afour)*p*p*p*p*p*p/6)-(afive*p*p*p*p*p*p*p/7));
	
	return (13.6*13.6*pow((1+0.038*log((sep-position)/X_0)),2)*sigtz2/X_0);
}
