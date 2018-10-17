#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <fstream>
#include <TBranch.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TStyle.h>
#include <Math/Vector3D.h>

#define azero 7.457e-6
#define aone 4.548e-7
#define atwo (-5.777e-8)
#define athree 1.301e-8
#define afour (-9.228e-10)
#define afive 2.687e-11
#define X_0 36.1

using namespace std;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > XYZVector;

float Sigmat1(float);
float Sigmaz1(float);
float Sigmatz1(float);
float Sigmat2(float, float);
float Sigmaz2(float, float);
float Sigmatz2(float, float);

XYZVector SplineMLP(Double_t t, XYZVector X0, XYZVector X1, XYZVector P0, XYZVector P1, Double_t Lambda0, Double_t Lambda1) {
   XYZVector P0Lambda, P1Lambda, X1mX0, S;
   Float_t tt = pow(t, 2), ttt = pow(t, 3);

   X1mX0 = X1 - X0;
   P0Lambda = P0 * Lambda0 * sqrt(X1mX0.Mag2());
   P1Lambda = P1 * Lambda1 * sqrt(X1mX0.Mag2());

   S = (2*ttt - 3*tt + 1) * X0 + (ttt - 2*tt + t) * P0Lambda + (-2*ttt + 3*tt) * X1 + (ttt - tt) * P1Lambda;

   return S;
}


void findMLP(Float_t phantomSize = 200, Float_t rotation = -1, Float_t spotsize = -1, Float_t initialEnergy = 200) {
   TFile *f = nullptr;

   if (rotation >= 0) {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%03.0fmm_rotation%02.0fmrad.root", initialEnergy, phantomSize, rotation));
   }

   else if (spotsize >= 0)  {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%03.0fmm_spotsize%04.1fmm.root", initialEnergy, phantomSize, spotsize));
   }
   
   else {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%03.0fmm.root", initialEnergy, phantomSize));
   }

   TTree *tree = (TTree*) f->Get("Hits");

   if (!tree) exit(0);

   Int_t      printed = 0;
   const Int_t eventsToUse = 10000;
   Float_t     x, y, z, edep, sum_edep = 0, residualEnergy = 0;
   Int_t       eventID, parentID, lastEID = -1;
   XYZVector   Xp0, Xp1, Xp2, Xp3, X0, X1, X0est, X0err, X0NoTrk, P0, P0NoTrk, P1, P0hat, P1hat, S, P1Rotated; // Xp are the plane coordinates, X are the tracker coordinates (X1 = (Xp1 + Xp2) / 2)
   Float_t     wepl, wet, Lambda0, Lambda1;
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
   Double_t    aPosMLPestx[10000], aPosMLPesty[10000];
   Int_t       volumeID[10];
   TBranch   * b_volumeID;
   Int_t       aIdxMC = 0;
   Int_t       aIdxMLP = 0;
   Bool_t      stop = false, stopacc = false;

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
   TH1F *hP1x = new TH1F("hP1x", "Outgoing angle histgram;#theta_{x} [mrad];Frequency", 200, 0, 250);
   TH1F *hP1y = new TH1F("hP1y", "Outgoing angle histgram;#theta_{y} [mrad];Frequency", 200, 0, 250);
   TH1F *hP2x = new TH1F("hP2x", "Outgoing angle histgram;#theta_{x} [mrad];Frequency", 200, 0, 250);
   TH1F *hP2y = new TH1F("hP2y", "Outgoing angle histgram;#theta_{y} [mrad];Frequency", 200, 0, 250);
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
   
   if (rotation < 0) {
      P0NoTrk.SetCoordinates(0, 0, 1); // Normalized
      X0NoTrk.SetCoordinates(0, 0,-phantomSize/2);
   }

   else {
      P0NoTrk.SetCoordinates(0, -sin(rotation / 1000), cos(rotation / 1000));
      X0NoTrk.SetCoordinates(0, -sin(rotation / 1000) * (sourceToX0dist + 15), -phantomSize/2);
   }

   Int_t maxAcc = 10;

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);

      if (eventID > eventsToUse) stop = true;
      if (eventID > 5) stopacc = true;

      if (lastEID < 0) {
         lastEID = eventID;
         sum_edep = edep;
         residualEnergy = 0;
      }

      if (lastEID != eventID) {
         X0 = (Xp0 + Xp1) / 2;
         X1 = (Xp2 + Xp3) / 2;
         P0 = Xp1 - Xp0;
         P1 = Xp3 - Xp2;
         P0hat = P0.Unit();
         P1hat = P1.Unit();

         P1Rotated = P1hat - P0NoTrk;
         P1Rotated.SetZ(P1hat.Z());
         P1Rotated = P1Rotated.Unit();
            
         wepl = splineWater->Eval(initialEnergy);
         wet = wepl - splineWater->Eval(residualEnergy);
         Float_t w = wet / wepl;
         Float_t w2 = pow(wet / wepl, 2);
            
         Float_t AX = 1 - 0.185 * w + 0.372 * pow(w,2) - 0.916 * pow(w,3);
         Float_t AP = -3.93 - 82.23 * w - 185.6 * pow(w,2) + 273.9 * pow(w,3);
        
         if (spotsize >= 0) {
            AX = 0.2333 * log(spotsize) + 0.4886;
            AP = -17.67 * log(spotsize) - 31.793;
         }

         Float_t dxy = sqrt(pow(P1Rotated.X(), 2) + pow(P1Rotated.Y(), 2));
         Float_t angle = fabs(atan2(dxy,1)) * 1000;
         hP2x->Fill(atan2(P0hat.X(), 1)*1000);
         hP2y->Fill(atan2(P0hat.Y(), 1)*1000);

         sigmaFilter = 139.15 * w + 13.991;
         printf("sigmaFilter = %.2f mrad.\n", sigmaFilter);

         if (residualEnergy > energyFilter && angle < sigmaFilter) {
            hResEnergy->Fill(residualEnergy);
         
            double trackerDist = 1;
            double trackerDistPhantom = 1;

            /*
            double px0 = P0hat.X();
            double py0 = P0hat.Y(); 
            double pz0 = P0hat.Z(); 
            */

            double px0 = P0NoTrk.X();
            double py0 = P0NoTrk.Y(); 
            double pz0 = P0NoTrk.Z(); 

            double x0 = X0NoTrk.X();
            double y0 = X0NoTrk.Y();
            double z0 = X0NoTrk.Z();
            
            double angleX2rad = atan2(px0, pz0);
            double angleY2rad = atan2(py0, pz0);
            double angleZ2rad = atan2(pz0, pz0);
            
            double angleXout = atan(P1hat.X());
            double angleYout = atan(P1hat.Y());
            double angleZout = atan(P1hat.Z());

            TVector3 m0(x0/10, y0/10, z0/10 + trackerDistPhantom); // in cm
            TVector3 p0(angleX2rad, angleY2rad, angleZ2rad);

            TVector3 m1(Xp2.X()/10 - trackerDistPhantom * tan(angleXout),
                        Xp2.Y()/10 - trackerDistPhantom * tan(angleYout),
                        Xp2.Z()/10 - trackerDistPhantom);
            TVector3 p1(angleXout, angleYout, angleZout);

            double step_length = (m1.z() - m0.z()) / 512;
            double posz = m0.z() + step_length;

            float sz1, sz2, st1, st2, stz1, stz2;
            float determinant_1, determinant_2, determinant_C12;
            float d_source = 10; // assume to first tracker layer
            float s_pos = 0.3; // cm
            float s_angle = 0.0025;

            if (spotsize >= 0) s_pos = spotsize / 10;
      
            // init MLP
            float X_mlp, Y_mlp;
            TVector3 p;

            int a; // iterator for matrix operations
            float scatter_1[4] = {0};
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

            sz1 = Sigmaz1(posz - m0.z());
            st1 = Sigmat1(posz - m0.z());
            stz1 = Sigmatz1(posz - m0.z());

            sz2 = Sigmaz2(m1.z() - m0.z(), posz - m0.z());
            stz2 = Sigmatz2(m1.z() - m0.z(), posz - m0.z());
            st2 = Sigmat2(m1.z() - m0.z(), posz - m0.z());

            R_0[0] = 1;
            R_0[1] = posz - m0.z();
            R_0[2] = 0;
            R_0[3] = 1;

            R_0_transpose[0] = 1;
            R_0_transpose[1] = 0;
            R_0_transpose[2] = posz - m0.z();
            R_0_transpose[3] = 1;

            R_1_inverse[0] = 1;
            R_1_inverse[1] = -(m1.z() - posz); // ok
            R_1_inverse[2] = 0;
            R_1_inverse[3] = 1;

            R_1_inverse_transpose[0] = 1;
            R_1_inverse_transpose[1] = 0;
            R_1_inverse_transpose[2] = -(m1.z() - posz);
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

            for (a=0; a<4; a++) {
               C1[a] = C1_2[a] + scatter_1[a];
            }

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

            y_0[0] = m0.x();
            y_0[1] = angleX2rad;
            y_2[0] = m1.x();
            y_2[1] = angleXout;
   
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

            // now do the y value
            y_0[0] = m0.y();
            y_0[1] = angleY2rad;

            y_2[0] = m1.y();
            y_2[1] = angleYout;
            
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
            hP1x->Fill(theta_X_mlp * 1000);
            hP1y->Fill(theta_Y_mlp * 1000);

            XYZVector X0krah(X_mlp*10, Y_mlp*10, -phantomSize/2);
            XYZVector P0krah(tan(theta_X_mlp), tan(theta_Y_mlp), 1);
            P0krah = P0krah.Unit();

            XYZVector projectedPathX0;
            projectedPathX0.SetCoordinates(15 * P0NoTrk.X(), 15 * P0NoTrk.Y(), 15);

            XYZVector projectedPath;
            projectedPath.SetCoordinates(X0NoTrk.X() + phantomSize * P0NoTrk.X(), X0NoTrk.Y() + phantomSize * P0NoTrk.Y(), 0);
         
            X0 += projectedPathX0;

            X1.SetCoordinates(X1.X() - 15 * P1hat.X(), X1.Y() - 15 * P1hat.Y(), X1.Z() - 15);

            X0est = X1 * AX + P1Rotated * AP;
            X0est.SetZ(-phantomSize/2);

            XYZVector X0errNaive, X0errKrah;
            X0err = X0est - X0;
            X0errNaive = X0NoTrk - X0;
            X0errKrah = X0krah - X0;

            hErrorSigmaScale->Fill(X0err.X(), X0err.Y());
            hErrorNaive->Fill(X0errNaive.X(), X0errNaive.Y());
            hErrorKrah->Fill(X0errKrah.X(), X0errKrah.Y());
            hErrorKrah1D->Fill(X0errKrah.X() / 10 ) ;

            // Find lambda values using polynomial in Fig. 4 in Fekete et al. 2015
            Lambda0 = 1.01 + 0.43 * w2;
            Lambda1 = 0.99 - 0.46 * w2;

            for (Float_t t=0; t<1; t += 0.01) {
               S = SplineMLP(t, X0, X1, P0hat, P1hat, Lambda0, Lambda1); // perfect info
               if (lastEID < maxAcc) {
                  aPosMLPx[aIdxMLP] = S.X();
                  aPosMLPy[aIdxMLP] = S.Y();
                  aPosMLPz[aIdxMLP] = S.Z();
               }
               
               arSplineMLPx[idxSplineMLP] = S.X(); 
               arSplineMLPy[idxSplineMLP] = S.Y();
               arSplineMLPz[idxSplineMLP] = S.Z();

               S = SplineMLP(t, X0NoTrk, X1, P0NoTrk, P1hat, Lambda0, Lambda1); // (0,0)
               arSplineMLPNoTrkx[idxSplineMLP] = S.X();
               arSplineMLPNoTrky[idxSplineMLP] = S.Y();
               
               S = SplineMLP(t, X0est, X1, P0NoTrk, P1hat, Lambda0, Lambda1); // mine
               arSplineMLPestx[idxSplineMLP] = S.X();
               arSplineMLPesty[idxSplineMLP] = S.Y();
               
               if (lastEID < maxAcc) {
                  aPosMLPesty[aIdxMLP] = S.Y();
                  aPosMLPestx[aIdxMLP] = S.X();
               }

               S = SplineMLP(t, X0krah, X1, P0krah, P1hat, Lambda0, Lambda1); // niels'
               arSplineMLPKrahx[idxSplineMLP] = S.X();
               arSplineMLPKrahy[idxSplineMLP++] = S.Y();

               if (lastEID < maxAcc) {
                  aPosMLPNoTrkx[aIdxMLP] = S.X();
                  aPosMLPNoTrky[aIdxMLP++] = S.Y();
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
               for (Double_t zSweep = -phantomSize/2; zSweep <= phantomSize/2; zSweep += 2) { // Keep Sweep increment high to speed up!
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
         if       (volumeID[2] == 0) Xp0.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 1) Xp1.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 2) Xp2.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 3) Xp3.SetCoordinates(x,y,z);
         else if  (volumeID[2] == 5) residualEnergy += edep;


         if  (volumeID[2] < 5) {
            arSplineMCx[idxSplineMC] = x;
            arSplineMCy[idxSplineMC] = y;
            arSplineMCz[idxSplineMC++] = z;
            
            if (eventID < maxAcc) {
               aPosMCx[aIdxMC] = x;
               aPosMCy[aIdxMC] = y;
               aPosMCz[aIdxMC++] = z;
            }
         }
      }
   }

   TCanvas *c0 = new TCanvas("c0", "Residual Energy", 1500, 600);
   c0->Divide(3,1,0.001,0.001);
   c0->cd(1);
   hResEnergy->Draw();
   c0->cd(2);
   hP1x->Draw();
   hP1x->SetFillColor(kRed);
   hP2x->Draw("same");
   c0->cd(3);
   hP1y->Draw();
   hP1y->SetFillColor(kRed);
   hP2y->Draw("same");

   TCanvas *c1 = new TCanvas("c1", "MLP estimation", 1500, 1200);
   c1->Divide(3, 2, 0.001, 0.001);
   c1->cd(1);
   TGraph *gMC = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMCy = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLP = new TGraph(aIdxMLP, aPosMLPz, aPosMLPx);
   TGraph *gMLPy = new TGraph(aIdxMLP, aPosMLPz, aPosMLPy);
   gMC->SetMarkerStyle(7);
   gMLP->SetMarkerStyle(7);
   gMC->SetMarkerColor(kRed);
   gMLP->SetMarkerColor(kBlue);
   gMC->Draw("AP");
   gMLP->Draw("P");
   gMC->SetTitle("Perfect knowledge;Depth [mm];X position [mm]");

   c1->cd(4);
   gMCy->SetMarkerStyle(7);
   gMLPy->SetMarkerStyle(7);
   gMCy->SetMarkerColor(kRed);
   gMLPy->SetMarkerColor(kBlue);
   gMCy->Draw("AP");
   gMLPy->Draw("P");
   gMCy->SetTitle("Perfect knowledge;Depth [mm];Y position [mm]");

   c1->cd(2);
   TGraph *gMCNoTrk = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPNoTrk = new TGraph(aIdxMLP, aPosMLPz, aPosMLPNoTrkx);
   gMLPNoTrk->SetMarkerStyle(7);
   gMLPNoTrk->SetMarkerColor(kBlue);
   gMCNoTrk->SetMarkerStyle(7);
   gMCNoTrk->SetMarkerColor(kRed);
   gMCNoTrk->Draw("AP");
   gMLPNoTrk->Draw("P");
   gMCNoTrk->SetTitle("No Front Tracker - use X_{0}^{est} = (0,0,z_{0});Depth [mm];X position [mm]");
   
   c1->cd(5);
   TGraph *gMCNoTrky = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPNoTrky = new TGraph(aIdxMLP, aPosMLPz, aPosMLPNoTrky);
   gMLPNoTrky->SetMarkerStyle(7);
   gMLPNoTrky->SetMarkerColor(kBlue);
   gMCNoTrky->SetMarkerStyle(7);
   gMCNoTrky->SetMarkerColor(kRed);
   gMCNoTrky->Draw("AP");
   gMLPNoTrky->Draw("P");
   gMCNoTrky->SetTitle("No Front Tracker - use Bayesian MLP;Depth [mm];Y position [mm]");

   c1->cd(3);
   TGraph *gMCest = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPest = new TGraph(aIdxMLP, aPosMLPz, aPosMLPestx);
   gMLPest->SetMarkerStyle(7);
   gMLPest->SetMarkerColor(kBlue);
   gMCest->SetMarkerStyle(7);
   gMCest->SetMarkerColor(kRed);
   gMCest->Draw("AP");
   gMLPest->Draw("P");
   gMCest->SetTitle("No Front Tracker - Projection Model;Depth [mm];X position [mm]");
   
   c1->cd(6);
   TGraph *gMCesty = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPesty = new TGraph(aIdxMLP, aPosMLPz, aPosMLPesty);
   gMLPesty->SetMarkerStyle(7);
   gMLPesty->SetMarkerColor(kBlue);
   gMCesty->SetMarkerStyle(7);
   gMCesty->SetMarkerColor(kRed);
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

   printf("Normalizing the %d indexes of differenceArray from %d.\n", idxDifferenceArray, nInDifferenceArray);
   for (Int_t i=0; i<idxDifferenceArray; i++) {
      differenceArrayDiff[i] /= nInDifferenceArray;
      differenceArrayDiffNoTrk[i]  /= nInDifferenceArray;
      differenceArrayDiffKrah[i]  /= nInDifferenceArray;
      differenceArrayDiffest[i] /= nInDifferenceArray;
      differenceArrayZ[i] += phantomSize/2;
   }

   TCanvas *c2 = new TCanvas("c2", "MC vs MLP difference", 1500, 600);
   c2->Divide(4,1,0.001,0.001);
   c2->cd(1);
   TGraph *gDifferencePerfect = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiff);
   gPad->SetGridy();
   gDifferencePerfect->SetTitle("Perfect knowledge;Depth in phantom [mm];Error MLP - MC [mm]");
   gDifferencePerfect->SetLineColor(kRed);
   gDifferencePerfect->SetLineWidth(2);
   gDifferencePerfect->Draw("AL");
   gDifferencePerfect->GetYaxis()->SetRangeUser(0, 6);
   
   c2->cd(2);
   TGraph *gDifferenceNoTrk = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffNoTrk);
   gPad->SetGridy();
   gDifferenceNoTrk->SetTitle("No Front Tracker - use X_{0}^{est} = (0,0,z_{0});Depth in phantom [mm]; Error MLP - MC [mm]");
   gDifferenceNoTrk->SetLineColor(kRed);
   gDifferenceNoTrk->SetLineWidth(2);
   gDifferenceNoTrk->Draw("AL");
   gDifferenceNoTrk->GetYaxis()->SetRangeUser(0, 6);

   c2->cd(3);
   TGraph *gDifferenceKrah = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffKrah);
   gPad->SetGridy();
   gDifferenceKrah->SetTitle("No Front Tracker - use Bayesian MLP;Depth in phantom [mm]; Error MLP - MC [mm]");
   gDifferenceKrah->SetLineColor(kRed);
   gDifferenceKrah->SetLineWidth(2);
   gDifferenceKrah->Draw("AL");
   gDifferenceKrah->GetYaxis()->SetRangeUser(0, 6);

   c2->cd(4);
   TGraph *gDifferenceest = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffest);
   gPad->SetGridy();
   gDifferenceest->SetTitle("No Front Tracker - use Projection Model;Depth in phantom [mm]; Error MLP - MC [mm]");
   gDifferenceest->SetLineColor(kRed);
   gDifferenceest->SetLineWidth(2);
   gDifferenceest->Draw("AL");
   gDifferenceest->GetYaxis()->SetRangeUser(0, 6);

   Float_t sigmaNoTrk = (hErrorNaive->GetStdDev(1) + hErrorNaive->GetStdDev(2)) / 2;
   Float_t sigmaEst = (hErrorSigmaScale->GetStdDev(1) + hErrorSigmaScale->GetStdDev(2)) / 2;

   if (spotsize >= 0) {
      ofstream file(Form("Output/MLPerror_energy%.0fMeV_Water_spotsize_krah.csv", initialEnergy), ofstream::out | ofstream::app);
//      file << phantomSize << " " << spotsize << " " <<  gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << sigmaNoTrk << " " << sigmaEst << endl;
      file << phantomSize << " " <<  spotsize << " " << gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << gDifferenceKrah->Eval(0) << " " << gDifferenceKrah->Eval(phantomSize/2) << " " << hResEnergy->GetMean() << endl;
      file.close();
   }
   else if (rotation >= 0) {
      ofstream file(Form("Output/MLPerror_energy%.0fMeV_Water_rotation.csv", initialEnergy), ofstream::out | ofstream::app);
      file << phantomSize << " " << rotation << " " <<  gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << sigmaNoTrk << " " << sigmaEst << endl;
      file.close();
   }
   else {
      ofstream file(Form("Output/MLPerror_energy%.0fMeV_Water_krah.csv", initialEnergy), ofstream::out | ofstream::app);
      file << phantomSize << " " <<  gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << gDifferenceKrah->Eval(0) << " " << gDifferenceKrah->Eval(phantomSize/2) << " " << hResEnergy->GetMean() << endl;
      file.close();
   }


   /*
   delete gDifferencePerfect;
   delete gDifferenceNoTrk;
   delete gDifferenceest;

   delete c1;
   delete c0;
   delete c2;

   delete gMC;
   delete gMLP;
   delete gMCNoTrk;
   delete gMLPNoTrk;
   delete gMCest;
   delete gMLPest;

   f->Close();

   */
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
