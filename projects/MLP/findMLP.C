#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <fstream>
#include <TBranch.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TStyle.h>
#include <Math/Vector3D.h>

using namespace std;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > XYZVector;

XYZVector SplineMLP(Double_t t, XYZVector X0, XYZVector X1, XYZVector P0, XYZVector P1, Double_t Lambda0, Double_t Lambda1) {
   XYZVector P0Lambda, P1Lambda, X1mX0, S;
   Float_t tt = pow(t, 2), ttt = pow(t, 3);

   X1mX0 = X1 - X0;
   P0Lambda = P0 * Lambda0 * sqrt(X1mX0.Mag2());
   P1Lambda = P1 * Lambda1 * sqrt(X1mX0.Mag2());

   S = (2*ttt - 3*tt + 1) * X0 + (ttt - 2*tt + t) * P0Lambda + (-2*ttt + 3*tt) * X1 + (ttt - tt) * P1Lambda;

   return S;
}

void findMLP(Float_t divergence = -1, Float_t phantomSize = 160, Float_t initialEnergy = 230, Float_t AX = 0.955, Float_t AP = -9.975) {
   TFile *f = nullptr;

   if (divergence < 0) {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%.0fmm.root", initialEnergy, phantomSize));
   }

   else {
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%.0fmm_rotation%.0fmrad.root", initialEnergy, phantomSize, divergence));
   }

   TTree *tree = (TTree*) f->Get("Hits");

   Int_t      printed = 0;
   const Int_t eventsToUse = 10000;
   Float_t     x, y, z, edep, sum_edep = 0, residualEnergy = 0;
   Int_t       eventID, parentID, lastEID = -1;
   XYZVector   Xp0, Xp1, Xp2, Xp3, X0, X1, X0est, X0err, X0NoTrk, P0, P0NoTrk, P1, P0hat, P1hat, S; // Xp are the plane coordinates, X are the tracker coordinates (X1 = (Xp1 + Xp2) / 2)
   Float_t     wepl, wet, Lambda0, Lambda1;
   ifstream    in;
   TSpline3  * splineMCx, * splineMCy, * splineMLPx, * splineMLPy;
   TSpline3  * splineMLPNoTrkx, * splineMLPNoTrky, * splineMLPestx, * splineMLPesty;
   Double_t    arSplineMCx[1000], arSplineMCy[1000], arSplineMCz[1000], arSplineMLPx[1000], arSplineMLPy[1000], arSplineMLPz[1000];
   Double_t    arSplineMLPNoTrkx[1000], arSplineMLPNoTrky[1000], arSplineMLPestx[1000], arSplineMLPesty[1000];
   Int_t       idxSplineMC = 0, idxSplineMLP = 0;
   Double_t    differenceArrayZ[1000];
   Double_t    differenceArrayDiff[1000] = {};
   Double_t    differenceArrayDiffNoTrk[1000] = {};
   Double_t    differenceArrayDiffest[1000] = {};
   Int_t       nInDifferenceArray = 0;
   Int_t       idxDifferenceArray = 0;
   Int_t       idxWater = 0;
   Double_t    energiesWater[500], rangesWater[500];
   Float_t     energy, range;
//   Float_t     sigmaFilter = 128.47;
   Float_t     sigmaFilter =5.0;
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

   TH1F *hResEnergy = new TH1F("hResEnergy", "Residual energy in calorimeter;Energy [MeV];Entries", 300, 0, 230);
   TH2F *hErrorNaive = new TH2F("hErrorNaive", "Beamspot uncertainty (assume point beam);X position [mm];Y position [mm]", 100, -25, 25, 100, -25, 25);
   TH2F *hErrorSigmaScale = new TH2F("hErrorSigmaScale", Form("Beamspot uncertainty (est);X position [mm];Y position [mm]"), 100, -25, 25, 100, -25, 25);

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("volumeID", volumeID, &b_volumeID);
   
   X0NoTrk.SetCoordinates(0,-divergence*1.015, -phantomSize/2);
   if (divergence < 1) {
      P0NoTrk.SetCoordinates(0, 0, 1); // Normalized
   }
   else {
      P0NoTrk.SetCoordinates(0, -sin(divergence / 1000), cos(divergence / 1000));
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
         // New particle, store last values
         if (residualEnergy > sigmaFilter) {
            
            // calculate Spline
            wepl = splineWater->Eval(initialEnergy);
            wet = wepl - splineWater->Eval(residualEnergy);

            Float_t w = pow(wet / wepl, 2);
            AX = 1.02 - 0.38 * w;
            AP = -13.74 * (wet/wepl);

            if (printed < 10) {
               printf("w = %.2f -> using AX = %.2f, and AP = %.2f.\n", w, AX, AP);
               printed++;
            }

            hResEnergy->Fill(residualEnergy);

            // Find vectors as defined in paper
            X0 = (Xp0 + Xp1) / 2;
            X1 = (Xp2 + Xp3) / 2;
            P0 = Xp1 - Xp0;
            P1 = Xp3 - Xp2;
            P0hat = P0.Unit();
            P1hat = P1.Unit();

            X0.SetZ(X0.Z() + 15);
            X1.SetCoordinates(X1.X() - 15 * P1hat.X(), X1.Y() - 15 * P1hat.Y(), X1.Z() - 15);

            X0est = X1 * AX + P1 * AP;
            X0est.SetZ(-phantomSize/2);
            if (divergence) {
               X0est.SetX(X0est.X() - phantomSize * P0NoTrk.X());
               X0est.SetY(X0est.Y() - phantomSize * P0NoTrk.Y());
            }
            
            X0err = X1 * AX + P1 * AP - X0;

            hErrorNaive->Fill(X0.X(), X0.Y());
            hErrorSigmaScale->Fill(X0err.X(), X0err.Y());

            // Find lambda values using polynomial in Fig. 4 in Fekete et al. 2015
            Lambda0 = 1.01 + 0.43 * w;
            Lambda1 = 0.99 - 0.46 * w;

            for (Float_t t=0; t<1; t += 0.01) {
               S = SplineMLP(t, X0, X1, P0hat, P1hat, Lambda0, Lambda1);
               if (lastEID < maxAcc) {
                  aPosMLPx[aIdxMLP] = S.X(); // Collect many tracks 10 times
                  aPosMLPy[aIdxMLP] = S.Y();
                  aPosMLPz[aIdxMLP] = S.Z();
               }

               arSplineMLPx[idxSplineMLP] = S.X(); // Collect a single track for all eventIDs
               arSplineMLPy[idxSplineMLP] = S.Y();
               arSplineMLPz[idxSplineMLP] = S.Z();

               S = SplineMLP(t, X0NoTrk, X1, P0NoTrk, P1hat, Lambda0, Lambda1);
               arSplineMLPNoTrkx[idxSplineMLP] = S.X();
               arSplineMLPNoTrky[idxSplineMLP] = S.Y();


               if (lastEID < maxAcc) {
                  aPosMLPNoTrky[aIdxMLP] = S.Y();
                  aPosMLPNoTrkx[aIdxMLP] = S.X();
               }

               S = SplineMLP(t, X0est, X1, P0NoTrk, P1hat, Lambda0, Lambda1);
               arSplineMLPestx[idxSplineMLP] = S.X();
               arSplineMLPesty[idxSplineMLP++] = S.Y();
               
               if (lastEID < maxAcc) {
                  aPosMLPesty[aIdxMLP] = S.Y();
                  aPosMLPestx[aIdxMLP++] = S.X();
               }
            }
         
            // Compare MC and MLP here
            // 1) Make splines to compare at same Z
            splineMCx  = new TSpline3("splineMCx", arSplineMCz, arSplineMCx, idxSplineMC);
            splineMCy  = new TSpline3("splineMCy", arSplineMCz, arSplineMCy, idxSplineMC);
            splineMLPx = new TSpline3("splineMLPx", arSplineMLPz, arSplineMLPx, idxSplineMLP);
            splineMLPy = new TSpline3("splineMLPy", arSplineMLPz, arSplineMLPy, idxSplineMLP);
            splineMLPNoTrkx = new TSpline3("splineMLPNoTrkx", arSplineMLPz, arSplineMLPNoTrkx, idxSplineMLP);
            splineMLPNoTrky = new TSpline3("splineMLPNoTrky", arSplineMLPz, arSplineMLPNoTrky, idxSplineMLP);
            splineMLPestx = new TSpline3("splineMLPestx", arSplineMLPz, arSplineMLPestx, idxSplineMLP);
            splineMLPesty = new TSpline3("splineMLPesty", arSplineMLPz, arSplineMLPesty, idxSplineMLP);

            // 2) Sweep Z values and add the absolute difference to an array
            Double_t diff_x, diff_y;
            idxDifferenceArray = 0; // Sweep each time
            nInDifferenceArray++; // To average at end
            for (Double_t zSweep = -phantomSize/2; zSweep <= phantomSize/2; zSweep += 2) { // Keep Sweep increment high to speed up!
               differenceArrayZ[idxDifferenceArray] = zSweep;
            
               diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPx->Eval(zSweep));
               diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPy->Eval(zSweep));
               differenceArrayDiff[idxDifferenceArray] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));

               diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPNoTrkx->Eval(zSweep));
               diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPNoTrky->Eval(zSweep));
               differenceArrayDiffNoTrk[idxDifferenceArray] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));

               diff_x = fabs(splineMCx->Eval(zSweep) - splineMLPestx->Eval(zSweep));
               diff_y = fabs(splineMCy->Eval(zSweep) - splineMLPesty->Eval(zSweep));
               differenceArrayDiffest[idxDifferenceArray++] += sqrt(pow(diff_x, 2) + pow(diff_y, 2));
            }

            delete splineMCx;
            delete splineMCy;
            delete splineMLPx;
            delete splineMLPy;
            delete splineMLPNoTrkx;
            delete splineMLPNoTrky;
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
   hResEnergy->Draw();

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
   gMCNoTrk->SetTitle("No Front Tracker - use (0,0);Depth [mm];X position [mm]");
   
   c1->cd(5);
   TGraph *gMCNoTrky = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPNoTrky = new TGraph(aIdxMLP, aPosMLPz, aPosMLPNoTrky);
   gMLPNoTrky->SetMarkerStyle(7);
   gMLPNoTrky->SetMarkerColor(kBlue);
   gMCNoTrky->SetMarkerStyle(7);
   gMCNoTrky->SetMarkerColor(kRed);
   gMCNoTrky->Draw("AP");
   gMLPNoTrky->Draw("P");
   gMCNoTrky->SetTitle("No Front Tracker - use (0,0);Depth [mm];Y position [mm]");

   c1->cd(3);
   TGraph *gMCest = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLPest = new TGraph(aIdxMLP, aPosMLPz, aPosMLPestx);
   gMLPest->SetMarkerStyle(7);
   gMLPest->SetMarkerColor(kBlue);
   gMCest->SetMarkerStyle(7);
   gMCest->SetMarkerColor(kRed);
   gMCest->Draw("AP");
   gMLPest->Draw("P");
   gMCest->SetTitle("No Front Tracker - estimate PB;Depth [mm];X position [mm]");
   
   c1->cd(6);
   TGraph *gMCesty = new TGraph(aIdxMC, aPosMCz, aPosMCy);
   TGraph *gMLPesty = new TGraph(aIdxMLP, aPosMLPz, aPosMLPesty);
   gMLPesty->SetMarkerStyle(7);
   gMLPesty->SetMarkerColor(kBlue);
   gMCesty->SetMarkerStyle(7);
   gMCesty->SetMarkerColor(kRed);
   gMCesty->Draw("AP");
   gMLPesty->Draw("P");
   gMCesty->SetTitle("No Front Tracker - estimate PB;Depth [mm];Y position [mm]");

   TCanvas *c = new TCanvas("c", "Beam spot estimation", 1000, 600);
   c->Divide(2, 1, 0.001, 0.001);

   c->cd(1);
   hErrorNaive->Draw("COLZ");
   c->cd(2);
   hErrorSigmaScale->Draw("COLZ");


   printf("Normalizing the %d indexes of differenceArray from %d.\n", idxDifferenceArray, nInDifferenceArray);
   for (Int_t i=0; i<idxDifferenceArray; i++) {
      differenceArrayDiff[i] /= nInDifferenceArray;
      differenceArrayDiffNoTrk[i]  /= nInDifferenceArray;
      differenceArrayDiffest[i] /= nInDifferenceArray;
      differenceArrayZ[i] += phantomSize/2 + 10;
   }

   TCanvas *c2 = new TCanvas("c2", "MC vs MLP difference", 1500, 600);
   c2->Divide(3,1,0.001,0.001);
   c2->cd(1);
   TGraph *gDifferencePerfect = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiff);
   gPad->SetGridy();
   gDifferencePerfect->SetTitle("Perfect knowledge;Depth in phantom [mm];Error MLP - MC [mm]");
   gDifferencePerfect->SetLineColor(kRed);
   gDifferencePerfect->SetLineWidth(2);
   gDifferencePerfect->Draw("AL");
   gDifferencePerfect->GetYaxis()->SetRangeUser(0, 4);

   c2->cd(2);
   TGraph *gDifferenceNoTrk = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffNoTrk);
   gPad->SetGridy();
   gDifferenceNoTrk->SetTitle("No Front Tracker - use (0,0);Depth in phantom [mm]; Error MLP - MC [mm]");
   gDifferenceNoTrk->SetLineColor(kRed);
   gDifferenceNoTrk->SetLineWidth(2);
   gDifferenceNoTrk->Draw("AL");
   gDifferenceNoTrk->GetYaxis()->SetRangeUser(0, 4);

   c2->cd(3);
   TGraph *gDifferenceest = new TGraph(idxDifferenceArray, differenceArrayZ, differenceArrayDiffest);
   gPad->SetGridy();
   gDifferenceest->SetTitle("No Front Tracker - estimate PB;Depth in phantom [mm]; Error MLP - MC [mm]");
   gDifferenceest->SetLineColor(kRed);
   gDifferenceest->SetLineWidth(2);
   gDifferenceest->Draw("AL");
   gDifferenceest->GetYaxis()->SetRangeUser(0, 4);

//   if (divergence > 0) {

      Float_t sigmaNoTrk = (hErrorNaive->GetStdDev(1) + hErrorNaive->GetStdDev(2)) / 2;
      Float_t sigmaEst = (hErrorSigmaScale->GetStdDev(1) + hErrorSigmaScale->GetStdDev(2)) / 2;

      ofstream file(Form("Output/MLPerror_energy%.0fMeV_water_rotation%.0fmrad.csv", initialEnergy, divergence), ofstream::out | ofstream::app);
      file << phantomSize << " " << divergence << " " <<  gDifferenceNoTrk->Eval(0) << " " << gDifferenceNoTrk->Eval(phantomSize/2) << " " << gDifferenceest->Eval(0) << " " << gDifferenceest->Eval(phantomSize/2) << " " << sigmaNoTrk << " " << sigmaEst << endl;
      file.close();
  // }

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
