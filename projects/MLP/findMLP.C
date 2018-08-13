#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <fstream>
#include <TGraph.h>
#include <TSpline.h>
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

void findMLP(Float_t sigmaScaleFactor = 1, Float_t vectorScaleFactor = 0) {
   TFile *f = new TFile("MC/simpleScanner.root");
   TTree *tree = (TTree*) f->Get("Hits");

   Float_t     initialEnergy = 230;
   Float_t     x, y, z, edep, sum_edep = 0, residualEnergy = 0;
   Int_t       eventID, parentID, lastEID = -1;
   XYZVector   Xp0, Xp1, Xp2, Xp3, X0, X1, X0est, P0, P1, P0hat, P1hat, S; // Xp are the plane coordinates, X are the tracker coordinates (X1 = (Xp1 + Xp2) / 2)
   Float_t     wepl, wet, Lambda0, Lambda1;
   ifstream    in;
   Int_t       idxWater = 0;
   Double_t    energiesWater[500], rangesWater[500];
   Float_t     energy, range;
   Float_t     sigmaFilter = 128.47;
   Double_t    aPosMCx[100000];
   Double_t    aPosMCz[100000];
   Double_t    aPosMLPx[100000];
   Double_t    aPosMLPz[100000];
   Int_t       aIdxMC = 0;
   Int_t       aIdxMLP = 0;
   Bool_t      stop = false;

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

   TH1F *hResEnergy = new TH1F("hResEnergy", "Residual energy in calorimeter;Energy [MeV];Entries", 300, 0, 150);
   TH2F *hErrorNaive = new TH2F("hErrorNaive", "Beamspot uncertainty (assume point beam);X position [mm];Y position [mm]", 200, -10, 10, 200, -10, 10);
   TH2F *hErrorSigmaScale = new TH2F("hErrorSigmaScale", Form("Beamspot uncertainty (use X1 * %.2f + P1 * %.2f);X position [mm];Y position [mm]", sigmaScaleFactor, vectorScaleFactor), 200, -10, 10, 200, -10, 10);

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
//      if (eventID > 1000) stop = true;

      tree->GetEntry(i);

      if (lastEID < 0) {
         lastEID = eventID;
         sum_edep = edep;
         residualEnergy = 0;
      }

      if (lastEID != eventID) {
         // New particle, store last values
         if (residualEnergy > sigmaFilter) {
            hResEnergy->Fill(residualEnergy);

            // Find vectors as defined in paper
            X0 = (Xp0 + Xp1) / 2;
            X1 = (Xp2 + Xp3) / 2;
            P0 = Xp1 - Xp0;
            P1 = Xp3 - Xp2;
            P0hat = P0.Unit();
            P1hat = P1.Unit();

            X0est = X1 * sigmaScaleFactor + P1 * vectorScaleFactor - X0;

            hErrorNaive->Fill(X0.X(), X0.Y());
            hErrorSigmaScale->Fill(X0est.X(), X0est.Y());
            
            // calculate Spline
            wepl = splineWater->Eval(initialEnergy);
            wet = splineWater->Eval(initialEnergy) - splineWater->Eval(residualEnergy);

            // Find lambda values using polynomial in Fig. 4 in Fekete et al. 2015
            Lambda0 = 1.01 + 0.43 * pow(wet/wepl, 2);
            Lambda1 = 0.99 - 0.46 * pow(wet/wepl, 2);

            for (Float_t t=0; t<1; t += 0.01) {
               S = SplineMLP(t, X0, X1, P0hat, P1hat, Lambda0, Lambda1);
               if (eventID < 10) {
                  aPosMLPx[aIdxMLP] = S.X();
                  aPosMLPz[aIdxMLP++] = S.Z();
               }
            }
         }
         
         if (stop) break;

         // Reset counters for next primary
         sum_edep = edep;
         residualEnergy = 0;
         lastEID = eventID;
      }

      else { // Still following the same particle
         sum_edep += edep;
      }


      if (parentID == 0) {
         if       (z < -119 && z > -121) { Xp0.SetCoordinates(x,y,z); }
         else if  (z < -109 && z > -111) { Xp1.SetCoordinates(x,y,z); }
         else if  (z >  109 && z <  111) { Xp2.SetCoordinates(x,y,z); }
         else if  (z >  119 && z <  121) { Xp3.SetCoordinates(x,y,z); }
         else if  (z >  130)             { residualEnergy += edep; }
         if (z < 130 && eventID < 10) {
            aPosMCx[aIdxMC] = x;
            aPosMCz[aIdxMC++] = z;
         }
      }
   }

   hResEnergy->Draw();

   TGraph *gMC = new TGraph(aIdxMC, aPosMCz, aPosMCx);
   TGraph *gMLP = new TGraph(aIdxMLP, aPosMLPz, aPosMLPx);
   gMC->SetMarkerStyle(7);
   gMLP->SetMarkerStyle(7);
   gMC->SetMarkerColor(kRed);
   gMLP->SetMarkerColor(kBlue);
   gMC->Draw("AP");
   gMLP->Draw("P");
   gMC->SetTitle("Optimized Spline MLP estimation (no nucl. interactions);Depth [mm];Lateral position [mm]");

   TCanvas *c = new TCanvas("c", "Beam spot estimation", 800, 600);
   c->Divide(2, 1, 0.001, 0.001);

   c->cd(1);
   hErrorNaive->Draw("COLZ");
   c->cd(2);
   hErrorSigmaScale->Draw("COLZ");
}
