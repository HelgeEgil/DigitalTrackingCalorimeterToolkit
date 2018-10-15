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

void findMLPLoop(Float_t phantomSize = 200, Float_t spotSize = -1);

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

void findMLPLoop(Float_t phantomSize, Float_t spotSize) {
   Float_t     initialEnergy = 200;
   Float_t     differenceArrayDZ = 3;
   const Int_t eventsToUse = 300;
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
   Double_t    differenceArrayDiff[1000] = {};
   Double_t    differenceArrayDiffNoTrk[1000] = {};
   Double_t    differenceArrayDiffest[1000] = {};
   Int_t       idxWater = 0;
   Double_t    energiesWater[500], rangesWater[500];
   Float_t     energy, range;
   //   Float_t     sigmaFilter = 128.47;
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
      f = new TFile(Form("MC/Output/simpleScanner_energy%.0fMeV_Water_phantom%03.0fmm_spotsize%.3fmm.root", initialEnergy, phantomSize, spotSize)); 
   }

   TTree *tree = (TTree*) f->Get("Hits");

   if (!tree) exit(0);

   in.open("Data/WaterPSTAR.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;
      rangesWater[idxWater] = range*10; // [mm]
      energiesWater[idxWater++] = energy;
   }
   in.close();
   
   TSpline3 *splineWater = new TSpline3("splineWater", energiesWater, rangesWater, idxWater);

   Float_t  AXlow = 0;
   Float_t  AXhigh = 1.05;
   Float_t  APlow = -110;
   Float_t  APhigh = 0;
   
   Float_t  AXdelta = 0.01;
   Float_t  APdelta = 1;

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
   
   P0NoTrk.SetCoordinates(0, 0, 1); // Normalized

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);

      if (eventID > eventsToUse) stop = true;

      if (lastEID < 0) {
         lastEID = eventID;
         sum_edep = edep;
         residualEnergy = 0;
      }

      if (lastEID != eventID) {
         if (lastEID % 10 == 0) printf("Particle %d/%d.\n", lastEID, eventsToUse);
         // New particle, store last values
         if (residualEnergy > sigmaFilter) {
            hResidualEnergy->Fill(residualEnergy);
            // Find vectors as defined in paper
            X0 = (Xp0 + Xp1) / 2;
            X1 = (Xp2 + Xp3) / 2;
            P0 = Xp1 - Xp0;
            P1 = Xp3 - Xp2;
            P0hat = P0.Unit();
            P1hat = P1.Unit();

            X0.SetZ(X0.Z() + 15);
            X1.SetCoordinates(X1.X() - 15 * P1hat.X(), X1.Y() - 15 * P1hat.Y(), X1.Z() - 15);

            // calculate Spline
            wepl = splineWater->Eval(initialEnergy);
            wet  = wepl - splineWater->Eval(residualEnergy);

            // Find lambda values using polynomial in Fig. 4 in Fekete et al. 2015
            Lambda0 = 1.01 + 0.43 * pow(wet/wepl, 2);
            Lambda1 = 0.99 - 0.46 * pow(wet/wepl, 2);
            
            splineMCx  = new TSpline3("splineMCx", arSplineMCz, arSplineMCx, idxSplineMC);
            splineMCy  = new TSpline3("splineMCy", arSplineMCz, arSplineMCy, idxSplineMC);
         
            // INNER MINIMIZATION LOOP
            for (Float_t AX = AXlow; AX <= AXhigh; AX += AXdelta) {
               for (Float_t AP = APlow; AP <= APhigh; AP += APdelta) {
                  X0est = X1 * AX + AP * P1hat; 
                  X0est.SetZ(-phantomSize/2); // We know the Z coordinate...
                  
                  Float_t lastZ = -1;
                  for (Float_t t=0; t<=1; t+= 0.05) { // was 0.01
                     S = SplineMLP(t, X0est, X1, P0NoTrk, P1hat, Lambda0, Lambda1);
                     if (lastZ<1) lastZ = S.Z();
                     diff_x = fabs(splineMCx->Eval(S.Z()) - S.X());
                     diff_y = fabs(splineMCy->Eval(S.Z()) - S.Y());
                     hErrorMatrix->Fill(AX, AP, fabs(S.Z()-lastZ)*sqrt(pow(diff_x, 2) + pow(diff_y, 2)));
                     hIdxMatrix->Fill(AX, AP);
                  }
               }
            }

            delete splineMCx;
            delete splineMCy;
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

         if  (volumeID[2] < 5) {
            arSplineMCx[idxSplineMC] = x;
            arSplineMCy[idxSplineMC] = y;
            arSplineMCz[idxSplineMC++] = z;
         }
      }
   }

   gStyle->SetOptStat(0);
   gStyle->SetNumberContours(99);

   hErrorMatrix->Divide(hIdxMatrix);
  
   Int_t nbins = hErrorMatrix->GetNcells();
   Float_t mincont = 1e5;
   Int_t mincell;
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
