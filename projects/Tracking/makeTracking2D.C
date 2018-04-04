#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <algorithm>
#include <TPolyLine.h>
#include <TSpline.h>
#include <iostream>
#include <fstream>

using namespace std;

Int_t  findNextPoint(Float_t *arrayx, Float_t *arrayy, Float_t *arrayz, Int_t thisIdx, Int_t maxIdx) {
   Float_t  x = arrayx[thisIdx];
   Float_t  y = arrayy[thisIdx];
   Float_t  z = arrayz[thisIdx];
   Float_t  x2, z2, y2;
   Float_t  delta, minDelta = 1e3;
   Int_t    minIdx = -1;
   Bool_t   use3DInformation = false;

   for (int idx = 0; idx < maxIdx; idx++) {
      if (idx == thisIdx) continue;

      x2 = arrayx[idx];
      y2 = arrayy[idx];
      z2 = arrayz[idx];
   
      if (z2 > z && z2 - z < 4) { // next layer
         if (use3DInformation) {
            delta = sqrt(pow(x-x2, 2) + pow(y-y2, 2));
         }
         else {
            delta = fabs(x-x2);
         }
         if (delta < minDelta) {
            minDelta = delta;
            minIdx = idx;
         }
      }
   }

   return minIdx;
}

void makeTracking2D() {
//   TFile    *f = new TFile("Data/DTC_Aluminium_Absorber3mm_Degrader250mm_250MeV.root");
   TFile    *f = new TFile("Data/RealisticPencilBeamThrough20cmWaterPhantom.root");
   TTree    *tree = (TTree*) f->Get("Hits");
   Float_t  x,y,z,edep;
   Int_t    parentID, eventID;
      
   const Int_t arraySize = 10000;
   Int_t    numberOfPoints = 0;
   const Int_t    numberOfTracks = 15;
   Float_t  trackPoints_x[arraySize];
   Float_t  trackPoints_y[arraySize];
   Float_t  trackPoints_z[arraySize];
   Int_t    eventIDs[arraySize];

   TCanvas *canvas = new TCanvas("canvas", "Dose profile", 1200, 600);
   
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("parentID", &parentID); // primary (0) or daughter (>0) particle?
   tree->SetBranchAddress("eventID", &eventID); // unique primary ID

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if (parentID == 0) { // primary particle, we ignore secondaries for now
         if (eventID > numberOfTracks) break;
         eventIDs[numberOfPoints] = eventID;
         trackPoints_x[numberOfPoints] = x;
         trackPoints_y[numberOfPoints] = y;
         trackPoints_z[numberOfPoints++] = z;
      }
   }

   // TGraph (nPoints, x_array, y_array)
   TGraph *tracks2D = new TGraph(numberOfPoints, trackPoints_z, trackPoints_x);
   tracks2D->SetTitle(Form("%d primary tracks, no secondaries;Z depth [mm];X position [mm]", numberOfTracks));
   tracks2D->SetMarkerStyle(21);
   tracks2D->SetMarkerSize(0.8);
   tracks2D->SetMarkerColor(kRed);
   gPad->SetGridx();
   gPad->SetGridy();
   Int_t firstID;
   Int_t correctTracks[numberOfTracks];
   for (Int_t i=0; i<numberOfTracks; i++) {
      correctTracks[i] = 1;
   }

   tracks2D->Draw("AP"); // A = Create axis, P = draw as points (L = lines, * = stars)

   // Tracking
   vector<Int_t> vSeen; // Hits that are already searched
   Int_t nextHit, nLinePoints = 0;

   for (Int_t idx=0; idx<numberOfPoints; idx++) {
      if (trackPoints_z[idx] < 1) { // seed layer
         if ( std::find(vSeen.begin(), vSeen.end(), idx) != vSeen.end() ) { // Have seen this hit before
            continue;
         }
         else { // Have NOT seen this hit before
            firstID = eventIDs[idx];
            TPolyLine *l = new TPolyLine(numberOfPoints);
            nLinePoints = 0;
            l->SetLineWidth(2);
            l->SetPoint(nLinePoints++, trackPoints_z[idx], trackPoints_x[idx]);
            vSeen.push_back(idx);
            nextHit = findNextPoint(trackPoints_x, trackPoints_y, trackPoints_z, idx, numberOfPoints);
            while (nextHit > 0) {
               if ( std::find(vSeen.begin(), vSeen.end(), nextHit) == vSeen.end() ) { // Have NOT seen this hit
                  vSeen.push_back(nextHit);
                  if (eventIDs[nextHit] != firstID) {
                     correctTracks[firstID] = 0;
                     l->SetLineColor(kRed);
                  }
                  l->SetPoint(nLinePoints++, trackPoints_z[nextHit], trackPoints_x[nextHit]);
                  nextHit = findNextPoint(trackPoints_x, trackPoints_y, trackPoints_z, nextHit, numberOfPoints);
               }
               else {
                  nextHit = -1;
               }
            }
         l->Draw();
         }
      }
   }

   Int_t nCorrectTracks = 0;
   for (Int_t i=0; i<numberOfTracks; i++) {
      nCorrectTracks += correctTracks[i];
   }
   printf("Fraction of correct tracks = %.2f %%.\n", 100 * nCorrectTracks / numberOfTracks);
}
