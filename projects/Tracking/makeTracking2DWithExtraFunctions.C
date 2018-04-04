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

Float_t initialEnergy = 160.57;
TSpline3 *rangeSpline;
TSpline3 *energySpline;

void Init() {
   Double_t  range, energy;
   Double_t  ranges[500];
   Double_t  energies[500];
   Int_t    idx = 0;
   ifstream in;

   in.open("Data/energy_range_lut.csv");
   while (1) {
      in >> energy >> range;
      if (!in.good()) break;
      ranges[idx] = range;
      energies[idx++] = energy;
   }
   in.close();
   if (idx>0) printf("Successfully read LUT from CSV file.\n");
   
   rangeSpline = new TSpline3("rangeSpline", energies, ranges, idx);
   energySpline = new TSpline3("energySpline", ranges, energies, idx);
}

Float_t getEnergy(Float_t range) {
   return energySpline->Eval(range);
}

Float_t getEnergyAtDepth(Float_t depth) {
   Float_t range = rangeSpline->Eval(initialEnergy);
   Float_t energy = energySpline->Eval(range - depth);
   return energy;
}

Float_t getRange(Float_t energy) {
   return rangeSpline->Eval(energy);
}


Float_t getRMSScatteringAngle(Float_t atDepth, Float_t slabThickness) {
   Float_t  X0 = 8.897;
   Float_t  proton_mass = 938.27;
   Float_t  a = 0.00184;
   Float_t  p = 1.663;
   Float_t  sumIntegral = 0;
   Int_t    nTerms = 1000;
   Float_t  range = a * pow(initialEnergy, p);
   if (slabThickness > 0)  range = slabThickness;

   Float_t  binThickness = (range) / float(nTerms);
   Float_t  depth, energy, gamma, beta, momentum, pv, RMSangle;
   
   for (Int_t i=1; i<nTerms; i++) {
      depth = atDepth + i * range / nTerms;
      energy = pow((range - depth) / a, 1/p);
      if (energy < 0) { 
         printf("energy < 0!\n");
         sumIntegral *= range / (float(i) * binThickness); // correct for this..?
         break;
      }
      gamma = (initialEnergy + proton_mass) / proton_mass;
      beta = sqrt(1 - pow(gamma, -2));
      momentum = gamma * beta * proton_mass;
      pv = beta * momentum;

      sumIntegral += pow(14.1/pv, 2) * binThickness / X0;
   }

   RMSangle = (1 + 0.11111 * log(range / X0)) * sqrt(sumIntegral) * 1000; // mrad

   return RMSangle;
}      

Int_t  findNextPoint(Float_t *arrayx, Float_t *arrayz, Int_t thisIdx, Int_t maxIdx) {
   Float_t x = arrayx[thisIdx];
   Float_t z = arrayz[thisIdx];
   Float_t x2, z2;
   Float_t delta, minDelta = 1e3;
   Int_t minIdx = -1;

//   printf("Searching for next point from (z,x) = (%.2f, %.2f)\n", z,x);

   for (int idx = 0; idx < maxIdx; idx++) {
      if (idx == thisIdx) continue;

      x2 = arrayx[idx];
      z2 = arrayz[idx];
   
      if (z2 > z && z2 - z < 4) { // next layer
         delta = fabs(x-x2);
         if (delta < minDelta) {
//            if (delta > 10) continue;
            minDelta = delta;
            minIdx = idx;
         }
      }
   }
//   printf("Using best candidate with at (z,x) = (%.2f, %.2f) with delta = %.2f mm and idx = %d.\n", z,x, minDelta, minIdx);

   return minIdx;
}

void makeTracking2D() {
   Init();
//   TFile    *f = new TFile("Data/DTC_Aluminium_Absorber3mm_Degrader250mm_250MeV.root");
   TFile    *f = new TFile("Data/RealisticPencilBeamThrough20cmWaterPhantom.root");
   TTree    *tree = (TTree*) f->Get("Hits");
   Float_t  x,y,z,edep;
   Int_t    parentID, eventID;
      
   const Int_t arraySize = 10000;
   Int_t    numberOfPoints = 0;
   Int_t    numberOfTracks = 10;
   Float_t  trackPoints_x[arraySize];
   Float_t  trackPoints_z[arraySize];

   printf("SOME DEBUG:\n");
   printf("The range of a 160 MeV proton is %.2f mm in Al\n", getRange(160));
   printf("The energy of a 100 mm proton is %.2f MeV\n", getEnergy(100));
   Float_t rmsAngle = getRMSScatteringAngle(0, -1);
   printf("The RMS scattering angle through the complete DTC is %.2f mrad (lateral scattering r*alpha = %.2f mm).\n", rmsAngle, 0.001*getRange(initialEnergy) * rmsAngle);

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
         trackPoints_x[numberOfPoints] = x;
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

   tracks2D->Draw("AP"); // A = Create axis, P = draw as points (L = lines, * = stars)

   // Tracking
   vector<Int_t> vSeen; // Hits that are already searched
   Int_t nextHit, nLinePoints = 0;

   Float_t nextTrack_x[arraySize];
   Float_t nextTrack_z[arraySize];

   for (Int_t idx=0; idx<numberOfPoints; idx++) {
      if (trackPoints_z[idx] < 1) { // seed layer
         if ( std::find(vSeen.begin(), vSeen.end(), idx) != vSeen.end() ) { // Have seen this hit before
            continue;
         }
         else { // Have NOT seen this hit before
            TPolyLine *l = new TPolyLine(numberOfPoints);
            nLinePoints = 0;
            l->SetLineWidth(2);
            l->SetPoint(nLinePoints++, trackPoints_z[idx], trackPoints_x[idx]);
            vSeen.push_back(idx);
            nextHit = findNextPoint(trackPoints_x, trackPoints_z, idx, numberOfPoints);
            while (nextHit > 0) {
               if ( std::find(vSeen.begin(), vSeen.end(), nextHit) == vSeen.end() ) { // Have NOT seen this hit
                  vSeen.push_back(nextHit);
                  l->SetPoint(nLinePoints++, trackPoints_z[nextHit], trackPoints_x[nextHit]);
                  nextHit = findNextPoint(trackPoints_x, trackPoints_z, nextHit, numberOfPoints);
               }
               else {
                  nextHit = -1;
               }
            }
         l->Draw();
         }
      }
   }
}
