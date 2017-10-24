#define DataInterface_cxx

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>

#include <TEllipse.h>
#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TView.h>
#include <TLeaf.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/RangeAndEnergyCalculations.h"
#include "Classes/Hit/Hit.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
// #include "Classes/Track/conversionFunctions.h"
#include "HelperFunctions/Tools.h"
#include "Classes/Layer/Layer.h"
#include "Classes/DataInterface/DataInterface.h"
#include "Classes/CalorimeterFrame/CalorimeterFrame.h"
#include "RootFiles/LinkDef.h"

using namespace DTC;

DataInterface::DataInterface(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("Hits",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
     char *materialChar = getMaterialChar();

     Float_t readoutAbsorber = kAbsorberThickness;
     if (roundf(readoutAbsorber) != readoutAbsorber) readoutAbsorber *= 10; // To avoid having int / float mismatches in filenames

     //
      TChain * chain = new TChain("Hits","");
      if (!useDegrader) {
         chain->Add(Form("Data/MonteCarlo/DTC_%s_%.0fMeV_%.0fmm.root/Hits", materialChar, run_energy, kAbsorberThickness));
      }
      else {
         printf("Opening file with degrader thickness %.0f mm, material is %s and abs. thicknesss %.0f mm.\n", run_degraderThickness, materialChar, readoutAbsorber);
         chain->Add(Form("Data/MonteCarlo/DTC_%s_Absorber%.0fmm_Degrader%.0fmm_250MeV.root/Hits", materialChar, readoutAbsorber, run_degraderThickness)); // Fix if original run_energy is != 250
      }
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

DataInterface::~DataInterface()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DataInterface::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DataInterface::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DataInterface::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).


   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PDGEncoding", &PDGEncoding, &b_PDGEncoding);
   fChain->SetBranchAddress("trackID", &trackID, &b_trackID);
   fChain->SetBranchAddress("parentID", &parentID, &b_parentID);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("edep", &edep, &b_edep);
   fChain->SetBranchAddress("stepLength", &stepLength, &b_stepLength);
   fChain->SetBranchAddress("posX", &posX, &b_posX);
   fChain->SetBranchAddress("posY", &posY, &b_posY);
   fChain->SetBranchAddress("posZ", &posZ, &b_posZ);
   fChain->SetBranchAddress("localPosX", &localPosX, &b_localPosX);
   fChain->SetBranchAddress("localPosY", &localPosY, &b_localPosY);
   fChain->SetBranchAddress("localPosZ", &localPosZ, &b_localPosZ);
   fChain->SetBranchAddress("baseID", &baseID, &b_baseID);
   fChain->SetBranchAddress("level1ID", &level1ID, &b_level1ID);
   fChain->SetBranchAddress("level2ID", &level2ID, &b_level2ID);
   fChain->SetBranchAddress("level3ID", &level3ID, &b_level3ID);
   fChain->SetBranchAddress("level4ID", &level4ID, &b_level4ID);
   fChain->SetBranchAddress("layerID", &layerID, &b_layerID);
   fChain->SetBranchAddress("photonID", &photonID, &b_photonID);
   fChain->SetBranchAddress("nPhantomCompton", &nPhantomCompton, &b_nPhantomCompton);
   fChain->SetBranchAddress("nCrystalCompton", &nCrystalCompton, &b_nCrystalCompton);
   fChain->SetBranchAddress("nPhantomRayleigh", &nPhantomRayleigh, &b_nPhantomRayleigh);
   fChain->SetBranchAddress("nCrystalRayleigh", &nCrystalRayleigh, &b_nCrystalRayleigh);
   fChain->SetBranchAddress("primaryID", &primaryID, &b_primaryID);
   fChain->SetBranchAddress("sourcePosX", &sourcePosX, &b_sourcePosX);
   fChain->SetBranchAddress("sourcePosY", &sourcePosY, &b_sourcePosY);
   fChain->SetBranchAddress("sourcePosZ", &sourcePosZ, &b_sourcePosZ);
   fChain->SetBranchAddress("sourceID", &sourceID, &b_sourceID);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("axialPos", &axialPos, &b_axialPos);
   fChain->SetBranchAddress("rotationAngle", &rotationAngle, &b_rotationAngle);
   fChain->SetBranchAddress("volumeID", volumeID, &b_volumeID);
   fChain->SetBranchAddress("processName", processName, &b_processName);
   fChain->SetBranchAddress("comptVolName", comptVolName, &b_comptVolName);
   fChain->SetBranchAddress("RayleighVolName", RayleighVolName, &b_RayleighVolName);
   Notify();
}

Bool_t DataInterface::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DataInterface::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DataInterface::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void DataInterface::getMCData(Int_t runNo, TH3F* Frame3D) {
   if (fChain==0) return;

   Int_t eventIdFrom = runNo * kEventsPerRun/10;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun/10;

   Float_t offsetX = (nx/2+2) * dx;
   Float_t offsetY = (ny/2) * dy;
   Float_t x,y;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nb = 0;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      nb = fChain->GetEntry(jentry);

      if (eventID < eventIdFrom) continue;
      else if (eventID >= eventIdTo) break;

      x = (posX + offsetX) * nx/2 / (offsetX);
      y = (posY + offsetY) * ny/2 / (offsetY);

      Frame3D->Fill(posZ, x, y, edep*1000/14); // keV / um
   }
   
} // end function GetData3D

void DataInterface::getEventIDs(Int_t runNo, Hits * hits) {
   if (fChain==0) return;

   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

   Int_t lastZ = -1;

   Float_t offsetX = (nx/2+2) * dx;
   Float_t offsetY = (ny/2) * dy;
   Float_t x,y,z;
   Float_t xS = 0, yS = 0, edepS = 0;
   Int_t n = 0;
   Int_t lastEventID = -1;

   Float_t xAvg, yAvg;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nb = 0;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      nb = fChain->GetEntry(jentry);

      if (eventID < eventIdFrom) continue;
      else if (eventID >= eventIdTo) break;

      if (parentID > 0) continue;

      x = (posX + offsetX) * nx/2 / (offsetX);
      y = (posY + offsetY) * ny/2 / (offsetY);
      z = level1ID + baseID; // makeGeometry.py

      // often the layer has more than one hit
      // and they must be accumulated before storing
      // (sum edep, average position)

      if (lastZ < 0 || (lastZ == z && lastEventID == eventID )) {
         // First hit in list OR hit in same layer
         xS += x;
         yS += y;
         edepS += edep;
         n++;
         lastZ = z;
         lastEventID = eventID;
//       cout << "First/same layer. x = " << x << ", y = " << y << ", edepS = " << edepS << ", z = " << z << ", eventID = " << eventID << ", n = " << n << endl;
      }

      else {
         // new layer
         xAvg = xS / n;
         yAvg = yS / n;
         hits->appendPoint(xAvg, yAvg, lastZ, lastEventID, edepS/14);
//       cout << "Saving old layer: x = " << xS << "/" << n << " = " << xAvg << ", y = " << yAvg << ", edep = " << edepS << ", z = " << lastZ << ", eventID = " << lastEventID << endl;

         xS = x;
         yS = y;
         edepS = edep;
         lastZ = z;
         lastEventID = eventID;
         n = 1;
//       cout << "New layer: x = " << x << ", y = " << y << " Z = " << z << ", eventID = " << eventID << endl;
      }
   }
   // last layer
   xAvg = xS / n;
   yAvg = yS / n;
   hits->appendPoint(xAvg, yAvg, lastZ, lastEventID, edepS/14);
// cout << "Saving last layer: x = " << xS << "/" << n << " = " << xAvg << ", y = " << yAvg << ", edep = " << edepS << ", z = " << lastZ << ", eventID = " << lastEventID << endl;

}

void DataInterface::getDataProfile(TH2F *hProfile, TH2F *hProjection, Int_t energy) {
   if (!existsEnergyFile(energy)) {
      cout << "There are no data files with energy " << energy << endl;
      return;
   }

   TString fn = Form("Data/ExperimentalData/DataFrame_%i_MeV.root", energy);
   TFile *f = new TFile(fn);
   TTree *tree = (TTree*) f->Get("tree");

   Int_t nentries = tree->GetEntries();

   TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
   TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
   TLeaf *lLayer = tree->GetLeaf("fDataFrame.fLayer");

   Float_t x, y, layer;

   for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);

      for (Int_t j=0; j<lY->GetLen(); j++) {
         x = lX->GetValue(j)  + nx/2;
         y = lY->GetValue(j)  + ny/2;
         layer = lLayer->GetValue(j);
         
         hProfile->Fill(y, layer);
         hProjection->Fill(x, y);
      }
   }
}

void DataInterface::getDataFrame(Int_t runNo, CalorimeterFrame * cf, Int_t energy) {

   if (!existsEnergyFile(energy)) {
      cout << "There are no data files with energy " << energy << endl;
      return;
   }

   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

   TString fn = Form("Data/ExperimentalData/DataFrame_%i_MeV.root", energy);
   TFile *f = new TFile(fn);
   TTree *tree = (TTree*) f->Get("tree");

   Int_t nentries = tree->GetEntries();
   if (eventIdTo > nentries) {
      eventIdTo = nentries;
   }
   cout << "Found " << nentries << " frames in the DataFrame.\n";

   TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
   TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
   TLeaf *lLayer = tree->GetLeaf("fDataFrame.fLayer");

   Int_t counter = 0;
   for (Int_t i=eventIdFrom; i<eventIdTo; i++) {
      tree->GetEntry(i);

      for (Int_t j=0; j<lX->GetLen(); j++) {
         Int_t x = lX->GetValue(j) + nx/2;
         Int_t y = lY->GetValue(j) + ny/2;
         Int_t z = lLayer->GetValue(j);

         if ( x > nx || y > ny ) printf("POINT (x,y,z) = (%d\t%d\t%d) OUT OF BOUNDS!!\n", x,y,z);
         cf->fillAt(z, x, y);

      }
      counter++;
   }
   delete f;
}

void DataInterface::writeDataFrame(Int_t energy) {
   if (!existsEnergyFile(energy)) {
      cout << "There are no data files with energy " << energy << endl;
   }
   
   TString fn = Form("Data/ExperimentalData/DataFrame_%i_MeV.root", energy);
   TFile *f = new TFile(fn);
   TTree *tree = (TTree*) f->Get("tree");

   Int_t nentries = tree->GetEntries();
   cout << "Found " << nentries << " frames in the DataFrame.\n";
  
   ofstream iofile(Form("OutputFiles/DataFrameCSV_%i_MeV.csv", energy), ofstream::out);

   TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
   TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
   TLeaf *lLayer = tree->GetLeaf("fDataFrame.fLayer");

   Int_t counter = 0;
   for (Int_t i=0; i<nentries; i++) {
      tree->GetEntry(i);

      for (Int_t j=0; j<lX->GetLen(); j++) {
         Int_t x = lX->GetValue(j) + nx/2;
         Int_t y = lY->GetValue(j) + ny/2;
         Int_t z = lLayer->GetValue(j);

         iofile << i << " " << x << " " << y << " " << z << endl;
      }
   }
   delete f;
}


void  DataInterface::getMCClusters(Int_t runNo, Clusters *clusters) {
   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

   printf("eventIdFrom = %d, eventIdTo = %d.\n", eventIdFrom, eventIdTo);

   if (runNo == 0) lastJentry_ = 0;
   
   Float_t  sum_edep = 0;
   Int_t    lastID = 0;
   Int_t    lastLayer = 0;
   Float_t  lastZ = 0;
   Int_t    layer;
   Int_t    n = 0;
   Float_t  sumX = 0, sumY = 0;
   Float_t  x,y;
   Bool_t   isInElastic;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=lastJentry_; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry<0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

      layer = level1ID + baseID - 1;
      if (parentID != 0) continue;
      
      if (lastID != eventID || lastLayer != layer) {
         x = sumX/n / dx + nx/2;
         y = sumY/n / dy + ny/2;
         
         clusters->appendClusterEdep(x, y, lastLayer, sum_edep/14, lastID);

         sum_edep = 0;
         sumY = 0;
         sumX = 0;
         n = 0;
      }
      
         
      if (eventID < eventIdFrom) {
         printf("eventID (%d) < eventIdFrom (%d), continuing.\n", eventID, eventIdFrom);
         continue;
      }
      
      else if (eventID >= eventIdTo) {
         printf("eventID (%d) >= eventIdTo (%d), breaking.\n", eventID, eventIdTo);
         lastJentry_ = jentry;
         break;
      }

      sum_edep += edep*1000;
      sumX += posX;
      sumY += posY;
      lastZ = posZ;
      lastID = eventID;
      lastLayer = layer;
      n++;
   }

}

Int_t DataInterface::getMCFrame(Int_t runNo, CalorimeterFrame *cf, Float_t *x_energy, Float_t *y_energy) {
   // Retrieve kEventsPerRun events and put them into CalorimeterFrame*
   // if x_energy and y_energy is supplied, it implies that a full simulation is performed, and hits are recorded everywhere
   
   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

   if (runNo == 0) lastJentry_ = 0;

   Float_t offsetX = (nx/2+2) * dx;
   Float_t offsetY = (ny/2) * dy;
   Float_t x,y;
   Float_t sum_edep = 0;
   Int_t n = 0;
   Int_t calorimeterLayer = 0;
   Int_t blackListEventID = -1;
   Int_t whiteListEventID = -1;
   Int_t lastID = 0;

   if (fChain == 0) return -1;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=lastJentry_; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry<0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (lastID != eventID) {
         sum_edep = 0;
         n = 0;
      }
      
      sum_edep += edep/14;
      
      if  (x_energy && y_energy && parentID == 0) {
         x_energy[n + sizeOfEventID*eventID] = posZ + 1.5; // 1.5 for Al absorber
         y_energy[n + sizeOfEventID*eventID] = run_energy - sum_edep;
         n++;
      }
      
//    By using the 'old' geometry (used in the 2016 NIMA publicatioN), volumeID[4] == 4 in a full simulation designated a chip
//    In the geometry generated by makeGeometry.py, it is a bit more complicated... But for simulations only using readout in chips, it is of no consequence

//    if (volumeID[4] == 4 || !x_energy) { // hit in chip, if no x_energy provided do this
      if (!x_energy) { // hit in chip, if no x_energy provided do this
         
         if (eventID < eventIdFrom) {
            continue;
         }
         
         else if (eventID >= eventIdTo) {
            lastJentry_ = jentry;
            break;
         }
         
         if (eventID == blackListEventID) continue;

         if (posZ < -100) { // scintillator hits
            blackListEventID = eventID;
            whiteListEventID = eventID;
            continue;
         }

         // In the geometry generated by makeGeometry, baseID==0 + level1ID==0 implies first layer, and baseID==0 + level1ID==0 is 2nd layer
         calorimeterLayer = level1ID + baseID; 

         if (calorimeterLayer<0 || posZ < -20) {
            continue;
         }

         x = (posX + offsetX) * nx/2 / (offsetX);
         y = (posY + offsetY) * ny/2 / (offsetY);

         cf->fillAt(calorimeterLayer, x, y, edep*1000/14);
      }
      
      lastID = eventID;
   }

   if (kEventsPerRun > 1) return -1;
   return eventID;
}
