#ifndef DataInterface_cxx
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

DataInterface::DataInterface(TTree *tree) : fChain(0) {
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
      if (!kUseDegrader) {
//         chain->Add(Form("Data/MonteCarlo/DTC_%s_%.0fMeV_%.0fmm.root/Hits", materialChar, run_energy, kAbsorberThickness));
         chain->Add(Form("Data/MonteCarlo/phantom_%.0fMeV_highacc3.root/Hits", run_energy));
      }
      else if (kFinalDesign) { // FINAL DESIGN
         if (!kHelium) { // FINAL DESIGN : PROTONS
            if (kSpotScanning) { // FINAL DESIGN : PROTONS : SPOT SCANNING
               chain->Add(Form("Data/MonteCarlo/DTC_Final_HeadPhantom_rotation%ddeg.root/Hits", kRotation));
               printf("Opening PROTON phantom file with spotX %04.0f and %d degrees rotation\n", kSpotX, kRotation);
            }
            else { // FINAL DESIGN : PROTONS : SINGLE PENCIL BEAM
               printf("Opening PROTON file with degrader thickness %.0f mm and FINAL design (3.5 mm)\n", run_degraderThickness);
               chain->Add(Form("Data/MonteCarlo/DTC_Final_Degrader%03.0fmm_%dMeV.root/Hits", run_degraderThickness, kEnergy));
            }
         }
         else { // FINAL DESIGN : HELIUM : SINGLE PENCIL BEAM
            printf("Opening HELIUM file with degrader thickness %.0f mm and FINAL design (3.5 mm)\n", run_degraderThickness);
            chain->Add(Form("Data/MonteCarlo/DTC_Final_Helium_Degrader%03.0fmm_%dMeV.root/Hits", run_degraderThickness, kEnergy));
         }
      }

      else { // OLD DESIGN
         if (!kHelium) {
            if (!kSpotScanning) {
               printf("Opening PROTON file with degrader thickness %.0f mm, material is %s and abs. thickness %.0f mm...\n", run_degraderThickness, materialChar, readoutAbsorber);
               chain->Add(Form("Data/MonteCarlo/DTC_%s_Absorber%.0fmm_Degrader%03.0fmm_%dMeV.root/Hits", materialChar, readoutAbsorber, run_degraderThickness, kEnergy)); // Fix if original run_energy is != 250
            }
         }
         else {
            printf("Opening HELIUM file with degrader thickness %.0f mm, material is %s and abs. thickness %.0f mm...\n", run_degraderThickness, materialChar, readoutAbsorber);
            if (!kSpotScanning) {
               chain->Add(Form("Data/MonteCarlo/DTC_%s_Helium_Absorber%.0fmm_Degrader%03.0fmm_%dMeV.root/Hits", materialChar, readoutAbsorber, run_degraderThickness, kEnergy));
            }
            else { // Load phantom (used with spot scanning)
               chain->Add(Form("Data/MonteCarlo/DTC_%s_CarbonHelium_Absorber%.0fmm_Headphantom_%dMeV.root/Hits", materialChar, readoutAbsorber, kEnergy));
               primariesInSpot_ = 0;
            }

            printf("OK!\n");
         }
      }
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

DataInterface::~DataInterface() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DataInterface::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t DataInterface::LoadTree(Long64_t entry) {
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

void DataInterface::Init(TTree *tree) {
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
//   fChain->SetBranchAddress("trackID", &trackID, &b_trackID);
   fChain->SetBranchAddress("parentID", &parentID, &b_parentID);
//   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("edep", &edep, &b_edep);
//   fChain->SetBranchAddress("stepLength", &stepLength, &b_stepLength);
   fChain->SetBranchAddress("posX", &posX, &b_posX);
   fChain->SetBranchAddress("posY", &posY, &b_posY);
   fChain->SetBranchAddress("posZ", &posZ, &b_posZ);
//   fChain->SetBranchAddress("localPosX", &localPosX, &b_localPosX);
//   fChain->SetBranchAddress("localPosY", &localPosY, &b_localPosY);
//   fChain->SetBranchAddress("localPosZ", &localPosZ, &b_localPosZ);
   fChain->SetBranchAddress("baseID", &baseID, &b_baseID);
   fChain->SetBranchAddress("level1ID", &level1ID, &b_level1ID);
//   fChain->SetBranchAddress("level2ID", &level2ID, &b_level2ID);
//   fChain->SetBranchAddress("level3ID", &level3ID, &b_level3ID);
//   fChain->SetBranchAddress("level4ID", &level4ID, &b_level4ID);
//   fChain->SetBranchAddress("layerID", &layerID, &b_layerID);
//   fChain->SetBranchAddress("photonID", &photonID, &b_photonID);
//   fChain->SetBranchAddress("nPhantomCompton", &nPhantomCompton, &b_nPhantomCompton);
//   fChain->SetBranchAddress("nCrystalCompton", &nCrystalCompton, &b_nCrystalCompton);
//   fChain->SetBranchAddress("nPhantomRayleigh", &nPhantomRayleigh, &b_nPhantomRayleigh);
//   fChain->SetBranchAddress("nCrystalRayleigh", &nCrystalRayleigh, &b_nCrystalRayleigh);
//   fChain->SetBranchAddress("primaryID", &primaryID, &b_primaryID);
//   fChain->SetBranchAddress("sourcePosX", &sourcePosX, &b_sourcePosX);
//   fChain->SetBranchAddress("sourcePosY", &sourcePosY, &b_sourcePosY);
//   fChain->SetBranchAddress("sourcePosZ", &sourcePosZ, &b_sourcePosZ);
//   fChain->SetBranchAddress("sourceID", &sourceID, &b_sourceID);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   if (kSpotScanning) {
      fChain->SetBranchAddress("spotPosX", &spotPosX, &b_spotPosX);
      fChain->SetBranchAddress("spotPosY", &spotPosY, &b_spotPosY);
   }

//   fChain->SetBranchAddress("runID", &runID, &b_runID);
//   fChain->SetBranchAddress("axialPos", &axialPos, &b_axialPos);
//   fChain->SetBranchAddress("rotationAngle", &rotationAngle, &b_rotationAngle);
//   fChain->SetBranchAddress("volumeID", volumeID, &b_volumeID);
//   fChain->SetBranchAddress("processName", processName, &b_processName);
//   fChain->SetBranchAddress("comptVolName", comptVolName, &b_comptVolName);
//   fChain->SetBranchAddress("RayleighVolName", RayleighVolName, &b_RayleighVolName);
   Notify();
}

Bool_t DataInterface::Notify() {
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DataInterface::Show(Long64_t entry) {
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t DataInterface::Cut(Long64_t entry) {
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
         hits->appendPoint(xAvg, yAvg, lastZ, edepS/14, lastEventID);
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
   hits->appendPoint(xAvg, yAvg, lastZ, edepS/14, lastEventID);
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

void DataInterface::getDataFrame(Int_t runNo, Layer * l, Int_t energy) {

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

   for (Int_t i=eventIdFrom; i<eventIdTo; i++) {
      tree->GetEntry(i);

      for (Int_t j=0; j<lX->GetLen(); j++) {
         Int_t x = lX->GetValue(j) + nx/2;
         Int_t y = lY->GetValue(j) + ny/2;
         Int_t z = lLayer->GetValue(j);

         if ( x > nx || y > ny ) printf("POINT (x,y,z) = (%d\t%d\t%d) OUT OF BOUNDS!!\n", x,y,z);

         if (l->getLayer() == z) { // Only fill the single layer in question
            l->Fill(x, y);
         }
      }
   }
   delete f;
}

void DataInterface::getDataHits(Int_t runNo, Hits * hits, Int_t energy) {
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

   for (Int_t i=eventIdFrom; i<eventIdTo; i++) {
      tree->GetEntry(i);

      for (Int_t j=0; j<lX->GetLen(); j++) {
         Int_t x = lX->GetValue(j) + nx/2;
         Int_t y = lY->GetValue(j) + ny/2;
         Int_t z = lLayer->GetValue(j);

         if ( x > nx || y > ny ) {
            printf("POINT (x,y,z) = (%d\t%d\t%d) OUT OF BOUNDS!!\n", x,y,z);
         }

         hits->appendPoint(x,y,z,i);
      }
   }
   delete f;
}

void  DataInterface::getMCClustersThreshold(Int_t runNo, Clusters *clusters, Hits *hits, Float_t useSpotX, Float_t useSpotY) {
   Int_t eventIdFrom = runNo * kEventsPerRun + kSkipTracks;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun + kSkipTracks;
   Int_t lastPropagated = -1;
   Float_t alpideThickness = 25; // was 14

   // let = (n/4.23)^1/0.65 [Pettersen et al. Phys Med 2019]
   // n < 2 -> let < (2/4.23)^1/0.65 = 0.316 keV/um
   // edep = let * 14 um = 0.316 kev/um * 14 um = 4.4e-3 MeV ~ 4e-3 MeV
   // new design: let * 50 um = 0.316 keV/um * 50 um = 15.8e-3 MeV ~ 1.5e-2 MeV
   // new design 2: let * 25 um = 0.316 keV/um * 25 um = 7.9e-3 MeV ~ 8e-3 MeV

   Float_t threshold = 8e-3;
   Int_t particlesBelowThreshold = 0;

   if (runNo == 0 && !kSpotScanning) lastJentry_ = 0;

// Split the file instead (to save memory when loading 40x copies...)
   if (lastJentry_ == 0 && kSpotScanning) {
      lastJentry_ = findSpotIndex(useSpotX);
   }


   Int_t    layer = 0, lastEventID=-1;
   Float_t  x,y;
   Bool_t   isSecondary = false;
   Long64_t nentries = fChain->GetEntriesFast();

   for (Long64_t jentry=lastJentry_; jentry<nentries; jentry++) { // new interaction
      Long64_t ientry = LoadTree(jentry);
      if (ientry<0) {
         lastJentry_ = jentry;
         break;
      }

      fChain->GetEntry(jentry);
      if (lastEventID != eventID) {
         primariesInSpot_++;
         isSecondary = false;
      }

      if (!kSpotScanning) {
         if (eventID < eventIdFrom) {
            continue;
         }
         
         else if (eventID >= eventIdTo || jentry == nentries-1) { // LAST ENTRY, EXIT AND APPEND THIS
            lastJentry_ = jentry;
            break;
         }
      }

      else {
         if (useSpotX > spotPosX || useSpotY > spotPosY) {
            continue;
         }

         else if (useSpotX < spotPosX || useSpotY < spotPosY || jentry == nentries-1) { // new spot in data -> Store & exit
            lastJentry_ = jentry;
            break;
         }
         else {
            if (primariesInSpot_ >= kEventsPerRun) { // same spot, new Run
               primariesInSpot_ = 0;
               lastJentry_ = jentry;
               break;
            }
         }
      }
      
      layer = level1ID + baseID - 2;

      if (kFinalDesign) {
         if (!kPhantom) {
            layer = baseID*2 + level1ID - 2;
         }

         else { // Re-ordering of baseID and level1ID due to phantom
//            layer = 2 * level1ID + baseID;
            layer = baseID*2 + level1ID;
         }
      }
      
      if (posZ < 0) { // Inside degrader, check if nuclear interactions
         
         if (TString(processName) == "xxhadElastic" || TString(processName) == "alphaInelastic" || TString(processName) == "protonInelastic") {
            isSecondary = true;
         }

         lastEventID = eventID;
         continue;
      }

      if (parentID != 0 || TString(processName) == "alphaInelastic" || TString(processName) == "protonInelastic") { // Secondary track
         isSecondary = true;

         /*
         if (clusters) {
            clusters->propagateSecondaryStatusFromTop();
            }

         if (hits) {
            hits->propagateSecondaryStatusFromTop();
         }
         */

         // If secondary particle is not electron or gamma, propagate secondary status to eventID-particle
         if (clusters && lastPropagated != eventID) {
            if (clusters->Last()) { // This would give segfault sometimes
               if (clusters->Last()->getEventID() == eventID && PDGEncoding > 1000) {
                  clusters->Last()->setSecondary(true);
//                  clusters->propagateSecondaryStatusFromTop(eventID);
                  lastPropagated = eventID;
               }
            }
         }

         if (hits && lastPropagated != eventID && false) {
            if (hits->Last()) { // This would give segfault sometimes
               if (hits->Last()->getEventID() == eventID && PDGEncoding > 1000) {
                  hits->Last()->setSecondary(true);
//                  hits->propagateSecondaryStatusFromTop(eventID);
                  lastPropagated = eventID;
               }
            }
         }
      }
      
      if (edep < threshold) {
         particlesBelowThreshold++;
         lastEventID = eventID;
         continue;
      }

      x = posX / dx + nx/2;
      y = posY / dy + ny/2;

/*      
      printf("VolumeIDs layer %d posz %.3f: ", layer, posZ);
      for (Int_t i=0; i<10; i++) {
         printf("%d = %d; ", i, volumeID[i]);
      }
      printf("\n");
  */    
//      printf("posz %.3f baseID %d level1ID %d -> layer %d. edep %.2f keV/um, PDG %d, parentID %d\n", posZ, baseID, level1ID, layer, edep/alpideThickness*1000, PDGEncoding, parentID);


//      if (volumeID[3] == 0 && volumeID[4] == 0 && volumeID[5] == -1) {
         if (layer < nLayers) {
            if (hits)      hits->appendPoint(x,y,layer,edep*1000/alpideThickness,eventID,isSecondary,PDGEncoding);
            if (clusters)  clusters->appendClusterEdep(x,y,layer,edep*1000/alpideThickness,eventID,isSecondary,PDGEncoding);
         }
//      }

      lastEventID = eventID;
   }
}

Long64_t DataInterface::findSpotIndex(Float_t findSpotX) {
   // Binary search of input ROOT file to file the first occurence of spot X
   // The file is of course sorted with increasing spotX values
   
   Long64_t nentries = fChain->GetEntriesFast(), firstIndex, maxIndex, nextIndex;
   nextIndex = nentries / 2;
   maxIndex = nentries;
   firstIndex = 0;

   fChain->GetEntry(0);
   if (spotPosX == findSpotX) return 0;

   Int_t nTries = 0;
   while (nTries < 200) {
      nTries++;

      if (maxIndex - firstIndex == 1) break;

      Int_t ientry = fChain->LoadTree(nextIndex);
      fChain->GetEntry(nextIndex);

      if (spotPosX >= findSpotX || ientry<0) {
         maxIndex = nextIndex;
         nextIndex = (firstIndex + maxIndex) / 2; 
      }

      else if (spotPosX < findSpotX) {
         firstIndex = nextIndex;
         nextIndex = (firstIndex + maxIndex) / 2;
      }
   }
   return maxIndex;
}

void  DataInterface::getMCClusters(Int_t runNo, Clusters *clusters, Hits * hits, Float_t useSpotX, Float_t useSpotY) {
   Int_t eventIdFrom = runNo * kEventsPerRun + kSkipTracks;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun + kSkipTracks;

   if (runNo == 0 && !kSpotScanning) lastJentry_ = 0;

   /* 
   if (kSpotScanning) {
      printf("New run... Searching for (%.0f,%.0f), starting with jentry = %lld\n", useSpotX, useSpotY, lastJentry_);
   }
   */

   Float_t  sum_edep = 0;
   Int_t    lastEventID = -1, lastParentID = -1, lastTrackID = -1;
   Int_t    lastLayer = 0;
   Float_t  lastZ = 0;
   Int_t    layer = 0;
   Int_t    n = 0, nAdded = 0;
   Float_t  sumX = 0, sumY = 0;
   Float_t  x,y;
   Bool_t   isInElastic, isSecondary = false;
   Char_t   lastProcessName[25];
   Int_t    lastPDG = 0;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=lastJentry_; jentry<nentries; jentry++) { // new interaction
      Long64_t ientry = LoadTree(jentry);
      if (ientry<0) {
         lastJentry_ = jentry;
         break;
      }

      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (kSpotScanning) {
         if (useSpotX > spotPosX || useSpotY > spotPosY) {
//            if (jentry%1000==0) printf("jentry = %lld\n", jentry)   
            lastEventID = -1;
            continue;
         }
         else if (useSpotX < spotPosX || useSpotY < spotPosY || jentry == nentries-1) { // new spot in data -> Store & exit
            primariesInSpot_ = 0;
            lastJentry_ = jentry;
            break;
         }
         else {
            if (primariesInSpot_ >= kEventsPerRun) { // same spot, new Run
               primariesInSpot_ = 0;
               lastJentry_ = jentry;
               break;
            }
         }
      }

      if (lastEventID < 0) {
         lastZ = posZ;
         lastLayer = layer;
         lastParentID = parentID;
         lastTrackID = trackID;
         lastPDG = PDGEncoding;
         strcpy(lastProcessName, processName);
         sumX = 0;
         sumY = 0;
         sum_edep = 0;
         n = 0;
      }

      layer = level1ID + baseID - 1;
//      if (kPhantom) layer++;

      if (kFilterNuclearInteractions == true && parentID != 0 && PDGEncoding != 11) {
         isSecondary = true;
         continue;
      }

      if (lastEventID == eventID && trackID > 1) {// && PDGEncoding != 100002004) {
         isSecondary = true;
//         cout << "Particle " << PDGEncoding << " generated from " << lastProcessName[0] << " at " << lastZ << endl;
      }

      if (TString(lastProcessName) == "alphaInelastic" || lastParentID > 0) {
         isSecondary = true;
      }


      if (lastEventID != eventID || lastLayer != layer) { // new layer -- store summed information from LAST layer now 
         x = sumX/n / dx + nx/2;
         y = sumY/n / dy + ny/2;

         if (lastLayer < nLayers) {
            if (hits)      hits->appendPoint(x, y, lastLayer, sum_edep/14, lastEventID, isSecondary, lastPDG);
            if (clusters)  clusters->appendClusterEdep(x,y, lastLayer, sum_edep/14, lastEventID, isSecondary, lastPDG);
         }

         sum_edep = 0;
         sumY = 0;
         sumX = 0;
         n = 0;
      }


      if (!kSpotScanning) {
         if (eventID < eventIdFrom) {
            continue;
         }
         
         else if (eventID >= eventIdTo || jentry == nentries-1) { // LAST ENTRY, EXIT AND APPEND THIS
            lastJentry_ = jentry;
            break;
         }
      }

      if (eventID != lastEventID) { // new track
         primariesInSpot_++;
         isSecondary = false;
      }

      sum_edep += edep*1000;
      sumX += posX;
      sumY += posY;
      lastZ = posZ;
      lastEventID = eventID;
      lastParentID = parentID;
      lastTrackID = trackID;
      lastLayer = layer;
      lastPDG = PDGEncoding;
      strcpy(lastProcessName, processName);
      n++;
   }
   
   if ((lastJentry_ == nentries-1) && (lastEventID != eventID || lastLayer != layer || lastParentID != parentID)) { // append last FIX THIS 
      fChain->GetEntry(lastJentry_);
      x = sumX/n / dx + nx/2;
      y = sumY/n / dy + ny/2;

      if (lastLayer < nLayers && !std::isnan(x+y)) {
         if (hits)      hits->appendPoint(x, y, lastLayer, sum_edep/14, lastEventID, isSecondary, lastPDG);
         if (clusters)  clusters->appendClusterEdep(x, y, lastLayer, sum_edep/14, lastEventID, isSecondary, lastPDG);
      }

      sum_edep = 0;
      sumY = 0;
      sumX = 0;
      n = 0;
   }
}

Int_t DataInterface::getMCFrame(Int_t runNo, CalorimeterFrame *cf) {
   // Retrieve kEventsPerRun events and put them into CalorimeterFrame*
   
   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

   if (runNo == 0) lastJentry_ = 0;

   Float_t offsetX = (nx/2+2) * dx;
   Float_t offsetY = (ny/2) * dy;
   Float_t x,y;
   Int_t calorimeterLayer = -1;

   if (fChain == 0) return -1;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   showDebug("Nentries = " << fChain->GetEntries() << endl);

   for (Long64_t jentry=lastJentry_; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry<0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (eventID < eventIdFrom) continue;
      else if (eventID >= eventIdTo) {
         lastJentry_ = jentry;
         break;
      }

      calorimeterLayer = level1ID + baseID - 1; 

      x = (posX + offsetX) * nx/2 / (offsetX);
      y = (posY + offsetY) * ny/2 / (offsetY);

      if (calorimeterLayer < nLayers) {
         cf->fillAt(calorimeterLayer, x, y, edep*1000/14);
      }
   }

   if (kEventsPerRun > 1) return -1;
   return eventID;
}

Int_t DataInterface::getMCFrame(Int_t runNo, Layer *l) {
   // Retrieve kEventsPerRun events and put them into Layer*
   
   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

   if (runNo == 0) lastJentry_ = 0;

   Float_t offsetX = (nx/2+2) * dx;
   Float_t offsetY = (ny/2) * dy;
   Float_t x,y;
   Int_t calorimeterLayer = -1;

   if (fChain == 0) return -1;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   showDebug("Nentries = " << fChain->GetEntries() << endl);

   for (Long64_t jentry=lastJentry_; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry<0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      if (eventID < eventIdFrom) continue;
      else if (eventID >= eventIdTo) {
         lastJentry_ = jentry;
         break;
      }

      calorimeterLayer = level1ID + baseID - 1; 

      x = (posX + offsetX) * nx/2 / (offsetX);
      y = (posY + offsetY) * ny/2 / (offsetY);

      if (l->getLayer() == calorimeterLayer) {
         l->Fill(x, y, edep*1000/14);
      }
   }

   if (kEventsPerRun > 1) return -1;
   return eventID;
}
#endif
