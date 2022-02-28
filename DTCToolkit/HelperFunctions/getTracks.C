#ifndef getTracks_cxx
#define getTracks_cxx

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <list>

#include <TRandom3.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TView.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeIndex.h>
#include <TFile.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "GlobalConstants/RangeAndEnergyCalculations.h"
#include "GlobalConstants/Misalign.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/DataInterface/DataInterface.h"
#include "HelperFunctions/Tools.h"
#include "HelperFunctions/getTracks.h"

#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace std;
using namespace DTC;

void makeTracks(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   tracks->extrapolateToLayer0();
}

void saveTracks(Tracks *tracks, Clusters * clusters, Int_t dataType, Float_t energy) {
   TString sDataType = (dataType == 0) ? "_MC_" : "_data_";

   tracks->CompressCWT();

   Float_t readoutAbsorber = (kAbsorberThickness == roundf(kAbsorberThickness)) ? kAbsorberThickness : kAbsorberThickness*10;

   TString sEnergy = Form("_%.2fMeV", energy);
   TString fileName = "Data/Tracks/tracks";
   TString sAbsThick = Form("_%.1fmm", readoutAbsorber);
   TString sDegThick = Form("_%.0fmm", run_degraderThickness);
   TString sMaterial = getMaterialChar();
   fileName.Append(sDataType);
   fileName.Append(sMaterial);
   fileName.Append(sAbsThick);
   fileName.Append(sEnergy);
   fileName.Append(sDegThick);
   fileName.Append(".root");
   
   TFile f(fileName, "recreate");
   f.SetCompressionLevel(1);
   TTree T("T", "tracks");
   T.Branch("tracks", &tracks, 256000, 1);
   T.Branch("clusters", &clusters, 256000, 1);
   T.Fill();
   T.Write();
   f.Close();
}

Tracks * loadTracks(Int_t Runs, Int_t dataType, Float_t energy, Clusters * clusters) {
   Float_t readoutAbsorber = (kAbsorberThickness == roundf(kAbsorberThickness)) ? kAbsorberThickness : kAbsorberThickness*10;

   TString sDataType = (dataType == 0) ? "_MC_" : "_data_";
   TString sAbsThick = Form("_%.1fmm", readoutAbsorber);
   TString sEnergy = Form("_%.2fMeV", energy);
   TString sDegThick = Form("_%.0fmm", run_degraderThickness);
   TString fileName = "Data/Tracks/tracks";
   TString sMaterial = getMaterialChar();
   fileName.Append(sDataType);
   fileName.Append(sMaterial);
   fileName.Append(sAbsThick);
   fileName.Append(sEnergy);
   fileName.Append(sDegThick);
   fileName.Append(".root");
   
   TFile *f = new TFile(fileName);
   if (!f) return 0;
   TTree *T = (TTree*) f->Get("T");
   Tracks * tracks = new Tracks();

   T->GetBranch("tracks")->SetAutoDelete(kFALSE);
   T->SetBranchAddress("tracks",&tracks);
   T->SetBranchAddress("clusters", &clusters);

   T->GetEntry(0);
   
   cout << "There are " << tracks->GetEntriesFast() << " tracks in " << fileName << ".\n";
   
   return tracks;
}

Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy, Float_t spotx, Float_t spoty, Clusters * clusters, Int_t startRun) {
   Tracks * tracks = nullptr;
   
   if (recreate) {
      tracks = getTracksFromClusters(Runs, dataType, kCalorimeter, energy, spotx, spoty, clusters, startRun);
      if (tracks) saveTracks(tracks, clusters, dataType, energy);
   }

   else {
      tracks = loadTracks(Runs, dataType, energy, clusters);
      if (!tracks) {
         cout << "!tracks, creating new file\n";
         tracks = getTracksFromClusters(Runs, dataType, kCalorimeter, energy, spotx, spoty, clusters, startRun);
         saveTracks(tracks, clusters, dataType, energy);
      }
   }
   return tracks;
}

Tracks * getTracksFromClusters(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy, Float_t spotx, Float_t spoty, Clusters * saveClusters, Int_t startRun) {
   DataInterface   * di = new DataInterface();
   Int_t             nClusters = kEventsPerRun * 10 * nLayers;
   Int_t             nTracks = kEventsPerRun * 2;
   Bool_t            breakSignal = false;
   Clusters        * clusters = nullptr;
   Clusters        * clustersCloned = nullptr;
   Tracks          * tracks = nullptr;
   Tracks          * allTracks = new Tracks(nTracks * Runs);
   TRandom3        * gRandom = new TRandom3(0);

   allTracks->SetOwner(kTRUE);

   TStopwatch t1, t2, t3, t4, t6, t7, t8, t9, t10;
   t1.Reset(); // get MC clusters
   t2.Reset(); // Stop
   t3.Reset();
   t4.Reset();
   t6.Reset();
   t7.Reset();
   t8.Reset();
   t9.Reset();
   t10.Reset();

   t4.Start(false);
   Cluster * thisCluster = nullptr;
   for (Int_t i=startRun; i<Runs+startRun; i++) {
      clusters = new Clusters(nClusters);
      showDebug("Start getMCClusters\n");
      t1.Start(false);
      if (kDoDiffusion) {
         Hits * hits = new Hits();
         Hits * diffusedHits = nullptr;
         t8.Start(false); di->getMCClustersThreshold(i, nullptr, hits, spotx, spoty); t8.Stop();
         hits->sortHits();

         // DETECTOR RESPONSE THRESHOLD
         for (Int_t j=0; j<=hits->GetEntriesFast(); j++) {
            if (hits->At(j)) {
               if (hits->At(j)->getEdep() < 0.1) {
                  hits->removeHitAt(j);
               }
            }
         }

         hits->Compress();
         t6.Start(false); diffusedHits = diffuseHits(gRandom, hits); t6.Stop();
//         diffusedHits->sortHits(); /// THIS IS THE CULPRIT TIMEWISE!!! DO WE ACTUALLY NEED IT? SHOULD BE SORTED ESP BY LAYER
         diffusedHits->makeLayerIndex();
         t7.Start(false); diffusedHits->findClustersFromHits(clusters); t7.Stop();

         // EDEPish CUT
         if (kHelium) {
            clusters->removeSmallClusters(5);
         }
         else {
            clusters->removeSmallClusters(1); // VALIDATE !!!!
         }

//         Int_t nRem = clusters->removeClustersInGap(1, 0);

         delete hits;
         delete diffusedHits;
      }
      else {
         di->getMCClustersThreshold(i, clusters, nullptr, spotx, spoty);
         showDebug("Found " << clusters->GetEntries() << " clusters...removing halo...");
         showDebug("ok\n");
//       Int_t nRem = clusters->removeClustersInGap(1, 0); 
      }
      if (saveClusters) {
         for (Int_t j=0; j<=clusters->GetEntriesFast(); j++) {
            if (!clusters->At(j)) continue;
            saveClusters->appendCluster(clusters->At(j));
         }
      }

      if (!clusters->GetEntriesFast()) continue;

      t1.Stop();
   
      showDebug("Finding calorimeter tracks\n");
      t2.Start(false);
      if (kDoTracking) {
         showDebug("sortClusters...");
         clusters->sortClusters();
         showDebug("ok!\n Start tracking (kWeightedRecursive)\n");
         tracks = clusters->findTracksWithRecursiveWeighting();
         if (i%10 == 0) printf("Found %d tracks from %d Clusters in run %.2f-%d.\n", tracks->GetEntriesFast(), clusters->GetEntriesFast(), run_energy, i);
      }
      else {
         showDebug("sortClusters...");
         clusters->sortClusters();
         tracks = clusters->findCalorimeterTracksWithMCTruth();
      }
      t2.Stop();
    
      showDebug("propagateSecondaryStatus..."); 
//      tracks->propagateSecondaryStatus();
      showDebug("ok!\n");

      if (tracks->GetEntriesFast() == 0) breakSignal = true; // to stop running

      showDebug("Found all the tracks ( " << tracks->GetEntriesFast() << ")\n");

      // Track improvements
      Int_t nTracksBefore = 0, nTracksAfter = 0;
      Int_t nIsInelastic = 0, nIsNotInelastic = 0;
     
      tracks->removeNANs();
      tracks->sortTracks(); // reverse order from retrograde reconstruction
      tracks->appendClustersWithoutTrack(clusters->getClustersWithoutTrack());
      
      t3.Start(false);
      tracks->removeEmptyTracks();
//      tracks->removeShortTracks(4); // rem [0, 4]
//      tracks->removeTracksLeavingDetector();
      tracks->fillOutIncompleteTracks(0.4);
      t3.Stop();
      
      // tracks->removeNuclearInteractions();
//      tracks->removeTracksWithMinWEPL(100);

      tracks->Compress();
      tracks->CompressClusters();
      showDebug("ok\n");

      // Marking incomplete tracks
      Track * thisTrack = nullptr;
      for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
         thisTrack = tracks->At(j);
         if (tracks->isTrackIncomplete(thisTrack)) {
            thisTrack->setIncomplete(true);
         }
      }
     
      showDebug("append tracks to alltracks...");
      showDebug(" (" << tracks->GetEntriesFast() << ") tracks\n");
      for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
         if (!tracks->At(j)) continue;
         allTracks->appendTrack(tracks->At(j));
      }
      showDebug("ok\n");

      showDebug("appendClustersWithoutTrack...");
//      if (kSaveCWT) allTracks->appendClustersWithoutTrack(tracks->getClustersWithoutTrack());
      if (kSaveCWT) allTracks->appendClustersWithoutTrack(clusters->getClustersWithoutTrack());
      showDebug("ok\n");

      delete clusters;
      delete tracks;

      if (breakSignal) break;
   }
   t4.Stop();

   printf("Timing: Total: %.3f s. Cluster retrieval: %.3f s. Tracking: %.3f s. Track improvements: %.3f s. Diffuse hits: %.3f. Find clusters %.3f Readout %.3fs,\n", t4.CpuTime(), t1.CpuTime(), t2.CpuTime(), t3.CpuTime(), t6.CpuTime(), t7.CpuTime(), t8.CpuTime());

   delete di;
   delete gRandom;

   return allTracks;
}

Hits * diffuseHits(TRandom3 *gRandom, Hits * hits) {
   Int_t          nHits = hits->GetEntriesFast();
   Hits         * hitsOut = new Hits();
   Int_t          x, y, outX, outY, layer, cs, eventID, idx_x, binPos, PDG;
   Float_t        edep;
   Bool_t         isSecondary;
   Int_t          randomClusterIdx;
   Int_t          nBefore, nAfter;
   TStopwatch t1;
   t1.Reset();

   showDebug("Diffusing hits (ALPIDE-Heidelberg). Number of hits = " << nHits << endl);

   Int_t binPosLUT[12] = {};
   for (Int_t i=0; i<12; i++) { binPosLUT[i] = pow(2, i); }
   
   Short_t circleX[70] = {0,1,0,-1,0,1,-1,-1,1,0,-2,0,2,1,-2,-1,2,-1,-2,1,2,-2,-2,2,2,0,-3,0,3,-1,-3,1,3,1,-3,-1,3,0,-4,0,4,2,-3,-2,3,-8,-2,-3,2,4,-1,-4,1,4,1,7,-1,3,3,-3,-3,4,2,-4,-2,4,-2,2,5,0};
   Short_t circleY[70] = {0,0,-1,0,1,-1,-1,1,1,-2,0,2,0,-2,-1,2,1,-2,1,2,-1,-2,2,2,-2,-3,0,3,0,-3,1,3,-1,-3,-1,3,1,-4,0,4,0,-3,-2,3,2,15,-3,2,3,-1,-4,1,4,1,-4,-18,4,3,-3,-3,3,2,-4,-2,4,-2,-4,4,0,5};

   for (Int_t h=0; h<nHits; h++) {
      x = hits->getX(h);
      y = hits->getY(h);
      layer = hits->getLayer(h);
      edep = hits->getEdep(h);
      eventID = hits->getEventID(h);
      isSecondary = hits->isSecondary(h);
      cs = getCSFromEdep(edep);
      PDG = hits->getPDG(h);

      if (cs < 27 && kUseExperimentalClusterPainting) {
         randomClusterIdx = gRandom->Integer(CDB_sortIndex[cs+1] - CDB_sortIndex[cs]) + CDB_sortIndex[cs];
         CDB_treeCluster->GetEntry(randomClusterIdx);
         idx_x = 0;
         for (Int_t i=0; i<10; i++) {
            for (Int_t binPosPow = 0; binPosPow < 10; binPosPow++) {
               binPos = binPosLUT[binPosPow];
               if (binPos & CDB_hit_array[i]) {
                  outX = x + (idx_x - CDB_x_mean) + 0.5;
                  outY = y + (binPosPow - CDB_y_mean) + 0.5;
                  hitsOut->appendPoint(outX, outY, layer, edep/CDB_clusterSize, eventID, isSecondary, PDG);
               }
            }
         idx_x++;
         }
      }
      else {
         int circleSize = min(cs, 70);
         for (Int_t i=0; i<circleSize; i++) {
            outX = x + circleX[i] + 0.5;
            outY = y + circleY[i] + 0.5;
            hitsOut->appendPoint(outX, outY, layer, edep/circleSize, eventID, isSecondary, PDG);
         }
      }
   }
   
   return hitsOut;
}

#endif
