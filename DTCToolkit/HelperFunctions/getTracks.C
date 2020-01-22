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
#include <TLeaf.h>
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

Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy, Float_t spotx, Float_t spoty, Clusters * clusters) {
   Tracks * tracks = nullptr;
   
   if (recreate) {
      tracks = getTracksFromClusters(Runs, dataType, kCalorimeter, energy, spotx, spoty, clusters);
      saveTracks(tracks, clusters, dataType, energy);
   }

   else {
      tracks = loadTracks(Runs, dataType, energy, clusters);
      if (!tracks) {
         cout << "!tracks, creating new file\n";
         tracks = getTracksFromClusters(Runs, dataType, kCalorimeter, energy, spotx, spoty, clusters);
         saveTracks(tracks, clusters, dataType, energy);
      }
   }
   return tracks;
}

Tracks * getTracksFromClusters(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy, Float_t spotx, Float_t spoty, Clusters * saveClusters) {
   DataInterface   * di = new DataInterface();
   Int_t             nClusters = kEventsPerRun * 5 * nLayers;
   Int_t             nTracks = kEventsPerRun * 2;
   Bool_t            breakSignal = false;
   Clusters        * clusters = nullptr;
   Tracks          * tracks = nullptr;
   Tracks          * allTracks = new Tracks(nTracks * Runs);
   TRandom3        * gRandom = new TRandom3(0);

   allTracks->SetOwner(kTRUE);

   TStopwatch t1, t2, t3, t4, t5, t6;
   t1.Reset(); // get MC clusters
   t2.Reset(); // Stop
   t3.Reset();
   t4.Reset();
   t5.Reset();
   t6.Reset();

   for (Int_t i=0; i<Runs; i++) {
      clusters = new Clusters(nClusters);
      showDebug("Start getMCClusters\n");
      t1.Start(false);
      if (kDoDiffusion) {
         Hits * hits = new Hits();
         Hits * diffusedHits = nullptr;
         di->getMCClustersThreshold(i, nullptr, hits, spotx, spoty);
         hits->sortHits(); 
         diffusedHits = diffuseHits(gRandom, hits);
         diffusedHits->sortHits();
         diffusedHits->makeLayerIndex();
         diffusedHits->findClustersFromHits(clusters);
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

//      printf("Found %d clusters...\n", clusters->GetEntriesFast());
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
         tracks = clusters->findCalorimeterTracksWithMCTruth();
      }
      t2.Stop();
    
      showDebug("propagateSecondaryStatus..."); 
//      tracks->propagateSecondaryStatus();
      showDebug("ok!\n");

      if (tracks->GetEntriesFast() == 0) breakSignal = true; // to stop running

      showDebug("Found all the tracks\n");

      // Track improvements
      Int_t nTracksBefore = 0, nTracksAfter = 0;
      Int_t nIsInelastic = 0, nIsNotInelastic = 0;
      showDebug("removeNANs...");
      tracks->removeNANs();
      showDebug("ok!\nsortTracks...");
      tracks->sortTracks(); // reverse order from retrograde reconstruction
      tracks->removeEmptyTracks();
      
      t3.Start(false);
      showDebug("ok!\nRemoving tracks leaving detector...");
      tracks->removeTracksLeavingDetector();
      showDebug("ok\ncompress tracks and clusters...");      
      t3.Stop();
      t4.Start(false);
      tracks->Compress();
      tracks->CompressClusters();
      showDebug("ok\n");
      
      // Flag tracks as too short
      Track * thisTrack;
      tracks->appendClustersWithoutTrack(clusters->getClustersWithoutTrack());
      for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
         thisTrack = tracks->At(i);
         if (tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), i)) {
            thisTrack->setIncomplete(true);
         }
      }

      showDebug("append tracks to alltracks...");
      for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
         if (!tracks->At(j)) continue;
         allTracks->appendTrack(tracks->At(j));
      }
      showDebug("ok\n");

      showDebug("appendClustersWithoutTrack...");
      allTracks->appendClustersWithoutTrack(clusters->getClustersWithoutTrack());
      t4.Stop();
      showDebug("ok\n");

      delete clusters;
      delete tracks;

      if (breakSignal) break;
   }

   printf("Timing: Cluster retrieval: %.3f s. Tracking: %.3f s. Track improvements: %.3f s. Logistics: %.3f s. Sort TCA by layer: %.3f s\n", t1.CpuTime(), t2.CpuTime(), t3.CpuTime(), t4.CpuTime(), t6.CpuTime());

   delete di;
   delete gRandom;

   return allTracks;
}

Tracks * getTracksFromClustersMT(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy, Float_t spotx, Float_t spoty, Clusters * saveClusters) {
   DataInterface   * di = new DataInterface();
   Int_t             nClusters = kEventsPerRun * 5 * nLayers;
   Int_t             nTracks = kEventsPerRun * 10;
   Bool_t            breakSignal = false;
   Clusters        * clusters = nullptr;
   TRandom3        * gRandom = new TRandom3(0);
   Int_t             nThreads = 8; // Move to Constants.h
   Int_t             runID = 0;
   
   TH2I            * hOriginal = new TH2I("hOriginal", "Original", 201, 4399.5, 4600.5, 201, 2149.5, 2350.5);
   TH2I            * hDiffused = new TH2I("hDiffused", "Diffused", 201, 4399.5, 4600.5, 201, 2149.5, 2350.5);
   TH2I            * hReconstructed = new TH2I("hReconstructed", "Reconstructed",201, 4399.5, 4600.5, 201, 2149.5, 2350.5);

   TStopwatch t1, t2, t3, t4, t5;
   t1.Reset(); 
   t2.Reset(); 
   t3.Reset();
   t4.Reset();
   t5.Reset();
  
   t1.Start(false); 
   TFile clusterFile("OutputFiles/clusterFile2.root", "recreate");
   TTree * clusterTree = new TTree("Clusters", "clusters");

   clusterTree->Branch("clusters", &clusters, 256000, 0);
   clusterTree->Branch("runID", &runID, 256000, 0);

   for (runID=0; runID<Runs; runID++) {
      clusters = new Clusters(nClusters);
      if (kDoDiffusion) {
         Hits * hits = new Hits();
         Hits * diffusedHits = nullptr;
         t3.Start(false);
         di->getMCClustersThreshold(runID, nullptr, hits, spotx, spoty);
         hits->removeHaloAtRadius(15);
         hits->sortHits();
         t3.Stop();
         t4.Start(false);
         diffusedHits = diffuseHitsMT(hits);
         diffusedHits->sortHits();
         t4.Stop();
         t5.Start(false);
         diffusedHits->makeLayerIndex();
         diffusedHits->findClustersFromHits(clusters);
         t5.Stop();
         if (runID%10 == 0) { 
            printf("Run %d of %d: Found %d clusters from %d diffused Hits\n", runID, Runs, clusters->GetEntriesFast(), hits->GetEntriesFast());
         }

         Int_t bin;
         for (Int_t i=0; i<hits->GetEntriesFast(); i++) {
            bin = hOriginal->FindBin(hits->At(i)->getX(), hits->At(i)->getY());

            if (hits->getLayer(i) == 20) hOriginal->SetBinContent(bin, 1);
         }
         for (Int_t i=0; i<diffusedHits->GetEntriesFast(); i++) {
            bin = hDiffused->FindBin(diffusedHits->At(i)->getX(), diffusedHits->At(i)->getY());
            if (diffusedHits->getLayer(i) == 20) hDiffused->SetBinContent(bin, 1);
         }
         for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
            bin = hReconstructed->FindBin(clusters->At(i)->getX(), clusters->At(i)->getY());
            if (clusters->getLayer(i) == 20) hReconstructed->SetBinContent(bin, 1);
         }

         delete hits;
         delete diffusedHits;
      }

      else {
         di->getMCClustersThreshold(runID, clusters, nullptr, spotx, spoty);
         clusters->removeHaloAtRadius(15);
      }
      
      if (saveClusters) {
         for (Int_t j=0; j<=clusters->GetEntriesFast(); j++) {
            if (!clusters->At(j)) continue;
            saveClusters->appendCluster(clusters->At(j));
         }
      }
      
      if (!clusters->GetEntriesFast()) continue;
      clusters->sortClusters();

//      printf("Filling %d clusters\n", clusters->GetEntriesFast());
      clusterTree->Fill();
      delete clusters;
   }

   TCanvas *cRecon = new TCanvas("cRecon", "cRecon", 1500, 800);
   cRecon->Divide(3,1);

   cRecon->cd(1);
   hOriginal->Draw("colz");
   cRecon->cd(2);
   hDiffused->Draw("colz");
   cRecon->cd(3);
   hReconstructed->Draw("colz");

   clusterFile.ReOpen("update"); // To bump the file to the top, otherwise the ...closed... tfile in diffuseHitsMT is used

   clusterTree->Write();
   clusterFile.Write();
   t1.Stop();

   t2.Start(false);
   ROOT::EnableImplicitMT(nThreads);
   ROOT::TThreadedObject<Tracks> tracksTTO(nTracks);
   ROOT::TTreeProcessorMT tp(*clusterTree);

   tp.SetMaxTasksPerFilePerWorker(2);
   
   auto trackingFunction = [&] ( TTreeReader &reader ) {
      TTreeReaderValue<int> runIDRV(reader, "runID");
      TTreeReaderValue<Clusters> clustersRV(reader, "clusters");
      
      auto tracksMT = tracksTTO.Get();

      while (reader.Next()) {
         auto clusterBatch = *clustersRV;
         auto runID = *runIDRV;

         if (runID%10 == 0) printf("Running MT tracking with runID = %d\n", runID);
         Tracks * localTracks = clusterBatch.findTracksWithRecursiveWeighting();

         for (Int_t trackID=0; trackID<localTracks->GetEntriesFast(); trackID++) {
            tracksMT->appendTrack(localTracks->At(trackID));
         }
         tracksMT->appendClustersWithoutTrack(clusterBatch.getClustersWithoutTrack());
         delete localTracks;
      }
   };

   tp.Process(trackingFunction);
   ROOT::DisableImplicitMT();
   
   auto allTracksShared = tracksTTO.Merge();
   auto allTracks = new Tracks();
   for (Int_t i=0; i<allTracksShared->GetEntriesFast(); i++) {
      allTracks->appendTrack(allTracksShared->At(i));
   }
   allTracks->appendClustersWithoutTrack(allTracksShared->getClustersWithoutTrack());
      
   // Flag tracks as too short
   Track * thisTrack;
   for (Int_t i=0; i<allTracks->GetEntriesFast(); i++) {
      thisTrack = allTracks->At(i);
      if (allTracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), i)) {
         thisTrack->setIncomplete(true);
      }
   }

   allTracks->SetOwner(kTRUE);
   clusterFile.Close();

   t2.Stop();

   // DO THIS AFTER MT
   //
   // Track improvements
   allTracks->removeNANs();
   allTracks->sortTracks(); // reverse order from retrograde reconstruction
   allTracks->removeTracksLeavingDetector();
   
   allTracks->Compress();
   allTracks->CompressClusters();
   allTracks->removeEmptyTracks();
   

   printf("Timing: Cluster retrieval: %.3f s. Tracking: %.3f s. hits log %.3f s, diffusion %.3f s, findClusterFromHits %.3f s\n", t1.CpuTime(), t2.CpuTime(), t3.CpuTime(), t4.CpuTime(), t5.CpuTime());

   delete di;

   printf("allTracks has %d tracks\n", allTracks->GetEntriesFast());

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

   showDebug("Diffusing hits (ALPIDE-Heidelberg). Number of hits = " << nHits << endl);

   Int_t binPosLUT[12] = {};
   for (Int_t i=0; i<12; i++) { binPosLUT[i] = pow(2, i); }
   
   int circleX[70] = {0,1,0,-1,0,1,-1,-1,1,0,-2,0,2,1,-2,-1,2,-1,-2,1,2,-2,-2,2,2,0,-3,0,3,-1,-3,1,3,1,-3,-1,3,0,-4,0,4,2,-3,-2,3,-8,-2,-3,2,4,-1,-4,1,4,1,7,-1,3,3,-3,-3,4,2,-4,-2,4,-2,2,5,0};

   int circleY[70] = {0,0,-1,0,1,-1,-1,1,1,-2,0,2,0,-2,-1,2,1,-2,1,2,-1,-2,2,2,-2,-3,0,3,0,-3,1,3,-1,-3,-1,3,1,-4,0,4,0,-3,-2,3,2,15,-3,2,3,-1,-4,1,4,1,-4,-18,4,3,-3,-3,3,2,-4,-2,4,-2,-4,4,0,5};

   for (Int_t h=0; h<nHits; h++) {
      x = hits->getX(h);
      y = hits->getY(h);
      layer = hits->getLayer(h);
      edep = hits->getEdep(h);
      eventID = hits->getEventID(h);
      isSecondary = hits->isSecondary(h);
      cs = getCSFromEdep(edep);
      PDG = hits->getPDG(h);
       
      randomClusterIdx = gRandom->Integer(CDB_sortIndex[cs+1] - CDB_sortIndex[cs]) + CDB_sortIndex[cs];
      CDB_treeCluster->GetEntry(randomClusterIdx);

      if (cs < 27 && kUseExperimentalClusterPainting) {
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

Hits * diffuseHitsMT(Hits * hits) {
   Int_t          nHits = hits->GetEntriesFast();
   Hits         * hitsOut = new Hits();
   Int_t          x, y, outX, outY, layer, cs, eventID, idx_x, binPos, PDG;
   Float_t        edep;
   Bool_t         isSecondary;
   Int_t          randomClusterIdx;
   Int_t          nBefore, nAfter;
   Int_t          nThreads = 8;

   showDebug("Diffusing hits (ALPIDE-Heidelberg). Number of hits = " << nHits << endl);

   // Make TTree out of the Hits
   TFile hitsFile("OutputFiles/hitsFile2.root", "recreate");
   TTree * hitsTree = new TTree("Hits", "hits");
   hitsTree->Branch("x", &x, 256000, 1);
   hitsTree->Branch("y", &y, 256000, 1);
   hitsTree->Branch("layer", &layer, 256000, 1);
   hitsTree->Branch("edep", &edep, 256000, 1);
   hitsTree->Branch("cs", &cs, 256000, 1);
   hitsTree->Branch("eventID", &eventID, 256000, 1);
   hitsTree->Branch("isSecondary", &isSecondary, 256000, 1);
   hitsTree->Branch("PDG", &PDG, 256000, 1);

   for (Int_t h=0; h<nHits; h++) {
      x = hits->getX(h);
      y = hits->getY(h);
      layer = hits->getLayer(h);
      eventID = hits->getEventID(h);
      isSecondary = hits->isSecondary(h);
      edep = hits->getEdep(h);
      cs = getCSFromEdep(edep);
      PDG = hits->getPDG(h);
      if (cs<2) continue;
//      if (cs>=27) cs=26;

      hitsTree->Fill();
   }
   hitsTree->Write();
   hitsFile.Write();

   ROOT::EnableImplicitMT(nThreads);
   ROOT::TThreadedObject<Hits> hitsTTO;
   ROOT::TTreeProcessorMT tp(*hitsTree);

   // manually drawn circle; starting from center growing spirally ca symmetrically
   int circleX[70] = {0,1,0,-1,0,1,-1,-1,1,0,-2,0,2,1,-2,-1,2,-1,-2,1,2,-2,-2,2,2,0,-3,0,3,-1,-3,1,3,1,-3,-1,3,0,-4,0,4,2,-3,-2,3,-8,-2,-3,2,4,-1,-4,1,4,1,7,-1,3,3,-3,-3,4,2,-4,-2,4,-2,2,5,0};

   int circleY[70] = {0,0,-1,0,1,-1,-1,1,1,-2,0,2,0,-2,-1,2,1,-2,1,2,-1,-2,2,2,-2,-3,0,3,0,-3,1,3,-1,-3,-1,3,1,-4,0,4,0,-3,-2,3,2,15,-3,2,3,-1,-4,1,4,1,-4,-18,4,3,-3,-3,3,2,-4,-2,4,-2,-4,4,0,5};
   
   auto diffuseFunction = [&] ( TTreeReader &reader ) {
      TTreeReaderValue<int> xRV(reader, "x");
      TTreeReaderValue<int> yRV(reader, "y");
      TTreeReaderValue<int> layerRV(reader, "layer");
      TTreeReaderValue<int> csRV(reader, "cs");
      TTreeReaderValue<int> eventIDRV(reader, "eventID");
      TTreeReaderValue<bool> isSecondaryRV(reader, "isSecondary");
      TTreeReaderValue<float> edepRV(reader, "edep");
      TTreeReaderValue<int> pdgRV(reader, "PDG");
      
      auto hitsMT = hitsTTO.Get();
      TRandom3 *gRandomHere = new TRandom3(0);
      
      Int_t binPosLUT[12] = {};
      for (Int_t i=0; i<12; i++) { binPosLUT[i] = pow(2, i); }

      while (reader.Next()) {
         auto thisX = *xRV;
         auto thisY = *yRV;
         auto thisLayer = *layerRV;
         auto thisCS = *csRV;
         auto thisEventID = *eventIDRV;
         auto thisIsSecondary = *isSecondaryRV;
         auto thisEdep = *edepRV;
         auto thisPDG = *pdgRV;

         if (thisCS < 27 && kUseExperimentalClusterPainting) {
            randomClusterIdx = gRandomHere->Integer(CDB_sortIndex[thisCS+1] - CDB_sortIndex[thisCS]) + CDB_sortIndex[thisCS];
            CDB_treeCluster->GetEntry(randomClusterIdx);

            idx_x = 0;
            for (Int_t i=0; i<10; i++) {
               for (Int_t binPosPow = 0; binPosPow < 10; binPosPow++) {
                  binPos = binPosLUT[binPosPow];
                  if (binPos & CDB_hit_array[i]) {
                     outX = thisX + (idx_x - CDB_x_mean) + 0.5;
                     outY = thisY + (binPosPow - CDB_y_mean) + 0.5;
                     hitsMT->appendPoint(outX, outY, thisLayer, thisEdep/CDB_clusterSize, thisEventID, thisIsSecondary, thisPDG);
                  }
               }
            idx_x++;
            }
         }
         else { // NO DATA EXISTS, MAKE CIRCLE
            int circleSize = min(thisCS, 70);
            for (Int_t i=0; i<circleSize; i++) {
               outX = thisX + circleX[i] + 0.5;
               outY = thisY + circleY[i] + 0.5;
               hitsMT->appendPoint(outX, outY, thisLayer, thisEdep/circleSize, thisEventID, thisIsSecondary, thisPDG);
            }
         }
      }
   };

   tp.Process(diffuseFunction);
   ROOT::DisableImplicitMT();
   
   auto hitsShared = hitsTTO.Merge();
   for (Int_t i=0; i<hitsShared->GetEntriesFast(); i++) {
      hitsOut->appendHit(hitsShared->At(i));
   }

   hitsFile.Close(); 

   return hitsOut;
   
}

#endif
