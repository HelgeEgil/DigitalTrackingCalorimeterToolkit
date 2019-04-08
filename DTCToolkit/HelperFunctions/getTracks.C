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

using namespace std;
using namespace DTC;

void makeTracks(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   tracks->extrapolateToLayer0();
}

void saveTracks(Tracks *tracks, Int_t dataType, Float_t energy) {
   TString sDataType = (dataType == 0) ? "_MC_" : "_data_";

   tracks->CompressCWT();

   Float_t readoutAbsorber = (kAbsorberThickness == roundf(kAbsorberThickness)) ? kAbsorberThickness : kAbsorberThickness*10;

   TString sEnergy = Form("_%.2fMeV", energy);
   TString fileName = "Data/Tracks/tracks";
   TString sAbsThick = Form("_%.0fmm", readoutAbsorber);
   TString sMaterial = getMaterialChar();
   fileName.Append(sDataType);
   fileName.Append(sMaterial);
   fileName.Append(sAbsThick);
   fileName.Append(sEnergy);
   fileName.Append(".root");
   
   TFile f(fileName, "recreate");
   f.SetCompressionLevel(1);
   TTree T("T", "tracks");
   T.Branch("tracks", &tracks, 256000, 1);
   T.Fill();
   T.Write();
   f.Close();
}

Tracks * loadTracks(Int_t Runs, Int_t dataType, Float_t energy) {
   Float_t readoutAbsorber = (kAbsorberThickness == roundf(kAbsorberThickness)) ? kAbsorberThickness : kAbsorberThickness*10;

   TString sDataType = (dataType == 0) ? "_MC_" : "_data_";
   TString sAbsThick = Form("_%.0fmm", readoutAbsorber);
   TString sEnergy = Form("_%.2fMeV", energy);
   TString fileName = "Data/Tracks/tracks";
   TString sMaterial = getMaterialChar();
   fileName.Append(sDataType);
   fileName.Append(sMaterial);
   fileName.Append(sAbsThick);
   fileName.Append(sEnergy);
   fileName.Append(".root");
   
   TFile *f = new TFile(fileName);
   if (!f) return 0;
   TTree *T = (TTree*) f->Get("T");
   Tracks * tracks = new Tracks();

   T->GetBranch("tracks")->SetAutoDelete(kFALSE);
   T->SetBranchAddress("tracks",&tracks);


   T->GetEntry(0);
   
   cout << "There are " << tracks->GetEntriesFast() << " tracks in " << fileName << ".\n";
   
   return tracks;
}

Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy) {
   Tracks * tracks = nullptr;
   
   if (recreate) {
      tracks = getTracksFromClusters(Runs, dataType, kCalorimeter, energy);
   }

   else {
      tracks = loadTracks(Runs, dataType, energy);
      if (!tracks) {
         cout << "!tracks, creating new file\n";
         tracks = getTracksFromClusters(Runs, dataType, kCalorimeter, energy);
         saveTracks(tracks, dataType, energy);
      }
   }
   return tracks;
}

Clusters * getClusters(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy) {
   DataInterface   * di = new DataInterface();
   Int_t             nClusters = kEventsPerRun * 5 * nLayers;
   Int_t             nHits = kEventsPerRun * 50;
   Int_t             nTracks = kEventsPerRun * 2;
   Bool_t            breakSignal = false;
   CalorimeterFrame *cf = new CalorimeterFrame();
   Clusters        * clusters = nullptr;
   Clusters        * trackerClusters = new Clusters(nClusters);
   Clusters        * allClusters = new Clusters(nClusters * Runs);
   Hits            * hits = new Hits(nHits);
   Hits            * eventIDs = new Hits(kEventsPerRun * sizeOfEventID);
   Int_t             eventID = -1;
   Hits            * trackerHits = new Hits(nHits);
   TRandom3        * gRandom = new TRandom3(0);

   for (Int_t i=0; i<Runs; i++) {

      cout << "Finding clusters " << i*kEventsPerRun << "->" << (i+1)*kEventsPerRun << " of " << Runs * kEventsPerRun << endl;

      if (dataType == kMC) {
         eventID = di->getMCFrame(i, cf);
         di->getEventIDs(i, eventIDs);
         cf->diffuseFrame(gRandom);
         hits = cf->findHits(eventID);
         clusters = hits->findClustersFromHits(); // badly optimized
         clusters->removeSmallClusters(2);

         clusters->matchWithEventIDs(eventIDs);
         eventIDs->Clear("C");
      }
      
      else if (dataType == kData) {
         di->getDataFrame(i, cf, energy);
         hits = cf->findHits();
         clusters = hits->findClustersFromHits();
         clusters->removeSmallClusters(2);
         clusters->removeAllClustersAfterLayer(8); // bad data in layer 10 and 11
      }
      
      clusters->Compress();
      
      if (clusters->GetEntriesFast() == 0) breakSignal = kTRUE; // to stop running

      for (Int_t j=0; j<clusters->GetEntriesFast(); j++) {
         allClusters->appendCluster(clusters->At(j));
      }

      cf->Reset();
      cf->Clear("C");
      hits->Clear("C");
      trackerHits->Clear("C");
      delete clusters;
      trackerClusters->Clear("C");
      
      if (breakSignal) break;
   }


   delete cf;
   delete trackerClusters;
   delete hits;
   delete trackerHits;
   delete di;

   return allClusters;
}

Hits * diffuseHits(TRandom3 *gRandom, Hits * hits) {
   Int_t          nHits = hits->GetEntriesFast();
   Hits         * hitsOut = new Hits();
   Int_t          x, y, outX, outY, layer, cs, eventID, idx_x, binPos;
   Float_t        edep;
   Bool_t         isSecondary;
   Int_t          randomClusterIdx;
   Int_t          nBefore, nAfter;

   showDebug("Diffusing hits (ALPIDE-Heidelberg). Number of hits = " << nHits << endl);

   for (Int_t h=0; h<nHits; h++) {
      x = hits->getX(h);
      y = hits->getY(h);
      layer = hits->getLayer(h);
      edep = hits->getEdep(h);
      eventID = hits->getEventID(h);
      isSecondary = hits->isSecondary(h);
      cs = getCSFromEdep(edep);
    
      if (cs<2) cs=2;
      if (cs>=27) cs=26;

      randomClusterIdx = gRandom->Integer(CDB_sortIndex[cs+1] - CDB_sortIndex[cs]) + CDB_sortIndex[cs];
      CDB_treeCluster->GetEntry(CDB_treeCluster->LoadTree(CDB_index->GetIndex()[randomClusterIdx]));

      idx_x = 0;
      for (Int_t n : *CDB_hit_array) {
         for (Int_t binPosPow = 0; binPosPow < 10; binPosPow++) {
            binPos = pow(2, binPosPow);
            if (binPos & n) {
               outX = x + (idx_x - CDB_x_mean) + 0.5;
               outY = y + (binPosPow - CDB_y_mean) + 0.5;
               hitsOut->appendPoint(outX, outY, layer, edep/CDB_clusterSize, eventID, isSecondary); 
            }
         }
      idx_x++;
      }
   }
   
   return hitsOut;
}

Tracks * getTracksFromClusters(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy) {
   DataInterface   * di = new DataInterface();
   Int_t             nClusters = kEventsPerRun * 5 * nLayers;
   Int_t             nTracks = kEventsPerRun * 2;
   Bool_t            breakSignal = false;
   Clusters        * clusters = nullptr;
   Tracks          * tracks = nullptr;
   Tracks          * allTracks = new Tracks(nTracks * Runs);
   TRandom3        * gRandom = new TRandom3(0);

   TH2C            * hHits = new TH2C("hHits", "", 500, 4250, 4750, 500, 2000, 2500);
   TH2C            * hDiffusedHits = new TH2C("hDiffusedHits", "", 500, 4250, 4750, 500, 2000, 2500);
   TH2C            * hClusters = new TH2C("hClusters", "", 500, 4250, 4750, 500, 2000, 2500);

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
         di->getMCClusters(i, nullptr, hits);
         hits->sortHits(); 
         diffusedHits = diffuseHits(gRandom, hits);
         diffusedHits->sortHits();
         diffusedHits->makeLayerIndex();
         clusters = diffusedHits->findClustersFromHits();

         for (int j=0; j<hits->GetEntriesFast(); j++) {
            if (hits->getLayer(j) == 7) {
               cout << *hits->At(j) << endl;
               hHits->Fill(hits->getX(j), hits->getY(j));
            }
         }
         
         for (int j=0; j<diffusedHits->GetEntriesFast(); j++) {
            if (diffusedHits->getLayer(j) == 7) {
               hDiffusedHits->Fill(diffusedHits->getX(j), diffusedHits->getY(j));
            }
         }
         
         for (int j=0; j<clusters->GetEntriesFast(); j++) {
            if (clusters->getLayer(j) == 7) {
               hClusters->Fill(clusters->getX(j), clusters->getY(j));
            }
         }

         delete hits;
         delete diffusedHits;
      }
      else {
         di->getMCClusters(i, clusters);
      }

      t1.Stop();
   
      showDebug("Finding calorimeter tracks\n");
      t2.Start(false);
      if (kDoTracking) {
         showDebug("sortClusters...");
         clusters->sortClusters();
         showDebug("ok!\n Start tracking (kWeightedRecursive)\n");
         tracks = clusters->findTracksWithRecursiveWeighting();
         if (i%10 == 0)    printf("Found %d tracks from %d Clusters in run %.2f-%d.\n", tracks->GetEntriesFast(), clusters->GetEntriesFast(), run_energy, i);
      }
      else {
         tracks = clusters->findCalorimeterTracksWithMCTruth();
      }
      t2.Stop();
    
      showDebug("propagateSecondaryStatus..."); 
      tracks->propagateSecondaryStatus();
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
      
      t3.Start(false);
      showDebug("ok!\nRemoving tracks leaving detector...");
      tracks->removeTracksLeavingDetector();
      showDebug("ok\ncompress tracks and clusters...");      
      t3.Stop();
      t4.Start(false);
      tracks->Compress();
      tracks->CompressClusters();
      showDebug("ok\n");
      
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

   TCanvas *c = new TCanvas;
   gStyle->SetPalette(56);
   gStyle->SetOptStat(0);

   c->Divide(3,1,1e-5,1e-5);
   c->cd(1);
   int yfrom = 2180;
   int yto = 2280;
   int xfrom = 4400;
   int xto = 4500;
   hHits->SetTitle("Monte Carlo pixel hits;X position [# of pixels];Y position [# of pixels]");
   hHits->GetXaxis()->SetRangeUser(xfrom, xto);
   hHits->GetYaxis()->SetRangeUser(yfrom, yto);
   hHits->Draw("col");
   hDiffusedHits->SetTitle("Diffused pixel hits;X position [# of pixels];Y position [# of pixels]");
   hDiffusedHits->GetXaxis()->SetRangeUser(xfrom, xto);
   hDiffusedHits->GetYaxis()->SetRangeUser(yfrom, yto);
   c->cd(2);
   hDiffusedHits->Draw("col");
   c->cd(3);
   hClusters->SetTitle("Identified cluster positions;X position [# of pixels];Y position [# of pixels]");
   hClusters->GetXaxis()->SetRangeUser(xfrom, xto);
   hClusters->GetYaxis()->SetRangeUser(yfrom, yto);
   hClusters->Draw("col");

   printf("Timing: Cluster retrieval: %.3f s. Tracking: %.3f s. Track improvements: %.3f s. Logistics: %.3f s. Sort TCA by layer: %.3f s\n", t1.CpuTime(), t2.CpuTime(), t3.CpuTime(), t4.CpuTime(), t6.CpuTime());

   delete di;

   return allTracks;
}

#endif
