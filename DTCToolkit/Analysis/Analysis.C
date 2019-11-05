#ifndef Analysis_cxx
#define Analysis_cxx

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TH2.h>
#include <TH3.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TPaveStats.h>
#include <TView.h>
#include <TLeaf.h>
#include <TArrow.h>
#include <TF1.h>
#include <Math/ProbFunc.h>

#include "Analysis/Analysis.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "GlobalConstants/RangeAndEnergyCalculations.h"
#include "GlobalConstants/Misalign.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/DataInterface/DataInterface.h"
#include "HelperFunctions/Tools.h"
#include "HelperFunctions/getTracks.h"

using namespace std;
using namespace DTC;


void drawTracksDeltaTheta(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness) {
   // For each track (found using event ID information), find the change in angle
   // at layer 1, 2, 3, 4, 5, ..., nLayers.

   Int_t    layers = 50;
   Int_t    lastActivatedLayer = 0;
   Bool_t   drawIndividualLayers = true;
   Track   *thisTrack = nullptr;
   Float_t  angle, totalAngle, y2, y1, y0, x2, x1, x0;
   Float_t  entering[3] = {};
   Float_t  leaving[3] = {};
   Float_t  dotproduct, scalarproduct;
   Double_t mus[50] = {};
   Double_t sigmas[50] = {};
   TCanvas *cSum = new TCanvas("cSum", "Total angular spread at increasing depth", 1200, 900);
   TH2F    *hAngleSumInAllLayers = new TH2F("hAngleSumInAllLayers", ";Layer number;Total angular spread [mrad]", layers, 0, layers-1, 100, 0.18, 2.8);
   
   kDoTracking = false;
   kEventsPerRun = 50000;
   run_degraderThickness = degraderThickness;
   run_energy = energy;

   if (kUseDegrader) run_energy = getEnergyAtWEPL(energy, degraderThickness);
   BinLogY(hAngleSumInAllLayers);

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, run_energy);

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      totalAngle = 0;

      for (Int_t layer=0; layer<layers; layer++) {
         if (thisTrack->GetEntriesFast() - 2 <= layer) {
            lastActivatedLayer = fmax(layer, lastActivatedLayer);
            break;
         }

         x2 = thisTrack->getXmm(layer+1);
         x1 = thisTrack->getXmm(layer);
         x0 = (layer > 0) ? thisTrack->getXmm(layer-1) : thisTrack->getXmm(layer); // parallel projection if layer=0

         y2 = thisTrack->getYmm(layer+1);
         y1 = thisTrack->getYmm(layer);
         y0 = (layer > 0) ? thisTrack->getYmm(layer-1) : thisTrack->getYmm(layer); // parallel projection if layer=0
         
         entering[0] = x1-x0;
         entering[1] = y1-y0;
         entering[2] = dz;

         leaving[0] = x2-x1;
         leaving[1] = y2-y1;
         leaving[2] = dz;

         // DOT PRODUCT RULE
         // dot(a,b) = |a| |b| cos theta
         
         scalarproduct = sqrt(pow(entering[0], 2) + pow(entering[1], 2) + pow(entering[2], 2)) * sqrt(pow(leaving[0], 2) + pow(leaving[1], 2) + pow(leaving[2], 2));
         dotproduct = entering[0] * leaving[0] + entering[1] * leaving[1] + entering[2] * leaving[2];
         angle = acos(dotproduct / scalarproduct);
         totalAngle = sqrt(pow(1000 * angle, 2) + pow(totalAngle, 2));

         hAngleSumInAllLayers->Fill(layer, totalAngle);
      }
   }

   gStyle->SetOptStat(0);

   hAngleSumInAllLayers->GetXaxis()->SetTitleFont(22);
   hAngleSumInAllLayers->GetYaxis()->SetTitleFont(22);
   hAngleSumInAllLayers->GetZaxis()->SetTitleFont(22);
   hAngleSumInAllLayers->GetXaxis()->SetTitleSize(0.05);
   hAngleSumInAllLayers->GetYaxis()->SetTitleSize(0.05);
   hAngleSumInAllLayers->GetZaxis()->SetTitleSize(0.05);
   hAngleSumInAllLayers->GetXaxis()->SetLabelFont(22);
   hAngleSumInAllLayers->GetYaxis()->SetLabelFont(22);
   hAngleSumInAllLayers->GetXaxis()->SetLabelSize(0.05);
   hAngleSumInAllLayers->GetYaxis()->SetLabelSize(0.05);

   hAngleSumInAllLayers->Draw("COLZ");

   Double_t layersArray[50] = {};
   Double_t threeSigmaArray[50] = {};
   for (int i=0; i<50; i++) {
      layersArray[i] = i;
      threeSigmaArray[i] = sigmas[i];
   }

   TPolyLine *pline = new TPolyLine(lastActivatedLayer-2, layersArray, threeSigmaArray);
   pline->SetLineColor(kRed);
   pline->SetLineWidth(3);
   pline->Draw();

   TLine *l190 = new TLine(0, 190, 45, 190);
   TLine *l300 = new TLine(0, 300, 45, 300);

   l190->SetLineColor(kBlue-4);
   l300->SetLineColor(kBlue-4);
   l190->SetLineWidth(3);
   l300->SetLineWidth(3);
   l190->Draw();
   l300->Draw();

}
   
void drawTracksDeltaThetaEachLayer(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness) {
   // Make one histogram for each layer (starting at layer 1)

   Int_t layers = 50;
   Int_t lastActivatedLayer = 0;
   Bool_t drawIndividualLayers = true;
   Track *thisTrack = nullptr;
   Float_t angle, totalAngle, y2, y1, y0, x2, x1, x0;
   Float_t entering[3] = {};
   Float_t leaving[3] = {};
   Float_t dotproduct, scalarproduct;
   Double_t mus[50] = {};
   Double_t sigmas[50] = {};

   kDoTracking = false;
   kEventsPerRun = 50000;
   run_degraderThickness = degraderThickness;
   run_energy = energy;
   if (kUseDegrader) run_energy = getEnergyAtWEPL(energy, degraderThickness);

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, run_energy);

   vector<TH1F*> *hAngleDifference = new vector<TH1F*>;
   vector<TCanvas*> *cCanvases = new vector<TCanvas*>;
   hAngleDifference->reserve(layers);
   cCanvases->reserve(layers);

   for (Int_t layer=0; layer<layers; layer++) {
      hAngleDifference->push_back(new TH1F(Form("hAngleDifference_layer_%i",layer), Form("Angular spread in layer %d for %.0f mm absorbator;Angular spread [rad];Entries",layer, degraderThickness), 1000, 0, 500));
   }


   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      totalAngle = 0;

      for (Int_t layer=0; layer<layers; layer++) {
         if (thisTrack->GetEntriesFast() - 2 <= layer) {
            lastActivatedLayer = fmax(layer, lastActivatedLayer);
            break;
         }

         x2 = thisTrack->getXmm(layer+1);
         x1 = thisTrack->getXmm(layer);
         x0 = (layer > 0) ? thisTrack->getXmm(layer-1) : thisTrack->getXmm(layer); // parallel projection if layer=0

         y2 = thisTrack->getYmm(layer+1);
         y1 = thisTrack->getYmm(layer);
         y0 = (layer > 0) ? thisTrack->getYmm(layer-1) : thisTrack->getYmm(layer); // parallel projection if layer=0
         
         entering[0] = x1-x0;
         entering[1] = y1-y0;
         entering[2] = dz;

         leaving[0] = x2-x1;
         leaving[1] = y2-y1;
         leaving[2] = dz;

         // DOT PRODUCT RULE
         // dot(a,b) = |a| |b| cos theta
         
         scalarproduct = sqrt(pow(entering[0], 2) + pow(entering[1], 2) + pow(entering[2], 2)) * sqrt(pow(leaving[0], 2) + pow(leaving[1], 2) + pow(leaving[2], 2));
         dotproduct = entering[0] * leaving[0] + entering[1] * leaving[1] + entering[2] * leaving[2];
         angle = acos(dotproduct / scalarproduct);
         totalAngle = sqrt(pow(1000 * angle, 2) + pow(totalAngle, 2));

         hAngleDifference->at(layer)->Fill(totalAngle);
      }
   }

   for (Int_t layer=0; layer<lastActivatedLayer; layer++) {
      cCanvases->push_back(new TCanvas(Form("canvas_%d", layer), Form("LAYER %d", layer), 1200, 900));
   }

   TH1F *h = nullptr;
   printf("Layer mean_angle sigma_angle\n");
   TLine *indLine1 = nullptr;
   TLine *indLine2 = nullptr;
   for (Int_t layer=0; layer<lastActivatedLayer; layer++) {
      cCanvases->at(layer)->cd();
      h = hAngleDifference->at(layer);
      h->SetFillColor(kBlue-4);
      h->Draw();
      TF1 *fit = new TF1("fit", "gaus");
      h->Fit(fit, "Q");

      printf("%d %.3f %.3f\n", layer, fit->GetParameter(1), fit->GetParameter(2));
      
      mus[layer] = h->GetBinCenter(h->GetMaximumBin());
      
      Float_t totalArea = h->Integral();
      Int_t binAt = 0;
      for (int i=h->GetNbinsX(); i>0; i--) {
         if (h->Integral(i, h->GetNbinsX()) > (1-0.97752) * totalArea) { // 2 sigma
            binAt = i;
            break;
         }
      }
      
      h->GetXaxis()->SetRange(h->GetMaximumBin(), 1000);
      sigmas[layer] = h->GetXaxis()->GetBinCenter(binAt);
      
      indLine1 = new TLine(mus[layer], 0, mus[layer], 5000);
      indLine2 = new TLine(sigmas[layer], 0, sigmas[layer], 5000);
      indLine1->SetLineColor(kRed);
      indLine2->SetLineColor(kBlack);
      indLine1->Draw();
      indLine2->Draw();

      delete fit;
   }
}

void drawRadiograph(Int_t nparticles, Float_t energy) {
   run_energy = energy;
   run_degraderThickness = 160;
   Float_t X0, Y0, Xplane, Yplane, P0x, P0y, WEPL;
   Track * thisTrack = nullptr;
   TGraphErrors *tge = nullptr;
   Float_t resolution = 0.25; // mm
//   TH2F * radiograph = new TH2F("radiograph", "27x13.5 cm^2 radiograph of 16 cm ball;X position [mm];Y position[mm];WEPL", 270/resolution, -135, 135, 135/resolution, -67.5, 67.5);
//   TH2I * radiographNorm = new TH2I("radiographNorm", "27x13.5 cm^2 radiograph of 16 cm ball;X position [mm];Y position[mm];WEPL", 270/resolution, -135, 135, 135/resolution, -67.5, 67.5);
   TH2F * radiograph = new TH2F("radiograph", "16 cm water - #sigma=3 mm - projected on hull;X position [mm];Y position[mm];WEPL", 30/resolution, -15, 15, 30/resolution, -15, 15);
   TH2I * radiographNorm = new TH2I("radiographNorm", "27x13.5 cm^2 radiograph of 16 cm ball;X position [mm];Y position[mm];WEPL", 30/resolution, -15, 15, 30/resolution, -15, 15);

   TH1F * oneDradiograph = new TH1F("1D_radiograph", "16 cm water - 3 mm sigma beam - projected on hull;Y position [mm]", 30/resolution, -15, 15);
   TH1F * oneDradiographNorm = new TH1F("1D_radiographNorm", "title;Y position [mm]", 30/resolution, -15, 15);
   
   Tracks * tracks = loadOrCreateTracks(1, int(nparticles/50), kMC, energy);

   printf("Found %d tracks.\n", tracks->GetEntriesFast());
//   tracks->removeTracksEndingInGapRegion();
   tracks->removeHighAngleTracks(75);
//   tracks->removeThreeSigmaShortTracks();
//   tracks->removeNuclearInteractions();
   tracks->Compress();

   Float_t d_plane = 20;

   for (int i=0; i<=tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (!thisTrack->At(0) || !thisTrack->At(1)) continue;

      tge = (TGraphErrors*) thisTrack->doTrackFit(false, kUseCSDA);
      WEPL = getWEPLFromTL(thisTrack->getFitParameterRange());


      X0 = thisTrack->getXmm(0);
      Y0 = thisTrack->getYmm(0);
      P0x = (thisTrack->getXmm(1) - thisTrack->getXmm(0)) / dz;
      P0y = (thisTrack->getYmm(1) - thisTrack->getYmm(1)) / dz;

      Xplane = X0 - P0x * d_plane;
      Yplane = Y0 - P0y * d_plane;

      radiograph->Fill(Xplane, Yplane, getWEPLFromEnergy(energy) - WEPL);
      radiographNorm->Fill(Xplane, Yplane);

      oneDradiograph->Fill(Yplane, getWEPLFromEnergy(energy) - WEPL);
      oneDradiographNorm->Fill(Yplane);
   }

   radiograph->Divide(radiographNorm);
   oneDradiograph->Divide(oneDradiographNorm);

   radiograph->Draw("COLZ");
   gStyle->SetOptStat(0);

   TCanvas *c2 = new TCanvas();
   oneDradiograph->Draw();
}

void getRangeFromRawBeam(Float_t energy) {
   run_energy = energy;
   run_degraderThickness = 0;
   Float_t  ranges[1000];
   Float_t  rangesBootstrap[5000];
   Float_t  errorsBootstrap[5000];
   Float_t  newRanges[1000];
   Int_t    nTracks;
   Track   *thisTrack = nullptr;
   Float_t  fitRange;
   Float_t  sum = 0;
   Float_t  sum2 = 0;
   Float_t  mean = 0;
   Float_t  x;
   TGraphErrors *tge = nullptr;

   Tracks * tracks = loadOrCreateTracks(1, 5, kMC, energy);
   printf("Found %d tracks.\n", tracks->GetEntriesFast());
   tracks->removeHighAngleTracks(75);
   tracks->removeThreeSigmaShortTracks();
   tracks->removeNuclearInteractions();
   tracks->Compress();

   nTracks = tracks->GetEntries();
   Int_t nTracksSeen = 0;
   for (int i=0; i<=nTracks; i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) {
         printf("!thisTrack %d\n", i);
         continue;
      }

      tge = (TGraphErrors*) thisTrack->doTrackFit(false, kUseCSDA); // (bool isScaleVariable, bool useTrackLength (~ CSDA))
      fitRange = getWEPLFromTL(thisTrack->getFitParameterRange());
      ranges[i] = fitRange;
      sum += fitRange;
      nTracksSeen++;
   }
   
   mean = sum/nTracksSeen;

   TH1I    *hBootstrap  = new TH1I("hBootstrap", "Bootstrap mean; Ranges [mm WEPL];Frequency", 300, mean-50, mean+50);
   TH1I    *hBootstrapStdDev  = new TH1I("hBootstrapStdDev", "Bootstrap standard deviation; Ranges [mm WEPL]", 300, 0, 50);
   TH1I    *hRanges     = new TH1I("hRanges", "Original range histogram; Ranges [mm WEPL];Frequency", 300, mean-50, mean+50);

   for (int i=0; i<=nTracks; i++) {
      hRanges->Fill(ranges[i]);
   }

   TRandom3 *gRand = new TRandom3(0);
   Int_t nextIdx, num_tracks_in_list;
   Float_t bs_mean, bs_mean2, bs_error, bs_stddev, bs_stddeverror;
   for (int i=0; i<5000; i++) {
      sum = 0; sum2 = 0;
      num_tracks_in_list = 0;

      while (num_tracks_in_list < nTracksSeen) { // make bootstrapped sample
         nextIdx = gRand->Integer(nTracksSeen);
         newRanges[num_tracks_in_list] = ranges[nextIdx];
         num_tracks_in_list++;
      }

      for (int i=0; i<nTracksSeen; i++) { // calculate statistics from bootstrapped sample
         x = newRanges[i];
         sum += x;
         sum2 += x*x;
      }

      bs_mean = sum/nTracksSeen;
      bs_mean2 = sum2/nTracksSeen;
      bs_error = sqrt((bs_mean2 - pow(bs_mean,2)));

      rangesBootstrap[i] = bs_mean;
      errorsBootstrap[i] = bs_error;

      hBootstrap->Fill(bs_mean);
      hBootstrapStdDev->Fill(bs_error);
   }

   TCanvas *c = new TCanvas();
   hBootstrap->Draw();

   TCanvas *c2 = new TCanvas();
   hRanges->Draw();

   Float_t bootstrapStdDev = hBootstrap->GetStdDev();
   Float_t bootstrapMean = hBootstrap->GetMean();
   printf("Mean range is %.2f +- %.2f mm.\n", mean, bootstrapStdDev);

   // 2x sample median - bootstrap mean
   Double_t sampleMedian;
   Double_t quantile = 0.5;
   hRanges->GetQuantiles(1, &sampleMedian, &quantile); // WHAT a bothersome method to get the median
   printf("The sample median is %.2f mm.\n", sampleMedian);
   printf("The corrected mean range is %.2f +- %.2f mm.\n", 2*sampleMedian - hBootstrap->GetMean(), bootstrapStdDev);
/*
   for (int i=0; i<=5000; i++) {
      int num_tracks_in_list = 0;
      sum = 0;
      while (num_tracks_in_list++ < nTracks) {
         nextIdx = gRand->Integer(nTracks);
         sum += pow(ranges[nextIdx] - bootstrapMean, 2);
      }
      sum /= nTracksSeen;
      hBootstrapStdDev->Fill(sqrt(sum));
   }
*/
   TCanvas *c3 = new TCanvas();
   hBootstrapStdDev->Draw();

   hBootstrapStdDev->Fill(sqrt(sum));
   Float_t bootstrapStdDevMean = hBootstrapStdDev->GetMean();
   Float_t bootstrapStdDevStdDev = hBootstrapStdDev->GetStdDev();
   printf("The bootstrapped range uncertainty is %.2f +- %.2f mm.\n", bootstrapStdDevMean, bootstrapStdDevStdDev);
   printf("The corrected range uncertanity is %.2f +- %.2f mm.\n", 2*hRanges->GetStdDev() - bootstrapStdDevMean, bootstrapStdDevStdDev);

}

void getTracksReconstructionEfficiency(Int_t dataType, Float_t energy, Float_t degraderThickness) {
   Int_t nRuns = 0;

   run_energy = energy;
   run_degraderThickness = degraderThickness;
   
   if (kUseDegrader) {
      run_energy = getEnergyAtWEPL(energy, degraderThickness);
   }

   Int_t nRunArray[12] = {3,4,5,8,16,32,64,128,181,256,512,1024};

   for (Int_t i=0; i<12; i++) { // 1 -> 30
      nRuns = nRunArray[i];

      kEventsPerRun = nRuns;
      Float_t factor = 2;

      Int_t totalNumberOfRuns = 10000 / kEventsPerRun;
      if (totalNumberOfRuns < 1) totalNumberOfRuns = 1;
      if (totalNumberOfRuns > 10000) totalNumberOfRuns = 10000;

      Tracks * tracks = loadOrCreateTracks(1, totalNumberOfRuns, dataType, energy);
      tracks->removeHighAngleTracks(75);
      tracks->removeThreeSigmaShortTracks();
      tracks->removeNuclearInteractions();
      tracks->extrapolateToLayer0();

      Track *thisTrack;
      Int_t EID, thisEID;
      Int_t nTotal = 0;
      Int_t nFirstAndLastAllTracks = 0;

      for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
         thisTrack = tracks->At(j);
         if (!thisTrack) continue;
         nTotal++;

         if (thisTrack->isFirstAndLastEventIDEqual() && tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), 1) == 0) {
            nFirstAndLastAllTracks++;
         }
      }

      Float_t ratioFirstAndLastAllTracks = (float) nFirstAndLastAllTracks / nTotal;
      Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;

      ofstream file2(Form("OutputFiles/lastLayerCorrect_different_nRuns_%.0f.csv", kAbsorberThickness), ofstream::out | ofstream::app);
      file2 << readoutAbsorber << " " << nRuns << " " << " " << ratioFirstAndLastAllTracks << endl;
      file2.close();
      
      delete tracks;
   }
}

void drawClusterShapes(Int_t Runs, Bool_t dataType, Float_t energy, Float_t degraderThickness) {
   // get vector of TH2F's, each with a hits distribution and cluster size
   // dataType = kMC (0) or kData (1)

   // See also ../Classes/Hit/findClusters.C

   showDebug("Initializing function...\n");
   run_energy = energy;
   run_degraderThickness = degraderThickness;

   Int_t useDataHits = false; // cpu saver?
   kDataType = dataType;

   Int_t nRows = 6;
   Int_t nRepeats = 10;
   Int_t nN = nRows * nRepeats;
   vector<TH2C> *hClusterMaps = new vector<TH2C>;
   hClusterMaps->reserve(nN);

   showDebug("Creating TH2Fs...");
   for (Int_t i=0; i<nN; i++) {
      hClusterMaps->push_back(TH2C(Form("hClusterMap_%i",i), "", 11, 0, 11, 11, 0, 11));
   }
   showDebug("OK!\n");

   Int_t nClusters = kEventsPerRun * 5;
   Int_t nHits = kEventsPerRun * 50;
   Int_t nTracks = kEventsPerRun * 2;

   showDebug("Creaing vectors...");
   DataInterface *di = new DataInterface();
   CalorimeterFrame *cf = new CalorimeterFrame();
   Hits *hits = new Hits(nHits);
   vector<Hits*> * tempClusterHitMap;
   vector<Hits*> * clusterHitMap = new vector<Hits*>;
   clusterHitMap->reserve(Runs*15*kEventsPerRun);
   showDebug("OK\n");
   Hits * dataHits = nullptr;

   if (useDataHits) {
      showDebug("Reserving " << Runs * 25 * nLayers * kEventsPerRun << " Hit objects in Hits\n");
      dataHits = new Hits(Runs*25*nLayers*kEventsPerRun);
   }

   for (Int_t i=0; i<Runs; i++) {
      if (dataType == kMC) { // Use Monte Carlo data
         showDebug("getMCFrame...");
         di->getMCFrame(i, cf);
         showDebug("OK!\nDiffuseFrame...");
         cf->diffuseFrame(new TRandom3(0)); // Model the cluster diffusion process
         showDebug("OK!\nfindHits...");
         hits = cf->findHits();
         showDebug("OK!\nfindClustersHitMap()...");
         tempClusterHitMap = hits->findClustersHitMap();
         showDebug("OK!\n");
      }

      else if (dataType == kData && useDataHits == false) { // Use experimental data (122, 140, 150, 160, 170, 180, 188 MeV)
         showDebug("getDataFrame...");
         di->getDataFrame(i, cf, energy);
         showDebug("OK!\nfindHits...");
         hits = cf->findHits();
         showDebug("OK!\nfindClustersHitMap()...");
         tempClusterHitMap = hits->findClustersHitMap();
         showDebug("OK!\n");
      }

      else if (dataType == kData && useDataHits == true) {
         di->getDataHits(i, dataHits, energy);
         showDebug("Found " << dataHits->GetEntriesFast() << " hits in exp. data\n");
         tempClusterHitMap = dataHits->findClustersHitMap();
      }  

      else {
        std::cerr << "Please choose between dataType = kMC (0) or kData (1).\n" << endl;
        exit(1);
      }
  
      showDebug("Add to clusterHitMap...");
      for (UInt_t j=0; j<tempClusterHitMap->size(); j++) {
         clusterHitMap->push_back( tempClusterHitMap->at(j) );
      }
      showDebug("OK!\n");

      cf->Reset();
   }
   
   delete hits;
   delete cf;
   delete di;
   
   // Here it is possible to access and modify the cluster shapes
   // Each cluster is stored as a Hits (Hit collection) pointer in the vector collection clusterHitMap.
   // To loop through each cluster use for (i=0; i<clusterHitMap->size(); i++) { Hits * myCluster = clusterHitMap->at(i); myCluster->....; }
   // E.g. to count the number of hits in each cluster:
   //
   // for (int i=0; i<clusterHitMap->size(); i++) {
   //    Hits *myCluster = clusterHitMap->at(i);
   // if (!myCluster) continue;
   // int counter = 0;
   // for (int j=0; j<myCluster->GetEntriesFast(); j++) {
   //    Hit *myHit = myCluster->At(j);
   //    if (!myHit) continue;
   //    counter++;
   // }
   // }
   //
   
   // fill hClusterMaps with cluster shapes from clusterHitMap
   // sizes 3-5 in first row
   // 6-8 in second row
   // eg. up to 33-35 in 11th row
   
   Int_t size_from, size_to, nInRow, x, y;

   Int_t nTotal = 0;
   Int_t kColor = 1;

   for (Int_t i=0; i<nRows; i++) {
      size_from = i*9;
      size_to = size_from + 8;
      nInRow = 0;

      for (UInt_t j=0; j<clusterHitMap->size(); j++) {
         Int_t csize = clusterHitMap->at(j)->GetEntriesFast();
         if (j == 284) continue;
         if (csize >= size_from && csize <= size_to) {
            // plot the cluster in the j'th row
            for (Int_t k=0; k<clusterHitMap->at(j)->GetEntriesFast(); k++) {
               x = clusterHitMap->at(j)->getX(k);
               y = clusterHitMap->at(j)->getY(k);
               if (useDataHits) {
                  // set color according to eventID of k; be smart
               }
               hClusterMaps->at(nTotal).SetBinContent(x,y, kColor);
            } // end loop through all points
            hClusterMaps->at(nTotal).SetTitle(Form("%i", csize));
            printf("row = %d, repeat = %d, j = %d\n", i, nInRow, j);
            nInRow++; nTotal++;
            if (nInRow >= nRepeats) break; // stop looping over clusters now
         }
      }
      if (nInRow < nRepeats) {
         cout << "Only " << nInRow << " clusters with size [" << i*3 << "," << i*3+2 << "] with i = " << i << ". Setting nTotal from " << nTotal << " to " << (i+1)*nRepeats << ".\n";
         nTotal = (i+1)*nRepeats;
      }
   }

   // draw the canvas
   
   TCanvas *c = new TCanvas("c", "c", 1000, 800);
   c->Divide(nRepeats, nRows, 0, 0);
   for (Int_t i=0; i<nN; i++) {
      c->cd(i+1);
      gStyle->SetTitleY(0.93);
      gStyle->SetTitleX(0.85);
      hClusterMaps->at(i).Draw("same, COL,ah,fb,bb");
      gStyle->SetOptStat(0);
      gPad->Update();
      TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
      if (title) {
         title->SetTextFont(22);
         title->SetTextSize(0.18);
      }
      gPad->Modified();
   }
}

void drawFitScale(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   TCanvas *cScale = new TCanvas("cScale", "Scale histogram", 1400, 1000);
   run_energy = energy;
   kDataType = dataType;
   
   TH1F *hScale = new TH1F("hScale", "Scale histogram", 800, 0, 50);
   TGraphErrors *outputGraph;

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   
   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;
      if (thisTrack->doesTrackEndAbruptly()) continue;

      outputGraph = (TGraphErrors*) thisTrack->doTrackFit(true, kUseCSDA);
      if (!outputGraph) continue;
      
      Float_t fitScale = thisTrack->getFitParameterScale();
      
      hScale->Fill(fitScale);
      delete outputGraph;
   }

   cScale->cd();
   hScale->SetYTitle("Number of protons");
   hScale->SetXTitle("Parameter 1 of fit (SCALE)");
   hScale->Draw();
}  

void drawTracksDepthDose(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness, Int_t eventsPerRun, Bool_t doTracking, Bool_t excludeNuclearInteractions) {
   run_degraderThickness = degraderThickness;
   run_energy = energy;
   kDataType = dataType;
   kEventsPerRun = eventsPerRun;
   
   if (kUseDegrader) {
      run_energy = getEnergyFromDegraderThickness(degraderThickness);
      printf("Using degrader, expecting nominal residual energy %.2f MeV\n", run_energy);
   }
   
   Bool_t         removeHighAngleTracks = true;
   Bool_t         removeNuclearInteractions = false;
   Bool_t         removeShortTracks = false;
   Float_t        fitRange, fitScale, fitError;
   Int_t          nCutDueToTrackEndingAbruptly = 0;
   Int_t          nPlotX = 3, nPlotY = 3;
   Int_t          skipPlot = 0;
   Int_t          fitIdx = 0, plotSize = nPlotX*nPlotY;
   TGraphErrors * outputGraph;
   char         * sDataType = getDataTypeChar(dataType);
   char         * sMaterial = getMaterialChar();
   TCanvas      * cGraph = new TCanvas("cGraph", "Fitted data points", nPlotX*500, nPlotY*500);

   cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
   cGraph->cd();
   gPad->SetBorderMode(0); gStyle->SetFrameBorderMode(0);
   gPad->SetTickx(1); gPad->SetTicky(1);
   gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.05);
   gPad->SetBottomMargin(0.05);
   gPad->SetLeftMargin(0.15);

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, run_energy,15,15);
   tracks->removeEmptyTracks();

   if (removeHighAngleTracks) {
      tracks->removeHighAngleTracks(100);
   }

   if (removeNuclearInteractions) {
      tracks->removeNuclearInteractions();
   }
   
   if (removeShortTracks) {
      tracks->removeThreeSigmaShortTracks();
   }

   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      if (j < skipPlot) continue;

      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;

      // Do track fit, extract all parameters for this track
      outputGraph = (TGraphErrors*) thisTrack->doTrackFit(false, kUseCSDA); // (bool isScaleVariable, bool useTrackLength (~ CSDA))
      if (!outputGraph) continue;

      fitRange = thisTrack->getFitParameterRange();
      fitScale = thisTrack->getFitParameterScale();
      fitError = quadratureAdd(thisTrack->getFitParameterError(), dz*0.28867); // latter term from error on layer position

      if (fitIdx < plotSize) {
         drawIndividualGraphs(cGraph, outputGraph, fitRange, fitScale, fitError, fitIdx++);
         if (thisTrack->Last()->isSecondary()) {
            outputGraph->SetTitle("Secondary track");
         }
         else {
            outputGraph->SetTitle("Primary track");
         }
         printf("Drawing plot number %d.\n", fitIdx);
      }
   
      else break;
   }
}

void drawTracksRangeHistogram(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness, Int_t eventsPerRun, Int_t outputFileIdx, Bool_t drawFitResults, Bool_t doTracking, Bool_t excludeNuclearInteractions, Int_t skipTracks) {
   run_degraderThickness = degraderThickness;
   run_energy = energy;
   kEventsPerRun = eventsPerRun;
   kDoTracking = doTracking;
   kFilterNuclearInteractions = excludeNuclearInteractions;
   
   if (kUseDegrader) {
      run_energy = getEnergyFromDegraderThickness(degraderThickness);
   }

   printf("Using water degrader of thickness %.0f mm, the initial energy of %.0f MeV is reduced to %.1f MeV.\n", degraderThickness, energy, run_energy);

   kDataType = dataType;
   Bool_t         removeHighAngleTracks = true;
   Bool_t         removeNuclearInteractions = true;
   Bool_t         drawVerticalLayerLines = false;

   Float_t        finalEnergy = 0;
   Float_t        fitRange, fitScale, fitError;
   Int_t          nCutDueToTrackEndingAbruptly = 0;
   TGraphErrors * outputGraph;
   char         * sDataType = getDataTypeChar(dataType);
   char         * sMaterial = getMaterialChar();
   char         * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", run_energy, sMaterial, sDataType);

   Int_t nEnergyBins = getUnitFromEnergy(run_energy);
   gStyle->SetOptTitle(0);

   if (kUseDegrader) {
      hTitle = Form("Fitted energy of a %.0f MeV nominal beam on %s DTC w/%.1f mm water degrader", energy, sMaterial, degraderThickness);
   }

   Float_t lowHistogramLimit = getUnitFromEnergy(0);
   Float_t highHistogramLimit = getUnitFromEnergy(run_energy)*1.2 + 10;
   if (isnan(highHistogramLimit)) highHistogramLimit = getUnitFromEnergy(run_energy) + 30;
   TH1F * hFitResults = new TH1F("fitResult", hTitle, fmax(nEnergyBins,100), lowHistogramLimit, highHistogramLimit);
   TH1F * hLastLayer = new TH1F("hLastLayer", "Last Layer;Layer;Entries", 50, 0, 50);
 
   printf("Using material: %s\n", sMaterial);
   printf("Histogram limits: %.2f to %.2f.\n", lowHistogramLimit, highHistogramLimit);
   printf("At energy %.0f, expecting range %.2f mm and WEPL %.2f mm.\n", run_energy, getTLFromEnergy(run_energy), getWEPLFromEnergy(run_energy));
   printf("This corresponds to a WEPL factor of %.2f.\n", getWEPLFactorFromEnergy(run_energy));
   cout << "Correcting for aluminum plate: " << kIsAluminumPlate << endl;
   cout << "Correcting for scintillators: " << kIsScintillator << endl;

   hFitResults->SetLineColor(kBlack); hFitResults->SetFillColor(kGreen-5);

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, run_energy);
   if (removeHighAngleTracks) {
      tracks->removeHighAngleTracks(75);
   }
   if (removeNuclearInteractions) {
      tracks->removeNuclearInteractions();
   }

   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;
      if (j < skipTracks) continue;
    
      if (thisTrack->doesTrackEndAbruptly()) {
         nCutDueToTrackEndingAbruptly++;
      }

      // Do track fit, extract all parameters for this track
      outputGraph = (TGraphErrors*) thisTrack->doTrackFit(false, kUseCSDA); // (bool isScaleVariable, bool useTrackLength (~ CSDA))
      if (!outputGraph) continue;
   
      delete outputGraph;

      fitRange = thisTrack->getFitParameterRange();
      hFitResults->Fill(getUnitFromTL(fitRange));
      hLastLayer->Fill(thisTrack->Last()->getLayer());
   }
  
   // Draw expected gaussian distribution of results from initial energy
   Float_t expectedStraggling = 0, expectedMean = 0, dlayer_down = 0, dlayer = 0;
   Float_t separationFactor = 0.9, nullStraggling = 0;
   Float_t sigma_energy = getSigmaEnergy(run_energy);
   
   expectedMean = getUnitFromEnergy(run_energy);
   expectedStraggling = getUnitStragglingFromEnergy(run_energy, sigma_energy);
   nullStraggling = getUnitStragglingFromEnergy(run_energy, 0);
   
   cout << "OutputUnit is " << kOutputUnit << " and the expected mean value is " << expectedMean 
       << ". The straggling including / excluding energy variation is " << expectedStraggling << " / " << nullStraggling << ".\n";
       
   Float_t means[10] = {};
   Float_t sigmas[10] = {};

   TF1 *gauss = doSimpleGaussianFit(hFitResults, means, sigmas, outputFileIdx);
   Float_t empiricalMean = means[9];
   Float_t empiricalSigma = sigmas[9];
   
   Float_t energySigma = getEnergyFromUnit(empiricalMean +  empiricalSigma/ 2) - getEnergyFromUnit(empiricalMean - empiricalSigma / 2);

   if (drawFitResults) {
      TCanvas * cFitResults = new TCanvas("cFitResults", hTitle, 1000, 1000);
      
      if       (kOutputUnit == kPhysical) hFitResults->SetXTitle("Physical range [mm]");
      else if  (kOutputUnit == kWEPL)     hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
      else if  (kOutputUnit == kUnitEnergy)   hFitResults->SetXTitle("Energy [MeV]");

      hFitResults->SetYTitle("Number of protons");
      hFitResults->GetXaxis()->SetTitleFont(22);
      hFitResults->GetXaxis()->SetLabelFont(22);
      hFitResults->GetYaxis()->SetTitleFont(22);
      hFitResults->GetYaxis()->SetLabelFont(22);
      hFitResults->GetYaxis()->SetTitleOffset(1.5);

      hFitResults->SetTitle("");   
      hFitResults->Draw();
      gPad->Update();

      if (drawVerticalLayerLines) {
         TLine *l = nullptr;
         Float_t line_z = 0;
         for (Int_t i=0; i<65; i++) {
            line_z = getUnitFromTL(getLayerPositionmm(i));
            if (line_z > gPad->GetUxmax()) break;
            l = new TLine(line_z, 0, line_z, hFitResults->GetMaximum()*1.05);
            l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw();
         }
      }

      gPad->Update();
      gStyle->SetOptStat(11);
      TPaveStats *ps = (TPaveStats*) cFitResults->GetPrimitive("stats");
      hFitResults->SetBit(TH1::kNoStats);
      ps->SetX1NDC(0.4176); ps->SetX2NDC(0.9257);
      ps->SetY1NDC(0.7415); ps->SetY2NDC(0.9712);
      ps->SetTextFont(22);
      ps->AddText(Form("Nominal WEPL = %.2f #pm %.2f", expectedMean, expectedStraggling));
      ps->AddText(Form("Calculated WEPL = %.2f #pm %.2f", empiricalMean, empiricalSigma));
      ps->AddText(Form("WEPL deviation = %.2f #pm %.2f", empiricalMean - expectedMean, sqrt(pow(empiricalSigma, 2) - pow(expectedStraggling, 2))));
      cFitResults->Modified();

      cFitResults->SaveAs(Form("OutputFiles/RangeHistogram/%.0f_%.0f.png", kAbsorberThickness, degraderThickness));

      TCanvas *c2 = new TCanvas();
      hLastLayer->Draw();
   }

   delete tracks;
}

void findTracksRangeAccuracy(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness, Int_t eventsPerRun, Int_t outputFileIdx, Bool_t doTracking, Bool_t excludeNuclearInteractions) {
   run_degraderThickness = degraderThickness;
   run_energy = energy;
   kEventsPerRun = eventsPerRun;
   kDoTracking = doTracking;
   kFilterNuclearInteractions = excludeNuclearInteractions;
         
   gROOT->SetBatch(kTRUE);

   if (kUseDegrader) {
      run_energy = getEnergyFromDegraderThickness(degraderThickness);
   }

   printf("Using water degrader of thickness %.0f mm, the initial energy of %.0f MeV is reduced to %.1f MeV.\n", degraderThickness, energy, run_energy);
   printf("-> The expected range is %.2f mm.\n", getTLFromEnergy(run_energy));

   kDataType = dataType;
   Bool_t         removeHighAngleTracks = true;
   Bool_t         removeNuclearInteractions = true;

   Float_t        finalEnergy = 0;
   Float_t        fitRange, fitScale, fitError;
   Int_t          nCutDueToTrackEndingAbruptly = 0;
   Int_t          skipPlot = 0;
   TGraphErrors * outputGraph;
   Float_t        expectedMean = getUnitFromEnergy(run_energy);
   Float_t        expectedStraggling = getUnitStragglingFromEnergy(run_energy, getSigmaEnergy(run_energy));
   Float_t        empiricalMean, empiricalSigma;
   Float_t        means[10] = {};
   Float_t        sigmas[10] = {};
   Float_t        listOfErrorValues[10000] = {};
   Int_t          idxOfErrorValues = 0;

   Int_t          nEnergyBins = getUnitFromEnergy(run_energy);
   Float_t        lowHistogramLimit = getUnitFromEnergy(0);
   Float_t        highHistogramLimit = getUnitFromEnergy(run_energy)*1.4 + 10;
   TH1F         * hFitResults = new TH1F("fitResult", "noDrawHistogram", fmax(nEnergyBins,100), lowHistogramLimit, highHistogramLimit);
 
   Tracks * tracks = loadOrCreateTracks(recreate, Runs*1.7, dataType, run_energy);

   if (removeHighAngleTracks) {
      tracks->removeHighAngleTracks(100);
   }
   if (removeNuclearInteractions) {
      tracks->removeNuclearInteractions();
   }

   Int_t nTracks = 0;
   printf("ntracks = %d\n", tracks->GetEntriesFast());
   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;

      outputGraph = (TGraphErrors*) thisTrack->doTrackFit(false, kUseCSDA); // (bool isScaleVariable, bool useTrackLength (~ CSDA))
      if (!outputGraph) continue; 
      delete outputGraph;
      
      nTracks++;
      if (nTracks > Runs * eventsPerRun) break;

      fitRange = thisTrack->getFitParameterRange();
      hFitResults->Fill(getUnitFromTL(fitRange));
      
      if (nTracks % eventsPerRun == 0) { // DO ANALYSIS
         TF1 *gauss = doSimpleGaussianFit(hFitResults, means, sigmas, outputFileIdx);
         delete gauss;

         empiricalMean = means[9];
         empiricalSigma = sigmas[9];
         
         printf("Run %d of %d ... Range error is %.2f mm WET, with %d protons.\n", nTracks / eventsPerRun, Runs, empiricalMean - expectedMean, eventsPerRun);
         if (!isnan(empiricalMean)) {
            listOfErrorValues[idxOfErrorValues++] = empiricalMean - expectedMean;
         }
         TCanvas *c = new TCanvas();
         hFitResults->Draw();
         c->SaveAs("tcanvas.png");
         delete c;
         hFitResults->Reset();
      }
   }

   delete hFitResults;
   delete tracks;

   Float_t meanError = 0;
   Float_t sigmaError = 0;

   for (Int_t i=0; i<idxOfErrorValues; i++) {
      meanError += listOfErrorValues[i];
   }
   
   meanError /= (idxOfErrorValues);

   for (Int_t i=0; i<idxOfErrorValues; i++) {
      sigmaError += pow(listOfErrorValues[i] - meanError, 2);
   }

   if (idxOfErrorValues == 1) sigmaError = 1e5;
   else sigmaError = sqrt(sigmaError / (idxOfErrorValues));

   printf("Mean error from %d runs is %.2f mm +- %.2f mm. Expected error is %.2f.\n", Runs, meanError, sigmaError, expectedStraggling / sqrt(eventsPerRun));
   ofstream file("OutputFiles/tracksRangeAccuracy.csv", ofstream::out | ofstream::app);
   file << run_energy << " " << Runs << " " << eventsPerRun<< " " << meanError << " " << sigmaError << endl;
   file.close();
}

void draw2DProjection(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   
   Tracks * tracks1 = loadOrCreateTracks(recreate, Runs, dataType, 170);
   Tracks * tracks2 = loadOrCreateTracks(recreate, Runs, dataType, 188);
   Tracks * tracks3 = loadOrCreateTracks(recreate, Runs, dataType, 190);
   Tracks * tracks4 = loadOrCreateTracks(recreate, Runs, dataType, 150);
   Tracks * tracks5 = loadOrCreateTracks(recreate, Runs, dataType, 160);
   tracks1->extrapolateToLayer0();
   tracks2->extrapolateToLayer0();
   tracks3->extrapolateToLayer0();
   tracks4->extrapolateToLayer0();
   tracks5->extrapolateToLayer0();
   
   Int_t hSizeX = nx/12;
   Int_t hSizeY = ny/12;

   Float_t angleCut = 5.;
   Float_t x0, y0, xL, yL, theta0;
   Int_t nPoints = 0;

   Float_t fit_energy;

   TCanvas *c1 = new TCanvas("c1", "Entry points", 1500, 1000);
   c1->Divide(2, 1, 0.001, 0.001);
   TCanvas *c2 = new TCanvas("c2", "1D entry points", 1200, 1200);
   c2->Divide(2,1,0.001,0.001);

   TH2F *normalizeFrame = new TH2F("normalizeFrame", "2D Projection of entry position from multiple proton beams on FOCAL;X position [mm];Y position [mm]", hSizeX, -30, 30, hSizeY, -30, 30);
   TH2F *normalizeFrameLast = new TH2F("normalizeFrameLast", "2D Projection of stopping position from multiple proton beams on FOCAL;X position [mm]; Y posision [mm]", hSizeX, -30, 30, hSizeY, -30, 30);

   TH1F *beamSpec = new TH1F("beamSpec", "1D projection of entry positions from multiple proton beams on FOCAL;X position [mm];Entries", hSizeX, -30, 30);
   TH1F *beamSpecLast = new TH1F("beamSpecLast", "1D projection of exit positions from multiple proton beams on FOCAL;X position [mm];Entries", hSizeX, -30, 30);

   for (Int_t i=0; i<tracks1->GetEntriesFast(); i++) {
      Track *thisTrack = tracks1->At(i);
      if (!thisTrack) continue;
      x0 = thisTrack->getXmm(0);
      y0 = thisTrack->getYmm(0);
      xL = thisTrack->Last()->getXmm();
      yL = thisTrack->Last()->getYmm();
      normalizeFrame->Fill(x0, y0);
      normalizeFrameLast->Fill(xL, yL);
      beamSpec->Fill(x0);
      beamSpecLast->Fill(xL);
   }
   
   for (Int_t i=0; i<tracks2->GetEntriesFast(); i++) {
      Track *thisTrack = tracks2->At(i);
      if (!thisTrack) continue;
      x0 = thisTrack->getXmm(0);
      y0 = thisTrack->getYmm(0);
      xL = thisTrack->Last()->getXmm();
      yL = thisTrack->Last()->getYmm();
      normalizeFrame->Fill(x0, y0);
      normalizeFrameLast->Fill(xL, yL);
      beamSpec->Fill(x0);
      beamSpecLast->Fill(xL);
   }
   
   for (Int_t i=0; i<tracks3->GetEntriesFast(); i++) {
      Track *thisTrack = tracks3->At(i);
      if (!thisTrack) continue;
      x0 = thisTrack->getXmm(0);
      y0 = thisTrack->getYmm(0);
      xL = thisTrack->Last()->getXmm();
      yL = thisTrack->Last()->getYmm();
      normalizeFrame->Fill(x0, y0);
      normalizeFrameLast->Fill(xL, yL);
      beamSpec->Fill(x0);
      beamSpecLast->Fill(xL);
   }
   
   for (Int_t i=0; i<tracks4->GetEntriesFast(); i++) {
      Track *thisTrack = tracks4->At(i);
      if (!thisTrack) continue;
      x0 = thisTrack->getXmm(0);
      y0 = thisTrack->getYmm(0);
      xL = thisTrack->Last()->getXmm();
      yL = thisTrack->Last()->getYmm();
      normalizeFrame->Fill(x0, y0);
      normalizeFrameLast->Fill(xL, yL);
      beamSpec->Fill(x0);
      beamSpecLast->Fill(xL);
   }
   
   for (Int_t i=0; i<tracks5->GetEntriesFast(); i++) {
      Track *thisTrack = tracks5->At(i);
      if (!thisTrack) continue;
      x0 = thisTrack->getXmm(0);
      y0 = thisTrack->getYmm(0);
      xL = thisTrack->Last()->getXmm();
      yL = thisTrack->Last()->getYmm();
      normalizeFrame->Fill(x0, y0);
      normalizeFrameLast->Fill(xL, yL);
      beamSpec->Fill(x0);
      beamSpecLast->Fill(xL);
   }

   delete tracks1;
   delete tracks2;
   delete tracks3;
   delete tracks4;
   delete tracks5;
   
   normalizeFrame->GetXaxis()->SetTitleFont(22);
   normalizeFrame->GetXaxis()->SetLabelFont(22);
   normalizeFrame->GetYaxis()->SetTitleFont(22);
   normalizeFrame->GetYaxis()->SetLabelFont(22);
   
   normalizeFrameLast->GetXaxis()->SetTitleFont(22);
   normalizeFrameLast->GetXaxis()->SetLabelFont(22);
   normalizeFrameLast->GetYaxis()->SetTitleFont(22);
   normalizeFrameLast->GetYaxis()->SetLabelFont(22);
  
   beamSpec->GetXaxis()->SetTitleFont(22);
   beamSpec->GetXaxis()->SetLabelFont(22);
   beamSpec->GetYaxis()->SetTitleFont(22);
   beamSpec->GetYaxis()->SetLabelFont(22);
  
   beamSpecLast->GetXaxis()->SetTitleFont(22);
   beamSpecLast->GetXaxis()->SetLabelFont(22);
   beamSpecLast->GetYaxis()->SetTitleFont(22);
   beamSpecLast->GetYaxis()->SetLabelFont(22);

   gStyle->SetOptStat(0);
   
   c1->cd(1);
   normalizeFrame->Draw("colz");
   c1->cd(2);
   normalizeFrameLast->Draw("colz");

   c2->cd(1);
   beamSpec->Draw();

   c2->cd(2);
   beamSpecLast->Draw();
}

void drawClusterSizeDistribution(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   DataInterface *di = new DataInterface();

   Int_t nClusters = kEventsPerRun * 5 * nLayers;
   Int_t nHits = kEventsPerRun * 50;
   Int_t nLayersToUse = 6;

   CalorimeterFrame *cf = new CalorimeterFrame();
   Clusters * clusters = new Clusters(nClusters);
   Hits * hits = new Hits(nHits);
   TH2F * hClusterSizes = new TH2F("hClusterSizes", "Cluster sizes vs layer", nLayers, 0, nLayers-1, 50, 0, 50);
   
   vector<TH1F*> *hCSVector = new vector<TH1F*>;
   for (Int_t i=0; i<nLayersToUse; i++) {
      hCSVector->push_back(new TH1F(Form("Layer %d", i), Form(";Cluster Size;Entries"), 50, 0, 50));
   }

   for (Int_t i=0; i<Runs; i++) {
      if (dataType == kMC) {
         di->getMCFrame(i, cf);
         cf->diffuseFrame(new TRandom3(0));
         hits = cf->findHits();
         clusters = hits->findClustersFromHits();
      }

      else if (dataType == kData) {
         di->getDataFrame(i, cf, energy);
         hits = cf->findHits();
         printf("Found %d hits\n", hits->GetEntriesFast());
         clusters = hits->findClustersFromHits();
         printf("Found %d clusters\n", clusters->GetEntriesFast());
         clusters->removeSmallClusters(2);
         clusters->removeAllClustersAfterLayer(8);
      }
      
      clusters->Compress();

      for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
         hClusterSizes->Fill(clusters->getLayer(i), clusters->getSize(i));
         if (clusters->getLayer(i) < nLayersToUse) hCSVector->at(clusters->getLayer(i))->Fill(clusters->getSize(i));
      }
      cf->Reset();
   }

   TCanvas *c1 = new TCanvas("c1", "2D CS distribution", 1200, 800);
   hClusterSizes->Draw("COLZ");
   
   gStyle->SetPadRightMargin(0.03);
   gStyle->SetPadTopMargin(0.02);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.15);
   TCanvas *c2 = new TCanvas("c2", "1D CS distributions", 1200, 800);
   c2->Divide(3,2,0.0001,0.0001);

   Int_t colors[7] = {3, 2, 1, -4, -7, -9, -10};
   gStyle->SetOptStat(1111);
   for (Int_t i=0; i<nLayersToUse; i++) {
      c2->cd(i+1);
      hCSVector->at(i)->GetXaxis()->SetLabelSize(0.06);
      hCSVector->at(i)->GetXaxis()->SetTitleSize(0.06);
      hCSVector->at(i)->GetYaxis()->SetLabelSize(0.06);
      hCSVector->at(i)->GetYaxis()->SetTitleSize(0.06);
      hCSVector->at(i)->GetXaxis()->SetLabelFont(22);
      hCSVector->at(i)->GetXaxis()->SetTitleFont(22);
      hCSVector->at(i)->GetYaxis()->SetLabelFont(22);
      hCSVector->at(i)->GetYaxis()->SetTitleFont(22);
      hCSVector->at(i)->GetYaxis()->SetTitleOffset(1.3);

      hCSVector->at(i)->SetLineColor(kBlack);
      hCSVector->at(i)->SetFillColor(kRed + colors[i]);

      hCSVector->at(i)->Draw();
/*   
      TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
      title->SetTextFont(22);
      */
      gPad->Update();
      TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
      s->SetTextSize(0.05);
      s->SetTextFont(22);
      s->SetX1NDC(0.59);
      s->SetY1NDC(0.75);
      s->SetX2NDC(0.95);
      s->SetY2NDC(0.95);
      gPad->Modified();
   }

}

void compareChargeDiffusionModels(Int_t Runs, Bool_t recreate, Float_t energy) {
   run_energy = energy;

   DataInterface *di = new DataInterface();

   Int_t nClusters = kEventsPerRun * 5 * nLayers;
   Int_t nHits = kEventsPerRun * 50;
   Int_t nLayersToUse = 6;
   
   TCanvas *c = new TCanvas("c", "Compare charge diffusion models", 1200, 800);
   TCanvas *c4 = new TCanvas("c4", "Compare charge diffusion models", 1200, 800);
   
   TCanvas *c2 = new TCanvas("c2", "Compare charge diffusion models", 1200, 800);
   c2->Divide(3,2,1e-6,1e-6);

   TCanvas *c3 = new TCanvas("c3", "Compare charge diffusion models", 1200, 800);
   c3->Divide(2,1,1e-6,1e-6);

   int binx = 40;

   Float_t  edeps[1000];
   Float_t  css[1000];
   Float_t  edeps_error[1000];
   Float_t  css_error_high[1000];
   Float_t  css_error_low[1000];
   Int_t    graphIdx = 0;

   gStyle->SetTitleFont(22);
   gStyle->SetLabelFont(22);
   gStyle->SetTextFont(22);
   gStyle->SetLabelFont(22, "Y");
   gStyle->SetTitleFont(22, "Y");

   TH2F *hEdepVSCS = new TH2F("hEdepVSCS", ";E_{dep} [keV/#mum]; Calibrated cluster size", binx, 0.5, 5, 31, 0, 30);
   TH2F *hEdepVSCS1 = new TH2F("hEdepVSCS1", "122 MeV;E_{dep} [keV/#mum]; Calibrated cluster size", binx, 0.5, 5, 31, 0, 30);
   TH2F *hEdepVSCS2 = new TH2F("hEdepVSCS2", "140 MeV;E_{dep} [keV/#mum]; Calibrated cluster size", binx, 0.5, 5, 31, 0, 30);
   TH2F *hEdepVSCS3 = new TH2F("hEdepVSCS3", "150 MeV;E_{dep} [keV/#mum]; Calibrated cluster size", binx, 0.5, 5, 31, 0, 30);
   TH2F *hEdepVSCS4 = new TH2F("hEdepVSCS4", "170 MeV;E_{dep} [keV/#mum]; Calibrated cluster size", binx, 0.5, 5, 31, 0, 30);
   TH2F *hEdepVSCS5 = new TH2F("hEdepVSCS5", "180 MeV;E_{dep} [keV/#mum]; Calibrated cluster size", binx, 0.5, 5, 31, 0, 30);
   TH2F *hEdepVSCS6 = new TH2F("hEdepVSCS6", "188 MeV;E_{dep} [keV/#mum]; Calibrated cluster size", binx, 0.5, 5, 31, 0, 30);
   TH2F *hEnergyVSCSCal = new TH2F("hEnergyVSCSCal", "Calibrated;Residual energy [MeV];Cluster Size", 200, 0, 250, 31, 0, 30);
   TH2F *hEnergyVSCSUncal = new TH2F("hEnergyVSCSUncal", "Uncalibrated;Residual energy [MeV];Cluster Size", 200, 0, 250, 31, 0, 30);
   TF1 *fPheno2 = new TF1("fPheno2", "6.6177 * x ** 0.9383", 0, 50);
   TF1 *fPheno3 = new TF1("fPheno3", "7.77457 * x ** 0.7674", 0, 50);
   TF1 *fAnaly_30 = new TF1("fAnaly_30", "[0] + 6.16421 * x ** 3.81925e-1", 0, 50);
   TF1 *fAnaly_60 = new TF1("fAnaly_60", "[0] + 8.72351 * x ** 4.47309e-1", 0, 50);
   TF1 *fAnaly_65 = new TF1("fAnaly_65", "[0] + 8.99963 * x ** 4.54506e-1", 0, 50);
   TF1 *fAnaly_infty = new TF1("fAnaly_infty", "[0] + 8.99963 * x ** 4.54506e-1", 0, 50);
   TF1 *fAnaly_45 = new TF1("fAnaly_45", "[0] + 1.67988 * x ** 5.66446e-1", 0, 50);
   TF1 *fFinck = new TF1("fFinck", "2*3.1415 * [0] * log(14*x / (2*3.1415 * 0.0036 * [1]))", 0, 50);
   Cluster * cluster = nullptr;
   Float_t inst_wepl, inst_dedx, range, calcEnergy;
   Int_t inst_cs;

   fAnaly_30->SetParameter(0, 0);
   fAnaly_45->SetParameter(0, 0);
   fAnaly_65->SetParameter(0, 0);
   fAnaly_60->SetParameter(0, 0);
   fAnaly_infty->SetParameter(0, 0);

   fAnaly_30->SetParLimits(0, 0, 0.001);
   fAnaly_45->SetParLimits(0, 0, 0.001);
   fAnaly_60->SetParLimits(0, 0, 0.001);
   fAnaly_65->SetParLimits(0, 0, 0.001);
   fAnaly_infty->SetParLimits(0, 0, 0.001);

   fFinck->SetParameters(1.202, 221.58);

   fPheno2->SetLineWidth(5);
   fPheno3->SetLineWidth(5);
   fAnaly_30->SetLineWidth(5);
   fAnaly_60->SetLineWidth(5);
   fAnaly_45->SetLineWidth(5);
   fAnaly_65->SetLineWidth(5);
   fAnaly_infty->SetLineWidth(5);
   fFinck->SetLineWidth(5);

   fFinck->SetLineColor(kBlack);
   fPheno2->SetLineColor(kBlue-4);
   fPheno3->SetLineColor(kBlack);
   fAnaly_30->SetLineColor(kBlack);
   fAnaly_65->SetLineColor(kBlack);
   fAnaly_infty->SetLineColor(kBlack);

   CalorimeterFrame *cf = new CalorimeterFrame();
   Clusters * clusters = new Clusters(nClusters);
   Hits * hits = new Hits(nHits);
   Tracks *tracks = nullptr;

   for (Int_t e=0; e<nEnergies; e++) {
      energy = energies[e];
      run_energy = energy;
      tracks = loadOrCreateTracks(recreate, Runs, 1, energy);

      for (Int_t k=0; k<tracks->GetEntriesFast(); k++) {
         if (!tracks->At(k)) continue;

         TGraphErrors * graph = (TGraphErrors*) tracks->At(k)->doTrackFit(false, kUseCSDA);
         if (!graph) continue;

         range = tracks->At(k)->getFitParameterRange();
         calcEnergy = getEnergyFromTL(range);

         for (Int_t j=0; j<tracks->At(k)->GetEntriesFast(); j++) {
            if (!tracks->At(k)->At(j)) continue;
            cluster = tracks->At(k)->At(j); 

            inst_wepl = getWEPLFromTL(cluster->getLayermm() + tracks->At(k)->getPreTL());
            inst_dedx = 43.95 * pow(getEnergyAtWEPL(energy, inst_wepl), -0.748);

            inst_cs =   cluster->getCalibratedSize();

            if (inst_cs < 0) continue;
            if (inst_dedx > 3) {
            }
               hEnergyVSCSCal->Fill(getEnergyAtWEPL(energy, inst_wepl), inst_cs);
               hEnergyVSCSUncal->Fill(getEnergyAtWEPL(energy, inst_wepl), cluster->getSize());

            hEdepVSCS->Fill(inst_dedx, inst_cs);

            if       (energy == 122) hEdepVSCS1->Fill(inst_dedx, inst_cs);
            else if  (energy == 140) hEdepVSCS2->Fill(inst_dedx, inst_cs);
            else if  (energy == 150) hEdepVSCS3->Fill(inst_dedx, inst_cs);
            else if  (energy == 170) hEdepVSCS4->Fill(inst_dedx, inst_cs);
            else if  (energy == 180) hEdepVSCS5->Fill(inst_dedx, inst_cs);
            else if  (energy == 188) hEdepVSCS6->Fill(inst_dedx, inst_cs);
         }
      }
   }

   hEdepVSCS->GetXaxis()->SetTitleFont(22);
   hEdepVSCS->GetXaxis()->SetLabelFont(22);
   hEdepVSCS->GetXaxis()->SetTitleSize(0.05);
   hEdepVSCS->GetXaxis()->SetLabelSize(0.05);
   hEdepVSCS->GetYaxis()->SetTitleFont(22);
   hEdepVSCS->GetYaxis()->SetLabelFont(22);
   hEdepVSCS->GetYaxis()->SetTitleSize(0.05);
   hEdepVSCS->GetYaxis()->SetLabelSize(0.05);
   hEdepVSCS->GetZaxis()->SetLabelFont(22);
   hEdepVSCS->GetZaxis()->SetLabelSize(0.05);

   gStyle->SetOptStat(0);

   c->cd();
   gPad->SetBottomMargin(0.12);
   gPad->SetRightMargin(0.2);
   gPad->SetTopMargin(0.05);
   gPad->SetLogz();
   hEdepVSCS->GetYaxis()->SetTitleOffset(0.8);
   hEdepVSCS->Draw("COL");
   gPad->Update();

   TLatex *t1 = new TLatex(5, 25, "Power fit");
   TLatex *t4 = new TLatex(5, 20, "Gaussian fit");
   TLatex *t2 = new TLatex(5, 15, "#splitline{Analytical}{#lambda = 30 #mum}");
   TLatex *t3 = new TLatex(5, 10, "Analytical model");
   t1->SetTextFont(22);
   t2->SetTextFont(22);
   t3->SetTextFont(22);
   t4->SetTextFont(22);
   t1->SetTextSize(0.05);
   t2->SetTextSize(0.05);
   t3->SetTextSize(0.05);
   t4->SetTextSize(0.05);
   t1->Draw();
//   t2->Draw();
   t3->Draw();
   t4->Draw();

   c4->cd();
   TF1 * fitFunc = new TF1("fitFunc", "[0] * x ** [1]", 0, 5);
   fitFunc->SetParameters(6.6177, 0.9383); // pheno2
   
   // CALCULATE STATISTICS FROM hEdepVSCS:
   Float_t  sum_mean;
   Float_t  sum_variance;
   Float_t  sum_variance_high;
   Float_t  sum_variance_low;
   Float_t  mean;
   Float_t  standarddeviation;
   Float_t  standarddeviation_high;
   Float_t  standarddeviation_low;
   Float_t  fScale_high, fScale_low;
   Int_t    nEntries = 0;
   Int_t    nEntries_high = 0;
   Int_t    nEntries_low = 0;
   TAxis *xaxis = hEdepVSCS->GetXaxis();
   TAxis *yaxis = hEdepVSCS->GetYaxis();

   Int_t    ntotal = hEdepVSCS->Integral();
   Int_t    nEntries0 = 0;

   for (Int_t horiz=0; horiz<hEdepVSCS->GetNbinsX(); horiz++) {
      sum_mean = 0;
      sum_variance = 0;
      sum_variance_high = 0;
      sum_variance_low = 0;
      nEntries = 0;
      nEntries_high = 0;
      nEntries_low = 0;

      // calculate mean
      for (Int_t vert=0; vert<hEdepVSCS->GetNbinsY(); vert++) {
         sum_mean += yaxis->GetBinCenter(vert) * hEdepVSCS->GetBinContent(horiz, vert);
         nEntries += hEdepVSCS->GetBinContent(horiz, vert);
      }
      
      if (nEntries > 0) {
         if (!nEntries0) nEntries0 = nEntries;
         mean = sum_mean / nEntries;

         printf("Mean value at %.1f kev/um is %.1f\n", xaxis->GetBinCenter(horiz), mean);
      
         for (Int_t vert=0; vert<hEdepVSCS->GetNbinsY(); vert++) {
            if (yaxis->GetBinCenter(vert) < mean) {
               sum_variance_low += hEdepVSCS->GetBinContent(horiz, vert) * pow( yaxis->GetBinCenter(vert) - mean, 2);
               nEntries_low += hEdepVSCS->GetBinContent(horiz,vert);
            }
            else {
               sum_variance_high += hEdepVSCS->GetBinContent(horiz, vert) * pow( yaxis->GetBinCenter(vert) - mean, 2);
               nEntries_high += hEdepVSCS->GetBinContent(horiz,vert);
            }
         }

         // rescale  nentries to the total number of entries, so that the nEntries number also acts as a point WEIGHT
         // so that nEntries_high + nEntries_low = ntotal

         if (!nEntries_high || !nEntries_low) continue;
         
         fScale_high = (nEntries0 / nEntries_high) * 2;
         fScale_low = (nEntries0 / nEntries_low) * 2;

         printf("At %.1f keV/um, fScale l/h is %.2f / %.2f.\n", xaxis->GetBinCenter(horiz), fScale_low, fScale_high);

         if (nEntries_high == 1) nEntries_high = 2;
         if (nEntries_low == 1) nEntries_low = 2;

         standarddeviation_high = sqrt(sum_variance_high / (nEntries_low - 1));
         standarddeviation_low = sqrt(sum_variance_low / (nEntries_high - 1));
         
         css[graphIdx] = mean;
         css_error_high[graphIdx] = standarddeviation_high;
         css_error_low[graphIdx] = standarddeviation_low;
         edeps[graphIdx] = xaxis->GetBinCenter(horiz);
         edeps_error[graphIdx++] = 0;
      }
   }

   TGraphAsymmErrors *gEdepVSCS = new TGraphAsymmErrors(graphIdx, edeps, css, edeps_error, edeps_error, css_error_low, css_error_high);
   gEdepVSCS->SetMarkerStyle(21);
   gEdepVSCS->Draw("AP");
   fitFunc->SetParameters(7.85, 0.7265);
   
   fitFunc->Draw("same");

   printf("Fit function: %.3f * edep ^ %.3f\n", fitFunc->GetParameter(0), fitFunc->GetParameter(1));
   printf("Finck function: rs = %.3f, Ts = %.3f\n", fFinck->GetParameter(0), fFinck->GetParameter(1));

   fFinck->Draw("same");
   fAnaly_infty->Draw("same");
   
   TLatex *tt1 = new TLatex(5, 25, "Power Fit");
   TLatex *tt4 = new TLatex(5, 20, "Gaussian Intensity");
   TLatex *tt3 = new TLatex(5, 10, "First Principles");
   tt1->SetTextFont(22);
   tt3->SetTextFont(22);
   tt4->SetTextFont(22);
   tt1->SetTextSize(0.05);
   tt3->SetTextSize(0.05);
   tt4->SetTextSize(0.05);
   tt1->Draw();
   tt3->Draw();
   tt4->Draw();


   c2->cd(1);
   gPad->SetLogz();
   hEdepVSCS1->Draw("COLZ");
   c2->cd(2);
   gPad->SetLogz();
   hEdepVSCS2->Draw("COLZ");
   c2->cd(3);
   gPad->SetLogz();
   hEdepVSCS3->Draw("COLZ");
   c2->cd(4);
   gPad->SetLogz();
   hEdepVSCS4->Draw("COLZ");
   c2->cd(5);
   gPad->SetLogz();
   hEdepVSCS5->Draw("COLZ");
   c2->cd(6);
   gPad->SetLogz();
   hEdepVSCS6->Draw("COLZ");
   c->cd();
   fPheno3->Draw("same");
   fAnaly_infty->Draw("same");
   fFinck->Draw("same");

   TLegend *leg = new TLegend(0.55, 0.17, 0.85, 0.34);
   leg->AddEntry(fPheno3, "Phenomenological", "L");
   leg->AddEntry(fAnaly_infty, "Analytic #lambda = #infty #mum", "L");
   leg->SetTextFont(22);


   c3->cd(1);
   hEnergyVSCSCal->Draw("COLZ");
   c3->cd(2);
   hEnergyVSCSUncal->Draw("COLZ");

}

void secondaryAnalysis(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t switchLayer, Float_t energy, Float_t degraderThickness, Int_t tracksperrun, Bool_t doTracking) {
   run_energy = energy;
   run_degraderThickness = degraderThickness;
   kEventsPerRun = tracksperrun;
   kDoTracking = doTracking;

   Track  * thisTrack = nullptr;
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   Int_t    nTracks = tracks->GetEntries();
   Int_t    nTracksNuclear =0;
   Int_t    nTracksNuclearAfterAngle = 0;
   Int_t    nTracksNuclearAfterAngleAndEdep = 0;
   Int_t    nTracksNuclearAfterAngleAndEdepAndRange = 0;
   Float_t  angle;

   Cluster *a;
   Cluster *b;

   tracks->doTrackFit();

   TH1F   * hAngleAnalysisP = new TH1F("hAngleAnalysisP", "Angular distribution of incoming tracks (primary);Incoming angle [mrad];Frequency", 100, 0, 500);
   TH1F   * hAngleAnalysisS = new TH1F("hAngleAnalysisS", "Angular distribution of incoming tracks (scndary);Incoming angle [mrad];Frequency", 100, 0, 500);
   TH1F   * hAngleAnalysisAfter = new TH1F("hAngleAnalysisAfter", "Angular distribution of incoming tracks (After);Incoming angle [mrad];Frequency", 100, 0, 500);
   TH1F   * hRangeAnalysisP = new TH1F("hRangeAnalysisP", "Range distribution of tracks (primary);Range [mm WEPL];Frequency", 100, 0, 300);
   TH1F   * hRangeAnalysisS = new TH1F("hRangeAnalysisS", "Range distribution of tracks (scndary);Range [mm WEPL];Frequency", 100, 0, 300);
   TH1F   * hRangeAnalysisAfter = new TH1F("hRangeAnalysisAfter", "Range distribution of tracks (after);Range [mm WEPL];Frequency", 100, 0, 300);

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) {
         printf("!thisTrack\n");
         continue;
      }

      a = thisTrack->At(0);
      b = thisTrack->At(1);
      if (!a || !b) continue;

      angle = getDotProductAngle(a,a,b);

      if (thisTrack->Last()->isSecondary()) {
         hRangeAnalysisS->Fill(thisTrack->getFitParameterRange());
         hAngleAnalysisS->Fill(angle*1000);
      }
      else {
         hRangeAnalysisP->Fill(thisTrack->getFitParameterRange());
         hAngleAnalysisP->Fill(angle*1000);
      }
   }

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) {
         printf("!thisTrack\n");
         continue;
      }

      if (thisTrack->Last()->isSecondary()) nTracksNuclear++;
   }

   printf("There are %d tracks in total. %d of them are secondary tracks (%.2f%%).\n", nTracks, nTracksNuclear, 100*float(nTracksNuclear)/nTracks);

   tracks->removeHighAngleTracks(20); // 75 in protons
   nTracks = tracks->GetEntries();
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (thisTrack->Last()->isSecondary()) nTracksNuclearAfterAngle++;
   }

   printf("After ANGLE cut: There are %d tracks in total. %d of them are secondary tracks (%.2f%%).\n", nTracks, nTracksNuclearAfterAngle, 100*float(nTracksNuclearAfterAngle)/nTracks);
   
   tracks->removeNuclearInteractions();
   nTracks = tracks->GetEntries();
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (thisTrack->Last()->isSecondary()) nTracksNuclearAfterAngleAndEdep++;
   }

   printf("After EDEP cut: There are %d tracks in total. %d of them are secondary tracks (%.2f%%).\n", nTracks, nTracksNuclearAfterAngleAndEdep, 100*float(nTracksNuclearAfterAngleAndEdep)/nTracks);

   tracks->removeThreeSigmaShortTracks();
   nTracks = tracks->GetEntries();
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (thisTrack->Last()->isSecondary()) nTracksNuclearAfterAngleAndEdepAndRange++;

      a = thisTrack->At(0);
      b = thisTrack->At(1);
      if (!a || !b) continue;
      angle = getDotProductAngle(a,a,b);
      hRangeAnalysisAfter->Fill(thisTrack->getFitParameterRange());
      hAngleAnalysisAfter->Fill(angle*1000);
   }

   

   printf("After RANGE cut: There are %d tracks in total. %d of them are secondary tracks (%.2f%%).\n", nTracks, nTracksNuclearAfterAngleAndEdepAndRange, 100*float(nTracksNuclearAfterAngleAndEdepAndRange)/nTracks);

   TCanvas *c1 = new TCanvas();
   hAngleAnalysisP->SetLineColor(kRed);
   hAngleAnalysisS->SetLineColor(kGreen-3);
   hAngleAnalysisAfter->SetLineColor(kBlack);
   hAngleAnalysisP->SetLineWidth(3);
   hAngleAnalysisS->SetLineWidth(3);
   hAngleAnalysisAfter->SetLineWidth(3);
   hAngleAnalysisP->Draw();
   hAngleAnalysisS->Draw("same");
   hAngleAnalysisAfter->Draw("same");

   TCanvas *c2 = new TCanvas();
   hRangeAnalysisP->SetLineColor(kRed);
   hRangeAnalysisS->SetLineColor(kGreen-3);
   hRangeAnalysisAfter->SetLineColor(kBlack);
   hRangeAnalysisP->SetLineWidth(3);
   hRangeAnalysisS->SetLineWidth(3);
   hRangeAnalysisAfter->SetLineWidth(3);
   hRangeAnalysisP->Draw();
   hRangeAnalysisS->Draw("same");
   hRangeAnalysisAfter->Draw("same");

}


void drawEdep(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t switchLayer, Float_t energy, Float_t degraderThickness, Int_t tracksperrun, Bool_t doTracking) {
   run_energy = energy;
   run_degraderThickness = degraderThickness;
   kEventsPerRun = tracksperrun;
   kDoTracking = doTracking;
   
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   tracks->removeEmptyTracks();

   Track *thisTrack = nullptr;

   TH1F * edepP = new TH1F("edepP", "primary;edep;freq;", 400, 0, 200);
   TH1F * edepS = new TH1F("edepS", "secondary;edep;freq;", 400, 0, 200);

   for (int i=0; i<=tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;

      if (thisTrack->Last()->isSecondary()) {
         edepP->Fill(thisTrack->Last()->getDepositedEnergy());
      }
      else {
         edepS->Fill(thisTrack->Last()->getDepositedEnergy());
      }

   }

   edepP->Draw("COLZ");
   edepS->Draw("COLZ same");

}

void analyseSecondaryFilters(Int_t Runs, Int_t tracksPerRun, Float_t degraderThickness, Bool_t doTracking) {
   Float_t energy = 917;
   Float_t recreate = 1;
   Float_t dataType = 0;

   run_energy = energy;
   run_degraderThickness = degraderThickness;
   kEventsPerRun = tracksPerRun;
   kDoTracking = doTracking;

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   tracks->removeHighAngleTracks(30);
   tracks->removeNuclearInteractions();
   tracks->fillOutIncompleteTracks();
   tracks->doTrackFit();
   tracks->removeThreeSigmaShortTracks();
   tracks->removeEmptyTracks();
   
   Track        * thisTrack = nullptr;
   Cluster      * thisCluster = nullptr;
   TGraphErrors * outputGraph = nullptr;
   TGraph       * outputGraphTest = nullptr;
   Int_t          nGraph = 1;
   TCanvas      * c = new TCanvas("analyseSecondaryFilters", "analyseSecondaryFilters", 800,600);
   c->Divide(6,2);

   Int_t used = 0;
   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;
      if (nGraph>=6) break;
      if (j<10) continue;

      // Do track fit, extract all parameters for this track
      outputGraph = (TGraphErrors*) thisTrack->doTrackFit(false, kUseCSDA); // (bool isScaleVariable, bool useTrackLength (~ CSDA))
      if (!outputGraph) continue;

      c->cd(nGraph);
      outputGraph->SetTitle(Form("%d isSecondary: %d", j, thisTrack->Last()->isSecondary()));
      outputGraph->Draw("AP");

      c->cd(nGraph+6);
      Float_t x[50] = {};
      Float_t y[50] = {};
      Int_t idx = 0, n = thisTrack->GetEntriesFast();
      Float_t lastTwo, nextTwo;
      for (Int_t i=1; i<n-1; i++) {
         if (!thisTrack->At(i)) continue;

         lastTwo = thisTrack->getAverageDepositedEnergy(i-2,i);
         nextTwo = thisTrack->getAverageDepositedEnergy(i,i+2);

         x[idx] = thisTrack->At(i)->getLayer();
         y[idx++] = lastTwo - nextTwo;
      }
      outputGraphTest = new TGraph(idx,x,y);
      outputGraphTest->SetTitle(Form("%d;Layer;lastTwo-nextTwo",j));
      outputGraphTest->SetMarkerStyle(7);
      outputGraphTest->Draw("AP");

      nGraph++;
   }
}

void analyseHelium(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t switchLayer, Float_t energy, Float_t degraderThickness, Int_t tracksperrun, Bool_t doTracking) {
   run_energy = energy;
   run_degraderThickness = degraderThickness;
   kEventsPerRun = tracksperrun;
   kDoTracking = doTracking;
   
   Clusters * savedClusters = new Clusters();

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy,0,0, savedClusters);
   tracks->removeEmptyTracks();

   Int_t nSecondariesEnd = 0, nPrimariesEnd = 0, nTotal = 0;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      nSecondariesEnd += tracks->At(i)->Last()->isSecondary();
      nPrimariesEnd += (!tracks->At(i)->Last()->isSecondary());
   }
   nTotal = tracks->GetEntries();

   // From a no-halo filter run
   nTotal = 126881;
   nPrimariesEnd = 50265;
   nSecondariesEnd = 76616;

   Float_t  cutHalo = 12;
   Float_t  cutMaxAngle = 60;
   Float_t  cutAngle = 40;
   Float_t  cutEdep = 10;

   printf("Found in total %d tracks. %d primaries, %d secondaries\n", nTotal, nPrimariesEnd, nSecondariesEnd);


// Don't use
//   track->removeHighChiSquare(210);

// Use
//   tracks->removeHighAngleTracks(cutAngle); // mrad
//   tracks->removeHighAngularChangeTracks(cutMaxAngle); // mrad
//   tracks->doTrackFit();
//   tracks->removeNuclearInteractions();
//   tracks->removeThreeSigmaShortTracks();
   
   Int_t nfSecondariesEnd = 0, nfPrimariesEnd = 0, nfTotal = 0;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      nfSecondariesEnd += tracks->At(i)->Last()->isSecondary();
      nfPrimariesEnd += (!tracks->At(i)->Last()->isSecondary());
   }
   nfTotal = tracks->GetEntries();

   printf("After filtering: %d tracks. %d primaries, %d secondaries\n", nfTotal, nfPrimariesEnd, nfSecondariesEnd);

   printf("True positives (filtered secondaries): %.1f %%\n", float(nSecondariesEnd-nfSecondariesEnd)/nTotal*100);
   printf("False positives (filtered primaries): %.1f %%\n", float(nPrimariesEnd-nfPrimariesEnd)/nTotal*100);
   printf("True negatives (unfiltered primaries): %.1f %%\n", (nfPrimariesEnd)/float(nTotal) * 100);
   printf("False negatives (unfiltered secondaries): %.1f %%\n", (nfSecondariesEnd)/float(nTotal) * 100);

   /*
   tracks->removeEmptyTracks();
   tracks->fillOutIncompleteTracks(0.05);
   */

   Track *thisTrack = nullptr;
   Cluster *a, *b;
   Float_t angle, range, edep, bragg;

   TCanvas *cEdep = new TCanvas("cEdep","cEdep",1200,600);
   cEdep->Divide(2,1);
   TH1F * edepP = new TH1F("edepP", "Primary at end;E_{dep} last layer [kev/#mum];Entries", 100, 0, 100);
   TH1F * edepS = new TH1F("edepS", "Secondary at end;E_{dep} last layer [kev/#mum];Entries", 100, 0, 100);

   TCanvas *cEdep2D = new TCanvas("cEdep2D","cEdep2D",1200,600);
   cEdep2D->Divide(2,1);
   TH2F *edep2DP = new TH2F("edep2DP", "Primary at end;Edep last layer [kev/#mum];Edep next-to-last layer", 60, 0, 60, 60, 0, 1000);
   TH2F *edep2DS = new TH2F("edep2DS", "Secondary at end;#Delta E detector mid [kev/#mum]; residual energy after mid [MeV]", 60, 0, 60,60, 0, 1000);

   TCanvas *cAllEdep = new TCanvas("cAllEdep","cAllEdep",1200,600);
   cAllEdep->Divide(2,1);
   TH1F * cAllEdepP = new TH1F("allEdepP", "Primary at end;E_{dep} all layers [keV/#mum];Entries",100,0,100);
   TH1F * cAllEdepS = new TH1F("allEdepS", "Secondary at end;E_{dep} all layers [kev/#mum];Entries",100,0,100);
   
   TCanvas *cAngle = new TCanvas("cAngle","cAngle",1200,600);
   cAngle->Divide(2,1);
   TH1F * angleP = new TH1F("angleP", "Primary at end;Incoming angle [mrad];Entries;", 100, 0, 200);
   TH1F * angleS = new TH1F("angleS", "Secondary at end;Incoming angle [mrad];Entries;", 100, 0, 200);
   
   TCanvas *cRange = new TCanvas("cRange","cRange",1200,600);
   cRange->Divide(2,1);
   TH1F * rangeP = new TH1F("rangeP", "Primary at end;Residual range [WEPL mm];Entries;", 100, 0, 300);
   TH1F * rangeS = new TH1F("rangeS", "Secondary at end;Residual range [WEPL mm];Entries;", 100, 0, 300);

   TCanvas *cRange2 = new TCanvas("cRange2", "cRange", 800, 600);
   TH1F * rangePS = new TH1F("rangePS", "All;Residual range [WEPL mm];Entries;", 100, 0, 300);

   TCanvas *cBragg = new TCanvas("cBragg","cBragg",1200,600);
   cBragg->Divide(2,1);
   TH1F * braggP = new TH1F("rangeP", "Primary at end;Depth-dose fit log_{10} #chi^{2};Entries;", 100, 1, 2000);
   TH1F * braggS = new TH1F("rangeS", "Secondary at end;Depth-dose fit log_{10} #chi^{2};Entries;", 100, 0, 2000);
   
   TCanvas *cSecondary = new TCanvas();
   TH1F *hSec = new TH1F("hSec", "Secondary prouction depth;Detector layer;Entries", 46, -0.5, 45.5);
   TH1F *hPrim = new TH1F("hPrim", "Primary beam;layer;freq", 50, 0, 50);

   TCanvas *cMaxDeltaTheta = new TCanvas("cMaxDeltaTheta","cMaxDeltaTheta", 1200,600);
   cMaxDeltaTheta->Divide(2,1);
   TH1F *hMaxDeltaThetaP = new TH1F("hMaxDeltaThetaP", "Primary at end;Max layer-wise angular change [mrad];Entries",100,0,300);
   TH1F *hMaxDeltaThetaS = new TH1F("hMaxDeltaThetaS", "Secondary at end;Max layer-wise angular change [mrad];Entries",100,0,300);

   TCanvas *cHits = new TCanvas("cHits", "cHits", 1200,600);
   cHits->Divide(2,1);
   TH1F *cHitsP = new TH1F("cHitsP", "Primary at end;Beam profile at entrance layer #sqrt{x^{2}+y^{2}} [mm];Entries", 100, 0, 100);
   TH1F *cHitsS = new TH1F("cHitsS", "Secondary at end;Beam profile at entrance layer #sqrt{x^{2}+y^{2}} [mm];Entries", 100, 0, 100);

   for (int i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
   
      bool was = false;
      float biggestChange = 0;
      for (int j=0; j<thisTrack->GetEntriesFast(); j++) {
         if (!was) {
            if (thisTrack->At(j)->isSecondary()) {
               hSec->Fill(thisTrack->getLayer(j));
               was = true;
            }
         }
      }
      int last = thisTrack->GetEntriesFast() - 1;
      if (last>=3) {
         float edeplast = thisTrack->Last()->getDepositedEnergy();
         float edepnlast = thisTrack->At(last-1)->getDepositedEnergy();
         float edepnnlast = thisTrack->At(last-2)->getDepositedEnergy();
         float edepnnnlast = thisTrack->At(last-3)->getDepositedEnergy();
         if (!thisTrack->Last()->isSecondary()) {
               cAllEdepP->Fill(edeplast+edepnlast);
         }
         else {
               cAllEdepS->Fill(edeplast+edepnlast);
         }
      }
      
      Cluster *a = thisTrack->At(0);
      Cluster *b = thisTrack->At(1);
      if (!a || !b) continue;
      angle = getDotProductAngle(a, a, b) * 1000;
      if (isnan(thisTrack->getFitParameterRange())) continue;

//      range = getUnitFromTL(thisTrack->getFitParameterRange());
      range = getUnitFromTL(thisTrack->getRangemm());
      bragg = thisTrack->getFitParameterChiSquare();
      edep = 0;
      last = thisTrack->GetEntriesFast() - 1;
      if (thisTrack->GetEntriesFast() > 3) {
//         edep = thisTrack->getDepositedEnergy(last) + thisTrack->getDepositedEnergy(last-1);
         edep = thisTrack->getDepositedEnergy(last);
      }

      rangePS->Fill(range);
      if (!thisTrack->Last()->isSecondary()) {
         edepP->Fill(edep);
         angleP->Fill(angle);
         rangeP->Fill(range);
         braggP->Fill(bragg);

         if (last>5) { // thisTrack->At(last-5) && thisTrack->At(last-4)) {
            float edep_mid = thisTrack->getDepositedEnergy(last-5);
            float residual_energy = getEnergyAtWEPL(917,thisTrack->getRangeWEPLAt(last-4));
            edep2DP->Fill(edep_mid, residual_energy);
         }


         hMaxDeltaThetaP->Fill(thisTrack->getMaximumSlopeAngleChange() * 3.1415 / 180 * 1000);
      }
      else {
         edepS->Fill(edep);
         angleS->Fill(angle);
         rangeS->Fill(range);
         braggS->Fill(bragg);
         hMaxDeltaThetaS->Fill(thisTrack->getMaximumSlopeAngleChange() * 3.1415 / 180 * 1000);
         if (last>5) {
            float edep_mid = thisTrack->getDepositedEnergy(last-5);
            float residual_energy = getEnergyAtWEPL(917,thisTrack->getRangeWEPLAt(last-4));
            edep2DS->Fill(edep_mid, residual_energy);
         }
      }
   }

   float x,y,z;
   Cluster * cl;
   int sgn;
   for (Int_t i=0; i<=savedClusters->GetEntriesFast(); i++) {
      sgn = 1;
      cl = savedClusters->At(i);
      if (!cl) continue;
      x = cl->getXmm();
      y = cl->getYmm();
      z = cl->getLayer();
      
      if (y<0) sgn = -1;

      if (cl->isSecondary()) {
        cHitsS->Fill(sgn*sqrt(x*x+y*y));
      }
      else {
         cHitsP->Fill(sgn*sqrt(x*x+y*y));
      }
   }
   delete savedClusters;

   rangePS->SetFillColor(kRed);
   edepP->SetFillColor(kGray);
   edepS->SetFillColor(kGreen-3);
   angleP->SetFillColor(kGray);
   angleS->SetFillColor(kGreen-3);
   rangeP->SetFillColor(kGray);
   rangeS->SetFillColor(kGreen-3);
   braggP->SetFillColor(kGray);
   braggS->SetFillColor(kGreen-3);
   cAllEdepP->SetFillColor(kGray);
   cAllEdepS->SetFillColor(kGreen-3);
   hSec->SetFillColor(kGreen-3);
   hMaxDeltaThetaP->SetFillColor(kGray);
   hMaxDeltaThetaS->SetFillColor(kGreen-3);
   cHitsP->SetFillColor(kGray);
   cHitsS->SetFillColor(kGreen-3);

   // Please tell me how to do this properly with gStyle
   rangePS->SetLineColor(kBlack);
   edepP->SetLineColor(kBlack);
   edepS->SetLineColor(kBlack);
   angleP->SetLineColor(kBlack);
   angleS->SetLineColor(kBlack);
   rangeP->SetLineColor(kBlack);
   rangeS->SetLineColor(kBlack);
   braggP->SetLineColor(kBlack);
   braggS->SetLineColor(kBlack);
   cAllEdepP->SetLineColor(kBlack);
   cAllEdepS->SetLineColor(kBlack);
   hSec->SetLineColor(kBlack);
   hMaxDeltaThetaP->SetLineColor(kBlack);
   hMaxDeltaThetaS->SetLineColor(kBlack);
   cHitsP->SetLineColor(kBlack);
   cHitsS->SetLineColor(kBlack);

   TLine *l1, *l2;

   cEdep->cd(1);
   edepP->Draw();
   l1 = new TLine(cutEdep, 0, cutEdep, edepP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cEdep->cd(2);
   edepS->Draw();
   l1 = new TLine(cutEdep, 0, cutEdep, edepS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   cAngle->cd(1);
   angleP->Draw();
   l1 = new TLine(cutAngle, 0, cutAngle, angleP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cAngle->cd(2);
   angleS->Draw();
   l1 = new TLine(cutAngle, 0, cutAngle, angleS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   Float_t mu = rangePS->GetXaxis()->GetBinCenter(rangePS->GetMaximumBin());
   Float_t sigma = 10; // as expected from ~homogenous area // WEPL 

   cRange->cd(1);
   rangeP->Draw();
   l1 = new TLine(mu-3*sigma, 0, mu-3*sigma, rangeP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   l2 = new TLine(mu+1.5*sigma, 0, mu+1.5*sigma, rangeP->GetMaximum()*1.05);
   l2->SetLineStyle(7);
   l2->SetLineColor(kBlack);
   l2->Draw();

   cRange->cd(2);
   rangeS->Draw();
   l1 = new TLine(mu-3*sigma, 0, mu-3*sigma, rangeS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   l2 = new TLine(mu+1.5*sigma, 0, mu+1.5*sigma, rangeS->GetMaximum()*1.05);
   l2->SetLineStyle(7);
   l2->SetLineColor(kBlack);
   l2->Draw();

   cRange2->cd();
   rangePS->Draw();
   l1 = new TLine(mu-3*sigma, 0, mu-3*sigma, rangePS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   l2 = new TLine(mu+1.5*sigma, 0, mu+1.5*sigma, rangePS->GetMaximum()*1.05);
   l2->SetLineStyle(7);
   l2->SetLineColor(kBlack);
   l2->Draw();

   cBragg->cd(1);
   braggP->Draw();
   l1 = new TLine(450, 0, 450, braggP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cBragg->cd(2);
   braggS->Draw();
   l1 = new TLine(450, 0, 450, braggS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   cAllEdep->cd(1);
   cAllEdepP->Draw();
   cAllEdep->cd(2);
   cAllEdepS->Draw();

   cEdep2D->cd(1);
   edep2DP->Draw("colz");
   cEdep2D->cd(2);
   edep2DS->Draw("colz");

   cHits->cd(1);
   cHitsP->Draw();
   l1 = new TLine(12, 0, 12, cHitsP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cHits->cd(2);
   cHitsS->Draw();
   l1 = new TLine(12, 0, 12, cHitsS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   // True positive: Secondaries removed by filter
   // True negative: Primaries not removed filter
   // False positive: Primaries removed by filter
   // False negative: Secondaries not removed by filter
   // Normalize by the total number of particles

   int bmin, bmax;
   int ntracks = tracks->GetEntries();

   bmin = cHitsP->GetXaxis()->FindBin(12);
   int totalHalo = cHitsS->Integral() + cHitsP->Integral();
   int removedByHalo = cHitsS->Integral(bmin,100) + cHitsP->Integral(bmin,100);
   int removedByHalo2nd = cHitsS->Integral(bmin,100);
   
   float haloTP = cHitsS->Integral(bmin,100) / float(totalHalo);
   float haloTN = (cHitsP->Integral() - cHitsP->Integral(bmin,100)) / float(totalHalo);
   float haloFP = cHitsP->Integral(bmin,100) / float(totalHalo);
   float haloFN = (cHitsS->Integral() - cHitsS->Integral(bmin,100)) / float(totalHalo);

   bmin = angleS->GetXaxis()->FindBin(40);
   int removedByAngle = angleS->Integral(bmin,100) + angleP->Integral(bmin,100);
   int removedByAngle2nd = angleS->Integral(bmin,100);
   
   float angleTP = angleS->Integral(bmin,100) / float(ntracks);
   float angleTN = (angleP->Integral() - angleP->Integral(bmin,100)) / float(ntracks);
   float angleFP = angleP->Integral(bmin,100) / float(ntracks);
   float angleFN = (angleS->Integral() - angleS->Integral(bmin,100)) / float(ntracks);

   bmin = hMaxDeltaThetaS->GetXaxis()->FindBin(50);
   int removedByMaxAngle = hMaxDeltaThetaS->Integral(bmin,100) + hMaxDeltaThetaP->Integral(bmin,100);
   int removedByMaxAngle2nd = hMaxDeltaThetaS->Integral(bmin,100);
   
   float maxAngleTP = hMaxDeltaThetaS->Integral(bmin,100) / float(ntracks);
   float maxAngleTN = (hMaxDeltaThetaP->Integral() - hMaxDeltaThetaP->Integral(bmin,100)) / float(ntracks);
   float maxAngleFP = hMaxDeltaThetaP->Integral(bmin,100) / float(ntracks);
   float maxAngleFN = (hMaxDeltaThetaS->Integral() - hMaxDeltaThetaS->Integral(bmin,100)) / float(ntracks);

   bmin = rangeS->GetXaxis()->FindBin(mu-3*sigma);
   bmax = rangeS->GetXaxis()->FindBin(mu+3*sigma);
   int removedByR = rangeS->Integral(0,bmin-1) + rangeS->Integral(bmax+1,100) + rangeP->Integral(0,bmin-1) + rangeP->Integral(bmax+1,100);
   int removedByR2nd = rangeS->Integral(0,bmin-1) + rangeS->Integral(bmax+1,100);
   
   float rangeTP = (rangeS->Integral() - rangeS->Integral(bmin-1,bmax+1)) / float(ntracks);
   float rangeTN = (rangeP->Integral(bmin,bmax)) / float(ntracks);
   float rangeFP = (rangeP->Integral() - rangeP->Integral(bmin-1,bmax+1)) / float(ntracks);
   float rangeFN = (rangeS->Integral(bmin,bmax)) / float(ntracks);

   bmin = edepP->GetXaxis()->FindBin(10);
   int removedByED = edepS->Integral(0,bmin) + edepP->Integral(0,bmin);
   int removedByED2nd = edepS->Integral(0,bmin);
   
   float edepTP = edepS->Integral(bmin,100) / float(ntracks);
   float edepTN = (edepP->Integral() - edepP->Integral(bmin,100)) / float(ntracks);
   float edepFP = edepP->Integral(bmin,100) / float(ntracks);
   float edepFN = (edepS->Integral() - edepS->Integral(bmin,100)) / float(ntracks);

   bmin = braggP->GetXaxis()->FindBin(450);
   int removedByChi2 = braggS->Integral(bmin,100) + braggP->Integral(bmin,100);
   int removedByChi22nd = braggS->Integral(bmin,100);
   
   float chi2TP = braggS->Integral(bmin,100) / float(ntracks);
   float chi2TN = (braggP->Integral() - braggP->Integral(bmin,100)) / float(ntracks);
   float chi2FP = braggP->Integral(bmin,100) / float(ntracks);
   float chi2FN = (braggS->Integral() - braggS->Integral(bmin,100)) / float(ntracks);

   bmin = cAllEdepP->GetXaxis()->FindBin(20);
   int removedByAllEdep = cAllEdepS->Integral(0,bmin) + cAllEdepP->Integral(bmin,100);
   int removedByAllEdep2nd = cAllEdepS->Integral(0,bmin);
   
   float allEdepTP = cAllEdepS->Integral(bmin,100) / float(ntracks);
   float allEdepTN = (cAllEdepP->Integral() - cAllEdepP->Integral(bmin,100)) / float(ntracks);
   float allEdepFP = cAllEdepP->Integral(bmin,100) / float(ntracks);
   float allEdepFN = (cAllEdepS->Integral() - cAllEdepS->Integral(bmin,100)) / float(ntracks);
      
   printf("Total number of tracks = %d\n", ntracks);
   printf("Tracks removed by Halo: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByHalo, removedByHalo2nd, 100.*removedByHalo2nd/removedByHalo, 100.*removedByHalo/totalHalo);
   printf("Tracks removed by Angle: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByAngle, removedByAngle2nd, 100.*removedByAngle2nd/removedByAngle, 100.*removedByAngle/ntracks);
   printf("Tracks removed by MaxAngle: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByMaxAngle, removedByMaxAngle2nd, 100.*removedByMaxAngle2nd/removedByMaxAngle, 100.*removedByMaxAngle/ntracks);
   printf("Tracks removed by R: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByR, removedByR2nd, 100.*removedByR2nd/removedByR, 100.*removedByR/ntracks);
   printf("Tracks removed by last ED: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByED, removedByED2nd, 100.*removedByED2nd/removedByED, 100.*removedByED/ntracks);
   printf("Tracks removed by last2 ED: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByAllEdep, removedByAllEdep2nd, 100.*removedByAllEdep2nd/removedByAllEdep, 100.*removedByAllEdep/ntracks);
   printf("Tracks removed by Chi2: %d (%d secondaries) (%.1f%% purity) (%.1f%% removed)\n", removedByChi2, removedByChi22nd, 100.*removedByChi22nd/removedByChi2, 100.*removedByChi2/ntracks);

   printf("True/False Positive/Negative: \n\n");
   printf("Halo: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", haloTP*100, haloFP*100, haloTN*100, haloFN*100);
   printf("Angle: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", angleTP*100, angleFP*100, angleTN*100, angleFN*100);
   printf("MaxAngle: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", maxAngleTP*100, maxAngleFP*100, maxAngleTN*100, maxAngleFN*100);
   printf("Range: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", rangeTP*100, rangeFP*100, rangeTN*100, rangeFN*100);
   printf("Edep: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", edepTP*100, edepFP*100, edepTN*100, edepFN*100);
   printf("Chi2: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", chi2TP*100, chi2FP*100, chi2TN*100, chi2FN*100);
   printf("AllEdep: TP = %.1f, FP = %.1f, TN = %.1f, FN = %.1f\n", allEdepTP*100, allEdepFP*100, allEdepTN*100, allEdepFN*100);

   cSecondary->cd();
//   hPrim->Draw();
   hSec->Draw();

   cMaxDeltaTheta->cd(1);
   hMaxDeltaThetaP->Draw();
   l1 = new TLine(cutMaxAngle, 0, cutMaxAngle, hMaxDeltaThetaP->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();
   cMaxDeltaTheta->cd(2);
   hMaxDeltaThetaS->Draw();
   l1 = new TLine(cutMaxAngle, 0, cutMaxAngle, hMaxDeltaThetaS->GetMaximum()*1.05);
   l1->SetLineStyle(7);
   l1->SetLineColor(kBlack);
   l1->Draw();

   Int_t particleListIndex[100] = {};
   Int_t particleListSort[100] = {};
   Int_t particleList[100] = {};
   Int_t idxPDG = 0, pdg;
   Bool_t inList = false;

   // last particle
   Cluster * c = nullptr;
   for (Int_t i=0; i<=tracks->GetEntriesFast(); i++) {
      Track *thisTrack = tracks->At(i);
      if (!thisTrack) continue;

      c = thisTrack->Last();
      if (!c) continue;
      pdg = c->getPDG();
      if (!c->isSecondary()) continue;
      
      inList = false;
      for (int j=0; j<idxPDG; j++) {
         if (particleListIndex[j] == pdg) {
            inList = true;
            particleList[j]++;
            break;
         }
      }
      if (!inList) {
         particleListIndex[idxPDG] = pdg;
         particleList[idxPDG++] = 1;
      }
   }

   // Make indexed sorting list
   Int_t currentMax = 0;
   Int_t maxIdx = 0;
   Int_t sortIdx = 0;
   while (sortIdx < idxPDG) {
      currentMax = 0;
      for (Int_t i=0; i<idxPDG; i++) {
         if (particleList[i] >= currentMax) {
            inList = false;
            for (Int_t j=0; j<sortIdx; j++) inList += (particleListSort[j] == i);
            if (inList) continue;

            currentMax = particleList[i];
            maxIdx = i;
         }
      }
      particleListSort[sortIdx++] = maxIdx;
   }

   string name;
   Int_t p;

   printf("Secondary particle types: \n");
   for (Int_t i=0; i<idxPDG; i++) {
      p = particleListIndex[particleListSort[i]];
      if (p == 11) name = "Electron";
      else if (p == -11) name = "Positron";
      else if (p == 2212) name = "Proton";
      else if (p == 22) name = "Gamma";
      else if (p == 1000030040) name = "Litium4"; 
      else if (p == 1000030060) name = "Litium";
      else if (p == 1000020040) name = "Helium";
      else if (p == 1000020030) name = "Helium3";
      else if (p == 1000010030) name = "Tritium";
      else if (p == 1000010020) name = "Deuterium";
      else if (p == 1000130270) name = "Aluminium";
      else if (p == 1000140280) name = "Silicon";
      else if (p == 1000140260) name = "Silicon26";
      else if (p == 1000140299) name = "Silicon29";
      else if (p == 1000120240) name = "Magnesium";
      else name = Form("%d", p);

      cout << name << ": " << particleList[particleListSort[i]] << " tracks\n";
   }

   TCanvas *cPDG = new TCanvas("cPDG", "PDG(z)");
   cPDG->cd();
   float xfrom = -0.5;
   float xto = 50.5;
   int xlen = 51;
   
   gStyle->SetOptStat(0);
/*
   gStyle->SetOptStat(0);
   gStyle->SetTitleFont(22);
   gStyle->SetLabelFont(22);
   gStyle->SetTitleSize(0.05);
   gStyle->SetLabelSize(0.05);
   gStyle->SetTextFont(22);
   gStyle->SetLabelFont(22, "Y");
   gStyle->SetLabelSize(0.05, "Y");
   gStyle->SetTitleFont(22, "Y");
   gStyle->SetTitleSize(0.05, "Y");
   gStyle->SetTitleOffset(0.92);
   gStyle->SetTitleOffset(0.92, "Y");
*/

   TH1F *hPDGe = new TH1F("hPDGe", "Electrons", xlen, xfrom, xto);
   TH1F *hPDGpos = new TH1F("hPDGpos", "Positrons", xlen, xfrom, xto);
   TH1F *hPDGpro = new TH1F("hPDGpro", ";Depth in detector [layer number];Fraction of secondary particles", xlen, xfrom, xto);
   TH1F *hPDGneu = new TH1F("hPDGneu", ";Depth in detector [layer number];Fraction of secondary particles", xlen, xfrom, xto);
   TH1F *hPDGhe = new TH1F("hPDGhe", "Helium", xlen, xfrom, xto);
   TH1F *hPDGhe3 = new TH1F("hPDGhe3", ";Depth in detector [layer number];Fraction of secondary particles", xlen, xfrom, xto);
   TH1F *hPDGgam = new TH1F("hPDGgam", "Gamma", xlen, xfrom, xto);
   TH1F *hPDGli = new TH1F("hPDGli", "Litium", xlen, xfrom, xto);
   TH1F *hPDGdeu = new TH1F("hPDGdeu", "Deuterium", xlen, xfrom, xto);
   TH1F *hPDGtri = new TH1F("hPDGtri", "Tritium", xlen, xfrom, xto);
   TH1F *hPDGsi = new TH1F("hPDGsi", "Silicon", xlen, xfrom, xto);
   TH1F *hPDGmg = new TH1F("hPDGmg", "Magnesium", xlen, xfrom, xto);
   TH1F *hPDGal = new TH1F("hPDGal", "Aluminum", xlen, xfrom, xto);

   Int_t l;
   for (Int_t i=0; i<=tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      for (Int_t j=0; j<=thisTrack->GetEntriesFast(); j++) {
         if (!thisTrack->At(j)) continue;
         l = thisTrack->getLayer(j);
         p = thisTrack->At(j)->getPDG();
         if (!thisTrack->At(j)->isSecondary()) continue;

         if (p == 11) hPDGe->Fill(l);
         else if (p == -11) hPDGpos->Fill(l);
         else if (p == 2212) hPDGpro->Fill(l);
         else if (p == 2112) hPDGneu->Fill(l);
         else if (p == 22) hPDGgam->Fill(l);
         else if (p == 1000030040) hPDGli->Fill(l);
         else if (p == 1000030060) hPDGli->Fill(l);
         else if (p == 1000020040) hPDGhe->Fill(l);
         else if (p == 1000020030) hPDGhe3->Fill(l);
         else if (p == 1000010030) hPDGtri->Fill(l);
         else if (p == 1000010020) hPDGdeu->Fill(l);
         else if (p == 1000130270) hPDGal->Fill(l);
         else if (p == 1000140280) hPDGsi->Fill(l);
         else if (p == 1000140260) hPDGsi->Fill(l);
         else if (p == 1000140299) hPDGsi->Fill(l);
         else if (p == 1000120240) hPDGmg->Fill(l);
      }
   }

   hPDGhe->SetLineColor(8); hPDGpro->SetLineColor(1); hPDGdeu->SetLineColor(2); hPDGhe3->SetLineColor(3);
   hPDGtri->SetLineColor(4); hPDGe->SetLineColor(9); hPDGgam->SetLineColor(6); hPDGpos->SetLineColor(7);
   hPDGal->SetLineColor(14); hPDGsi->SetLineColor(16); hPDGmg->SetLineColor(18); hPDGli->SetLineColor(12);
   hPDGneu->SetLineColor(28); 

   hPDGhe->SetLineWidth(2); hPDGpro->SetLineWidth(2); hPDGdeu->SetLineWidth(2); hPDGhe3->SetLineWidth(2);
   hPDGtri->SetLineWidth(2); hPDGe->SetLineWidth(2); hPDGli->SetLineWidth(2); hPDGgam->SetLineWidth(2);
   hPDGal->SetLineWidth(2); hPDGsi->SetLineWidth(2); hPDGmg->SetLineWidth(2); hPDGpos->SetLineWidth(2);
   hPDGneu->SetLineWidth(2);

   // NORMALIZE
   Float_t totalParticles = tracks->GetEntries();

   hPDGhe->Scale(100 / totalParticles);
   hPDGpro->Scale(100 / totalParticles);
   hPDGdeu->Scale(100 / totalParticles);
   hPDGhe3->Scale(100 / totalParticles);
   hPDGtri->Scale(100 / totalParticles);
   hPDGe->Scale(100 / totalParticles);
   hPDGgam->Scale(100 / totalParticles);
   hPDGpos->Scale(100 / totalParticles);
   hPDGal->Scale(100 / totalParticles);
   hPDGsi->Scale(100 / totalParticles);
   hPDGmg->Scale(100 / totalParticles);
   hPDGli->Scale(100 / totalParticles);
   hPDGneu->Scale(100 / totalParticles);

   TLegend *leg = new TLegend(0.8,0.4,0.95,0.95);

   hPDGpro->Draw("hist");
   leg->AddEntry(hPDGpro, "Protons", "L");
   hPDGhe3->Draw("same hist"); 
   leg->AddEntry(hPDGhe3, "Helium3", "L");

   hPDGpro->GetYaxis()->SetRangeUser(0, 1.2*max(hPDGhe->GetMaximum(), max(hPDGpro->GetMaximum(), hPDGhe3->GetMaximum())));

   if (hPDGe->Integral()) {
      hPDGe->Draw("same hist"); 
      leg->AddEntry(hPDGe, "Electrons", "L");
   }

   if (hPDGneu->Integral()) {
      hPDGneu->Draw("same hist");
      leg->AddEntry(hPDGneu, "Neutrons", "L");
   }

   if (hPDGpos->Integral()) {
      hPDGpos->Draw("same hist");
      leg->AddEntry(hPDGpos, "Positrons", "L");
   }


   if (hPDGhe->Integral()) {
      hPDGhe->Draw("same hist"); 
      leg->AddEntry(hPDGhe, "Helium", "L"); 
   }


   if (hPDGgam->Integral()) { 
      hPDGgam->Draw("same hist");
      leg->AddEntry(hPDGgam, "Gamma", "L");
   }

   if (hPDGli->Integral()) { 
      hPDGli->Draw("same hist"); 
      leg->AddEntry(hPDGli, "Litium", "L");
   }

   if (hPDGdeu->Integral()) { 
      hPDGdeu->Draw("same hist"); 
      leg->AddEntry(hPDGdeu, "Deuterium", "L");
   }

   if (hPDGtri->Integral()) { 
      hPDGtri->Draw("same hist"); 
      leg->AddEntry(hPDGtri, "Tritium", "L");
   }

   if (hPDGsi->Integral()) { 
      hPDGsi->Draw("same hist"); 
      leg->AddEntry(hPDGsi, "Silicon", "L");
   }

   if (hPDGmg->Integral()) { 
      hPDGmg->Draw("same hist"); 
      leg->AddEntry(hPDGmg, "Magnesium", "L");
   }

   if (hPDGal->Integral()) { 
      hPDGal->Draw("same hist"); 
      leg->AddEntry(hPDGal, "Aluminum", "L");
   }

   leg->SetTextFont(22);
   leg->SetTextSize(0.04);
   leg->Draw();
 
   // ADD % TO AXIS
   TText *t = new TText();
   hPDGpro->GetYaxis()->SetLabelOffset(5);
   t->SetTextAlign(32);
   t->SetTextSize(0.04);
   t->SetTextFont(22);
   Float_t at;
   for (Int_t i=0; i<8;i++) {
      // max > i*6
      // i < max/6
      
      at = i/8. * 20; // max(hPDGhe->GetMaximum(), max(hPDGpro->GetMaximum(), hPDGhe3->GetMaximum()));
      t->DrawText(-0.42, at, Form("%.1f%%", at));
   }

}

void showOutputFileForImageReconstruction() {
   TH2F * hSimpleImagePrim = new TH2F("hSimpleImagePrim", "CHeCT projection (only primaries);X [mm];Y [mm]", 150, -75, 75, 150, -75,75);
   TH2F * hSimpleImageNormPrim = new TH2F("hSimpleImageNormPrim", "CHeCT projection (only primaries);X [mm];Y [mm]", 150, -75, 75, 150, -75,75);
   TH2F * hSimpleImage = new TH2F("hSimpleImage", "CHeCT projection;X [mm];Y [mm]", 150, -75, 75, 150, -75,75);
   TH2F * hSimpleImageNorm = new TH2F("hSimpleImageNorm", "CHeCT projection;X [mm];Y [mm]", 150, -75, 75, 150, -75,75);

   TFile *fIn = new TFile("OutputFiles/CHeCT_PedHead_190MeVu_15x15cm2.root");
   TTree *tIn = (TTree*) fIn->Get("WEPLData");

   Float_t x2x, x2y, wepl;
   Bool_t isSecondary;

   tIn->SetBranchAddress("X2x", &x2x);
   tIn->SetBranchAddress("X2y", &x2y);
   tIn->SetBranchAddress("wepl", &wepl);
   tIn->SetBranchAddress("isSecondary", &isSecondary);

   for (Int_t i=0; i<tIn->GetEntriesFast(); i++) {
      hSimpleImage->Fill(x2x, x2y, wepl);
      hSimpleImageNorm->Fill(x2x, x2y);

      if (!isSecondary) {
         hSimpleImagePrim->Fill(x2x, x2y, wepl);
         hSimpleImageNormPrim->Fill(x2x, x2y);
      }
   }

   hSimpleImage->Divide(hSimpleImageNorm);
   hSimpleImagePrim->Divide(hSimpleImageNormPrim);

   TCanvas *c;
   c->Divide(1,2,1e-5,1e-5);
   c->cd(1);
   hSimpleImage->Draw("COLZ");
   c->cd(2);
   hSimpleImagePrim->Draw("COLZ");

}

void makeOutputFileForImageReconstruction(Int_t Runs, Int_t tracksperrun, Int_t spotPosY) {
   run_energy = 760;
   run_degraderThickness = 180; // This number is useless I hope
   kEventsPerRun = tracksperrun;
   kDoTracking = true;

   Bool_t  kDraw = true;
   Float_t spotXFrom = -75;
   Float_t spotXTo = 75;
   Float_t spotXSpacing = 5;
   Float_t spotYFrom = -75;
   Float_t spotYTo = 75;
   Float_t spotYSpacing = 5;

   TH2F *hSimpleImage = nullptr;
   TH2F *hSimpleImageNorm = nullptr;

   if (kDraw) {
      hSimpleImage = new TH2F("hSimpleImage", "Simple Image;X [mm];Y [mm]", 150, -75, 75, 150, -75,75);
      hSimpleImageNorm = new TH2F("hSimpleImageNorm", "Simple Image;X [mm];Y [mm]", 150, -75, 75, 150, -75,75);
   }

   Float_t outSpotX, outSpotY, outWEPL, outX2x, outX2y, outP2x, outP2y, residualRange;
   Track * thisTrack = nullptr;

   TFile *fOut = new TFile(Form("OutputFiles/PedHeadCarbonHelium_spotsy_%03d.root", spotPosY), "recreate");
   TTree *tOut = new TTree("WEPLData", "Carbon+Helium beam");

   tOut->Branch("spotX", &outSpotX, "spotX/F");
   tOut->Branch("spotY", &outSpotY, "spotY/F");
   tOut->Branch("WEPL", &outWEPL, "WEPL/F");
   tOut->Branch("X2x", &outX2x, "X2x/F");
   tOut->Branch("X2y", &outX2y, "X2y/F");
   tOut->Branch("P2x", &outP2x, "P2x/F");
   tOut->Branch("P2y", &outP2y, "P2y/F");

   Int_t nruns = 0;

   for (float spotX = spotXFrom; spotX <= spotXTo; spotX += spotXSpacing) {
//      for (float spotY = spotYFrom; spotY <= spotYTo; spotY += spotYSpacing) {
//         printf("Spot (%.0f,%.0f)\n", spotX, spotY);
      Tracks * tracks = loadOrCreateTracks(true, Runs, false, run_energy, spotX, spotPosY);
      tracks->removeEmptyTracks();

      tracks->removeHighAngleTracks(30);
      tracks->removeNuclearInteractions();
      tracks->fillOutIncompleteTracks();
      tracks->doTrackFit();
      tracks->removeThreeSigmaShortTracks();
      tracks->removeEmptyTracks();

      for (Int_t i=0; i<=tracks->GetEntriesFast(); i++) {
         thisTrack = tracks->At(i);
         if (!thisTrack) continue;

         outSpotX = spotX;
         outSpotY = spotPosY;
         residualRange = getWEPLFromTL(thisTrack->getFitParameterRange());
         if (isnan(residualRange)) continue;

         outWEPL = 238.8 - residualRange; // 190 MeV/u Helium
         outX2x = thisTrack->At(0)->getXmm();
         outX2y = thisTrack->At(0)->getYmm();
         outP2x = (thisTrack->At(1)->getXmm() - outX2x) / dz;
         outP2y = (thisTrack->At(1)->getYmm() - outX2y) / dz;

         tOut->Fill();

         if (kDraw) {
            hSimpleImage->Fill(outX2x, outX2y, outWEPL);
            hSimpleImageNorm->Fill(outX2x, outX2y);
         }
      }
      nruns++;
//      }
   }
   
   fOut->Write();
//   printf("Image Reconstruction Input data written to Output/PedHeadCarbonHelium.root.\n");
   
   if (kDraw) {
      hSimpleImage->Divide(hSimpleImageNorm);
      hSimpleImage->Draw("COLZ");
   }
}

void drawTracks3D(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t switchLayer, Float_t energy, Float_t degraderThickness, Int_t tracksperrun, Bool_t doTracking) {
   run_energy = energy;
   run_degraderThickness = degraderThickness;
   kEventsPerRun = tracksperrun;
   kDoTracking = doTracking;
   
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy, 0, 0);
   tracks->removeEmptyTracks();

   printf("Found %d tracks before filtering.\n", tracks->GetEntries());

   Int_t numberOfPrimaries = 0;
   Int_t numberOfSecondaries = 0;
   for (int i=0; i<tracks->GetEntriesFast(); i++) {
      if (!tracks->At(i)) continue;
      if (!tracks->At(i)->Last()->isSecondary()) numberOfPrimaries++;
      else numberOfSecondaries++;
   }
   cout << "Number of primaries = " << numberOfPrimaries << ", number of secondaries = " << numberOfSecondaries << endl;
 
   tracks->removeHighAngleTracks(30); // mrad -> 30 He
   tracks->removeNuclearInteractions();
   tracks->removeEmptyTracks();

   tracks->fillOutIncompleteTracks(0.05);
   tracks->doTrackFit();

//   tracks->removeHighChiSquare(2500);
   tracks->removeThreeSigmaShortTracks();
   tracks->removeEmptyTracks();

   Bool_t   kDraw = true;

   TH1I  *hWEPLCorrect = new TH1I("hWEPLCorrect", ";Range in detector [mm WEPL];Frequency", 400, 0, 250);
   TH1I  *hWEPLSecondary = new TH1I("hWEPLSecondary", "", 400, 0, 250);
   TH1I  *hWEPLConfused = new TH1I("hWEPLConfused", "", 400, 0, 250);
   TH1I  *hProjConfused = new TH1I("hProjConfused", ";Projected error on phantom 10 cm from tracker;Frequency", 100, 0, 50);

   Float_t means[10];
   Float_t sigmas[10];

   TCanvas *c1 = nullptr;
   if (kDraw) {
      c1 = new TCanvas("c1", "c1", 1000, 1000);
      c1->SetTitle(Form("Tracks from %.2f MeV protons on %s", energy, getMaterialChar()));
   }
  
   TView *view = nullptr; 
   if (kDraw) view = TView::CreateView(1);
   float fromx = 0.1 * nx;
   float tox = 0.9 * nx;
   float fromy = 0.1 * ny;
   float toy = 0.9 * ny;

   /*
   fromy = 0, toy = ny;
   fromx = 0, tox = nx;
   */

   Int_t zoom = 750; // 750

   fromx = nx/2 - zoom*2;
   fromy = ny/2 - zoom*2;
   tox = nx/2 + zoom*2;
   toy = ny/2 + zoom*2;

   Int_t iret;
   Float_t theta = -3.5;
   Float_t phi = 89.5;

   if (kDraw) {
      view->SetRange(fromx, 0, fromy, tox, 35, toy);
      view->SetView(theta, phi, 0, iret);
   }

   TClonesArray *restPoints = tracks->getClustersWithoutTrack();
   Clusters * conflictClusters = nullptr;

   Int_t nClusters = 0;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      if (!tracks->At(i)) continue;
      nClusters += tracks->GetEntriesFast(i);
   }

   TPolyMarker3D *pMarker = new TPolyMarker3D(restPoints->GetEntriesFast(), 7);
   TPolyMarker3D *EIDMarker = new TPolyMarker3D(nClusters, 7);
   TPolyMarker3D *conflictMarker = new TPolyMarker3D(nClusters, 7);
   pMarker->SetMarkerColor(kBlue); // Missing cluster
//   pMarker->SetMarkerStyle(15);
   EIDMarker->SetMarkerColor(kGreen);
   conflictMarker->SetMarkerColor(kRed); // Conflicting cluster
   
   for (Int_t i=0; i<restPoints->GetEntriesFast(); i++) {
      if (!restPoints->At(i))
         continue;

      Cluster *thisCluster = (Cluster*) restPoints->At(i);
      Float_t x = thisCluster->getX();
      Float_t z = thisCluster->getY();
      Float_t y = thisCluster->getLayer();
   
      if (thisCluster->isSecondary()) {
         EIDMarker->SetPoint(i,x,y,z);
      }
      else {
         pMarker->SetPoint(i, x, y, z);
      }
   }

   printf("There are %d unused clusters.\n", restPoints->GetEntries());

   if (kDraw) pMarker->Draw();
   
   Int_t ntracks = tracks->GetEntriesFast();
   Int_t EIDidx = 0;
   Int_t conflictIdx = 0;

   Int_t medianEventID = -1;

   Int_t nTrueTracks = 0;
   Int_t nOKTracks = 0;
   Int_t nOKTracksAllClusters = 0;
   Int_t nOKTracksAllClustersOK2nd = 0;
   Int_t nOKMinusTracks = 0;
   Int_t nOKLastLayers = 0;
   Int_t nOneWrong = 0;
   Int_t nMissingEID;
   Int_t nOkPrimary = 0;
   Int_t nPrimaryIncomplete = 0;
   Int_t nPrimaryConfused = 0;
   Int_t nPrimaryConfusedIncomplete = 0;
   Int_t nOkSecondary = 0;
   Int_t nSecondaryIncomplete = 0;
   Int_t nSecondaryConfused = 0;
   Int_t nSecondaryConfusedIncomplete = 0;
   numberOfPrimaries = 0;
   numberOfSecondaries = 0;
   Float_t correctPosx;
   Float_t thisPosx;
   Float_t correctPosy;
   Float_t thisPosy;
   Float_t distToPhantom = 100;
   Cluster *a = nullptr;
   Cluster *b = nullptr;
   Cluster *aTrue = nullptr;
   Cluster *bTrue = nullptr;

   Track * thisTrack = nullptr;
//   tracks->createEIDSortList();

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      if (isnan(thisTrack->getFitParameterRange())) continue;

      nMissingEID = tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), i); 

      if (!thisTrack->Last()->isSecondary() && !thisTrack->isSecondary(0)) { // primary
         numberOfPrimaries++;

         if (thisTrack->isFirstAndLastEventIDEqual()) { // No confused tracks
            if (nMissingEID == 0) { // NO missing clusters
               nOkPrimary++;
               hWEPLCorrect->Fill(getUnitFromTL(thisTrack->getFitParameterRange()));
            }
            
            else { // Missing clusters
               nPrimaryIncomplete++;
               hWEPLConfused->Fill(getUnitFromTL(thisTrack->getFitParameterRange()));
            }

         }

         else { // Confused tracks
            a = thisTrack->At(0);
            b = thisTrack->At(1);
            int eid = thisTrack->Last()->getEventID();
            // Track *corTrack = tracks->getTrackWithEID(eid);
            Track * corTrack = nullptr;
            if (corTrack) {
               aTrue = corTrack->At(0);
               bTrue = corTrack->At(1);

               if (a && b && aTrue && bTrue) {
                  float deltax = b->getXmm() - a->getXmm();
                  float deltay = b->getYmm() - a->getYmm();
                  float deltaz = b->getLayermm() - a->getLayermm();
                  thisPosx = a->getXmm() - (deltax/deltaz) * distToPhantom;
                  thisPosy = a->getYmm() - (deltay/deltaz) * distToPhantom;

                  deltax = bTrue->getXmm() - aTrue->getXmm();
                  deltay = bTrue->getYmm() - aTrue->getYmm();
                  deltaz = bTrue->getLayermm() - aTrue->getLayermm();
                  correctPosx = aTrue->getXmm() - (deltax/deltaz) * distToPhantom;
                  correctPosy = aTrue->getYmm() - (deltay/deltaz) * distToPhantom;

                  float delta = sqrt(pow(thisPosx - correctPosx, 2) + pow(thisPosy - correctPosy, 2));
                  hProjConfused->Fill(delta);
               }
            }
            if (nMissingEID == 0) { // No missing clusters
               nPrimaryConfused++;
            }

            else { // Confused and missing clusters
               nPrimaryConfusedIncomplete++;
            }
         }
      }

      else { // secondary
         numberOfSecondaries++;
         hWEPLSecondary->Fill(getUnitFromTL(thisTrack->getFitParameterRange()));

         if (thisTrack->isFirstAndLastEventIDEqual()) { // No confused tracks
            if (nMissingEID == 0) { // No missing clusters
               nOkSecondary++;
            }

            else { // Missing clusters
               nSecondaryIncomplete++;
            }
         }

         else { // Confused tracks
            if (nMissingEID == 0) { // No missing clusters
               nSecondaryConfused++;
            }
            else { // Confused AND Missing clusters
               nSecondaryConfusedIncomplete++;
            }
         }
      }

      if (thisTrack->isOneEventID()) nTrueTracks++;
      
      
      if (thisTrack->isFirstAndLastEventIDEqual()) nOKTracks++;
      else {
         if (!thisTrack->Last()) continue;
         Int_t lastEID = thisTrack->Last()->getEventID();   
         Bool_t lastSecondary = thisTrack->Last()->isSecondary();
         if (lastSecondary) continue;

         Int_t trackID = tracks->getTrackIdxFromFirstLayerEID(lastEID);

         if (trackID < 0) continue;
         if (!tracks->At(trackID)->At(0)) continue;
   
         Float_t delta = diffmmXY(tracks->At(trackID)->At(0), thisTrack->At(0));
         Float_t phi0 = thisTrack->getSlopeAngleBetweenLayers(1);
         Float_t phi1 = tracks->At(trackID)->getSlopeAngleBetweenLayers(1);

         Float_t deltaphi = fabs(phi0 - phi1);

         if (delta < 0.5 && deltaphi < 1) {
            nOKMinusTracks++;
            nOKLastLayers++;
         }
         else if (thisTrack->getWEPL() < 0.2 * getWEPLFromEnergy(run_energy)) {
            // Bad track ends early. OK...
            nOKLastLayers++;
         }
      }
   }

   nOKMinusTracks += nOKTracks;
   nOKLastLayers += nOKTracks;

   Int_t numberOfTracks = tracks->GetEntries();

   cout << endl << "Total number of tracks = " << numberOfPrimaries + numberOfSecondaries << endl;
   cout << "Number of primaries = " << numberOfPrimaries << ", number of secondaries = " << numberOfSecondaries << endl;
   cout << "Number of primaries tracked correctly = " << nOkPrimary << " (" << 100 * (float) nOkPrimary / numberOfPrimaries << "%)\n";
   cout << "Number of primaries confused = " << nPrimaryConfused << " (" << 100 * (float) nPrimaryConfused / numberOfPrimaries << "%)\n";
   cout << "Number of primaries incompletely tracked = " << nPrimaryIncomplete << " (" << 100 * (float) nPrimaryIncomplete / numberOfPrimaries << "%)\n";
   cout << "Number of primaries confused + incompletely tracked = " << nPrimaryConfusedIncomplete << " (" << 100 * (float) nPrimaryConfusedIncomplete / numberOfPrimaries << "%)\n";
   
   cout << endl << endl;
   cout << "Total number of primaries = " << numberOfPrimaries << endl;
   cout << "Primaries tracked correctly = " << nOkPrimary << " (" << 100 * (float) nOkPrimary / numberOfPrimaries;
   cout << "%),  Secondaries tracked correctly = " << nOkSecondary << " (" << 100 * (float) nOkSecondary / numberOfSecondaries << "%)\n";

   /*
   cout << endl << endl;
   cout << "Total number of tracks = " << numberOfPrimaries + numberOfSecondaries << endl;
   cout << "Tracks tracked correctly = " << nOkPrimary + nOkSecondary << " (" << 100 * (float) ( nOkPrimary + nOkSecondary) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nOkPrimary << " (" << 100 * (float) nOkPrimary / (nOkPrimary + nOkSecondary);
   cout << "%) are primaries, and " << nOkSecondary << " (" << 100 * (float) nOkSecondary / (nOkSecondary + nOkPrimary) << "%) are secondaries.\n";
   cout << "Tracks confused = " << nPrimaryConfused + nSecondaryConfused << " (" << 100 * (float) (nPrimaryConfused + nSecondaryConfused) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nPrimaryConfused << " (" << 100 * (float) nPrimaryConfused / (nPrimaryConfused + nSecondaryConfused);
   cout << "%) are primaries and " << nSecondaryConfused << " (" << 100 * (float) nSecondaryConfused / (nPrimaryConfused + nSecondaryConfused);
   cout << "%) are secondaries.\n";
   cout << "Tracks incompletely tracked = " << nPrimaryIncomplete + nSecondaryIncomplete << " (" << 100 * (float) (nPrimaryIncomplete + nSecondaryIncomplete) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nPrimaryIncomplete << " (" << 100 * (float) nPrimaryIncomplete / (nPrimaryIncomplete + nSecondaryIncomplete);
   cout << "%) are primaries and " << nSecondaryIncomplete << " (" << 100 * (float) nSecondaryIncomplete/ (nPrimaryIncomplete+ nSecondaryIncomplete);
   cout << "%) are secondaries.\n";
   cout << "Tracks confused and incompletely tracked = " << nPrimaryConfusedIncomplete + nSecondaryConfusedIncomplete << " (" << 100 * (float) (nPrimaryConfusedIncomplete + nSecondaryConfusedIncomplete) / (numberOfPrimaries + numberOfSecondaries) << "%), of which " << nPrimaryConfusedIncomplete << " (" << 100 * (float) nPrimaryConfusedIncomplete / (nPrimaryConfusedIncomplete + nSecondaryConfusedIncomplete);
   cout << "%) are primaries and " << nSecondaryConfusedIncomplete << " (" << 100 * (float) nSecondaryConfusedIncomplete/ (nPrimaryConfusedIncomplete+ nSecondaryConfusedIncomplete);
   cout << "%) are secondaries.\n";

   Float_t factorEIDOK = 100 * ((float) nOKTracks / numberOfTracks);
   Float_t factorEIDOKAllClusters = 100 * ((float) nOKTracksAllClusters / numberOfTracks);
   Float_t factorEIDOKAllClustersOK2nd = 100 * ((float) nOkPrimary / numberOfTracks);
   Float_t factorEIDOKMinus = 100 * ((float) nOKMinusTracks / numberOfTracks);
   Float_t factorLastLayers = 100 * ((float) nOKLastLayers / numberOfTracks);

   cout << endl << endl;
   cout << nOKTracks << " of total " << numberOfTracks << " tracks has the same first/last ID (" << factorEIDOK << "%)\n";
   cout << nOKTracksAllClustersOK2nd << " of total " << numberOfTracks << " track has first/last event ID + no missing clusters (" << factorEIDOKAllClustersOK2nd << "%)\n";
   cout << nOKMinusTracks << " of total " << numberOfTracks << " tracks has a close match (0.5 mm, 1 degree) on first / last cluster (" << factorEIDOKMinus << "%)\n";
   cout << nOKLastLayers << " of total " << numberOfTracks << " tracks has a close match (0.5 mm, 1 degree) or is a very short track (" << factorLastLayers << "%)\n";
   */

   Int_t badSecondary = 0;
   Int_t badPrimary = 0;
   Int_t okPrimary = 0;
   Int_t badShort = 0;
   Int_t incompletePrimary = 0;

   if (kDraw) {
      for (Int_t i=0; i<ntracks; i++) {
         Track *thisTrack = tracks->At(i);
         if (!thisTrack) continue;

         Int_t n = thisTrack->GetEntriesFast();

         TPolyLine3D *l = new TPolyLine3D(n);
         l->SetLineWidth(2);
         TPolyMarker3D *trackPoints = new TPolyMarker3D(nClusters, 7);
         nMissingEID = tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer(), i);
         if (!thisTrack->isFirstAndLastEventIDEqual()) {
            l->SetLineColor(kRed);
         }

         else if (nMissingEID>0) {
            l->SetLineColor(kGray);
         }
        
         if (thisTrack->isSecondary(0) || thisTrack->Last()->isSecondary()) {
            l->SetLineColor(kGreen);
         }

         Int_t lineElementNumber = 0;
         Int_t pointNumber = 0;
         for (Int_t j=0; j<n; j++) {
            if (!thisTrack->At(j)) continue;

            Float_t x = thisTrack->getX(j);
            Float_t z = thisTrack->getY(j);
            Float_t y = thisTrack->getLayer(j);
            
            if (thisTrack->getLayer(j) < switchLayer) {
               l->SetPoint(lineElementNumber++,x,y,z);
            }
            else {
               trackPoints->SetPoint(pointNumber++, x, y, z);
            }
         }

         conflictClusters = (Clusters*) thisTrack->getConflictClusters();
         for (Int_t j=0; j<conflictClusters->GetEntriesFast(); j++) {
            if (!conflictClusters->At(j)) continue;
            
            Float_t x = conflictClusters->getX(j);
            Float_t z = conflictClusters->getY(j);
            Float_t y = conflictClusters->getLayer(j);
            
            conflictMarker->SetPoint(conflictIdx++, x,y,z);
         }
         if (kDraw) {
   //       l->SetLineColor(kBlack);
            l->SetLineWidth(3);
            if (l->GetLineColor() == kRed) l->Draw();
            if (l->GetLineColor() == kGreen) {
               l->SetLineWidth(2);
               l->Draw();
            }
            if (l->GetLineColor() == kGray) l->Draw();
            l->Draw();
         }

         if (l->GetLineColor() == kGreen) badSecondary++;
         if (l->GetLineColor() == kRed)   badPrimary++;
         else if (l->GetLineColor() == kGray) incompletePrimary++;
         else okPrimary++;
         
         if (kDraw) {
            trackPoints->Draw();
   //         EIDMarker->Draw();
   //      conflictMarker->Draw();
         }
      }
      
      view->ShowAxis(); // comment for pure display
      c1->Update();

      TAxis3D *axis = TAxis3D::GetPadAxis();
      axis->SetTitle("3D view of tracks and clusters");
      axis->SetLabelColor(kBlack);
      axis->SetAxisColor(kBlack);
      axis->SetXTitle("Pixels in X");
      axis->SetYTitle("Layer number");
      axis->SetZTitle("Pixels in Y");
      axis->SetLabelSize(0.025);
      axis->SetTitleOffset(2);

   }

   TCanvas *c2 = new TCanvas();
   TH1I *changeEID = new TH1I("changeEID", "Change EID;Layer;EID changes", 100, 0, 100);
   for (int i=0; i<tracks->GetEntriesFast(); i++) {
      Track * thisTrack = tracks->At(i);
      if (!thisTrack) continue;
      int firstEID = thisTrack->getEventID(0);
      for (int j=0; j<thisTrack->GetEntriesFast(); j++) {
         if (thisTrack->getEventID(j) != firstEID) {
            changeEID->Fill(thisTrack->getLayer(j));
            break;
         }
      }
   }
   changeEID->Draw();

   Int_t badTotal = badPrimary + badSecondary;
   printf("Of %d bad tracks, %d (%.2f %%) are primaries and %d (%.2f %%) are secondaries.\n", badTotal, badPrimary, 100 * float(badPrimary) / badTotal, badSecondary, 100 * float(badSecondary) / badTotal);
   printf("okPrimary = %d\n", okPrimary);
   printf("Incomplete primary = %d\n", incompletePrimary);

   /*
   vector<Int_t> * conflictTracks = tracks->getTracksWithConflictClusters();
   vector<Int_t> * oneConflictPair = nullptr;
   vector<Int_t> * allConflictPairs = new vector<Int_t>;
   
   for (UInt_t i=0; i<conflictTracks->size(); i++) {
      oneConflictPair = tracks->getConflictingTracksFromTrack(conflictTracks->at(i));

      Int_t idx0 = oneConflictPair->at(0);
      Int_t idx1 = -1;
      if (oneConflictPair->size() > 1) {
         idx1 = oneConflictPair->at(1);
      }

      allConflictPairs->push_back(idx0);
      allConflictPairs->push_back(idx1);
   }

   // PRINTING
   
   cout << "Found the following tracks with conflicting clusters: ";
   for (UInt_t i=0; i<conflictTracks->size(); i++) {
      cout << conflictTracks->at(i) << " (eventID " << tracks->At(conflictTracks->at(i))->getEventID(0) << "), ";
   }
   cout << "\n";

   for (UInt_t i=0; i<allConflictPairs->size() / 2; i++) {
      if (allConflictPairs->at(2*i+1) < 0) continue;
      
      Track * trackA = tracks->At(allConflictPairs->at(2*i));
      Track * trackB = tracks->At(allConflictPairs->at(2*i+1));

      cout << "Track pair number " << i+1 << " found is: \n\tTRACK A: ";
      for (Int_t j=0; j<trackA->GetEntriesFast(); j++) { 
         if ( ! trackA->At(j) ) continue;
         cout << *trackA->At(j) << ", ";
      }

      cout << "\n\tTRACK B: ";
      for (Int_t j=0; j<trackB->GetEntriesFast(); j++) { 
         if ( ! trackB->At(j) ) continue;
         cout << *trackB->At(j) << ", ";
      }
      cout << endl;
   }
   */

   if (kDraw) c1->SaveAs(Form("OutputFiles/figures/testOutput_switchLayer%d.png", switchLayer));
/*
   TCanvas *c2 = new TCanvas();
   c2->Divide(2,1,1e-5,1e-5);
   c2->cd(1);
   hWEPLCorrect->SetLineColor(kBlack);
   hWEPLCorrect->Draw();
   hWEPLConfused->SetLineColor(kRed);
   hWEPLConfused->Draw("same");
   hWEPLSecondary->SetLineColor(kGreen);
   hWEPLSecondary->Draw("same");

   TLegend *l = new TLegend(0.3,0.3,0.4,0.4);
   l->AddEntry(hWEPLCorrect, "Good tracks", "L");
   l->AddEntry(hWEPLSecondary, "Nuclear interacting tracks");
   l->AddEntry(hWEPLConfused, "Tracks confused during recon");
   l->Draw();

   c2->cd(2);
   hProjConfused->SetLineColor(kRed);
   hProjConfused->Draw();
*/
   delete tracks;
//   delete conflictTracks;
//   delete allConflictPairs;
}

void drawAlignmentCheck(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   kDataType = dataType;

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);

   TCanvas   * c1 = new TCanvas("c1", "Alignment check", 1800, 1500);
   c1->Divide(2, 4, 0.01, 0.01);
   gStyle->SetOptStat(0);

   TH1F      * alignmentX1 = new TH1F("alignmentX1", "Alignment check in X1 direction;Layer;Mean alignment in x1 [mm]", 8, 0, 8);
   TH1F      * alignmentY1 = new TH1F("alignmentY1", "Alignment check in Y1 direction;Layer;Mean alignment in y1 [mm]", 8, 0, 8);
   TH1F      * alignmentX2 = new TH1F("alignmentX2", "Alignment check in X2 direction;Layer;Mean alignment in x2 [mm]", 8, 0, 8);
   TH1F      * alignmentY2 = new TH1F("alignmentY2", "Alignment check in Y2 direction;Layer;Mean alignment in y2 [mm]", 8, 0, 8);
   TH1F      * alignmentX3 = new TH1F("alignmentX3", "Alignment check in X3 direction;Layer;Mean alignment in x3 [mm]", 8, 0, 8);
   TH1F      * alignmentY3 = new TH1F("alignmentY3", "Alignment check in Y3 direction;Layer;Mean alignment in y3 [mm]", 8, 0, 8);
   TH1F      * alignmentX4 = new TH1F("alignmentX4", "Alignment check in X4 direction;Layer;Mean alignment in x4 [mm]", 8, 0, 8);
   TH1F      * alignmentY4 = new TH1F("alignmentY4", "Alignment check in Y4 direction;Layer;Mean alignment in y4 [mm]", 8, 0, 8);
   TH1F      * normX1 = new TH1F("normX1", "Alignment check in X1 direction", 8, 0, 8);
   TH1F      * normY1 = new TH1F("normY1", "Alignment check in Y1 direction", 8, 0, 8);
   TH1F      * normX2 = new TH1F("normX2", "Alignment check in X2 direction", 8, 0, 8);
   TH1F      * normY2 = new TH1F("normY2", "Alignment check in Y2 direction", 8, 0, 8);
   TH1F      * normX3 = new TH1F("normX3", "Alignment check in X3 direction", 8, 0, 8);
   TH1F      * normY3 = new TH1F("normY3", "Alignment check in Y3 direction", 8, 0, 8);
   TH1F      * normX4 = new TH1F("normX4", "Alignment check in X4 direction", 8, 0, 8);
   TH1F      * normY4 = new TH1F("normY4", "Alignment check in Y4 direction", 8, 0, 8);

   Track     * thisTrack = nullptr;
   Cluster   * thisCluster = nullptr;
   Cluster   * nextCluster = nullptr;
   Int_t       thisLayer;
   Float_t     x, y, xp, yp, diffx, diffy;
   Int_t       chip = 0;
   vector<Float_t> deflection;

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);

      for (Int_t j=0; j<thisTrack->GetEntriesFast()-1; j++) {
         thisCluster = thisTrack->At(j);
         nextCluster = thisTrack->At(j+1);

         if (thisCluster && nextCluster) {
            x = thisCluster->getXmm();
            y = thisCluster->getYmm();
            xp = nextCluster->getXmm();
            yp = nextCluster->getYmm();

            diffx = xp - x;
            diffy = yp - y;
            thisLayer = nextCluster->getLayer();

            deflection = thisTrack->getLateralDeflectionFromExtrapolatedPosition(thisLayer);
            diffx = deflection.at(0);
            diffy = deflection.at(1);

            if       (x >  0 && y >  0) chip = 0;
            else if  (x <= 0 && y >  0) chip = 3;
            else if  (x <= 0 && y <= 0) chip = 2;
            else if  (x >  0 && y <= 0) chip = 1;

            if (thisLayer>9) continue;

            if (chip == 0) {
               alignmentX1->Fill(thisLayer, diffx);
               alignmentY1->Fill(thisLayer, diffy);
               normX1->Fill(thisLayer);
               normY1->Fill(thisLayer);
            }
            else if (chip == 1) {
               alignmentX2->Fill(thisLayer, diffx);
               alignmentY2->Fill(thisLayer, diffy);
               normX2->Fill(thisLayer);
               normY2->Fill(thisLayer);
            }

            else if (chip == 2) {
               alignmentX3->Fill(thisLayer, diffx);
               alignmentY3->Fill(thisLayer, diffy);
               normX3->Fill(thisLayer);
               normY3->Fill(thisLayer);
            }

            else if (chip == 3) {
               alignmentX4->Fill(thisLayer, diffx);
               alignmentY4->Fill(thisLayer, diffy);
               normX4->Fill(thisLayer);
               normY4->Fill(thisLayer);
            }
         }
      }
   }

   alignmentX1->Divide(normX1);
   alignmentY1->Divide(normY1);
   alignmentX2->Divide(normX2);
   alignmentY2->Divide(normY2);
   alignmentX3->Divide(normX3);
   alignmentY3->Divide(normY3);
   alignmentX4->Divide(normX4);
   alignmentY4->Divide(normY4);

   c1->cd(1);
   alignmentX1->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentX1->GetYaxis()->SetTitleOffset(1.2);
   alignmentX1->GetYaxis()->SetTitleFont(22);
   alignmentX1->GetYaxis()->SetLabelFont(22);
   alignmentX1->GetXaxis()->SetTitleFont(22);
   alignmentX1->GetXaxis()->SetLabelFont(22);
   alignmentX1->SetFillColor(kGreen+3);
   alignmentX1->SetLineColor(kBlack);
   alignmentX1->Draw();
   gPad->Update();
   TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();

   c1->cd(2);
   alignmentY1->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentY1->GetYaxis()->SetTitleOffset(1.2);
   alignmentY1->GetYaxis()->SetTitleFont(22);
   alignmentY1->GetYaxis()->SetLabelFont(22);
   alignmentY1->GetXaxis()->SetTitleFont(22);
   alignmentY1->GetXaxis()->SetLabelFont(22);
   alignmentY1->SetFillColor(kGreen+3);
   alignmentY1->SetLineColor(kBlack);
   alignmentY1->Draw();
   gPad->Update();
   TPaveText *title2 = (TPaveText*) gPad->GetPrimitive("title");
   title2->SetTextFont(22);
   gPad->Modified();

   c1->cd(3);
   alignmentX2->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentX2->GetYaxis()->SetTitleOffset(1.2);
   alignmentX2->GetYaxis()->SetTitleFont(22);
   alignmentX2->GetYaxis()->SetLabelFont(22);
   alignmentX2->GetXaxis()->SetTitleFont(22);
   alignmentX2->GetXaxis()->SetLabelFont(22);
   alignmentX2->SetFillColor(kGreen+3);
   alignmentX2->SetLineColor(kBlack);
   alignmentX2->Draw();
   gPad->Update();
   TPaveText *title3 = (TPaveText*) gPad->GetPrimitive("title");
   title3->SetTextFont(22);
   gPad->Modified();

   c1->cd(4);
   alignmentY2->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentY2->GetYaxis()->SetTitleOffset(1.2);
   alignmentY2->GetYaxis()->SetTitleFont(22);
   alignmentY2->GetYaxis()->SetLabelFont(22);
   alignmentY2->GetXaxis()->SetTitleFont(22);
   alignmentY2->GetXaxis()->SetLabelFont(22);
   alignmentY2->SetFillColor(kGreen+3);
   alignmentY2->SetLineColor(kBlack);
   alignmentY2->Draw();
   gPad->Update();
   TPaveText *title4 = (TPaveText*) gPad->GetPrimitive("title");
   title4->SetTextFont(22);
   gPad->Modified();

   c1->cd(5);
   alignmentX3->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentX3->GetYaxis()->SetTitleOffset(1.2);
   alignmentX3->GetYaxis()->SetTitleFont(22);
   alignmentX3->GetYaxis()->SetLabelFont(22);
   alignmentX3->GetXaxis()->SetTitleFont(22);
   alignmentX3->GetXaxis()->SetLabelFont(22);
   alignmentX3->SetFillColor(kGreen+3);
   alignmentX3->SetLineColor(kBlack);
   alignmentX3->Draw();
   gPad->Update();
   TPaveText *title5 = (TPaveText*) gPad->GetPrimitive("title");
   title5->SetTextFont(22);
   gPad->Modified();

   c1->cd(6);
   alignmentY3->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentY3->GetYaxis()->SetTitleOffset(1.2);
   alignmentY3->GetYaxis()->SetTitleFont(22);
   alignmentY3->GetYaxis()->SetLabelFont(22);
   alignmentY3->GetXaxis()->SetTitleFont(22);
   alignmentY3->GetXaxis()->SetLabelFont(22);
   alignmentY3->SetFillColor(kGreen+3);
   alignmentY3->SetLineColor(kBlack);
   alignmentY3->Draw();
   gPad->Update();
   TPaveText *title6 = (TPaveText*) gPad->GetPrimitive("title");
   title6->SetTextFont(22);
   gPad->Modified();
   
   c1->cd(7);
   alignmentX4->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentX4->GetYaxis()->SetTitleOffset(1.2);
   alignmentX4->GetYaxis()->SetTitleFont(22);
   alignmentX4->GetYaxis()->SetLabelFont(22);
   alignmentX4->GetXaxis()->SetTitleFont(22);
   alignmentX4->GetXaxis()->SetLabelFont(22);
   alignmentX4->SetFillColor(kGreen+3);
   alignmentX4->SetLineColor(kBlack);
   alignmentX4->Draw();
   gPad->Update();
   TPaveText *title7 = (TPaveText*) gPad->GetPrimitive("title");
   title7->SetTextFont(22);
   gPad->Modified();
   
   c1->cd(8);
   alignmentY4->GetYaxis()->SetRangeUser(-0.4, 0.4);
   alignmentY4->GetYaxis()->SetTitleOffset(1.2);
   alignmentY4->GetYaxis()->SetTitleFont(22);
   alignmentY4->GetYaxis()->SetLabelFont(22);
   alignmentY4->GetXaxis()->SetTitleFont(22);
   alignmentY4->GetXaxis()->SetLabelFont(22);
   alignmentY4->SetFillColor(kGreen+3);
   alignmentY4->SetLineColor(kBlack);
   alignmentY4->Draw();
   gPad->Update();
   TPaveText *title8 = (TPaveText*) gPad->GetPrimitive("title");
   title8->SetTextFont(22);
   gPad->Modified();

   if (dataType == kMC) {
      c1->SaveAs(Form("OutputFiles/figures/AlignmentCheck_MC_%.0fMeV.pdf", run_energy));
   }
   else {
      c1->SaveAs(Form("OutputFiles/figures/AlignmentCheck_EXP_Corrected_%.0fMeV.pdf", run_energy));
   }

   Float_t rmsX1 = 0;
   Float_t rmsX2 = 0;
   Float_t rmsY1 = 0;
   Float_t rmsY2 = 0;
   Float_t rmsX3 = 0;
   Float_t rmsX4 = 0;
   Float_t rmsY3 = 0;
   Float_t rmsY4 = 0;
   for (Int_t i=0; i<7; i++) {
      rmsX1 += pow(alignmentX1->GetBinContent(i+1), 2);
      rmsX2 += pow(alignmentX2->GetBinContent(i+1), 2);
      rmsX3 += pow(alignmentX3->GetBinContent(i+1), 2);
      rmsX4 += pow(alignmentX4->GetBinContent(i+1), 2);
      rmsY1 += pow(alignmentY1->GetBinContent(i+1), 2);
      rmsY2 += pow(alignmentY2->GetBinContent(i+1), 2);
      rmsY3 += pow(alignmentY3->GetBinContent(i+1), 2);
      rmsY4 += pow(alignmentY4->GetBinContent(i+1), 2);
   }

   rmsX1 = sqrt(rmsX1); rmsX2 = sqrt(rmsX2);
   rmsX3 = sqrt(rmsX3); rmsX4 = sqrt(rmsX4);
   rmsY1 = sqrt(rmsY1); rmsY2 = sqrt(rmsY2);
   rmsY3 = sqrt(rmsY3); rmsY4 = sqrt(rmsY4);

   cout << "The RMS values for the distributions are: \n";
   cout << "Total RMS = " << sqrt(pow(rmsX1, 2) + pow(rmsX2, 2) + pow(rmsX3, 2) + pow(rmsX4, 2) + pow(rmsY1, 2) + pow(rmsY2, 2) + pow(rmsY3, 2) + pow(rmsY4, 2)) << endl;

   ifstream in("Data/ExperimentalData/Alignment_mine.txt");

   chipAlignment chipAlignmentArray[96];
   chipAlignment chipFile;

   while (1) {
      in >> chipFile.idx >> chipFile.deltaX >> chipFile.deltaY >> chipFile.deltaTheta;
      if (!in.good()) break;
      chipFile.deltaX *= 10;
      chipFile.deltaY *= 10;

      chipAlignmentArray[chipFile.idx] = chipFile;
   }
   in.close();

   ofstream file("Data/ExperimentalData/Alignment_mine.txt");

   Float_t theta[32] = {0.0032, -0.00017, 0.0015, 0.0021, -0.0059, -0.007, -0.0037, -0.0032, -0.0056, -0.0060, 0.00016, -0.0007, -0.0003, -0.0014, -0.0006, -0.00015, 0.0008, 0.0007, -0.002, -0.005, -0.007, -0.049, -0.0038, -0.006, 0, 0, 0, 0, -0.0004, -0.00032, -0.007, -0.008};

   Float_t deltaX, deltaY;
   for (Int_t i=0; i<7; i++) {
      deltaX = chipAlignmentArray[i*4+0].deltaX;
      deltaY = chipAlignmentArray[i*4+0].deltaY;
      cout << Form("In layer %d, old value was (%.4f, %.4f), with corrections (%.4f, %.4f) new value is (%.4f, %.4f)\n", i, deltaX, deltaY, alignmentX4->GetBinContent(i+1), alignmentY4->GetBinContent(i+1), deltaX - alignmentX1->GetBinContent(i+1), deltaY - alignmentY1->GetBinContent(i+1));
      file << i*4 + 0 << " " << -0.1 * (alignmentX1->GetBinContent(i+1) - deltaX) << " " << -0.1 * (alignmentY1->GetBinContent(i+1) - deltaY) << " " << theta[i*4+0] << endl;
      deltaX = chipAlignmentArray[i*4+1].deltaX;
      deltaY = chipAlignmentArray[i*4+1].deltaY;
      file << i*4 + 1 << " " << -0.1 * (alignmentX2->GetBinContent(i+1) - deltaX) << " " << -0.1 * (alignmentY2->GetBinContent(i+1) - deltaY) << " " << theta[i*4+1] << endl;
      deltaX = chipAlignmentArray[i*4+2].deltaX;
      deltaY = chipAlignmentArray[i*4+2].deltaY;
      file << i*4 + 2 << " " << -0.1 * (alignmentX3->GetBinContent(i+1) - deltaX) << " " << -0.1 * (alignmentY3->GetBinContent(i+1) - deltaY) << " " << theta[i*4+2] << endl;
      deltaX = chipAlignmentArray[i*4+3].deltaX;
      deltaY = chipAlignmentArray[i*4+3].deltaY;
      file << i*4 + 3 << " " << -0.1 * (alignmentX4->GetBinContent(i+1) - deltaX) << " " << -0.1 * (alignmentY4->GetBinContent(i+1) - deltaY) << " " << theta[i*4+3] << endl;
   }
   file.close();
}

void drawDiffusionCheck(Int_t Runs, Int_t Layer, Float_t energy) {
   run_energy = energy;
   DataInterface *di = new DataInterface();
   CalorimeterFrame *cf = new CalorimeterFrame();
   
   for (Int_t i=0; i<=Runs; i++) {
      di->getMCFrame(i, cf); // Remember to have MC data available at ./test.root
   }

   delete di;

   TCanvas *c1 = new TCanvas("c1", "multipads", 1400, 900);
   gStyle->SetOptStat(0);
   c1->Divide(2,1,0.01,0.01,0);
   
   TH2F *undiffusedTH2F = (TH2F*) cf->getTH2F(Layer)->Clone();
   undiffusedTH2F->SetName("undiffusedTH2F");
   
   cf->diffuseFrame(new TRandom3(0));
   
   TH2F *diffusedTH2F = (TH2F*) cf->getTH2F(Layer)->Clone();
   diffusedTH2F->SetName("diffusedTH2F");

   c1->cd(1);
   undiffusedTH2F->Draw("COLZ");
   
   c1->cd(2);
   diffusedTH2F->Draw("COLZ");
   
   c1->Update();
}

void drawFrame2D(Int_t dataType, Int_t layerNo, Float_t energy, Float_t degraderThickness) {
   run_energy = energy;
   run_degraderThickness = degraderThickness;

   printf("Creating a 2D histogram of layer %d using %s and %d number of primary protons\n", layerNo, getDataTypeChar(dataType), kEventsPerRun); 

   DataInterface *di = new DataInterface();
   Layer *l = new Layer(layerNo);
   TCanvas *c1 = new TCanvas("c1", Form("Hit distribution in layer %d", layerNo), 1200, 800);

/*
 * // Draw several layers at once
 *
   vector<TCanvas*> *cvec = new vector<TCanvas*>;
   for (Int_t i=0; i<nLayers; i++) {
      cvec->push_back(new TCanvas(Form("c%d", i), Form("Hit distribution in layer %d", i), 800, 800));
   }
*/
   
   if (dataType == kData) {
      di->getDataFrame(0, l, energy);
   }

   else {
      showDebug("Get MC frame...");
      di->getMCFrame(0, l);
      showDebug("OK!\nDiffuseLayer...");
      l->diffuseLayer(new TRandom3(0)); // Model the cluster diffusion process
      showDebug("OK!\n");
   }

   delete di;

   gStyle->SetOptStat(0);
   TH2F *Frame2D = l->getTH2F();
   Frame2D->Draw("COLZ");

/* 
 * // Draw several layers at once
 *
   for (Int_t i=0; i<6; i++) {
      cvec->at(i)->cd();
      Frame2D = cf->getTH2F(i);
      Frame2D->Draw("COLZ");
      Frame2D->GetXaxis()->SetRangeUser(600,900);
      Frame2D->GetYaxis()->SetRangeUser(400,700);
   }
*/

}

void drawDataProfile(Float_t energy) {
   run_energy = energy;

   TCanvas *c1 = new TCanvas("c1", "Beam profiles", 1200, 700);
   c1->Divide(2,1,0.0001,0.0001);

   DataInterface *di = new DataInterface();

   gStyle->SetOptStat(0);

   Float_t scaleFactor = 15;

   TH2F *hProfile = new TH2F("h2", "Beam profile in detector;Y posizion;Layer number", ny/scaleFactor,0,ny,20,0,20);
   TH2F *hProjection = new TH2F("hProjection", "Beam projection in detector;X position;Y position", nx/scaleFactor,0,nx, ny/scaleFactor,0,ny);
   di->getDataProfile(hProfile, hProjection, energy);

   c1->cd(1);
   hProfile->Draw("colz");

   c1->cd(2);
   hProjection->Draw("colz");
}

/*
void drawTrackingError(Int_t Runs, Int_t dataType = kMC, Bool_t recreate, Float_t energy, Float_t degraderThickness) {
   // Draw a histogram on the error of different track recon scenarios ...
   // Error ^ 2 = error_xy ^ 2 + ( 15 cm * error_rads ) ^ 2 
   // Check for different track recon efficiencies, and find RMS values :)
   
   run_energy = energy;
   run_degraderThickness = degraderThickness;

   Track *  thisTrack = nullptr;
   Track *  originalTrack = nullptr;
   Int_t    finalEID;
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);

   printf("Create EID index list ... ");
   tracks->createEIDSortList();
   printf("OK!\n");

   TCanvas *c = new TCanvas("c", "Total error from tracking errors", 1200, 900);
   c->Divide(3,1,0,0);

   TH1F *hPosError = new TH1F("hPosError", "Positional Error;Error [mm];Entries", 0, 100);
   TH1F *hAngError = new TH1F("hAngError", "Angular Error;Error [mrad];Entries", 0, 100);
   TH1F *hTotError = new TH1F("hTotError", "Total Projected Error;Error [mm];Entries", 0, 100);
      
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;

      if (thisTrack->isFirstAndLastEventIDEqual()) continue;

      finalEid = thisTrack->Last()->getEventID();
      originalTrack = (Track*) tracks->getTrackWithEID(finalEid);

      Float_t x0, y0, x1, y1, x0p, y0p, x1p, y0p;
      Float_t theta0, theta1;
      Float_t scalar0, dot0, scalar1, dot1;

      x0 = thisTrack->getXmm(0);
      y0 = thisTrack->getYmm(0);
      x1 = thisTrack->getXmm(1);
      y1 = thisTrack->getYmm(1);

      x0p = originalTrack->getXmm(0);
      y0p = originalTrack->getYmm(0);
      x1p = originalTrack->getXmm(1);
      y1p = originalTrack->getYmm(1);

      scalar1 = sqrt(pow(x0, 2) + pow(y0, 2) + pow(-150, 2)) * sqrt(pow(x1, 2) + pow(y1, 2) + pow(0, 2)); // angle = 15 cm in front of first layer
      dot1 = x0 * x1 + y0 * y1;

      angle = fabs((acos(dot1 / scalar1) - acos(dot0 / scalar0) )* 1000;
*/

#endif
