#define Analysis_cxx

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
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
//#include "Classes/Track/conversionFunctions.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/DataInterface/DataInterface.h"
#include "HelperFunctions/Tools.h"
#include "HelperFunctions/getTracks.h"

using namespace std;
using namespace DTC;

void writeDataFrame(Int_t energy) {
   DataInterface *di = new DataInterface();
   di->writeDataFrame(energy);
}

void findMCSAngles(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness) {
   // Make one histogram for each layer (starting at layer 1)
   // For each track (found using event ID information), find the change in angle
   // at layer 1, 2, 3, 4, 5, ..., nLayers.
   // The angle at layer 1 is defined as 
   // THETA1 = atan2(sqrt((x2-x1)^2 + (y2-y1)^2), dz)/sqrt(2) - atan2(sqrt(x1-x0)^2 + (y1-y0)^2, dz)/sqrt(2).
   // Use MM units for all distances.

   Int_t layers = 150;
   Int_t lastActivatedLayer = 0;

   run_degraderThickness = degraderThickness;
   run_energy = energy;
   if (useDegrader) run_energy = getEnergyAtWEPL(energy, degraderThickness);

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, run_energy);

   vector<TH1F*> *hAngleDifference = new vector<TH1F*>;
   vector<TCanvas*> *cCanvases = new vector<TCanvas*>;
   hAngleDifference->reserve(layers);
   cCanvases->reserve(layers);

   for (Int_t layer=0; layer<layers; layer++) {
      hAngleDifference->push_back(new TH1F(Form("hAngleDifference_layer_%i",layer), Form("Angular spread in layer %d for %.0f mm absorbator;Angular spread [rad];Entries",layer, degraderThickness), 100, 0, 0.1));
   }

   Track *thisTrack = nullptr;
   Float_t angle, y2, y1, y0, x2, x1, x0;
   Float_t entering[3] = {};
   Float_t leaving[3] = {};
   Float_t dotproduct, scalarproduct;

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;

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

         hAngleDifference->at(layer)->Fill(angle);
      }
   }
   
   for (Int_t layer=0; layer<lastActivatedLayer; layer++) {
      cCanvases->push_back(new TCanvas(Form("canvas_%d", layer), Form("LAYER %d", layer), 1200, 900));
   }

   TH1F *h = nullptr;
   printf("Layer mean_angle sigma_angle\n");
   for (Int_t layer=0; layer<lastActivatedLayer; layer++) {
      cCanvases->at(layer)->cd();
      h = hAngleDifference->at(layer);
      h->SetFillColor(kBlue-4);
      h->Draw();
      TF1 *fit = new TF1("fit", "gaus");
      h->Fit(fit, "Q");

      printf("%d %.3f %.3f\n", layer, fit->GetParameter(1), fit->GetParameter(2));

      delete fit;
   }
}

void drawTrackAngleAtVaryingRunNumbers(Int_t dataType, Float_t energy, Float_t degraderThickness) {
   Int_t nRuns = 0;

   run_energy = energy;
   run_degraderThickness = degraderThickness;
   
   if (useDegrader) {
      run_energy = getEnergyAtWEPL(energy, degraderThickness);
   }

   for (Int_t i=1; i<30; i++) { // 1 -> 30
      nRuns = pow(2, 4 + 0.25 * i) + 0.5;

      kEventsPerRun = nRuns;
      Float_t factor = 2;

      Int_t totalNumberOfRuns = 2600 / kEventsPerRun;
      if (totalNumberOfRuns < 1) totalNumberOfRuns = 1;
      if (totalNumberOfRuns > 75) totalNumberOfRuns = 75;

      Tracks * tracks = loadOrCreateTracks(1, totalNumberOfRuns, dataType, energy);

      char * sDataType = getDataTypeChar(dataType);
      TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
      TH1F *hAngles = new TH1F("hAngles", Form("Proton angle plot with %d protons in frame (%s)", nRuns, sDataType), 500, 0, 30);
      TCanvas *c2 = new TCanvas("c2", "Number of correct proton tracks with depth", 1200, 800);
      TH1F *hCorrectTracks = new TH1F("hCorrectTracks", "Number of correct proton tracks with depth", 10, 0, 10);
      TH1F *normCorrectTracks = new TH1F("normCorrectTracks", "Normalisation histogram", 10,0,10);

      hAngles->SetXTitle("Protons angle from initial measurement to layer 2");
      hAngles->SetYTitle("Number of protons");
      hAngles->SetFillColor(kCyan-8);
      hAngles->SetLineColor(kBlack);
      gStyle->SetOptStat(0);

      hCorrectTracks->SetTitle(Form("Tracks with same eventID using %d protons/run", kEventsPerRun));
      hCorrectTracks->SetXTitle("Layer number");
      hCorrectTracks->SetYTitle("Number of protons");
      hCorrectTracks->SetFillColor(kBlue-7);
      hCorrectTracks->SetLineColor(kBlack);

      Track *thisTrack;
      Int_t EID, thisEID;
      Int_t nTotal = tracks->GetEntries();
      Int_t nTotal2 = 0;
      Int_t nFirstAndLast = 0;
      Int_t nFirstAndLastAllTracks = 0;
      Int_t nLastCloseToFirst = 0;
      Int_t nCorrect = 0;

      for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
         thisTrack = tracks->At(j);
         if (!thisTrack) continue;
         hAngles->Fill(thisTrack->getSlopeAngleChangeBetweenLayers(2));

         EID = thisTrack->getEventID(0);
         if (EID > -1) { hCorrectTracks->Fill(0); }
         normCorrectTracks->Fill(0);
         nCorrect += (int) thisTrack->isOneEventID();
         nFirstAndLast += (int) thisTrack->isFirstAndLastEventIDEqual();
         nLastCloseToFirst += (Int_t) tracks->isLastEventIDCloseToFirst(j);
         nFirstAndLastAllTracks += (int) (thisTrack->isFirstAndLastEventIDEqual() && tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer()) == 0);
      
         for (Int_t k=1; k<thisTrack->GetEntriesFast(); k++) {
            if (!thisTrack->At(k)) continue;
            normCorrectTracks->Fill(thisTrack->getLayer(k));
            thisEID = thisTrack->getEventID(k);
            if (thisEID == EID || EID == -1) {
               hCorrectTracks->Fill(thisTrack->getLayer(k));
            }
         }
      }

      Float_t ratioCorrect = (float) nCorrect / nTotal;
      Float_t ratioFirstAndLast = (float) nFirstAndLast / nTotal;
      Float_t ratioFirstAndLastAllTracks = (float) nFirstAndLastAllTracks / nTotal;
      Float_t ratioLastCloseToFirst = (float) nLastCloseToFirst / nTotal;

      hCorrectTracks->Divide(normCorrectTracks);
      
      Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;

      ofstream file2("OutputFiles/lastLayerCorrect_different_nRuns.csv", ofstream::out | ofstream::app);
      file2 << readoutAbsorber << " " << nRuns << " " << ratioCorrect << " " << ratioFirstAndLast << " " << ratioLastCloseToFirst << " " << ratioFirstAndLastAllTracks << endl;
      file2.close();

      c1->cd();
      hAngles->Draw();
      c1->SaveAs(Form("OutputFiles/figures/angles/angles_layer%.1f_with_nRuns-%d.png", factor, nRuns));
      c1->SaveAs(Form("OutputFiles/figures/angles/angles_layer%.1f_with_nRuns-%d.root", factor, nRuns));

      c2->cd();
      hCorrectTracks->Draw();

      c2->SaveAs(Form("OutputFiles/figures/angles/correctTracks_factor%.1f_nruns%d.png", factor, nRuns));

      Float_t rms = hAngles->GetRMS();
      Float_t mean = hAngles->GetMean();
      Int_t binmax = hAngles->FindLastBinAbove(1);
      Float_t maximum = hAngles->GetXaxis()->GetBinCenter(binmax);

      ofstream file("OutputFiles/angles_different_nRuns.csv", ofstream::out | ofstream::app);
      file << factor << ";" << nRuns << ";" << rms << ";" << mean << ";"
           << maximum << endl;

      file.close();

      delete tracks;
      delete hAngles;
      delete c1;
      delete c2;
   }
}

void getTrackStatistics(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Int_t epr) {
   run_energy = energy;
   kDataType = dataType;
   
   if (epr>0) {
      kEventsPerRun = epr;
   }

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
   
   Int_t nTracksToPlot = 25;
   Int_t nTracksToPlot1D = 5;

   char * sDataType = getDataTypeChar(dataType);

   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
   TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
   TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
   TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);
   TCanvas *c5 = new TCanvas("c5", "c5", 800, 600);
   TCanvas *c6 = new TCanvas("c6", "c6", 1000, 800);
   TCanvas *c7 = new TCanvas("c7", "c7", 1000, 800);
   TCanvas *c8 = new TCanvas("c8", "c8", 1000, 800);

   c5->Divide(2, 2, 0.01, 0.01, 0);
   c6->Divide(3, 3, 0.01, 0.01, 0);
   c7->Divide(nTracksToPlot1D, nTracksToPlot1D, 0.001, 0.001, 0);
   c8->Divide(3,3,0.01,0.01,0);

   TH1F *hTrackLengths = new TH1F("hTrackLengths", Form("Track Lengths (%s)", sDataType), 100, 0, 120);
   TH2F *hClusterSizeAlongTrack = new TH2F("hClusterSizeAlongTrack",
            Form("Cluster size along track length (%s)", sDataType), 1.5*nLayers, 0, 1.5*nLayers*dz, 50, 0, 50);
   TH1F *hStraightness = new TH1F("hStraightness", Form("Sinuosity plot (%s)", sDataType), 500, 1, 1.01);
   TH1F *hSlope = new TH1F("hSlope", Form("Proton angle plot (%s)", sDataType), 500, 0, 20);

   // Average cluster size
   vector<TH1F*> *hAvgCS = new vector<TH1F*>;
   hAvgCS->reserve(4);
   for (Int_t chip=0; chip<4; chip++) {
      hAvgCS->push_back(new TH1F(Form("hAvgCS_chip_%i",chip),
            Form("Average Cluster Size vs Track Length for chip %i (%s)",chip, sDataType), 50, 0, 50));
   }

   // Cluster size for individual layers
   vector<TH1F*> *hCSLayer = new vector<TH1F*>;
   hCSLayer->reserve(9);
   for (Int_t layer=0; layer<9; layer++) {
      hCSLayer->push_back(new TH1F(Form("hCSLayer_%i", layer),
            Form("Cluster size for layer %i (%s)", layer, sDataType), 50, 0, 50));
   }

   // Cluster size along track length for a single track
   vector<TH1F*> *hFollowTrack = new vector<TH1F*>;
   hFollowTrack->reserve(nTracksToPlot);
   for (Int_t track=0; track<nTracksToPlot; track++) {
      hFollowTrack->push_back(new TH1F(Form("hFollowTrack_%i", track),
            Form("Cluster size along track length for a single track (%s)", sDataType), 50, 0, 50));
      hFollowTrack->at(track)->SetXTitle("Track Length [mm]");
      hFollowTrack->at(track)->SetYTitle("Cluster size [# of pixels]");
   }

   // Proton angle distribution in layer
   vector<TH1F*> *hAngles = new vector<TH1F*>;
   hAngles->reserve(9);
   for (Int_t layer=0; layer<9; layer++) {
      hAngles->push_back(new TH1F(Form("hAngles_%i", layer),
            Form("Proton angle distribution in layer %i (%s)", layer, sDataType), 50, 0, 20));
      hAngles->at(layer)->SetXTitle("Proton angle [deg]");
      hAngles->at(layer)->SetYTitle("Cluster size [# of pixels]");
   }
   
   hTrackLengths->SetXTitle("Track length [mm]");
   hClusterSizeAlongTrack->SetXTitle("Track length [mm]");
   hClusterSizeAlongTrack->SetYTitle("Cluster size [# of pixels]");
   hStraightness->SetXTitle("Sinuosity parameter");
   hSlope->SetXTitle("Total track angle (degree)");

   Float_t trackLengthSoFar = 0;
   Int_t trackNum = 0;
   Int_t chip = 0; // the quadrant
   Bool_t cutTL = false;
   Bool_t cutCHIP = false;
   Int_t okTL = 0;
   Int_t okCHIP = 0;
   Float_t ang = 0;

   Track *thisTrack;
   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);

      Float_t TL = thisTrack->getTrackLengthmm();
      Int_t x0 = thisTrack->getX(0);
      Int_t y0 = thisTrack->getY(0);

      cutTL = (TL > kMinimumTracklength) ? true : false;

      chip = (x0 >= nx/2) + 2 * (y0 < ny/2);
      cutCHIP = (chip<2 || dataType == kMC) ? true : false;

      hTrackLengths->Fill(TL);
      hStraightness->Fill(thisTrack->getSinuosity());
      hSlope->Fill(thisTrack->getSlopeAngle());

      for (Int_t j=0; j<tracks->GetEntriesFast(i); j++) {

         trackLengthSoFar += thisTrack->getTrackLengthmmAt(j);
         if (cutTL && cutCHIP)
            hClusterSizeAlongTrack->Fill(trackLengthSoFar, thisTrack->getSize(j));

         if (cutTL)
            hAvgCS->at(chip)->Fill(trackLengthSoFar, thisTrack->getSize(j));

         Int_t layer = thisTrack->getLayer(j);
         if (layer<9) {
            hCSLayer->at(layer)->Fill(thisTrack->getSize(j));
//          hAngles->at(layer)->Fill(thisTrack->getSlopeAngleAtLayer(j));
            ang = thisTrack->getSlopeAngleChangeBetweenLayers(j);
            if (ang>=0) {
               hAngles->at(layer)->Fill(ang);
            }
         }

         if (trackNum < nTracksToPlot && cutTL && cutCHIP) {
            hFollowTrack->at(trackNum)->Fill(trackLengthSoFar, thisTrack->getSize(j));
            hFollowTrack->at(trackNum)->SetTitle(Form("Track length histogram for run %i (%s)", i, sDataType));
         }
      }
      trackLengthSoFar = 0;

      if (cutTL) okTL++;
      if (cutCHIP) okCHIP++;
      if (cutTL && cutCHIP) trackNum++;
   }

   cout << "Total number of tracks: " << tracks->GetEntriesFast() << endl;
   cout << "Passed track length: " << okTL << " (" << 100 * okTL / tracks->GetEntriesFast() << ")\n";
   cout << "Passed chip #: " << okCHIP << " (" << 100 * okCHIP / tracks->GetEntriesFast() << ")\n";

   c1->cd();
      hTrackLengths->Draw();
   c2->cd();
      gStyle->SetOptStat(0);
      hClusterSizeAlongTrack->Draw("COLZ");
   c3->cd();
      hStraightness->Draw();
   c4->cd();
      hSlope->Draw();

   for (Int_t chip=0; chip<4; chip++) {
      c5->cd(chip+1);
      hAvgCS->at(chip)->SetFillColor(kRed-chip*2);
      hAvgCS->at(chip)->Draw();
   }
   for (Int_t layer=0; layer<9; layer++) {
      c6->cd(layer+1);
      hCSLayer->at(layer)->SetFillColor(kRed-9+layer);
      hCSLayer->at(layer)->Draw("same");

//    cout << "Average cluster size and RMS for layer " << layer << " is \033[1m " << hCSLayer->at(layer)->GetMean() << " pixels \033[0m and \033[1m " << hCSLayer->at(layer)->GetRMS() << "pixels \033[0m\n";
   }

   for (Int_t track=0; track<nTracksToPlot; track++) {
      c7->cd(track+1);
      gPad->DrawFrame(0, 0, 50, 35);
      hFollowTrack->at(track)->SetFillColor(kBlue-2);
      hFollowTrack->at(track)->Draw("same");
   }

   fillMCSRadiusList(1);
   for (Int_t layer=0; layer<9; layer++) {
      c8->cd(layer+1);


      Float_t meanAngleAtLayer = findMCSAtLayerRad(layer, run_energy);
      Float_t mcs = getMCSAngleForLayer(layer) / cos(meanAngleAtLayer);

//    cout << "The added MCS factor due to inclined crossing is " << 1/(cos(meanAngleAtLayer)) << ".\n";

      TF1 *mcsGauss = new TF1("mcsGauss", "gaus(0)", 0, 25);
      if (hAngles->at(layer)->Integral()>0) {
         mcsGauss->SetParameters(10, 0, mcs);
         mcsGauss->SetParLimits(0, 1, 1000);
         mcsGauss->SetParLimits(1, 0, 0);
         mcsGauss->SetParLimits(2, mcs, mcs);

         hAngles->at(layer)->Fit("mcsGauss", "M,B,Q");
      }

      int maxHeight = hAngles->at(layer)->GetMaximum();
      TLine *line = new TLine(mcs, 0, mcs, maxHeight * 1.05);
      TLine *line2 = new TLine(mcs  * 2, 0, mcs * 2, maxHeight * 1.05);
      TLine *line3 = new TLine(mcs  * 3, 0, mcs * 3, maxHeight * 1.05);

      hAngles->at(layer)->SetFillColor(kRed-9+layer);
      hAngles->at(layer)->Draw("same");
      mcsGauss->Draw("same");
      line->Draw("same");
      line2->Draw("same");
      line3->Draw("same");
   }

   delete tracks;
}
//
void drawClusterShapes(Int_t Runs, Bool_t dataType, Bool_t recreate, Float_t energy) {
   // get vector of TH2F's, each with a hits distribution and cluster size
   // dataType = kMC (0) or kData (1)

   run_energy = energy;
   kDataType = dataType;

   Int_t nRows = 15;
   Int_t nRepeats = 20;
   Int_t nN = nRows * nRepeats;
   vector<TH2C*> *hClusterMaps = new vector<TH2C*>;
   hClusterMaps->reserve(nN);
   for (Int_t i=0; i<nN; i++)
      hClusterMaps->push_back(new TH2C(Form("hClusterMap_%i",i), "", 11, 0, 11, 11, 0, 11));

   Int_t nClusters = kEventsPerRun * 5;
   Int_t nHits = kEventsPerRun * 50;
   Int_t nTracks = kEventsPerRun * 2;

   DataInterface *di = new DataInterface();
   CalorimeterFrame *cf = new CalorimeterFrame();
   Hits *hits = new Hits(nHits);
   vector<Hits*> * tempClusterHitMap;
   vector<Hits*> * clusterHitMap = new vector<Hits*>;
   clusterHitMap->reserve(Runs*500*kEventsPerRun);

   for (Int_t i=0; i<Runs; i++) {
      if (dataType == kMC) {
         di->getMCFrame(i, cf);
         cf->diffuseFrame(new TRandom3(0)); // THE MAGIC PART
         hits = cf->findHits();
         tempClusterHitMap = hits->findClustersHitMap();
      }

      else if (dataType == kData) {
         di->getDataFrame(i, cf, energy);
         hits = cf->findHits();
         tempClusterHitMap = hits->findClustersHitMap();
      }
      else {
        std::cerr << "Please choose between dataType = kMC (0) or kData (1).\n" << endl;
        exit(1);
      }
   
      for (UInt_t j=0; j<tempClusterHitMap->size(); j++)
         clusterHitMap->push_back( tempClusterHitMap->at(j) );
   }
   // delete tempClusterHitMap;
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
   for (Int_t i=0; i<nRows; i++) {
      size_from = i*3;
      size_to = size_from + 2;
      nInRow = 0;

      for (UInt_t j=0; j<clusterHitMap->size(); j++) {
         Int_t csize = clusterHitMap->at(j)->GetEntriesFast();
         if (csize >= size_from && csize <= size_to) {
            // plot the cluster in the j'th row
            for (Int_t k=0; k<clusterHitMap->at(j)->GetEntriesFast(); k++) {
               x = clusterHitMap->at(j)->getX(k);
               y = clusterHitMap->at(j)->getY(k);
               hClusterMaps->at(nTotal)->Fill(x,y);
            } // end loop through all points
            hClusterMaps->at(nTotal)->SetTitle(Form("Cluster size [%i, %i]", size_from, size_to));
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
   c->Divide(nRepeats, nRows, 0.0001, 0.0001);
   for (Int_t i=0; i<nN; i++) {
      c->cd(i+1);
      hClusterMaps->at(i)->Draw("same, COL,ah,fb,bb");
      gStyle->SetOptStat(0);
   }
}

// a function of the same signature is also defined in HelperFunctions/getTracks.C
// in the end it looks like it is not used anyway
/*
void makeTracks(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
}
*/

void drawTrackRanges(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   kDataType = dataType;
   
   char * sDataType = getDataTypeChar(dataType);
   char * sMaterial = getMaterialChar();
   char * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", energy, sMaterial, sDataType);
   TCanvas *cFitResults = new TCanvas("cFitResults", hTitle, 1400, 1000);
   
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetTitleH(0.06);
   gStyle->SetTitleYOffset(1);
   
   TGraphErrors *outputGraph;
   TH1F *hFitResults = new TH1F("fitResult", hTitle, 500, getWEPLFromEnergy(0), getWEPLFromEnergy(energy*1.2));
   hFitResults->SetLineColor(kBlack);
   hFitResults->SetFillColor(kGreen-5);
   hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
   hFitResults->SetYTitle("Number of protons");
   hFitResults->Draw();

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);

   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;
      Float_t preEnergyLoss = thisTrack->getPreEnergyLoss();

      hFitResults->Fill(getWEPLFromEnergy(thisTrack->getEnergy() + preEnergyLoss)); 
   }
}

void drawTungstenSpectrum(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   kDataType = dataType;
   
   char * sDataType = getDataTypeChar(dataType);
   char * sMaterial = getMaterialChar();
   char * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", energy, sMaterial, sDataType);
   TCanvas *cFitResults = new TCanvas("cFitResults", hTitle, 1400, 1000);
   
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetTitleH(0.06);
   gStyle->SetTitleYOffset(1);
   
   
   TGraphErrors *outputGraph;
   TH1F *hFitResults = new TH1F("fitResult", hTitle, 500, getWEPLFromEnergy(0), getWEPLFromEnergy(energy*1.2));
   hFitResults->SetLineColor(kBlack);
   hFitResults->SetFillColor(kGreen-5);
   hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
   hFitResults->SetYTitle("Number of protons");
   hFitResults->Draw();

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);

   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;
      
//       if (thisTrack->doesTrackEndAbruptly()) {
//          hFitResults->Fill(getWEPLFromEnergy(thisTrack->getEnergy()));
//          continue;
//       }
      
      outputGraph = (TGraphErrors*) thisTrack->doRangeFit();
      if (!outputGraph) continue;
      
      hFitResults->Fill(getUnitFromWEPL(thisTrack->getFitParameterRange()));
   }
}

void drawScintillatorStatistics(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   kDataType = dataType;
   run_energy = energy;
   
   char * sDataType = getDataTypeChar(dataType);
   char * sMaterial = getMaterialChar();
   char * hTitle = Form("Number of scintillators hit with a %.2f MeV beam in %s (%s)", energy, sMaterial, sDataType);
   TCanvas *cFitResults = new TCanvas("cFitResults", hTitle, 1400, 1000);
   
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetTitleH(0.06);
   gStyle->SetTitleYOffset(1);
   
   TGraphErrors *outputGraph;
   TH1I *hFitResults = new TH1I("fitResult", hTitle, 5, 0, 4);
   hFitResults->SetLineColor(kBlack);
   hFitResults->SetFillColor(kGreen-5);
   hFitResults->SetXTitle("Number of scintillators hit");
   hFitResults->SetYTitle("Number of protons");
   hFitResults->Draw();

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);

   Float_t nScintillators = 0;
   
   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;

      nScintillators = thisTrack->getNScintillators();   
      hFitResults->Fill(nScintillators);  
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

      outputGraph = (TGraphErrors*) thisTrack->doRangeFit(true);
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

void drawTracksWithEnergyLoss(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   kDataType = dataType;
   
   Float_t x_energy[sizeOfEventID*400] = {};
   Float_t y_energy[sizeOfEventID*400] = {};
   
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy, x_energy, y_energy);
   
   cout << "Using aluminum plate: " << kIsAluminumPlate << endl;
   cout << "Using scintillators: " << kIsScintillator << endl;

   Int_t nPlotX = 2, nPlotY = 2;
   Int_t fitIdx = 0, plotSize = nPlotX*nPlotY;
   Int_t skipIdx = 0;

   TGraphErrors *outputGraph;
   TCanvas *cAccuracy = new TCanvas("cAccuracy", "Accuracy of fit", 1400, 1000);
   TCanvas *cGraph = new TCanvas("cGraph", "Fitted data points", 1400, 1000);
   TH2F *hAccuracy = new TH2F("hAccuracy", Form("Accuracy of fit in a %.0f MeV nominal beam", run_energy), 400, 0, 180, 400, 0, 180);
   cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
   gStyle->SetPadBorderMode(0); gStyle->SetFrameBorderMode(0);
   gStyle->SetTitleH(0.06); gStyle->SetTitleYOffset(1);
   hAccuracy->SetYTitle("Fit energy [MeV]"); hAccuracy->SetXTitle("Real energy [MeV]");

   Float_t finalEnergy, realEnergy, preTL = 0;
   Float_t fitEnergy, fitRange, fitScale, fitError = 0;
   Int_t eventID = -1;
   
   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;

      outputGraph = (TGraphErrors*) thisTrack->doRangeFit();
      if (!outputGraph) {
         continue;
      }

      eventID = thisTrack->At(0)->getEventID()-1;

      preTL = thisTrack->getPreTL() - 1.5;
//    if       (thisTrack->getYmm(0) > 0)    { preTL = thisTrack->getPreTL() + firstUpperLayerZ; }
//    else if  (thisTrack->getYmm(0) < 0)    { preTL = thisTrack->getPreTL() + firstLowerLayerZ; }
      for (Int_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) { x_energy[i] += preTL; }

      convertXYToWEPL(x_energy, y_energy, eventID);
      realEnergy = getEnergyFromXY(x_energy, y_energy, eventID);

      fitRange = thisTrack->getFitParameterRange();
      fitEnergy = getEnergyFromWEPL(fitRange);
      fitScale  = thisTrack->getFitParameterScale();
      fitError = quadratureAdd(thisTrack->getFitParameterError(), dz/sqrt(12));
      
      hAccuracy->Fill(realEnergy, fitEnergy);

      if (fitIdx < plotSize) {
         drawIndividualGraphs(cGraph, outputGraph, fitRange, fitScale, fitError, fitIdx++, eventID, x_energy, y_energy);
      }

      else delete outputGraph;
   }
   cAccuracy->cd();
   hAccuracy->Draw("COLZ");
}

Float_t drawTH2FRangeAccuracy() {
   TCanvas      * c = new TCanvas("c", "Range determination accuracy", 1200, 1200);
   TH2F         * hRangeAccuracy = new TH2F("hRangeAccuracy", "Range Determination Accuracy;Nominal range [WEPL mm];Calculated range [WEPL mm]", 400, 0, 400, 800, 0, 400);
   Float_t        nominalRange, calculatedRange;
   TGraphErrors * outputGraph = nullptr;
   Track        * thisTrack = nullptr;

   for (Int_t degraderThickness = 20; degraderThickness < 350; degraderThickness += 2) {
      run_degraderThickness = degraderThickness;
      run_energy = getEnergyAtWEPL(250, degraderThickness);
      nominalRange = getUnitFromEnergy(run_energy);

      Tracks * tracks = loadOrCreateTracks(1, 1, 0, run_energy);

      for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
         thisTrack = tracks->At(i);

         if (!thisTrack) continue;
         if (thisTrack->doesTrackEndAbruptly()) continue;

         outputGraph = (TGraphErrors*) thisTrack->doRangeFit(false);
         hRangeAccuracy->Fill(nominalRange, thisTrack->getFitParameterRange());
     }

      delete tracks;
   }
   
   gStyle->SetOptStat(0);
   gPad->SetLogz();
   hRangeAccuracy->Draw("COLZ");

   TLine *line = new TLine(0, 0, 400, 400);
   line->Draw();
}

Float_t drawBraggPeakGraphFit(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness) {
   run_degraderThickness = degraderThickness;
   run_energy = energy;

   printf("WEPL factors for different energies... \n");
   for (int i=50; i<260; i += 50) {
      printf("%d MeV -> factor = %.3f\n", i, getWEPLFactorFromEnergy(i));
   }

   if (useDegrader) {
//      run_energy = getEnergyAtWEPL(energy, degraderThickness);
      run_energy = getEnergyFromDegraderThickness(degraderThickness);
   }

   printf("Using water degrader of thickness %.0f mm, the initial energy of %.0f MeV is reduced to %.1f MeV.\n", degraderThickness, energy, run_energy);

   kDataType = dataType;
   Bool_t         kDrawHorizontalLines = false;
   Bool_t         kDrawVerticalLayerLines = false;
   Bool_t         kDrawIndividualGraphs = true;
   Bool_t         acceptAngle = true;
   Float_t        maxAngle, thisAngle;
   Float_t        cutAngle = (run_degraderThickness * 0.0219 + 0.556) * 0.8;
   Float_t        finalEnergy = 0;
   Float_t        fitRange, fitScale, fitError;
   Int_t          nCutDueToTrackEndingAbruptly = 0;
   Int_t          nPlotX = 3, nPlotY = 3;
   Int_t          fitIdx = 0, plotSize = nPlotX*nPlotY;
   Int_t          skipPlot = 10;
   TGraphErrors * outputGraph;
   char         * sDataType = getDataTypeChar(dataType);
   char         * sMaterial = getMaterialChar();
   char         * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", run_energy, sMaterial, sDataType);

   Int_t nEnergies = getUnitFromEnergy(run_energy);
   gStyle->SetOptTitle(0);

   if (useDegrader) {
      hTitle = Form("Fitted energy of a %.0f MeV nominal beam on %s DTC w/%.1f mm water degrader", energy, sMaterial, degraderThickness);
   }
   TCanvas      * cGraph = new TCanvas("cGraph", "Fitted data points", nPlotX*500, nPlotY*500);
   TCanvas      * cFitResults = new TCanvas("cFitResults", hTitle, 1000, 1000);
//   TCanvas      * cAngle = new TCanvas("cAngle", "Incoming angles", 800, 800);
   TH1F         * hFitResults = new TH1F("fitResult", hTitle, fmax(nEnergies*8,200), getUnitFromEnergy(0), getUnitFromEnergy(run_energy)*1.4+10);
   TH1F         * hFitResultsDroppedData = new TH1F("fitResultsDroppedData", hTitle, 200, getUnitFromEnergy(0), getUnitFromEnergy(run_energy)*1.4+10);
   TH1F         * hMaxAngle = new TH1F("hMaxAngle", "Maximum angle for proton track", 200, 0, 25);
 
   printf("Histogram limits: %.2f to %.2f.\n", getUnitFromEnergy(0), getUnitFromEnergy(run_energy)*1.2);

   printf("At energy %.0f, expecting range %.2f mm and WEPL %.2f mm.\n", run_energy, getTLFromEnergy(run_energy), getWEPLFromEnergy(run_energy));
   printf("This corresponds to a WEPL factor of %.2f.\n", getWEPLFactorFromEnergy(run_energy));
   cout << "Correcting for aluminum plate: " << kIsAluminumPlate << endl;
   cout << "Correcting for scintillators: " << kIsScintillator << endl;

   // Histogram options
   cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
   gStyle->SetPadBorderMode(0); gStyle->SetFrameBorderMode(0);
   gStyle->SetTitleH(0.06); gStyle->SetTitleYOffset(1);
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   gStyle->SetPadTopMargin(0.05); gStyle->SetPadRightMargin(0.05);
   gStyle->SetPadBottomMargin(0.05);
   gStyle->SetPadLeftMargin(0.15);
   hFitResults->SetLineColor(kBlack); hFitResults->SetFillColor(kGreen-5);
   hFitResultsDroppedData->SetLineColor(kBlack); hFitResultsDroppedData->SetFillColor(kYellow-2);
   hMaxAngle->SetLineColor(kBlack); hMaxAngle->SetFillColor(kGreen-5);

   // Create or load all tracks
   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, run_energy);
   
   for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
      Track *thisTrack = tracks->At(j);
      if (!thisTrack) continue;
    
      if (thisTrack->doesTrackEndAbruptly()) {
         hFitResultsDroppedData->Fill(getUnitFromEnergy(thisTrack->getEnergy()));
         nCutDueToTrackEndingAbruptly++;
         continue;
      }
   
      maxAngle = thisTrack->getSlopeAngleAtLayer(1);
      if (maxAngle > cutAngle && acceptAngle) continue;
      hMaxAngle->Fill(maxAngle);

      // Do track fit, extract all parameters for this track
      outputGraph = (TGraphErrors*) thisTrack->doRangeFit(false);
      if (!outputGraph) continue;

      fitRange = thisTrack->getFitParameterRange();
      fitScale = thisTrack->getFitParameterScale();
      fitError = quadratureAdd(thisTrack->getFitParameterError(), dz/sqrt(12)); // latter term from error on layer position

      cGraph->Update();

      hFitResults->Fill(getUnitFromTL(fitRange));

      if (fitIdx < plotSize && kDrawIndividualGraphs && j>=skipPlot) {
         drawIndividualGraphs(cGraph, outputGraph, fitRange, fitScale, fitError, fitIdx++);

         // Drawing of three panels with individual range-Edep graphs
         if (fitIdx == 1) {
            gPad->SetRightMargin(0.01);
            gPad->SetLeftMargin(0.15);
//            outputGraph->GetXaxis()->SetTitle("");
         }

         else if (fitIdx == 2) {
            gPad->SetLeftMargin(0.01);
            gPad->SetRightMargin(0.01);
            outputGraph->GetYaxis()->SetLabelOffset(2);
            outputGraph->SetTitle(Form("Bragg-Kleeman fit to exp. data at %.0f MeV", run_energy));
           /* 
            gPad->Update();
            TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
            title->SetTextFont(22);
            title->SetTextSize(0.06);
            gPad->Modified();
            */
            
         }

         else if (fitIdx == 3) {
            gPad->SetRightMargin(0.1);
            gPad->SetLeftMargin(0.01);
            outputGraph->GetXaxis()->SetTitle("");
            outputGraph->GetYaxis()->SetLabelOffset(2);
         }
      }
      
      else delete outputGraph;

   }
   
   if (!kDrawIndividualGraphs) delete cGraph;
   
   cout << 100 * float(nCutDueToTrackEndingAbruptly) / tracks->GetEntriesFast() << " % of the tracks were cut due to seemingly inelastic nuclear interactions.\n";

   TF1 *fMaxAngle = new TF1("fMaxAngle", "gaus(0)", 0, 25);
   fMaxAngle->SetParameters(100, 4, 6);
   hMaxAngle->Fit(fMaxAngle, "M, W, Q, N", "", 0, 25);

   Float_t angleTo = fMaxAngle->GetParameter(1) + 3 * fMaxAngle->GetParameter(2);

   cout << "3 sigma Confidence Limit for angular spread  = " << angleTo << endl;

   Int_t nAccepted = hMaxAngle->Integral(0,hMaxAngle->GetXaxis()->FindBin(angleTo));
   Float_t percentAccepted = 100 * nAccepted / hMaxAngle->Integral(0);

   cout << "Number of accepted events = " << nAccepted << " of total " << hMaxAngle->Integral()  << "(" <<  percentAccepted << " %) " << endl;

   cFitResults->cd();

   if       (kOutputUnit == kPhysical) hFitResultsDroppedData->SetXTitle("Physical range [mm]");
   else if  (kOutputUnit == kWEPL)     hFitResultsDroppedData->SetXTitle("Range in Water Equivalent Path Length [mm]");
   else if  (kOutputUnit == kEnergy)   hFitResultsDroppedData->SetXTitle("Energy [MeV]");

   hFitResultsDroppedData->SetYTitle("Number of protons");
   hFitResultsDroppedData->GetXaxis()->SetTitleFont(22);
   hFitResultsDroppedData->GetXaxis()->SetLabelFont(22);
   hFitResultsDroppedData->GetYaxis()->SetTitleFont(22);
   hFitResultsDroppedData->GetYaxis()->SetLabelFont(22);
   hFitResultsDroppedData->GetYaxis()->SetTitleOffset(1.5);
   
   if       (kOutputUnit == kPhysical) hFitResults->SetXTitle("Physical range [mm]");
   else if  (kOutputUnit == kWEPL)     hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
   else if  (kOutputUnit == kEnergy)   hFitResults->SetXTitle("Energy [MeV]");

   hFitResults->SetYTitle("Number of protons");
   hFitResults->GetXaxis()->SetTitleFont(22);
   hFitResults->GetXaxis()->SetLabelFont(22);
   hFitResults->GetYaxis()->SetTitleFont(22);
   hFitResults->GetYaxis()->SetLabelFont(22);
   hFitResults->GetYaxis()->SetTitleOffset(1.5);

   hFitResults->SetTitle("");
   
//   hFitResultsDroppedData->Draw();
   hFitResults->Draw(); // was same
//   hFitResultsDroppedData->GetYaxis()->SetRangeUser(0, max(hFitResultsDroppedData->GetMaximum()*1.05, hFitResults->GetMaximum()*1.05));

   /* 
    * Uncomment when finished making plots for optimization poster
    *
   gPad->Update();
   TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();
   */

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

   TF1 *gauss = doSimpleGaussianFit(hFitResults, means, sigmas);
   Float_t empiricalMean = means[9];
   Float_t empiricalSigma = sigmas[9];
   
   Float_t energySigma = getEnergyFromUnit(empiricalMean +  empiricalSigma/ 2) - getEnergyFromUnit(empiricalMean - empiricalSigma / 2);

   cFitResults->Update();

   gPad->Update();

   TLine *l = nullptr;
   if (kDrawVerticalLayerLines) {
      Float_t line_z = 0;
      for (Int_t i=0; i<65; i++) {
         line_z = getUnitFromTL(getLayerPositionmm(i));
         if (line_z > gPad->GetUxmax()) break;
         l = new TLine(line_z, 0, line_z, hFitResults->GetMaximum()*1.05);
         l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw();
      }
   }

   gPad->Update();
   Float_t bip_value = empiricalMean - 6 * empiricalSigma;
   if (bip_value == 0) bip_value = getWEPLFromEnergy(run_energy)*0.9;
   TLine *bip = new TLine(bip_value, gPad->GetUymax(), bip_value, 0);
   bip->Draw();

   Float_t bip_value2 = empiricalMean + 6 * empiricalSigma;
   TLine *bip2 = new TLine(bip_value2, 0, bip_value2, gPad->GetUymax());
   bip2->Draw();

   TLegend *legend = new TLegend(0.16, 0.78, 0.42, 0.88);
   legend->SetTextSize(0.028);
   legend->AddEntry(hFitResults, "Accepted tracks", "F");
   legend->AddEntry(hFitResultsDroppedData, "Tracks without BP rise", "F");
   legend->AddEntry(gauss, "Fitted Gaussians", "L");
   legend->AddEntry(bip, "x_{i'} bin", "L");
   legend->SetTextFont(22);
   if (kDrawVerticalLayerLines) legend->AddEntry(l, "Sensor layer positions", "L");
//   legend->Draw();

   gStyle->SetOptStat(11);
   TPaveStats *ps = (TPaveStats*) cFitResults->GetPrimitive("stats");
   hFitResultsDroppedData->SetBit(TH1::kNoStats);
   hFitResults->SetBit(TH1::kNoStats);
   ps->SetX1NDC(0.4176); ps->SetX2NDC(0.9257);
   ps->SetY1NDC(0.7415); ps->SetY2NDC(0.9712);
   ps->SetTextFont(22);
   ps->AddText(Form("Nominal WEPL = %.2f #pm %.2f", expectedMean, expectedStraggling));
//   ps->AddText(Form("Nominal straggling = %.2f", expectedStraggling));
  
   /*
   for (Int_t i=0; i<2; i++) {
      ps->AddText(Form("Fit %d WEPL = %.2f", i+1, means[i]));
      ps->AddText(Form("Fit %d #sigma_{WEPL} = %.2f", i+1, sigmas[i]));
      ps->AddText(Form("Fit %d energy = %.2f", i+1, getEnergyFromUnit(means[i])));
   }
   */

   ps->AddText(Form("Calculated WEPL = %.2f #pm %.2f", empiricalMean, empiricalSigma));
   ps->AddText(Form("WEPL deviation = %.2f #pm %.2f", empiricalMean - expectedMean, empiricalSigma - expectedStraggling));
//   ps->AddText(Form("Resulting energy = %.2f #pm %.2f", getEnergyFromUnit(empiricalMean), energySigma));
      
/*
   cAngle->cd();
   hMaxAngle->SetTitle("Proton angles;Incoming angle [AU];Entries");
   hMaxAngle->Draw();
*/
   
   if (kOutputUnit == kPhysical) {
      cFitResults->SaveAs(Form("OutputFiles/figures/Fitted_energies_%.2f_MeV_%s_%s_physical.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
   }
   else if (kOutputUnit == kEnergy) {
      cFitResults->SaveAs(Form("OutputFiles/figures/Fitted_energies_%.2f_MeV_%s_%s_energy.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
   }
   
   else if (kOutputUnit == kWEPL) {
      cFitResults->SaveAs(Form("OutputFiles/figures/Fitted_energies_%.2f_MeV_%s_%s_WEPL.png", energy, getMaterialChar(), getDataTypeChar(dataType)));
   }

   delete tracks;
   
   return empiricalMean;
}

void writeClusterFile(Int_t Runs, Int_t dataType, Float_t energy) {
   run_energy = energy;
   kDataType = dataType;
   
   Int_t nClusters = kEventsPerRun * 5 * nLayers;
   Int_t nHits = kEventsPerRun * 50;
   Bool_t kRemoveSmallClusters = true;
   
   DataInterface *di = new DataInterface();
   CalorimeterFrame *cf = new CalorimeterFrame();
   Hits * hits = new Hits(nHits);
   Clusters * clusters = new Clusters(nClusters);
   
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
         clusters = hits->findClustersFromHits();
         
         if (kRemoveSmallClusters) {
            Int_t maxRemoveSize = 2;
            clusters->removeSmallClusters(maxRemoveSize);
         }
      }
   }

   delete di;

   ofstream file("OutputFiles/output_all_layers.csv");
   file << "layer;x;y;clustersize" << endl;

   for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
      file << clusters->getLayer(i) << ";" 
           << clusters->getX(i) << ";"
          << clusters->getY(i) << ";" 
          << clusters->getSize(i) << endl;
   }
   
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

   /*
   gPad->Update();
   TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();
*/
   gStyle->SetOptStat(0);
   
   c1->cd(1);
   normalizeFrame->Draw("colz");
   c1->cd(2);
   normalizeFrameLast->Draw("colz");

   c2->cd(1);
   beamSpec->Draw();

   c2->cd(2);
   beamSpecLast->Draw();

   /*
   // draw lines from 474 to 806
   TLine *l1 = new TLine(474, 0, 474, 1280);
   TLine *l2 = new TLine(806, 0, 806, 1280);
   TLine *l3 = new TLine(0, 474, 1280, 474);
   TLine *l4 = new TLine(0, 806, 1280, 806);
   l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
   */
}

Hits * getEventIDs(Int_t Runs, Float_t energy) {
   run_energy = energy;
   DataInterface *di = new DataInterface();

   Hits * hits = new Hits(kEventsPerRun * sizeOfEventID * Runs);

   for (Int_t i=0; i<Runs; i++) {
      di->getEventIDs(i, hits);
   }

   return hits;
}

void drawClusterSizeDistribution(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   DataInterface *di = new DataInterface();

   Int_t nClusters = kEventsPerRun * 5 * nLayers;
   Int_t nHits = kEventsPerRun * 50;

   CalorimeterFrame *cf = new CalorimeterFrame();
   Clusters * clusters = new Clusters(nClusters);
   Hits * hits = new Hits(nHits);
   TH2F * hClusterSizes = new TH2F("hClusterSizes", "Cluster sizes vs layer", nLayers, 0, nLayers-1, 50, 0, 50);

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
         clusters = hits->findClustersFromHits();
         clusters->removeSmallClusters(2);
         clusters->removeAllClustersAfterLayer(8);
      }
      
      clusters->Compress();

      for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
         hClusterSizes->Fill(clusters->getLayer(i), clusters->getSize(i));
      }
   }

   hClusterSizes->Draw("COLZ");
}
      
void drawTracks3D(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t switchLayer, Float_t energy, Float_t degraderThickness) {
   run_energy = energy;
   run_degraderThickness = degraderThickness;

   Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);

   TCanvas *c1 = new TCanvas("c1", "c1", 1600, 800);
   c1->SetTitle(Form("Tracks from %.2f MeV protons on %s", energy, getMaterialChar()));
   TView *view = TView::CreateView(1);
   float fromx = 0.1 * nx;
   float tox = 0.9 * nx;
   float fromy = 0.1 * ny;
   float toy = 0.9 * ny;

   view->SetRange(fromx, 0, fromy, tox, 70, toy);
   Int_t iret;
   Float_t theta = 285;
   Float_t phi = 80;

   view->SetView(theta, phi, 0, iret);

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
   pMarker->SetMarkerStyle(15);
   EIDMarker->SetMarkerColor(kRed);
   conflictMarker->SetMarkerColor(kRed); // Conflicting cluster
   
   for (Int_t i=0; i<restPoints->GetEntriesFast(); i++) {
      if (!restPoints->At(i))
         continue;

      Cluster *thisCluster = (Cluster*) restPoints->At(i);
      Float_t x = thisCluster->getX();
      Float_t z = thisCluster->getY();
      Float_t y = thisCluster->getLayer();

      pMarker->SetPoint(i, x, y, z);
   }

   pMarker->Draw();
   
   Int_t ntracks = tracks->GetEntriesFast();
   Int_t firstEID;
   Int_t EIDidx = 0;
   Int_t conflictIdx = 0;

   Int_t medianEventID = -1;

   Int_t nTrueTracks = 0;
   Int_t nOKTracks = 0;
   Int_t nOKTracksAllClusters = 0;
   Int_t nOKMinusTracks = 0;
   Int_t nOKLastLayers = 0;
   Int_t nOneWrong = 0;
   Int_t nMissingEID;

   Track * thisTrack = nullptr;

   for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
      thisTrack = tracks->At(i);
      if (!thisTrack) continue;

      if (thisTrack->isOneEventID()) nTrueTracks++;

      if (!thisTrack->isOneEventID()) {
         medianEventID = thisTrack->getModeEventID();
      }
      
      nMissingEID = tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer()); // only count missing tracks after layer
      if (thisTrack->isFirstAndLastEventIDEqual() && nMissingEID == 0) nOKTracksAllClusters++;

      if (thisTrack->isFirstAndLastEventIDEqual()) nOKTracks++;

      else {
         if (!thisTrack->Last()) continue;
         if (!thisTrack->At(thisTrack->GetEntriesFast() - 2)) continue;

         Int_t lastEID = thisTrack->Last()->getEventID();   
         if (lastEID < 0) continue;

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
         else if (thisTrack->getWEPL() < 0.6 * getWEPLFromEnergy(run_energy)) {
            // Bad track ends early. OK...
            nOKLastLayers++;
         }
      }
   }

   nOKMinusTracks += nOKTracks;
   nOKLastLayers += nOKTracks;

   Int_t numberOfTracks = tracks->GetEntries();

   Float_t factorEID = 100 * ((float) nTrueTracks / numberOfTracks);
   Float_t factorEIDOK = 100 * ((float) nOKTracks / numberOfTracks);
   Float_t factorEIDOKAllClusters = 100 * ((float) nOKTracksAllClusters / numberOfTracks);
   Float_t factorEIDOKMinus = 100 * ((float) nOKMinusTracks / numberOfTracks);
   Float_t factorLastLayers = 100 * ((float) nOKLastLayers / numberOfTracks);

   cout << nTrueTracks << " of total " << numberOfTracks << " tracks has the same event ID (" << factorEID << "%)\n";
   cout << nOKTracks << " of total " << numberOfTracks << " tracks has the same first/last event ID (" << factorEIDOK << "%)\n";
   cout << nOKTracksAllClusters << " of total " << numberOfTracks << " tracks has the same first/last event ID and NO missing clusters in track (" << factorEIDOKAllClusters << "%)\n";
   cout << nOKMinusTracks << " of total " << numberOfTracks << " tracks has a close match (0.5 mm, 1 degree) on first / last cluster (" << factorEIDOKMinus << "%)\n";
   cout << nOKLastLayers << " of total " << numberOfTracks << " tracks has a close match (0.5 mm, 1 degree) or is a very short tracr (" << factorLastLayers << "%)\n";

   /*
   cout << "Tracks with no EID in first layer: ";
   for (Int_t i=0; i<ntracks; i++) {
      if (!tracks->At(i)) continue;
      if (!tracks->At(i)->At(0) || !tracks->At(i)->At(1)) continue;
      if (tracks->At(i)->getEventID(0) < 0) cout << tracks->At(i)->getEventID(1) << ", ";
   }
   cout << endl;
   */


   for (Int_t i=0; i<ntracks; i++) {
      Track *thisTrack = tracks->At(i);
      if (!thisTrack) continue;
//      if (thisTrack->getTrackLengthmm() < 2) continue;

      Int_t n = thisTrack->GetEntriesFast();

      TPolyLine3D *l = new TPolyLine3D(n);
      l->SetLineWidth(2);
      TPolyMarker3D *trackPoints = new TPolyMarker3D(nClusters, 7);
      
      if (!thisTrack->isFirstAndLastEventIDEqual()) {
         if (!thisTrack->Last()) continue;
         if (!thisTrack->At(thisTrack->GetEntriesFast() - 2)) continue;

         Int_t lastEID = thisTrack->Last()->getEventID();   
         if (lastEID < 0) continue;

         Int_t trackID = tracks->getTrackIdxFromFirstLayerEID(lastEID);

         if (trackID < 0) continue;
         if (!tracks->At(trackID)->At(0)) continue;
   
         Float_t delta = diffmmXY(tracks->At(trackID)->At(0), thisTrack->At(0));
         Float_t phi0 = thisTrack->getSlopeAngleBetweenLayers(1);
         Float_t phi1 = tracks->At(trackID)->getSlopeAngleBetweenLayers(1);
         Float_t deltaphi = fabs(phi0 - phi1);
/*
         if (delta < 0.5 && deltaphi < 1) {
         }

         else if (thisTrack->getWEPL() < 0.6 * getWEPLFromEnergy(run_energy)) {
            // Bad track ends early. OK...
         }
         
         nMissingEID = tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0));
         if (!thisTrack->isFirstAndLastEventIDEqual() || nMissingEID > 0) {
            l->SetLineColor(kRed);
         }
         */
      
         nMissingEID = tracks->getNMissingClustersWithEventID(thisTrack->getEventID(0), thisTrack->Last()->getLayer());
         if (thisTrack->isFirstAndLastEventIDEqual() && nMissingEID == 0) {}
         else {
            l->SetLineColor(kRed);
         }
      }
      /* Old method, draw all tracks with different event IDs red
      if (!thisTrack->isFirstAndLastEventIDEqual()) {
            l->SetLineColor(kRed);
      }
      */

      firstEID = thisTrack->getEventID(0);
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
//      l->SetLineColor(kBlack);
      l->SetLineWidth(3);
      l->Draw();
      trackPoints->Draw();
//    EIDMarker->Draw();
      conflictMarker->Draw();
   }
   view->ShowAxis();
   c1->Update();

   TAxis3D *axis = TAxis3D::GetPadAxis();
   axis->SetTitle("3D view of tracks and clusters");
   axis->SetLabelColor(kBlack);
   axis->SetAxisColor(kBlack);
   axis->SetXTitle("Pixels in X");
   axis->SetYTitle("Layer number");
   axis->SetZTitle("Pixels in Y");
   axis->SetTitleOffset(2);

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
   /*
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

   c1->SaveAs(Form("OutputFiles/figures/testOutput_switchLayer%d.png", switchLayer));

   delete tracks;
   delete conflictTracks;
   delete allConflictPairs;
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

void drawRegionOccupancy(Int_t Runs, Int_t Layer, Float_t energy) {
   run_energy = energy;
   DataInterface *di = new DataInterface();
   CalorimeterFrame *cf = new CalorimeterFrame();
  
   TCanvas *c2 = new TCanvas("c2", "Region dependent occupancy", 1200, 800);
   Int_t counterOverall = 0;
   Int_t fromX =0;
   Int_t toX = 0;
   Int_t counter;
   Int_t nPixels = ny*2;
   TH1F *hCount = new TH1F("hCount", "Priority Encoding (FOCAL) region occupancy @ 1000 protons/event;Occupancy (%);Counts",30,0,15);
   Float_t occupancy;
   
   for (Int_t k=0; k<Runs; k++) {
      di->getDataFrame(k, cf, energy);
      TH2F *Frame2D = cf->getTH2F(Layer);
      if (Frame2D->Integral() == 0) break;
      cout << "Hits in calorimeterframe = " << Frame2D->Integral();

      for (Int_t i=0; i<nx/2; i++) { 
         counter = 0;
         fromX = i*2;
         toX = (i+1)*2;

         for (Int_t j=0; j<ny; j++) {
            if (Frame2D->GetBinContent(Frame2D->FindBin(i, j)) > 0){
               counter++;
               counterOverall++;
            }
         }

         occupancy = (float) counter / nPixels * 100;
         hCount->Fill(occupancy);
      }
      cf->Reset();
   }
 
   delete di;

   c2->cd();
   hCount->SetFillColor(kBlue-7);
   hCount->Draw();
      
   printf("The overall occupancy is %.2f %%", 100*(float) counterOverall / nx / ny);

   gPad->Update();
   TLine *l = new TLine(12.5, 0, 12.5, 500);
   l->Draw();
}

void drawFrame2D(Int_t Runs, Int_t dataType, Int_t Layer, Float_t energy) {
   run_energy = energy;
   DataInterface *di = new DataInterface();
   CalorimeterFrame *cf = new CalorimeterFrame();
  
   TCanvas *c1 = new TCanvas("c1", "Hit distribution in layer", 1200, 800);

   TList *histogramList = new TList;
   for (Int_t k=0; k<Runs; k++) {
      if (dataType == kData) {
         di->getDataFrame(k, cf, energy);
      }
      else {
         di->getMCFrame(k, cf);
      }

      TH2F *Frame2D = cf->getTH2F(Layer);
      histogramList->Add(Frame2D);
   }

   delete di;
   TH2F *Frame2D = new TH2F("Frame2D", Form("Hit distribution in layer %i", Layer),
                              nx, 0, nx, ny, 0, ny);

   c1->cd();
   Frame2D->Merge(histogramList);
   Frame2D->Draw("COLZ");
   gStyle->SetOptStat(0);
}

void drawData3D(Int_t Runs, Float_t energy) {
   run_energy = energy;

   DataInterface *di = new DataInterface();

   TH3F *Frame3D = new TH3F("Frame3D", "3D map of energy deposition [keV]", 
                              100, -120, -50, 100, 0, nx, 100, 0, ny);

   Frame3D->SetXTitle("Z axis");
   Frame3D->SetYTitle("X axis"); // to get projection right (Z is depth, not up)
   Frame3D->SetZTitle("Y axis"); 

   for (Int_t run=0; run<Runs; run++) {
      di->getMCData(run, Frame3D);
   }
   
   delete di;

   Frame3D->Draw("LEGO");
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


void compareClusterSizes(Int_t Runs, Bool_t recreate, Float_t energy) {
   run_energy = energy;
   Tracks    * MCTracks = nullptr;
   Tracks    * DataTracks = nullptr;
   Cluster   * thisCluster = nullptr;
   Clusters  * MCClusters = nullptr;
   Clusters  * DataClusters = nullptr;
   Int_t       nTracksMC, nTracksData, thisLayer, thisSize;
   Bool_t      useChip = true;
   Int_t       nClustersInChip[30] = {0};
   Bool_t      useLowerMCEnergy = true;
   Float_t     altEnergy = energy;

   Int_t fChip = 1;
   if (useChip) fChip = 4;

   const Int_t nLayersToUse = 7*4;

   if (useLowerMCEnergy) {
      if (energy == 188) altEnergy = 184;
      if (energy == 180) altEnergy = 171;
      if (energy == 170) altEnergy = 166;
      if (energy == 160) altEnergy = 155;
      if (energy == 150) altEnergy = 150;
   }
   
   cout << "Finding MC tracks...\n";
   MCTracks = loadOrCreateTracks(recreate, Runs, kMC, altEnergy);

   cout << "Finding EXP tracks...\n";
   DataTracks = loadOrCreateTracks(recreate, Runs, kData, energy);

   MCClusters = getClusters(Runs, kMC, kCalorimeter, energy);
   DataClusters = getClusters(Runs, kData, kCalorimeter, energy);

   TCanvas *c2 = new TCanvas("c2", "Individual cluster size distributions MC", 1200, 800);
   c2->Divide(4, 2, 0.01, 0.01);
   TCanvas *c3 = new TCanvas("c3", "Individual cluster size distributions DATA", 1200, 800);
   c3->Divide(4, 2, 0.01, 0.01);
   TCanvas *c1 = new TCanvas("c1", "Cluster size distribution comparison", 1022, 645);
   c1->Divide(1,2,0.01, 0.01);

   vector<TH1F*> *hCSVectorMC = new vector<TH1F*>;
   for (Int_t i=0; i<nLayersToUse; i++) {
      hCSVectorMC->push_back(new TH1F(Form("hCSIndMC_%d", i), Form("CS histogram %d", i), 60, 0, 60));
   }
   
   vector<TH1F*> *hCSVectorData = new vector<TH1F*>;
   for (Int_t i=0; i<nLayersToUse; i++) {
      hCSVectorData->push_back(new TH1F(Form("hCSIndData_%d", i), Form("CS histogram %d", i), 60, 0, 60));
   }
   
   vector<TH1F*> *hCSVectorDataCorrected = new vector<TH1F*>;
   for (Int_t i=0; i<nLayersToUse; i++) {
      hCSVectorDataCorrected->push_back(new TH1F(Form("hCSIndDataCorrected_%d", i), Form("CS histogram %d", i), 60, 0, 60));
   }

   for (Int_t i=0; i<MCClusters->GetEntriesFast(); i++) {
      thisCluster = MCClusters->At(i);
        if (useChip)    thisLayer = thisCluster->getChip();
        else            thisLayer = thisCluster->getLayer();
      
        if (thisLayer >= nLayersToUse) continue;

      thisSize = thisCluster->getDepositedEnergy(false);
      hCSVectorMC->at(thisLayer)->Fill(thisSize);
   }
   
   for (Int_t i=0; i<DataClusters->GetEntriesFast(); i++) {
      thisCluster = DataClusters->At(i);

      if (useChip) {
        thisLayer = thisCluster->getChip();
        nClustersInChip[thisLayer]++;
      }
      else {
         thisLayer = thisCluster->getLayer();
      }

      if (thisLayer >= nLayersToUse) continue;
      thisSize = thisCluster->getDepositedEnergy(false);
      hCSVectorData->at(thisLayer)->Fill(thisSize);
      hCSVectorDataCorrected->at(thisLayer)->Fill(thisCluster->getDepositedEnergy(true));
   }

   Float_t layerMC[nLayersToUse]; //= {0, 1, 2, 3, 4, 5, 6, 7};
   Float_t layerData[nLayersToUse];// = {0, 1, 2, 3, 4, 5, 6, 7};
     
   for (Int_t i=0; i<nLayersToUse; i++) {
        layerMC[i] = i;
        layerData[i] = i;
   }

   Float_t errorLayer[nLayersToUse] = {0};
   Float_t clusterSizeMC[nLayersToUse] = {0};
   Float_t errorClusterSizeMC[nLayersToUse] = {0};
   Float_t clusterSizeData[nLayersToUse] = {0};
   Float_t clusterSizeDataCorrected[nLayersToUse] = {0};
   Float_t errorClusterSizeData[nLayersToUse] = {0};
   Float_t errorClusterSizeDataCorrected[nLayersToUse] = {0};
   Float_t clusterSizeRatio[nLayersToUse] = {0};
   Float_t errorClusterSizeRatio[nLayersToUse] = {0};


   for (Int_t i=0; i<nLayersToUse; i++) {
      layerMC[i] -= 0.07;
      layerData[i] += 0.07;

      clusterSizeMC[i] = hCSVectorMC->at(i)->GetMean();
      clusterSizeData[i] = hCSVectorData->at(i)->GetMean();
      clusterSizeDataCorrected[i] = hCSVectorDataCorrected->at(i)->GetMean();
      errorClusterSizeMC[i] = hCSVectorMC->at(i)->GetRMS();
      errorClusterSizeData[i] = hCSVectorData->at(i)->GetRMS();
      errorClusterSizeDataCorrected[i] = hCSVectorDataCorrected->at(i)->GetRMS();
      clusterSizeRatio[i] = clusterSizeData[i] / clusterSizeMC[i];
      Float_t e_high = (clusterSizeData[i] + errorClusterSizeData[i]/2) / (clusterSizeMC[i] - errorClusterSizeMC[i]/2);
      Float_t e_low  = (clusterSizeData[i] - errorClusterSizeData[i]/2) / (clusterSizeMC[i] + errorClusterSizeMC[i]/2);
      errorClusterSizeRatio[i] = (e_high - e_low);
   }

   for (Int_t i=0; i<nLayersToUse; i++) {
      cout << Form("Chip %d has a correction factor MC/data of %.3f with %d entries.\n", i, clusterSizeMC[i]/clusterSizeData[i], nClustersInChip[i]);
   }

   TGraphErrors * graphCSMC = new TGraphErrors(nLayersToUse, layerMC, clusterSizeMC, errorLayer, errorClusterSizeMC);
   TGraphErrors * graphCSMC2 = new TGraphErrors(nLayersToUse, layerMC, clusterSizeMC, errorLayer, errorClusterSizeMC);
   TGraphErrors * graphCSData = new TGraphErrors(nLayersToUse, layerData, clusterSizeData, errorLayer, errorClusterSizeData);
   TGraphErrors * graphCSDataCorrected = new TGraphErrors(nLayersToUse, layerData, clusterSizeDataCorrected, errorLayer, errorClusterSizeDataCorrected);
   TGraph       * graphRatios = new TGraphErrors(nLayersToUse, layerMC, clusterSizeRatio);
   
   graphCSMC->SetTitle(Form("Uncalibrated energy deposition distribution comparison at %.0f MeV;Chip number; Energy deposition [keV/#mum]", run_energy));
   graphCSMC2->SetTitle(Form("Calibrated energy deposition distribution comparison at %.0f MeV;Chip number; Energy deposition [keV/#mum]", run_energy));

   graphCSMC->SetMinimum(-2);
   graphCSMC->SetMaximum(10);
   graphCSMC->SetMarkerStyle(21);
   graphCSMC->SetMarkerColor(kBlue);
   graphCSMC2->SetMinimum(-2);
   graphCSMC2->SetMaximum(10);
   graphCSMC2->SetMarkerStyle(21);
   graphCSMC2->SetMarkerColor(kBlue);
   graphCSData->SetMarkerStyle(22);
   graphCSData->SetMarkerColor(kRed);
   graphCSData->SetMarkerSize(1.5);
   graphCSData->GetXaxis()->SetTitleFont(22);
   graphCSData->GetXaxis()->SetLabelFont(22);
   graphCSData->GetYaxis()->SetTitleFont(22);
   graphCSData->GetYaxis()->SetTitleSize(0.07);
   graphCSData->GetYaxis()->SetTitleOffset(0.3);
   graphCSData->GetYaxis()->SetLabelFont(22);
   graphCSData->GetXaxis()->SetNdivisions(54);
   graphCSDataCorrected->SetMarkerStyle(22);
   graphCSDataCorrected->SetMarkerColor(kRed);
   graphCSDataCorrected->SetMarkerSize(1.5);
   graphCSMC->SetMarkerSize(1.25);
   graphCSMC->GetXaxis()->SetTitleFont(22);
   graphCSMC->GetXaxis()->SetLabelFont(22);
   graphCSMC->GetYaxis()->SetTitleFont(22);
   graphCSMC->GetXaxis()->SetTitleSize(0.07);
   graphCSMC->GetXaxis()->SetLabelSize(0.07);
   graphCSMC->GetYaxis()->SetTitleSize(0.07);
   graphCSMC->GetYaxis()->SetTitleOffset(0.3);
   graphCSMC->GetYaxis()->SetLabelFont(22);
   graphCSMC->GetXaxis()->SetNdivisions(54);
   graphCSMC2->SetMarkerSize(1.25);
   graphCSMC2->GetXaxis()->SetTitleFont(22);
   graphCSMC2->GetXaxis()->SetTitleSize(0.07);
   graphCSMC2->GetXaxis()->SetLabelSize(0.07);
   graphCSMC2->GetXaxis()->SetLabelFont(22);
   graphCSMC2->GetYaxis()->SetTitleFont(22);
   graphCSMC2->GetYaxis()->SetTitleSize(0.07);
   graphCSMC2->GetYaxis()->SetTitleOffset(0.3);
   graphCSMC2->GetYaxis()->SetLabelFont(22);
   graphCSMC2->GetXaxis()->SetNdivisions(54);
   graphRatios->SetMinimum(0.5);
   graphRatios->SetMaximum(1.5);
   graphRatios->SetTitle("Cluster size distribution ratios (DATA / MC) at 188 MeV;Layer Number;Cluster Size Ratio (DATA/MC)");
   if (useChip) {
      graphRatios->GetXaxis()->SetTitle("Chip number");
      graphRatios->GetYaxis()->SetTitle("E_{dep} ratio (data / MC)");
   }
   graphRatios->SetMarkerStyle(21);
   graphRatios->SetMarkerColor(kBlue);
   graphRatios->SetMarkerSize(1.5);
   graphRatios->GetXaxis()->SetTitleFont(22);
   graphRatios->GetXaxis()->SetLabelFont(22);
   graphRatios->GetYaxis()->SetTitleFont(22);
   graphRatios->GetYaxis()->SetLabelFont(22);
   graphRatios->GetXaxis()->SetNdivisions(54);

   TLegend *leg = new TLegend(0.16, 0.70, 0.28, 0.86);
   leg->SetTextFont(22);
   leg->AddEntry(graphCSMC, "MC", "ep");
   leg->AddEntry(graphCSData, "Exp. data", "ep");

   c1->cd(1);
   graphCSMC->Draw("AP");
   graphCSData->Draw("P");
   leg->Draw();

   gPad->Update();
   TPaveText *title = (TPaveText*) gPad->GetPrimitive("title");
   title->SetTextFont(22);
   gPad->Modified();

   c1->cd(2);
   graphCSMC2->Draw("AP");
   graphCSDataCorrected->Draw("P");
   
   gPad->Update();
   TPaveText *title2 = (TPaveText*) gPad->GetPrimitive("title");
   title2->SetTextFont(22);
   gPad->Modified();

   for (Int_t i=0; i<nLayersToUse; i++) {
      c2->cd(i+1);
      hCSVectorMC->at(i)->Draw();
      c3->cd(i+1);
      hCSVectorData->at(i)->Draw();
   }
}

Bool_t getCutTrackLength(Float_t energy, Track *track) {
   Int_t minTL = getMinimumTrackLength(energy);
   Float_t TL = track->getTrackLengthmm();

   Bool_t cutTL = (TL > minTL) ? true : false;

   return cutTL;
}

Bool_t getCutWEPL(Track *track) {
   Float_t minTLWEPL = 150;
   Float_t WEPL = track->getWEPL();

   Bool_t cutTL = (WEPL > minTLWEPL) ? true : false;

   return cutTL;
}

Bool_t getCutChipNumber(Track *track) {
   Int_t x0 = track->getX(0);
   Int_t y0 = track->getY(0);
   Int_t chip = (x0 >= nx/2) + 2 * (y0 < ny/2);
   Bool_t cutChipNumber = (chip<2) ? true : false;

   return cutChipNumber;
}

Bool_t getCutBraggPeakInTrack(Track *track) {
   Float_t braggPeakRatio = 2.5;

   Int_t lastBin = track->GetEntriesFast() - 1;
   if (lastBin < 4) return false;

   Float_t lowStd = track->getStdSizeToIdx(lastBin-1);
   Int_t   nEmptyBins = track->getNMissingLayers();


   if (lowStd > 4) return false;
   if (nEmptyBins > 1) return false;

   Float_t rMean = track->getMeanSizeToIdx(lastBin - 1);
   Float_t rrMean = track->getMeanSizeToIdx(lastBin - 2);

   Float_t r = track->getSize(lastBin) / rMean;
   Float_t rr = track->getSize(lastBin-1) / rrMean;

   if (r > braggPeakRatio) return true;
   else {
      if (rr > braggPeakRatio) return true;
      else return false;
   }
}

void drawIndividualGraphs(TCanvas *cGraph, TGraphErrors* outputGraph, Float_t fitRange, Float_t fitScale, Float_t fitError, Int_t fitIdx, Int_t eventID, Float_t *x_energy, Float_t *y_energy) {
   cGraph->cd(fitIdx+1);
   Bool_t kDrawFit = true;
   Bool_t kDrawText = true;

   outputGraph->SetMinimum(0);
   outputGraph->SetMaximum(20);
   outputGraph->SetTitle("");
   
   if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
      outputGraph->GetXaxis()->SetTitle("Water Equivalent Thickness [mm]");
   }
   
   else if (kOutputUnit == kPhysical) {
      outputGraph->GetXaxis()->SetTitle("Projected range [mm]");
   }

   outputGraph->GetYaxis()->SetTitle("Energy deposition [keV/#mum]");
   outputGraph->GetYaxis()->SetTitleOffset(1);
   outputGraph->GetXaxis()->SetTitleSize(0.05);
   outputGraph->GetYaxis()->SetTitleSize(0.05);
   outputGraph->GetXaxis()->SetLabelSize(0.04);
   outputGraph->GetYaxis()->SetLabelSize(0.04);
   outputGraph->GetXaxis()->SetTitleFont(22);
   outputGraph->GetXaxis()->SetLabelFont(22);
   outputGraph->GetYaxis()->SetTitleFont(22);
   outputGraph->GetYaxis()->SetLabelFont(22);

   Float_t low = getUnitFromEnergy(0);
   Float_t high = getUnitFromEnergy(run_energy) * 1.7;

   outputGraph->GetXaxis()->SetLimits(low, high);

   outputGraph->SetMarkerColor(4);
   outputGraph->SetMarkerStyle(21);
   outputGraph->Draw("PA");
   
   TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
   func->SetParameters(fitRange, fitScale);
   func->SetNpx(500);

   func->Draw("same");

   if (x_energy) {
      Float_t WEPLFactor = getWEPLFactorFromEnergy(run_energy);
      Long64_t n=0;
      Long64_t j=0;
      
      for (Long64_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) {
         if (y_energy[i] == 0) {
            n = j;
            break;
         }
         y_energy[i] *=  2;
         
         j++;
      }
      
      if (n>1) {
         TGraph *energyLoss = new TGraph(n, (x_energy + eventID*sizeOfEventID), (y_energy + eventID*sizeOfEventID));
         energyLoss->SetLineColor(kBlack); energyLoss->SetLineWidth(2);
         energyLoss->Draw("same");
      }
      
      Float_t lastX = x_energy[eventID*sizeOfEventID + n-1];
      TLine *l = new TLine(lastX, 0, lastX, 600);
      l->Draw();
      
      Float_t realEnergy = getEnergyFromWEPL(x_energy[eventID*sizeOfEventID + n-1]);
      TLatex *text2 = new TLatex(15, 28, Form("'Real' energy: %.1f #pm MeV", realEnergy));
      text2->SetTextSize(0.045);
      text2->SetTextFont(22);
      text2->Draw();
   }

   if (kDrawText) {
      TLatex *text = new TLatex(15, 18, Form("Fitted range: %.1f #pm %.1f mm", fitRange, fitError));
      text->SetTextSize(0.05);
      text->SetTextFont(22);
      text->Draw();
   }
   
   outputGraph->GetXaxis()->SetRangeUser(0, 350);

   cGraph->Update();
}

TF1 *  doSimpleGaussianFit (TH1F *h, Float_t *means, Float_t *sigmas) {
   Float_t estimated_energy = 0, estimated_range = 0;
   Float_t estimated_energy_error = 0;
   Float_t sum_constant = 0, sumSigma = 0;

   Float_t nominalMean = getUnitFromEnergy(run_energy);
   Float_t nominalSigma = getUnitStragglingFromEnergy(run_energy, 0);

   TAxis *axis = h->GetXaxis();
   TF1 *gauss = new TF1("gauss", "gaus");
   gauss->SetParameter(1, nominalMean);
   gauss->SetParameter(2, nominalSigma);
   h->Fit("gauss", "B,W,M");

   Float_t mu, sigma;
   mu = gauss->GetParameter(1);
   sigma = fabs(gauss->GetParameter(2));

   // nominalSigma is too low at very low energies
   // fitted sigma may be too low at high energies, where fit can be ~= 1 layer
   Float_t sigmaToUse = fmax(sigma, nominalSigma);

   Float_t sigmaFrom = fmax(mu - 4 * sigmaToUse, 0);
   Float_t sigmaTo = fmin(mu + 4 * sigmaToUse, nominalMean * 1.4 + 10);

   Int_t binSigmaFrom = axis->FindBin(sigmaFrom);
   Int_t binSigmaTo = axis->FindBin(sigmaTo);
   
   if (binSigmaTo > h->GetNbinsX()) binSigmaTo = h->GetNbinsX();
   if (binSigmaTo < binSigmaFrom) binSigmaTo = h->GetNbinsX();

   printf("doSimpleGaussianFit: Searching from %.2f (bin %d)  to %.2f (bin %d).\n", sigmaFrom, binSigmaFrom, sigmaTo, binSigmaTo);

   Float_t squareMeanDifference = 0;
   Float_t empiricalMean = 0;
   Int_t N = 0;

   for (Int_t i=binSigmaFrom; i<=binSigmaTo; i++) {
      empiricalMean += h->GetBinContent(i) * axis->GetBinCenter(i);
      N += h->GetBinContent(i);
   } 
  
   empiricalMean /= N;

   for (Int_t i=binSigmaFrom; i<=binSigmaTo; i++) {
      squareMeanDifference += h->GetBinContent(i) * pow(axis->GetBinCenter(i) - empiricalMean, 2);
   }

   Float_t empiricalSigma = sqrt(abs(squareMeanDifference / N));
   cout << "The empirical estimated range is " << empiricalMean << " mm. (which is " << getEnergyFromUnit(empiricalMean) << " MeV).\n";
   cout << "The empirical standard deviation is " << empiricalSigma << " mm.\n";
   
   Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;
   
   ofstream file2("OutputFiles/result_makebraggpeakfit.csv", ofstream::out | ofstream::app);
   // absorber thickness; energy; nominal range; estimated range; range sigma
   if (run_degraderThickness == 0) {
      file2 << readoutAbsorber << " " << run_energy << " " << getUnitFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
   }
   else {
      file2 << readoutAbsorber << " " << run_degraderThickness << " " << getUnitFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
   }

   file2.close();
   
   means[9] = empiricalMean;
   sigmas[9] = empiricalSigma;

   return gauss;

}

Float_t doNGaussianFit ( TH1F *h, Float_t *means, Float_t *sigmas) {
   TF1 *gauss;
   cout << "Energy " << run_energy << endl;
   cout << "Layer;\t Constant;\t Mean;\t\t Energy;\t Sigma;\t\t Fits in layer;\t Chi2/n\n";
   
   Float_t constant, mean, lEnergy, sigma;
   Float_t meanValueFromSumming = 0;
   
   Float_t array_mean[3] = {};
   Int_t array_layer[3] = {};
   Float_t array_constant[3] = {};
   Float_t array_sigma[3] = {};
   Float_t array_sum[3] = {};
   Float_t array_energy[3] = {};
   Float_t array_f[3] = {};
   Float_t array_chi2n[3] = {};
   
   TAxis *axis = h->GetXaxis();
   Float_t fullIntegral = h->Integral();
   Int_t binSearchFrom = axis->FindBin(getWEPLFromEnergy(run_energy * 0.9));
   Int_t binSearchTo   = axis->FindBin(getWEPLFromEnergy(run_energy * 1.1));
   Float_t partialIntegral = h->Integral(binSearchFrom, binSearchTo);

   for (Int_t i=binSearchFrom; i<binSearchTo; i++) {
      meanValueFromSumming += axis->GetBinCenter(i) * h->GetBinContent(i) / partialIntegral;
   }

   Float_t error = abs(100 * (meanValueFromSumming - getWEPLFromEnergy(run_energy)) / getWEPLFromEnergy(run_energy));
   cout << "\033[1mUsing the sum method, the weighted mean value of the distribution is " << meanValueFromSumming << " mm. (" << error << ")\033[0m\n";
   
   Float_t maxBinHeight = h->GetMaximum();

   Bool_t isLastLayer, wasLastLayer = false;;
   

   Int_t j=0;
   for (Int_t i=0; i<nLayers; i++) {
      if (getWEPLFromTL(getLayerPositionmm(i)) > getUnitFromEnergy(run_energy*1.2)) continue;
      if (getWEPLFromTL(getLayerPositionmm(i)) < getWEPLFromEnergy(run_energy * 0.8)) continue;

      Float_t searchFrom = getWEPLFromTL(getLayerPositionmm(i))+10;
      Float_t searchTo = getWEPLFromTL(getLayerPositionmm(i+1))+10;
      
      Int_t bmin = axis->FindBin(searchFrom);
      Int_t bmax = axis->FindBin(searchTo);
      
      Float_t integral = h->Integral(bmin,bmax);
      Float_t ratio = integral / fullIntegral;
      
      isLastLayer = ((getWEPLFromTL(getLayerPositionmm(i)) > getUnitFromEnergy(run_energy-10)) && !wasLastLayer && ratio > 0.01) ;
   
      if (i<=3) continue;
      if (ratio < 0.1 && !isLastLayer) continue;
      
      gauss = new TF1(Form("Gaus_%d", i), "gaus(0)", searchFrom, searchTo);
   
      gauss->SetLineWidth(2);

      sigma = 3;
      if (isLastLayer && ratio < 0.05) sigma = 0.2;
      
      gauss->SetParameters(10, (searchFrom+searchTo)/2, sigma);
      gauss->SetParLimits(0, 0, maxBinHeight);
      gauss->SetParLimits(1, searchFrom+5, searchTo);
      gauss->SetParLimits(2, 2, 9);
      
      h->Fit(gauss, "M, B, WW, Q, 0", "", searchFrom, searchTo);
      
      sigma = gauss->GetParameter(2);
      constant = gauss->GetParameter(0);
      mean = gauss->GetParameter(1);
      lEnergy = getEnergyFromUnit(mean);
      
      Float_t integralSigma = h->Integral(axis->FindBin(mean - sigma), axis->FindBin(mean + sigma));
      Float_t chi2 = gauss->GetChisquare();
      Float_t chi2n = chi2 / integral;
      Float_t chi2nSigma = chi2 / integralSigma;
      
      cout << Form("Searching from %.1f to %.1f, with midpoint at %.1f. Found best fit @ %.1f with chi2 = %.2f and chi2/n = %.2f and chi2/nSIGMA = %.2f, ratio = %.2f.\n", searchFrom, searchTo,(searchTo+searchFrom)/2 , mean, chi2, chi2n, chi2nSigma, ratio);

   /*
      if (run_energy > 182 && run_energy < 193) {
         if (chi2n > 200 || sigma < 2.5 || constant < 3) {
            delete gauss;
            continue;
         }
      }
   */

      if (ratio > 0.1) {
         gauss->SetLineColor(kRed);
         gauss->SetLineWidth(3);
//         gauss->Draw("same");
         cout << Form("%d;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f;\t %.2f\n", i, constant, mean, lEnergy, sigma, ratio, chi2n);
         
         if (j<3) {
            array_mean[j] = mean;
            array_layer[j] = i;
            array_constant[j] = constant;
            array_energy[j] = lEnergy;
            array_sigma[j] = sigma;
            array_f[j] = ratio;
            array_sum[j] = integralSigma;
         }
         
         sigmas[j] = sigma;
         means[j++] = mean;
      }
      wasLastLayer = isLastLayer;
   }

   Float_t estimated_energy = 0, estimated_range = 0;
   Float_t estimated_energy_error = 0;
   Float_t sum_constant = 0, sumSigma = 0;
   for (Int_t i=0 ; i<3; i++) {
      estimated_range += array_sum[i] * array_mean[i];
      sum_constant += array_sum[i];
      sumSigma += pow(array_sigma[i], 2);
   }

   sumSigma = sqrt(sumSigma);

   estimated_range /= sum_constant;

   estimated_energy = getEnergyFromUnit(estimated_range);
   estimated_energy_error = getEnergyFromWEPL(estimated_range + sumSigma/2) - getEnergyFromWEPL(estimated_range - sumSigma/2);
   cout << "ESTIMATED ENERGY FROM RUN IS " << estimated_energy << " +- " << estimated_energy_error << endl;
   cout << "Estimated range = " << estimated_range << " +- " << sumSigma << endl;
  
   Float_t lastMean = array_mean[1];
   Float_t lastSigma = array_sigma[1];
   if (lastMean == 0 || axis->FindBin(lastMean) > h->GetNbinsX()) {
      lastMean = array_mean[0];
      lastSigma = array_sigma[0];
   }

   printf("lastMean = %.2f. lastSigma = %.2f\n", lastMean, lastSigma);
   printf("sigmaFrom = %.2f\n", lastMean - 9 * lastSigma);

   Int_t binSigmaFrom = axis->FindBin(lastMean - 9*lastSigma);

   if (lastMean == 0) {
      binSigmaFrom = axis->FindBin(getWEPLFromEnergy(run_energy)*0.9);
      cout << "Setting binSigmaFrom to " << binSigmaFrom;
   }

   Float_t squareMeanDifference = 0;
   Float_t empiricalMean = 0;
   Int_t N = 0;

   printf("Looping from %d to %d.\n", binSigmaFrom, h->GetNbinsX());

   for (Int_t i=binSigmaFrom; i<=h->GetNbinsX(); i++) {
      empiricalMean += h->GetBinContent(i) * axis->GetBinCenter(i);
      N += h->GetBinContent(i);
   } 
  
   empiricalMean /= N;

   printf("empiricalMean = %.2f\n", empiricalMean);

   for (Int_t i=binSigmaFrom; i<=h->GetNbinsX(); i++) {
      squareMeanDifference += h->GetBinContent(i) * pow(axis->GetBinCenter(i) - empiricalMean, 2);
   }
   printf("squareMeanDifference = %.2f\n", squareMeanDifference);

   Float_t empiricalSigma = sqrt(abs(squareMeanDifference / N));
   cout << "The empirical estimated range is " << empiricalMean << " mm. (which is " << getEnergyFromWEPL(empiricalMean) << " MeV).\n";
   cout << "The empirical standard deviation is " << empiricalSigma << " mm.\n";
   
   ofstream file("OutputFiles/output_gauss.csv", ofstream::out | ofstream::app);
   // run_energy; layer[i], constant[i], mpv[i], energy[i], sigma[i], ratio[i];  i 1->3
   
   file << run_energy << "; ";
   for (Int_t j=0; j<3; j++) {
      file << Form("%d; %.5f; %.5f; %.5f; %.5f; %.5f; %.5f",
                array_layer[j], array_constant[j], array_mean[j],
                array_energy[j], array_sigma[j], array_f[j], array_chi2n[j]);
   }

   file << endl;
   file.close();
   
   Float_t last_range = array_mean[0];
   if (array_constant[1] > 0.1) last_range = array_mean[1];
   if (array_constant[2] > 0.1) last_range = array_mean[2];

   Float_t nominalSigma = getWEPLStragglingFromEnergy(run_energy, 0);

   Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;

   ofstream file2("OutputFiles/result_makebraggpeakfit.csv", ofstream::out | ofstream::app);
   // absorber thickness; energy; nominal range; estimated range; range sigma
   if (run_degraderThickness == 0) {
      file2 << readoutAbsorber << " " << run_energy << " " << getWEPLFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
   }
   else {
      file2 << readoutAbsorber << " " << run_degraderThickness << " " << getWEPLFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
   }

   file2.close();

   means[9] = empiricalMean;
   sigmas[9] = empiricalSigma;

   return empiricalMean;
}

void generateWaterDegraderValues() {
   // Use the Water.csv file to generate residual energies after a varying water phantom thickness
   // Store the values in Data/Ranges/EnergyAfterDegrader.csv

   // getEnergyAtWEPL(E0, depth)

   ofstream file("Data/Ranges/EnergyAfterDegraderG4.csv", ofstream::out);
   for (Int_t depth = 0; depth < 380; depth++) {
      file << depth << " " << getEnergyAtWEPL(250, (float) depth) << endl;
      printf("Writing to file: %d %.4f\n", depth, getEnergyAtWEPL(250, (float) depth));
   }
}
