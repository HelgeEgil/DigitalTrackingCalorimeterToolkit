#define Analysis_cxx
#include "Analysis.h"
#include "Constants.h"
#include "MaterialConstants.h"
#include "Track_conversion.h"
#include "TFocal.h"
#include "Tracks.h"
#include "Tools.h"
#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
#include <TEllipse.h>
#include <vector>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TPaveStats.h>
#include <algorithm>
#include <ctime>
#include <TView.h>
#include <TLeaf.h>
#include <TArrow.h>
#include <TF1.h>
#include "Math/ProbFunc.h"
#include <ctime>

using namespace std;

/*
void GetTrackerStatistics(Int_t Events, Int_t Runs) {
	
	Focal f;

	TrackerCollection *allTrackerFrames = new TrackerCollection(2*Events*Runs);
	Clusters *restPoints = new Clusters(Runs*200);


	allTrackerFrames->SetOwner(true);

	for (Int_t i=0; i<Runs; i++) {
		// events should be low in this case !!
		// Not more than 50 maybe?

		Int_t EventFrom = Runs * i;
		Int_t EventTo = EventFrom + Events;

		// could use absorbobjects(tracks), but i don't get it to work
		// Track objects are small though.
		
		TrackerCollection *trackers = f.FindTrackerFrame(EventFrom, EventTo, restPoints);
		trackers->SetOwner(true);

		cout << "--- Run " << i+1 << " of " << Runs << " complete --- \n";

		for (Int_t j=0; j<trackers->GetEntriesFast(); j++) {
			allTrackerFrames->AddTracker(trackers->At(j));
		}
		delete trackers;
	}

	cout << "From " << Runs << " runs with " <<Events << " events in each, we have collected"
			<< allTrackerFrames->GetEntriesFast() << " tracks.\n";

	cout << "Here is info about the first 10 tracker frames: \n\n\n";
	for (Int_t i=0; i<10; i++) {
		if (!allTrackerFrames->At(i)) continue;

		cout << "Tracker frame " << i << " : \n";
		cout << allTrackerFrames->At(i) << endl;
	}

	// start comment here ...

	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
	TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);

	TH1F *hTrackLengths = new TH1F("hTrackLengths", "Track Lengths (MC)", 100, 0, 120);
	TH2F *hClusterSizeAlongTrack = new TH2F("hClusterSizeAlongTrack",
				"Cluster size along track length (MC)", 1.5*nlayers, 0, 1.5*nlayers*dz, 50, 0, 50);
	TH1F *hStraightness = new TH1F("hStraightness", "Sinuosity plot (MC)", 500, 1, 1.01);
	TH1F *hSlope = new TH1F("hSlope", "Proton angle plot (MC)", 500, 0, 20);
	
	hTrackLengths->SetXTitle("Track length [mm]");
	hClusterSizeAlongTrack->SetXTitle("Track length [mm]");
	hClusterSizeAlongTrack->SetYTitle("Cluster size [# of pixels]");
	hStraightness->SetXTitle("Sinuosity parameter");
	hSlope->SetXTitle("Total track angle (degree)");

	Float_t trackLengthSoFar = 0;
	for (Int_t i=0; i<allTracks->GetEntriesFast(); i++) {
		hTrackLengths->Fill(allTracks->GetTrackLengthmm(i));
		hStraightness->Fill(allTracks->GetSinuosity(i));
		hSlope->Fill(allTracks->GetSlopeAngle(i));
		for (Int_t j=0; j<allTracks->GetEntriesFast(i); j++) {
			trackLengthSoFar += allTracks->At(i)->GetTrackLengthAtmm(j);
			hClusterSizeAlongTrack->Fill(trackLengthSoFar, allTracks->At(i)->GetSize(j));
			cout << "GetSize for cluster at " << *allTracks->At(i)->At(j) << ": " << allTracks->At(i)->GetSize(j) << endl;
		} // end loop over clusters
		trackLengthSoFar = 0;
	} // end loop over tracks

	c1->cd();
		hTrackLengths->Draw();
	c2->cd();
		gStyle->SetOptStat(0);
		hClusterSizeAlongTrack->Draw("COLZ");
	c3->cd();
		hStraightness->Draw();
	c4->cd();
		hSlope->Draw();

	// end comment here ...

	delete allTrackerFrames;
}

*/

void drawTrackAngleAtVaryingRunNumbers(Int_t dataType, Float_t energy) {
	Int_t nRuns = 0;

	for (Int_t i=26; i<34; i++) {
		nRuns = pow(2, 4 + 0.25 * i) + 0.5;

		cout << "Reconstructing tracks with " << nRuns << " protons per frame.";

		kEventsPerRun = nRuns;

		Tracks * tracks = loadOrCreateTracks(1, 100, dataType, energy);
		tracks->extrapolateToLayer0();

		char * sDataType = getDataTypeChar(dataType);
		TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
		TH1F *hAngles = new TH1F("hAngles", Form("Proton angle plot with %d protons in frame (%s)", nRuns, sDataType), 500, 0, 30);
		hAngles->SetXTitle("Protons angle from initial measurement to layer 1");
		hAngles->SetXTitle("Number of protons");
		hAngles->SetFillColor(kCyan-8);
		hAngles->SetLineColor(kBlack);
		gStyle->SetOptStat(0);

		Track *thisTrack;
		for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
			thisTrack = tracks->At(i);
			if (!thisTrack) continue;
			hAngles->Fill(thisTrack->getSlopeAngleAtLayer(1));
		}

		hAngles->Draw();
		c1->SaveAs(Form("figures/angles/angles_layer1_with_nRuns-%d.png", nRuns));
//		break;

		delete tracks;
		delete hAngles;
		delete c1;
	}
}

void getTrackStatistics(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	Focal f;

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
	
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
		hAngles->at(layer)->SetXTitle("Track Length [mm]");
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

	Track *thisTrack;
	for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
		thisTrack = tracks->At(i);

		Float_t TL = thisTrack->getTrackLengthmm();
		Int_t x0 = thisTrack->getX(0);
		Int_t y0 = thisTrack->getY(0);

		cutTL = (TL > kMinimumTracklength) ? true : false;

		chip = (x0 >= nx) + 2 * (y0 < ny);
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
				hAngles->at(layer)->Fill(thisTrack->getSlopeAngleAtLayer(j));
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
	}

	for (Int_t track=0; track<nTracksToPlot; track++) {
		c7->cd(track+1);
		gPad->DrawFrame(0, 0, 50, 35);
		hFollowTrack->at(track)->SetFillColor(kBlue-2);
		hFollowTrack->at(track)->Draw("same");
	}

	for (Int_t layer=0; layer<9; layer++) {
		c8->cd(layer+1);
		hAngles->at(layer)->SetFillColor(kRed-9+layer);
		hAngles->at(layer)->Draw("same");
	}

	delete tracks;
}
//
void drawClusterShapes(Int_t Runs, Bool_t dataType, Bool_t recreate, Float_t energy) {
	// get vector of TH2F's, each with a hits distribution and cluster size
	// dataType = kMC (0) or kData (1)

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

	Focal *f = new Focal();
	CalorimeterFrame *cf = new CalorimeterFrame();
	Hits *hits = new Hits(nHits);
	vector<Hits*> * tempClusterHitMap;
	vector<Hits*> * clusterHitMap = new vector<Hits*>;
	clusterHitMap->reserve(Runs*500*kEventsPerRun);

	for (Int_t i=0; i<Runs; i++) {
		if (dataType == kMC) {
			f->getMCFrame(i, cf);
			cf->diffuseFrame(new TRandom3(0)); // THE MAGIC PART
			hits = cf->findHits();
			tempClusterHitMap = hits->findClustersHitMap();
		}

		else if (dataType == kData) {
			f->getDataFrame(i, cf, energy);
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
	delete f;
	
	// Here it is possible to access and modify the cluster shapes
	// Each cluster is stored as a Hits (Hit collection) pointer in the vector collection clusterHitMap.
	// To loop through each cluster use for (i=0; i<clusterHitMap->size(); i++) { Hits * myCluster = clusterHitMap->at(i); myCluster->....; }
	// E.g. to count the number of hits in each cluster:
	//
	// for (int i=0; i<clusterHitMap->size(); i++) {
	// 	Hits *myCluster = clusterHitMap->at(i);
	//	if (!myCluster) continue;
	//	int counter = 0;
	//	for (int j=0; j<myCluster->GetEntriesFast(); j++) {
	//		Hit *myHit = myCluster->At(j);
	//		if (!myHit) continue;
	//		counter++;
	//	}
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

void makeTracks(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
}

void drawTrackRanges(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	run_energy = energy;
	setWaterPValues(kHigh);
	
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
	tracks->extrapolateToLayer0();

	for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
		Track *thisTrack = tracks->At(j);
		if (!thisTrack) continue;
		Float_t preEnergyLoss = thisTrack->getPreEnergyLoss();

		hFitResults->Fill(getWEPLFromEnergy(thisTrack->getEnergy() + preEnergyLoss));	
	}
}

void drawTungstenSpectrum(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	run_energy = energy;
	setWaterPValues(kHigh);
	
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
	tracks->extrapolateToLayer0();

	for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
		Track *thisTrack = tracks->At(j);
		if (!thisTrack) continue;
		
// 		if (thisTrack->doesTrackEndAbruptly()) {
// 			hFitResults->Fill(getWEPLFromEnergy(thisTrack->getEnergy()));
// 			continue;
// 		}
		
		outputGraph = (TGraphErrors*) thisTrack->doRangeFit();
		if (!outputGraph) continue;
		
		hFitResults->Fill(getWEPLFromEnergy(thisTrack->getFitParameterEnergy()));	
	}
}

void drawScintillatorStatistics(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	run_energy = energy;
	setWaterPValues(kHigh);
	
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
	tracks->extrapolateToLayer0();

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
	if (run_energy < 150) setWaterPValues(kLow);
	else setWaterPValues(kHigh);
	TH1F *hScale = new TH1F("hScale", "Scale histogram", 800, 0, 800);
	TGraphErrors *outputGraph;

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
	
	for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
		Track *thisTrack = tracks->At(j);
		if (!thisTrack) continue;
		
		outputGraph = (TGraphErrors*) thisTrack->doRangeFit();
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
	
// 	setCalibratedTungstenPValues();
	
	Float_t x_energy[sizeOfEventID*400] = {};
	Float_t y_energy[sizeOfEventID*400] = {};
	
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy, x_energy, y_energy);
	tracks->extrapolateToLayer0();
	
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
	Float_t fitEnergy, fitScale = 0;
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
//		if 		(thisTrack->getYmm(0) > 0)		{ preTL = thisTrack->getPreTL() + firstUpperLayerZ; }
//		else if 	(thisTrack->getYmm(0) < 0) 	{ preTL = thisTrack->getPreTL() + firstLowerLayerZ; }
		for (Int_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) { x_energy[i] += preTL; }

		convertXYToWEPL(x_energy, y_energy, eventID);
		realEnergy = getEnergyFromXY(x_energy, y_energy, eventID);

		fitEnergy = thisTrack->getFitParameterEnergy();
		fitScale  = thisTrack->getFitParameterScale();
		
		cout << "realEnergy = " << realEnergy << ", fitEnergy = " << fitEnergy << endl;

		hAccuracy->Fill(realEnergy, fitEnergy);

		if (fitIdx < plotSize) {
			drawIndividualGraphs(cGraph, outputGraph, fitEnergy, fitScale, fitIdx++, eventID, x_energy, y_energy);
		}

		else delete outputGraph;
	}
	cAccuracy->cd();
	hAccuracy->Draw("COLZ");
}

Float_t drawBraggPeakGraphFit(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	run_energy = energy;
	
	if (run_energy < 150) setWaterPValues(kLow);
	else setWaterPValues(kHigh);
	
 	setCalibratedTungstenPValues();
	
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();

	cout << "Using aluminum plate: " << kIsAluminumPlate << endl;
	cout << "Using scintillators: " << kIsScintillator << endl;
	
	Bool_t kDrawHorizontalLines = false;
	Bool_t kDrawVerticalLayerLines = false;
	Bool_t kDrawIndividualGraphs = true;
	Bool_t kUseTrackLength = false;
	
	Bool_t kOutsideScintillatorCross = false;
	Bool_t kInsideScintillatorCross = false;

	char * sDataType = getDataTypeChar(dataType); char * sMaterial = getMaterialChar();
	char * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", energy, sMaterial, sDataType);

	Int_t nPlotX = 4, nPlotY = 4;
	Int_t fitIdx = 0, plotSize = nPlotX*nPlotY;

	TGraphErrors *outputGraph;
	TCanvas *cGraph = new TCanvas("cGraph", "Fitted data points", 1400, 1000);
	TCanvas *cFitResults = new TCanvas("cFitResults", hTitle, 1400, 1000);
	cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
	gStyle->SetPadBorderMode(0); gStyle->SetFrameBorderMode(0);
	gStyle->SetTitleH(0.06); gStyle->SetTitleYOffset(1);
	
	TH1F *hFitResults = new TH1F("fitResult", hTitle, 250, getUnitFromEnergy(0), getUnitFromEnergy(energy*1.2));
	hFitResults->SetLineColor(kBlack); hFitResults->SetFillColor(kGreen-5);

	Float_t finalEnergy = 0;
	Int_t nCutDueToTrackEndingAbruptly = 0;
	
	for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
		Track *thisTrack = tracks->At(j);
		if (!thisTrack) continue;
		
		if (thisTrack->doesTrackEndAbruptly()) {
			hFitResults->Fill(getUnitFromEnergy(thisTrack->getEnergy()));
			nCutDueToTrackEndingAbruptly++;
			continue;
		}
		
		if (kOutsideScintillatorCross) {
			if (thisTrack->isHitOnScintillators()) continue;
		}
		
		if (kInsideScintillatorCross) {
			if (!thisTrack->isHitOnScintillators()) continue;
		}
		
		if (kUseTrackLength)
			outputGraph = (TGraphErrors*) thisTrack->doFit();
		else
			outputGraph = (TGraphErrors*) thisTrack->doRangeFit();
		
		if (!outputGraph) continue;
			
		Float_t fitEnergy = thisTrack->getFitParameterEnergy();
		Float_t fitScale = thisTrack->getFitParameterScale();
		
		hFitResults->Fill(getUnitFromEnergy(fitEnergy));

		if (fitIdx < plotSize && kDrawIndividualGraphs ) {
			drawIndividualGraphs(cGraph, outputGraph, fitEnergy, fitScale, fitIdx++);
		}
		
		else delete outputGraph;
	}
	
	if (!kDrawIndividualGraphs) delete cGraph;
	
	cout << 100 * float(nCutDueToTrackEndingAbruptly) / tracks->GetEntriesFast() << " % of the tracks were cut due to seemingly inelastic nuclear interactions.\n";
	
	cFitResults->cd();

	if (kOutputUnit == kPhysical) hFitResults->SetXTitle("Physical range [mm]");
	else if (kOutputUnit == kWEPL) hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
	else if (kOutputUnit == kEnergy) { hFitResults->SetXTitle("Energy [MeV]"); }
	hFitResults->SetYTitle("Number of protons");

	hFitResults->Draw();
	
	// Draw expected gaussian distribution of results from initial energy
	
	Float_t expectedStraggling = 0, expectedMean = 0, dlayer_down = 0, dlayer = 0;
	Float_t separationFactor = 0.9, nullStraggling = 0;
	
	Float_t sigma_energy = getSigmaEnergy(energy);
	
	expectedMean = getUnitFromEnergy(energy);
	expectedStraggling = getUnitStragglingFromEnergy(energy, sigma_energy);
	nullStraggling = getUnitStragglingFromEnergy(energy, 0);
	
	cout << "OutputUnit is " << kOutputUnit << " and the expected mean value is " << expectedMean 
		 << ". The straggling including / excluding energy variation is " << expectedStraggling << " / " << nullStraggling << ".\n";

		 
	Float_t means[10] = {};
	Float_t sigmas[10] = {};
	
// 	doNLandauFit(hFitResults, means);
	Float_t nGaussianFitRange = doNGaussianFit(hFitResults, means, sigmas);
		
	Int_t nMean = 0;
	for (Int_t i=0; i<10; i++) {
		if (means[i]) nMean++;
	}
	
	cFitResults->Update();
	
	TLine *l;
	if (kDrawVerticalLayerLines) {
		Float_t line_z = 0;
		for (Int_t i=0; i<10; i++) {
			line_z = getWEPLFromTL(getLayerPositionmm(i));
			l = new TLine(line_z, 0, line_z, hFitResults->GetMaximum()*1.05);
			l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw();
		}
	}
	
	TLegend *legend = new TLegend(0.15, 0.7, 0.48, 0.85);
	legend->SetTextSize(0.02);
	legend->AddEntry(hFitResults, "Results from individual track fits", "F");
//  	legend->AddEntry(landau, Form("Fit with E = %.1f MeV and #sigma = %.1f mm ", landau_energy, landau->GetParameter(2)*1.7), "F");
	if (kDrawVerticalLayerLines) legend->AddEntry(l, "Sensor layer positions", "L");
// 	legend->Draw();

	Float_t lowerRange = expectedMean - expectedStraggling;
	Float_t higherRange = expectedMean + expectedStraggling;
	
	Float_t sigma_fwhm = getFWxMInRange(hFitResults, lowerRange, higherRange, 2);
	Float_t sigma_fwtm = getFWxMInRange(hFitResults, lowerRange, higherRange, 3);
	Float_t sigma_fwqm = getFWxMInRange(hFitResults, lowerRange, higherRange, 4);
	
	TPaveStats *ps = (TPaveStats*) cFitResults->GetPrimitive("stats");
	hFitResults->SetBit(TH1::kNoStats);
	ps->AddText(Form("Nominal mean = %.2f", expectedMean));
	ps->AddText(Form("Nominal straggling = %.2f", expectedStraggling));
	
	for (Int_t i=0; i<nMean; i++) {
		ps->AddText(Form("WEPL fit %d = %.2f", i+1, means[i]));
		ps->AddText(Form("Energy fit %d = %.2f", i+1, getEnergyFromUnit(means[i])));
		ps->AddText(Form("Sigma fit %d = %.2f", i+1, sigmas[i]));
	}
		
	if (kOutputUnit == kPhysical) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_physical.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	else if (kOutputUnit == kEnergy) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_energy.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	
	else if (kOutputUnit == kWEPL) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_WEPL.png", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	
  	delete tracks;
	
	return nGaussianFitRange;
}

void writeClusterFile(Int_t Runs, Int_t dataType, Float_t energy) {
	run_energy = energy;
	
	Int_t nClusters = kEventsPerRun * 5 * nLayers;
	Int_t nHits = kEventsPerRun * 50;
	Bool_t kRemoveSmallClusters = true;
	
	Focal *f = new Focal();
	CalorimeterFrame *cf = new CalorimeterFrame();
	Hits * hits = new Hits(nHits);
	Clusters * clusters = new Clusters(nClusters);
	
	for (Int_t i=0; i<Runs; i++) {
		if (dataType == kMC) {
			f->getMCFrame(i, cf);	
			cf->diffuseFrame(new TRandom3(0));
			hits = cf->findHits();
			clusters = hits->findClustersFromHits();
		}
		
		else if (dataType == kData) {
			f->getDataFrame(i, cf, energy);
			hits = cf->findHits();
			clusters = hits->findClustersFromHits();
			
			if (kRemoveSmallClusters) {
				Int_t maxRemoveSize = 2;
				clusters->removeSmallClusters(maxRemoveSize);
			}
		}
	}

	ofstream file("output_all_layers.csv");
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

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
	
	Int_t hSizeX = nx/8;
	Int_t hSizeY = ny/8;

	Float_t angleCut = 5.;
	Float_t x0, y0, theta0;
	Int_t nPoints = 0;

	Float_t fit_energy;

	Int_t ntracks = tracks->GetEntriesFast();
	Int_t nScintillatorHits = 0;

	char * sDataType = getDataTypeChar(dataType); char * sMaterial = getMaterialChar();
	char *title = (char*) Form("2D Projection of data from %.2f MeV proton beam in %s (%s)", energy, sMaterial, sDataType);

	TCanvas *c1 = new TCanvas(title);
	
	TH2F *Frame2D = new TH2F("Frame2D", title, hSizeX, 0, nx*2, hSizeY, 0, ny*2);
	TH2F *normalizeFrame = new TH2F("normalizeFrame", "title", hSizeX, 0, nx*2, hSizeY, 0, ny*2);

	for (Int_t i=0; i<ntracks; i++) {
		Track *thisTrack = tracks->At(i);

		if (!thisTrack) continue;
		
		x0 = thisTrack->getX(0);
		y0 = thisTrack->getY(0);
		theta0 = thisTrack->getSlopeAngle();

		fit_energy = thisTrack->getEnergy();
		if (fit_energy > 180) continue;

		if (thisTrack->isHitOnScintillators())	{
			nScintillatorHits ++;

	//		if (theta0 > angleCut) continue;
		Frame2D->Fill(x0, y0, fit_energy);
		normalizeFrame->Fill(x0, y0);
		}
		nPoints++;

	}

	Frame2D->Divide(normalizeFrame);

	delete tracks;

	cout << "nPoints = " << nPoints << endl;
	cout << "Estimated scintillator hits = " << nScintillatorHits << endl;
	
	Frame2D->Draw("COLZ");
	gStyle->SetOptStat(0);

	// draw lines from 474 to 806
	TLine *l1 = new TLine(474, 0, 474, 1280);
	TLine *l2 = new TLine(806, 0, 806, 1280);
	TLine *l3 = new TLine(0, 474, 1280, 474);
	TLine *l4 = new TLine(0, 806, 1280, 806);
	l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();

}

void saveTracks(Tracks *tracks, Int_t dataType, Float_t energy) {
  	
	// C++ / ROOT has something to learn from Python... ;)
	TString sDataType = (dataType == 0) ? "_MC_" : "_data_";

	TString sEnergy = Form("_%.2fMeV", energy);
	TString fileName = "tracks/tracks";
	TString sMaterial = getMaterialChar();
	fileName.Append(sDataType);
	fileName.Append(sMaterial);
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
  	TString sDataType = (dataType == 0) ? "_MC_" : "_data_";
	TString sEnergy = Form("_%.2fMeV", energy);
	TString fileName = "tracks/tracks";
	TString sMaterial = getMaterialChar();
	fileName.Append(sDataType);
	fileName.Append(sMaterial);
	fileName.Append(sEnergy);
	fileName.Append(".root");
	
	TFile *f = new TFile(fileName);
	if (!f) return 0;
	TTree *T = (TTree*) f->Get("T");
	Tracks * tracks = new Tracks();

	T->GetBranch("tracks")->SetAutoDelete(kFALSE);
	T->SetBranchAddress("tracks",&tracks);

//	cout << "There are " << T->GetEntriesFast() << " entries in TTree.\n";

	T->GetEntry(0);
	
	cout << "There are " << tracks->GetEntriesFast() << " tracks in " << fileName << ".\n";
	
	return tracks;
}

Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy, Float_t *x, Float_t *y) {
	Tracks * tracks;
	
	if (recreate) {
		tracks = getTracks(Runs, dataType, kCalorimeter, energy, x, y);
		saveTracks(tracks, dataType, energy);
	}

	else {
		tracks = loadTracks(Runs, dataType, energy);
	
		if (!tracks) {
			cout << "!tracks, creating new file\n";
			tracks = getTracks(Runs, dataType, kCalorimeter, energy);
			saveTracks(tracks, dataType, energy);
		}
	}
	return tracks;
}

Tracks * getTracks(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy, Float_t *x, Float_t *y) {
	Focal *f = new Focal();

	Int_t nClusters = kEventsPerRun * 5 * nLayers;
	Int_t nHits = kEventsPerRun * 50;
	Int_t nTracks = kEventsPerRun * 2;

	Bool_t breakSignal = false;

	CalorimeterFrame *cf = new CalorimeterFrame();
	TrackerFrame *tf = new TrackerFrame();
	Clusters * clusters = new Clusters(nClusters);
	Clusters * trackerClusters = new Clusters(nClusters);
	Hits *hits = new Hits(nHits);
	Hits *trackerHits = new Hits(nHits);
	Tracks *calorimeterTracks = new Tracks(nTracks);
	Tracks *trackerTracks = new Tracks(nTracks);
	Tracks *allTracks = new Tracks(nTracks * Runs);
	TRandom3 *gRandom = new TRandom3(0);
	
	Int_t eventID = -1;

	TStopwatch t1, t2, t3, t4, t5, t6;
	
	for (Int_t i=0; i<Runs; i++) {

		cout << "Finding track " << (i+1)*kEventsPerRun << " of " << Runs*kEventsPerRun << "... ";
		if (dataType == kMC) {
		
			t1.Start();
			eventID = f->getMCFrame(i, cf, x, y);
			if (kDebug) cout << "Sum frame layer 2 = " << cf->getTH2F(2)->GetSum() << endl;
			t1.Stop(); t2.Start();
			if (kDebug) cout << "Start diffuseFrame\n";
			cf->diffuseFrame(gRandom);
			if (kDebug) cout << "End diffuseFrame, start findHits\n";
			t2.Stop(); t3.Start();
			hits = cf->findHits(eventID);
			if (kDebug) cout << "Number of hits in frame: " << hits->GetEntriesFast() << endl;
			t3.Stop(); t4.Start();
			clusters = hits->findClustersFromHits(); // badly optimized
			clusters->removeSmallClusters(2);
			if (kDebug) cout << "Number of clusters in frame: " << clusters->GetEntriesFast() << endl;
			t4.Stop();
		}
		
		else if (dataType == kData) {
			f->getDataFrame(i, cf, energy);
			hits = cf->findHits();
			clusters = hits->findClustersFromHits();
			clusters->removeSmallClusters(2);
		}
		
		t5.Start();
		calorimeterTracks = clusters->findCalorimeterTracks();
		t5.Stop();

		if (calorimeterTracks->GetEntriesFast() == 0) breakSignal = kTRUE; // to stop running

		if (kUseTrackSplitting) {
			calorimeterTracks->splitSharedClusters();
		}

		// should do track matching here
		// and append calorimeterTracks to trackerTracks...

		for (Int_t j=0; j<calorimeterTracks->GetEntriesFast(); j++) {
			if (!calorimeterTracks->At(j)) continue;

			allTracks->appendTrack(calorimeterTracks->At(j));
		}

		allTracks->appendClustersWithoutTrack(clusters->getClustersWithoutTrack());
		cout << Form("Timing: getMCframe (%.2f sec), diffuseFrame (%.2f sec), findHits (%.2f sec), findClustersFromHits (%.2f sec), findTracks (%.2f sec)\n",
			     t1.RealTime(), t2.RealTime(), t3.RealTime(), t4.RealTime(), t5.RealTime());

		cf->Reset();
//		tf->Reset();
		hits->clearHits();
		trackerHits->clearHits();
		clusters->clearClusters();
		trackerClusters->clearClusters();
		calorimeterTracks->clearTracks();
		trackerTracks->clearTracks();

		if (breakSignal) break;
	}

	delete cf;
//	delete tf;
	delete clusters;
	delete trackerClusters;
	delete hits;
	delete trackerHits;
	delete calorimeterTracks;
	delete trackerTracks;
	delete f;

	return allTracks;
}

void drawTracks3D(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	// FIXME Add kTracker to loadOrCreateTracks
  
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();

	TCanvas *c1 = new TCanvas("c1");
	c1->SetTitle(Form("Tracks from %.2f MeV protons on %s", energy, getMaterialChar()));
	TView *view = TView::CreateView(1);
	view->SetRange(0, 0, 0, 2*nx, nLayers, 2*ny);

	TClonesArray *restPoints = tracks->getClustersWithoutTrack();

	TPolyMarker3D *pMarker = new TPolyMarker3D(restPoints->GetEntriesFast(), 7);
   pMarker->SetMarkerColor(kRed);

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

	for (Int_t i=0; i<ntracks; i++) {
		Track *thisTrack = tracks->At(i);
		if (thisTrack->getTrackLengthmm() < 2) continue;
		
		Int_t n = thisTrack->GetEntriesFast();

		TPolyLine3D *l = new TPolyLine3D(n);
		l->SetLineWidth(1);
		
		// for points
		Int_t pointNumber = 0;

		for (Int_t j=0; j<n; j++) {
			// since we now keep null points here...
			if (!thisTrack->At(j)) continue;

			Float_t x = thisTrack->getX(j);
			Float_t z = thisTrack->getY(j);
			Float_t y = thisTrack->getLayer(j);
			l->SetPoint(pointNumber++,x,y,z);
		}
		l->Draw();
	}
	view->ShowAxis();
   c1->Update();

   TAxis3D *axis = TAxis3D::GetPadAxis();
   axis->SetLabelColor(kBlack);
   axis->SetAxisColor(kBlack);

   tracks->checkLayerOrientation();

   delete tracks;
}

void drawDiffusionCheck(Int_t Runs, Int_t Layer) {
   Focal *f = new Focal();
   CalorimeterFrame *cf = new CalorimeterFrame();
   
   for (Int_t i=0; i<=Runs; i++) {
   	f->getMCFrame(i, cf); // Remember to have MC data available at ./test.root
   }

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

void drawFrame2D(Int_t Runs, Int_t Layer) {
   // Draw one layer using GetFrame2D

   Focal *f = new Focal();
   CalorimeterFrame *cf = new CalorimeterFrame();
   
   TList *histogramList = new TList;
   
   for (Int_t i=0; i<Runs; i++) {
	   f->getMCFrame(i, cf);
	   histogramList->Add(cf->getTH2F(Layer));
   }
 
   TH2F *Frame2D = new TH2F("Frame2D", Form("Hit distribution in layer %i", Layer),
                              nx*2, 0, nx*2, ny*2, 0, ny*2);

   Frame2D->Merge(histogramList);
   Frame2D->Draw("COLZ");
   gStyle->SetOptStat(0);
}

void drawData3D(Int_t Runs) {
   Focal *f = new Focal();

   TH3F *Frame3D = new TH3F("Frame3D", "3D map of energy deposition [keV]", 
                              100, -120, -50, 100, 0, 2*nx, 100, 0, 2*ny);

   Frame3D->SetXTitle("Z axis");
   Frame3D->SetYTitle("X axis"); // to get projection right (Z is depth, not up)
   Frame3D->SetZTitle("Y axis"); 

   for (Int_t run=0; run<Runs; run++) {
   	f->getMCData(run, Frame3D);
   }

   Frame3D->Draw("LEGO");
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
	Int_t chip = (x0 >= nx) + 2 * (y0 < ny);
	Bool_t cutChipNumber = (chip<2) ? true : false;

	return cutChipNumber;
}

Bool_t getCutBraggPeakInTrack(Track *track) {
	Float_t braggPeakRatio = 2.5;

	Int_t lastBin = track->GetEntriesFast() - 1;
	if (lastBin < 4) return false;

	Float_t lowStd = track->getStdSizeToIdx(lastBin-1);
	Int_t nEmptyBins = track->getNMissingLayers();

//	cout << "The standard deviation for the first " << lastBin << " bins is " << lowStd << ".\n";
//	cout << "Number of missing layers = " << nEmptyBins << endl;

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

void getPValues() {
	// create list with energies vs range

	TCanvas *cCustom = new TCanvas("cCustom", "Range fit for custom range-energy list", 1200, 900);
	Float_t energies[6] = {150, 160, 170, 180, 190, 200};
	Float_t ranges[6] = {18.92, 20.95, 23.42, 25.61, 28.27, 30.97};

	TGraph *graph_custom = new TGraph(6, energies, ranges);
	graph_custom->SetTitle("Range fit for custom range-energy list");
	graph_custom->Draw("A*");

	TF1 *fCustom = new TF1("fCustom", "[0]*pow(x, [1])",90,200);
	fCustom->SetParameter(0, 0.01); fCustom->SetParameter(1, 1.7);
	graph_custom->Fit("fCustom");

	Float_t a = fCustom->GetParameter(0);
	Float_t p = fCustom->GetParameter(1);

	cout << Form("For custom, the range fit is R = %.6f * E ^ %.6f\n", a, p);

}

void drawIndividualGraphs(TCanvas *cGraph, TGraphErrors* outputGraph, Float_t fitEnergy, Float_t fitScale, Int_t fitIdx, Int_t eventID, Float_t *x_energy, Float_t *y_energy) {
	cGraph->cd(fitIdx+1);
	Bool_t kDrawFit = true;
	Bool_t kDrawText = true;

	outputGraph->SetMinimum(0);
	outputGraph->SetMaximum(600);
	outputGraph->SetTitle("");
	
	if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		outputGraph->GetXaxis()->SetTitle("Water Equivalent Path Length [mm]");
	}
	
	else if (kOutputUnit == kPhysical) {
		outputGraph->GetXaxis()->SetTitle("Physical path length [mm]");
	}

	outputGraph->GetYaxis()->SetTitle("Deposited energy per layer [keV]");
	outputGraph->GetXaxis()->SetTitleSize(0.05);
	outputGraph->GetYaxis()->SetTitleSize(0.05);
	outputGraph->GetXaxis()->SetLabelSize(0.04);
	outputGraph->GetYaxis()->SetLabelSize(0.04);

	Float_t low = getUnitFromEnergy(0);
	Float_t high = getUnitFromEnergy(run_energy * 1.2);

	outputGraph->GetXaxis()->SetLimits(low, high);

	outputGraph->SetMarkerColor(4);
	outputGraph->SetMarkerStyle(21);
	outputGraph->Draw("AP");

	TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
	func->SetParameters(fitEnergy, fitScale);

	if (kDrawFit) {
		func->Draw("same");
	}
	
	if (x_energy) {
		Float_t WEPLFactor = getWEPLFactorFromEnergy(run_energy);
		Long64_t n=0;
		Long64_t j=0;
		
		for (Long64_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) {
			
			cout << i << ", " << x_energy[i] << ", " << y_energy[i] << endl;

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
		TLatex *text2 = new TLatex(13, 450, Form("'Real' energy: %.1f MeV", realEnergy));
		text2->SetTextSize(0.06);
		text2->Draw();

	}

	if (kDrawText) {
		TLatex *text = new TLatex(10, 500, Form("Fitted energy: %.1f MeV", fitEnergy));
		text->SetTextSize(0.06);
		text->Draw();
	}

	
	cGraph->Update();
}

Float_t  doNGaussianFit ( TH1F *h, Float_t *means, Float_t *sigmas) {
	TF1 *gauss;
	cout << "Energy " << run_energy << endl;
	cout << "Layer;\t Constant;\t Mean;\t\t Energy;\t Sigma;\t\t Fits in layer;\t Chi2/n\n";
	
	Float_t constant, mean, lEnergy, sigma;
	
	Float_t array_mean[3] = {};
	Int_t array_layer[3] = {};
	Float_t array_constant[3] = {};
	Float_t array_sigma[3] = {};
	Float_t array_energy[3] = {};
	Float_t array_f[3] = {};
	Float_t array_chi2n[3] = {};
	
	TAxis *axis = h->GetXaxis();
	Float_t fullIntegral = h->Integral();
	
	Float_t maxBinHeight = h->GetMaximum();

	Bool_t isLastLayer, wasLastLayer = false;;
	
	Int_t j=0;
	for (Int_t i=0; i<15; i++) {
 		if (getWEPLFromTL(getLayerPositionmm(i)) > getUnitFromEnergy(run_energy*1.1)) continue;

		Float_t searchFrom = getWEPLFromTL(getLayerPositionmm(i))+10;
		Float_t searchTo = getWEPLFromTL(getLayerPositionmm(i+1))+10;
		
		Int_t bmin = axis->FindBin(searchFrom);
		Int_t bmax = axis->FindBin(searchTo);
		
		TLine *l1 = new TLine(searchFrom, 0, searchFrom, 1000); l1->SetLineColor(kGreen); l1->Draw();
		TLine *l2 = new TLine(searchTo, 0, searchTo, 1000); l2->SetLineColor(kRed); l2->Draw();
		
		Float_t integral = h->Integral(bmin, bmax);
		Float_t ratio = integral / fullIntegral;
		
		isLastLayer = ((getWEPLFromTL(getLayerPositionmm(i)) > getUnitFromEnergy(run_energy-10)) && !wasLastLayer && ratio > 0.01) ;
		
 		if (i<=3) continue;
 		if (ratio < 0.05 && !isLastLayer) continue;
		
		gauss = new TF1(Form("Gaus_%d", i), "gaus(0)", searchFrom, searchTo);
		
		sigma = 3;
		if (isLastLayer && ratio < 0.05) sigma = 0.2;
		
		gauss->SetParameters(10, (searchFrom+searchTo)/2, sigma);
		gauss->SetParLimits(0, 1, maxBinHeight);
		gauss->SetParLimits(1, searchFrom, searchTo);
		gauss->SetParLimits(2, 2, 15);
		
		h->Fit(gauss, "M, B, WW, Q", "", searchFrom, searchTo);
		
		Float_t chi2 = gauss->GetChisquare();
		Float_t chi2n = chi2 / integral;
		
		sigma = gauss->GetParameter(2);
		constant = gauss->GetParameter(0);
		mean = gauss->GetParameter(1);
		lEnergy = getEnergyFromUnit(mean);
		
 		cout << Form("Searching from %.1f to %.1f, with midpoint at %.1f. Found best fit @ %.1f with chi2 = %.2f and chi2/n = %.2f, ratio = %.2f.\n", searchFrom, searchTo,(searchTo+searchFrom)/2 , mean, chi2, chi2n, ratio);
		
  		if (chi2n > 12) {
  			delete gauss;
  			continue;
		}


 		if (ratio > 0.05 || isLastLayer) {
 			gauss->Draw("same");
			cout << Form("%d;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f;\t %.2f\n", i, constant, mean, lEnergy, sigma, ratio, chi2n);
			
			if (j<3) {
				array_mean[j] = mean;
				array_layer[j] = i;
				array_constant[j] = constant;
				array_energy[j] = lEnergy;
				array_sigma[j] = sigma;
				array_f[j] = ratio;
			}
			
			sigmas[j] = sigma;
			means[j++] = mean;

 		}
 		wasLastLayer = isLastLayer;
	}
	
	Float_t estimated_energy = 0, estimated_range = 0;
	Float_t sum_constant = 0;
	for (Int_t i=0 ; i<3; i++) {
		estimated_range += array_constant[i] * array_mean[i];
		sum_constant += array_constant[i];
	}
	estimated_range /= sum_constant;
	estimated_energy = getEnergyFromUnit(estimated_range);
	cout << "ESTIMATED ENERGY FROM RUN IS " << estimated_energy << endl;
	cout << "Estimated range = " << estimated_range << endl;
	
	if (true) {
		ofstream file("output_gauss.csv", ofstream::out | ofstream::app);
		// run_energy; layer[i], constant[i], mpv[i], energy[i], sigma[i], ratio[i];  i 1->3
		
		file << run_energy << "; ";
		for (Int_t j=0; j<3; j++) {
			file << Form("%d; %.5f; %.5f; %.5f; %.5f; %.5f; %.5f",
						 array_layer[j], array_constant[j], array_mean[j],
						 array_energy[j], array_sigma[j], array_f[j], array_chi2n[j]);
		}
		
		file << endl;
		file.close();
	}
	return estimated_range;
}

void doNLandauFit(TH1F *h, Float_t *mpvs) {
	TF1 *landau;
	cout << "Energy " << run_energy << endl;
	cout << "Layer;\t Constant;\t MPV;\t\t Energy;\t Sigma;\t\t Fits in layer\n";
	Float_t constant, mpv, lEnergy, sigma;
	
	Bool_t kWriteFile = true;
	Float_t array_mpv[3] = {};
	Int_t array_layer[3] = {};
	Float_t array_constant[3] = {};
	Float_t array_sigma[3] = {};
	Float_t array_energy[3] = {};
	Float_t array_f[3] = {};
	
	TAxis *axis = h->GetXaxis();
	Float_t fullIntegral = h->Integral();
	
	Int_t j=0;
	for (Int_t i=0; i<15; i++) {
		// do Landau fit for all possible layers
		
 		if (getWEPLFromTL(getLayerPositionmm(i)) > getUnitFromEnergy(run_energy*1.2)) continue;

		Float_t searchFrom = getWEPLFromTL(getLayerPositionmm(i)) + 10;
		Float_t searchTo = getWEPLFromTL(getLayerPositionmm(i+1)) + 10;
		
		Int_t bmin = axis->FindBin(searchFrom);
		Int_t bmax = axis->FindBin(searchTo);
		
		Float_t integral = h->Integral(bmin, bmax);
		Float_t ratio = integral / fullIntegral;
		
 		if (ratio < 0.05 || i<=4) continue;

		landau = new TF1(Form("Landau_%d", i), "landau(0)", searchFrom, searchTo);
		
  		landau->SetParameters(500, searchFrom+2, 0.5);
  		landau->SetParLimits(0, 75, 2000);
  		landau->SetParLimits(1, searchFrom, searchTo);
  		landau->SetParLimits(2, 0, 2);
		landau->SetNpx(1000);
		
		h->Fit(landau, "Q,M,B", "", searchFrom, searchTo);
		
		sigma = landau->GetParameter(2);
		constant = landau->GetParameter(0);
		mpv = landau->GetParameter(1);
		lEnergy = getEnergyFromUnit(mpv);
		
		cout << Form("Searching from %.1f to %.1f, with midpoint at %.1f. Found best fit @ %.1f.\n", searchFrom, searchTo, searchTo, mpv);
		
 		if (sigma<5 && constant != 250 and ratio > 0.05) {
			landau->Draw("same");
			cout << Form("%d;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f\n", i, constant, mpv, lEnergy, sigma, ratio);
			
			if (j<3) {
				array_mpv[j] = mpv;
				array_layer[j] = i;
				array_constant[j] = constant;
				array_energy[j] = lEnergy;
				array_sigma[j] = sigma;
				array_f[j] = ratio;
			}
			
			mpvs[j++] = mpv;
 		}
	}
	
	if (kWriteFile) {
		ofstream file("output_landau.csv", ofstream::out | ofstream::app);
		// run_energy; layer[i], constant[i], mpv[i], energy[i], sigma[i], ratio[i];  i 1->3
		
		file << run_energy << "; ";
		for (Int_t j=0; j<3; j++) {
			file << Form("%d; %.5f; %.5f; %.5f; %.5f; %.5f; ",
						 array_layer[j], array_constant[j], array_mpv[j],
						 array_energy[j], array_sigma[j], array_f[j]);
		}
		
		file << endl;
		file.close();
	}
}
