#define Analysis_cxx
#include "Analysis.h"
#include "Constants.h"
#include "TFocal.h"
#include "Tracks.h"
#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TEllipse.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <TView.h>
#include <TLeaf.h>
#include <TArrow.h>
#include <TF1.h>
#include "Math/ProbFunc.h"

using namespace std;

/*
void CompareDataAndMC(Int_t MCRuns, Int_t DataRuns) {

	Focal f;

	// Get tracks from MC
	Int_t nEvents = 200;
	Int_t energy = 190;

	Tracks *allTracksMC = new Tracks(2*nEvents * MCRuns);
	allTracksMC->SetOwner(kTRUE);

	for (Int_t i=0; i<MCRuns; i++) {
		Int_t EventFrom = MCRuns * i;
		Int_t EventTo = EventFrom + nEvents;
		Tracks *tracks = f.FindTracks(EventFrom, EventTo);
		tracks->SetOwner(kTRUE);
		cout << "--- Run " << i+1 << " of " << MCRuns << " complete --- \n";
		for (Int_t j=0; j<tracks->GetEntriesFast(); j++) allTracksMC->appendTrack(tracks->At(j));
		delete tracks;
	}
	cout << "From " << MCRuns << " runs with " << nEvents << " events in each, we have collected "
			<< allTracksMC->GetEntriesFast() << " tracks.\n";

	// Get tracks from data
	Clusters *restPoints = new Clusters(DataRuns * 20);
	Tracks *allTracksData = f.FindRealTracks(DataRuns, restPoints, 190);
	allTracksData->SetOwner(kTRUE);

	// Plot statistics
	// First: Cluster sizes per layer
	
	Float_t small = 1e-5;
	
	TCanvas *c1 = new TCanvas("c1", "c1", 1200, 900);
	gStyle->SetStatFont(22);
	gStyle->SetTitleFont(22);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	c1->Divide(7, 2, small, small); // nlayers & Data/MC

	// Then: Proton angles
	TCanvas *c2 = new TCanvas("c2", "c2", 1200, 900);
	gStyle->SetStatFont(22);
	gStyle->SetTitleFont(22);
	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	c2->Divide(7, 2, small, small); // x: nlayers, y: data / MC

	vector<TH1F*> *hCSMC = new vector<TH1F*>;
	vector<TH1F*> *hCSData = new vector<TH1F*>;

	vector<TH1F*> *hAngleMC = new vector<TH1F*>;
	vector<TH1F*> *hAngleData = new vector<TH1F*>;

	hCSMC->reserve(7); hCSData->reserve(7);
	hAngleMC->reserve(7); hAngleData->reserve(7);

	for (Int_t layer=0; layer<7; layer++) {
		hCSMC->push_back(new TH1F(Form("hCSMC_%i", layer), 
					Form("Cluster sizes layer %i (MC)", layer), 60, 0, 60));
		hCSData->push_back(new TH1F(Form("hCSData_%i", layer), 
					Form("Cluster sizes layer %i (Data)", layer), 60, 0, 60));

		hCSMC->back()->SetXTitle("Cluster size [# of pixels]");
		hCSMC->back()->SetYTitle("Number of hits");
		hCSMC->back()->GetYaxis()->SetTitleOffset(1.5);
		hCSData->back()->SetXTitle("Cluster size [# of pixels]");
//		hCSData->back()->SetYTitle("Number of hits");
//		hCSData->back()->GetYaxis()->SetTitleOffset(1.4);

		hAngleMC->push_back(new TH1F(Form("hAngleMC_%i", layer), 
					Form("Proton angles from layer 0 to %i (MC)", layer), 100, 0, 25));
		hAngleData->push_back(new TH1F(Form("hAngleData_%i", layer), 
					Form("Proton angles from layer 0 to %i (Data)", layer), 100, 0, 25));

		hAngleMC->back()->SetXTitle("Proton angle [deg]");
		hAngleMC->back()->SetYTitle("Number of hits");
		hAngleMC->back()->GetYaxis()->SetTitleOffset(1.5);
		hAngleData->back()->SetXTitle("Proton angle [deg]");
//		hAngleData->back()->SetYTitle("Number of hits");
//		hAngleData->back()->GetYaxis()->SetTitleOffset(1.4);
	}
	
	// fill histogram with data
	Bool_t cut;
	Int_t minTL;

	Bool_t kCutTL = kFALSE;
	Bool_t kCutChip = kFALSE;
	Int_t layer, cs;
	Float_t TL, x0, y0;
	Float_t angle;

	if (energy == 190) minTL = 20;
	if (energy == 180 || energy == 170) minTL = 15;
	if (energy == 160 || energy == 150) minTL = 10;
	if (energy == 140 || energy == 130 || energy == 122) minTL = 5;

	// loop through tracks
	for (Int_t i=0; i<allTracksData->GetEntriesFast(); i++) {

		if (!allTracksData->At(i)) continue;

		cut = kTRUE;
		TL = allTracksData->getTrackLengthmm(i);
		x0 = allTracksData->At(i)->getX(0);
		y0 = allTracksData->At(i)->getY(0);

		if (kCutTL)		cut *= (TL>minTL); // cut on minimum track length
		if (kCutChip)	cut *= (y0>=ny); // cut on chip NW or NE

		// loop through clusters
		for (Int_t j=0; j<allTracksData->GetEntriesFast(i); j++) {
			layer = allTracksData->At(i)->getLayer(j);
			if (layer<7) {
				cs = allTracksData->At(i)->getSize(j);
				angle = allTracksData->At(i)->getSlopeAngleAtLayer(j);
				
				hCSData->at(layer)->Fill(cs);
				hAngleData->at(layer)->Fill(angle);
			}
		}
	}

	// fill histogram with MC
	
	for (Int_t i=0; i<allTracksMC->GetEntriesFast(); i++) {
		cut = kTRUE;
		if (!allTracksMC->At(i)) continue;
		TL = allTracksMC->getTrackLengthmm(i);

		if (kCutTL)		cut *= (TL>minTL); // cut on minimum track length

		// loop through clusters
		for (Int_t j=0; j<allTracksMC->GetEntriesFast(i); j++) {
			if (!allTracksMC->At(i)) continue;

			layer = allTracksMC->At(i)->getLayer(j);
			if (layer<7) {
				cs = allTracksMC->At(i)->getSize(j);
				angle = allTracksMC->At(i)->getSlopeAngleAtLayer(j);

				hCSMC->at(layer)->Fill(cs);
				hAngleMC->at(layer)->Fill(angle);
			}
		}
	}


	// plot everything
	TH1F *h;
	Float_t left = 2e-1;
	for (Int_t layer=0; layer<7; layer++) {

		h = hCSMC->at(layer);
		c1->cd(layer+1);
			h->SetFillColor(kGreen-5+layer);
			h->SetLineColor(kBlack);
			gPad->SetBottomMargin(small);
			gPad->SetRightMargin(small);
			if (layer>0) gPad->SetLeftMargin(small);
			if (layer==0) gPad->SetLeftMargin(left);
			h->Draw();
		
		h = hCSData->at(layer);
		c1->cd(layer+8);
			h->SetFillColor(kRed-5+layer);
			h->SetLineColor(kBlack);
			gPad->SetTopMargin(small);
			gPad->SetRightMargin(small);
			if (layer>0) gPad->SetLeftMargin(small);
			if (layer==0) gPad->SetLeftMargin(left);
			h->Draw();

		h = hAngleMC->at(layer);
		c2->cd(layer+1);
			h->SetFillColor(kGreen-5+layer);
			h->SetLineColor(kBlack);
			gPad->SetBottomMargin(small);
			gPad->SetRightMargin(small);
			if (layer>0) gPad->SetLeftMargin(small);
			if (layer==0) gPad->SetLeftMargin(left);
			h->Draw();

		h = hAngleData->at(layer);
		c2->cd(layer+8);
			h->SetFillColor(kRed-5+layer);
			h->SetLineColor(kBlack);
			gPad->SetTopMargin(small);
			gPad->SetRightMargin(small);
			if (layer>0) gPad->SetLeftMargin(small);
			if (layer==0) gPad->SetLeftMargin(left);
			h->Draw();
	}
}


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

void GetTrackStatistics(Int_t Events, Int_t Runs) {
	
	Focal f;

	Tracks *allTracks = new Tracks(2*Events*Runs);

	allTracks->SetOwner(true);

	for (Int_t i=0; i<Runs; i++) {
		Int_t EventFrom = Runs * i;
		Int_t EventTo = EventFrom + Events;

		// could use absorbobjects(tracks), but i don't get it to work
		// Track objects are small though.
		
		Tracks *tracks = f.FindTracks(EventFrom, EventTo);
		tracks->SetOwner(true);

		cout << "--- Run " << i+1 << " of " << Runs << " complete --- \n";

		for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
			allTracks->appendTrack(tracks->At(j));
		}
		delete tracks;
	}

	cout << "From " << Runs << " runs with " <<Events << " events in each, we have collected"
			<< allTracks->GetEntriesFast() << " tracks.\n";

	// now we can calulate some statistics, such as
	//  - average track length
	//  - position of bragg peak ( with track -> cluster -> size )
	//  - make histogram of naive track lengths (number of particles per layer)
	//			versus realistic track lengths
	//	 - make histogram of naive size distributions along track vs realistic
	
	TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
	TCanvas *c4 = new TCanvas("c4", "c4", 800, 600);

	TH1F *hTrackLengths = new TH1F("hTrackLengths", "Track Lengths (MC)", 100, 0, 120);
	TH2F *hClusterSizeAlongTrack = new TH2F("hClusterSizeAlongTrack",
				"Cluster size along track length (MC)", 1.5*nLayers, 0, 1.5*nLayers*dz, 50, 0, 50);
	TH1F *hStraightness = new TH1F("hStraightness", "Sinuosity plot (MC)", 500, 1, 1.01);
	TH1F *hSlope = new TH1F("hSlope", "Proton angle plot (MC)", 500, 0, 20);
	
	hTrackLengths->SetXTitle("Track length [mm]");
	hClusterSizeAlongTrack->SetXTitle("Track length [mm]");
	hClusterSizeAlongTrack->SetYTitle("Cluster size [# of pixels]");
	hStraightness->SetXTitle("Sinuosity parameter");
	hSlope->SetXTitle("Total track angle (degree)");

	Float_t trackLengthSoFar = 0;
	for (Int_t i=0; i<allTracks->GetEntriesFast(); i++) {
		hTrackLengths->Fill(allTracks->getTrackLengthmm(i));
		hStraightness->Fill(allTracks->getSinuosity(i));
		hSlope->Fill(allTracks->getSlopeAngle(i));
		for (Int_t j=0; j<allTracks->GetEntriesFast(i); j++) {
			trackLengthSoFar += allTracks->At(i)->getTrackLengthmmAt(j);
			hClusterSizeAlongTrack->Fill(trackLengthSoFar, allTracks->At(i)->getSize(j));
			cout << "GetSize for cluster at " << *allTracks->At(i)->At(j) << ": " << allTracks->At(i)->getSize(j) << endl;
		} // end loop over clusters
		trackLengthSoFar = 0;
	} // end loop over tracks

	c1->cd();
		hTrackLengths->Draw();
	c2->cd();
		gStyle->SetOptStat(0);
		hClusterSizeAlongTrack->Draw("COLZ");
	c3->cd();f
		hStraightness->Draw();
	c4->cd();
		hSlope->Draw();

	delete allTracks;
}



// Kristian
// - range resolution for eksperimentelle data
// - antall partikler vs rekkevidde

void DrawClusterShapes() {
	// get vector of TH2F's, each with a hits distribution and cluster size
	
	Focal f;
	
	TCanvas *c = new TCanvas("c", "c", 1000, 800);

	Int_t nRows = 12;
	Int_t nRepeats = 20;
	Int_t nN = nRows * nRepeats;

	c->Divide(nRepeats, nRows, 0.0001, 0.0001);

	vector<TH2F*> *hClusterMaps = new vector<TH2F*>;
	hClusterMaps->reserve(nN);
	for (Int_t i=0; i<nN; i++) {
		hClusterMaps->push_back(new TH2F(Form("hClusterMap_%i",i),
								"", 11, 0, 11, 11, 0, 11));
	}

	vector<Hits*> *clusterHitMap = new vector<Hits*>;
	clusterHitMap->reserve(10000);

	Bool_t MonteCarlo = kFALSE;

	f.GetClusterFrames(clusterHitMap, MonteCarlo);
	
	// fill hClusterMaps with cluster shapes from clusterHitMap
	// sizes 3-5 in first row
	// 6-8 in second row
	// 9-11 in third row
	// 12-14 in fourth row
	// 15-17 in fifth row
	// 18-20 in sixth row
	// 21-23 in seventh row
	// 24-26 in eighth row
	// 27-29 in ninth row
	// 30-32 in tenth row
	// 33-35 in eleventh row
	
	Int_t size_from, size_to, nInRow, x, y;

	Int_t nTotal = 0;
	for (Int_t i=0; i<nRows; i++) {
		size_from = i*3;
		size_to = size_from + 2;
		nInRow = 0;

		for (UInt_t j=0; j<clusterHitMap->size(); j++) {
			Int_t csize = clusterHitMap->at(j)->GetEntriesFast();
			if (csize >= size_from && csize <= size_to) {
				// plot it in the j'th row
				for (Int_t k=0; k<clusterHitMap->at(j)->GetEntriesFast(); k++) {
					x = clusterHitMap->at(j)->getX(k);
					y = clusterHitMap->at(j)->getY(k);
					hClusterMaps->at(nTotal)->Fill(x,y);
				} // end loop through all points
				hClusterMaps->at(nTotal)->SetTitle(Form("Cluster size [%i, %i]", size_from, size_to));
				nInRow++; nTotal++;
				if (nInRow >= nRepeats) break; // stop looping over clusters now
			} // end plot point
		} // end loop over all clusters
		if (nInRow < nRepeats) {
			cout << "Only " << nInRow << " clusters with size [" << i*3 << "," << i*3+2 << "] with i = " << i << ". Setting nTotal from " << nTotal << " to " << (i+1)*nRepeats << ".\n";
			nTotal = (i+1)*nRepeats;
		}
	} // end loop over row with predefined cluster size window

	// draw the canvas
	
	for (Int_t i=0; i<nN; i++) {
		c->cd(i+1);
		hClusterMaps->at(i)->Draw("same, COL,ah,fb,bb");
		gStyle->SetOptStat(0);
	}
}

void DrawBraggPeakFit(Int_t Runs, Int_t energy) {

	Focal f;

//	Int_t energy = 190;

	Clusters *restPoints = new Clusters(Runs*200);
	Tracks *allTracks = f.FindRealTracks(Runs, restPoints, energy);
	allTracks->SetOwner(true);

//	TCanvas *cFollowTrack = new TCanvas("c1", "c1", 1000, 800);
//	cFollowTrack->Divide(4, 4, 0.01, 0.01, 0);

	vector<TH1F*> *hFollowTrackCut = new vector<TH1F*>;
	hFollowTrackCut->reserve(allTracks->GetEntriesFast());

	vector<TH1F*> *hFollowTrack = new vector<TH1F*>;
	hFollowTrack->reserve(allTracks->GetEntriesFast());
	for (Int_t track=0; track<allTracks->GetEntriesFast(); track++) {
		// why doesn't the titles appear?
		hFollowTrack->push_back(new TH1F(Form("hFollowTrack_%i", track), "Cluster size along track length for a single track", nLayers, 0, nLayers*dz));
		hFollowTrack->at(track)->SetXTitle("Track Length [mm]");
		hFollowTrack->at(track)->SetYTitle("Cluster size [# of pixels]");
		
		hFollowTrackCut->push_back(new TH1F(Form("hFollowTrackCut_%i", track), "Cluster size along track length for a single track", nLayers, 0, nLayers*dz));
		hFollowTrackCut->at(track)->SetXTitle("Track Length [mm]");
		hFollowTrackCut->at(track)->SetYTitle("Cluster size [# of pixels]");
	}

	Float_t trackLengthSoFar = 0;
	Int_t trackNum = 0;
	Int_t chip = 0; // quadrant to plot
	Bool_t cutTL;
	Bool_t cutCHIP;
	Bool_t cutBP;
	Bool_t cutNTBP; // bin next-to-bp must be filled
	Bool_t cutCombined;
	
	Int_t okCHIP = 0;
	Int_t okTL = 0;
	Int_t okBP = 0;
	Int_t okALL = 0;

	vector<TArrow*> arrows;
	arrows.reserve(100);

	for (Int_t i=0; i<allTracks->GetEntriesFast(); i++) {

		Float_t TL = allTracks->getTrackLengthmm(i);
		Int_t x0 = allTracks->At(i)->getX(0);
		Int_t y0 = allTracks->At(i)->getY(0);

		// Find the largest bin of the last three
		// This will probably be where the bragg peak is
		Int_t BP = 0;
		Int_t BPidx = -1;
		Int_t BPidxFrom = allTracks->GetEntriesFast(i) - 4;
		if (BPidxFrom < 0) BPidxFrom = 0;
		
		cutNTBP = kTRUE;

		for (Int_t j=BPidxFrom; j<allTracks->GetEntriesFast(i); j++) {
			if (BP <= allTracks->At(i)->getSize(j)) {
				BP = allTracks->At(i)->getSize(j);
				BPidx = j;
			}
			if (allTracks->At(i)->getSize(j) < 0.5) cutNTBP = kFALSE;
		}

		if (!cutNTBP) {
			cout << "!cutNTBP for i = " << i << ". GetSizes follows.\n";
			for (Int_t j=0; j<allTracks->GetEntriesFast(i); j++) {
				cout << allTracks->At(i)->getSize(j) << endl;
			}
		}


		if (BP > 20) {
			cutBP = kTRUE;
			okBP++;
		}
		else cutBP = kFALSE;

		// Is the bragg peak located in the last three bins?
		Int_t nextToMaxBP = 0;
		for (Int_t j=BPidxFrom; j<allTracks->GetEntriesFast(i); j++) {
			if (j==BPidx) continue;
			if (nextToMaxBP <= allTracks->At(i)->getSize(j)) {
				nextToMaxBP = allTracks->At(i)->getSize(j);
			}
		}
		
		Bool_t BPSampled;
		if ((Float_t) BP / nextToMaxBP >= 1.2 && BP > 20) {
			BPSampled = kTRUE;
		}
		else BPSampled = kFALSE;

		Int_t minTL;

		if (energy == 190) minTL = 20;
		if (energy == 180) minTL = 15;
		if (energy == 170) minTL = 15;
		if (energy == 160) minTL = 10;
		if (energy == 150) minTL = 10;
		if (energy == 140) minTL = 5;
		if (energy == 130) minTL = 5;
		if (energy == 122) minTL = 5;

		if (TL>minTL) {
			cutTL = kTRUE;
			okTL++;
		}
		else cutTL = kFALSE;
				
		if (x0<nx && y0>=ny) chip = 0;
		else if (x0>=nx && y0>=ny) chip = 1;
		else if (x0<nx && y0<ny) chip = 2;
		else if (x0>=nx && y0<ny) chip = 3;

		if (chip == 0 || chip == 1) {
			cutCHIP = kTRUE;
			okCHIP++;
		}
		
		else cutCHIP = kFALSE;
		
		cutCombined = cutCHIP * cutTL;// * BPSampled;
		if (cutCombined) okALL++;

		// start comment here

		for (Int_t j=0; j<allTracks->GetEntriesFast(i); j++) {
			trackLengthSoFar += allTracks->At(i)->GetTrackLengthAtmm(j);
			if (trackNum < 16 && cutCombined) {
				hFollowTrack->at(trackNum)->Fill(trackLengthSoFar, allTracks->At(i)->GetSize(j));
	//				hFollowTrack->at(trackNum)->SetTitle(Form("Track length histogram for run %i", i));
				if (BPSampled && BPidx == j) {
					// draw arrow on top of bin
					// containing the bragg peak
					arrows.push_back(new TArrow(trackLengthSoFar+dz/2, allTracks->At(i)->GetSize(j) * 1.3, 
										  trackLengthSoFar+dz/2, allTracks->At(i)->GetSize(j) * 1.1, 0.005, "|>"));
					arrows.back()->SetLineColor(kRed);
					arrows.back()->SetLineWidth(2.);
					arrows.back()->SetFillColor(kRed);
					cout << "BP ok for trackNum " << trackNum << endl;
				}
				else if (!BPSampled && BPidx == j) {
					arrows.push_back(NULL);
				}
			}
			else if (cutCombined) {
				// add to hFollowTrack anyway
				// but don't plot it
				hFollowTrack->at(trackNum)->Fill(trackLengthSoFar, allTracks->At(i)->GetSize(j));
			}

		} // end loop over clusters

		// end comment here

		for (Int_t j=0; j<allTracks->GetEntriesFast(i); j++) {
			trackLengthSoFar += allTracks->At(i)->getTrackLengthmmAt(j);
			if (cutCombined) {
				hFollowTrack->at(trackNum)->Fill(trackLengthSoFar, allTracks->At(i)->getSize(j));

				if (BPSampled) {
					hFollowTrackCut->at(trackNum)->Fill(trackLengthSoFar, allTracks->At(i)->getSize(j));
				}
			}
		}

		trackLengthSoFar = 0;
		if (cutCombined) trackNum++;
	} // end loop over tracks

	cout << "Total tracks: " << allTracks->GetEntriesFast() << endl;
	cout << "Total okCHIP: " << okCHIP << endl;
	cout << "Total okTL: " << okTL << endl;
	cout << "Total okBP: " << okBP << endl;
	cout << "Total ok combined: " << okALL << endl;

	// start comment here
	// Fill histogram hFollowTrack with 16 first tracks
	for (Int_t track=0; track<16; track++) {
		cFollowTrack->cd(track+1);
		gPad->DrawFrame(0, 0, 50, 45);
		hFollowTrack->at(track)->SetFillColor(kBlue-2);
		hFollowTrack->at(track)->Draw("same");
		if (track < arrows.size()) {
			if (arrows.at(track)) {
				arrows.at(track)->Draw();
			}
		}
	}
	// end comment here

	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	TH1F *hMu = new TH1F("hMu", "Mean fit value", 400, 0, 30);

	TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
	TH1F *hTotal = new TH1F("hTotal", "Sum of clustersizes in last two layers", 150, 0, 150);
	TH1F *hTotalCut = new TH1F("hTotalBPcut", "Size of last two bins, with cuts", 150, 0, 150);

	Float_t cog;

	for (UInt_t track=0; track<hFollowTrack->size(); track++) {
// start comment here
 *   BRAGG PEAK GAUSSIAN FIT ALGORITHM
 *   NOT SUITED FOR THIS PROBLEM ...
 
		
 TF1 *g1 = new TF1("m1", "gaus", 15, 35);
		
		Float_t BPpos = hFollowTrack->at(track)->GetMaximumBin() * dz;
		Float_t BP = hFollowTrack->at(track)->GetBinContent(BPpos/dz);
		Float_t BPwidth = 2;

		g1->SetParameters(BP, BPpos, BPwidth); // BP, BPpos, 2.5
		g1->SetParLimits(0, BP*0.8, BP*1.2);
		g1->SetParLimits(1, BPpos*0.9, BPpos*1.1);
		g1->SetParLimits(2, BPwidth, BPwidth);

		cFollowTrack->cd(track+1);

		hFollowTrack->at(track)->Fit("m1", "R, WW");
	// end comment here


		// CoG for last two bins...
		cog = 0;
		Int_t hLen = hFollowTrack->at(track)->FindLastBinAbove();
		if (hLen<3) continue;

		// pen-penultimate, penultimate and ultimate bin
		Float_t ppu = hFollowTrack->at(track)->GetBinContent(hLen-2);
		Float_t pu = hFollowTrack->at(track)->GetBinContent(hLen-1);
		Float_t u = hFollowTrack->at(track)->GetBinContent(hLen);

		// sometimes penultimate channel is bad, in that case fix it
		if (pu<1 && ppu > pu*5) pu = ppu;
//		if (pu<1) continue;
	
		// if the BP is not located in last layer, "push" the cog a bit
//		if (u<21) u += (ppu+pu)/4;

		// cog for last two bins
		cog = (pu * (hLen-1)*dz + u * (hLen)*dz) / (u + pu);
		
		// using cumulative Gaussian as LUT
		Float_t cumGaus = 0;
		Float_t sumLimit = 45; // sum under whole gaussian
		Float_t mu = 30;
		Float_t sigma = dz;

		Float_t limit;
		if (u<25) { // no bragg peak
			limit = u;
		}

		else {
			limit = pu + u;
		}

		while (sumLimit*cumGaus < limit) {
			cumGaus = ROOT::Math::normal_cdf(u, sigma, mu);
			mu -= 0.05;
			if (mu<10) break;
		}

//		if (cog<22) continue;
		hMu->Fill(cog);
		hTotal->Fill(u);

	}

	// loop over track with extra cuts
	for (UInt_t track=0; track<hFollowTrackCut->size(); track++) {
		Int_t hLen = hFollowTrackCut->at(track)->FindLastBinAbove();
		if (hLen<3) continue;

		Float_t ppu = hFollowTrackCut->at(track)->GetBinContent(hLen-2);
		Float_t pu = hFollowTrackCut->at(track)->GetBinContent(hLen-1);
		Float_t u = hFollowTrackCut->at(track)->GetBinContent(hLen);
		if (pu<1 && ppu > pu*5) pu = ppu;
		hTotalCut->Fill(u);
	}

	c2->cd();
	hMu->SetFillColor(kGray);
	hMu->SetTitle(Form("Bragg peak positions for %i MeV protons", energy));
	hMu->Draw();

	// save hMu to root file
	c2->Print(Form("bp_root/%iMeV.root", energy));

	c3->cd();
	hTotal->SetFillColor(kBlue-3);
	hTotal->SetFillStyle(3002);
	hTotal->Draw();
	hTotalCut->SetFillColor(kGreen-5);
//	hTotalCut->Draw("SAME");
	hTotalCut->SetFillStyle(3002);

	delete allTracks;
}

void GetRealTrackStatistics(Int_t Runs) {
	
	Focal f;

	Clusters *restPoints = new Clusters(10000);
	Tracks *allTracks = f.FindRealTracks(Runs, restPoints);
	allTracks->SetOwner(true);

	// now we can calulate some statistics, such as
	//  - average track length
	//  - position of bragg peak ( with track -> cluster -> size )
	//  - make histogram of naive track lengths (number of particles per layer)
	//			versus realistic track lengths
	//	 - make histogram of naive size distributions along track vs realistic

	Int_t nTracksToPlot = 25;
	Int_t nTracksToPlot1D = 5;

	TCanvas *c1d = new TCanvas("c1d", "c1", 800, 600);
	TCanvas *c2d = new TCanvas("c2d", "c2", 800, 600);
	TCanvas *c3d = new TCanvas("c3d", "c3", 800, 600);
	TCanvas *c4d = new TCanvas("c4d", "c4", 800, 600);
	TCanvas *c5d = new TCanvas("c5d", "c5", 800, 600);
	c5d->Divide(2, 2, 0.01, 0.01, 0);

	TCanvas *c6d = new TCanvas("c6d", "c6", 1000, 800);
	c6d->Divide(3, 3, 0.01, 0.01, 0);

	TCanvas *c7d = new TCanvas("c7d", "c7", 1000, 800);
	c7d->Divide(nTracksToPlot1D, nTracksToPlot1D, 0.01, 0.01, 0);

	TCanvas *c8d = new TCanvas("c8d", "c8", 1000, 800);
	c8d->Divide(3,3,0.01,0.01,0);

	TH1F *hTrackLengths = new TH1F("hTrackLengthsd", "Track Lengths (data)", nLayers, 0, nLayers*dz);
	TH2F *hClusterSizeAlongTrack = new TH2F("hClusterSizeAlongTrackd",
				"Cluster size along track length (data)", nLayers, 0, nLayers*dz, 50, 0, 50);
	TH1F *hStraightness = new TH1F("hStraightnessd", "Sinuosity plot (data)", 50, 1, 1.11);
	TH1F *hSlope = new TH1F("hSloped", "Proton angle plot (data)", 50, 0, 20);

	vector<TH1F*> *hAvgCS = new vector<TH1F*>;
	hAvgCS->reserve(4);
	for (Int_t chip=0; chip<4; chip++) {
		hAvgCS->push_back(new TH1F("hAvgCS", Form("Average Cluster Size vs Track Length for chip %i (data)",chip), 50, 0, 50));
	}

	vector<TH1F*> *hCSLayer = new vector<TH1F*>;
	hCSLayer->reserve(9);
	for (Int_t layer=0; layer<9; layer++) {
		hCSLayer->push_back(new TH1F(Form("hCSLayer%id", layer), Form("Cluster Size for layer %i", layer), 50, 0, 50));
	}

	vector<TH1F*> *hFollowTrack = new vector<TH1F*>;
	hFollowTrack->reserve(nTracksToPlot);
	for (Int_t track=0; track<nTracksToPlot; track++) {
		hFollowTrack->push_back(new TH1F(Form("hFollowTrack_%i", track), "Cluster size along track length for a single track", 50, 0, 50));
		hFollowTrack->at(track)->SetXTitle("Track Length [mm]");
		hFollowTrack->at(track)->SetYTitle("Cluster size [# of pixels]");
	Int_t trackNum = 0;
	}
	
	vector<TH1F*> *hAngles = new vector<TH1F*>;
	hAngles->reserve(9);
	for (Int_t layer=0; layer<9; layer++) {
		hAngles->push_back(new TH1F(Form("hANgles_%i", layer), Form("Proton angle distribution in layer %i", layer), 50, 0, 20));
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
	Int_t chip = 0; // quadrant to plot
	Bool_t cutTL;
	Bool_t cutCHIP;
	
	Int_t okCHIP = 0;
	Int_t okTL = 0;
	Int_t okBOTH = 0;

	for (Int_t i=0; i<allTracks->GetEntriesFast(); i++) {
		Float_t TL = allTracks->getTrackLengthmm(i);
		Int_t x0 = allTracks->At(i)->getX(0);
		Int_t y0 = allTracks->At(i)->getY(0);
	
		if (TL>22) {
			cutTL = kTRUE;
			okTL++;
		}

		else cutTL = kFALSE;
				
		if (x0<nx && y0>=ny) chip = 0;
		else if (x0>=nx && y0>=ny) chip = 1;
		else if (x0<nx && y0<ny) chip = 2;
		else if (x0>=nx && y0<ny) chip = 3;

		if (chip == 0 || chip == 1) {
			cutCHIP = kTRUE;
			okCHIP++;
		}
		else cutCHIP = kFALSE;
		
		if (cutTL && cutCHIP) {
			okBOTH++;
		}
		
		hTrackLengths->Fill(TL);
		hStraightness->Fill(allTracks->getSinuosity(i));
		hSlope->Fill(allTracks->getSlopeAngle(i));

		for (Int_t j=0; j<allTracks->GetEntriesFast(i); j++) {
			trackLengthSoFar += allTracks->At(i)->getTrackLengthmmAt(j);
			if (cutTL && cutCHIP) hClusterSizeAlongTrack->Fill(trackLengthSoFar, allTracks->At(i)->getSize(j));
			if (cutTL) {
				// find quadrant to plot from
				hAvgCS->at(chip)->Fill(trackLengthSoFar, allTracks->At(i)->getSize(j));
			}
			Int_t layer = allTracks->At(i)->getLayer(j);
			if (layer < 9) {
				hCSLayer->at(layer)->Fill(allTracks->At(i)->getSize(j));
			}
			if (layer<9) {
				hAngles->at(layer)->Fill(allTracks->At(i)->getSlopeAngleAtLayer(j));
				// slope angle between track at layer 0 and layer j
			}

			if (trackNum < nTracksToPlot && cutTL && cutCHIP) { // don't plot the shortest tracks
				hFollowTrack->at(trackNum)->Fill(trackLengthSoFar, allTracks->At(i)->getSize(j));
				hFollowTrack->at(trackNum)->SetTitle(Form("Track length histogram for run %i", i));
//				cout << "hFollowTrack_" << trackNum << ": Run " << i << endl;
			}
		} // end loop over clusters
		trackLengthSoFar = 0;
		if (cutTL && cutCHIP) trackNum++;
	} // end loop over tracks

	cout << "Total tracks: " << allTracks->GetEntriesFast() << endl;
	cout << "Total okCHIP: " << okCHIP << endl;
	cout << "Total okTL: " << okTL << endl;
	cout << "Total ok(CHIP+TL): " << okBOTH << endl;

	c1d->cd();
		hTrackLengths->Draw();
	c2d->cd();
		gStyle->SetOptStat(0);
		hClusterSizeAlongTrack->Draw("COLZ");
	c3d->cd();
		hStraightness->Draw();
	c4d->cd();
		hSlope->Draw();
	for (Int_t chip=0; chip<4; chip++) {
		c5d->cd(chip+1);
		hAvgCS->at(chip)->SetFillColor(kRed-chip*2);
		hAvgCS->at(chip)->Draw();
	}
	for (Int_t layer=0; layer<9; layer++) {
		c6d->cd(layer+1);
		hCSLayer->at(layer)->SetFillColor(kRed-9+layer);
		hCSLayer->at(layer)->Draw("same");
	}

	for (Int_t track=0; track<nTracksToPlot; track++) {
		c7d->cd(track+1);
		gPad->DrawFrame(0, 0, 50, 35);
		hFollowTrack->at(track)->SetFillColor(kBlue-2);
		hFollowTrack->at(track)->Draw("same");
	}

	for (Int_t layer=0; layer<9; layer++) {
		c8d->cd(layer+1);
		hAngles->at(layer)->SetFillColor(kRed-9+layer);
		hAngles->at(layer)->Draw("same");
	}

	delete allTracks;
}
	
void WriteClusterFile(Int_t Runs) {
	// loop through hits in histogram and make a TTree branch of them
	// Find all hits with neighbours and return a TTree of clusters 
	// with size + center of gravity

	Focal f;

	vector<TH2F*> *Frame3D = new vector<TH2F*>;
	f.GetFrame3D(0, Runs, Frame3D);

	ofstream file("output_all_layers.csv");
	file << "layer;x;y;clustersize" << endl;
	// make list of hits
	for (Int_t L = 0; L < nLayers; L++) { // loop over layers
		Hits *hits = new Hits(Runs*nLayers*100);
		Clusters *clusters = new Clusters(Runs*nLayers*20);

		f.FindHits(Frame3D->at(L), hits, L);
		Int_t nhits = hits->GetEntriesFast();

		cout << "Found " << nhits << " hits." << endl;

		clock_t t;
		t = clock();
		cout << "Finding clusters for layer " << L << "....";
		
		// in ClusterFinder
		f. FindClustersFromHits(hits, clusters, Frame3D);

		t = clock() - t;
		cout << "...done (" << ((float)t)/CLOCKS_PER_SEC << " s)" << endl;

		for (Int_t i = 0; i < clusters->GetEntriesFast(); i++) { // write clusters to file
			file  << L << ";"
					<< clusters->getX(i) << ";"
					<< clusters->getY(i) << ";"
					<< clusters->getSize(i) << ";" << endl;
		} // end loop over clusters
	} // end loop over layers

} // end function FindClustersToLayers

void Draw2DProjection(Int_t Runs, Int_t Events) {
	
	// first get tracks

	Focal f;
	Tracks *allTracks = new Tracks(2*Events*Runs);

	allTracks->SetOwner(true);

	for (Int_t i=0; i<Runs; i++) {
		Int_t EventFrom = Runs * i;
		Int_t EventTo = EventFrom + Events;

		// could use absorbobjects(tracks), but i don't get it to work
		// Track objects are small though.
		
		Tracks *tracks = f.FindTracks(EventFrom, EventTo);
		tracks->SetOwner(true);

		cout << "--- Run " << i+1 << " of " << Runs << " complete --- \n";

		for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
			allTracks->appendTrack(tracks->At(j));
		}
		delete tracks;
	}

	Float_t angleCut = 5.;
	Float_t TLCut = 15.;
	Float_t x0, y0, theta0, TL;
	Int_t nPoints = 0;

	allTracks->extrapolateToLayer0();

	Int_t ntracks = allTracks->GetEntriesFast();

	TCanvas *c1 = new TCanvas("2D Projection of data");
	
	cout << TLCut;
	
	TH2F *Frame2D = new TH2F("Frame2D", "2D Projection of data, at first detector frame",
			nx/8, 0, nx*2, ny/8, 0, ny*2);

//   gPad->Update();

//   TPaveText *pt = (TPaveText*) gPad->FindObject("title");
//   pt->InsertText(Form("Angle threshold at %d degrees", angleCut));
//   pt->InsertText(Form("Track length treshold at %d mm", TLCut));

	for (Int_t i=0; i<ntracks; i++) {
		// loop through tracks to get the relevant information:
		//    - track starting point
		//    - total track angle (until end point)
		//    - track length (somewhat \propto WEPL)
		//
		// Then cut at the angle, and plot a value \propto the TL at the track starting point
		x0 = allTracks->At(i)->getX(0);
		y0 = allTracks->At(i)->getY(0);
		theta0 = allTracks->At(i)->getSlopeAngle();
		TL = allTracks->At(i)->getTrackLengthmm();
		
		if (theta0 > angleCut) continue;
		if (TL < TLCut) {
			cout << "TL = " << TL << ", continuing\n";
			continue;
		}

		Frame2D->Fill(x0, y0, TL);
		nPoints++;
	}

	delete allTracks;

	cout << "Plotted " << nPoints << " data points\n";

	Frame2D->Draw("COLZ");
	gStyle->SetOptStat(0);

}
*/

Tracks * getTracks(Int_t Runs) {

	Focal *f = new Focal();

	Int_t nClusters = kEventsPerRun * 5;
	Int_t nHits = kEventsPerRun * 50;
	Int_t nTracks = kEventsPerRun * 2;

	CalorimeterFrame *cf = new CalorimeterFrame();
	Clusters * clusters = new Clusters(nClusters);
	Hits *hits = new Hits(nHits);
	Tracks *tracks = new Tracks(nTracks);
	Tracks *allTracks = new Tracks(nTracks * Runs);

	for (Int_t i=0; i<Runs; i++) {
		f->getFrame3D(i, cf);
		cf->diffuseFrame();
		hits = cf->findHits();
		clusters = hits->findClustersFromHits();
		tracks = clusters->findTracks();

		for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
			allTracks->appendTrack(tracks->At(j));
		}

		allTracks->appendClustersWithoutTrack(clusters->getClustersWithoutTrack());

		cf->Reset();
		hits->clearHits();
		clusters->clearClusters();
		tracks->clearTracks();
	}
	return allTracks;
}

void DrawTracks3D(Int_t Runs) {
	Tracks *tracks = getTracks(Runs);
	tracks->extrapolateToLayer0();

	TCanvas *c1 = new TCanvas("c1");
	TView *view = TView::CreateView(1);
	Int_t xf = 0;
	Int_t yf = 0;
	Int_t zf = 0;
	Int_t xt = 2*nx;
	Int_t yt = 2*ny;
	Int_t zt = nLayers;

	view->SetRange(xf, zf, yf, xt, zt, yt);

	TClonesArray *restPoints = tracks->getClustersWithoutTrack();

	TPolyMarker3D *p = new TPolyMarker3D(restPoints->GetEntriesFast(), 7);
   p->SetMarkerColor(kRed);

   for (Int_t i=0; i<restPoints->GetEntriesFast(); i++) {

      if (!restPoints->At(i))
      	continue;

      Cluster *thisCluster = (Cluster*) restPoints->At(i);
      Float_t x = thisCluster->getX();

      // reverse to get correct viewpoint
      Float_t z = thisCluster->getY();
      Float_t y = thisCluster->getLayer();

      p->SetPoint(i, x, y, z);
   }

   p->Draw();
   
	Int_t ntracks = tracks->GetEntriesFast();

	// for tracks
	for (Int_t i=0; i<ntracks; i++) {
		Track *thisTrack = tracks->At(i);
		if (thisTrack->getTrackLengthmm() < 2) continue;
		
		Int_t n = thisTrack->GetEntriesFast();

		TPolyLine3D *l = new TPolyLine3D(n);
		l->SetLineWidth(2);
		
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
}

/*

void DrawRealTracks3D(Int_t Runs) {

   // Draw a vector representation of the tracks
   // found in TrackCollections

   Focal f;
   Clusters *restPoints = new Clusters(10000);
   Tracks *tracks = f.FindRealTracks(Runs, restPoints);

	TCanvas *c1 = new TCanvas("c1");
	TView *view = TView::CreateView(1);
	Int_t xf = 0.2*nx*dx;
	Int_t yf = 0.2*ny*dy;
	Int_t zf = 0.1*dz*nLayers/3;
	Int_t xt = 1.8*nx*dx;
	Int_t yt = 1.8*ny*dy;
	Int_t zt = 0.9*dz*nLayers/3;

	view->SetRange(xf, zf, yf, xt, zt, yt);

   // for Clusters

   TPolyMarker3D *p = new TPolyMarker3D(restPoints->GetEntriesFast(), 7);
   p->SetMarkerColor(kRed);

   for (Int_t i=0; i<restPoints->GetEntriesFast(); i++) {
      Cluster *thisCluster = restPoints->At(i);
      Float_t x = thisCluster->getXmm();

      // reverse to get correct viewpoint
      Float_t z = thisCluster->getYmm();
      Float_t y = thisCluster->getLayermm();

      p->SetPoint(i, x, y, z);
   }

   p->Draw();

	Int_t ntracks = tracks->GetEntriesFast();
	// for tracks
	for (Int_t i=0; i<ntracks; i++) {
		Track *thisTrack = tracks->At(i);
		
		Int_t n = thisTrack->GetEntriesFast();

		TPolyLine3D *l = new TPolyLine3D(n);
		
		// for lines
		for (Int_t j=0; j<n; j++) {
			Float_t x = thisTrack->getXmm(j);
			Float_t z = thisTrack->getYmm(j);
			Float_t y = thisTrack->getLayermm(j);
			l->SetPoint(j,x,y,z);
		}
		l->Draw();
	} // end loop over tracks
	view->ShowAxis();
   c1->Update();

} // end function DrawRealTracks3D


void FindClusters(Int_t Runs, Int_t Layer) {
   // loop through hits in histogram and make a TTree branch of them
   // Find all hits with neighbours and return a TTree of clusters with size
   //			+ center of gravity

   Focal f; 

   TCanvas *c1 = new TCanvas("c1", "multipads", 1400, 900);
   gStyle->SetOptStat(0);
   c1->Divide(2,1,0.01,0.01, 0);

   TH2F *Frame2D = new TH2F("Frame2D", Form("Energy deposition in layer %i [keV]", Layer), 
                              nx*2, 0, 2*nx, ny*2, 0, 2*ny);
  
   f.GetFrame2D(Runs, Layer, Frame2D);
   f.DiffuseFrame(Frame2D);

   
   // Until layer geometry is flipped
//   Layer = (23 - Layer);
//
//   Long64_t nbytes = 0, nb = 0;
//   for (Long64_t jentry=0; jentry<nentries;jentry++) {
//      Long64_t ientry = LoadTree(jentry);
//      if (ientry < 0) break;
//      nb = fChain->GetEntry(jentry);   nbytes += nb;
//
//      if (eventID > Runs) break;
//
//      if ((Layer > 0 && level1ID == Layer) || (Layer < 0)) {
//         Int_t H = level2ID; // up or down
//         Int_t C = level3ID; // left or right
//         Int_t px = level4ID; // pixel number
//         Int_t X = px % nx + nx * C;
//         if (H == 1) px = ny*nx - px; // if H=1 (down) flip y values
//         Int_t Y = px / nx + ny * (1 - H);
//
//         DiffuseTH2F(Frame2D, X, Y, edep);
//
//      } // end if layer is correct
//   } // end loop over entries

   

   // make list of hits
   Hits *hits = new Hits(Runs*100);
   Clusters *clusters = new Clusters(Runs*20);

   f.FindHits(Frame2D, hits, Layer);
   Int_t nhits = hits->GetEntriesFast();

   cout << "Found " << nhits << " hits." << endl;

   clock_t t;
   t = clock();
   cout << "Finding clusters....";

   f.FindClustersFromHits(hits, clusters);

   t = clock() - t;
   cout << "...done (" << ((float)t)/CLOCKS_PER_SEC << " s)" << endl;

   cout << "Some statistics...." << endl;
   cout << "From " << Runs << " primary particles, " << clusters->GetEntriesFast() 
         << " clusters were found." << endl;
   cout << "The efficiency is " << float(clusters->GetEntriesFast()) / float(Runs) * 100. 
         << "%." << endl;

   c1->cd(1);
   Frame2D->Draw("COLZ");

   // Draw ellipses on each cluster:

   for (Int_t i = 0; i < clusters->GetEntriesFast(); i++) {
      Float_t clusterRadius = sqrt(clusters->getSize(i) / 3.14159265);
      TEllipse *ellipse = new TEllipse(clusters->getX(i), clusters->getY(i),
                                       clusterRadius, clusterRadius);
      ellipse->SetLineWidth(1);
      ellipse->SetLineColor(2);
      ellipse->SetFillStyle(0);
      ellipse->Draw();
   }

   c1->Update();

   TH1I *cSizeDistribution = new TH1I("cSizeDistribution", 
         "Cluster Size Distribution using connected neighbours", 70, 0, 70);

   for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
      cSizeDistribution->Fill(clusters->getSize(i));
   }

   c1->cd(2);
   cSizeDistribution->Draw();
   c1->Update();
   
   // Write information to CSV file

   ofstream csvfile("pixels.csv");
   csvfile << "Size; x; y" << endl;
   for (Int_t i = 0; i < clusters->GetEntriesFast();  i++) {
      csvfile << clusters->getSize(i) << ";" <<
               clusters->getX(i) << ";" <<
               clusters->getY(i) << endl;
   } // end draw first nout clusters
   csvfile.close();

} // end function FindClusters

void DrawDiffusionCheck(Int_t Runs, Int_t Layer) {
   Focal f;

   TCanvas *c1 = new TCanvas("c1", "multipads", 1400, 900);
   gStyle->SetOptStat(0);
   c1->Divide(2,1,0.01,0.01,0);

   c1->cd(2);
   TH2F *Frame2D = new TH2F("Frame2D", "Diffused hitsmap in all layers", 
                              nx*2, 0, nx*2, ny*2, 0, ny*2);
   Frame2D->SetXTitle("Pixel number");
   Frame2D->SetYTitle("Pixel number");
   
   c1->cd(1);

   f.GetFrame2D(Runs, Layer, Frame2D);
   TH2F *DiffusedFrame2D = (TH2F*) Frame2D->Clone();
   
   f.DiffuseFrame(DiffusedFrame2D);
   
   DiffusedFrame2D->SetName("DiffusedFrame2D");
   DiffusedFrame2D->SetTitle("Original hitsmap in all layers");

   c1->cd(1);
   Frame2D->Draw("COLZ");
   c1->cd(2);
   DiffusedFrame2D->Draw("COLZ");

   c1->Update();
} // end function Frame2DWithDiffusion

void DrawRealFrame2D(Int_t Runs, Int_t Layer) {
   Focal f;

   TH2F *Frame2D = new TH2F("Frame2D", "Hitsmap in all layers", 
                              nx*2, 0, nx*2, ny*2, 0, ny*2);
   Frame2D->SetXTitle("Pixel number");
   Frame2D->SetYTitle("Pixel number");
  	
   f.GetRealFrame2D(Runs, Layer, Frame2D);

	Frame2D->Draw("COLZ");
	gStyle->SetOptStat(0);

} // end function DrawRealFrame2D

void DrawFrame2D(Int_t Runs, Int_t Layer) {
   // Draw one layer using GetFrame2D

   Focal f;

   TH2F *Frame2D = new TH2F("Frame2D", "Energy deposition in all layers [keV]",
                              nx*2, 0, nx*2, ny*2, 0, ny*2);

   f.GetFrame2D(Runs, Layer, Frame2D);

   Frame2D->Draw("COLZ");
   gStyle->SetOptStat(0);
}

void DrawData3D(Int_t RunFrom, Int_t RunTo) {
   // Get hits directly from posX / posY / posZ to show MC data
   // Plots runs from RunFrom until RunTo
   // no histogramming
   // no diffusion
   // NO NOTHING

   Focal f;
   TH3F *Frame3D = new TH3F("Frame3D", "3D map of energy deposition [keV]", 
                              100, -250, 50, 100, -20, 20, 100, -20, 20);

   Frame3D->SetXTitle("Z axis");
   Frame3D->SetYTitle("X axis"); // to get projection right (Z is depth, not up)
   Frame3D->SetZTitle("Y axis"); 

   f.GetData3D(RunFrom, RunTo, Frame3D);

   Frame3D->Draw("LEGO");
}

void DrawRealData3D(Int_t RunFrom, Int_t RunTo) {
   // Get hits directly from posX / posY / posZ to show MC data
   // Plots runs from RunFrom until RunTo
   // no histogramming
   // no diffusion
   // NO NOTHING

   Focal f;
   TH3F *Frame3D = new TH3F("Frame3D", "3D map of energy deposition [keV]", 
                              25, 0, 25, 100, -650, 650, 100, -650, 650);

   Frame3D->SetXTitle("Z axis");
   Frame3D->SetYTitle("X axis"); // to get projection right (Z is depth, not up)
   Frame3D->SetZTitle("Y axis"); 

   f.GetRealData3D(RunFrom, RunTo, Frame3D);

   Frame3D->Draw("LEGO");
}

void MakeLayerPNGs(Int_t Runs, Int_t Diffusion) {
   // write PNG files from all layers
   // diffusion = 1 or = 0 to activate / deactivate
   
   Focal f;
   Int_t nlayers = 23;

   TCanvas *c = new TCanvas;

   vector<TH2F*> *Frame3D = new vector<TH2F*>;
   f.GetFrame3D(0, Runs, Frame3D);

   if (Diffusion > 0) {
      f.DiffuseFrame(Frame3D);
   }

   for (Int_t Layer = 0; Layer < nlayers; Layer++) {
   
      // Save TH2F to Diffused_Frame2D_LayerXX.png
      Frame3D->at(Layer)->Draw("COLZ");
      if (Diffusion > 0) {
         c->Print(Form("Frames/Diffused_Frame2D_Layer%i.png",Layer));
      }
      else {
         c->Print(Form("Frames/Undiffused_Frame2D_Layer%i.png", Layer));
      }
   } // end loop over layers
   delete Frame3D;
   delete c;
} // end function MakeLayerPNGs

// start comment here
 * Must be added in Focal:: namespace before
 * it can work properly, but we don't need it for now...
 * 
void EdepHistogram()
{
   if (fChain == 0) return;

Long64_t nentries = fChain->GetEntriesFast();
   TH1I *eDepHisto = new TH1I("eDepHisto", "Energy Deposition on pixels [keV]", 
                              1000, 0, 3000);
   eDepHisto->SetXTitle("Deposited Energy [keV]");
   eDepHisto->SetYTitle("Number of pixels");

   Long64_t nb = 0;
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry<0) break;
      nb = fChain->GetEntry(jentry);
      eDepHisto->Fill(edep * 1000);
   }
   eDepHisto->Draw();
}
// end comment here

 */
