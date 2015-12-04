#define Analysis_cxx
#include "Analysis.h"
#include "Constants.h"
#include "MaterialConstants.h"
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
#include <TCanvas.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <TEllipse.h>
#include <vector>
#include <TStopwatch.h>
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

void getTrackStatistics(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {
	Focal f;

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
	
	Int_t nTracksToPlot = 25;
	Int_t nTracksToPlot1D = 5;

	TString sDataType = getDataTypeString(dataType);

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

	TH1F *hTrackLengths = new TH1F("hTrackLengths", Form("Track Lengths (%s)", sDataType.Data()), 100, 0, 120);
	TH2F *hClusterSizeAlongTrack = new TH2F("hClusterSizeAlongTrack",
				Form("Cluster size along track length (%s)", sDataType.Data()), 1.5*nLayers, 0, 1.5*nLayers*dz, 50, 0, 50);
	TH1F *hStraightness = new TH1F("hStraightness", Form("Sinuosity plot (%s)", sDataType.Data()), 500, 1, 1.01);
	TH1F *hSlope = new TH1F("hSlope", Form("Proton angle plot (%s)", sDataType.Data()), 500, 0, 20);

	// Average cluster size
	vector<TH1F*> *hAvgCS = new vector<TH1F*>;
	hAvgCS->reserve(4);
	for (Int_t chip=0; chip<4; chip++) {
		hAvgCS->push_back(new TH1F(Form("hAvgCS_chip_%i",chip),
				Form("Average Cluster Size vs Track Length for chip %i (%s)",chip, sDataType.Data()), 50, 0, 50));
	}

	// Cluster size for individual layers
	vector<TH1F*> *hCSLayer = new vector<TH1F*>;
	hCSLayer->reserve(9);
	for (Int_t layer=0; layer<9; layer++) {
		hCSLayer->push_back(new TH1F(Form("hCSLayer_%i", layer),
				Form("Cluster size for layer %i (%s)", layer, sDataType.Data()), 50, 0, 50));
	}

	// Cluster size along track length for a single track
	vector<TH1F*> *hFollowTrack = new vector<TH1F*>;
	hFollowTrack->reserve(nTracksToPlot);
	for (Int_t track=0; track<nTracksToPlot; track++) {
		hFollowTrack->push_back(new TH1F(Form("hFollowTrack_%i", track),
				Form("Cluster size along track length for a single track (%s)", sDataType.Data()), 50, 0, 50));
		hFollowTrack->at(track)->SetXTitle("Track Length [mm]");
		hFollowTrack->at(track)->SetYTitle("Cluster size [# of pixels]");
	}

	// Proton angle distribution in layer
	vector<TH1F*> *hAngles = new vector<TH1F*>;
	hAngles->reserve(9);
	for (Int_t layer=0; layer<9; layer++) {
		hAngles->push_back(new TH1F(Form("hAngles_%i", layer),
				Form("Proton angle distribution in layer %i (%s)", layer, sDataType.Data()), 50, 0, 20));
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
				hFollowTrack->at(trackNum)->SetTitle(Form("Track length histogram for run %i (%s)", i, sDataType.Data()));
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

void drawClusterShapes(Int_t Runs, Bool_t dataType, Bool_t recreate, Int_t energy) {
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
			cf->diffuseFrame(); // THE MAGIC PART
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

TString getDataTypeString(Int_t dataType) {
	TString sDataType;
	if (dataType == kMC)
		sDataType = "MC";
	else if (dataType == kData)
		sDataType = "Exp. data";
	else
		sDataType = "Unknown source";

	return sDataType;
}

Int_t getMinimumTrackLength(Int_t energy) {
	Int_t minTL = 0;

	if (energy < 150) minTL = 5;
	else if (energy < 170) minTL = 10;
	else if (energy < 190) minTL = 15;
	else if (energy < 200) minTL = 20;
	else if (energy < 230) minTL = 23;

	return minTL;
}

Double_t fitfunc_DBP(Double_t *v, Double_t *par) {
	// we convert to water equivalent path length before fitting
	// so the numbers are given in water
	// Based on Bortfeld and Schlegel 1996
	// v[0] is depth and par[0] is initial energy E0. par[1] is a scaling.

	Float_t R_mm = 0.0019 * pow(par[0],1.8) * 10;
	Double_t fitval = par[1] / ( 0.0554 * pow((R_mm - v[0]), 0.444) );
	if (isnan(fitval)) fitval = 0;
	return fitval;
}

void makeTracks(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
}

void drawBraggPeakGraphFit(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {

	run_energy = energy;
	TString sDataType = getDataTypeString(dataType);
	
	Float_t trackLengthSoFar;
	Bool_t cut;

	TCanvas *cGraphAll = new TCanvas("cGraphAll", "Fitting data points", 1000, 800);
	TCanvas *cGraph = new TCanvas("cGraph", "Fitted data points", 1400, 1000);
	TCanvas *c3 = new TCanvas("c3", "Fit results");
	cGraph->Divide(4,4, 0.000001, 0.000001, 0);

	char *sMaterial;
	if (kMaterial == kTungsten) { sMaterial = "W"; }
	else if (kMaterial == kAluminum) { sMaterial = "Al"; }
	else if (kMaterial == kPMMA) { sMaterial = "PMMA";}
	
	TH1F *fitResult = new TH1F("fitResult", Form("Fitted energy of a %d MeV beam in %s", energy, sMaterial), 500, 100, 325);

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();

	vector<TGraph*> vGraph;
	vGraph.reserve(tracks->GetEntriesFast());

	Int_t m = 0;
	Track *thisTrack = 0;
	Float_t xx[nLayers*tracks->GetEntriesFast()];
	Float_t yy[nLayers*tracks->GetEntriesFast()];

	Int_t cutWEPL = 0, cutBPinT = 0;
	
	for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
		thisTrack = tracks->At(j);
		if (!thisTrack) continue;
		
		if (getCutWEPL(thisTrack)) cutWEPL++;
		if (getCutBraggPeakInTrack(thisTrack)) cutBPinT++;
		
		cut = getCutWEPL(thisTrack);//* getCutBraggPeakInTrack(thisTrack);
		//if (dataType == kData) cut *= getCutChipNumber(thisTrack);

		if (!cut) continue;
		Int_t n = thisTrack->GetEntriesFast();
		Float_t x[n], y[n];

		Float_t trackLengthSoFar = 0;
		Float_t trackLengthSoFarWEPL = 0;
		for (Int_t k=0; k<n; k++) {
			if (!thisTrack->At(k)) continue;
			trackLengthSoFar += thisTrack->getTrackLengthmmAt(k);
			trackLengthSoFarWEPL += thisTrack->getTrackLengthWEPLmmAt(k);
			x[k] = trackLengthSoFarWEPL;
			y[k] = thisTrack->getDepositedEnergy(k);
			xx[m+k] = trackLengthSoFarWEPL;
			yy[m+k] = thisTrack->getDepositedEnergy(k);
		}
		m += n;
		vGraph.push_back(new TGraph(n, x, y));
	}
	
	cout << "Of " << tracks->GetEntriesFast() << ", " << cutWEPL << " made the WEPL cut and " << cutBPinT << " made the BP cut.\n";

	vGraph.push_back(new TGraph(m, xx, yy));
	
	Int_t gsize = vGraph.size();
	Int_t plotSize = 16;
	//if (gsize>25) gsize = 25;
	
	for (Int_t i=0; i<gsize; i++) {
		if (!vGraph.at(i)) continue;
		TGraph *gr = vGraph.at(i);
		if (i < plotSize) {
			cGraph->cd(i+1);
		   gStyle->SetPadBorderMode(0);
		   gStyle->SetFrameBorderMode(0);
		   gStyle->SetTitleH(0.06);
		   gStyle->SetTitleYOffset(1);
			gr->SetMaximum(500);
			gr->SetMinimum(0);
			gr->SetTitle("");
			gr->GetXaxis()->SetTitle("Water Equivalent Path Length [mm]");
			gr->GetXaxis()->SetTitleSize(0.05);
			gr->GetYaxis()->SetTitleSize(0.05);
			gr->GetYaxis()->SetTitle("Deposited energy [keV/#mum]");
			gr->GetXaxis()->SetLabelSize(0.04);
			gr->GetYaxis()->SetLabelSize(0.04);
			gr->Draw("A*");
			if (i<gsize-1)
				gr->GetXaxis()->SetRangeUser(0, 500);
			gr->Draw("A*");
			cGraph->Update();
		}

		// fitting
		Int_t n = gr->GetN();
		Float_t maxVal = 0;
		Int_t maxIdx = 0;
		
		for (Int_t j=0.5*n; j<n; j++) {
		    if (gr->GetY()[j] > maxVal) {
		    maxVal = gr->GetY()[j];
		    maxIdx = j;
		   }
		}
				
		Float_t constHeight = 7;
		Float_t BPpos = maxIdx * dz;
		Float_t BP = maxVal - constHeight;
		Float_t BPwidth = 1.7;

		/*
		TF1 *g1 = new TF1("m1", "gaus(0) + [3]", 0, BPpos*2);
		g1->SetParameters(BP, BPpos, BPwidth, constHeight); // BP, BPpos, 2.5
		g1->SetParLimits(0, BP*0.8, BP*2);
		g1->SetParLimits(1, BPpos*0.9, BPpos*1.3);
		g1->SetParLimits(2, BPwidth, BPwidth);
		gr->Fit("m1", "B, W, Q", "", 0, BPpos*1.5);
		Float_t fit_tl = g1->GetParameter(1);
		//cout << Form("The fitted energy is %d MeV.\n", thisTrack->getEnergyFromWEPL(fit_tl));
		 */

		TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
		func->SetParName(0,"Initial energy [MeV]");
		func->SetParName(1, "Factor");
		func->SetParameter(0,energy);
		func->SetParameter(1, 40);
		func->SetParLimits(0, 10, energy*1.25);
		func->SetParLimits(1, 30,50);
		gr->Fit("fit_BP", "B, W, Q", "", 0, 500);
		Float_t fit_t = func->GetParameter(0);
		if (fit_t < energy*1.24) fitResult->Fill(fit_t);
// 		cout << Form("The fitted energy is %f MeV... Jesus christ that wasn't hard was it\n", fit_t);

		if (i<plotSize) {
			TLatex *myl = new TLatex(20,400,Form("Fitted energy: %.1f MeV", fit_t));
			myl->SetTextSize(0.06);
			myl->Draw();
		}
			

	}

	cGraphAll->cd();
	TGraph *g = vGraph.at(vGraph.size() - 1);
	g->GetXaxis()->SetTitle("Water Equivalent Path Length");
	g->Draw("A*");
	
	//Track *t = new Track;
	// fit main function
	TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 300, 2);
	func->SetParName(0,"Initial energy [MeV]");
	func->SetParName(1, "Factor");
	func->SetParameter(0,200.);
	func->SetParameter(1, 40);
	func->SetParLimits(0, 10, 215);
	func->SetParLimits(1, 30,50);
	g->Fit("fit_BP", "B, W, Q", "", 0, 300);
	Float_t fit_t = func->GetParameter(0);
	//cout << Form("The fitted energy is %f MeV... Jesus christ that wasn't hard was it\n", fit_t);

	c3->cd();
	fitResult->SetXTitle("Energy [MeV]");
	fitResult->SetYTitle("Number of protons");
	fitResult->Draw();

  	delete tracks;
}

void drawBraggPeakFit(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {

	Focal f;
	Float_t trackLengthSoFar = 0;
	Int_t trackNum = 0, chip = 0, minTL = 0;
	Bool_t cutTrackLength, cutChipNumber, cutBraggPeakInTrack;
	Bool_t cutCombined, BPSampled;
	Int_t okCHIP = 0, okTL = 0, okBP = 0, okALL = 0, okNTBP = 0;
	vector<TArrow*> arrows;
	arrows.reserve(100);

	vector<TArrow*> cogArrows;
	cogArrows.reserve(100);

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, kCalorimeter, energy);
	tracks->extrapolateToLayer0();
	
	TString sDataType = getDataTypeString(dataType);

	TCanvas *cAverageTrack = new TCanvas("cAverageTrack", "Average Cluster Size along track", 1000, 800);
	TH1F *hAverageTrack = new TH1F("hAverageTrack", "Average Cluster Size along track", nLayers+1, -dz/2, nLayers*dz+dz/2);

	TCanvas *cBPResult = new TCanvas("cBPResult", "Bragg Peak fit mean values", 1000, 800);
	TH1F *hBPResult = new TH1F("hBPResult", "Bragg peak fit mean values", 400, 0, 30);

	TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
	TH1F *hMu = new TH1F("hMu", "Weighted CoG BP position", 400, 0, 30);

	TCanvas *cFollowTrack = new TCanvas("c1", "c1", 1000, 800);
	cFollowTrack->Divide(4, 4, 0.01, 0.01, 0);

	vector<TH1F*> *hFollowTrack = new vector<TH1F*>;
	hFollowTrack->reserve(tracks->GetEntriesFast());
	for (Int_t track=0; track<tracks->GetEntriesFast(); track++) {
		hFollowTrack->push_back(new TH1F(Form("hFollowTrack_%i", track), "Cluster size along track length for a single track", nLayers+1, -dz/2, nLayers*dz+dz/2));
		hFollowTrack->at(track)->SetXTitle("Track Length [mm]");
		hFollowTrack->at(track)->SetYTitle("Cluster size [# of pixels]");
	}

	Track *thisTrack;
	for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
		thisTrack = tracks->At(i);

		cutTrackLength = getCutTrackLength(energy, thisTrack);
		cutBraggPeakInTrack = getCutBraggPeakInTrack(thisTrack);

		if (dataType == kData)
			cutChipNumber = getCutChipNumber(thisTrack);
		else
			cutChipNumber = true;

		
		cutCombined = cutBraggPeakInTrack * cutTrackLength;
		if (cutCombined) okALL++;

		trackLengthSoFar = 0;
		for (Int_t j=0; j<thisTrack->GetEntriesFast(); j++) {
			trackLengthSoFar += thisTrack->getTrackLengthmmAt(j);
			Float_t trackLengthLayer = thisTrack->getLayermm(j);

//			Float_t arrowPos = trackLengthSoFar + dz/2;

			if (cutCombined) {
				hAverageTrack->Fill(trackLengthLayer, thisTrack->getSize(j));
				hFollowTrack->at(trackNum)->Fill(trackLengthLayer, thisTrack->getSize(j));
			}

			/*
			if (trackNum < 16 && cutCombined) {
				if (BPSampled && BPidx == j) {
					arrows.push_back(new TArrow(arrowPos, thisTrack->getSize(j) * 1.3,
										arrowPos, thisTrack->getSize(j) * 1.1, 0.005, "|>"));
					arrows.back()->SetLineColor(kRed);
					arrows.back()->SetLineWidth(2.);
					arrows.back()->SetFillColor(kRed);
				}

				else if (!BPSampled && BPidx == j)
					arrows.push_back(NULL);
			}
			*/
		}

		if (cutCombined) {
			Float_t BPpos = hFollowTrack->at(trackNum)->GetMaximumBin() * dz;
			Float_t BP = hFollowTrack->at(trackNum)->GetBinContent(BPpos/dz);
			Float_t BPwidth = 1.7;

			TF1 *g1 = new TF1("m1", "gaus", 0, BPpos*2);

			g1->SetParameters(BP, BPpos, BPwidth); // BP, BPpos, 2.5
			g1->SetParLimits(0, BP*0.8, BP*1.2);
			g1->SetParLimits(1, BPpos*0.85, BPpos*1.15);
			g1->SetParLimits(2, BPwidth, BPwidth);

			if (trackNum<16)
				cFollowTrack->cd(trackNum+1);

			hFollowTrack->at(trackNum)->Fit("m1", "B, W, Q", "", BPpos*0.7, BPpos*1.3);

			cout << "BPpos = " << BPpos << ", mean value after fit: " << g1->GetParameter(1) << endl;
			hBPResult->Fill(g1->GetParameter(1));
		}

		if (cutCombined) trackNum++;
		if (cutTrackLength) okTL++;
		if (cutChipNumber) okCHIP++;
		if (cutBraggPeakInTrack) okBP++;
	}

	hAverageTrack->Scale(1/(float) trackNum);

	TF1 *g1 = new TF1("m1", "gaus", 15, 35);

	Float_t BPpos = hAverageTrack->GetMaximumBin() * dz;
	Float_t BP = hAverageTrack->GetBinContent(BPpos/dz);
	Float_t BPwidth = 1.7;

	g1->SetParameters(BP, BPpos, BPwidth); // BP, BPpos, 2.5
	g1->SetParLimits(0, BP*0.8, BP*1.2);
	g1->SetParLimits(1, BPpos*0.85, BPpos*1.15);
	g1->SetParLimits(2, BPwidth, BPwidth);

	cAverageTrack->cd();

	hAverageTrack->Fit("m1", "R, WW", "", BPpos*0.7, BPpos*1.3);

	cout << "Total tracks: " << tracks->GetEntriesFast() << endl;
	cout << "Total okCHIP: " << okCHIP << endl;
	cout << "Total okTL: " << okTL << endl;
	cout << "Total okBP: " << okBP << endl;
	cout << "Total ok combined: " << okALL << endl;

	cAverageTrack->cd();

	hAverageTrack->SetFillColor(kOrange+7);
	hAverageTrack->Draw();

	//TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);

	Float_t cog;

	for (UInt_t track=0; track<hFollowTrack->size(); track++) {
		cog = 0;
		Int_t hLen = hFollowTrack->at(track)->FindLastBinAbove();
		if (hLen<3) continue;

		// pen-penultimate, penultimate and ultimate bin
		Float_t ppu = hFollowTrack->at(track)->GetBinContent(hLen-2);
		Float_t pu = hFollowTrack->at(track)->GetBinContent(hLen-1);
		Float_t u = hFollowTrack->at(track)->GetBinContent(hLen);

		// sometimes penultimate channel is bad, in that case fix it
		if (pu<1 && ppu > pu*5) pu = ppu;

		// cog for last two bins
		cog = (pu * (hLen-1)*dz + u * (hLen)*dz) / (u + pu);

		hMu->Fill(cog);
		
		/*
		cogArrows.push_back(new TArrow(cog, max(u,pu) * 1.3,
							cog, max(u,pu) * 1.1, 0.005, "|>"));
		cogArrows.back()->SetLineColor(kGreen);
		cogArrows.back()->SetLineWidth(2.);
		cogArrows.back()->SetFillColor(kGreen);
		*/
	}

	for (UInt_t trackNo=0; trackNo<16; trackNo++) {
		cFollowTrack->cd(trackNo+1);
		gPad->DrawFrame(0, 0, 50, 45);
		hFollowTrack->at(trackNo)->SetFillColor(kBlue-2);
		hFollowTrack->at(trackNo)->Draw("same");

		/*
		if (trackNo < arrows.size()) {
			if (arrows.at(trackNo)) {
				arrows.at(trackNo)->Draw();
			}
		}
		if (trackNo < cogArrows.size()) {
			if (cogArrows.at(trackNo)) {
				cogArrows.at(trackNo)->Draw();
			}
		}
		*/
	}

	cBPResult->cd();
	hBPResult->Draw();

	c2->cd();
	hMu->SetFillColor(kGray);
	hMu->SetTitle(Form("Bragg peak positions for %i MeV protons", energy));
	hMu->Draw();

	delete tracks;
}

/*
	
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

void saveTracks(Tracks *tracks, Int_t dataType, Int_t energy) {
  	
	// C++ / ROOT has something to learn from Python... ;)
	TString sDataType = (dataType == 0) ? "_MC" : "_data";
	TString sEnergy = Form("_%dMeV", energy);
	TString fileName = "tracks";
	fileName.Append(sDataType);
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

Tracks * loadTracks(Int_t Runs, Int_t dataType, Int_t energy) {
  	TString sDataType = (dataType == 0) ? "_MC" : "_data";
	TString sEnergy = Form("_%dMeV", energy);
	TString fileName = "tracks";
	fileName.Append(sDataType);
	fileName.Append(sEnergy);
	fileName.Append(".root");
	
	TFile *f = new TFile(fileName);
	if (!f) return 0;
	TTree *T = (TTree*) f->Get("T");
	Tracks * tracks = new Tracks();

	T->GetBranch("tracks")->SetAutoDelete(kFALSE);
	T->SetBranchAddress("tracks",&tracks);

	cout << "There are " << T->GetEntriesFast() << " entries in TTree.\n";

	T->GetEntry(0);
	
	cout << "There are " << tracks->GetEntriesFast() << " entries in tracks.\n";
	
	return tracks;
}

Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Int_t energy) {
	Tracks * tracks;
	
	if (recreate) {
		tracks = getTracks(Runs, dataType, kCalorimeter, energy);
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

Tracks * getTracks(Int_t Runs, Int_t dataType, Int_t frameType, Int_t energy) {
	Focal *f = new Focal();

	Int_t nClusters = kEventsPerRun * 5 * nLayers;
	Int_t nHits = kEventsPerRun * 50;
	Int_t nTracks = kEventsPerRun * 2;

	Bool_t breakSignal = kFALSE;

	CalorimeterFrame *cf = new CalorimeterFrame();
	TrackerFrame *tf = new TrackerFrame();
	Clusters * clusters = new Clusters(nClusters);
	Clusters * trackerClusters = new Clusters(nClusters);
	Hits *hits = new Hits(nHits);
	Hits *trackerHits = new Hits(nHits);
	Tracks *calorimeterTracks = new Tracks(nTracks);
	Tracks *trackerTracks = new Tracks(nTracks);
	Tracks *allTracks = new Tracks(nTracks * Runs);
	
	TStopwatch t1, t2, t3, t4, t5, t6;
	
	for (Int_t i=0; i<Runs; i++) {

		cout << "Finding track " << (i+1)*kEventsPerRun << " of " << Runs*kEventsPerRun << "...\n";
		
		if (dataType == kMC) {
			t1.Start();
			f->getMCFrame(i, cf);
			t1.Stop(); t2.Start();
			cf->diffuseFrame();
			t2.Stop(); t3.Start();
			hits = cf->findHits();
			t3.Stop(); t4.Start();
			clusters = hits->findClustersFromHits(); // badly optimized
			t4.Stop();

/*			f->getMCTrackerFrame(i, tf);
			tf->diffuseFrame();
			trackerHits = tf->findHits();
			trackerClusters = trackerHits->findClustersFromHits();*/
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
//		trackerTracks = trackerClusters->findTrackerTracks();

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

void drawTracks3D(Int_t Runs, Int_t dataType, Int_t frameType, Int_t energy) {
	// FIXME Add kTracker to loadOrCreateTracks
  
	Tracks *tracks = getTracks(Runs, dataType, frameType, energy);
	tracks->extrapolateToLayer0();

	TCanvas *c1 = new TCanvas("c1");
	TView *view = TView::CreateView(1);
	view->SetRange(0, 0, 0, 2*nx, nLayers, 2*ny);

	TClonesArray *restPoints = tracks->getClustersWithoutTrack();

	TPolyMarker3D *p = new TPolyMarker3D(restPoints->GetEntriesFast(), 7);
   p->SetMarkerColor(kRed);

   for (Int_t i=0; i<restPoints->GetEntriesFast(); i++) {

      if (!restPoints->At(i))
      	continue;

      Cluster *thisCluster = (Cluster*) restPoints->At(i);
      Float_t x = thisCluster->getX();

      Float_t z = thisCluster->getY();
      Float_t y = thisCluster->getLayer();

      p->SetPoint(i, x, y, z);
   }

   p->Draw();
   
	Int_t ntracks = tracks->GetEntriesFast();

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


void DrawFrame2D(Int_t Runs, Int_t Layer) {
   // Draw one layer using GetFrame2D

   Focal f;

   TH2F *Frame2D = new TH2F("Frame2D", "Energy deposition in all layers [keV]",
                              nx*2, 0, nx*2, ny*2, 0, ny*2);

   f.GetFrame2D(Runs, Layer, Frame2D);

   Frame2D->Draw("COLZ");
   gStyle->SetOptStat(0);
}

*/

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

/*

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

Bool_t getCutTrackLength(int energy, Track *track) {
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



