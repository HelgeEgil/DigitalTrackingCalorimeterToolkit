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
#include <TAxis3D.h>
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


void makeTracks(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
}

void drawBraggPeakFit(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {
	run_energy = energy;
	TStopwatch t1, t2;

	char *title = Form("Fitted energy of a %d MeV beam in %s (%s)", energy, getMaterialChar(), getDataTypeChar(dataType));
	TCanvas *c = new TCanvas("c", title, 2000, 1400);
	TH1F *fitResult = new TH1F("fitResult", title, 100, 80, 240);

	fitResult->SetLineColor(kBlack);
	if (dataType == kMC) { fitResult->SetFillColor(kGreen-5); }
	else if (dataType == kData) { fitResult->SetFillColor(kOrange+4); }

	cout << "Loading tracks...\n";
	t1.Start(); Tracks *tracks = loadOrCreateTracks(recreate, Runs, dataType, energy); t1.Stop();
	cout << Form("Done (%.2f s).\n", t1.RealTime());
	tracks->extrapolateToLayer0();
	cout << "Fitting tracks...\n";
	t2.Start(); tracks->doFit(); t2.Stop();
	cout << Form("Done (%.2f s).\n", t2.RealTime());

	for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
		Track *thisTrack = tracks->At(i);
		if (!thisTrack) {	continue; }

		Float_t res = thisTrack->getFitParameterEnergy();
		if (res) { fitResult->Fill(res); }
	}

	fitResult->Draw();
	
	c->SaveAs(Form("Fitted_energies_%d_MeV_%s_material_%s_datatype.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
}

void drawBraggPeakGraphFit(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {
	run_energy = energy;
	
	char * sDataType = getDataTypeChar(dataType);
	char * sMaterial = getMaterialChar();
	char * hTitle = Form("Fitted energy of a %d MeV beam in %s (%s)", energy, sMaterial, sDataType);

	Int_t nPlotX = 4, nPlotY = 4;
	Int_t fitIdx = 0, plotSize = nPlotX*nPlotY;
	Bool_t isFitOk = false;

	TCanvas *cGraph = new TCanvas("cGraph", "Fitted data points", 1400, 1000);
	TCanvas *cFitResults = new TCanvas("cFitResults", hTitle);
	cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
	TH1F *hFitResults = new TH1F("fitResult", hTitle, 1000, 100, 325);

	gStyle->SetPadBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetTitleH(0.06);
	gStyle->SetTitleYOffset(1);
	
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();

	for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
		Track *thisTrack = tracks->At(j);
		if (!thisTrack) continue;
	
		isFitOk = thisTrack->doFit();
		if (!isFitOk) continue;
			
		Float_t fitEnergy = thisTrack->getFitParameterEnergy();
		Float_t fitScale = thisTrack->getFitParameterScale();

		hFitResults->Fill(fitEnergy);

		Int_t n = thisTrack->GetEntriesFast();
		Float_t x[n], y[n];

		if (fitIdx < plotSize) {		
			cGraph->cd(fitIdx+1);
			
			Float_t trackLengthWEPL = 0;
			for (Int_t k=0; k<n; k++) {
				if (!thisTrack->At(k)) continue;
				trackLengthWEPL += thisTrack->getTrackLengthWEPLmmAt(k);
				x[k] = trackLengthWEPL;
				y[k] = thisTrack->getDepositedEnergy(k);
			}
			
			TGraph *graph = new TGraph(n,x,y);

			graph->SetMinimum(0);
			graph->SetMaximum(600);
			graph->SetTitle("");
			graph->GetXaxis()->SetTitle("Water Equivalent Path Length [mm]");
			graph->GetYaxis()->SetTitle("Deposited energy per layer [keV]");
			graph->GetXaxis()->SetTitleSize(0.05);
			graph->GetYaxis()->SetTitleSize(0.05);
			graph->GetXaxis()->SetLabelSize(0.04);
			graph->GetYaxis()->SetLabelSize(0.04);
			graph->GetXaxis()->SetLimits(0, 300);
			
			graph->Draw("A*");
			

			TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
			func->SetParameters(fitEnergy, fitScale);
			func->Draw("same");

			TLatex *text = new TLatex(20, 400, Form("Fitted energy: %.1f MeV", fitEnergy));
			text->SetTextSize(0.06);
			text->Draw();
		
			cGraph->Update();
			fitIdx++;
		}
	}
	
	cFitResults->cd();
	hFitResults->SetXTitle("Energy [MeV]");
	hFitResults->SetYTitle("Number of protons");
	hFitResults->Draw();

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

*/

void Draw2DProjection(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, kCalorimeter, energy);
	tracks->extrapolateToLayer0();
	
	Int_t hSizeX = nx/8;
	Int_t hSizeY = ny/8;

	Float_t angleCut = 5.;
	Float_t x0, y0, theta0;
	Int_t nPoints = 0;

	Float_t fit_energy;

	Int_t ntracks = tracks->GetEntriesFast();

	char *title = (char*) Form("2D Projection of data from %d MeV proton beam on %s type FOCAL", energy, getMaterialChar());

	TCanvas *c1 = new TCanvas(title);
	
	TH2F *Frame2D = new TH2F("Frame2D", title, hSizeX, 0, nx*2, hSizeY, 0, ny*2);
	TH2F *normalizeFrame = new TH2F("normalizeFrame", "title", hSizeX, 0, nx*2, hSizeY, 0, ny*2);

//   gPad->Update();

//   TPaveText *pt = (TPaveText*) gPad->FindObject("title");
//   pt->InsertText(Form("Angle threshold at %d degrees", angleCut));
//   pt->InsertText(Form("Track length treshold at %d mm", TLCut));

	for (Int_t i=0; i<ntracks; i++) {
		Track *thisTrack = tracks->At(i);

		x0 = thisTrack->getX(0);
		y0 = thisTrack->getY(0);
		theta0 = thisTrack->getSlopeAngle();
		fit_energy = thisTrack->getEnergy();

		if (theta0 > angleCut) continue;

		Frame2D->Fill(x0, y0, energy);
		normalizeFrame->Fill(x0, y0);

		nPoints++;
	}

	Frame2D->Divide(normalizeFrame);

	delete tracks;

	Frame2D->Draw("COLZ");
	gStyle->SetOptStat(0);

}

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

//	cout << "There are " << T->GetEntriesFast() << " entries in TTree.\n";

	T->GetEntry(0);
	
	cout << "There are " << tracks->GetEntriesFast() << " tracks in " << fileName << ".\n";
	
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
	
		cout << "Finished loading tracks\n";

		if (!tracks) {
			cout << "!tracks, creating new file\n";
			tracks = getTracks(Runs, dataType, kCalorimeter, energy);
			saveTracks(tracks, dataType, energy);
		}
	}
	cout << "returning tracks\n";
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
			cf->diffuseFrame(new TRandom3(0));
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

void drawTracks3D(Int_t Runs, Int_t dataType, Bool_t recreate, Int_t energy) {
	// FIXME Add kTracker to loadOrCreateTracks
  
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, kCalorimeter, energy);
	tracks->extrapolateToLayer0();

	TCanvas *c1 = new TCanvas("c1");
	c1->SetTitle(Form("Tracks from %d MeV protons on %s", energy, getMaterialChar()));
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



