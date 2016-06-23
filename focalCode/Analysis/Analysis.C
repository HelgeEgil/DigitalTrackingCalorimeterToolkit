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
#include "Classes/Track/conversionFunctions.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/DataInterface/DataInterface.h"
#include "HelperFunctions/Tools.h"

using namespace std;

void drawTrackAngleAtVaryingRunNumbers(Int_t dataType, Float_t energy) {
	Int_t nRuns = 0;
	Hits * eventIDs = nullptr;

	for (Int_t i=2; i<32; i++) {
		nRuns = pow(2, 4 + 0.3 * i) + 0.5;

		kEventsPerRun = nRuns;
		Float_t factor = 2;

		Int_t totalNumberOfRuns = 5000 / kEventsPerRun;
		if (totalNumberOfRuns < 1) totalNumberOfRuns = 1;
		if (totalNumberOfRuns > 250) totalNumberOfRuns = 250;

		Tracks * tracks = loadOrCreateTracks(1, totalNumberOfRuns, dataType, energy);
		tracks->extrapolateToLayer0();

//		eventIDs = getEventIDs(totalNumberOfRuns, energy);
//		tracks->matchWithEventIDs(eventIDs);

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

		hCorrectTracks->SetTitle(Form("Tracks with same eventID using search cone size of MCS * %.1f and %d protons/run", factor, kEventsPerRun));
		hCorrectTracks->SetXTitle("Layer number");
		hCorrectTracks->SetYTitle("Number of protons");
		hCorrectTracks->SetFillColor(kBlue-7);
		hCorrectTracks->SetLineColor(kBlack);

		Track *thisTrack;
		Int_t EID, thisEID;
		for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
			thisTrack = tracks->At(j);
			if (!thisTrack) continue;
			hAngles->Fill(thisTrack->getSlopeAngleChangeBetweenLayers(2));

			EID = thisTrack->getEventID(0);
			if (EID > -1) { hCorrectTracks->Fill(0); }
			normCorrectTracks->Fill(0);
			for (Int_t k=1; k<thisTrack->GetEntriesFast(); k++) {
				if (!thisTrack->At(k)) continue;
				normCorrectTracks->Fill(thisTrack->getLayer(k));
				thisEID = thisTrack->getEventID(k);
				if (thisEID == EID) {
					hCorrectTracks->Fill(thisTrack->getLayer(k));
				}
			}
		}

		hCorrectTracks->Divide(normCorrectTracks);

		ofstream file2("OutputFiles/lastLayerCorrect_different_nRuns.csv", ofstream::out | ofstream::app);
		file2 << factor << ";" << nRuns << ";" << hCorrectTracks->GetBinContent(5) << endl;
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
	
	if (epr>0) {
		kEventsPerRun = epr;
	}

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
//				hAngles->at(layer)->Fill(thisTrack->getSlopeAngleAtLayer(j));
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

//		cout << "Average cluster size and RMS for layer " << layer << " is \033[1m " << hCSLayer->at(layer)->GetMean() << " pixels \033[0m and \033[1m " << hCSLayer->at(layer)->GetRMS() << "pixels \033[0m\n";
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

//		cout << "The added MCS factor due to inclined crossing is " << 1/(cos(meanAngleAtLayer)) << ".\n";

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
	Float_t fitEnergy, fitScale, fitError = 0;
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
		fitEnergy = quadratureAdd(thisTrack->getFitParameterError(), dz/sqrt(12));
		
		cout << "realEnergy = " << realEnergy << ", fitEnergy = " << fitEnergy << endl;

		hAccuracy->Fill(realEnergy, fitEnergy);

		if (fitIdx < plotSize) {
			drawIndividualGraphs(cGraph, outputGraph, fitEnergy, fitScale, fitError, fitIdx++, eventID, x_energy, y_energy);
		}

		else delete outputGraph;
	}
	cAccuracy->cd();
	hAccuracy->Draw("COLZ");
}

Float_t drawBraggPeakGraphFit(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	run_energy = energy;
	
 	cout << "At energy " << run_energy << ", expecting TL = " << getTLFromEnergy(run_energy) << " and WEPL = " << getWEPLFromEnergy(run_energy) << endl;

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
//	TCanvas *cMaxAngle = new TCanvas("cMaxAngle", "Maxium angle for proton track", 1400, 1000);
	cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
	gStyle->SetPadBorderMode(0); gStyle->SetFrameBorderMode(0);
	gStyle->SetTitleH(0.06); gStyle->SetTitleYOffset(1);
	
	TH1F *hFitResults = new TH1F("fitResult", hTitle, 250, getUnitFromEnergy(0), getUnitFromEnergy(energy*1.2));
	hFitResults->SetLineColor(kBlack); hFitResults->SetFillColor(kGreen-5);

	TH1F *hMaxAngle = new TH1F("hMaxAngle", "Maximum angle for proton track", 200, 0, 25);
	hMaxAngle->SetLineColor(kBlack); hMaxAngle->SetFillColor(kGreen-5);

	Bool_t acceptAngle = false;

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
	
		Float_t maxAngle = 0, thisAngle = 0;
		for (Int_t i=0; i<thisTrack->GetEntriesFast(); i++) {
			thisAngle = thisTrack->getSlopeAngleAtLayer(i);
			maxAngle = max(thisAngle, maxAngle);
		}

		acceptAngle = (maxAngle < 13);

		hMaxAngle->Fill(maxAngle);

		if (kUseTrackLength)
			outputGraph = (TGraphErrors*) thisTrack->doFit();
		else
			outputGraph = (TGraphErrors*) thisTrack->doRangeFit();
		
		if (!outputGraph) continue;

		Float_t fitEnergy = thisTrack->getFitParameterEnergy();
		Float_t fitScale = thisTrack->getFitParameterScale();
		
		Float_t fitError = quadratureAdd(thisTrack->getFitParameterError(), dz/sqrt(12));

		if (acceptAngle || true) {
			hFitResults->Fill(getUnitFromEnergy(fitEnergy));

			if (fitIdx < plotSize && kDrawIndividualGraphs) {
				drawIndividualGraphs(cGraph, outputGraph, fitEnergy, fitScale, fitError, fitIdx++);
			}
		}
		
		else delete outputGraph;
	}
	
	if (!kDrawIndividualGraphs) delete cGraph;
	
	cout << 100 * float(nCutDueToTrackEndingAbruptly) / tracks->GetEntriesFast() << " % of the tracks were cut due to seemingly inelastic nuclear interactions.\n";
	
//	cMaxAngle->cd();
//	hMaxAngle->Draw();
	TF1 *fMaxAngle = new TF1("fMaxAngle", "gaus(0)", 0, 25);
	fMaxAngle->SetParameters(100, 4, 6);
	hMaxAngle->Fit(fMaxAngle, "M, W, Q, N", "", 0, 25);

	Float_t angleTo = fMaxAngle->GetParameter(1) + 3 * fMaxAngle->GetParameter(2);

	cout << "3 sigma CL = " << angleTo << endl;

	Int_t nAccepted = hMaxAngle->Integral(0,hMaxAngle->GetXaxis()->FindBin(angleTo));
	Float_t percentAccepted = 100 * nAccepted / hMaxAngle->Integral(0);

	cout << "Number of accepted events = " << nAccepted << " of total " << hMaxAngle->Integral()  << "(" <<  percentAccepted << " %) " << endl;

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
	
	Float_t nGaussianFitRange = doNGaussianFit(hFitResults, means, sigmas);
		
	Int_t nMean = 0;
	for (Int_t i=0; i<10; i++) {
		if (means[i]) nMean++;
	}
	
	Float_t rangeSigma = 0;
	for (Int_t i=0; i<nMean; i++) {
		rangeSigma += pow(sigmas[i], 2);
	}
	rangeSigma = sqrt(rangeSigma);

	Float_t energySigma = getEnergyFromUnit(nGaussianFitRange + rangeSigma / 2) - getEnergyFromUnit(nGaussianFitRange - rangeSigma / 2);

	cFitResults->Update();
	
	TLine *l = nullptr;
	if (kDrawVerticalLayerLines) {
		Float_t line_z = 0;
		for (Int_t i=0; i<10; i++) {
			line_z = getWEPLFromTL(getLayerPositionmm(i));
			l = new TLine(line_z, 0, line_z, hFitResults->GetMaximum()*1.05);
			l->SetLineColor(kBlack); l->SetLineWidth(2); l->Draw();
		}
	}
	
	TLegend *legend = new TLegend(0.15, 0.6, 0.48, 0.85);
	legend->SetTextSize(0.02);
	legend->AddEntry(hFitResults, "Results from individual track fits", "F");
//  	legend->AddEntry(landau, Form("Fit with E = %.1f MeV and #sigma = %.1f mm ", landau_energy, landau->GetParameter(2)*1.7), "F");
	if (kDrawVerticalLayerLines) legend->AddEntry(l, "Sensor layer positions", "L");
// 	legend->Draw();

	TPaveStats *ps = (TPaveStats*) cFitResults->GetPrimitive("stats");
	hFitResults->SetBit(TH1::kNoStats);
	ps->AddText(Form("Nominal mean = %.2f", expectedMean));
	ps->AddText(Form("Nominal straggling = %.2f", expectedStraggling));
	
	for (Int_t i=0; i<nMean; i++) {
		ps->AddText(Form("WEPL fit %d = %.2f", i+1, means[i]));
		ps->AddText(Form("Energy fit %d = %.2f", i+1, getEnergyFromUnit(means[i])));
		ps->AddText(Form("Sigma fit %d = %.2f", i+1, sigmas[i]));
	}
	ps->AddText(Form("Resulting range = %.2f #pm %.2f", nGaussianFitRange, rangeSigma));
	ps->AddText(Form("Resulting energy = %.2f #pm %.2f", getEnergyFromUnit(nGaussianFitRange), energySigma));
		
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
	
	return nGaussianFitRange;
}

void writeClusterFile(Int_t Runs, Int_t dataType, Float_t energy) {
	run_energy = energy;
	
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
	TString fileName = "Data/Tracks/tracks";
	TString sMaterial = getMaterialChar();
	fileName.Append(sDataType);
	fileName.Append(sMaterial);
	fileName.Append(sEnergy);
	fileName.Append(".root");
	
	TFile f(fileName, "recreate");
	f.SetCompressionLevel(1);
	TTree T("T", "tracks");
	T.Branch("tracks", &tracks, 256000, 1);
	cout << "Length of CWOT: " << tracks->GetEntriesFastClustersWithoutTrack() << endl; 
	T.Fill();
	T.Write();
	f.Close();
}

Tracks * loadTracks(Int_t Runs, Int_t dataType, Float_t energy) {
  	TString sDataType = (dataType == 0) ? "_MC_" : "_data_";
	TString sEnergy = Form("_%.2fMeV", energy);
	TString fileName = "Data/Tracks/tracks";
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
	Tracks * tracks = nullptr;
	
	if (recreate) {
		tracks = getTracks(Runs, dataType, kCalorimeter, energy, x, y);
		if (tracks->GetEntries()) {
			cout << "Saving " << tracks->GetEntries() << " tracks.\n";
//			saveTracks(tracks, dataType, energy);
		}
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

Hits * getEventIDs(Int_t Runs, Float_t energy) {
	run_energy = energy;
	DataInterface *di = new DataInterface();

	Hits * hits = new Hits(kEventsPerRun * sizeOfEventID * Runs);

	for (Int_t i=0; i<Runs; i++) {
		di->getEventIDs(i, hits);
	}

	return hits;
}

Tracks * getTracks(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy, Float_t *x, Float_t *y) {
	run_energy = energy;

	DataInterface *di = new DataInterface();

	Int_t nClusters = kEventsPerRun * 5 * nLayers;
	Int_t nHits = kEventsPerRun * 50;
	Int_t nTracks = kEventsPerRun * 2;

	Bool_t breakSignal = false;

	CalorimeterFrame *cf = new CalorimeterFrame();
//	TrackerFrame *tf = new TrackerFrame();
	Clusters * clusters = new Clusters(nClusters);
	Clusters * trackerClusters = new Clusters(nClusters);
	Hits *hits = new Hits(nHits);
	Hits *eventIDs = new Hits(kEventsPerRun * sizeOfEventID);
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
			eventID = di->getMCFrame(i, cf, x, y);
			di->getEventIDs(i, eventIDs);
			if (kDebug) cout << "Sum frame layer 2 = " << cf->getTH2F(2)->GetSum() << endl;
			t1.Stop(); t2.Start();
			if (kDebug) cout << "Start diffuseFrame\n";
			cf->diffuseFrame(gRandom);
			cout << "Occupancy = " << 100 * cf->getOccupancyLastLayer() << "%.\n";
			if (kDebug) cout << "End diffuseFrame, start findHits\n";
			t2.Stop(); t3.Start();
			hits = cf->findHits(eventID);
			if (kDebug) cout << "Number of hits in frame: " << hits->GetEntriesFast() << endl;
			t3.Stop(); t4.Start();
			clusters = hits->findClustersFromHits(); // badly optimized
			clusters->removeSmallClusters(2);
		//	cout << "Number of clusters in frame: " << clusters->GetEntriesFast() << endl;
			cout << "Number of clusters in last layer: " << clusters->GetEntriesFastLastLayer() << endl;
			t4.Stop();
			clusters->matchWithEventIDs(eventIDs);
			eventIDs->Clear();
		}
		
		else if (dataType == kData) {
			di->getDataFrame(i, cf, energy);
			cout << "Occupancy = " << 100 * cf->getOccupancyLastLayer() << "%.\n";
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

//		calorimeterTracks->matchWithEventIDs(eventIDs);

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
		hits->clearHits();
		trackerHits->clearHits();
		clusters->clearClusters();
		trackerClusters->clearClusters();
		calorimeterTracks->clearTracks();
		trackerTracks->clearTracks();

		if (breakSignal) break;
	}

	delete cf;
	delete clusters;
	delete trackerClusters;
	delete hits;
	delete trackerHits;
	delete calorimeterTracks;
	delete trackerTracks;
	delete di;

	return allTracks;
}

void drawTracks3D(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();

	TCanvas *c1 = new TCanvas("c1");
	c1->SetTitle(Form("Tracks from %.2f MeV protons on %s", energy, getMaterialChar()));
	TView *view = TView::CreateView(1);
	view->SetRange(0, 0, 0, 2*nx, 10, 2*ny);

	TClonesArray *restPoints = tracks->getClustersWithoutTrack();
//	TClonesArray *conflictClusters = nullptr; 
	Clusters * conflictClusters = nullptr;

	Int_t nClusters = 0;
	for (Int_t i=0; i<tracks->GetEntriesFast(); i++) {
		nClusters += tracks->GetEntriesFast(i);
	}

	TPolyMarker3D *pMarker = new TPolyMarker3D(restPoints->GetEntriesFast(), 7);
	TPolyMarker3D *EIDMarker = new TPolyMarker3D(nClusters, 7);
	TPolyMarker3D *conflictMarker = new TPolyMarker3D(nClusters, 7);
   pMarker->SetMarkerColor(kBlue); // Missing cluster
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

	for (Int_t i=0; i<ntracks; i++) {

		Track *thisTrack = tracks->At(i);
		if (thisTrack->getTrackLengthmm() < 2) continue;
		Int_t n = thisTrack->GetEntriesFast();

		TPolyLine3D *l = new TPolyLine3D(n);
		l->SetLineWidth(1);

		firstEID = thisTrack->getEventID(0);
		Int_t pointNumber = 0;
		for (Int_t j=0; j<n; j++) {
			if (!thisTrack->At(j)) continue;

			Float_t x = thisTrack->getX(j);
			Float_t z = thisTrack->getY(j);
			Float_t y = thisTrack->getLayer(j);
			l->SetPoint(pointNumber++,x,y,z);

			if (thisTrack->getEventID(j) != firstEID) {
				EIDMarker->SetPoint(EIDidx++, x, y, z);
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

		l->Draw();
//		EIDMarker->Draw();
		conflictMarker->Draw();
	}
	view->ShowAxis();
   c1->Update();

   TAxis3D *axis = TAxis3D::GetPadAxis();
   axis->SetLabelColor(kBlack);
   axis->SetAxisColor(kBlack);

//   tracks->checkLayerOrientation();

   delete tracks;
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

void drawFrame2D(Int_t Runs, Int_t Layer, Float_t energy) {
   run_energy = energy;
	DataInterface *di = new DataInterface();
   CalorimeterFrame *cf = new CalorimeterFrame();
   
   TList *histogramList = new TList;
   
   for (Int_t i=0; i<Runs; i++) {
	   di->getMCFrame(i, cf);
	   histogramList->Add(cf->getTH2F(Layer));
   }
 
	delete di;

   TH2F *Frame2D = new TH2F("Frame2D", Form("Hit distribution in layer %i", Layer),
                              nx*2, 0, nx*2, ny*2, 0, ny*2);

   Frame2D->Merge(histogramList);
   Frame2D->Draw("COLZ");
   gStyle->SetOptStat(0);
}

void drawData3D(Int_t Runs, Float_t energy) {
	run_energy = energy;

   DataInterface *di = new DataInterface();

   TH3F *Frame3D = new TH3F("Frame3D", "3D map of energy deposition [keV]", 
                              100, -120, -50, 100, 0, 2*nx, 100, 0, 2*ny);

   Frame3D->SetXTitle("Z axis");
   Frame3D->SetYTitle("X axis"); // to get projection right (Z is depth, not up)
   Frame3D->SetZTitle("Y axis"); 

   for (Int_t run=0; run<Runs; run++) {
   	di->getMCData(run, Frame3D);
   }
	
   delete di;

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
	TCanvas *cCustomInv = new TCanvas("cCustomInv", "Range fit for custom range-energy list", 1200, 900);
	
	TCanvas *cCustomW = new TCanvas("cCustomW", "Range fit for custom range-energy list tungsten", 1200, 900);
	TCanvas *cCustomInvW = new TCanvas("cCustomInvW", "Range fit for custom range-energy list tungsten", 1200, 900);

	Float_t energies[19] = {10, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400};
	Float_t ranges[13] = {1.047, 4.061, 8.639, 14.655, 22.022, 45.868, 76.789, 114.074, 157.148, 205.502, 258.75, 316.419, 378.225};
	Float_t rangesW[15] = {3.027, 5.989, 9.691, 13.983, 18.899, 24.586, 30.541, 37.219, 43.945, 51.257, 59.058, 67.172, 75.522, 84.338, 93.62};

	cCustom->cd();
	TGraph *graph_custom = new TGraph(13, energies, ranges);
	graph_custom->SetTitle("Range fit for custom range-energy list");
	graph_custom->Draw("A*");

	cCustomInv->cd();
	TGraph *graph_inv = new TGraph(13, ranges, energies);
	graph_inv->SetTitle("Energy fit for custom energy-range list");
	graph_inv->Draw("A*");

	cCustomW->cd();
	TGraph *graph_custom_w = new TGraph(15, energies, rangesW);
	graph_custom_w->SetTitle("Range fit for W");
	graph_custom_w->Draw("A*");

	cCustomInvW->cd();
	TGraph *graph_inv_w = new TGraph(15, rangesW, energies);
	graph_inv_w->SetTitle("Energy for for W");
	graph_inv_w->Draw("A*");

	TF1 *fCustom = new TF1("fCustom", "[0] * x * (1 + ([1] - [1]* exp(-[2] * x)) + ([3] - [3] * exp(-[4]*x)))", 10, 250);
	TF1 *fCustomInv = new TF1("fCustomInv", "x * ([0] * exp ( - [1] * x) + [2] * exp( - [3] * x) + [4] * exp( - [5] * x) + [6] * exp(-[7] * x) + [8] * exp(-[9] * x))", 1, 400);
	TF1 *fCustomW = new TF1("fCustomW", "[0] * x * (1 + ([1] - [1]* exp(-[2] * x)) + ([3] - [3] * exp(-[4]*x)))", 10, 250);
	TF1 *fCustomInvW = new TF1("fCustomInvW", "x * ([0] * exp ( - 1./[1] * x) + [2] * exp( - 1./[3] * x) + [4] * exp( - 1./[5] * x) + [6] * exp(-1./[7] * x) + [8] * exp(-1./[9] * x))", 1, 400);
	fCustom->SetParameters(6.94656e-2, 15.14450027, 0.001260021, 29.84400076, 0.003260031);
	fCustomInv->SetParameters(9.663872, 1/0.975, 2.50472, 1/12.4999, 0.880745, 1/57.001, 0.419001, 1/106.501, 0.92732, 1/1067.2784);
	
	fCustomW->SetParameters(6.94656e-2/9, 15.14450027, 0.001260021, 29.84400076, 0.003260031);
	fCustomInvW->SetParameters(9.663872*9, 0.975/9, 2.50472*9, 12.4999/9, 0.880745*9, 57.001/9, 0.419001*9, 106.501/9, 0.92732*9, 1067.2784/9);
	
	fCustom->SetParLimits(0, 0.02, 0.1);
	fCustom->SetParLimits(1, 5, 45);
	fCustom->SetParLimits(2, 0.0005, 0.0045);
	fCustom->SetParLimits(3, 10, 90);
	fCustom->SetParLimits(4, 0.001, 0.01);
	
	fCustomInv->SetParLimits(0, 0.1, 100);
	fCustomInv->SetParLimits(1, 0.001, 5);
	fCustomInv->SetParLimits(2, 0.1, 100);
	fCustomInv->SetParLimits(3, 0.001, 5);
	fCustomInv->SetParLimits(4, 0.1, 100);
	fCustomInv->SetParLimits(5, 0.001, 5);
	fCustomInv->SetParLimits(6, 0.1, 100);
	fCustomInv->SetParLimits(7, 0.001, 5);
	fCustomInv->SetParLimits(8, 0.1, 100);
	fCustomInv->SetParLimits(9, 0.001, 5);

	fCustomInvW->SetParLimits(0, 1, 50);
	fCustomInvW->SetParLimits(1, 0.01, 0.2);
	fCustomInvW->SetParLimits(2, 5, 30);
	fCustomInvW->SetParLimits(3, 1, 8);
	fCustomInvW->SetParLimits(4, 1, 15);
	fCustomInvW->SetParLimits(5, 5, 50);
	fCustomInvW->SetParLimits(6, 0.1, 10);
	fCustomInvW->SetParLimits(7, 0.5, 30);
	fCustomInvW->SetParLimits(8, 1, 20);
	fCustomInvW->SetParLimits(9, 30, 500);

	cCustom->cd();
	fCustom->SetNpx(500);
	graph_custom->Fit("fCustom", "M, B");
	cCustomInv->cd();
	fCustomInv->SetNpx(500);
	graph_inv->Fit("fCustomInv", "M, B");
	cCustomW->cd();
	fCustomW->SetNpx(500);
	graph_custom_w->Fit("fCustomW", "M, B");
	cCustomInvW->cd();
	fCustomInvW->SetNpx(500);
	graph_inv_w->Fit("fCustomInvW", "M, B");

	Float_t sqrt1 = 0, sqrt2 = 0, sqrt3 = 0, sqrt4 = 0; 
	Float_t avgRangeSum = 0, avgRangeSumW = 0;
	for (Int_t i=0; i<9; i++) {
		sqrt1 += pow(ranges[i] - fCustom->Eval(energies[i]), 2);
		sqrt2 += pow(energies[i] - fCustomInv->Eval(ranges[i]), 2);
		sqrt3 += pow(rangesW[i] - fCustomW->Eval(energies[i]), 2);
		sqrt4 += pow(energies[i] - fCustomInvW->Eval(rangesW[i]), 2);
		avgRangeSum += fabs(ranges[i] - fCustom->Eval(energies[i]));
		avgRangeSumW += fabs(rangesW[i] - fCustomW->Eval(energies[i]));
	}

	avgRangeSum /= 9;
	avgRangeSumW /= 15;

	sqrt1 = sqrt(sqrt1);
	sqrt2 = sqrt(sqrt2);
	sqrt3 = sqrt(sqrt3);
	sqrt4 = sqrt(sqrt4);

	cout << "RSQ of Gompertz-function = " << sqrt1 << endl;
	cout << "RSQ of INV fitted Gompertz-function = " << sqrt2 << endl;
	cout << "RSQ of Gompertz-function in tungsten = " << sqrt3 << endl;
	cout << "RSQ of INV fitted Gompertz-function in tungsten = " << sqrt4 << endl;
	cout << endl << "Average range error for fitted Gompertz is " << avgRangeSum << " mm.\n";
	cout << "Average range error for fitted Gompertz in tungsten is " << avgRangeSumW << " mm.\n";

	Float_t E1 = fCustomInv->Eval(fCustom->Eval(10));
	Float_t E2 = fCustomInv->Eval(fCustom->Eval(30));
	Float_t E3 = fCustomInv->Eval(fCustom->Eval(50));

	Float_t E1w = fCustomInvW->Eval(fCustomW->Eval(10));
	Float_t E2w = fCustomInvW->Eval(fCustomW->Eval(30));
	Float_t E3w = fCustomInvW->Eval(fCustomW->Eval(50));

	cout << Form("From energy of 10, 30, 50 MeV, the calc-invcalc values are %.10f, %.10f, %.10f.\n", E1, E2, E3);
	cout << Form("TUNGSTEN - From energy of 10, 30, 50 MeV, the calc-invcalc values are %.10f, %.10f, %.10f.\n", E1w, E2w, E3w);

	cout << " --- \033[1m VALUES FOR WATER \033[0m --- " << endl;
	cout << "a1      = " << fCustom->GetParameter(0) << endl;
	cout << "b1      = " << fCustom->GetParameter(1) << endl;
	cout << "g1      = " << fCustom->GetParameter(2) << endl;
	cout << "b2      = " << fCustom->GetParameter(3) << endl;
	cout << "g2      = " << fCustom->GetParameter(4) << endl;
	cout << endl;
	cout << "c1      = " << fCustomInv->GetParameter(0) << endl;
	cout << "lambda1 = " << 1/fCustomInv->GetParameter(1) << endl;
	cout << "c2      = " << fCustomInv->GetParameter(2) << endl;
	cout << "lambda2 = " << 1/fCustomInv->GetParameter(3) << endl;
	cout << "c3      = " << fCustomInv->GetParameter(4) << endl;
	cout << "lambda3 = " << 1/fCustomInv->GetParameter(5) << endl;
	cout << "c4      = " << fCustomInv->GetParameter(6) << endl;
	cout << "lambda4 = " << 1/fCustomInv->GetParameter(7) << endl;
	cout << "c5      = " << fCustomInv->GetParameter(8) << endl;
	cout << "lambda5 = " << 1/fCustomInv->GetParameter(9) << endl;

	cout << endl << endl;

	cout << " --- \033[1m VALUES FOR TUNGSTEN \033[0m --- " << endl;
	cout << "a1      = " << fCustomW->GetParameter(0) << endl;
	cout << "b1      = " << fCustomW->GetParameter(1) << endl;
	cout << "g1      = " << fCustomW->GetParameter(2) << endl;
	cout << "b2      = " << fCustomW->GetParameter(3) << endl;
	cout << "g2      = " << fCustomW->GetParameter(4) << endl;
	cout << endl;
	cout << "c1      = " << fCustomInvW->GetParameter(0) << endl;
	cout << "lambda1 = " << 1/fCustomInvW->GetParameter(1) << endl;
	cout << "c2      = " << fCustomInvW->GetParameter(2) << endl;
	cout << "lambda2 = " << 1/fCustomInvW->GetParameter(3) << endl;
	cout << "c3      = " << fCustomInvW->GetParameter(4) << endl;
	cout << "lambda3 = " << 1/fCustomInvW->GetParameter(5) << endl;
	cout << "c4      = " << fCustomInvW->GetParameter(6) << endl;
	cout << "lambda4 = " << 1/fCustomInvW->GetParameter(7) << endl;
	cout << "c5      = " << fCustomInvW->GetParameter(8) << endl;
	cout << "lambda5 = " << 1/fCustomInvW->GetParameter(9) << endl;

}

void drawIndividualGraphs(TCanvas *cGraph, TGraphErrors* outputGraph, Float_t fitEnergy, Float_t fitScale, Float_t fitError, Int_t fitIdx, Int_t eventID, Float_t *x_energy, Float_t *y_energy) {
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
		TLatex *text2 = new TLatex(13, 450, Form("'Real' energy: %.1f #pm MeV", realEnergy));
		text2->SetTextSize(0.06);
		text2->Draw();
	}

	if (kDrawText) {
		TLatex *text = new TLatex(10, 500, Form("Fitted energy: %.1f #pm %.1f MeV", fitEnergy, fitError));
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
		
//		TLine *l1 = new TLine(searchFrom, 0, searchFrom, 1000); l1->SetLineColor(kGreen); l1->Draw();
//		TLine *l2 = new TLine(searchTo, 0, searchTo, 1000); l2->SetLineColor(kRed); l2->Draw();
		
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
		gauss->SetParLimits(1, searchFrom+8, searchTo-8);
		gauss->SetParLimits(2, 2, 12);
		
		h->Fit(gauss, "M, B, WW, Q", "", searchFrom, searchTo);
		
		Float_t chi2 = gauss->GetChisquare();
		Float_t chi2n = chi2 / integral;
		
		sigma = gauss->GetParameter(2);
		constant = gauss->GetParameter(0);
		mean = gauss->GetParameter(1);
		lEnergy = getEnergyFromUnit(mean);
		
 		cout << Form("Searching from %.1f to %.1f, with midpoint at %.1f. Found best fit @ %.1f with chi2 = %.2f and chi2/n = %.2f, ratio = %.2f.\n", searchFrom, searchTo,(searchTo+searchFrom)/2 , mean, chi2, chi2n, ratio);
//
//  		if (chi2n > 7) {
//  			delete gauss;
//  			continue;
//		}

 		if (ratio > 0.05 || (isLastLayer && ratio>0.025)) {
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
	Float_t estimated_energy_error = 0;
	Float_t sum_constant = 0, sumSigma = 0;
	for (Int_t i=0 ; i<3; i++) {
		estimated_range += array_constant[i] * array_mean[i];
		sum_constant += array_constant[i];
		sumSigma += pow(array_sigma[i], 2);
	}

	sumSigma = sqrt(sumSigma);

	estimated_range /= sum_constant;

	estimated_energy = getEnergyFromUnit(estimated_range);
	estimated_energy_error = getEnergyFromWEPL(estimated_range + sumSigma/2) - getEnergyFromWEPL(estimated_range - sumSigma/2);
	cout << "ESTIMATED ENERGY FROM RUN IS " << estimated_energy << " +- " << estimated_energy_error << endl;
	cout << "Estimated range = " << estimated_range << " +- " << sumSigma << endl;
	
	if (true) {
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
	}
	return estimated_range;
}
