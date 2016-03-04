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

void drawBraggPeakFit(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	run_energy = energy;
	kOutputUnit = kPhysical;
	
	char * sDataType = getDataTypeChar(dataType);
	char * sMaterial = getMaterialChar();
	char * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", energy, sMaterial, sDataType);

	Bool_t isFitOk = false;

	TCanvas *cFitResults = new TCanvas("cFitResults", hTitle, 1400, 1000);
	TH1F *hFitResults = 0;
	
	if (kOutputUnit == kPhysical) {
		hFitResults = new TH1F("fitResult", hTitle, 1000, 0, 60);
	}
	else if (kOutputUnit == kWEPL) {
		hFitResults = new TH1F("fitResult", hTitle, 1000, 50, 350);
	}
	else if (kOutputUnit == kEnergy) {
		hFitResults = new TH1F("fitResult", hTitle, 1000, 100, 300);
	}
		
	hFitResults->SetLineColor(kBlack);
	hFitResults->SetFillColor(kGreen-5);
	
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();

	for (Int_t j=0; j<tracks->GetEntriesFast(); j++) {
		Track *thisTrack = tracks->At(j);
		if (!thisTrack) continue;
	
		isFitOk = thisTrack->doFit();
		if (!isFitOk) continue;
			
		Float_t fitEnergy = thisTrack->getFitParameterEnergy();
		Float_t fitScale = thisTrack->getFitParameterScale();

		if (fitEnergy < run_energy*1.245) {
			if (kOutputUnit == kPhysical) {
				hFitResults->Fill(getTLFromEnergy(fitEnergy));
			}
			else if (kOutputUnit == kWEPL) {
				hFitResults->Fill(getWEPLFromEnergy(fitEnergy));
			}
			else if (kOutputUnit == kEnergy) {
				hFitResults->Fill(fitEnergy);
			}
		}
	}
	
	if (kOutputUnit == kPhysical) hFitResults->SetXTitle("Physical range [mm]");
	else if (kOutputUnit == kWEPL) hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
	else if (kOutputUnit == kEnergy) { hFitResults->SetXTitle("Energy [MeV]"); }
		
	hFitResults->SetYTitle("Number of protons");
	hFitResults->Draw();
	
	Float_t maxBin = hFitResults->GetMaximum();
	Int_t nBins = hFitResults->GetSize() - 2;
	Int_t nextMaxBin = 0;
	for (Int_t i=0; i<nBins; i++) {
		Float_t bin = hFitResults->GetBinContent(i);
		if (bin > nextMaxBin && bin < maxBin) {
			nextMaxBin = bin;
		}
	}
	
	// Draw expected gaussian distribution of results from initial energy
	Float_t expectedStraggling, expectedMean;
	
	if (kOutputUnit == kPhysical) {
		expectedMean = getTLFromEnergy(energy);
		expectedStraggling = getRangeStragglingFromEnergy(energy);
	}
	else if (kOutputUnit == kEnergy) {
		expectedMean = energy;
		expectedStraggling = getEnergyStragglingFromEnergy(expectedMean);
	}
	else if (kOutputUnit == kWEPL) {
		expectedMean = getWEPLFromEnergy(energy);
		expectedStraggling = getWEPLStragglingFromEnergy(energy);
	}
	
	cout << "Expected mean = " << expectedMean << endl;
	cout << "Expected stragglig = " << expectedStraggling << endl;
	
	
	TF1 *ES = new TF1("ES", "gaus(0)", expectedMean*0.85, expectedMean*1.15);
	ES->SetParameters(maxBin,expectedMean,expectedStraggling);
 	ES->SetParLimits(0,maxBin,maxBin);
 	ES->SetParLimits(1,expectedMean*0.85,expectedMean*1.15);
 	ES->SetParLimits(2,expectedStraggling,expectedStraggling);
 	hFitResults->Fit(ES,"B,Q", "", expectedMean*0.85, expectedMean*1.15);
	ES->Draw("same");
	cFitResults->Update();

	TLegend *legend = new TLegend(0.65, 0.6, 0.98, 0.75);
	legend->SetTextSize(0.03);
	legend->AddEntry(hFitResults, "Energy fits to tracks", "F");
	legend->AddEntry(ES, "Expected range straggling", "F");
	legend->Draw();
	
	Float_t lowerRange = expectedMean - 2*expectedStraggling;
	Float_t higherRange = expectedMean + 2*expectedStraggling;
	
	Float_t sigma_fwhm = getFWxMInRange(hFitResults, lowerRange, higherRange, 2);
	Float_t sigma_fwtm = getFWxMInRange(hFitResults, lowerRange, higherRange, 3);
	Float_t sigma_fwqm = getFWxMInRange(hFitResults, lowerRange, higherRange, 4);
	
	TPaveStats *ps = (TPaveStats*) cFitResults->GetPrimitive("stats");
	hFitResults->SetBit(TH1::kNoStats);
	ps->AddText(Form("Nominal mean = %.2f", expectedMean));
	ps->AddText(Form("Nominal straggling = %.2f", expectedStraggling));
	ps->AddText(Form("Straggling (from FWHM) = %.2f", sigma_fwhm));
	ps->AddText(Form("Straggling (from FWTM) = %.2f", sigma_fwtm));
	ps->AddText(Form("Straggling (from FWQM) = %.2f", sigma_fwqm));
	
	// draw lines on HM, TM, QM
	Float_t len_hm = expectedStraggling * 2.355/2;
	Float_t len_tm = expectedStraggling * 2.9646/2;
	Float_t len_qm = expectedStraggling * 3.330/2;
	
	Float_t width = hFitResults->GetXaxis()->GetBinWidth(1)/2;
	TLine *hm = new TLine(expectedMean-sigma_fwhm/2-width, maxBin/2, expectedMean + sigma_fwhm/2 + width, maxBin/2);
	TLine *tm = new TLine(expectedMean-sigma_fwtm/2-width, maxBin/3, expectedMean + sigma_fwtm/2 + width, maxBin/3);
	TLine *qm = new TLine(expectedMean-sigma_fwqm/2-width, maxBin/4, expectedMean + sigma_fwqm/2 + width, maxBin/4);
	
	hm->SetLineColor(kRed); tm->SetLineColor(kRed); qm->SetLineColor(kRed);
	hm->SetLineWidth(3); tm->SetLineWidth(3); qm->SetLineWidth(3);
	hm->Draw(); tm->Draw(); qm->Draw();

	Bool_t kOutput = false;
	
	if (kOutput) {
		ofstream outputFile;
		outputFile.open("fit_result.txt", ios::out | ios::app);
		outputFile << energy << "; " << (int) kOutputUnit << "; " << expectedMean << "; " << hFitResults->GetMean() << "; " 
				<< expectedStraggling << "; " << sigma_fwhm << "; " << sigma_fwtm << "; "
				<< sigma_fwqm << endl;
		outputFile.close();
	}
	
	cFitResults->Modified();
	
	if (kOutputUnit == kPhysical) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_physical.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	else if (kOutputUnit == kEnergy) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_energy.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	
	else if (kOutputUnit == kWEPL) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_WEPL.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
  	delete tracks;
}

void drawBraggPeakGraphFit(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	run_energy = energy;
//	kOutputUnit = kEnergy;
	
	char * sDataType = getDataTypeChar(dataType);
	char * sMaterial = getMaterialChar();
	char * hTitle = Form("Fitted energy of a %.2f MeV beam in %s (%s)", energy, sMaterial, sDataType);

	Int_t nPlotX = 4, nPlotY = 4;
	Int_t fitIdx = 0, plotSize = nPlotX*nPlotY;
	Bool_t isFitOk = false;

	TCanvas *cGraph = new TCanvas("cGraph", "Fitted data points", 1400, 1000);
	TCanvas *cFitResults = new TCanvas("cFitResults", hTitle, 1400, 1000);
	TCanvas *cScale = new TCanvas("cScale", "Scale histogram", 1400, 1000);
	cGraph->Divide(nPlotX,nPlotY, 0.000001, 0.000001, 0);
	
	TH1F *hFitResults = 0;
	if (kOutputUnit == kPhysical) {
		if (kMaterial == kTungsten) 
			hFitResults = new TH1F("fitResult", hTitle, 1000, 0, 60);
		else if (kMaterial == kAluminum)
			hFitResults = new TH1F("fitResult", hTitle, 1000, 0, 200);
		else
			hFitResults = new TH1F("fitResult", hTitle, 1000, 0, 400);

	}
	else if (kOutputUnit == kWEPL) {
		hFitResults = new TH1F("fitResult", hTitle, 1000, 50, 350);
	}
	else if (kOutputUnit == kEnergy) {
		hFitResults = new TH1F("fitResult", hTitle, 1000, 100, 300);
	}
	else {
		cout << "Unknown output unit!\n";
	}

	TH1F *hScale = new TH1F("hScale", "Scale histogram", 400, 0, 400);
	
	hFitResults->SetLineColor(kBlack);
	hFitResults->SetFillColor(kGreen-5); 
	
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
		
		hScale->Fill(fitScale);

		if (fitEnergy < run_energy*1.245) {
			if (kOutputUnit == kPhysical) {
				hFitResults->Fill(getTLFromEnergy(fitEnergy));
			}
			else if (kOutputUnit == kWEPL) {
				hFitResults->Fill(getWEPLFromEnergy(fitEnergy));
			}
			else if (kOutputUnit == kEnergy) {
				hFitResults->Fill(fitEnergy);
			}
		}

		Int_t n = thisTrack->GetEntriesFast();
		Float_t x[n], y[n];

		if (fitIdx < plotSize) {
			cGraph->cd(fitIdx+1);
			
			if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
				Float_t trackLengthWEPL = 0;
				for (Int_t k=0; k<n; k++) {
					if (!thisTrack->At(k)) continue;
					trackLengthWEPL += thisTrack->getTrackLengthWEPLmmAt(k);
					x[k] = trackLengthWEPL;
					y[k] = thisTrack->getDepositedEnergy(k);
				}
			}
			
			else if (kOutputUnit == kPhysical) {
				Float_t trackLength = 0;
				for (Int_t k=0; k<n; k++) {
					if (!thisTrack->At(k)) continue;
					trackLength += thisTrack->getTrackLengthmmAt(k);
					x[k] = trackLength;
					y[k] = thisTrack->getDepositedEnergy(k);
				}
			}
			
			TGraph *graph = new TGraph(n,x,y);

			graph->SetMinimum(0);
			graph->SetMaximum(700);
			graph->SetTitle("");
			graph->GetXaxis()->SetTitle("Water Equivalent Path Length [mm]");
			graph->GetYaxis()->SetTitle("Deposited energy per layer [keV]");
			graph->GetXaxis()->SetTitleSize(0.05);
			graph->GetYaxis()->SetTitleSize(0.05);
			graph->GetXaxis()->SetLabelSize(0.04);
			graph->GetYaxis()->SetLabelSize(0.04);
			graph->GetXaxis()->SetLimits(0, 350);
			
			graph->Draw("A*");
			
			TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
			func->SetParameters(fitEnergy, fitScale);
			func->Draw("same");

			TLatex *text = new TLatex(20, 600, Form("Fitted energy: %.1f MeV", fitEnergy));
			text->SetTextSize(0.06);
			text->Draw();
		
			cGraph->Update();
			fitIdx++;
		}
	}
	
	cScale->cd();
	hScale->SetYTitle("Number of protons");
	hScale->SetXTitle("Parameter 1 of fit (SCALE)");
	hScale->Draw();
	
	cFitResults->cd();
	if (kOutputUnit == kPhysical) hFitResults->SetXTitle("Physical range [mm]");
	else if (kOutputUnit == kWEPL) hFitResults->SetXTitle("Range in Water Equivalent Path Length [mm]");
	else if (kOutputUnit == kEnergy) { hFitResults->SetXTitle("Energy [MeV]"); }
		
	hFitResults->SetYTitle("Number of protons");
	hFitResults->Draw();
	
	Float_t maxBin = hFitResults->GetMaximum();
	Int_t nBins = hFitResults->GetSize() - 2;
	Int_t nextMaxBin = 0;
	for (Int_t i=0; i<nBins; i++) {
		Float_t bin = hFitResults->GetBinContent(i);
		if (bin > nextMaxBin && bin < maxBin) {
			nextMaxBin = bin;
		}
	}
	
	// Draw expected gaussian distribution of results from initial energy
	
	Float_t expectedStraggling = 0, expectedMean = 0;
	
	if (kOutputUnit == kPhysical) {
		expectedMean = getTLFromEnergy(energy);
		expectedStraggling = getRangeStragglingFromEnergy(energy);
	}
	else if (kOutputUnit == kEnergy) {
		expectedMean = energy;
		expectedStraggling = getEnergyStragglingFromEnergy(expectedMean);
	}
	else if (kOutputUnit == kWEPL) {
		expectedMean = getWEPLFromEnergy(energy);
		expectedStraggling = getWEPLStragglingFromEnergy(energy);
	}
	
	TF1 *ES = new TF1("ES", "gaus(0)", expectedMean*0.85, expectedMean*1.15);
	ES->SetParameters(maxBin,expectedMean,expectedStraggling);
 	ES->SetParLimits(0,maxBin,maxBin);
 	ES->SetParLimits(1,expectedMean,expectedMean);
 	ES->SetParLimits(2,expectedStraggling,expectedStraggling);
 	hFitResults->Fit(ES,"B,Q", "", expectedMean*0.85, expectedMean*1.15);
	cout << "Tryna draw Gauss(" << maxBin << "," << expectedMean << "," << expectedStraggling << ").\n";
	ES->Draw("same");
	cFitResults->Update();

	TLegend *legend = new TLegend(0.65, 0.6, 0.98, 0.75);
	legend->SetTextSize(0.03);
	legend->AddEntry(hFitResults, "Energy fits to tracks", "F");
	legend->AddEntry(ES, "Expected range straggling", "F");
	legend->Draw();
	
	// find Full Width Half Maximum and draw it in the TPaveStats
	Int_t bin1 = hFitResults->FindFirstBinAbove(maxBin/2);
	Int_t bin2 = hFitResults->FindLastBinAbove(maxBin/2);
	Float_t fwhm = hFitResults->GetBinCenter(bin2) - hFitResults->GetBinCenter(bin1);
	
	Int_t bin3 = hFitResults->FindFirstBinAbove(maxBin/4);
	Int_t bin4 = hFitResults->FindLastBinAbove(maxBin/4);
	Float_t fwqm = hFitResults->GetBinCenter(bin4) - hFitResults->GetBinCenter(bin3);
	
	Int_t bin5 = hFitResults->FindFirstBinAbove(maxBin/3);
	Int_t bin6 = hFitResults->FindLastBinAbove(maxBin/3);
	Float_t fwtm = hFitResults->GetBinCenter(bin6) - hFitResults->GetBinCenter(bin5);
	
	TPaveStats *ps = (TPaveStats*) cFitResults->GetPrimitive("stats");
	hFitResults->SetBit(TH1::kNoStats);
	ps->AddText(Form("Nominal mean = %.2f", expectedMean));
	ps->AddText(Form("Nominal straggling = %.2f", expectedStraggling));
	//ps->AddText(Form("FWHM = %.2f", fwhm));
	ps->AddText(Form("Straggling (from FWHM) = %.2f", fwhm/2.355));
	ps->AddText(Form("Straggling (from FWTM) = %.2f", fwtm/2.9646));
	ps->AddText(Form("Straggling (from FWQM) = %.2f", fwqm/3.33));
	
	// draw lines on HM, TM, QM
	Float_t len_hm = expectedStraggling * 2.355/2;
	Float_t len_tm = expectedStraggling * 2.9646/2;
	Float_t len_qm = expectedStraggling * 3.330/2;
	
	//TLine *hm = new TLine(expectedMean - len_hm, maxBin/2, expectedMean + len_hm, maxBin/2);
	TLine *hm = new TLine(hFitResults->GetXaxis()->GetBinLowEdge(bin1), maxBin/2, hFitResults->GetXaxis()->GetBinUpEdge(bin2), maxBin/2);
	TLine *tm = new TLine(hFitResults->GetXaxis()->GetBinLowEdge(bin5), maxBin/3, hFitResults->GetXaxis()->GetBinUpEdge(bin6), maxBin/3);
	TLine *qm = new TLine(hFitResults->GetXaxis()->GetBinLowEdge(bin3), maxBin/4, hFitResults->GetXaxis()->GetBinUpEdge(bin4), maxBin/4);
	
	hm->SetLineColor(kRed); tm->SetLineColor(kRed); qm->SetLineColor(kRed);
	hm->SetLineWidth(3); tm->SetLineWidth(3); qm->SetLineWidth(3);
	hm->Draw(); tm->Draw(); qm->Draw();
	
	cFitResults->Modified();
	
	if (kOutputUnit == kPhysical) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_physical.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	else if (kOutputUnit == kEnergy) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_energy.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	
	else if (kOutputUnit == kWEPL) {
		cFitResults->SaveAs(Form("figures/Fitted_energies_%.2f_MeV_%s_%s_WEPL.pdf", energy, getMaterialChar(), getDataTypeChar(dataType)));
	}
	
  	delete tracks;
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

void Draw2DProjection(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {

	Tracks * tracks = loadOrCreateTracks(recreate, Runs, kCalorimeter, energy);
	tracks->extrapolateToLayer0();
	
	Int_t hSizeX = nx/8;
	Int_t hSizeY = ny/8;

	Float_t angleCut = 5.;
	Float_t x0, y0, theta0;
	Int_t nPoints = 0;

	Float_t fit_energy;

	Int_t ntracks = tracks->GetEntriesFast();

	char *title = (char*) Form("2D Projection of data from %.2f MeV proton beam on %s type FOCAL", energy, getMaterialChar());

	TCanvas *c1 = new TCanvas(title);
	
	TH2F *Frame2D = new TH2F("Frame2D", title, hSizeX, 0, nx*2, hSizeY, 0, ny*2);
	TH2F *normalizeFrame = new TH2F("normalizeFrame", "title", hSizeX, 0, nx*2, hSizeY, 0, ny*2);

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

Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy) {
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

Tracks * getTracks(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy) {
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
			if (kDebug) cout << "Sum frame layer 2 = " << cf->getTH2F(2)->GetSum() << endl;
			t1.Stop(); t2.Start();
			if (kDebug) cout << "Start diffuseFrame\n";
			cf->diffuseFrame(new TRandom3(0));
			if (kDebug) cout << "End diffuseFrame, start findHits\n";
			t2.Stop(); t3.Start();
			hits = cf->findHits();
			if (kDebug) cout << "Number of hits in frame: " << hits->GetEntriesFast() << endl;
			t3.Stop(); t4.Start();
			clusters = hits->findClustersFromHits(); // badly optimized
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
  
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, kCalorimeter, energy);
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
   	f->GetMCFrame(i, cf); // Remember to have MC data available at ./test.root
   }

   TCanvas *c1 = new TCanvas("c1", "multipads", 1400, 900);
   gStyle->SetOptStat(0);
   c1->Divide(2,1,0.01,0.01,0);
   
   TH2F *undiffusedTH2F = (TH2F*) cf->GetTH2F(Layer)->Clone();
   undiffusedTH2F->SetName("undiffusedTH2F");
   
   cf->DiffuseFrame(new TRandom3(0));
   
   TH2F *diffusedTH2F = (TH2F*) cf->GetTH2F(Layer)->Clone();
   diffusedTH2F->SetName("diffusedTH2F");

   c1->cd(1);
   undiffusedTH2F->Draw("COLZ");
   
   c1->cd(1);
   diffusedTH2F->Draw("COLZ");
   
   c1->Update();
}

void DrawFrame2D(Int_t Runs, Int_t Layer) {
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
	Int_t n = nPLEnergies;

	TCanvas *cW = new TCanvas("cW", "Range fit for tungsten", 1200, 900);
	TCanvas *cAl = new TCanvas("cAl", "Range fit for aluminum", 1200, 900);
	TCanvas *cPMMA = new TCanvas("cPMMA", "Range fit for PMMA", 1200, 900);
	TCanvas *cwater = new TCanvas("cwater", "Range fit for water", 1200, 900);
	TCanvas *cCustom = new TCanvas("cCustom", "Range fit for custom range-energy list", 1200, 900);

	Float_t kPLFocal_H2O[n];
	for (Int_t i=0; i<n; i++) {
		kPLFocal_H2O[i] = kPLFocal_Al[i] * kWEPLRatio_Al[i];
	}

	const Int_t n2 = 16;
	Float_t energies[n2] =
				{70, 90, 110, 130, 150, 170, 180, 200,
				 225, 250, 275, 300, 325, 350, 375, 400};

// G4 ranges for pure Tungsten
//	Float_t ranges[n2] = 
//				{4.42, 6.83, 9.62, 12.8, 16.3, 20.1,
//				 22.1, 26.3, 32.0, 38.0, 44.4, 51.0, 58.2, 65.3, 72.9, 80.6};

// G4 ranges for pure Aluminium
//	Float_t ranges[n2] = 
//				{19.74, 30.788, 43.7312, 58.546, 75.122, 93.277, 102.893, 123.255,
//					150.405, 180.101, 210.362, 242.602, 276.825, 312.74, 350.071, 387.685};

	// PSTAR CSDA ranges water
	Float_t ranges[n2] = 
		{40.8, 63.98, 91.4, 122.8, 157.7, 196.1, 216.5, 259.6, 317.4, 379.4, 445.2, 514.5, 587.1, 662.8, 741.3, 822.5};

	// G4 ranges for pure Aluminium
//	Float_t ranges[n2] = 
//		{40.48, 65.8, 91.22, 122.74, 157.334, 195.919, 216.203, 259.057, 316.699,
//		 381.18, 446.793, 513.254, 586.785, 656.357, 739.504, 827.741};

// PSTAR CSDA ranges for pure Tungsten
//	Float_t ranges[n2] = 
//				{0.18, 1.08, 2.56, 4.52, 6.92, 9.72, 12.9, 16.4, 20.2,
//				 22.3, 26.5, 32.2, 38.3, 44.7, 51.5, 58.5, 65.9, 72.4, 81.2};


	TGraph *graph_W = new TGraph(n,kPLEnergies,kPLFocal_W);
	TGraph *graph_Al = new TGraph(n,kPLEnergies,kPLFocal_Al);
	TGraph *graph_PMMA = new TGraph(n,kPLEnergies,kPLFocal_PMMA);
	TGraph *graph_water = new TGraph(n,kPLEnergies, kPLFocal_H2O);
	TGraph *graph_custom = new TGraph(7, energies, ranges);

	graph_W->SetTitle("Range fit for tungsten");
	graph_Al->SetTitle("Range fit for aluminum");
	graph_PMMA->SetTitle("Range fit for PMMA");
	graph_water->SetTitle("Range fit for water");
	graph_custom->SetTitle("Range fit for custom range-energy list");


	cW->cd();
		graph_W->Draw("A*");

	cAl->cd();
		graph_Al->Draw("A*");

	cPMMA->cd();
		graph_PMMA->Draw("A*");

	cwater->cd();
		graph_water->Draw("A*");

	cCustom->cd();
		graph_custom->Draw("A*");

	TF1 *fW = new TF1("fW", "[0]*pow(x,[1])",0,400);
	TF1 *fAl = new TF1("fAl", "[0]*pow(x,[1])",0,400);
	TF1 *fPMMA = new TF1("fPMMA", "[0]*pow(x,[1])",0,400);
	TF1 *fwater = new TF1("fwater", "[0]*pow(x,[1])",0,400);
	TF1 *fCustom = new TF1("fCustom", "[0]*pow(x, [1])",0,400);

	fW->SetParameter(0, 0.01); fW->SetParameter(1, 1.7);
	fAl->SetParameter(0, 0.01); fAl->SetParameter(1, 1.7);
	fPMMA->SetParameter(0, 0.01); fPMMA->SetParameter(1, 1.7);
	fwater->SetParameter(0, 0.01); fwater->SetParameter(1, 1.7);
	fCustom->SetParameter(0, 0.01); fCustom->SetParameter(1, 1.7);

	graph_W->Fit("fW");
	graph_Al->Fit("fAl");
	graph_PMMA->Fit("fPMMA");
	graph_water->Fit("fwater");
	graph_custom->Fit("fCustom");

	Float_t a1 = fW->GetParameter(0);
	Float_t p1 = fW->GetParameter(1);
	Float_t a2 = fAl->GetParameter(0);
	Float_t p2 = fAl->GetParameter(1);
	Float_t a3 = fPMMA->GetParameter(0);
	Float_t p3 = fPMMA->GetParameter(1);
	Float_t a4 = fwater->GetParameter(0);
	Float_t p4 = fwater->GetParameter(1);
	Float_t a5 = fCustom->GetParameter(0);
	Float_t p5 = fCustom->GetParameter(1);

	cout << Form("For tungsten, the range fit is R = %.6f * E ^ %.6f\n", a1, p1);
	cout << Form("For aluminum, the range fit is R = %.6f * E ^ %.6f\n", a2, p2);
	cout << Form("For PMMA, the range fit is R = %.6f * E ^ %.6f\n", a3, p3);
	cout << Form("For water, the range fit is R = %.6f * E ^ %.6f\n", a4, p4);
	cout << Form("For custom, the range fit is R = %.6f * E ^ %.6f\n", a5, p5);

}

