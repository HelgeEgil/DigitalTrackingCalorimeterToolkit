#define getTracks_cxx

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TRandom3.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TView.h>
#include <TLeaf.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "Classes/Track/conversionFunctions.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/DataInterface/DataInterface.h"
#include "HelperFunctions/Tools.h"
#include "HelperFunctions/getTracks.h"

#ifdef USEDEBUG
#define showDebug(x) std::cout << x
#else
#define showDebug(x)
#endif

using namespace std;

void makeTracks(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy) {
	Tracks * tracks = loadOrCreateTracks(recreate, Runs, dataType, energy);
	tracks->extrapolateToLayer0();
}

void saveTracks(Tracks *tracks, Int_t dataType, Float_t energy) {
	TString sDataType = (dataType == 0) ? "_MC_" : "_data_";

	tracks->CompressCWT();

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
			saveTracks(tracks, dataType, energy);
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

Tracks * getTracks(Int_t Runs, Int_t dataType, Int_t frameType, Float_t energy, Float_t *x, Float_t *y) {
	run_energy = energy;

	DataInterface *di = new DataInterface();

	Int_t nClusters = kEventsPerRun * 5 * nLayers;
	Int_t nHits = kEventsPerRun * 50;
	Int_t nTracks = kEventsPerRun * 2;

	Bool_t breakSignal = false;

	CalorimeterFrame *cf = new CalorimeterFrame();
	Clusters * clusters = new Clusters(nClusters);
	Clusters * trackerClusters = new Clusters(nClusters);
	Hits *hits = new Hits(nHits);
	Hits *eventIDs = new Hits(kEventsPerRun * sizeOfEventID);
	Int_t eventID = -1;
	Hits *trackerHits = new Hits(nHits);
	Tracks *calorimeterTracks = new Tracks(nTracks);
	Tracks *trackerTracks = new Tracks(nTracks);
	Tracks *allTracks = new Tracks(nTracks * Runs);
	TRandom3 *gRandom = new TRandom3(0);

	TStopwatch t1, t2, t3, t4, t5, t6;
	
	for (Int_t i=0; i<Runs; i++) {

		cout << "Finding track " << (i+1)*kEventsPerRun << " of " << Runs*kEventsPerRun << "... ";
		
		if (dataType == kMC) {
			t1.Start();

			eventID = di->getMCFrame(i, cf, x, y);
			di->getEventIDs(i, eventIDs);

			t1.Stop(); t2.Start();
			showDebug("Start diffuseFrame\n");

			cf->diffuseFrame(gRandom);

			showDebug("End diffuseFrame, start findHits\n");
			t2.Stop(); t3.Start();

			hits = cf->findHits(eventID);

			showDebug("Number of hits in frame: " << hits->GetEntriesFast() << endl);
			t3.Stop(); t4.Start();

			clusters = hits->findClustersFromHits(); // badly optimized
			clusters->removeSmallClusters(2);

			t4.Stop();

			clusters->matchWithEventIDs(eventIDs);
			eventIDs->Clear();
		}
		
		else if (dataType == kData) {
			di->getDataFrame(i, cf, energy);
			hits = cf->findHits();
			clusters = hits->findClustersFromHits();
			clusters->removeSmallClusters(2);
		}
		
		t5.Start();
		calorimeterTracks = clusters->findCalorimeterTracks();
		t5.Stop();
		
		if (calorimeterTracks->GetEntriesFast() == 0) breakSignal = kTRUE; // to stop running

		// Track improvements
		calorimeterTracks->extrapolateToLayer0();
		calorimeterTracks->splitSharedClusters();
		calorimeterTracks->removeTracksLeavingDetector();
		calorimeterTracks->removeTrackCollisions();
//		calorimeterTracks->retrogradeTrackImprovement(clusters);

		calorimeterTracks->Compress();
		calorimeterTracks->CompressClusters();
		
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
