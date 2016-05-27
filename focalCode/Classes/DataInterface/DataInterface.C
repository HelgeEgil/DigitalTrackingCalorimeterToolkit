#define DataInterface_cxx

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>

#include <TEllipse.h>
#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TView.h>
#include <TLeaf.h>

#include "GlobalConstants/Constants.h"
#include "Classes/Hit/Hit.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Track/conversionFunctions.h"
#include "HelperFunctions/Tools.h"
#include "Classes/Layer/Layer.h"
#include "Classes/DataInterface/DataInterface.h"
#include "Classes/CalorimeterFrame/CalorimeterFrame.h"
#include "RootFiles/LinkDef.h"

void DataInterface::getMCData(Int_t runNo, TH3F* Frame3D) {
	if (fChain==0) return;

   Int_t eventIdFrom = runNo * kEventsPerRun/10;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun/10;

   Float_t offsetX = (nx+2) * dx;
   Float_t offsetY = (ny) * dy;
   Float_t x,y;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nb = 0;

	for (Long64_t jentry=0; jentry<nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		nb = fChain->GetEntry(jentry);

		if (eventID < eventIdFrom) continue;
		else if (eventID >= eventIdTo) break;

      x = (posX + offsetX) * nx / (offsetX);
      y = (posY + offsetY) * ny / (offsetY);

		Frame3D->Fill(posZ, x, y, edep*1000);
	}
	
	if (kDebug) cout << "From runNo " << runNo << ", nentries = " << nentries << endl;
	
} // end function GetData3D

void DataInterface::getDataFrame(Int_t runNo, CalorimeterFrame * cf, Int_t energy) {

	if (!existsEnergyFile(energy)) {
		cout << "There are no data files with energy " << energy << endl;
		return;
	}

   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

	TString fn = Form("../../Data/ExperimentalData/DataFrame_%i_MeV.root", energy);
	TFile *f = new TFile(fn);
	TTree *tree = (TTree*) f->Get("tree");

	Int_t nentries = tree->GetEntries();
	if (eventIdTo > nentries) {
		eventIdTo = nentries;
	}
	cout << "Found " << nentries << " frames in the DataFrame.\n";

	TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
	TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
	TLeaf *lLayer = tree->GetLeaf("fDataFrame.fLayer");

	Int_t counter = 0;
	for (Int_t i=eventIdFrom; i<eventIdTo; i++) {
		tree->GetEntry(i);

		for (Int_t j=0; j<lX->GetLen(); j++) {
			Int_t x = lX->GetValue(j) + nx;
			Int_t y = lY->GetValue(j) + ny;
			Int_t z = lLayer->GetValue(j);

			cf->fillAt(z, x, y);
		}
		counter++;
	}
	delete f;
}

Int_t DataInterface::getMCFrame(Int_t runNo, CalorimeterFrame *cf, Float_t *x_energy, Float_t *y_energy) {
	// Retrieve kEventsPerRun events and put them into CalorimeterFrame*
	
	Int_t eventIdFrom = runNo * kEventsPerRun;
	Int_t eventIdTo = eventIdFrom + kEventsPerRun;

	if (runNo == 0) lastJentry_ = 0;

	Float_t offsetX = (nx+2) * dx;
	Float_t offsetY = (ny) * dy;
	Float_t x,y;
	Float_t sum_edep = 0;
	Int_t n = 0;
	Int_t calorimeterLayer = 0;
	Int_t blackListEventID = -1;
	Int_t whiteListEventID = -1;
	Int_t lastID = 0;

	if (fChain == 0) return -1;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	for (Long64_t jentry=lastJentry_; jentry<nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry<0) break;
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		
		if (lastID != eventID) {
			sum_edep = 0;
			n = 0;
		}
		
		sum_edep += edep;
		
		if  (x_energy && y_energy && parentID == 0) {
			x_energy[n + sizeOfEventID*eventID] = posZ + 1.5; // 1.5 for Al absorber
			y_energy[n + sizeOfEventID*eventID] = run_energy - sum_edep;
			n++;
		}
		
		if (volumeID[4] == 4 || !x_energy) { // hit in chip, if no x_energy provided do this
			
			if (eventID < eventIdFrom) {
				continue;
			}
			
			else if (eventID >= eventIdTo) {
				lastJentry_ = jentry;
				break;
			}
			
			if (eventID == blackListEventID) continue;

			if (posZ < -100) {
				blackListEventID = eventID;
				whiteListEventID = eventID;
				continue;
			}
			
			calorimeterLayer = level1ID;

			if (calorimeterLayer<0 || posZ < -20) {
				continue;
			}

			x = (posX + offsetX) * nx / (offsetX);
			y = (posY + offsetY) * ny / (offsetY);

			cf->fillAt(calorimeterLayer, x, y, edep*1000);
		}
		
		lastID = eventID;
	}

	if (kEventsPerRun > 1) return -1;
	return eventID;
}
