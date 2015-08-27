#define Focal_cxx
#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
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

// mine
#include "Constants.h"
#include "TFocal.h"
#include "Hit.h"
#include "Hits.h"
#include "Cluster.h"
#include "Clusters.h"
#include "Track.h"
#include "Tracks.h"
#include "Layer.h"
#include "CalorimeterFrame.h"
#include "TrackerFrame.h"
#include "LinkDef.h"

//using namespace std;

/*
void Focal::GetData3D(Int_t RunFrom, Int_t RunTo, TH3F* Frame3D) {
	if (fChain==0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nb = 0;

	for (Long64_t jentry=0; jentry<nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
		nb = fChain->GetEntry(jentry);

		if (eventID < RunFrom) continue;
		if (eventID > RunTo) break;

		Frame3D->Fill(posZ, posX, posY, edep*1000);
	}
} // end function GetData3D

void Focal::GetTrackerData(Int_t EventFrom, Int_t EventTo, vector<TH2F*> *Frame3D) {
   // contains a list of clusters

   if (fChain == 0) return;
   TRandom3 *gRandom = new TRandom3(0);

   Long64_t nentries = fChain->GetEntriesFast();

   Frame3D->reserve((int) nTrackers);
   for (Int_t i=0; i<nTrackers; i++) {
      Frame3D->push_back(new TH2F(Form("Frame3D.%i", i), 
                           Form("Energy deposition in tracker %i", i),
                           nx*2, 0, nx*2, ny*2, 0, ny*2));
   } // end fill vector with TH2F layers

   Long64_t nbytes = 0, nb = 0;
	Float_t X = 0, Y = 0;
	Int_t Tracker = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) { // loop over entries in tree
		Long64_t ientry = LoadTree(jentry);
      if (ientry<0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

		if (eventID <= EventFrom) continue;
		if (eventID > EventTo) break;
      if (level1ID > 3) continue; // regular layers

      Tracker = level1ID - 4;
      
		// with integration of four extra trackers
      // -4, -3, -2, -1 are the trackers
      
      X = (posX + 19.25) * nx / (19.25);
      Y = (posY + 19.2) * ny / (19.2);
     
      DiffuseAndFill(Frame3D, gRandom, X, Y, Tracker, edep*1000);
   } // end loop over entries
} // end function GetTrackerData

void Focal::GetRealData3D(Int_t RunFrom, Int_t RunTo, TH3F* Frame3D) {
	// Get data from DataFrame root files
	// TODO: Write a function to do data reduction (maybe later when doing clustering...?)
	TFile *f = new TFile("DataFrame/DataFrame_170_MeV.root");
	TTree *tree = (TTree*) f->Get("tree");

	Int_t nentries = tree->GetEntries();

	// Read information about
	// X position (-644 -> 644)
	// Y position (-644 -> 644)
	// layer (0 -> 24)

	TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
	TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
	TLeaf *lLayer = tree->GetLeaf("fDataFrame.fLayer");

	Int_t counter = 0;
	for (Int_t i=0; i<nentries; i++) {
		if (counter < RunFrom) continue;
		if (counter > RunTo) break;
		tree->GetEntry(i);

		Int_t nleafs = lX->GetLen();
		cout << " ------ " << nleafs << " leafs ------ \n";
		for (Int_t j=0; j<nleafs; j++) {
			Int_t x = lX->GetValue(j);
			Int_t y = lY->GetValue(j);
			Int_t z = lLayer->GetValue(j);

			cout << "(" << x << "," << y << "," << z << ")\n";

			Frame3D->Fill(z, x, y);
			counter++;
		}
	}
}

void Focal::GetClusterFrames(vector<Hits*> *clusterHitMap, Bool_t MC) {
	// return TH2F's with different cluster shapes
	
	// Initialize vector layerHits
	vector<Hits*> *layerHits = new vector<Hits*>;
	layerHits->reserve((int) nLayers);
	for (Int_t layer=0; layer<nLayers; layer++) {
		layerHits->push_back(new Hits(500));
	}

	if (!MC) {
		TFile *f = new TFile("DataFrame/DataFrame_190_MeV.root");
		TTree *tree = (TTree*) f->Get("tree");

		Int_t nFrames = tree->GetEntries();

		cout << "There are " << nFrames << " spills in Run.\n";

		TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
		TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
		TLeaf *lZ = tree->GetLeaf("fDataFrame.fLayer");

		Int_t nLeafs;


		// Loop over frames
		for (Int_t f=0; f<nFrames; f++) {
			tree->GetEntry(f);
			nLeafs = lX->GetLen();

			// Fill *hits and find layers for this Frame
			for (Int_t h=0; h<nLeafs; h++) {
				Int_t x = lX->GetValue(h) + nx;
				Int_t y = lY->GetValue(h) + ny;
				Int_t z = lZ->GetValue(h);

				layerHits->at(z)->appendPoint(x, y, z, f);
			}

			for (Int_t layer=0; layer<nLayers; layer++) {
				FindClustersHitMap(layerHits->at(layer), clusterHitMap);
			}

			if (clusterHitMap->size() > 1000) break;
		} // end loop over frame
	}

	else { // Monte Carlo
	
		Int_t f=0;

		if (fChain == 0) return;

		Long64_t nentries = fChain->GetEntriesFast();
		Long64_t nbytes = 0, nb = 0;

		for (Long64_t jentry=0; jentry<nentries; jentry++) {
			Long64_t ientry = LoadTree(jentry);
			if (ientry<0) break;
			nb = fChain->GetEntry(jentry); nbytes += nb;

			Float_t X = (posX + 19.25) * nx / (19.25);
			Float_t Y = (posY + 19.2) * ny / (19.2);
			Float_t Z = level1ID;

			layerHits->at(Z)->appendPOint(X,Y,Z,f++);
		}

		// FIXME Do this based on eventID;

		for (Int_t layer=0; layer<nLayers; layer++) {
			FindClustersHitMap(layerHits->at(layer), clusterHitMap);
		}
//		if (clusterHitMap->size() > 1000) break;
	}

}

void Focal::GetRealFrame2D(Int_t Runs, Int_t Layer, TH2F *Frame2D) {

	TFile *f = new TFile("DataFrame/DataFrame_190_MeV.root");
	TTree *tree = (TTree*) f->Get("tree");

	Int_t nentries = tree->GetEntries();

	TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
	TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
	TLeaf *lLayer = tree->GetLeaf("fDataFrame.fLayer");

	Double_t x;
	Double_t y;

	for (Int_t i=0; i<nentries; i++) {
		tree->GetEntry(i);

		Int_t nleafs = lX->GetLen();
		for (Int_t j=0; j<nleafs; j++) {
			if (i>=Runs) break;
			if (lLayer->GetValue(j) == Layer) {
				x = lX->GetValue(j) + nx;
				y = lY->GetValue(j) + ny;
				Frame2D->Fill(x, y);
			}
		}
   }
}

void Focal::GetFrame2D(Int_t Runs, Int_t Layer, TH2F *Frame2D) {
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   // Until layer geometry is flipped
//   Layer = (23 - Layer);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { // jentry<nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
     if (eventID > Runs) break;

      // calculate X and Y
     // Layer < 0 projects all layers
   
   if ((Layer > 0 && level1ID == Layer) || (Layer < 0)) {

      // ok, let's try to draw without 2ndies
      if (parentID==0) continue;
      
      Double_t X = (posX + 19.25) * nx / (19.25);
      Double_t Y = (posY + 19.2) * ny / (19.2);

      Frame2D->Fill(X, Y, edep*1000);
      }
   }
}
*/

/*
// TOOL
TrackerCollection* Focal::FindTrackerFrame(Int_t EventFrom, Int_t EventTo, Clusters *restPoints) {
	Int_t Runs = EventTo - EventFrom;	

	cout << "Before GetTrackerData\n";
   vector<TH2F*> *Frame3D = new vector<TH2F*>;
   GetTrackerData(EventFrom, EventTo, Frame3D);
	cout << "After GetTrackerData\n";

   // working cluster hit list
   // to be emptied throughout run
	Clusters *preClusters = new Clusters(Runs*ntrackers*10); // layer 0-1
   Clusters *postClusters = new Clusters(Runs*nTrackers*10); // layer 2-3
   TrackerCollection *trackerFrames = new TrackerCollection(Runs*nTrackers*20);
   TrackerCollection *preTrackerFrames = new TrackerCollection(Runs*nTrackers*10);
   TrackerCollection *postTrackerFrames = new TrackerCollection(Runs*nTrackers*10);
   Cluster *nullCluster = new Cluster();

	Tracker *tracker = new Tracker();
   Tracker *preTracker = new Tracker();
   Tracker *postTracker = new Tracker();

   // Loop over layers to find connected clusters
   // Put them in vector clusters with information about

   for (Int_t l=-4; l<-4+ntrackers; l++) {
      
      Hits *hits = new Hits(Runs*200);
      Clusters *clustersThisLayer = new Clusters(Runs*20);
		
		cout << "Before FindHits\n";
      FindHits(Frame3D->at(l + nTrackers), hits, l);
		cout << "After FindHits, before FindClustersFromHits\n";
      FindClustersFromHits(hits, clustersThisLayer, Frame3D); 
		cout << "After FindClustersFromHits\n";

      Int_t nClusters = 0;

      for (Int_t i=0; i<clustersThisLayer->GetEntriesFast(); i++) {
         if (l==-4 || l==-3) preClusters->appendCluster(clustersThisLayer->At(i));
         if (l==-2 || l==-1) postClusters->appendCluster(clustersThisLayer->At(i));
         nClusters++;
      } // end loop over clusters for layer l

      delete clustersThisLayer;
      delete hits;
   } // end loop over all layers

	for (UInt_t f=0; f<Frame3D->size(); f++) {
		delete Frame3D->at(f);
	}

	delete Frame3D;

   cout << "Found " << preClusters->GetEntriesFast() << " clusters in pre-tracker frames.\n";
   cout << "Found " << postClusters->GetEntriesFast() << " clusters in post-tracker frames.\n";

	cout << "Before optimization\n";
   preClusters->makeLayerIndex(); // optimization step
   postClusters->makeLayerIndex(); // optimization step
	cout << "After optimization\n";


   Clusters *seeds = new Clusters(1000);
   Track* bestTrack = new Track();
  
	cout << "Before FindSeeds\n";	
   FindSeeds(preClusters, seeds, -4);
	cout << "After FindSeeds\n";

   cout << "Found " << seeds->GetEntriesFast() << " seeds in pre-tracker layer.\n";

   Int_t searchToLayer = -3;

   for (Int_t i=0; i<seeds->GetEntriesFast(); i++) { // loop through layer -4 seeds
		if (!seeds->At(i)) continue;
		cout << "pre followseeds\n";
		FollowSeeds(preClusters, seeds->At(i), bestTrack, searchToLayer);
		cout << "post followseeds\n";

      if (bestTrack->GetEntriesFast() > 0) {
         preTracker->SetPreTracker1(bestTrack->At(0));
         preTracker->SetPreTracker2(bestTrack->At(1));
         
         preTrackerFrames->AddTracker(preTracker);
         preClusters->removeAllClustersInTrack(bestTrack);
      } // end if bestTrack
      bestTrack->Clear();
   } // end loop through layer 0 seeds

   seeds->Clear();
   FindSeeds(postClusters, seeds, -2);

   searchToLayer = -1;
   
   cout << "Found " << seeds->GetEntriesFast() << " seeds in post-tracker layer.\n";

   for (Int_t i=0; i<seeds->GetEntriesFast(); i++) { // loop through layer -2 seeds
		if (!seeds->At(i)) continue;
      FollowSeeds(postClusters, seeds->At(i), bestTrack, searchToLayer);

      if (bestTrack->GetEntriesFast() > 0) {
         // found four here...?
			cout << "Found " << bestTrack->GetEntriesFast() << " entries in bestTrack for post-tracker.\n";
			for (Int_t j=0; j<bestTrack->GetEntriesFast(); j++) {
				cout << bestTrack->At(j) << endl;
				if (bestTrack->At(j)) cout << *bestTrack->At(j) << endl;
			}
         
         postTracker->SetPostTracker1(bestTrack->At(0));
         postTracker->SetPostTracker2(bestTrack->At(1));
         
         postTrackerFrames->AddTracker(postTracker);
         postClusters->removeAllClustersInTrack(bestTrack);
      } // end if bestTrack
		bestTrack->Clear();
   } // end loop through layer 1 seeds


	// Try to connect tracks from the pre-tracker to the post-tracker...
	// In this part the event count should be LOW to get a good match!
	
	Int_t nPre = preTrackerFrames->GetEntriesFast();
	Int_t nPost = postTrackerFrames->GetEntriesFast();

	for (Int_t i=0; i<nPre; i++) {
		// 1) Find entry angle & position
		// 2) Find projected (straight line) path to hit on post-tracker layer
		// 3) For each post-track hit, calculate distance between expected hit and hit
		// 4) If the smallest distance is < delta, connect the two events:
		//		4a) Add them to trackerFrames object
		//		4b) Remove them from PreTrackerFrame and PostTrackerFrame

		if (!preTrackerFrames->At(i)) continue;


		Cluster preDerivative = preTrackerFrames->At(i)->GetPreDerivative();
		
		Cluster *prePosition = preTrackerFrames->At(i)->GetPreTracker2();
		Cluster *postPosition = new Cluster(0, 0, -2);
		
		cout << "Treating tracker hit at " << *prePosition << endl;

		Float_t delta_z = postPosition->getLayermm() - prePosition->getLayermm();

		postPosition->setXmm(prePosition->getXmm() + delta_z * preDerivative.getX());
		postPosition->setYmm(prePosition->getYmm() + delta_z * preDerivative.getY());

		cout << "Searching for layer in " << *postPosition << "... ";

		Int_t minidx = -1;
		Float_t mindelta = 1e10;

		for (Int_t j=0; j<nPost; j++) {
			if (!postTrackerFrames->At(j)) continue;

			Float_t delta = sqrt(pow(postPosition->getXmm() - postTrackerFrames->At(j)->GetPostPositionXmm(), 2) +
										pow(postPosition->getYmm() - postTrackerFrames->At(j)->GetPostPositionYmm(), 2));

			if (delta < mindelta) {
				// new candidate for tracker connection
				mindelta = delta;
				minidx = j;
			}
		}

		cout << " ... done. The closest connection at actual point " << postTrackerFrames->At(minidx) << " with mindelta = " << mindelta << endl;
	
		tracker->SetPreTracker1(preTrackerFrames->At(i)->GetPreTracker1());
		tracker->SetPreTracker2(preTrackerFrames->At(i)->GetPreTracker2());
		tracker->SetPostTracker1(postTrackerFrames->At(minidx)->GetPostTracker1());
		tracker->SetPostTracker2(postTrackerFrames->At(minidx)->GetPostTracker2());

		trackerFrames->AddTracker(tracker);
	}

	cout << "After FindClustersFromHits\n";
	cout << "There are now " << trackerFrames->GetEntriesFast() << " entries in trackerFrames.\n";

	Int_t clustersLeft = 0;
	for (Int_t i=0; i<preClusters->GetEntriesFast(); i++) {
		if (preClusters->At(i)) {
		   clustersLeft++;
		   if (restPoints)  restPoints->appendCluster(postClusters->At(i));
      }
	}
	for (Int_t i=0; i<postClusters->GetEntriesFast(); i++) {
		if (postClusters->At(i)) {
		   clustersLeft++;
		   if (restPoints)  restPoints->appendCluster(postClusters->At(i));
      }
	}

   cout << "There are " << clustersLeft << " of " << preClusters->GetEntriesFast() + postClusters->GetEntriesFast() << " remaining hits not assigned to track!\n";

   delete seeds;
   delete bestTrack;
   delete preClusters;
   delete postClusters;
	delete tracker;
	delete preTracker;
	delete postTracker;
	delete nullCluster;
	delete preClusters;
	delete postClusters;
	delete preTrackerFrames;
	delete postTrackerFrames;
   
   return trackerFrames;
   
} // end function FindTracks

// TOOL

// TOOL
Tracks* Focal::FindRealTracks(Int_t Runs, Clusters *restPoints, Int_t energy) {
	// energy is default = -1

	// number of objects is total number of frames
	// Frames in spill O(1000) * # spills O(50) * # Runs

	TFile *f;

	if (energy < 0) f = new TFile("DataFrame/DataFrame_190_MeV.root");
	else				 f = new TFile(Form("DataFrame/DataFrame_%i_MeV.root", energy));

	TTree *tree = (TTree*) f->Get("tree");

	Int_t nFrames = tree->GetEntries();

	if (nFrames > Runs) {
		nFrames = Runs;
	}

	cout << "There are " << nFrames << " spills in Run.\n";

	TLeaf *lX = tree->GetLeaf("fDataFrame.fX");
	TLeaf *lY = tree->GetLeaf("fDataFrame.fY");
	TLeaf *lZ = tree->GetLeaf("fDataFrame.fLayer");

	Int_t nLeafs;

	// Make some containers
   Tracks *tracks = new Tracks(1000);
	vector<Hits*> *layerHits = new vector<Hits*>;
	layerHits->reserve(nLayers);

//	HitCollection *hits = new HitCollection(2000);
   Clusters *clustersThisLayer = new Clusters(500);
	TClonesArray *Frames = new TClonesArray("ClusterCollection", nFrames);
	Frames->SetOwner(true);

	// Initialize vector
	for (Int_t layer=0; layer<nLayers; layer++) {
		layerHits->push_back(new Hits(500));
	}

	// Loop over frames
	for (Int_t f=0; f<nFrames; f++) {
//		cout << "Finding hits + clusters for frame number " << f << endl;

		// make new ClusterCollection* entry in Frames
		Int_t next = Frames->GetEntriesFast();
		Clusters *c = (Clusters*) Frames->ConstructedAt(next);

		// Get data
		tree->GetEntry(f);
		nLeafs = lX->GetLen();

//		cout << "Found " << nLeafs << " hits in Frame.\n";

		// Fill *hits and find layers for this Frame
		for (Int_t h=0; h<nLeafs; h++) {
			Int_t x = lX->GetValue(h) + nx;
			Int_t y = lY->GetValue(h) + ny;
			Int_t z = lZ->GetValue(h);

			layerHits->at(z)->appendPOint(x, y, z, f);
		}

		for (Int_t layer=0; layer<nlayers; layer++) {
			FindClustersFromHits(layerHits->at(layer), clustersThisLayer);
		}

//		cout << "Number of clusters in frame before removal: " << clustersThisLayer->GetEntries();
		RemoveSmallClusters(clustersThisLayer, 2); // 2nd arg is remove clusters up to this size
//		cout << ", and after: " << clustersThisLayer->GetEntries() << endl;

		Int_t nClusters = 0;
		for (Int_t i=0; i<clustersThisLayer->GetEntriesFast(); i++) {
			if (!clustersThisLayer->At(i)) continue;
			c->appendCluster(clustersThisLayer->At(i));
			nClusters++;
		} // end loop over all clusters in layer

		for (Int_t layer=0; layer<nLayers; layer++) layerHits->at(layer)->Clear();
		clustersThisLayer->Clear();
	} // end loop over frame

	while (!layerHits->empty()) {
		delete layerHits->back();
		layerHits->pop_back();
	}

	delete layerHits;
	delete clustersThisLayer;

	// TODO:
	// Here we can introduce a test:
	// If there are very few hits in clustersThisLayer and
	// they are not succesive, continue

	for (Int_t i=0; i<nFrames; i++) { // loop over frames

		Clusters *clusters = ((Clusters*) Frames->At(i));

		clusters->makeLayerIndex(); // optimization step
		Clusters *seeds = new Clusters(1000);
		Track *bestTrack = new Track();

		FindSeeds(clusters, seeds, 0);

		for (Int_t i=0; i<seeds->GetEntriesFast(); i++) { // loop through layer 0 seeds
			if (!seeds->At(i)) continue;
			FollowSeeds(clusters, seeds->At(i), bestTrack);

			if (bestTrack->GetEntriesFast() > 0) {
				tracks->appendTrack(bestTrack); // copy objects
				clusters->removeAllClustersInTrack(bestTrack);
			} // end if bestTrack
			bestTrack->Clear();
		} // end loop through layer 0 seeds

		seeds->Clear();
		FindSeeds(clusters, seeds, 1);

		for (Int_t i=0; i<seeds->GetEntriesFast(); i++) { // loop through layer 1 seeds
			if (!seeds->At(i)) continue;
			FollowSeeds(clusters, seeds->At(i), bestTrack);

			if (bestTrack->GetEntriesFast() > 0) {
				tracks->appendTrack(bestTrack);
				clusters->removeAllClustersInTrack(bestTrack);
			} // end if bestTrack
			bestTrack->Clear();
		} // end loop through layer 1 seeds

//		cout << "Found " << tracks->GetEntriesFast() << " tracks from " << Runs << " runs.\n";
		Int_t clustersLeft = 0;
		for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
			if (clusters->At(i)) {
			   clustersLeft++;
            // add rest cluster to restPoints object
            restPoints->appendCluster(clusters->At(i));
         }
		}
//		cout << "There are " << clustersLeft << " of " << clusters->GetEntriesFast() 
//			  << " remaining hits not assigned to track!\n";
		delete seeds;
		delete bestTrack;
	} // end loop over frames

	delete Frames;

   cout << endl;
   return tracks;
} // end function FindTracks

*/

// TOOL

void Focal::getFrame3D(Int_t runNo, CalorimeterFrame *cf) {
	// Retrieve kEventsPerRun events and put them into CalorimeterFrame*
	
   Int_t eventIdFrom = runNo * kEventsPerRun;
   Int_t eventIdTo = eventIdFrom + kEventsPerRun;

   Float_t offsetX = (nx+2) * dx;
   Float_t offsetY = (ny) * dy;
   Float_t x,y;
   Int_t calorimeterLayer;

	if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
		Long64_t ientry = LoadTree(jentry);
      if (ientry<0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;

		if (eventID < eventIdFrom) continue;
		if (eventID >= eventIdTo) break;

      calorimeterLayer = level1ID - 4;
      if (calorimeterLayer<0) continue;

      x = (posX + offsetX) * nx / (offsetX);
      y = (posY + offsetY) * ny / (offsetY);

	   cf->fillAt(calorimeterLayer, x, y, edep*1000);
   }
}
