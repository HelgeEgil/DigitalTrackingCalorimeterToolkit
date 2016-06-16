#include <iostream>

#include <TClonesArray.h>
#include <TObject.h>
#include <TCanvas.h>
#include <TView.h>
#include <TAttMarker.h>
#include <TAttLine.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>


#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

using namespace std;

Tracks::~Tracks() {
   tracks_.Delete();
   clustersWithoutTrack_.Delete();
}

void Tracks::Clear(Option_t *) {
	tracks_.Clear("C");
	clustersWithoutTrack_.Clear("C");
}

void Tracks::appendTrack(Track *copyTrack, Int_t startOffset /* default 0 */) {
   Int_t newIdx = tracks_.GetEntriesFast();
   Track *track = (Track*) tracks_.ConstructedAt(newIdx);

   for (Int_t i=0; i<copyTrack->GetEntriesFast(); i++) {
   	if(!copyTrack->At(i))
   		continue;

      track->appendCluster(copyTrack->At(i), startOffset);
   }
}

void Tracks::appendClustersWithoutTrack(TClonesArray *clustersWithoutTrack) {
	Int_t idxFrom = clustersWithoutTrack_.GetEntriesFast();
	Cluster *newCluster = 0;

	for (Int_t i=0; i<clustersWithoutTrack->GetEntriesFast(); i++) {
		if (!clustersWithoutTrack->At(i))
			continue;

		newCluster = (Cluster*) clustersWithoutTrack_.ConstructedAt(idxFrom + i);
		newCluster->set((Cluster*) clustersWithoutTrack->At(i));
	}
}

void Tracks::extrapolateToLayer0() {
   // do this for all tracks in collection
   for (Int_t i=0; i<tracks_.GetEntriesFast(); i++) {
   	if (!At(i)) continue;
      if (!At(i)->At(0)) At(i)->extrapolateToLayer0();
   }
}

void Tracks::splitSharedClusters() {
	// to be run as early as possible before joining Tracks* objects

	Cluster *interpolatedCluster;
	Clusters* interpolatedClusters = new Clusters();
	Clusters* interpolatedClosestClusters = new Clusters();

	Float_t dist, minDist = 1e5;
	Float_t x,y,sizeFactor, size;
	Int_t minIdx = 0, idx = 0, clusterIdx, nSplit = 0, nMissing = 0, nInterpolated = 0;
	Bool_t isMissing;

	for (Int_t layer=1; layer<nLayers; layer++) {
		
		vector<trackCluster> clustersThisLayer;
		vector<trackCluster> missingClustersThisLayer;
		
		clustersThisLayer.reserve(kEventsPerRun * 20);
		missingClustersThisLayer.reserve(kEventsPerRun * 5);

		for (Int_t i=0; i<GetEntriesFast(); i++) {
			Track *thisTrack = At(i);

			if (!thisTrack) continue;
			if (thisTrack->getLastLayer() < layer) continue;

			trackCluster thisCluster;
			thisCluster.track = i;

			if (thisTrack->hasLayer(layer)) {
				thisCluster.cluster = thisTrack->getClusterFromLayer(layer);
				clustersThisLayer.push_back(thisCluster);
			}

			else {
				missingClustersThisLayer.push_back(thisCluster);
				nMissing++;
			}
		}
		
		for (UInt_t i=0; i<missingClustersThisLayer.size(); i++) {
			Int_t trackIdx = missingClustersThisLayer.at(i).track;
			Track *missingTrack = At(trackIdx);
			if (!missingTrack ) continue;
		
			interpolatedCluster = missingTrack->getInterpolatedClusterAt(layer);
			if (interpolatedCluster) {
				minIdx = getClosestCluster(clustersThisLayer, interpolatedCluster);
				if (minIdx<0) continue; // no clusters this layer? What else could provoke this error?
				nInterpolated++;

				interpolatedClusters->appendCluster(interpolatedCluster);
				trackCluster closestTC= clustersThisLayer.at(minIdx); // THE CULPRIT
				Track *closestTrack = At(closestTC.track);
				Cluster *closestCluster = closestTrack->At(closestTC.cluster);
				
				minDist = diffmmXY(closestCluster, interpolatedCluster);
				Float_t clusterRadius = closestCluster->getRadiusmm();
				if (minDist < clusterRadius * 2) {
					
					Cluster *interpolatedClosestCluster = closestTrack->getInterpolatedClusterAt(layer);
					if (interpolatedClosestCluster) {
						interpolatedClosestClusters->appendCluster(interpolatedClosestCluster);
						
						// cluster size stays the same if minDist = 0, and increases to 2x when they the projected clusters barely touch (minDist > clusterRadius)
						sizeFactor = minDist / closestCluster->getRadiusmm() + 1;
						if (sizeFactor > 2) sizeFactor = 2;
						x = ( closestCluster->getX() + interpolatedCluster->getX() ) / 2;
						y = ( closestCluster->getY() + interpolatedCluster->getY() ) / 2;
						size = closestCluster->getSize() / sizeFactor;
						missingTrack->appendCluster(new Cluster(x, y, layer, size));
						nSplit++;
						
						// set size of closest cluster as well!
						closestCluster->setSize(size);

						sortTrackByLayer(trackIdx);
					}
				}
			}
		}
	}
	if (nMissing>0)
		cout << "From " << nMissing << " missing clusters, " << nInterpolated << " interpolated clusters were tested and " << nSplit << " clusters were split.\n";
}

Int_t Tracks::getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster) {
	Float_t minDist = 1e5, dist;
	Int_t minIdx = -1;

	for (UInt_t i=0; i<clusters.size(); i++) {
		trackCluster thisTrackCluster = clusters.at(i);
		Cluster *thisCluster = At(thisTrackCluster.track)->At(thisTrackCluster.cluster);
		
		dist = diffmmXY(thisCluster, interpolatedCluster);
		if (dist < minDist) {
			minDist = dist;
			minIdx = i;
		}
	}
	return minIdx;
}

void Tracks::sortTrackByLayer(Int_t trackIdx) {
	Int_t lastLayer = 0;
	Bool_t isSorted = true;
	Cluster *lastCluster;
	
	Int_t n = GetEntriesFast( trackIdx );
	
	Int_t firstLayer = 0;
	for (Int_t i=0; i<n; i++) {
		if (!At(i)) continue;
		firstLayer = i;
		break;
	}

	for (Int_t i=0; i<n; i++) {
		if (!At( trackIdx )->At(i)) continue;

		if (lastLayer > At( trackIdx )->getLayer(i))
			isSorted = false;

		lastLayer = At( trackIdx )->getLayer(i);
	}
	
	if (!isSorted) {
		vector<Int_t> sortList;
		sortList.resize(nLayers);
		for (Int_t i=0; i<nLayers; i++) sortList.at(i) = -1;
		for (Int_t i=0; i<n; i++) {
			if (!At(trackIdx)->At(i)) continue;
			sortList.at(At(trackIdx)->getLayer(i)) = i;
		}
		
		Track *newTrack = new Track();
		Cluster *nextCluster = 0;
		for (UInt_t i=0; i<sortList.size(); i++) {
			if (sortList.at(i) <0) continue;
			nextCluster = At(trackIdx)->At(sortList.at(i));
			newTrack->appendCluster(nextCluster, firstLayer);
		}

		if (newTrack->GetEntriesFast()) {
			appendTrack(newTrack);
			removeTrackAt(trackIdx);
		}
	}
}

void Tracks::checkLayerOrientation() {
	Float_t n[nLayers-1], x[nLayers-1], y[nLayers-1];
	Float_t dx, dy;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		Track *track = At(i);
		Int_t first = 0;
		if (!At(first)) first = 1;
		if (!At(first)) continue;
		Cluster *lastCluster = track->At(first);
		for (Int_t j=first+1; j<track->GetEntriesFast(); j++) {
			Cluster *thisCluster = track->At(j);
			dx = thisCluster->getX() - lastCluster->getX();
			dy = thisCluster->getY() - lastCluster->getY();

			x[j] += dx;
			y[j] += dy;
			n[j]++;

			lastCluster = thisCluster;
		}
	}

	cout << "The average track movement between layers is: ";
	for (Int_t i=0; i<nLayers-1; i++) {
		if (n[i])
			cout << Form("Layer %d: (%.1f, %.1f), ", i, x[i]/n[i], y[i]/n[i]);
	}
	cout << ".\n";
}

void Tracks::doFit() {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		At(i)->doFit();
	}
}

void Tracks::matchWithEventIDs(Hits * eventIDs) {
	// WRITE THIS TOMORROW
	// Check each point in each track against the hits in Hits * (and remove when found...)
	// N^2, but only done once and on MC
	// IF DIFFXY ( thisHit , trackHit ) on same layer is identical, assign eventID to cluster
	// is eventID same over whole track...? On which percent level?
	
	// ... and delete the sumEventID or something in Tools and in Hits.
	
	TStopwatch t1;
	Float_t minDist = 1e5; // px
	Float_t thisDist = 0;
	Track *thisTrack = nullptr;
	Cluster *thisCluster = nullptr;
	Hit *thisHit = nullptr;
	Int_t layer = -1;
	Int_t minIdx = 0;
	Bool_t doLoop = true;
	Int_t nHits = eventIDs->GetEntriesFast();

	Float_t cX, cY;

	TCanvas *c3 = new TCanvas("c3", "Distribution of closest Hit*", 1200, 800);
	TH1F *hist = new TH1F("hist", "Distribution of closest Hit*", 100, 0, 10);

	TCanvas *c4 = new TCanvas("c4", "HIT vs Cluster", 1200, 800);
	TH2F *hist2 = new TH2F("hist2", "HIT vs Cluster", 1280, 0, 1280, 1280, 0, 1280);

	hist->SetXTitle("Distance between cluster and hit (pixels)");
	hist->SetYTitle("Number of matchings");
	hist->SetFillColor(kBlue-7);
	hist->SetLineColor(kBlack);

	t1.Start();

	Int_t nClusters = 0;
	for (Int_t t=0; t<GetEntriesFast(); t++) {
		nClusters += At(t)->GetEntriesFast();
	}

	cout << "Number of clusters is " << nClusters << ", number of Hits is " << eventIDs->GetEntriesFast() << endl;


	for (Int_t t=0; t<GetEntriesFast(); t++) {
		thisTrack = At(t);

		for (Int_t c=0; c<thisTrack->GetEntriesFast(); c++) {
			thisCluster = thisTrack->At(c);
			layer = thisCluster->getLayer();

			cX = thisCluster->getX();
			cY = thisCluster->getY();

			minDist = diffXY(thisCluster, eventIDs->At(minIdx+1));

			if (minDist < 10) {
				minIdx++;
				doLoop = false;
			}

			else {
				doLoop = true;
				minDist = 1e5;
				minIdx = -1;
			}

			if (doLoop) {
				for (Int_t h=0; h<eventIDs->GetEntriesFast(); h++) {
					thisHit = eventIDs->At(h);
					
					if (!thisHit) continue;
					if (thisHit->getLayer() != layer) continue;

					if (fabs(cX - thisHit->getX()) < 10) {
						if (fabs(cY - thisHit->getY()) < 10) {

							thisDist = diffXY(thisCluster, thisHit);
							if (thisDist < minDist) {
								minDist = thisDist;
								minIdx = h;
							}
						}
					}
				}
			}
			
			if (minIdx >= 0 && minDist < 14) {
				thisCluster->setEventID(eventIDs->getEventID(minIdx));
				eventIDs->removeHitAt(minIdx);
				hist->Fill(minDist);
			}
		}
	}
	t1.Stop();

	c3->cd();
	hist->Draw();
	c3->SaveAs("OutputFiles/figures/EventIDMatching.png");

	cout << "Total time for function: " << t1.RealTime() << " seconds.\n";
	cout << "Number of Hits not matched (out of " << nHits << "): " << eventIDs->GetEntries() << endl;

	Int_t cWithoutEventID = 0;
	Bool_t printTrack = false;
	for (Int_t t=0; t<GetEntriesFast(); t++) {
		for (Int_t c=0; c<At(t)->GetEntriesFast(); c++) {
			if (At(t)->getEventID(c) < 0) {
				cWithoutEventID++;
				hist2->Fill(At(t)->getX(c), At(t)->getY(c), 2);
			}
		}
	}

	for (Int_t i=0; i<eventIDs->GetEntriesFast(); i++) {
		if (!eventIDs->At(i)) continue;
		hist2->Fill(eventIDs->getX(i), eventIDs->getY(i));
	}

	c4->cd();
	hist2->Draw("COLZ");
	c4->SaveAs("OutputFiles/figures/EventIDMatching2D.png");
	cout << "Number of clusters without eventID: " << cWithoutEventID << "( " << (float) cWithoutEventID / nClusters * 100 << "%)\n";
}
