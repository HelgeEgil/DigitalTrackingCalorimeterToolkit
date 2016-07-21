#include <iostream>
#include <vector>

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
						missingTrack->appendCluster(new Cluster(x, y, layer, size, -1));
						nSplit++;
						
						// set size of closest cluster as well!
						closestCluster->setSize(size);
						closestCluster->setEventID(-1);

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

	Int_t nClusters = 0;
	for (Int_t t=0; t<GetEntriesFast(); t++) {
		nClusters += At(t)->GetEntriesFast();
	}

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
	c3->cd();
	hist->Draw();
	c3->SaveAs("OutputFiles/figures/EventIDMatching.png");

//	cout << "Number of Hits not matched (out of " << nHits << "): " << eventIDs->GetEntries() << endl;

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

Int_t Tracks::getTrackIdxFromFirstLayerEID(Int_t eid) {
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		if (!At(i)->At(0)) continue;
		if (eid == At(i)->At(0)->getEventID()) return i;
	}

	return -1;
}

void Tracks::removeTracksLeavingDetector() {
	Int_t nTracksRemoved = 0;
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		if (!At(i)->Last()) continue;

		Int_t lastLayer = At(i)->Last()->getLayer();
		Cluster *nextPoint = getTrackExtrapolationToLayer(At(i), lastLayer + 1);
		
		if (isPointOutOfBounds(nextPoint, -15)) {
			removeTrack(At(i));
			nTracksRemoved++;
		}
	}
	cout << "Tracks::removeTracksLeavingDetector() has removed " << nTracksRemoved  << " tracks.\n";
}

vector<Int_t> * Tracks::getTracksWithConflictClusters() {
	const Int_t n = GetEntriesFast();
	Track * thisTrack = nullptr;

	vector<Int_t> * trackIdx = new vector<Int_t>;
	trackIdx->reserve(n/10.);

	for (Int_t i=0; i<n; i++) {
		thisTrack = At(i);
		if (!thisTrack) continue;
		if (thisTrack->isUsedClustersInTrack()) { trackIdx->push_back(i); }
	}

	return trackIdx;
}

vector<Int_t> * Tracks::getConflictingTracksFromTrack(Int_t trackIdx) {
	Cluster * conflictCluster = nullptr;
	vector<Int_t> * conflictingTracks = new vector<Int_t>;
	conflictingTracks->reserve(5);

	Track * trackA = At(trackIdx);

	if (!trackA) { return conflictingTracks; }

	Clusters * conflictClusters = trackA->getConflictClusters();
	Int_t nClusters = conflictClusters->GetEntriesFast();
	Int_t nClustersReal = conflictClusters->GetEntries();

	for (Int_t i=0; i<nClusters; i++) {
		conflictCluster = conflictClusters->At(i);
		if (!conflictCluster) continue;

		vector<Int_t> * possibleTracks = getTracksFromCluster(conflictCluster);

		for (UInt_t j=0; j<possibleTracks->size(); j++) {
			if (!isItemInVector(possibleTracks->at(j), conflictingTracks)) {
				conflictingTracks->push_back(possibleTracks->at(j));
			}
		}
	}

	return conflictingTracks;
}

Int_t Tracks::getTrackIdxFromCluster(Cluster * cluster) {
	Track * thisTrack = nullptr;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		thisTrack = At(i);
		if (!thisTrack) continue;

		if (thisTrack->isClusterInTrack(cluster)) {
			return i;
		}
	}
}

vector<Int_t> * Tracks::getTracksFromCluster(Cluster * cluster) {
	Track * thisTrack = nullptr;
	vector<Int_t> * tracksWithCluster = new vector<Int_t>;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		thisTrack = At(i);
		if (!thisTrack) continue;

		if (thisTrack->isClusterInTrack(cluster)) {
		   Int_t idx = thisTrack->getClusterIdx(cluster);
		   thisTrack->At(idx)->markUsed();
			tracksWithCluster->push_back(i);
		}
	}

	return tracksWithCluster;
}	

void Tracks::removeTrackCollisions() {
	// sometimes two tracks end in the same cluster
	// remove the last cluster from the smallest track IF THE SMALLEST TRACK IS <= 2.. delete the track
	// if the final length is 1
	
	vector<Int_t> * conflictPair = nullptr;
	Clusters * conflictClusters[2] = {};
	Track * track[2] = {};
	Int_t len[2] = {0};
	Int_t longTrack, shortTrack;

	Int_t nRemoved = 0;
	Int_t nShortRemoved = 0;

	vector<Int_t> * conflictTracks = getTracksWithConflictClusters();

	for (UInt_t i=0; i<conflictTracks->size(); i++) {
		conflictPair = getConflictingTracksFromTrack(conflictTracks->at(i)); // 2 tracks
      if (conflictPair->size() < 2) continue;

		for (Int_t j=0; j<2; j++) {
			track[j] = At(conflictPair->at(j));
			conflictClusters[j] = track[j]->getConflictClusters(); // 1 cluster
			cout << "Added j=" << j << " conflict cluster \033[1m " << conflictClusters[j]->At(0) << "\033[0m ";
			if (conflictClusters[j]->At(0)) cout << "-> \033[1m" << *conflictClusters[j]->At(0) << "\033[0m.\n";
			else cout << endl;
			len[j] = track[j]->GetEntriesFast();
		}
		
		longTrack = (len[0] > len[1]) ? 0 : 1;
		shortTrack = 1 - longTrack;

		Bool_t sameCluster = isSameCluster(conflictClusters[0]->At(0), conflictClusters[1]->At(0));
   	if (!sameCluster) { continue; }
		
		Cluster * conflictCluster = conflictClusters[0]->At(0); // only one conflicting cluster!

      if (!conflictCluster) continue;

		Int_t longClusterIdx = track[longTrack]->getClusterIdx(conflictCluster);
		Int_t shortClusterIdx = track[shortTrack]->getClusterIdx(conflictCluster);

      if (shortClusterIdx != track[shortTrack]->GetEntriesFast() - 1) {
         // not a terminating cluster!
         continue;
      }

		track[longTrack]->At(longClusterIdx)->markUnused();
	
		if (len[shortTrack] > 2) {
			track[shortTrack]->removeClusterAt(shortClusterIdx);
		}
		
		else {
			removeTrack(track[shortTrack]);
			nShortRemoved++;
		}
		nRemoved++;
	}

	cout << "Tracks::removeTrackCollisions removed " << nRemoved << " collisions, and removed " << nShortRemoved << " too short tracks.\n";
}

void Tracks::retrogradeTrackImprovement(Clusters * clusters) {
	Track * thisTrack = nullptr;
	Cluster * estimatedRetrogradeCluster = nullptr;
	Cluster * nearestNeighborAllTracks = nullptr;
	Cluster * actualCluster = nullptr;
	Float_t deltaThisTrack = 0;
	Float_t deltaBestTrack = 0;
	Bool_t anotherTrackIsCloserMatch;

	Track * switchTrack = nullptr;
	Int_t switchTrackIdx = -1;
	Int_t switchTrackAtLayer = -1;
	Int_t idxThis, idxSwitch;

	Tracks * newTracks = new Tracks();

	Bool_t ignoreTrack[GetEntriesFast()];
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		ignoreTrack[i] = false;
	}

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (ignoreTrack[i]) continue;

		switchTrack = nullptr;
		switchTrackIdx = -1;
		switchTrackAtLayer = -1;

		thisTrack = At(i);
		if (!thisTrack) continue;

		Int_t lastIdx = thisTrack->GetEntriesFast() - 1;
		Int_t lastLayer = thisTrack->getLayer(lastIdx);
		
		for (Int_t layer = lastLayer - 2; layer>=0; layer--) {
			if (thisTrack->getIdxFromLayer(layer) < 0) continue;

			estimatedRetrogradeCluster = getRetrogradeTrackExtrapolationToLayer(thisTrack, layer);
			if (!estimatedRetrogradeCluster) continue;
			
			actualCluster = thisTrack->At(thisTrack->getClusterFromLayer(layer));
			nearestNeighborAllTracks = clusters->findNearestNeighbour(estimatedRetrogradeCluster, false);	

			if (!nearestNeighborAllTracks) {
				delete estimatedRetrogradeCluster;
				continue;
			}

			deltaThisTrack = diffXY(estimatedRetrogradeCluster, actualCluster);
			deltaBestTrack = diffXY(estimatedRetrogradeCluster, nearestNeighborAllTracks);
			delete estimatedRetrogradeCluster;
			
			anotherTrackIsCloserMatch = (deltaBestTrack < deltaThisTrack);

			if (anotherTrackIsCloserMatch) {
				switchTrackIdx = getTrackIdxFromCluster(nearestNeighborAllTracks);
				switchTrack = At(switchTrackIdx);
				switchTrackAtLayer = layer;
				break;
			}
		}

		if (switchTrack) {
			Track * newTrackA = new Track();
			Track * newTrackB = new Track();
			
			Int_t longestTrack = max(thisTrack->Last()->getLayer(), switchTrack->Last()->getLayer());
			
			for (Int_t layer=0; layer<=longestTrack; layer++) {
				idxThis = thisTrack->getIdxFromLayer(layer);
				idxSwitch = switchTrack->getIdxFromLayer(layer);

				Cluster *clusterThis   = (idxThis>=0) ? thisTrack->At(idxThis) : nullptr;
				Cluster *clusterSwitch = (idxSwitch>=0) ? switchTrack->At(idxSwitch) : nullptr;
				
				if (layer <= switchTrackAtLayer) { // SWITCH TRACKS
					if (clusterThis)   newTrackB->appendCluster(clusterThis);
					if (clusterSwitch) newTrackA->appendCluster(clusterSwitch);
				}

				else { // KEEP ORIGINAL TRACKS
					if (clusterThis)   newTrackA->appendCluster(clusterThis);
					if (clusterSwitch) newTrackB->appendCluster(clusterSwitch);
				}
			}

			Float_t sin1 = thisTrack->getSinuosity();
			Float_t sin2 = switchTrack->getSinuosity();

			Float_t sin3 = newTrackA->getSinuosity();
			Float_t sin4 = newTrackB->getSinuosity();
			
			Float_t metricOld = quadratureAdd(sin1, sin2);
			Float_t metricNew = quadratureAdd(sin3, sin4);

			if (metricNew < metricOld) {
				newTracks->appendTrack(newTrackA);
				newTracks->appendTrack(newTrackB);
				removeTrackAt(switchTrackIdx);
				removeTrackAt(i);
			}

			ignoreTrack[switchTrackIdx] = true;
			delete newTrackA;
			delete newTrackB;
		}
	}

	cout << "Track::retrogradeTrackImprovement found " << newTracks->GetEntriesFast() << " improved tracks!\n";

	if (newTracks->GetEntriesFast()) {
		for (Int_t i=0; i<newTracks->GetEntriesFast(); i++) {
			appendTrack(newTracks->At(i));
		}
	}
}
