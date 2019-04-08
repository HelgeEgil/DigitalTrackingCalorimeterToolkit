#ifndef findTracks_cxx
#define findTracks_cxx

#include <vector>
#include<fstream>
#include<iostream>

#include <string>
#include <sstream>


//#include <TClonesArray.h>

#include "Classes/Cluster/Clusters.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Node.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "HelperFunctions/Tools.h"

#include <TH2.h>
#include <TH3.h>
#include <TPolyLine.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TPaveStats.h>
#include <TView.h>
#include <TLeaf.h>
#include <TArrow.h>
#include <TF1.h>
#include <Math/ProbFunc.h>


using namespace std;
using namespace DTC;

Tracks * Clusters::findTracksWithRecursiveWeighting() {

   Track     * track = nullptr;
   Tracks    * tracks = new Tracks(50 * 200 * kEventsPerRun);
  
   Int_t       spotSize = 33;
   // for pencil beams.
   // 22 = 2x2 mm, 33 = 3x3 mm etc.
   // See WoC publication for description of functions below

   if      (spotSize == 22) kMaxTrackScore = 0.433 * pow(kEventsPerRun, -0.172);
   else if (spotSize == 33) kMaxTrackScore = 0.469 * pow(kEventsPerRun, -0.176);
   else if (spotSize == 42) kMaxTrackScore = 0.397 * pow(kEventsPerRun, -0.150);
   else if (spotSize == 55) kMaxTrackScore = 0.455 * pow(kEventsPerRun, -0.162);
   else if (spotSize == 44) kMaxTrackScore = 0.460 * pow(kEventsPerRun, -0.168);
   showDebug("makeLayerIndex...");
   makeLayerIndex();
   showDebug("ok!\nfillMSCradiusList...");
   showDebug("ok!\n");
   kMCSFactor = 25;

   if (kAbsorberThickness == 2) kMaxTrackScore *= 2; // phenomenological factor ..? 

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) continue;
      appendClusterWithoutTrack(At(i));
   }

    Int_t FirstLayer = nLayers-1;
    Int_t LastLayer = 0;
 
    for (Int_t s=FirstLayer; s>=LastLayer; s--) {
    Cluster   * nextCluster = nullptr;
    Float_t     clusterScore;
    Node      * nextNode = nullptr;
    Node      * seedNode = nullptr;
    vector<Int_t> * nextClusters = new vector<Int_t>;
    vector<Int_t> * seeds = new vector<Int_t>;
    nextClusters->reserve(100);
   
    showDebug("findSeeds from " << GetEntriesFast() << " clusters." << endl);
    findSeeds(seeds, s, false); // starts in s, false-> does not use a cluster if it has already been used
    //cout<< "         size seeds in "<<s<<": "<< seeds->size()<< endl;
    showDebug("Found " << seeds->size() << " seeds in first go.\n");

    for (UInt_t i=0; i<seeds->size(); i++) {
        
       
        Cluster *seed = At(seeds->at(i)); //relates each seed to the cluster it belongs
        nextClusters->clear();

        findNearestClustersInNextLayer(seed, nextClusters);
    
        seedNode = new Node(nullptr, seed, 0); //   Node(Node *parent, Cluster *connectedCluster, Float_t score);
        
        showDebug("Found " << nextClusters->size() << " nextClusters first time.\n");
       // cout<< "size next clusters: "<< nextClusters->size() << endl;
       for (UInt_t j=0; j<nextClusters->size(); j++) {
           
            nextCluster = At(nextClusters->at(j)); // analyze cluster by cluster
            clusterScore = seedNode->getNextScore(nextCluster); // getParent-getChild should not depend on the starting point
            clusterScore /= 100;
            if (std::isnan(clusterScore)) cout << "clusterScore at " << *nextCluster << " isNan!\n";
            
            if (clusterScore < kMaxTrackScore) {
                nextNode = new Node(seedNode, nextCluster, clusterScore); // initial vector is 10 % weighted
                seedNode->addChild(nextNode);
                showDebug("Adding nextNode\n");
            }
        }
        seedNode->markExplored();
        vector<Node*> * endNodes = new vector<Node*>;
        //endNodes->reserve(kEventsPerRun * 5);
        endNodes->reserve(200* 50 * kEventsPerRun);
        seedNode->getUnexploredEndNodes(endNodes); //EndNodes=nodes you are studying
        doRecursiveWeightedTracking(seedNode, endNodes); // all other layers, the first part is only done in the first layer
        showDebug("..ok!\ngetBestTrack...");
        track = seedNode->getBestTrack();
        showDebug("ok!\n");
        
        if ((track->GetEntries() >= (s+1)) && (track->GetEntries() >1)){
            //cout << "value of s is: " << s <<endl;
            //cout << " size track: " << track->GetEntries() <<endl;
            if (track->Last()->getLayer() != 0) {
               printf("Track without last layer...\n");
               cout << *track << endl;
            }
            track->sortTrack();
            tracks->appendTrack(track);
            removeTrackFromClustersWithoutTrack(track);
            markUsedClusters(track);
        }
        delete track;
        seedNode->deleteNodeTree();
        delete seedNode;
        delete endNodes;
    }
    delete seeds;
    delete nextClusters;
}
    return tracks;
}

void Clusters::doRecursiveWeightedTracking(Node * seedNode, vector<Node*> * endNodes) {
   // Input: vector of nodes
   // Find all the next potential nodes segments w/acceptable score
   // Output: Vector of nodes (new end nodes)

   //cout<<"doRecursiveWeightedTracking from layer " << endNodes->at(0)->getCluster()->getLayer() << " with "  << seedNode->getNChildren() << " children" <<endl;
   Node    * thisNode;
   Cluster * nextCluster;
   Float_t   nextScore, nextAngle;

 
   for (UInt_t i=0; i<endNodes->size(); i++) { // All identified endpoints to the tree so far
      thisNode = endNodes->at(i);
      thisNode->markExplored();
      Int_t searchLayer = thisNode->getCluster()->getLayer() - 1;  // +1 because you are studying the next layer
      Int_t idxFrom = getFirstIndexOfLayer(searchLayer); //Only clusters in this layer-first cluster. Clusters contains all the cluster of the detector and depending on the index, these clusters corresponds to one layer or another.
      Int_t idxTo = getLastIndexOfLayer(searchLayer);// only clusters in this layer-last cluster
      if (idxFrom < 0) continue; // no more clusters in deeper layers

      for (Int_t j=idxFrom; j<idxTo; j++) { // Loop through possible additions to the track
         nextCluster = At(j);
         if (!nextCluster) continue;
        
         nextScore = thisNode->getNextScore(nextCluster);
         nextAngle = thisNode->getNodeAngle(nextCluster);

          if ((nextScore < kMaxTrackScore || (nextAngle < 0.03)) ){// || (nextAngle < kMaxTrackAngle)) {
              thisNode->addChild(new Node(thisNode, nextCluster, nextScore)); // it is either appended to the tree or deleted if thisNode is full
          }
      }
   }
   Int_t CurrentLayer =thisNode->getCluster()->getLayer();
   //cout << "Max angle in layer " << CurrentLayer << " is : " << MaxAngle << endl;
   endNodes->clear();
   seedNode->getUnexploredEndNodes(endNodes);
  if (endNodes->size() > 0){
      doRecursiveWeightedTracking(seedNode, endNodes);
   }
}

Tracks * Clusters::findCalorimeterTracksWithMCTruth() {
   // Now clusters is sorted by event ID
   // Loop through list, and fill to track when new event ID is encountered

   Cluster   * cluster = nullptr;
   Track     * track = new Track();
   Tracks    * tracks = new Tracks(kEventsPerRun * 5);
   Int_t       lastEventID, eventID;

   lastEventID = At(0)->getEventID();

   for (Int_t i=0; i<GetEntriesFast(); i++) {
      if (!At(i)) {
         cout << "Clusters::findCalorimeterTracksWithMCTruth could not find cluster at " << i << endl;
         continue;
      }

      cluster = At(i);
      eventID = cluster->getEventID();

      if (eventID != lastEventID) {
         // Cluster from new event ID encountered, store what we have already
         showDebug("Track at  is full, storing the pointer to tracks having currently " << tracks->GetEntriesFast() << " elements.\n");
         tracks->appendTrack(track);
         track->Clear("C");
      }

      showDebug("Appending cluster at " << *cluster << " with  number " << track->GetEntriesFast() + 1 << " to track.\n");
      track->appendCluster(cluster);

      lastEventID = eventID;
   }

   // Store last track
   tracks->appendTrack(track);

   delete track;
    
   return tracks;
}

void Clusters::findSeeds(vector<Int_t> * seeds, Int_t layer, Bool_t kUsedClustersInSeeds) {
   showDebug("finding layerIdxFrom...");
   Int_t layerIdxFrom = getFirstIndexOfLayer(layer);
   showDebug("to...");
   Int_t layerIdxTo = getLastIndexOfLayer(layer);
   showDebug("ok!\n");

   showDebug("findSeeds: layerIdxFrom = " << layerIdxFrom << ", layerIdxTO = " << layerIdxTo << endl);

   if (layerIdxFrom>=0) { 
      for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
         if (!At(i)) continue;
         if (!kUsedClustersInSeeds && isUsed(i)) continue;
         

         seeds->push_back(i);
      }
   }
}

void Clusters::findNearestClustersInNextLayer(Cluster *seed, vector<Int_t> * nextClusters) {
   for (Int_t skipLayers=0; skipLayers<2; skipLayers++) { // Loop to enable tracks skipping a layer without data
      Int_t nextLayer = seed->getLayer() - 1 - skipLayers;
      findClustersFromSeedInLayer(seed, nextLayer, nextClusters);
      if (nextClusters->size() > 0) break;
   }
}

void Clusters::findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer, vector<Int_t> * nextClusters) { 
   Int_t layerIdxFrom = getFirstIndexOfLayer(nextLayer);
   Int_t layerIdxTo = getLastIndexOfLayer(nextLayer);
   Float_t maxAngle, thisAngle;

   if (kUseEmpiricalMCS) {
      maxAngle = 3 * getEmpiricalMCSAngle(nextLayer - 1); // This function is only run once, to find seed candidates in the 2nd layer - here it's okay with a high angle
      if (layerIdxFrom >= 0) {
         for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
            if (!At(i)) { continue; }
            thisAngle = getDotProductAngle(seed, seed, At(i));
             if ((thisAngle < maxAngle)&& !isUsed(i)){
               nextClusters->push_back(i);
            }
         }
      }
   }

   else { // Use "old" MCS estimation method
      maxAngle = 3 * getSearchRadiusForLayer(nextLayer) * 0.75 * MCSMultiplicationFactor;
      if (layerIdxFrom >= 0) {
         for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
            if (!At(i)) { continue; }
            thisAngle = diffmmXY(seed, At(i));
            if (thisAngle < maxAngle) {
               nextClusters->push_back(i);
            }
         } 
      }
   }
}

Cluster * Clusters::findNearestNeighbour(Track *track, Cluster *projectedPoint, Bool_t rejectUsed) {
   // Finds the cluster in layer projectedPoint->getLayer() closest to the projectedPoint,
   // based on the distance from a projected vector from *track. (angle minimization)
   //
   // In some cases we don't know which track to connect the point to - in that case
   // Track * track = nullptr, and the legacy method (diffmmXY) is used
   Cluster *nearestNeighbour = nullptr;
   Float_t  thisAngle, maxAngle;
   Bool_t   kFoundNeighbour = false;
   Bool_t   reject = false;
   Int_t    searchLayer = projectedPoint->getLayer();
   Int_t    layerIdxFrom = getFirstIndexOfLayer(searchLayer);
   Int_t    layerIdxTo = getLastIndexOfLayer(searchLayer);

   if (kUseEmpiricalMCS && track)   maxAngle = getEmpiricalMCSAngle(searchLayer-1);
   else                             maxAngle = getSearchRadiusForLayer(searchLayer) * MCSMultiplicationFactor;


   if (layerIdxFrom < 0) return 0;

   for (Int_t i=layerIdxFrom; i<layerIdxTo; i++) {
      if (!At(i))
         continue;

      if (kUseEmpiricalMCS && track)   thisAngle = getDotProductAngle(track->At(track->GetEntriesFast()-2), track->Last(), At(i));
      else                             thisAngle = diffmmXY(projectedPoint, At(i));

      reject = (At(i)->isUsed() && rejectUsed);

      if (thisAngle < maxAngle && !reject) {
         nearestNeighbour = At(i);
         maxAngle = thisAngle;
         kFoundNeighbour = true;
      }
   }
   
   if (kFoundNeighbour)   return nearestNeighbour;
   else                   return nullptr;
}

#endif
