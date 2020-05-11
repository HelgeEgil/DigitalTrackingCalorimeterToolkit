#ifndef ClusterCollection_h
#define ClusterCollection_h

#include <vector>

#include <TClonesArray.h>

#include "GlobalConstants/Constants.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Node.h"
#include "Classes/Track/Tracks.h"

using namespace std;

namespace DTC {
class Track;
class Hits;

class Clusters : public TObject {

private:
   TClonesArray   clusters_;
   TClonesArray   clustersWithoutTrack_;
   vector<Int_t>  layerIndex_;
   vector<Int_t>  clustersPerEventID_;
   Bool_t         frameType_; // Not used yet... kCalorimeter or kTracker (for the 4 trackers)
   Float_t        MCSMultiplicationFactor;
   Float_t        kMCSFactorFirstPass;
   Float_t        kMCSFactorSecondPass;
   Float_t        kMCSFactorLastPass1;
   Float_t        kMCSFactorLastPass2;
   Float_t        kMCSFactorLastPass3;

public:
   
   Clusters(Bool_t frameType = kCalorimeter);
   virtual ~Clusters(); 

   // ROOT & I/O
   virtual Cluster * At(Int_t i)                { return ((Cluster*) clusters_.At(i)); }
   virtual Cluster * AtCWT(Int_t i)             { return ((Cluster*) clustersWithoutTrack_.At(i)); }
   virtual Cluster * Last()                     { return ((Cluster*) clusters_.Last()); }
   TClonesArray    * getClustersWithoutTrack()  { return (TClonesArray*) &clustersWithoutTrack_; }
   virtual Int_t     GetEntriesFast()           { return clusters_.GetEntriesFast(); }
   virtual Int_t     GetEntries()               { return clusters_.GetEntries(); }
   virtual Int_t     GetEntriesFastCWT()        { return clustersWithoutTrack_.GetEntriesFast(); }
   virtual Int_t     GetEntriesCWT()            { return clustersWithoutTrack_.GetEntries(); }
   virtual void      Compress()                 { clusters_.Compress(); }
   virtual Int_t     GetEntriesFastLastLayer();
   virtual Int_t     GetEntriesInLayer(Int_t i);
   virtual void      Clear(Option_t * option= "");
   void              sortClusters()             { clusters_.Sort(); }
   
   // Add and remove clusters
   virtual TObject * removeClusterWithoutTrackAt(Int_t i)   { return clustersWithoutTrack_.RemoveAt(i); }
   virtual void      removeCluster(Cluster *c);
   virtual TObject * removeClusterAt(Int_t i);
   void              removeAllClustersInTrack(Track *track);
   void              removeTrackFromClustersWithoutTrack(Track *track);
   void              removeSmallClusters(Int_t size);
   void              removeAllClustersAfterLayer(Int_t afterLayer);
   virtual void      appendCluster(Float_t x, Float_t y, Int_t layer = -1, Int_t size = -1, Int_t eventID = -1, Bool_t isSecondary = false, Int_t PDG = 0);
   virtual void      appendClusterEdep(Float_t x, Float_t y, Int_t layer = -1, Float_t edep = -1, Int_t eventID = -1, Bool_t isSecondary = false, Int_t PDG = 0);
   virtual void      appendCluster(Cluster *cluster);
   virtual void      appendClusterWithoutTrack(Cluster *cluster);
   virtual void      appendTrackToClustersWithoutTrack(Track *track);
   Float_t           removeClustersInGap(Float_t gapSizemm, Float_t gapPosmm);

   // Event ID operations
   Int_t             getClustersForEventID(Int_t eventID);
   void              findNumberOfClustersForEachEventID();
   void              propagateSecondaryStatusFromTop(Int_t eventID = -1);

   // Getters and setters
   virtual Float_t   getX(Int_t i)        { return At(i)->getX(); }
   virtual Float_t   getY(Int_t i)        { return At(i)->getY(); }
   virtual Int_t     getLayer(Int_t i)    { return At(i)->getLayer(); }
   virtual Int_t     getSize(Int_t i)     { return At(i)->getSize(); }
   virtual Int_t     getEventID(Int_t i)  { return At(i)->getEventID(); }
   virtual Int_t     getFrameType()       { return frameType_; }
   virtual Bool_t    isSecondary(Int_t i) { return At(i)->isSecondary(); }
   virtual Bool_t    isUsed(Int_t i)      { return At(i)->isUsed(); }
   virtual void      markUsed(Int_t i)    { At(i)->markUsed(); }
   virtual void      markUnused(Int_t i)  { At(i)->markUnused(); }
   void              markUsedClusters(Track *track);
   Int_t             getFirstIndexOfLayer(UInt_t layer);
   Int_t             getLastIndexOfLayer(UInt_t layer);
   Int_t             getClusterIdx(Float_t x, Float_t y, Int_t layer);
   Int_t             getLastActiveLayer();

   // Longer functions in Clusters.C
   virtual void makeLayerIndex();
   virtual void matchWithEventIDs(Hits * eventIDs);
   virtual void removeHaloAtRadius(Float_t radius);

   // in file findTracks.C 
   // This is the tracking algo used in 2018 WoC paper
   Tracks    * findCalorimeterTracksWithMCTruth();
   Tracks    * findTracksWithRecursiveWeighting();
   void        doRecursiveWeightedTracking(Node * seedNode, vector<Node*> * endNodes, Float_t thisMaxTrackScore);
   void        findSeeds(vector<Int_t> *seeds, Int_t layer, Bool_t kUsedClustersInSeeds = true);
   void        findNearestClustersInNextLayer(Cluster *seed, vector<Int_t> * nextClusters);
   void        findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer, vector<Int_t> * nextClusters);
   Cluster   * findNearestNeighbour(Track* track, Cluster *projectedPoint, Bool_t rejectUsed = true);

   // Work in progress, in file findTracksCandidateAlgorithms.C
//    Tracks    * findTracksWithReverseRecursiveWeighting(); (*)
//    void        doReverseRecursiveWeightedTracking(Node * seedNode, vector<Node*> * endNodes); (*)
//    void        findRemainingTracks(Tracks *tracks); (*)
//    Clusters  * findNearestClustersInLastLayer(Cluster *seed); (*)
//    void        growTrackFromLayer(Track *track, Int_t fromLayer); (*)
//    Track     * findLongestTrack(Tracks *seedTracks); (*)

// Pruned functions, but kept here for posteriority: Earlier versions of tracking algo, as used in 2017 NIMA paper
//    Tracks    * findCalorimeterTracks(); // find function in pre-2018-08 version on GH
//    Tracks    * findCalorimeterTracksAlpide(); // find function in pre-2018-08 version on GH
//    void        findTracksFromLayer(Tracks *tracks, Int_t layer, Bool_t kUsedClustersInSeeds = true); // find function in pre-2018-08 version on GH
//    Track     * trackPropagation(Cluster *seed); // find function in pre-2018-08 version on GH

   ClassDef(Clusters,6)
};
}
#endif
