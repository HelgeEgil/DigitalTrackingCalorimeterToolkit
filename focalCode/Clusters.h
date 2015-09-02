#ifndef ClusterCollection_h
#define ClusterCollection_h
#include <TClonesArray.h>
#include "Tracks.h"
#include "Constants.h"
#include "Cluster.h"
#include <vector>

using namespace std;

class Track;

class Clusters : public TObject {

   private:
      TClonesArray clusters_;
      TClonesArray clustersWithoutTrack_;
      vector<Int_t> layerIndex_;
      Bool_t frameType_; // FIXME Add frameType_ when creating Clusters

   public:
      Clusters(Bool_t frameType = kCalorimeter);

      virtual ~Clusters(); 

      virtual void appendCluster(Float_t x, Float_t y, Int_t layer = -1, Int_t size = -1);
      virtual void appendCluster(Cluster *cluster);
      virtual void appendClusterWithoutTrack(Cluster *cluster);


      virtual Int_t GetEntriesFast() { return clusters_.GetEntriesFast(); }
		virtual Int_t GetEntries() { return clusters_.GetEntries(); }
		virtual Cluster* At(Int_t i) { return ((Cluster*) clusters_.At(i)); }
      
		virtual void clearClusters();

      virtual Float_t getX(Int_t i) { return At(i)->getX(); }
      virtual Float_t getY(Int_t i) { return At(i)->getY(); }
		virtual Int_t getLayer(Int_t i) { return At(i)->getLayer(); }
      virtual Int_t getSize(Int_t i) { return At(i)->getSize(); }
      virtual Int_t getFrameType() { return frameType_; }
      virtual TClonesArray * getClustersWithoutTrack() { return (TClonesArray*) &clustersWithoutTrack_; }

      virtual void removeCluster(Cluster *c) { clusters_.Remove((TObject*) c); }
      virtual TObject* removeClusterAt(Int_t i) { return clusters_.RemoveAt(i); }
      virtual Bool_t removeClusterFromCoords(Float_t x, Float_t y, Int_t layer);
      virtual void removeAllClustersInTrack(Track *track);
      void removeSmallClusters(Int_t size);

      virtual void makeLayerIndex();
      virtual Int_t getFirstIndexOfLayer(UInt_t layer);
      virtual Int_t getLastIndexOfLayer(UInt_t layer);

   	Bool_t isPointOutOfBounds(Cluster *point);

      // in file Clusters_findTracks.C
   	Int_t getLastActiveLayer();
   	Tracks * findTracks();
   	void findTracksFromLayer(Tracks * tracks, Int_t layer);
   	Clusters * findSeeds(Int_t layer);
   	void doNearestClusterTrackPropagation(Track *track, Int_t lastHitLayer);
   	Track * recursiveTrackPropagation(Cluster *cluster, Track currentTrack);
   	Track * nearestClusterTrackPropagation(Cluster *seed);
   	Clusters * findNearestClustersInNextLayer(Cluster *seed);
   	Clusters * findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer);
   	Cluster * getTrackPropagationToLayer(Track *track, Int_t layer);
   	Cluster * findNearestNeighbour(Cluster *projectedPoint);
   	Clusters * findAllClosePointsInNextLayer(Cluster *projectedPoint);
   	Track * findLongestTrack(Tracks *seedTracks);

   ClassDef(Clusters,1);
};
#endif
