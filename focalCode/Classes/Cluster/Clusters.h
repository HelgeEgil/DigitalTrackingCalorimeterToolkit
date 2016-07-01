#ifndef ClusterCollection_h
#define ClusterCollection_h

#include <vector>

#include <TClonesArray.h>

#include "GlobalConstants/Constants.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Track/Tracks.h"

using namespace std;

class Track;
class Hits;

class Clusters : public TObject {

private:
	TClonesArray clusters_;
	TClonesArray clustersWithoutTrack_;
	TClonesArray conflictClusters_;
	vector<Int_t> layerIndex_;
	Bool_t frameType_; // FIXME Add frameType_ when creating Clusters

public:
	Clusters(Bool_t frameType = kCalorimeter);

	virtual ~Clusters(); 

	virtual void appendCluster(Float_t x, Float_t y, Int_t layer = -1, Int_t size = -1);
	virtual void appendCluster(Cluster *cluster);
	virtual void appendClusterWithoutTrack(Cluster *cluster);
	virtual void appendConflictClusters(Clusters * clusters);

	virtual Int_t GetEntriesFast() { return clusters_.GetEntriesFast(); }
	virtual Int_t GetEntries() { return clusters_.GetEntries(); }
	virtual Cluster* At(Int_t i) { return ((Cluster*) clusters_.At(i)); }
	
	virtual void clearClusters();
	virtual void Clear(Option_t * = "");

	virtual Int_t GetEntriesFastLastLayer();
	virtual void markUsed(Int_t i) { At(i)->markUsed(); }
	virtual void markUnused(Int_t i) { At(i)->markUnused(); }
	virtual Bool_t isUsed(Int_t i) { return At(i)->isUsed(); }

	virtual Float_t getX(Int_t i) { return At(i)->getX(); }
	virtual Float_t getY(Int_t i) { return At(i)->getY(); }
	virtual Int_t getLayer(Int_t i) { return At(i)->getLayer(); }
	virtual Int_t getSize(Int_t i) { return At(i)->getSize(); }
	virtual Int_t getEventID(Int_t i) { return At(i)->getEventID(); }
	virtual Int_t getFrameType() { return frameType_; }
	virtual TClonesArray * getClustersWithoutTrack() { return (TClonesArray*) &clustersWithoutTrack_; }
	virtual TClonesArray * getConflictClusters() { return (TClonesArray*) &conflictClusters_; }

	virtual void removeCluster(Cluster *c) { clusters_.Remove((TObject*) c); }
	virtual TObject* removeClusterAt(Int_t i) { return clusters_.RemoveAt(i); }
	virtual TObject* removeClusterWithoutTrackAt(Int_t i) { return clustersWithoutTrack_.RemoveAt(i); }
	virtual Bool_t removeClusterFromCoords(Float_t x, Float_t y, Int_t layer);
	virtual void removeAllClustersInTrack(Track *track);
	virtual void removeAllClustersInTrackFromClustersWithoutTrack(Track *track);
	void removeSmallClusters(Int_t size); 
	virtual void matchWithEventIDs(Hits * eventIDs);

	void markUsedClusters(Track *track);

	virtual void makeLayerIndex();
	virtual Int_t getFirstIndexOfLayer(UInt_t layer);
	virtual Int_t getLastIndexOfLayer(UInt_t layer);
	Int_t getClusterIdx(Float_t x, Float_t y, Int_t layer);

	// in file findTracks.C
	Int_t getLastActiveLayer();
	Tracks * findCalorimeterTracks();
	void findTracksFromLayer(Tracks *tracks, Int_t layer, Bool_t kUsedClustersInSeeds = true);
	Clusters * findSeeds(Int_t layer, Bool_t kUsedClustersInSeeds = true);
	void doNearestClusterTrackPropagation(Track *track, Int_t lastHitLayer);
	Track * nearestClusterTrackPropagation(Cluster *seed);
	Clusters * findNearestClustersInNextLayer(Cluster *seed);
	Clusters * findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer);
	Cluster * findNearestNeighbour(Cluster *projectedPoint, Bool_t rejectUsed = true);
	Track * findLongestTrack(Tracks *seedTracks);

	ClassDef(Clusters,4);
};
#endif
