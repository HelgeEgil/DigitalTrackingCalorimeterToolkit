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
	TClonesArray	clusters_;
	TClonesArray	clustersWithoutTrack_;
	vector<Int_t>	layerIndex_;
	Bool_t			frameType_; // Not used yet... kCalorimeter or kTracker (for the 4 trackers)
	Float_t			MCSMultiplicationFactor;

public:
	
	Clusters(Bool_t frameType = kCalorimeter);
	virtual ~Clusters(); 

	// ROOT & I/O
	virtual Cluster *	At(Int_t i)						{ return ((Cluster*) clusters_.At(i)); }
	TClonesArray	 *	getClustersWithoutTrack()	{ return (TClonesArray*) &clustersWithoutTrack_; }
	virtual Int_t		GetEntriesFast()				{ return clusters_.GetEntriesFast(); }
	virtual Int_t		GetEntries()					{ return clusters_.GetEntries(); }
	virtual Int_t		GetEntriesFastCWT()			{ return clustersWithoutTrack_.GetEntriesFast(); }
	virtual void		Compress()						{ clusters_.Compress(); }
	virtual Int_t		GetEntriesFastLastLayer();
	virtual Int_t		GetEntriesInLayer(Int_t i);
	virtual void		clearClusters();
	virtual void		Clear(Option_t * = "");
	
	// Add and remove clusters
	virtual void		removeCluster(Cluster *c)					{ clusters_.Remove((TObject*) c); }
	virtual TObject *	removeClusterAt(Int_t i)					{ return clusters_.RemoveAt(i); }
	virtual TObject *	removeClusterWithoutTrackAt(Int_t i)	{ return clustersWithoutTrack_.RemoveAt(i); }
	void					removeAllClustersInTrack(Track *track);
	void					removeTrackFromClustersWithoutTrack(Track *track);
	void					removeSmallClusters(Int_t size);
	void					removeAllClustersAfterLayer(Int_t afterLayer);
	virtual void		appendCluster(Float_t x, Float_t y, Int_t layer = -1, Int_t size = -1, Int_t eventID = -1);
	virtual void		appendClusterEdep(Float_t x, Float_t y, Int_t layer = -1, Float_t edep = -1, Int_t eventID = -1);
	virtual void		appendCluster(Cluster *cluster);
	virtual void		appendClusterWithoutTrack(Cluster *cluster);

	// Getters and setters
	virtual Float_t	getX(Int_t i)			{ return At(i)->getX(); }
	virtual Float_t	getY(Int_t i)			{ return At(i)->getY(); }
	virtual Int_t		getLayer(Int_t i)		{ return At(i)->getLayer(); }
	virtual Int_t		getSize(Int_t i)		{ return At(i)->getSize(); }
	virtual Int_t		getEventID(Int_t i)	{ return At(i)->getEventID(); }
	virtual Int_t		getFrameType()			{ return frameType_; }
	virtual Bool_t		isUsed(Int_t i)		{ return At(i)->isUsed(); }
	virtual void		markUsed(Int_t i)		{ At(i)->markUsed(); }
	virtual void		markUnused(Int_t i)	{ At(i)->markUnused(); }
	void					markUsedClusters(Track *track);
	Int_t					getFirstIndexOfLayer(UInt_t layer);
	Int_t					getLastIndexOfLayer(UInt_t layer);
	Int_t					getClusterIdx(Float_t x, Float_t y, Int_t layer);
	Int_t					getLastActiveLayer();

	// Longer functions in Clusters.C
	virtual void makeLayerIndex();
	virtual void matchWithEventIDs(Hits * eventIDs);

	// in file findTracks.C
	Tracks    * findCalorimeterTracks();
	void			findTracksFromLayer(Tracks *tracks, Int_t layer, Bool_t kUsedClustersInSeeds = true);
	Clusters	 * findSeeds(Int_t layer, Bool_t kUsedClustersInSeeds = true);
	Track		 *	trackPropagation(Cluster *seed);
	Clusters  *	findNearestClustersInNextLayer(Cluster *seed);
	Clusters  *	findClustersFromSeedInLayer(Cluster *seed, Int_t nextLayer);
	void			growTrackFromLayer(Track *track, Int_t fromLayer);
	Cluster   *	findNearestNeighbour(Cluster *projectedPoint, Bool_t rejectUsed = true);
	Track     *	findLongestTrack(Tracks *seedTracks);

	ClassDef(Clusters,5)
};
#endif
