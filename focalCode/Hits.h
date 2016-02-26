#ifndef HitCollection_h
#define HitCollection_h
#include "TClonesArray.h"
#include "Hit.h"
#include "Clusters.h"
#include "Constants.h"
#include <vector>

class Hits : public TObject {

private:
	TClonesArray hits_;
	vector<Int_t> layerIndex_;
	vector<Int_t> verticalIndexOfLayer_;

public:

	Hits(Int_t nHits) : hits_("Hit", nHits) {}
	Hits() : hits_("Hit", kEventsPerRun*200) {}
	virtual ~Hits(); 

	virtual void appendPoint(Int_t x, Int_t y, Int_t layer = -1, Int_t event = -1, Float_t edep = 0);
	virtual void appendHits(Hits *hits);

	virtual void clearHits() { hits_.Clear("C"); }
	virtual void Clear(Option_t * = "");

	virtual Int_t GetEntriesFast() { return hits_.GetEntriesFast(); }
	virtual Int_t GetEntries() { return hits_.GetEntries(); }
	
	virtual Int_t getX(Int_t i) { return At(i)->getX(); }
	virtual Int_t getY(Int_t i) { return At(i)->getY(); }
	virtual Int_t getLayer(Int_t i) { return At(i)->getLayer(); }
	virtual Int_t getEvent(Int_t i) { return At(i)->getEvent(); }
	virtual Float_t getEdep(Int_t i) { return At(i)->getEdep(); }

	virtual Hit* At(Int_t i) { return ((Hit*) hits_.At(i)); }
	virtual Hit* Last() { return ((Hit*) hits_.Last()); }

	virtual Int_t getI(Int_t x, Int_t y);

	Clusters * findClustersFromHits();
	vector<Int_t> findNeighbours(Int_t index);
	vector<Int_t> * findExpandedCluster(Int_t i, vector<Int_t> *checkedIndices);
	void appendExpandedClusterToClusters(vector<Int_t> *expandedCluster, Clusters *clusters);
	void appendExpandedClusterToClusterHitMap(vector<Int_t> *expandedCluster, vector<Hits*> *clusterHitMap);
	void checkAndAppendAllNextCandidates(vector<Int_t> nextCandidates, vector<Int_t> *checkedIndices,
			vector<Int_t> *toCheck, vector<Int_t> *expandedCluster);
	vector<Hits*> * findClustersHitMap();

	void makeLayerIndex();
	Int_t getFirstIndexOfLayer(UInt_t layer);
	Int_t getLastIndexOfLayer(UInt_t layer);
	Int_t getLastActiveLayer();

	void makeVerticalIndexOnLayer(Int_t layer);
	Int_t getFirstIndexBeforeY(Int_t y);
	Int_t getLastIndexAfterY(Int_t y);


ClassDef(Hits,1);
};
#endif
