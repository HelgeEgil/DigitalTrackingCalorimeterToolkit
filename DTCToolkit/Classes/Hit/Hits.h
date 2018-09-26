#ifndef HitCollection_h
#define HitCollection_h

#include <vector>

#include <TClonesArray.h>

#include "Classes/Hit/Hit.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"

namespace DTC {

class Hits : public TObject {

private:
   TClonesArray hits_;
   vector<Int_t> layerIndex_;
   vector<Int_t> verticalIndexOfLayer_;

public:
   Hits(Int_t nHits) : hits_("DTC::Hit", nHits) { hits_.SetOwner(kTRUE); }
   Hits() : hits_("DTC::Hit", kEventsPerRun*200) { hits_.SetOwner(kTRUE); }
   virtual ~Hits(); 

   // ROOT & I/O
   virtual Hit     * At(Int_t i)          { return ((Hit*) hits_.At(i)); }
   virtual Hit     * Last()               { return ((Hit*) hits_.Last()); }
   virtual Int_t     GetEntriesFast()     { return hits_.GetEntriesFast(); }
   virtual Int_t     GetEntries()         { return hits_.GetEntries(); }
   virtual void      clearHits()          { hits_.Clear("C"); }
   virtual void      Clear(Option_t * = "");
   
   // Add and remove hits
   TObject         * removeHitAt(Int_t i) { return hits_.RemoveAt(i); }
   void              appendPoint(Int_t x, Int_t y, Int_t layer = -1, Int_t event = -1, Float_t edep = 0);
   void              appendHits(Hits *hits);

   // Getters and setters
   virtual Int_t     getX(Int_t i)        { return At(i)->getX(); }
   virtual Int_t     getY(Int_t i)        { return At(i)->getY(); }
   virtual Int_t     getLayer(Int_t i)    { return At(i)->getLayer(); }
   virtual Int_t     getEventID(Int_t i)  { return At(i)->getEventID(); }
   virtual Float_t   getEdep(Int_t i)     { return At(i)->getEdep(); }
   virtual Int_t     getI(Int_t x, Int_t y);
   
   // Layer indexing - optimizstion
   void              makeLayerIndex();
   Int_t             getFirstIndexOfLayer(UInt_t layer);
   Int_t             getLastIndexOfLayer(UInt_t layer);
   Int_t             getLastActiveLayer();
   void              makeVerticalIndexOnLayer(Int_t layer);
   Int_t             getFirstIndexBeforeY(Int_t y);
   Int_t             getLastIndexAfterY(Int_t y);

   // In file findClusters.C
   // Methods to find clusters from the hits
   Clusters        * findClustersFromHits();
   vector<Int_t>   * findNeighbours(Int_t index);
   vector<Int_t>   * getAllNeighboursFromCluster(Int_t i, vector<Int_t> *checkedIndices);
   void              appendNeighboursToClusters(vector<Int_t> *expandedCluster, Clusters *clusters);
   void              checkAndAppendAllNextCandidates(vector<Int_t> *nextCandidates, vector<Int_t> *checkedIndices, vector<Int_t> *toCheck, vector<Int_t> *expandedCluster);
   
   // In find findClusters.C
   // Methods to find cluster hitmap (to visualize cluster shapes)
   vector<Hits*>   * findClustersHitMap();
   void              appendExpandedClusterToClusterHitMap(vector<Int_t> *expandedCluster, vector<Hits*> *clusterHitMap);

   ClassDef(Hits,2)
};
}
#endif
