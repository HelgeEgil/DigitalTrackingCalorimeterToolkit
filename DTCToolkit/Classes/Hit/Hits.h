#ifndef HitCollection_h
#define HitCollection_h

#include <vector>

#include <TClonesArray.h>
#include <TCollection.h>

#include "Classes/Hit/Hit.h"
#include "Classes/Cluster/Clusters.h"
#include "GlobalConstants/Constants.h"

namespace DTC {

class Hits : public TObject {

private:
   TClonesArray hits_;
   vector<Int_t> layerIndex_;
//   vector<Int_t> verticalIndexOfLayer_;

public:
   Hits(Int_t nHits) : hits_("DTC::Hit", nHits) { hits_.SetOwner(kTRUE); hits_.SetBit(kCanDelete); hits_.SetBit(kMustCleanup); }
   Hits() : hits_("DTC::Hit", kEventsPerRun*1000) { hits_.SetOwner(kTRUE); hits_.SetBit(kCanDelete); hits_.SetBit(kMustCleanup); }
   virtual ~Hits();

   // ROOT & I/O
   virtual Hit     * At(Int_t i)          { return ((Hit*) hits_.At(i)); }
   virtual Hit     * Last()               { return ((Hit*) hits_.Last()); }
   virtual Int_t     GetEntriesFast()     { return hits_.GetEntriesFast(); }
   virtual Int_t     GetEntries()         { return hits_.GetEntries(); }
   virtual void      clearHits()          { hits_.Clear("C"); }
   virtual void      Clear(Option_t * = "");
   virtual void      Compress()           { hits_.Compress(); }
   void              sortHits()           { hits_.Sort(); }
   Long64_t          Merge(TCollection *hlist);
   
   // Add and remove hits
   TObject         * removeHitAt(Int_t i) { return hits_.RemoveAt(i); }
   void              appendPoint(Int_t x, Int_t y, Int_t layer = -1, Float_t edep = 0, Int_t event = -1, Bool_t isSecondary = false, Int_t PDG = 0);
   void              appendHits(Hits *hits);
   void              appendHit(Hit *hit);
   void              removeHaloAtRadius(Float_t radius);

   // Getters and setters
   virtual Int_t     getX(Int_t i)        { return At(i)->getX(); }
   virtual Int_t     getY(Int_t i)        { return At(i)->getY(); }
   virtual Float_t   getXmm(Int_t i)      { return At(i)->getXmm(); }
   virtual Float_t   getYmm(Int_t i)      { return At(i)->getYmm(); }
   virtual Int_t     getLayer(Int_t i)    { return At(i)->getLayer(); }
   virtual Int_t     getEventID(Int_t i)  { return At(i)->getEventID(); }
   virtual Float_t   getEdep(Int_t i)     { return At(i)->getEdep(); }
   virtual Bool_t    isSecondary(Int_t i) { return At(i)->isSecondary(); }
   virtual Int_t     getPDG(Int_t i)      { return At(i)->getPDG(); }
   virtual Int_t     getI(Int_t x, Int_t y);
   void              propagateSecondaryStatusFromTop(Int_t eventID = -1);
   
   // Layer indexing - optimizstion
   Int_t             findLayerIndex(Int_t findLayer); // Without precaching
   void              makeLayerIndex();
   Int_t             getFirstIndexOfLayer(UInt_t layer);
   Int_t             getLastIndexOfLayer(UInt_t layer);
   Int_t             getLastActiveLayer();
//   void              makeVerticalIndexOnLayer(Int_t layer);
//   Int_t             getFirstIndexBeforeY(Int_t y);
//   Int_t             getLastIndexAfterY(Int_t y);

   // In file findClusters.C
   // Methods to find clusters from the hits
   void              findClustersFromHits(Clusters * clusters);
   vector<Int_t>   * findNeighbours(Int_t index);
   vector<Int_t>   * getAllNeighboursFromCluster(Int_t i, vector<Int_t> *checkedIndices);
   void              appendNeighboursToClusters(vector<Int_t> *expandedCluster, Clusters *clusters);
   void              checkAndAppendAllNextCandidates(vector<Int_t> *nextCandidates, vector<Int_t> *checkedIndices, vector<Int_t> *toCheck, vector<Int_t> *expandedCluster);
   
   // In find findClusters.C
   // Methods to find cluster hitmap (to visualize cluster shapes)
   vector<Hits*>   * findClustersHitMap();
   void              appendExpandedClusterToClusterHitMap(vector<Int_t> *expandedCluster, vector<Hits*> *clusterHitMap);

   ClassDef(Hits,3)
};
}
#endif
