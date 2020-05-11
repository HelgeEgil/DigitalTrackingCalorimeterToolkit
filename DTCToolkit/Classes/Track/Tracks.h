#ifndef TrackCollection_h
#define TrackCollection_h

#include <vector>

#include <TClonesArray.h>

#include "Classes/Track/Track.h"
#include <TCollection.h>

using namespace std;

namespace DTC {
class Hits;

struct trackCluster { 
  int track; 
  int cluster;
};

class Tracks : public TObject {

   private:
      TClonesArray tracks_;
      TClonesArray clustersWithoutTrack_;
//      Int_t        EIDindex_[100000];

   public:
      Tracks() : tracks_("DTC::Track", 2000), clustersWithoutTrack_("DTC::Cluster", 10000) { tracks_.SetOwner(kTRUE); clustersWithoutTrack_.SetOwner(kTRUE); tracks_.SetBit(kCanDelete); clustersWithoutTrack_.SetBit(kCanDelete); tracks_.SetBit(kMustCleanup); clustersWithoutTrack_.SetBit(kMustCleanup); }
//      Tracks(Int_t nTracks) : tracks_("DTC::Track", nTracks), clustersWithoutTrack_("DTC::Cluster", nTracks*5) { if (nTracks > 200000) cout << "Remember to increase size of EIDindex array!!!! (now = 100 000)\n"; tracks_.SetOwner(kTRUE); clustersWithoutTrack_.SetOwner(kTRUE); }
      Tracks(Int_t nTracks) : tracks_("DTC::Track", nTracks), clustersWithoutTrack_("DTC::Cluster", nTracks*5) { tracks_.SetOwner(kTRUE); clustersWithoutTrack_.SetOwner(kTRUE); tracks_.SetBit(kCanDelete); clustersWithoutTrack_.SetBit(kCanDelete); tracks_.SetBit(kMustCleanup); clustersWithoutTrack_.SetBit(kMustCleanup); }
      virtual ~Tracks(); 

      // ROOT & I/O     
      virtual void      clearTracks()           { tracks_.Clear("C"); }
      virtual Track   * At(Int_t i)             { return ((Track*) tracks_.At(i)); }
      virtual Cluster * AtCWT(Int_t i)          { return ((Cluster*) clustersWithoutTrack_.At(i)); }
      virtual Track   * Last()                  { return ((Track*) tracks_.At(GetEntriesFast()-1)); }
      virtual void      SetOwner(Bool_t val)    { tracks_.SetOwner(val); }
      virtual Int_t     GetEntriesFast()        { return tracks_.GetEntriesFast(); }
      virtual Int_t     GetEntriesFastCWT()     { return clustersWithoutTrack_.GetEntriesFast(); }
      virtual Int_t     GetEntries()            { return tracks_.GetEntries(); }
      virtual Int_t     GetEntriesFast(Int_t i) { return At(i)->GetEntriesFast(); }
      virtual void      Compress()              { tracks_.Compress(); }
      virtual void      CompressCWT()           { clustersWithoutTrack_.Compress(); }
      virtual void      CompressClusters();
      virtual void      Clear(Option_t * option = "");
      virtual void      sortTracks();
      void              sortTracksByLength() { tracks_.Sort(); }
      virtual Long64_t  Merge(TCollection *tlist);

      Bool_t            isSortable() const { return kTRUE; }

      // Add and remove tracks
      virtual void      removeTrack(Track *t); //   { tracks_.Remove((TObject*) t); }
      virtual TObject*  removeTrackAt(Int_t i); //  { return tracks_.RemoveAt(i); }
      virtual TObject*  removeCWTAt(Int_t i)    { return clustersWithoutTrack_.RemoveAt(i); }
      void              appendTrack(Track *copyTrack, Int_t startOffset = 0);
      void              appendClustersWithoutTrack(TClonesArray *clustersWithoutTrack);
      void              appendClusterWithoutTrack(Cluster * c);
      void              removeHighAngleTracks(Float_t mradLimit);
      void              removeHighAngleTracksRelativeToSpot(Float_t mradLimit, Float_t angleX, Float_t angleY);
      void              removeHighAngularChangeTracks(Float_t mradLimit);
      void              removeTracksStartingInDetector();
      void              removeShortTracks(Int_t len);

      // Retrieve tracks
      TClonesArray    * getClustersWithoutTrack() { return (TClonesArray*) &clustersWithoutTrack_; }
      vector<Int_t>   * getTracksWithConflictClusters();
      vector<Int_t>   * getConflictingTracksFromTrack(Int_t trackIdx);

      // Operations on single Track objects
      Float_t           getSinuosity(Int_t i)                  { return At(i)->getSinuosity(); }
      Float_t           getSlopeAngle(Int_t i)                 { return At(i)->getSlopeAngle(); }
      Float_t           getTrackLengthmm(Int_t i)              { return At(i)->getTrackLengthmm(); }
      Float_t           getTrackScore(Int_t i)                 { return At(i)->getTrackScore(); }
      virtual Bool_t    isUsedClustersInTrack(Int_t i)         { return At(i)->isUsedClustersInTrack(); }
      virtual Int_t     getNumberOfConflictClusters(Int_t i)   { return At(i)->getNumberOfConflictClusters(); }
      virtual void      extrapolateToLayer0();
      virtual void      doTrackFit(Bool_t freeScale = false);
      void              removeNANs();
      void              propagateSecondaryStatus();

      // Search for clusters
      Int_t             getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster);
      vector<Int_t>   * getTracksFromCluster(Cluster * cluster);
      Int_t             getTrackIdxFromFirstLayerEID(Int_t eventID);
      Int_t             getTrackIdxFromLastLayerEID(Int_t eventID);
      Int_t             getTrackIdxFromCluster(Cluster * cluster);
      Int_t             getNMissingClustersWithEventID(Int_t eventID, Int_t afterLayer = -1, Int_t trackID = -1);

      // Calculations on involving all tracks
      void              matchWithEventIDs(Hits *eventIDs);
      void              sortTrackByLayer(Int_t track);
      void              checkLayerOrientation();
      Bool_t            isLastEventIDCloseToFirst(Int_t trackIdx);
      Bool_t            isTrackIncomplete(Track * originalTrack);
//      void              createEIDSortList();
//      Track           * getTrackWithEID(Int_t eid);
      
      // Optimization of the tracking function
      // tracksOptimization.C
      void              splitSharedClusters();
      void              removeTracksLeavingDetector(); 
      void              removeTrackCollisions();
      void              retrogradeTrackImprovement(Clusters * clusters);
      void              removeTracksEndingInBadChannels();
      void              removeNuclearInteractions();
      void              removeThreeSigmaShortTracks();
      void              removeTracksWithMinWEPL(Float_t minWEPL);
      void              removeHighChiSquare(Float_t chi2limit = 210);
      void              fillOutIncompleteTracks(float angleLimit = 0.1);
      void              removeEmptyTracks();
      void              removeTracksEndingInGapRegion();
      void              removeTracksEndingInHalo(Float_t haloRadius = 15);

   ClassDef(Tracks,2)
};
}
#endif
