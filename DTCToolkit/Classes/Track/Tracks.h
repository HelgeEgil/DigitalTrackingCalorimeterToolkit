#ifndef TrackCollection_h
#define TrackCollection_h

#include <vector>

#include <TClonesArray.h>

#include "Classes/Track/Track.h"

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

   public:
      Tracks() : tracks_("Track", 1000), clustersWithoutTrack_("Cluster", 5000) {}
      Tracks(Int_t nTracks) : tracks_("Track", nTracks), clustersWithoutTrack_("Cluster", nTracks*5) {}
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
      virtual void      Clear(Option_t * = "");

      // Add and remove tracks
      virtual void      removeTrack(Track *t)   { tracks_.Remove((TObject*) t); }
      virtual TObject*  removeTrackAt(Int_t i)  { return tracks_.RemoveAt(i); }
      void              appendTrack(Track *copyTrack, Int_t startOffset = 0);
      void              appendClustersWithoutTrack(TClonesArray *clustersWithoutTrack);

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
      virtual void      doFit();

      // Search for clusters
      Int_t             getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster);
      vector<Int_t>   * getTracksFromCluster(Cluster * cluster);
      Int_t             getTrackIdxFromFirstLayerEID(Int_t eventID);
      Int_t             getTrackIdxFromLastLayerEID(Int_t eventID);
      Int_t             getTrackIdxFromCluster(Cluster * cluster);

      // Calculations on involving all tracks
      void              matchWithEventIDs(Hits *eventIDs);
      void              sortTrackByLayer(Int_t track);
      void              checkLayerOrientation();
      Bool_t            isLastEventIDCloseToFirst(Int_t trackIdx);
      
      // Optimization of the tracking function
      // tracksOptimization.C
      void              splitSharedClusters();
      void              removeTracksLeavingDetector(); 
      void              removeTrackCollisions();
      void              retrogradeTrackImprovement(Clusters * clusters);
      void              removeTracksEndingInBadChannels();

   ClassDef(Tracks,2)
};
}
#endif
