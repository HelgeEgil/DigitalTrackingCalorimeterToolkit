#ifndef TrackCollection_h
#define TrackCollection_h
#include "TClonesArray.h"
#include "Track.h"
// 	#include "Clusters.h"
#include <vector>

using namespace std;

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

      virtual void appendTrack(Track *copyTrack, Int_t startOffset = 0);
      virtual void appendClustersWithoutTrack(TClonesArray *clustersWithoutTrack);

      virtual void clearTracks() { tracks_.Clear("C"); }
      virtual void Clear(Option_t * = "");
      virtual void SetOwner(Bool_t val) { tracks_.SetOwner(val); }

      virtual Int_t GetEntriesFast() { return tracks_.GetEntriesFast(); }
      virtual Int_t GetEntries() { return tracks_.GetEntries(); }
      virtual Int_t GetEntriesFast(Int_t i) { return At(i)->GetEntriesFast(); }

      Float_t getSinuosity(Int_t i) { return At(i)->getSinuosity(); }
      Float_t getSlopeAngle(Int_t i) { return At(i)->getSlopeAngle(); }
      Float_t getTrackLengthmm(Int_t i) { return At(i)->getTrackLengthmm(); }
      Float_t getTrackScore(Int_t i) { return At(i)->getTrackScore(); }
      TClonesArray * getClustersWithoutTrack() { return (TClonesArray*) &clustersWithoutTrack_; }

      virtual void extrapolateToLayer0();
      virtual void splitSharedClusters();
      Int_t getClosestCluster(vector<trackCluster> clusters, Cluster* interpolatedCluster);


      virtual Track* At(Int_t i) { return ((Track*) tracks_.At(i)); }

      virtual void removeTrack(Track *t) { tracks_.Remove((TObject*) t); }
      virtual TObject* removeTrackAt(Int_t i) { return tracks_.RemoveAt(i); }
      void sortTrackByLayer(Int_t track);

      
//       virtual void drawClusters(Clusters *a, Clusters *b);
	  
   ClassDef(Tracks,1);
};
#endif
