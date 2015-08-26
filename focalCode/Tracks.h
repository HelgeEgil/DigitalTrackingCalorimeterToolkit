#ifndef TrackCollection_h
#define TrackCollection_h
#include "TClonesArray.h"
#include "Track.h"
#include <vector>

class Tracks : public TClonesArray {

   private:
      TClonesArray tracks_;
      TClonesArray clustersWithoutTrack_;

   public:

      Tracks() : tracks_("Track", 1000), clustersWithoutTrack_("Cluster", 5000) {}
      Tracks(Int_t nTracks) : tracks_("Track", nTracks), clustersWithoutTrack_("Cluster", nTracks*5) {}
      virtual ~Tracks(); 

      virtual void appendTrack(Track *copyTrack, Int_t startOffset = 0);
      virtual void appendClustersWithoutTrack(TClonesArray* clustersWithoutTrack);

      virtual void Clear() { tracks_.Clear("C"); }

      virtual Int_t GetEntriesFast() { return tracks_.GetEntriesFast(); }
		virtual Int_t GetEntriesFast(Int_t i) { return At(i)->GetEntriesFast(); }

		Float_t getSinuosity(Int_t i) { return At(i)->getSinuosity(); }
		Float_t getSlopeAngle(Int_t i) { return At(i)->getSlopeAngle(); }
		Float_t getTrackLengthmm(Int_t i) { return At(i)->getTrackLengthmm(); }
		TClonesArray * getClustersWithoutTrack() { return (TClonesArray*) &clustersWithoutTrack_; }


      virtual void extrapolateToLayer0();

      virtual Track* At(Int_t i) { return ((Track*) tracks_.At(i)); }

      virtual void removeCluster(Track *c) { tracks_.Remove((TObject*) c); }
      virtual TObject* removeClusterAt(Int_t i) { return tracks_.RemoveAt(i); }

   ClassDef(Tracks,1);
};
#endif
