#ifndef TrackerCollection_h
#define TrackerCollection_h
#include "TClonesArray.h"
#include "Cluster.h"
#include "Tracker.h"

#include <vector>
using namespace std;

class TrackerCollection : public TClonesArray {

   private:
      TClonesArray fTrackerCollection;

   public:

      TrackerCollection(Int_t nTrackers) : fTrackerCollection("Tracker", nTrackers) {}
		TrackerCollection() : fTrackerCollection("Tracker", 500) {}
      virtual ~TrackerCollection(); 

      virtual void Clear() { fTrackerCollection.Clear("C"); }

      virtual Int_t GetEntriesFast() { return fTrackerCollection.GetEntriesFast(); }
		virtual Int_t GetEntries() { return fTrackerCollection.GetEntries(); }
		virtual Tracker* At(Int_t i) { return ((Tracker*) fTrackerCollection.At(i)); }

      virtual void AddTracker(Tracker* tracker);
   
      virtual void Remove(Tracker *t) { fTrackerCollection.Remove((TObject*) t); }
      virtual TObject* RemoveAt(Int_t i) { return fTrackerCollection.removeClusterAt(i); }

      virtual Cluster* GetPreTracker1(Int_t i) { return At(i)->GetPreTracker1(); }
      virtual Cluster* GetPreTracker2(Int_t i) { return At(i)->GetPreTracker2(); }
      virtual Cluster* GetPostTracker1(Int_t i) { return At(i)->GetPostTracker1(); }
      virtual Cluster* GetPostTracker2(Int_t i) { return At(i)->GetPostTracker2(); }

   ClassDef(TrackerCollection,1);
};
#endif
