#include "Cluster.h"
#include "Tracker.h"
#include "TrackerCollection.h"
#include "Constants.h"
#include <iostream>
#include <TClonesArray.h>
#include "Clusters.h"

using namespace std;

// ClassImp(TrackerCollection)

TrackerCollection::~TrackerCollection() {
   // Destructor
   Clear();
}

void TrackerCollection::AddTracker(Tracker *tracker) {
   Int_t i = GetEntriesFast();
   Tracker *t = (Tracker*) fTrackerCollection.ConstructedAt(i);

   t->Set(tracker->GetPreTracker1(), tracker->GetPreTracker2(),
          tracker->GetPostTracker1(), tracker->GetPostTracker2());
}

