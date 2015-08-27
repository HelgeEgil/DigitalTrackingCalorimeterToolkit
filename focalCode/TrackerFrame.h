#ifndef TrackerFrame_h
#define TrackerFrame_h

#include "TClonesArray.h"
#include "Constants.h"
#include "Layer.h"
#include <vector>

using namespace std;

// A Frame contains many layers

class TrackerFrame : public TClonesArray {

   private:
      TClonesArray frameTracker3D_;

   public:
      TrackerFrame() : frameTracker3D_("Layer", nTrackers) { }

   ClassDef(TrackerFrame,1);
};
#endif
