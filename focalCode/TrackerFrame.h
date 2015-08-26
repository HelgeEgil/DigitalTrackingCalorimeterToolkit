#ifndef FrameTracker_h
#define FrameTracker_h
#include "TClonesArray.h"

#include <vector>
using namespace std;

// A Frame contains many layers

class TrackerFrame : public TClonesArray {

   private:
      TClonesArray frameTracker3D_;

   public:
      TrackerFrame() : frameTracker3D_("Layer", nLayers + nTrackers) { }

   ClassDef(TrackerFrame,1);
};
#endif
