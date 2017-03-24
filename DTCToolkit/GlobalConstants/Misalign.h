#ifndef Misalign_h
#define Misalign_h

#include <TObject.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"

namespace DTC {
class Cluster;
class Clusters;

struct chipAlignment {
   Float_t deltaX, deltaY, deltaTheta;
   Int_t idx;
};

class Misalign : public TObject {
  private:
//   Float_t misalignX[2][24];
//   Float_t misalignY[2][24];
     chipAlignment chipAlignmentArray_[96];

  public:
      Misalign();
      virtual ~Misalign(); 

      void           correctClusters(Clusters * clusters);
      chipAlignment  getMisalign(Cluster * cluster);

//      ClassDef(Misalign,1)
};
}
#endif
