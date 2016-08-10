#ifndef Misalign_h
#define Misalign_h

#include <TObject.h>

#include "GlobalConstants/Constants.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"

struct point {
	Float_t x, y;
}; 

class Misalign : public TObject {
  private:
	  Float_t misalignX[2][24];
	  Float_t misalignY[2][24];

  public:
      Misalign();
      virtual ~Misalign(); 

		void	correctClusters(Clusters * clusters);
		point getMisalign(Cluster * cluster);
		point getMisalign(Int_t assembly, Int_t unit);

//      ClassDef(Misalign,1)
};
#endif
