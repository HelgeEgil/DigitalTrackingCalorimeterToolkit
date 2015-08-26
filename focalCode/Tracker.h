#ifndef Tracker_h
#define Tracker_h
#include "Cluster.h"

#include <vector>
using namespace std;

// The tracking system consists of four layers
// Two in front of phantom / object to be reconstructed
//    pre-trackers
// Two behind (downstream)
//    post-trackers
// Both position and angle of individual protons can be
// recorded with such a system
//
// In this class, some simple reconstruction tasks can also be done

class Tracker : public TObject {

   private:
      // pointers since they could be null
      Cluster PreTracker1;
      Cluster PreTracker2;
      Cluster PostTracker1;
      Cluster PostTracker2;

   public:

      Tracker() {}
      Tracker(Cluster *pre1, Cluster *pre2, Cluster *post1, Cluster *post2);
      
      virtual ~Tracker(); 
		
		virtual void Set(Cluster *pre1, Cluster *pre2, Cluster *post1, Cluster *post2);

      virtual void SetPreTracker1(Cluster *cluster);
      virtual void SetPreTracker2(Cluster *cluster);
      virtual void SetPostTracker1(Cluster *cluster);
      virtual void SetPostTracker2(Cluster *cluster);

      virtual Cluster* GetPreTracker1() { return (Cluster*) &PreTracker1; }
      virtual Cluster* GetPreTracker2() { return (Cluster*) &PreTracker2; }
      virtual Cluster* GetPostTracker1() { return (Cluster*) &PostTracker1; }
      virtual Cluster* GetPostTracker2() { return (Cluster*) &PostTracker2; }

      virtual Float_t GetPrePositionX() { return PreTracker2.getX(); }
      virtual Float_t GetPrePositionY() { return PreTracker2.getY(); }
      virtual Float_t GetPostPositionX() { return PostTracker1.getX(); }
      virtual Float_t GetPostPositionY() { return PostTracker1.getY(); }
      
      virtual Float_t GetPrePositionXmm() { return PreTracker2.getXmm(); }
      virtual Float_t GetPrePositionYmm() { return PreTracker2.getYmm(); }
      virtual Float_t GetPostPositionXmm() { return PostTracker1.getXmm(); }
      virtual Float_t GetPostPositionYmm() { return PostTracker1.getYmm(); }

      virtual Float_t GetAngle(Cluster *c1, Cluster *c2);
      virtual Cluster GetPreDerivative();
      virtual Cluster GetPostDerivative();


      virtual Float_t GetPreAngle();
      virtual Float_t GetPostAngle();

      friend ostream& operator<<(ostream &os, Tracker &t);

   ClassDef(Tracker,1);
};
#endif
