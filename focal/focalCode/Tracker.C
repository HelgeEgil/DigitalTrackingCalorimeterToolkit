#include "Cluster.h"
#include "Tracker.h"
#include "Constants.h"
#include <iostream>
#include <cmath>

using namespace std;


Tracker::~Tracker() {
   // Destructor
   Clear();
}

Tracker::Tracker(Cluster *pre1, Cluster *pre2, Cluster *post1, Cluster *post2) {
   PreTracker1.set(pre1);
   PreTracker2.set(pre2);
   PostTracker1.set(post1);
   PostTracker2.set(post2);
}

void Tracker::Set(Cluster *pre1, Cluster *pre2, Cluster *post1, Cluster *post2) {
   PreTracker1.set(pre1);
   PreTracker2.set(pre2);
   PostTracker1.set(post1);
   PostTracker2.set(post2);
}


void Tracker::SetPreTracker1(Cluster *cluster) {
   PreTracker1.set(cluster->getX(), cluster->getY(), -4, cluster->getSize());
}

void Tracker::SetPreTracker2(Cluster *cluster) {
   PreTracker2.set(cluster->getX(), cluster->getY(), -3, cluster->getSize());
}

void Tracker::SetPostTracker1(Cluster *cluster) {
   PostTracker1.set(cluster->getX(), cluster->getY(), -2, cluster->getSize());
}

void Tracker::SetPostTracker2(Cluster *cluster) {
   PostTracker2.set(cluster->getX(), cluster->getY(), -1, cluster->getSize());
}

Float_t Tracker::GetAngle(Cluster *c1, Cluster *c2) {
   // in degrees

   Float_t straightLength = 
            sqrt(pow(c2->getXmm() - c1->getXmm(), 2) +
            pow(c2->getYmm() - c1->getYmm(), 2) +
            pow(c2->getLayermm() - c1->getLayermm(), 2));

   Float_t xyDist =
            sqrt(pow(c2->getXmm() - c1->getXmm(), 2) +
            pow(c2->getYmm() - c1->getYmm(), 2));

   Float_t angle = atan2(xyDist, straightLength) * 180 / 3.14159265;

   return angle;
}

Float_t Tracker::GetPreAngle() {
   return GetAngle(&PreTracker1, &PreTracker2);
}

Float_t Tracker::GetPostAngle() {
   return GetAngle(&PostTracker1, &PostTracker2);
}

Cluster Tracker::GetPreDerivative() {
	// in geometrical coordinates
   Float_t dx = PreTracker2.getXmm() - PreTracker1.getXmm();
	Float_t dy = PreTracker2.getYmm() - PreTracker1.getYmm();
   Float_t dz = PreTracker2.getLayermm() - PreTracker1.getLayermm();

   return Cluster(dx / dz, dy / dz);
}

Cluster Tracker::GetPostDerivative() {
	// in geometrical coordinates
   Float_t dx = PostTracker2.getXmm() - PostTracker1.getXmm();
   Float_t dy = PostTracker2.getYmm() - PostTracker1.getYmm();
	Float_t dz = PostTracker2.getLayermm() - PostTracker1.getLayermm();

   return Cluster(dx / dz, dy / dz);
}

ostream& operator<< (ostream &os, Tracker &t) {
   os << "Pre-tracker 1: (" << t.PreTracker1.getX()
      << "," << t.PreTracker1.getY() << ")\n"
      << "Pre-tracker 2: (" << t.PreTracker2.getX()
      << "," << t.PreTracker2.getY() << ")\n"
      << "Post-tracker 1: (" << t.PostTracker1.getX()
      << "," << t.PostTracker1.getY() << ")\n"
      << "Post-tracker 2: (" << t.PostTracker2.getX()
      << "," << t.PostTracker2.getY() << ")\n";
   return os;
}
