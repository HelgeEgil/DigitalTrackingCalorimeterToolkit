#include "Track.h"
#include "Cluster.h"
#include "Hit.h"
#include <iostream>
#include <cmath>
#include "Constants.h"
#include <TClonesArray.h>
// #include <TObject.h>

Track::~Track() {
   // Destructor
   clearTrack();
}

Float_t Track::diffmm(Cluster *p1, Cluster *p2) {
	return sqrt(pow(p2->getXmm() - p1->getXmm(), 2) +
					pow(p2->getYmm() - p1->getYmm(), 2) +
					pow(p2->getLayermm() - p1->getLayermm(), 2));
}


void Track::setTrack(Track *copyTrack, Int_t startOffset /* default 0 */) {
   clearTrack();
	for (Int_t i=0; i<copyTrack->GetEntriesFast(); i++) {
      appendCluster(copyTrack->At(i), startOffset);
   }
}

void Track::appendCluster(Cluster *copyCluster, Int_t startOffset /* default 0 */) {
   Int_t i = GetEntriesFast() + startOffset;

   // copy PROPERTIES from copycluster onto track
   // new object, no pointers are copied

   if (!copyCluster) return;

   Cluster *c = (Cluster*) track_.ConstructedAt(i);
   c->set(copyCluster->getX(), copyCluster->getY(), 
          copyCluster->getLayer(), copyCluster->getSize());
}

void Track::appendPoint(Float_t x, Float_t y, Int_t layer, Int_t size) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) track_.ConstructedAt(i);
   c->set(x,y,layer,size);
}

Float_t Track::getTrackLengthmm() {
   // return geometrical track length

   Float_t trackLength = 0;
   Float_t x = -1; Float_t xp = -1;
   Float_t y = -1; Float_t yp = -1;
   Float_t z = -1; Float_t zp = -1;

   xp = getXmm(0);
   yp = getYmm(0);
   zp = getLayermm(0);

   for (Int_t i=1; i<GetEntriesFast(); i++) {

		if (!At(i)) continue;

      x = getXmm(i);
      y = getYmm(i);
      z = getLayermm(i);
		
      trackLength += sqrt((xp-x)*(xp-x) + (yp-y)*(yp-y) + (zp-z)*(zp-z));

      xp = x; yp = y; zp = z;
   }
   return trackLength;
}

Float_t Track::getSinuosity() {
	Cluster *a = At(0);
	Cluster *b = Last();

	if (a == b) return 0;

	Float_t straightLength = diffmm(a,b);
	Float_t actualLength = getTrackLengthmm();

	if (straightLength == 0 || actualLength == 0) return 0;

	return actualLength / straightLength;
}

Float_t Track::getSlopeAngle() {
	Int_t lastIdx = GetEntriesFast() - 1;
	return getSlopeAngleAtLayer(lastIdx);
}

Float_t Track::getSlopeAngleAtLayer(Int_t i) {
	// return the slope angle of the track between
	// entry point in layer 0 and at layer i

	Cluster *a = At(0);
	Cluster *b = At(i);

	Float_t dx = b->getXmm() - a->getXmm();
	Float_t dy = b->getYmm() - a->getYmm();

	Float_t straightLength = diffmm(a,b);
	Float_t xyDist = sqrt(dx*dx + dy*dy);

	Float_t angle = atan2(xyDist, straightLength) * 180 / 3.14159265;
	return angle;
}

Float_t Track::getTrackLengthmmAt(Int_t i) {
	// returns geometrical track length at point i
	
	if (i==0) return 0;
	if (i>GetEntriesFast()) return 0;

	return diffmm(At(i-1), At(i));
}

Int_t Track::getNMissingLayers() {
   Int_t missingLayers = 0;
   Int_t lastLayer = getLayer(0);
   for (Int_t i=1; i<track_.GetEntriesFast(); i++) {
      if (getLayer(i) - lastLayer > 1) missingLayers++;
      lastLayer = getLayer(i);
   }
   return missingLayers;
}

void Track::extrapolateToLayer0() {
   Int_t nTracks = GetEntriesFast();
   if (!At(0)) {
      cout << "No pointer found for location (good)\n";
      // OK, no information in first layer
      // (stored as a 0x0 pointer)

      // now GetLayer(1) should be 1, and GetLayer(2) should be 2
      // If both layer 0 and 2 are skipped, some more extrapolation should be done
      // maybe set new x as [1.x - slope.x * (2.z - 1.z)]

      Hit *slope = new Hit();

      if (!At(1)) {
         cout << "No pointer for At(1) as well... Aborting this track extrapolation.\n";
         return;
      }
      if (At(2)) {
         slope->set(getX(2) - getX(1), getY(2) - getY(1), getLayer(2) - getLayer(1));
      }
      else if (!At(2) && At(3)) {
         slope->set(getX(3) - getX(1), getY(3) - getY(1), getLayer(3) - getLayer(1));
      }
      else {
         cout << "Too many holes (!At(0), At(1), !At(2), !At(3), .....), aborting this track extrapolation.\n";
         return;
      }

      // must create pointer!
      Cluster *newStart = (Cluster*) track_.ConstructedAt(0);

      newStart->set(getX(1) - slope->getLayer() * slope->getX(), // new x
                 getY(1) - slope->getLayer() * slope->getY(), // new y
                 0); // layer = 0

      cout << "Extrapolated to " << *At(0) << " from " << *At(1) << endl;
   }
}
