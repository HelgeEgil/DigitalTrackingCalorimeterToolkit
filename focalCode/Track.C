#include "Track.h"
#include "Cluster.h"
#include "Hit.h"
#include <iostream>
#include <cmath>
#include "Constants.h"
#include "MaterialConstants.h"
#include <TClonesArray.h>
// #include <TObject.h>

using namespace std;

Track::~Track() {
   // Destructor
//   ClearTrack();
	track_.Delete();
}

void Track::Clear(Option_t *) {
	track_.Clear("C");
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
   Int_t i = GetEntriesFast();
   if (i==0) i += startOffset;

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

   if (!At(0)) return 0;

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
	
	if (i==0) return 3.15; // pre-detector material
	if (i>GetEntriesFast()) return 0;

	return diffmm(At(i-1), At(i));
}

Float_t Track::getWEPL() {
	Float_t tl = getTrackLengthmm();
	if (!tl) return 0;
	return getWEPLFromTL(tl) + 3.09 + 1.65 * getWEPLFactorFromEnergy(run_energy);
}

Float_t Track::getTrackLengthWEPLmmAt(Int_t i) {
	// Returns geometrical track length in WEPL at point i
	// the WEPL conversion is a LUT calculated from the bethe equation
	// (see rangeCalculator.py)
	
	Float_t tl = getTrackLengthmmAt(i);
	
	if (i==0) {
		// WEPL of 1.5 mm Al + WEPL of 1.65 (avg) detector materal
		return 3.09 + 1.65 * getWEPLFactorFromEnergy(run_energy);
	}
	
	if (!tl) return 0;
	
	return getWEPLFromTL(tl);
}

Float_t Track::getWEPLFromTL(Float_t tl) {

	for (Int_t i=0; i<nPLEnergies; i++) {
		if (kPLFocal[i] < tl) continue;
		else {
			Float_t ratio = (tl - kPLFocal[i-1]) / (kPLFocal[i] - kPLFocal[i-1]);
			return tl * (kWEPLRatio[i-1] + ratio * (kWEPLRatio[i] - kWEPLRatio[i-1]));
		}
	}
}

Float_t Track::getEnergyFromTL(Float_t tl) {
	// TODO: optimize via stackoverflow.com/questions/11396860
	// (and create bigger table for tungsten....)

	for (Int_t i=0; i<nPLEnergies; i++) {
		if (kPLFocal[i] < tl) continue;
		else {
			Float_t ratio = (tl - kPLFocal[i-1]) / (kPLFocal[i] - kPLFocal[i-1]);
			return kPLEnergies[i-1] + ratio * (kPLEnergies[i] - kPLEnergies[i-1]);
		}
	}
	
}

Float_t Track::getEnergyFromWEPL(Float_t wepl) {

	for (Int_t i=0; i<nPLEnergies; i++) {
		if (kPLFocal[i] * kWEPLRatio[i] < wepl) continue;
		else {
			Float_t ratio = (wepl - kPLFocal[i-1] * kWEPLRatio[i-1]) / (kPLFocal[i]* kWEPLRatio[i] - kPLFocal[i-1]* kWEPLRatio[i-1]);
			return kPLEnergies[i-1] + ratio * (kPLEnergies[i] - kPLEnergies[i-1]);
		}
	}
}

Float_t Track::getWEPLFactorFromEnergy(Float_t energy) {
	for (Int_t i=0; i<nPLEnergies; i++) {
		if (kPLEnergies[i] < energy) { continue; }
		else { return kWEPLRatio[i]; }
	}
}

Float_t Track::getEnergy() {
	Float_t tl = getTrackLengthmm();
	if (!tl) return 0;
	return getEnergyFromTL(tl);
}

Float_t Track::getSnakeness() {
	// How much the track bends
	Int_t n = GetEntriesFast();
	Float_t snakeness = 0;

	if (n<3 || !At(0) || !At(1)) snakeness = 3;
	else {
		Float_t x, y, xp, yp, dXY, dX, dY;

		xp = At(1)->getXmm(); yp = At(1)->getYmm();

		Float_t slopeX = xp - At(0)->getXmm();
		Float_t slopeY = yp - At(0)->getYmm();

		for (Int_t i=2; i<n; i++) {
			if (!At(i)) continue;
			x = At(i)->getXmm(); y = At(i)->getYmm();

			dX = xp + slopeX - x;
			dY = yp + slopeY - y;
			dXY = pow(dX*dX + dY*dY, 0.5);

			slopeX = x - xp; slopeY = y - yp;

			snakeness += dXY;
			xp = x; yp = y;
		}
	}
	return snakeness;
}

Float_t Track::getTrackScore() {
	Float_t upperTrackLength = 35;
	Float_t upperSnakeness = 3;
	Int_t snakePoints = 5;
	Int_t trackLengthPoints = 25;

	Float_t trackLength = getTrackLengthmm();
	Float_t snakeness = getSnakeness();

	if (trackLength == 0) return 0;

	Float_t points = trackLength * (trackLengthPoints / upperTrackLength)
						+ (upperSnakeness - snakeness) * snakePoints/upperSnakeness;

	return points;
}

Int_t Track::getNMissingLayers() {
   Int_t missingLayers = 0;
   Int_t lastLayer = 0, firstLayer = 1e6;
   
   for (Int_t i=0; i<GetEntriesFast(); i++) {
     if (At(i)) {
       firstLayer = i;
       break;
     }
   }
   if (firstLayer > 1e5) cout << "Couldn't set first layer!\n";
   
   for (Int_t i=GetEntriesFast()-1; i>=0; i--) {
     if (At(i)) {
       lastLayer = i;
       break;
     }   
  }
  if (lastLayer == 0) cout << "Couldn't set last layer!\n";
  
  Int_t previousLayer = firstLayer - 1;
  Int_t diff = 0;
  
  for (Int_t i=firstLayer; i<=lastLayer; i++) {
    if (!At(i)) continue;
    
    diff = getLayer(i) - previousLayer - 1;
    if (diff>0) missingLayers += diff;
    previousLayer = getLayer(i);
  }
  
  return missingLayers;
}

Float_t Track::getMeanSizeToIdx(Int_t toIdx) {
	Int_t tempSum = 0;

	if (toIdx<0 || toIdx >= GetEntriesFast())
		toIdx = GetEntriesFast() - 1;

	for (Int_t i=0; i<toIdx; i++) {
		tempSum += getSize(i);
	}

	return (float) tempSum / (toIdx + 1 - getNMissingLayers());
}

Float_t Track::getStdSizeToIdx(Int_t toIdx) {
	Float_t tempSum = 0;

	if (toIdx<0 || toIdx >= GetEntriesFast())
		toIdx = GetEntriesFast() - 1;

	Int_t mean = getMeanSizeToIdx(toIdx);

	for (Int_t i=0; i<toIdx; i++) {
		tempSum += pow(getSize(i) - mean, 2);
	}

	Float_t std = sqrt((float) tempSum / (toIdx + 1));
	return std;
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

Int_t Track::getClusterFromLayer(Int_t layer) {
	// slow implementation first
	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;

		if (getLayer(i) == layer)
			return i;
		else if (getLayer(i) > layer) {
			cout << "GOT BREAK COMMAND. COULD NOT FIND LAYER " << layer << ": ";
			for (Int_t j=0; j<GetEntriesFast(); j++) cout << getLayer(j) << ", ";
			cout << endl;
		}
	}
	return -1;
}

Bool_t Track::hasLayer(Int_t layer) {
	for (Int_t i = 0; i < GetEntriesFast(); i++) {
		if (!At(i)) continue;
		if (getLayer(i) == layer) return true;
	}
	return false;
}

Int_t Track::getLastLayer() {
	for (Int_t i = GetEntriesFast() - 1; i >= 0; i--) {
		if (!At(i)) continue;
		return getLayer(i);
	}
	return -1;
}

Int_t Track::getFirstLayer() { 
	for (Int_t i = 0; i < GetEntriesFast(); i++) {
		if (!At(i)) continue;
		return getLayer(i);
	}
	return -1;
}

Cluster * Track::getInterpolatedClusterAt(Int_t layer) {
	Cluster *pre = 0;
	Cluster *post = 0;
	Int_t lastIdx = 0;
	
	Int_t firstLayer = getFirstLayer();
	Int_t lastLayer = getLastLayer();
	Int_t firstIdx = getClusterFromLayer(firstLayer);
	
	if (firstLayer == layer) return 0;
	if (lastLayer <= layer) return 0;
	
	for (Int_t i=firstIdx; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;

		if (getLayer(i) > layer) {
			post = At(i);
			lastIdx = i;
			break;
		}
	}

	for (Int_t i=lastIdx; i>=firstIdx; i--) {
		if (!At(i)) continue;

		if (getLayer(i) < layer) {
			pre = At(i);
			break;
		}
	}

	Float_t x = ( pre->getX() + post->getX() ) / 2.;
	Float_t y = ( pre->getY() + post->getY() ) / 2.;

	return new Cluster(x, y, layer);
}
