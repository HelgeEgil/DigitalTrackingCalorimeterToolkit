#ifndef trackAngleCalculations_cxx
#define trackAngleCalculations_cxx

#include <iostream>
#include <cmath>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TClonesArray.h>
#include <TF1.h>

#include "Classes/Track/Track.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Hit/Hit.h"
#include "HelperFunctions/Tools.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"

Float_t Track::getSlopeAngleAtLayer(Int_t i) {
   // Returns slope angle between entry point and point i
   // In degrees
   
   Cluster *a, *b;
   Float_t  diffx, diffy, diffz;
   Float_t  xyDist, angle;

   if (!At(i) || !At(0)) {
      return 0;
   }

   a = At(0); b = At(i);

   diffx = b->getXmm() - a->getXmm();
   diffy = b->getYmm() - a->getYmm();
   diffz = b->getLayermm() - a->getLayermm();
   
   xyDist = quadratureAdd(diffx, diffy);

   angle = atan2(xyDist, diffz) / kRad;

   return angle;
}

Float_t Track::getSlopeAngleBetweenLayers(Int_t i) {
   // return slope angle between i and i-1 (if i=0, return 0)
   // In degrees
   
   if (i == 0) return 0;

   Cluster *a, *b;
   Float_t diffx, diffy, diffz;
   Float_t xyDist, angle;

   a = At(i-1); b = At(i);

   if (!a || !b) return -1;
   
   diffx = b->getXmm() - a->getXmm();
   diffy = b->getYmm() - a->getYmm();
   diffz = b->getLayermm() - a->getLayermm();

   xyDist = quadratureAdd(diffx, diffy);

   angle = atan2(xyDist, diffz) / kRad;

   return angle;
}  

Float_t Track::getSlopeAngle() {
   // Returns slope angle between entry point and last layer
   // In degrees

   Int_t lastIdx = GetEntriesFast() - 1;
   return getSlopeAngleAtLayer(lastIdx);
}
Float_t Track::getSlopeAngleChangeBetweenLayers(Int_t i) {
   // The change in angle due to material between sensor layer i and i+1
   // In degrees

   if ((i+1) >= GetEntriesFast()) return 0;
   
   Float_t theta1 = 0, theta2 = 0;

   if (i>0) { theta1 = getSlopeAngleBetweenLayers(i); }
   theta2 = getSlopeAngleBetweenLayers(i+1);

   return theta2  - theta1;
}

Float_t Track::getSlopeAngleDifferenceSum() {
   // Returns the total slope angle difference
   // In degrees

   Float_t angleSum = 0;

   for (Int_t i=1; i<GetEntriesFast() - 1; i++) {
      angleSum += fabs(getSlopeAngleChangeBetweenLayers(i));
   }

   return angleSum;
}

Float_t Track::getSlopeAngleDifferenceSumInTheta0() {
   // Returns the total slope angle difference ROOT SUM SQUARE
   // In units of theta0 (expected multiple coulomb scattering)

   Float_t theta0sum = 0;
   Float_t dtheta = 0;
   Float_t expected_mcs = 0;

   for (Int_t i=1; i<GetEntriesFast() - 1; i++) {
      if (!At(i)) continue;
      dtheta = fabs(getSlopeAngleChangeBetweenLayers(i));
      expected_mcs = getMCSAngleForLayer(getLayer(i));

      theta0sum += pow(dtheta / expected_mcs, 2);
   }

   return sqrt(theta0sum);
}

Float_t Track::getSinuosity() {
   // The sinuosity is the total track length divided by the
   // straight line from start to end

   Cluster *a, *b;
   Float_t straightLength, actualLength;

   a = At(0); b = Last();

   if (a == b) return 0;

   straightLength = diffmmXYZ(a,b); 
   actualLength = getTrackLengthmm();

   if (straightLength == 0) return 0;

   return actualLength / straightLength;
}

Float_t Track::getProjectedRatio() {
   Float_t angle = getSlopeAngle() * kRad;
   
   Float_t actualLength = getTrackLengthmm();
   Float_t straightLength = actualLength / getSinuosity();
   Float_t projectedLength = straightLength * cos(angle);
   
   return actualLength / projectedLength;
}

Float_t Track::getMaximumSlopeAngleChange() {
   Float_t maxChange = 0;

   for (Int_t i=1; i<GetEntriesFast() - 1; i++) {
      maxChange = max(maxChange, fabs(getSlopeAngleChangeBetweenLayers(i)));
   }

   return maxChange;
}

Float_t Track::getAbsorberLength(Int_t i) {
   // returns the geometric absorber length between sensor i-1 and sensor i
   // If the (y<0) changes between two layers, add/remove a small difference

   if (i == 0) return dz/2;
   Float_t angle, distanceBetweenSensorLayers, absorberLength;
   
   angle = getSlopeAngleBetweenLayers(i) * kRad;
   distanceBetweenSensorLayers = dz;
   
   if (getYmm(i-1) > 0 && getYmm(i) < 0) {
      distanceBetweenSensorLayers += firstLowerLayerZ - firstUpperLayerZ;
   }

   else if (getYmm(i-1) < 0 && getYmm(i) > 0) {
      distanceBetweenSensorLayers += firstUpperLayerZ - firstLowerLayerZ;
   }
   
   absorberLength = distanceBetweenSensorLayers / cos(angle);
   
   if (kOutputUnit == kWEPL) {
      absorberLength = getWEPLFromTL(absorberLength);
   }
   
   return absorberLength;
}

#endif
