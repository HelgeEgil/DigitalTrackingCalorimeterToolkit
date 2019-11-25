#ifndef trackRangeCalculations_cxx
#define trackRangeCalculations_cxx

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

Float_t Track::getTrackLengthmm() {
   Float_t trackLength = 0;
   Float_t x = -1, xp = -1;
   Float_t y = -1, yp = -1;
   Float_t z = -1, zp = -1;

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

Float_t Track::getTrackLengthmmAt(Int_t i) {
   if (i==0)               return 0;
   if (i>GetEntriesFast()) return 0;

   return diffmmXYZ(At(i-1), At(i));
}

Float_t Track::getTrackLengthWEPLmmAt(Int_t i) {
   Float_t tl = getTrackLengthmmAt(i);
   if (!tl)    return 0;
   
   return getWEPLFromTL(tl);
}

Float_t Track::getRangemm() {
   Float_t range;
   if (!At(0)) return 0;

   range = Last()->getLayermm(); //  - getLayermm(0);

   return range;
}

Float_t Track::getRangemmAt(Int_t i) {
   // Changed this function, was Delta range at layer i!! Same as dz..?

   if (i==0) return 0;
   if (i>GetEntriesFast()) return 0;

   return getLayermm(i); //  - getLayermm(0);
}

Float_t Track::getWEPL() {
   Float_t tl = getTrackLengthmm();
   if (!tl) return 0;
   
   return getWEPLFromTL(tl);
}

Float_t Track::getRangeWEPLAt(Int_t i) {  
   Float_t range = getRangemmAt(i);
   if (!range) return 0;

   return getWEPLFromTL(range);
}

Float_t Track::getEnergy() {
   Float_t energy = 0;
   Float_t range = getRangemm();

   if (range>0) energy = getEnergyFromTL(range);
   
   return energy;
}

Float_t Track::getEnergyStraggling() {
   Float_t energy = getEnergy();
   if (!energy) return 0;
   
   Float_t tl = getTLFromEnergy(energy);
   if (!tl) return 0;
   
   return getEnergyStragglingFromTL(tl, 0);
}

#endif
