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

Cluster * Track::getInterpolatedClusterAt(Int_t layer) {
	Cluster   *	pre = 0;
	Cluster   *	post = 0;
	Int_t			lastIdx = 0;
	Int_t			firstLayer = getFirstLayer();
	Int_t			lastLayer = getLastLayer();
	Int_t			firstIdx = getClusterFromLayer(firstLayer);
	Float_t		x, y;
	
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

	x = ( pre->getX() + post->getX() ) / 2.;
	y = ( pre->getY() + post->getY() ) / 2.;

	return new Cluster(x, y, layer);
}

Cluster * Track::getExtrapolatedClusterAt(Float_t mmBeforeDetector) {
	Int_t		firstLayer = getFirstLayer();
	Int_t		firstIdx = getClusterFromLayer(firstLayer);
	Int_t		nextIdx = 0;
	Float_t	diffx, diffy, diffz;
	Float_t	extra_mm = getLayermm(firstLayer);

	for (Int_t i=firstIdx+1; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;

		nextIdx = i;
		break;
	}
	
	Cluster *first = At(firstIdx);
	Cluster *next = At(nextIdx);

	diffz = (next->getLayermm() - first->getLayermm());
	diffx = (first->getX() - next->getX()) / diffz;
	diffy = (first->getY() - next->getY()) / diffz;

	Float_t new_x = first->getX() + diffx * (mmBeforeDetector + extra_mm);
	Float_t new_y = first->getY() + diffy * (mmBeforeDetector + extra_mm);

	Cluster *extrapolated = new Cluster(new_x, new_y);	
	
	return extrapolated;
}

void Track::extrapolateToLayer0() {
   Int_t nTracks = GetEntriesFast();
   Int_t eventID = -1;

   if (!At(0)) {
      // OK, no information in first layer

      // now GetLayer(1) should be 1, and GetLayer(2) should be 2
      // If both layer 0 and 2 are skipped, some more extrapolation should be done
      // maybe set new x as [1.x - slope.x * (2.z - 1.z)]

      Hit *slope = new Hit();

      if (!At(1)) {
         cout << "No pointer for At(1) as well... Aborting this track extrapolation.\n";
         return;
      }
		else {
			eventID = At(1)->getEventID();
			if (getLayer(1) == 0) return; // hotfix... Why is layer 0 sometimes placed in idx 1? Must fix this.
		}

      if (At(2)) {
         slope->set(getX(2) - getX(1), getY(2) - getY(1), getLayer(2) - getLayer(1));
         if (eventID<0) eventID = At(2)->getEventID();
      }
      else if (!At(2) && At(3)) {
         slope->set(getX(3) - getX(1), getY(3) - getY(1), getLayer(3) - getLayer(1));
         if (eventID<0) eventID = At(3)->getEventID();
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

		newStart->setEventID(eventID);
   }
}

