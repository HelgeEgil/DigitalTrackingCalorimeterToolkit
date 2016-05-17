#include "Track.h"
#include "Cluster.h"
#include "Hit.h"
#include "Tools.h"
#include <iostream>
#include <cmath>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "Constants.h"
#include "MaterialConstants.h"
#include <TClonesArray.h>
#include <TF1.h>
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

   if (kEventsPerRun == 1) {
   	c->setEventID(copyCluster->getEventID());
   }
}

void Track::appendPoint(Float_t x, Float_t y, Int_t layer, Int_t size, Int_t eventID) {
   Int_t i = GetEntriesFast();
   Cluster *c = (Cluster*) track_.ConstructedAt(i);
   c->set(x,y,layer,size);
   c->setEventID(eventID);
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

Float_t Track::getRangemm() {
	// return range (z2 - z1), not track length
	Float_t range = 0;
	if (!At(0)) return 0;

	range = Last()->getLayermm() - At(0)->getLayermm();

	return range;
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

Float_t Track::getProjectedRatio() {
	Float_t angle = getSlopeAngle() * kRad;
	
	Float_t actualLength = getTrackLengthmm();
	Float_t straightLength = actualLength / getSinuosity();
	Float_t projectedLength = straightLength * cos(angle);
	
	return actualLength / projectedLength;
}

Float_t Track::getSlopeAngle() {
	Int_t lastIdx = GetEntriesFast() - 1;
	return getSlopeAngleAtLayer(lastIdx);
}



Float_t Track::getSlopeAngleAtLayer(Int_t i) {
	// return the slope angle of the track between
	// entry point in layer 0 and at layer i

	if (!At(i) || !At(0)) {
		cout << "Could not find slope angle...\n";
		return 0;
	}

	Cluster *a = At(0);
	Cluster *b = At(i);

	Float_t diffx = b->getXmm() - a->getXmm();
	Float_t diffy = b->getYmm() - a->getYmm();

	Float_t straightLength = diffmm(a,b);
	Float_t xyDist = sqrt(diffx*diffx + diffy*diffy);

	Float_t angle = atan2(xyDist, straightLength) / kRad;
	return angle;
}

Float_t Track::getSlopeAngleBetweenLayers(Int_t i) {
	// return slope angle between i and i-1 (if i=0, return 0)
	if (i == 0) return 0;
	
	Cluster *a = At(i-1);
	Cluster *b = At(i);
	
	Float_t diffx = b->getXmm() - a->getXmm();
	Float_t diffy = b->getYmm() - a->getYmm();
	
	Float_t straightLength = diffmm(a,b);
	Float_t xyDist = sqrt(diffx*diffx + diffy*diffy);

	Float_t angle = atan2(xyDist, straightLength) / kRad;
	return angle;
}	

Float_t Track::getAbsorberLength(Int_t i) {
	// returns the traversed absorber length before sensor i
	
	if (i == 0) {
		return dz/2; // due to the absorber being one-sided
	}
	
	Float_t angle = getSlopeAngleBetweenLayers(i) * kRad;
	Float_t distanceBetweenSensorLayers = dz;
	
	if (getYmm(i-1) > 0 && getYmm(i) < 0) {
		distanceBetweenSensorLayers += firstLowerLayerZ - firstUpperLayerZ;
	}

	else if (getYmm(i-1) < 0 && getYmm(i) > 0) {
		distanceBetweenSensorLayers += firstUpperLayerZ - firstLowerLayerZ;
	}
	
	Float_t absorberLength = distanceBetweenSensorLayers / cos(angle);
	
	if (kOutputUnit == kWEPL) {
		absorberLength = getWEPLFromTL(absorberLength);
	}
	
	return absorberLength;
}

Float_t Track::getTrackLengthmmAt(Int_t i) {
	// returns geometrical track length at point i
	
	if (i==0) return 0;
	if (i>GetEntriesFast()) return 0;

	return diffmm(At(i-1), At(i));
}

Float_t Track::getRangemmAt(Int_t i) {
	if (i==0) return 0;
	if (i>GetEntriesFast()) return 0;

	return diffmm(new Cluster(getX(i-1), getY(i-1), getLayer(i-1)),
				   new Cluster(getX(i-1), getY(i-1), getLayer(i)));
}

Float_t Track::getWEPL() {
	Float_t tl = getTrackLengthmm();
	if (!tl) return 0;
	
	return getWEPLFromTL(tl);
}

Float_t Track::getTrackLengthWEPLmmAt(Int_t i) {
	if (i==0) {
		return 0; 
	}

	Float_t tl = getTrackLengthmmAt(i);
	if (!tl) return 0;
	
	return getWEPLFromTL(tl);
}

Float_t Track::getRangeWEPLAt(Int_t i) {	
	Float_t range = getRangemmAt(i);
	if (i==0) {
		return 0;
	}

	if (!range) return 0;

	return getWEPLFromTL(range);
}

Float_t Track::getEnergyStraggling() {
	Float_t energy = getEnergy();
	if (!energy) return 0;
	
	Float_t tl = getTLFromEnergy(energy);
	if (!tl) return 0;
	
	return getEnergyStragglingFromTL(tl, 0);
}

Float_t Track::getEnergy() {
	Float_t energy = 0;
	
// 	Float_t tl = getTrackLengthmm();
	Float_t range = getRangemm();
	if (range>0) energy = getEnergyFromTL(range);
	
	return energy;
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

      if (kEventsPerRun == 1) {
      	newStart->setEventID(At(1)->getEventID());
      }

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

Cluster * Track::getExtrapolatedClusterAt(Float_t mmBeforeDetector) {
	Int_t firstLayer = getFirstLayer();
	Int_t firstIdx = getClusterFromLayer(firstLayer);
	Int_t nextIdx = 0;

	Float_t extra_mm = getLayermm(firstLayer);

	for (Int_t i=firstIdx+1; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;

		nextIdx = i;
		break;
	}
	
	Float_t softenFactor = 1;

	Cluster *first = At(firstIdx);
	Cluster *next = At(nextIdx);

	Float_t diff_z = (next->getLayermm() - first->getLayermm());
	Float_t diff_x = (first->getX() - next->getX()) / diff_z;
	Float_t diff_y = (first->getY() - next->getY()) / diff_z;

	Float_t new_x = first->getX() + softenFactor * diff_x * (mmBeforeDetector + extra_mm);
	Float_t new_y = first->getY() + softenFactor * diff_y * (mmBeforeDetector + extra_mm);

	Cluster *extrapolated = new Cluster(new_x, new_y);	
	
	return extrapolated;
}

Float_t Track::getAverageCS() {
	Int_t n = GetEntries();
	if (!n) return 0;

	Float_t sum=0, avg=0;

	for (Int_t i=0; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;

		sum += getSize(i);
	}
	avg = sum / n;

	return avg;
}


Float_t Track::getAverageCSLastN(Int_t last_n) {
	Int_t n=0;

	Float_t sum=0;
	Float_t avg=0;

	// find last index
	Int_t nFound = 0;
	Int_t lastIdx = 0;
	for (Int_t i=GetEntriesFast()-1; i>=0; i--) {
		if (!At(i)) continue;
		nFound++;
		if (nFound == last_n) {
			lastIdx = i;
			break;
		}
	}

	for (Int_t i=lastIdx; i<GetEntriesFast(); i++) {
		if (!At(i)) continue;
		n++;
		sum += getSize(i);
	}

	avg = sum / n;
	return avg;
}

TGraphErrors * Track::doRangeFit() {
	Bool_t newCutBraggPeak = (getAverageCSLastN(2) > getAverageCS()*kBPFactorAboveAverage);
	Bool_t cutNPointsInTrack = (GetEntries()>3);
	Bool_t cut = newCutBraggPeak * cutNPointsInTrack;
	if (!cut) return 0;

	Int_t n = GetEntriesFast();
	Float_t x[n], y[n];
	Float_t erx[n], ery[n];
	Float_t preTL = getPreTL();
	
	for (Int_t i=0; i<n; i++) {
		if (!At(i)) continue;
		x[i] = preTL + getLayermm(i);
		y[i] = getDepositedEnergy(i);
		ery[i] = getDepositedEnergyError(i);
		erx[i] = dz / sqrt(12);
	}
	
	Float_t max_energy = getEnergyFromTL(x[n-1] + 0.75*dz);
	Float_t estimated_energy = getEnergyFromTL(x[n-1] + 0.5*dz);

	if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		Float_t WEPLFactor = getWEPLFactorFromEnergy(estimated_energy);
		for (Int_t i=0; i<n; i++) {
			x[i] = x[i] * WEPLFactor;
			erx[i] = erx[i] * WEPLFactor;
		}
	}
	
	TGraphErrors *graph = new TGraphErrors(n, x, y, erx, ery);
	Float_t scaleParameter = 0;
	
	if (kOutputUnit == kPhysical) {
		if (kMaterial == kTungsten) scaleParameter = 14;
		if (kMaterial == kAluminum) scaleParameter = 65;
	}
	
	else if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		if (kMaterial == kTungsten) scaleParameter = 100;
		if (kMaterial == kAluminum) scaleParameter = 126;
	}

	TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, getWEPLFromEnergy(max_energy*1.2), 2);
	func->SetParameter(0, estimated_energy);
	func->SetParameter(1, scaleParameter);
	func->SetParLimits(0, 0, max_energy);
	
	func->SetNpx(500);
	
	func->SetParLimits(1, scaleParameter,scaleParameter);
	graph->Fit("fit_BP", "B, Q, N, W", "", 0, getWEPLFromEnergy(max_energy*1.2));
	
	fitEnergy_ = correctForEnergyParameterisation(func->GetParameter(0));
	fitScale_ = func->GetParameter(1);
	fitError_ = func->GetParError(0);

	return graph;
}

TGraphErrors * Track::doFit() {
	Bool_t newCutBraggPeak = (getAverageCSLastN(2) > getAverageCS()*kBPFactorAboveAverage);
	Bool_t cutNPointsInTrack = (GetEntries()>3);
	Bool_t cut = newCutBraggPeak * cutNPointsInTrack;
	if (!cut) return 0;

	Int_t n = GetEntriesFast();
	Float_t x[n], y[n];
	Float_t erx[n], ery[n];
//	Float_t preTL = getPreTLFromScintillatorAndAluminum();
	Float_t preTL = getPreTL();
	Float_t trackLength = preTL;
	
	for (Int_t i=0; i<n; i++) {
		if (!At(i)) continue;
		trackLength += getTrackLengthmmAt(i);
		x[i] = trackLength;
		y[i] = getDepositedEnergy(i);
		ery[i] = getDepositedEnergyError(i);
		erx[i] = getAbsorberLength(i);
	}
	
	Float_t max_energy = getEnergyFromTL(trackLength + 0.75*dz);
	Float_t estimated_energy = getEnergyFromTL(trackLength + 0.5*dz);

	if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		Float_t WEPLFactor = getWEPLFactorFromEnergy(estimated_energy);
		for (Int_t i=0; i<n; i++) {
			x[i] = x[i] * WEPLFactor;
		}
	}
	
	TGraphErrors *graph = new TGraphErrors(n, x, y, erx, ery);
	Float_t scaleParameter = 0;
	
	if (kOutputUnit == kPhysical) {
		if (kMaterial == kTungsten) scaleParameter = 14;
		if (kMaterial == kAluminum) scaleParameter = 65;
	}
	
	else if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		if (kMaterial == kTungsten) scaleParameter = 100;
		if (kMaterial == kAluminum) scaleParameter = 126;
	}
		
	TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
	func->SetParameter(0, estimated_energy);
	func->SetParameter(1, scaleParameter);
	func->SetParLimits(0, 0, max_energy);
	func->SetParLimits(1, scaleParameter,scaleParameter);
	graph->Fit("fit_BP", "B, Q, N, WW", "", 0, 500);
	
	fitEnergy_ = func->GetParameter(0);
	fitScale_ = func->GetParameter(1);

	return graph;
}

Float_t Track::getFitParameterEnergy() {
	if (!fitEnergy_) {
		if (!run_energy) {
			return 0;
		}
		else {
			doFit();
		}
	}

	return fitEnergy_;
}

Float_t Track::getFitParameterScale() {
	if (!fitScale_) {
		if (!run_energy) {
			return 0;
		}
		else {
			doFit();
		}
	}
	return fitScale_;
}

Float_t Track::getFitParameterError() {
	if (!fitScale_) {
		if (!run_energy) {
			return 0;
		}
		else {
			doFit();
		}
	}
	return fitError_;
}

Float_t Track::getPreEnergyLoss() {
	Float_t energyLoss = 0;
	Int_t nScintillators = 0;

	if (kIsScintillator) {
		nScintillators = 1 + isHitOnScintillatorH() + isHitOnScintillatorV();
 		energyLoss = getEnergyLossFromScintillators(run_energy, nScintillators);
	}

	if (kIsAluminumPlate) {
		energyLoss += getEnergyLossFromAluminumAbsorber(run_energy);
	}
	
	return energyLoss;
}

Float_t Track::getPreWEPL() {
	Float_t energyLoss = 0;
	Int_t nScintillators = 0;

	if (kIsScintillator) {
		nScintillators = 1 + isHitOnScintillatorH() + isHitOnScintillatorV();
		energyLoss = getEnergyLossFromScintillators(run_energy, nScintillators);
	}

	if (kIsAluminumPlate) {
		energyLoss += getEnergyLossFromAluminumAbsorber(run_energy);
	}
	
	Float_t wepl = getWEPLFromEnergy(run_energy) - getWEPLFromEnergy(run_energy - energyLoss);
	return wepl;
}

Float_t Track::getPreTL() {
	Float_t energyLoss = 0;
	Int_t nScintillators = 0;

	if (kIsScintillator) {
		nScintillators = getNScintillators();
		energyLoss = getEnergyLossFromScintillators(run_energy, nScintillators);
	}

	if (kIsAluminumPlate) {
		energyLoss += getEnergyLossFromAluminumAbsorber(run_energy);
	}
	
	Float_t tl = getTLFromEnergy(run_energy) - getTLFromEnergy(run_energy - energyLoss);
	return tl;
}

Int_t Track::getNScintillators() {
	return 1 + isHitOnScintillatorH() + isHitOnScintillatorV();
}

Float_t Track::getPreEnergyLossError() {
	Float_t energyLossError = 0;

	Bool_t scintillatorH = isHitOnScintillatorH();
	Bool_t scintillatorV = isHitOnScintillatorV();
	Bool_t scintillatorF = true; // 4x4 cm scintillator, full field

	Int_t nScintillators = scintillatorF + scintillatorH + scintillatorV;

	energyLossError = quadratureAdd(getEnergyLossErrorFromScintillators(nScintillators),
											  getEnergyLossErrorFromAluminumAbsorber());

	return energyLossError;
}

Bool_t Track::isHitOnScintillators() {
	return (isHitOnScintillatorH() || isHitOnScintillatorV());
}

Bool_t Track::isHitOnScintillatorH() {
	// (dx, dy, dz) = (4, 1, 0.5) cm @ -174 mm
	Float_t scintillatorMeanDepth = 174;
	Bool_t isOnScintillator = false;

	Cluster *extrapolatedCluster = getExtrapolatedClusterAt(0);
	Float_t y = extrapolatedCluster->getYmm();

	if (fabs(y)<5) {
		isOnScintillator = true;
	}

	return isOnScintillator;
}

Bool_t Track::isHitOnScintillatorV() {
	// (dx, dy, dz) = (1, 4, 0.5) cm @ -166 mm
	Float_t scintillatorMeanDepth = 166;
	Bool_t isOnScintillator = false;

	Cluster *extrapolatedCluster = getExtrapolatedClusterAt(0);
	Float_t x = extrapolatedCluster->getXmm();

	if (fabs(x)<5) {
		isOnScintillator = true;
	}

	return isOnScintillator;
}

Bool_t Track::doesTrackEndAbruptly() {
	Float_t expectedTL = getTLFromEnergy(run_energy);
	Float_t actualTL = getTLFromEnergy(getEnergy());
	
	Float_t riseFactor = getRiseFactor();
	
	if (actualTL < expectedTL * 0.6 && riseFactor < 1.5) {
// 		cout << "Track does end abruptly! Expected (actual) TL = " << expectedTL << " (" << actualTL << ", " << getEnergy() << " MeV), riseFactor = " << riseFactor << endl;
		return true;
	}
	
	else return false;
}

Float_t Track::getRiseFactor() {
	Float_t normalization = 0, riseFactor = 0;
	Int_t nNorm = 0, nNormActual = 0;
	Int_t n = GetEntriesFast();
	
	if (n==0) return 0;
	
	nNorm = n/2;
	
	for (Int_t i=0; i<nNorm; i++) {
		if (!At(i)) continue;
		normalization += getSize(i);
		nNormActual++;
	}
	
	normalization /= nNormActual;
	
	if (At(n-1) && At(n-2)) {
		riseFactor = (getSize(n-1) + getSize(n-2))/2 / normalization;
	}
	
	else if (At(n-1) && !At(n-2)) {
		riseFactor = getSize(n-1) / normalization;
	}
	
	else if (!At(n-1) && At(n-2)) {
		riseFactor = getSize(n-2) / normalization;
	}
	
	else return 0;
	
	return riseFactor;
}
	
