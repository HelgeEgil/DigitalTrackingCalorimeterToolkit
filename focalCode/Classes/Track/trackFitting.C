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

TGraphErrors * Track::doFit() {
	// Fit a bragg curve on the clusters in this track
	// See HelperFunctions/Tools.C::fitfunc_DBP for the BC function used
	
	Bool_t newCutBraggPeak = (getAverageCSLastN(2) > getAverageCS()*kBPFactorAboveAverage);
	Bool_t cutNPointsInTrack = (GetEntries()>3);
	Bool_t cut = newCutBraggPeak * cutNPointsInTrack;
	if (!cut) return 0;

	TGraphErrors * graph = nullptr;
	Int_t				n = GetEntriesFast();
	Float_t			x[n], y[n];
	Float_t			erx[n], ery[n];
	Float_t			preTL = getPreTL();
	Float_t			trackLength = preTL;
	Float_t			maxEnergy, estimatedEnergy;
	Float_t			scaleParameter = 0;
	Float_t			WEPLFactor;
	
	for (Int_t i=0; i<n; i++) {
		if (!At(i)) continue;
		trackLength += getTrackLengthmmAt(i);
		x[i] = trackLength;
		y[i] = getDepositedEnergy(i);
		ery[i] = getDepositedEnergyError(i);
		erx[i] = getAbsorberLength(i);
	}
	
	maxEnergy = getEnergyFromTL(trackLength + 0.75*dz);
	estimatedEnergy = getEnergyFromTL(trackLength + 0.5*dz);

	if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		WEPLFactor = getWEPLFactorFromEnergy(estimatedEnergy);
		for (Int_t i=0; i<n; i++) {
			x[i] = x[i] * WEPLFactor;
		}
	}
	
	graph = new TGraphErrors(n, x, y, erx, ery);
	
	if (kOutputUnit == kPhysical) {
		if (kMaterial == kTungsten) scaleParameter = 14;
		if (kMaterial == kAluminium) scaleParameter = 65;
	}
	
	else if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		if (kMaterial == kTungsten) scaleParameter = 100;
		if (kMaterial == kAluminium) scaleParameter = 126;
	}
	
	TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
	func->SetParameter(0, estimatedEnergy);
	func->SetParameter(1, scaleParameter);
	func->SetParLimits(0, 0, maxEnergy);
	func->SetParLimits(1, scaleParameter, scaleParameter);
	func->SetNpx(500);

	graph->Fit("fit_BP", "B, Q, N, WW", "", 0, 500);
	
 	fitEnergy_ = correctForEnergyParameterisation(func->GetParameter(0));
	fitScale_ = func->GetParameter(1);
	fitError_ = func->GetParError(0);

	return graph;
}

TGraphErrors * Track::doRangeFit(Bool_t isScaleVariable) {
	// Fit a bragg curve on the clusters in this track
	// See HelperFunctions/Tools.C::fitfunc_DBP for the BC function used
	//
	// This version uses the range (z2 - z1), and not the track length
	// The difference should be minimal! However this is how the conversion functions
	// are defined.


	TGraphErrors * graph = nullptr;
	Int_t				n = GetEntriesFast();
	Float_t			x[n], y[n];
	Float_t			erx[n], ery[n];
	Float_t			preTL = getPreTL();
	Float_t			maxEnergy, estimatedEnergy;
	Float_t			scaleParameter = 0;
	Float_t			WEPLFactor;
   Bool_t         checkResistivity = false;

   if (kDataType == kData) {
      checkResistivity = true;
   }
	
	for (Int_t i=0; i<n; i++) {
		if (!At(i)) continue;
		x[i] = preTL + getLayermm(i);
		y[i] = getDepositedEnergy(i, checkResistivity); // normalize to per um
		ery[i] = getDepositedEnergyError(i, checkResistivity);
		erx[i] = dz / sqrt(12);
	}
	
	maxEnergy = getEnergyFromTL(x[n-1] + 0.75*dz);
	estimatedEnergy = getEnergyFromTL(x[n-1] + 0.5*dz);

	if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		WEPLFactor = getWEPLFactorFromEnergy(estimatedEnergy);
		for (Int_t i=0; i<n; i++) {
			x[i] = x[i] * WEPLFactor;
			erx[i] = erx[i] * WEPLFactor;
		}
	}
	
	graph = new TGraphErrors(n, x, y, erx, ery);
	
	if (kOutputUnit == kPhysical) {
		if (kMaterial == kTungsten) scaleParameter = 14/14.;
		if (kMaterial == kAluminium) scaleParameter = 65/14.;
	}
	
	else if (kOutputUnit == kWEPL || kOutputUnit == kEnergy) {
		if (kMaterial == kTungsten) scaleParameter = 3.1; // validated through drawFitScale
		if (kMaterial == kAluminium) scaleParameter = 3.1;
	}

	if (kDataType == kData) scaleParameter = 2.7;

	TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, getWEPLFromEnergy(maxEnergy*1.2), 2);
	func->SetParameter(0, estimatedEnergy);
	func->SetParameter(1, scaleParameter);
	func->SetParLimits(0, 0, maxEnergy);
	func->SetParLimits(1, scaleParameter, scaleParameter);
   if (isScaleVariable) {
      func->SetParLimits(1, 0.01 * scaleParameter, 100 * scaleParameter);
   }
	func->SetNpx(500);

	graph->Fit("fit_BP", "B, N, Q, W", "", 0, getWEPLFromEnergy(maxEnergy*1.2));
	
	fitEnergy_ = correctForEnergyParameterisation(func->GetParameter(0));
//	fitEnergy_ = func->GetParameter(0);
	fitScale_ = func->GetParameter(1);
	fitError_ = func->GetParError(0);

	return graph;
}

Float_t Track::getFitParameterEnergy() {
	if (!fitEnergy_) {
		if (!run_energy) {
			return 0;
		}
		else {
			doRangeFit();
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
			doRangeFit();
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
			doRangeFit();
		}
	}
	return fitError_;
}

Float_t Track::getTrackScore() {
	// The scoring of the track quality
	//
	// If angularChange is 3, give 0 points
	// If its 0, give 10 points
	// If bragg peak is found, give 10 points
	
	Float_t	upperTrackLength = getTLFromEnergy(run_energy);
	Float_t	upperAngularChange = 3;
	Float_t	points = 0;
	Int_t		angularChangePoints = 10;
	Int_t		trackLengthPoints = 25;
	Int_t		braggPeakPoints = 10;
	Float_t	trackLength = getTrackLengthmm();
	Float_t	angularChange = getSlopeAngleDifferenceSumInTheta0();
	
	if (trackLength == 0) return 0;

	points = trackLength * (trackLengthPoints / upperTrackLength);
	points += (upperAngularChange - angularChange) * (angularChangePoints / upperAngularChange);

	if (getAverageCSLastN(2) > getAverageCS() * kBPFactorAboveAverage) {
		points += braggPeakPoints;
	}

	return points;
}
