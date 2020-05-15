#ifndef trackFitting_cxx
#define trackFitting_cxx


#include <iostream>
#include <cmath>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TStopwatch.h>

#include "Classes/Track/Track.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Hit/Hit.h"
#include "HelperFunctions/Tools.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"

TGraphErrors * Track::doTrackFit(Bool_t isScaleVariable, Bool_t useTrackLength) {
   // Fit a bragg curve on the clusters in this track
   // See HelperFunctions/Tools.C::fitfunc_DBP for the BC function used
   //
   // This version uses the range (z2 - z1), and not the track length
   // The difference should be minimal! However this is how the conversion functions
   // are defined.

   alpha = alpha_water;
   p = p_water;

   TGraphErrors * graph = nullptr;
   Int_t          n = GetEntriesFast();
   Float_t        x[n], y[n];
   Float_t        erx[n], ery[n];
   Float_t        trackLength = 0;
   Float_t        maxRange, minRange, estimatedRange;
   Float_t        scaleParameter = 0;
   Float_t        overFittingDistance, startFittingDistance;
   Bool_t         checkResistivity = false;

   if (kDataType == kData) {
      checkResistivity = true;
   }

   for (Int_t i=0; i<n; i++) {
      if (!At(i)) continue;

      x[i] = getWETFromLayer(i);
      y[i] = getDepositedEnergy(i, checkResistivity);
      ery[i] = getDepositedEnergyError(i, checkResistivity);
      erx[i] = dz * 0.28867; // 1/sqrt(12)
   }

   // how much beyond the last measurement the fit is allowed to go
   overFittingDistance = dz;
   startFittingDistance = dz * 0.5;

   minRange = x[n-2];
   maxRange = x[n-1] + overFittingDistance;
   estimatedRange = x[n-1] + startFittingDistance;

   graph = new TGraphErrors(n, x, y, erx, ery); // maybe a speedup here is possible
   
   scaleParameter = 0.73 / (p * pow(alpha, 1/p)); // Found from fits to all tracks

   if (kDataType == kData) {
      scaleParameter = 3.42;
   }

   scaleParameter = 7;
   if (kHelium) scaleParameter = 28; // was 22.5

   // scaleParameter *= 1.2; // Empirical tests to reduce range bias

   TF1 * func = new TF1("fit_BP", fitfunc_DBP, 0, maxRange, 2);
   func->SetParameter(0, estimatedRange);
   func->SetParameter(1, scaleParameter);
   func->SetParLimits(0, minRange, maxRange);
   func->SetParLimits(1, scaleParameter, scaleParameter);

   if (isScaleVariable) {
      func->SetParLimits(1, 0.01 * scaleParameter, 100 * scaleParameter);
   }

   func->SetNpx(750);

   graph->Fit("fit_BP", "B, N, Q, W, G", "", 0, maxRange*1.2);

   fitRange_ = func->GetParameter(0);
   fitScale_ = func->GetParameter(1);
   fitError_ = func->GetParError(0);
   fitChi2_  = func->GetChisquare();

   fitSigma_ = 0;
   Int_t ii;
   for (Int_t i=1; i<6; i++) {
      ii = GetEntriesFast()-i;
      fitSigma_ += pow(y[ii] - func->Eval(x[ii]), 2);
   }
   fitSigma_ = sqrt(fitSigma_ /(5-1));

   delete func;
   
   return graph;
}

Float_t Track::getFitParameterRange() {
   if (!fitRange_) {
      if (!run_energy) {
         return 0;
      }
      else {
         TGraphErrors *out = doTrackFit();
         delete out;
      }
   }

   return fitRange_;
}

Float_t Track::getFitParameterScale() {
   if (!fitScale_) {
      if (!run_energy) {
         return 0;
      }
      else {
         TGraphErrors *out = doTrackFit();
         delete out;
      }
   }
   return fitScale_;
}

Float_t Track::getFitParameterSigma() {
   return fitSigma_;
}

Float_t Track::getFitParameterError() {
   if (!fitScale_) {
      if (!run_energy) {
         return 0;
      }
      else {
         TGraphErrors *out = doTrackFit();
         delete out;
      }
   }
   return fitError_;
}

Float_t Track::getFitParameterChiSquare() {
   if (!fitChi2_) {
      if (!run_energy) {
         return 0;
      }
      else {
         TGraphErrors *out = doTrackFit();
         delete out;
      }
   }
   return fitChi2_;
}

Float_t Track::getTrackScore() {
   // The scoring of the track quality
   //
   // If angularChange is 3, give 0 points
   // If its 0, give 10 points
   // If bragg peak is found, give 10 points
  
   if (!kUseAlpide) {
      Float_t  upperTrackLength = getTLFromEnergy(run_energy);
      Float_t  upperAngularChange = 3;
      Float_t  points = 0;
      Int_t    angularChangePoints = 10;
      Int_t    trackLengthPoints = 25;
      Int_t    braggPeakPoints = 10;
      Float_t  trackLength = getTrackLengthmm();
      Float_t  angularChange = getSlopeAngleDifferenceSumInTheta0();
      
      if (trackLength == 0) return 0;

      points = trackLength * (trackLengthPoints / upperTrackLength);
      points += (upperAngularChange - angularChange) * (angularChangePoints / upperAngularChange);

      if (getAverageCSLastN(2) > getAverageCS() * kBPFactorAboveAverage) {
         points += braggPeakPoints;
      }

      return points;
   }

   else {
      Float_t  angularChange = 0;
      Int_t    idx1, idx2, idx3;

      for (Int_t layer=0; layer<nLayers; layer++) {
         idx1 = getClusterFromLayer((layer>0) ? layer-1 : layer);
         idx2 = getClusterFromLayer(layer);
         idx3 = getClusterFromLayer(layer+1);

         if (idx1<0 || idx2<0 || idx3<0) continue;
         
         angularChange += pow(getDotProductAngle(At(idx1), At(idx2), At(idx3)) / getEmpiricalMCSAngle(layer), 2);
      }

      angularChange = sqrt(angularChange);

      showDebug("Track::getTrackScore: Angular change chi2 = " << angularChange << " and tracklength = " << getTrackLengthmm() << endl);
      return getTrackLengthmm();
   }
   return 0;
}
#endif
