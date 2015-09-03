#ifndef Constants_h
#define Constants_h

#include <TObject.h>
#include <vector>

const Int_t nx = 640;
const Int_t ny = 640;

const Int_t nLayers = 24;
const Int_t nTrackers = 4;

const Float_t dx = 0.03; // mm
const Float_t dy = 0.03; // mm
const Float_t dz = 3.84; // mm

const Float_t initialSearchRadius = 50 * dx; // 20 if dimensionless
const Float_t searchRadius = 40 * dx; // 20 if dimensionless

// Diffusion constants

const Double_t SpreadNumber = 0.1; // number of iterations in gaussian spread
const Double_t SpreadSigma  = 0.7; // sigma (in pixels) to spread

const Bool_t kCalorimeter = 0;
const Bool_t kTracker = 1;

const Bool_t kMC = 0;
const Bool_t kData = 1;

const Int_t kEventsPerRun = 150;

Int_t energies[7] = {122, 140, 150, 160, 170, 180, 190};
std::vector<Int_t> kPossibleEnergies(&energies[0], &energies[0]+7);

/*
 * Track reconstruction algorithms
 *  recursive is excellent on small frames, really slow on large
 *  nearestCluster is good on everything
 */

const Int_t kRecursive = 0;
const Int_t kNearestCluster = 1;
const Int_t kTrackFindingAlgorithm = kNearestCluster;

// For the bragg peak analysis
const Int_t kMinimumTracklength = 22;

#endif
