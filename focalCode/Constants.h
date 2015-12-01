#ifndef Constants_h
#define Constants_h

#include <TObject.h>
#include <cstring>
#include <vector>

const Int_t nx = 640;
const Int_t ny = 640;

const Int_t nLayers = 41; // Al
// const Int_t nLayers = 24; // W
//const Int_t nLayers = 65; // PMMA
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
 * Choose material used in simulation
 * Comment out the irrelevant ones...
 */

// MATERIAL CONSTANTS CALCULATED USING SRIM VERSION 2013

#include "wepl_arrays.h"

/*
 * Track reconstruction algorithms
 *  recursive is excellent on small frames, really slow on large
 *  nearestCluster is good on everything
 */

const Int_t kRecursive = 0;
const Int_t kNearestCluster = 1;
const Int_t kTrackFindingAlgorithm = kNearestCluster;

/*
 * Track splitting
 *  Look through all tracks after track finding. If any tracks miss
 *  a cluster in the vicinity of a track crossing (looks for track connected
 *  clusters near the /interpolated/ missing cluster positing (between track cluster i-1 and i+1) )
 *  Split the cluster into two, and connect the other half to the other track.
 *
 *  The new cluster sizes range in size from 50 % to 100 % of the original cluster size, depending on
 *  how close the tracks collide.
 */

const Bool_t kUseTrackSplitting = kTRUE;

/*
 * Minimum track length (in mm) to account for the track in the bragg peak fitting analysis.
 */

const Int_t kMinimumTracklength = 22;


#endif
