#ifndef Constants_h
#define Constants_h

#include <TObject.h>
#include <cstring>
#include <vector>

enum eFrameType {kCalorimeter, kTracker};
enum eDataType {kMC, kData};

Float_t run_energy = 0;
Bool_t kIsAluminumPlate = false;

Bool_t kDebug = false;

// natural unit is mm
const Float_t cm = 0.1;
const Float_t um = 1000;


// Some general run parameters
const Int_t nx = 640;
const Int_t ny = 640;
const Int_t nTrackers = 4;

// nLayers are loaded in WEPL-LUT.C according to the detector geometry
const Float_t dx = 0.03; // mm
const Float_t dy = 0.03; // mm
const Float_t dz = 3.84; // mm
const Int_t kEventsPerRun = 100;

// Used for treatment of available experimental data files
const Int_t nEnergies = 8;
Int_t energies[nEnergies] = {122, 140, 150, 160, 170, 180, 188, 190};

/*
 * Material used for the conversion between the absorber material and the energy
 * and/or the water equivalent path length of the detector
 * 
 * The actual lists (kPLFocal for energies and kWEPLRatio for WEPL conversion) are
 * loaded in WEPL-LUT.C.
 * 
*/

enum eMaterial {kTungsten, kAluminum, kPMMA, kWater, kFocalTungsten, kFocalAluminum};
const Int_t kMaterial = kAluminum;

enum eOutputUnit {kPhysical, kWEPL, kEnergy};
Int_t kOutputUnit = kPhysical;

/*
 * Track reconstruction algorithms
 *  recursive is excellent on small frames, really slow on large
 *  nearestCluster is good on everything
 */

enum eTrackFindingAlgorithm {kRecursive, kNearestCluster};
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

// Tracking parameters
const Float_t initialSearchRadius = 50 * dx; // 20 if dimensionless
const Float_t searchRadius = 40 * dx; // 20 if dimensionless
const Double_t SpreadNumber = 0.1; // number of iterations in gaussian spread
const Double_t SpreadSigma  = 0.7; // sigma (in pixels) to spread


// Minimum track length (in mm) to account for the track in the bragg peak fitting analysis.
const Int_t kMinimumTracklength = 5;

// How much above the average edep must the bragg peak (last two layers) be?
const Float_t kBPFactorAboveAverage = 1.3;

#endif
