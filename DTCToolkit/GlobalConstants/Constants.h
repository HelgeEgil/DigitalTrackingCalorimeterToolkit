#ifndef Constants_h
#define Constants_h

#include <cstring>
#include <vector>
#include <TObject.h>

// Set these compiler directives to adjust the logic flow of this file

// #define USEALPIDE // Comment to use experimental data; uncomment to use Monte Carlo data
#define USEDEBUG // Uncomment to print more debug information
// #define ONECHIP // (MC) Limit the data area to a single chip (simplify visualization & memory usage)

#ifdef USEDEBUG
#define showDebug(x) std::cout << x
#else
#define showDebug(x)
#endif

#ifdef USEALPIDE // Monte Carlo
#define NX 9000
#define DX 0.030
#define DY 0.030
#define NY 4500
#define DZ 0.435
#else // exp. data
#define NX 1280
#define NY 1280
#define DX 0.03
#define DY 0.03
#define DZ 0.975
#endif

#ifdef ONECHIP // In case 9000 x 4500 x (number of layers) is too taxing
#undef NX
#undef NY
#define NX 1024
#define NY 512
#endif

// -------------------------------------

Bool_t   kIsAluminumPlate = false;
Bool_t   kIsScintillator = false;
Bool_t   kIsFirstLayerAir = false;
Bool_t   kDoTracking = true;
Bool_t   kUseEmpiricalMCS = true;
Bool_t   kFilterNuclearInteractions = false;
Int_t    kEventsPerRun = 500;

#ifdef USEALPIDE
Bool_t   kUseDegrader = true;
Bool_t   kUseAlpide = true;
#else
Bool_t   kUseAlpide = false;
Bool_t   kUseDegrader = false;
#endif

const Int_t sizeOfEventID = 25;
const Int_t nChildrenInNode = 2; // max concurrent track segments to follow
Float_t kMaxTrackScore = 0.3; // cumulative rad // was 0.19

// natural unit is mm
const Float_t cm = 0.1;
const Float_t um = 1000;
const Float_t kRad = 3.14159265/180.;

// Some general run parameters
const Int_t nx = NX;
const Int_t ny = NY;
const Int_t nTrackers = 4;

#ifdef USEALPIDE
const Float_t kAbsorberThickness = 3.5; // ALPIDE, CHANGE TO FIT MC DATA GEOMETRY
#else
const Float_t kAbsorberThickness = 3.3; // FOCAL EXPERIMENTAL DATA, DON'T CHANGE
#endif

// nLayers are loaded in MaterialConstants.C according to the detector geometry
const Float_t dx = DX; // mm
const Float_t dy = DY; // mm
const Float_t dz = DZ + kAbsorberThickness;

// Used for treatment of available experimental data files
const Int_t nEnergies = 6;
Int_t energies[nEnergies] = {122, 140, 150, 170, 180, 188};

enum eFrameType {kCalorimeter, kTracker};
enum eDataType {kMC, kData};
enum eMaterial {kTungsten, kAluminum, kPMMA, kWater, kCarbon};
enum eOutputUnit {kPhysical, kWEPL, kEnergy};

#ifdef USEALPIDE
const Int_t kMaterial = kAluminum;
#else
const Int_t kMaterial = kTungsten;
#endif

Int_t kOutputUnit = kPhysical;

/*
 * Track reconstruction algorithms
 *  recursive is excellent on small frames, really slow on large
 *  nearestCluster is good on everything
 */

enum eTrackFindingAlgorithm {kWeightedRecursive, kNearestCluster};
// const Int_t kTrackFindingAlgorithm = kNearestCluster;
const Int_t kTrackFindingAlgorithm = kWeightedRecursive;

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
const Bool_t kUseTrackSplitting = true;

// Use refined clustering model -- empirical model with updated parameters
const Bool_t kUseRefinedClustering = true;

// Tracking parameters
Float_t initialSearchRadius = 50 * dx; // 20 if dimensionless
Float_t searchRadius = 40 * dx; // 20 if dimensionless
Float_t kMCSFactor = 3;
const Double_t SpreadNumber = 0.1; // number of iterations in gaussian spread
const Double_t SpreadSigma  = 0.7; // sigma (in pixels) to spread

// Minimum track length (in mm) to account for the track in the bragg peak fitting analysis.
const Int_t kMinimumTracklength = 5;

// How much above the average edep must the bragg peak (last two layers) be?
const Float_t kBPFactorAboveAverage = 1.3;

Int_t    GlobalLayerID = 0;
Float_t  run_energy = 0;
Float_t  run_degraderThickness = 0;
Int_t    kDataType = kData;
#endif
