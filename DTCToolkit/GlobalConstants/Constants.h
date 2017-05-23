#ifndef Constants_h
#define Constants_h

#include <cstring>
#include <vector>
#include <TObject.h>

#define USEALPIDE

#ifdef USEDEBUG
#define showDebug(x) std::cout << x
#else
#define showDebug(x)
#endif


// Was 3072x1024 w/29.3 um px
// With 270x135 and 30 um px this should be 9000x4500
#ifdef USEALPIDE
#define NX 9000
#define DX 0.030
#define DY 0.030
#define NY 4500
#define DZ 0.435
#else
#define NX 1280
#define NY 1280
#define DX 0.03
#define DY 0.03
#define DZ 0.975
#endif

enum eFrameType {kCalorimeter, kTracker};
enum eDataType {kMC, kData};

Float_t  run_energy = 0;
Float_t  run_degraderThickness = 0;
Bool_t   kIsAluminumPlate = false;
Bool_t   kIsScintillator = false;
Bool_t   kIsFirstLayerAir = true;
Bool_t   kUseAlpide = true;
Bool_t   kDoTracking = false;
Bool_t   kUseEmpiricalMCS = true;
Bool_t   kFilterNuclearInteractions = true;
Bool_t   useDegrader = true;

const Int_t sizeOfEventID = 500;

// natural unit is mm
const Float_t cm = 0.1;
const Float_t um = 1000;
const Float_t kRad = 3.14159265/180.;

// Some general run parameters
const    Int_t nx = NX;
const    Int_t ny = NY;
const    Int_t nTrackers = 4;
const    Float_t kAbsorbatorThickness = 2;
// FOCAL IS 3 mm (2x absorbers = 3 mm)

// nLayers are loaded in MaterialConstants.C according to the detector geometry
const Float_t dx = DX; // mm
const Float_t dy = DY; // mm
const Float_t dz = DZ + kAbsorbatorThickness;
Int_t kEventsPerRun = 800000;

// Used for treatment of available experimental data files
const Int_t nEnergies = 8;
Int_t energies[nEnergies] = {122, 140, 150, 160, 170, 180, 188, 190};

enum eMaterial {kTungsten, kAluminum, kPMMA, kWater, kCarbon};
const Int_t kMaterial = kAluminum;

Int_t kDataType = kMC;

enum eOutputUnit {kPhysical, kWEPL, kEnergy};
Int_t kOutputUnit = kWEPL;

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
Float_t initialSearchRadius = 50 * dx; // 20 if dimensionless
Float_t searchRadius = 40 * dx; // 20 if dimensionless
Float_t kMCSFactor = 3;
const Double_t SpreadNumber = 0.1; // number of iterations in gaussian spread
const Double_t SpreadSigma  = 0.7; // sigma (in pixels) to spread

// Minimum track length (in mm) to account for the track in the bragg peak fitting analysis.
const Int_t kMinimumTracklength = 5;

// How much above the average edep must the bragg peak (last two layers) be?
const Float_t kBPFactorAboveAverage = 1.3;

#endif
