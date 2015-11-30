#ifndef Constants_h
#define Constants_h

#include <TObject.h>
#include <cstring>
#include <vector>

const Int_t nx = 640;
const Int_t ny = 640;

// const Int_t nLayers = 41; // Al
const Int_t nLayers = 24; // W
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

const Int_t nPLEnergies = 29;
Int_t kPLEnergies[nPLEnergies]
	= {1, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130,
	   140, 150, 160, 170, 180, 190, 200, 225, 250, 275, 300, 350, 400};

/*
 * Choose material used in simulation
 * Comment out the irrelevant ones...
 */

Float_t kPLFocal[nPLEnergies] // Tungsten
	= {0.0358, 0.0560, 0.0860, 0.1964, 0.5536, 1.0740, 1.7440, 2.5480, 3.4590, 4.5150, 5.6860,
		6.9670, 8.3520, 9.8390, 11.441, 13.117, 14.880, 16.729, 18.660, 20.670, 22.758, 24.919,
		27.147, 33.023, 39.294, 45.930, 52.901, 67.750, 83.668};

Float_t kWEPLRatio[nPLEnergies] // Tungsten
	= {0.7430, 2.6589, 4.1628, 6.2062, 7.6553, 8.2123, 8.5161, 8.7253, 8.9306, 9.0286, 9.1094,
	   9.1771, 9.2364, 9.2863, 9.3151, 9.3567, 9.3947, 9.4289, 9.4602, 9.4891, 9.5150, 9.5393,
	   9.5634, 9.3139, 9.6572, 9.6949, 9.7285, 9.7862, 9.8340};


// Float_t kPLFocal[nPLEnergies] // Aluminum
// 	= {0.0234, 0.0090, 0.2000, 0.6440, 2.1698, 4.4619, 7.4556, 11.104, 15.371, 20.224, 25.637,
// 		31.583, 38.040, 44.989, 52.408, 60.280, 68.583, 77.313, 86.449, 95.796, 105.88, 116.15,
// 		126.78, 154.81, 184.81, 216.63, 250.13, 321.72,	398.70};
//
// Float_t kWEPLRatio[nPLEnergies] // Aluminum
// 	= {1.137, 1.654, 1.790, 1.892, 1.953, 1.977, 1.992, 2.002, 2.009, 2.016, 2.020, 2.024, 2.028,
// 		2.031, 2.034, 2.036, 2.038, 2.040, 2.042, 2.047, 2.045, 2.047, 2.048, 2.051, 2.053, 2.056,
// 		2.057, 2.061, 2.064};

// const Float_t kPLFocal[nPLEnergies] // PMMA
// 	= {0.0247, 0.1270, 0.3010, 1.0145, 3.5084, 7.2895, 12.255, 18.327, 25.448, 33.564, 42.630, 52.605,
// 		63.450, 75.132, 87.618,	100.877, 114.882,	129.604,	145.022,	161.110,	177.846,	195.207,	213.174,
// 		260.614,	311.435,	365.383,	422.227,	543.793,	674.678};
//
// const Float_t kWEPLRatio[nPLEnergies] // PMMA
// 	= {1.0769, 1.1727, 1.1894,	1.2015, 1.2080, 1.2100,	1.2120, 1.2131, 1.2139,	1.2145, 1.2150, 1.2154,
// 		1.2158, 1.2161, 1.2164, 1.2166, 1.2168, 1.2171, 1.2173, 1.2174, 1.2176, 1.2177, 1.2179, 1.2182,
// 		1.2184, 1.2187, 1.2189, 1.2192, 1.2195};

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
