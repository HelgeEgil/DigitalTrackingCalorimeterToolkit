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

const Int_t nPLEnergies = 29;
Int_t kPLEnergies[nPLEnergies]
	= {1, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130,
	   140, 150, 160, 170, 180, 190, 200, 225, 250, 275, 300, 350, 400};

/*
 * Choose material used in simulation
 * Comment out the irrelevant ones...
 */


// MATERIAL CONSTANTS CALCULATED USING SRIM VERSION 2013

// Float_t kPLFocal[nPLEnergies] // Tungsten SRIM
// 	= {0.005, 0.026, 0.057, 0.169, 0.527, 1.040, 1.690, 2.470, 3.380, 4.400,
// 		5.520, 6.750, 8.080, 9.510, 11.030, 12.630, 14.320, 16.090, 17.940,
// 		19.870, 21.860, 23.950, 26.060, 31.680, 37.670, 44.010, 50.670, 64.870, 80.080};
// 	   
// Float_t kWEPLRatio[nPLEnergies] // Tungsten SRIM
// 	= {5.009434, 5.700960, 6.259096, 7.093456, 7.915417, 8.355769, 8.650888, 8.862348,
// 		9.002959, 9.127273, 9.244565, 9.333333, 9.409653, 9.470032, 9.525839, 9.581156,
// 		9.627095, 9.669360, 9.706243, 9.738299, 9.774016, 9.793319, 9.830775, 9.890152,
// 		9.942925, 9.988412, 10.028419, 10.094497, 10.150475};
		
Float_t kPLFocal[nPLEnergies] // Aluminium SRIM
	= {0.01532, 0.0823, 0.19169, 0.62863, 2.12, 4.36, 7.28, 10.84, 14.99, 19.72, 24.99,
		30.78, 37.07, 43.84, 51.07, 58.74, 66.84, 75.35, 84.25, 93.54, 103.2, 113.21, 123.57,
		150.89, 180.15, 211.2, 243.91, 313.83, 389.08};
		
Float_t kWEPLRatio[nPLEnergies] // Aluminium SRIM
	= {1.733029, 1.804496, 1.848766, 1.908913, 1.966981, 1.993119, 2.008242, 2.019373, 
		2.030020, 2.036511, 2.042017, 2.046784, 2.050985, 2.054288, 2.057372, 2.060095, 
		2.062537, 2.064764, 2.066825, 2.068634, 2.070349, 2.071813, 2.073238, 2.076480,
		2.079101, 2.081392, 2.083309, 2.086576, 2.089159};
		
// Float_t kPLFocal[nPLEnergies] // PMMA SRIM
// 	= {0.02684, 0.15544, 0.37317, 1.27, 4.41, 9.18, 15.47, 23.16, 32.19, 42.49, 54, 66.67,
// 		80.45, 95.32, 111.21, 128.09, 145.94, 164.71, 184.37, 204.89, 226.24, 248.38, 271.32,
// 		331.88, 396.82, 465.8, 538.54, 694.25, 862.09};
// 
// Float_t kWEPLRatio[nPLEnergies] // PMMA SRIM
// 	= {0.989195, 0.955417, 0.949674, 0.944882, 0.945578, 0.946623, 0.945055, 0.945164, 
// 		0.945325, 0.945164, 0.945000, 0.944953, 0.945059, 0.944817, 0.944789, 0.944726,
// 		0.944635, 0.944569, 0.944460, 0.944409, 0.944395, 0.944319, 0.944236, 0.944076,
// 		0.943879, 0.943731, 0.943551, 0.943219, 0.942883};

// THESE VALUES WERE CALCULATED USING THE OLD MATERIAL CONSTANTS

// Float_t kPLFocal[nPLEnergies] // Tungsten
// 	= {0.0358, 0.0560, 0.0860, 0.1964, 0.5536, 1.0740, 1.7440, 2.5480, 3.4590, 4.5150, 5.6860,
// 		6.9670, 8.3520, 9.8390, 11.441, 13.117, 14.880, 16.729, 18.660, 20.670, 22.758, 24.919,
// 		27.147, 33.023, 39.294, 45.930, 52.901, 67.750, 83.668};
// 
// Float_t kWEPLRatio[nPLEnergies] // Tungsten
// 	= {0.7430, 2.6589, 4.1628, 6.2062, 7.6553, 8.2123, 8.5161, 8.7253, 8.9306, 9.0286, 9.1094,
// 	   9.1771, 9.2364, 9.2863, 9.3151, 9.3567, 9.3947, 9.4289, 9.4602, 9.4891, 9.5150, 9.5393,
// 	   9.5634, 9.3139, 9.6572, 9.6949, 9.7285, 9.7862, 9.8340};


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
