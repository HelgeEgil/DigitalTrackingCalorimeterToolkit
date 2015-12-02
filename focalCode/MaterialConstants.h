#ifndef MaterialConstants_h
#define MaterialConstants_h

#include <TObject.h>

const Int_t nPLEnergies = 70;

Float_t kPLFocal[nPLEnergies];
Float_t kWEPLRatio[nPLEnergies];
Int_t nLayers;

Float_t kPLEnergies[nPLEnergies] =
		   {1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2, 2.25, 2.5, 2.75,
			3, 3.25, 3.5, 3.75, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14,
			15, 16, 17, 18, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40, 45, 50, 55,
			60, 65, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 225,
			250, 275, 300, 325, 350, 375, 400};

Float_t kPLFocal_W[nPLEnergies] = // Tungsten
		  {0.005770, 0.006580, 0.007430, 0.008320, 0.009240, 0.010190,
			0.011170, 0.012180, 0.013230, 0.015400, 0.018280, 0.021330,
			0.024550, 0.027940, 0.031500, 0.035210, 0.039080, 0.043110,
			0.051610, 0.060710, 0.070390, 0.080640, 0.091430, 0.102770,
			0.127030, 0.153350, 0.181680, 0.211980, 0.244190, 0.278280,
			0.314220, 0.351960, 0.391490, 0.432770, 0.475790, 0.566880,
			0.690070, 0.823410, 0.966620, 1.120000, 1.280000, 1.450000,
			1.630000, 1.820000, 2.230000, 2.670000, 3.140000, 3.640000,
			4.180000, 4.740000, 5.960000, 7.290000, 8.730000, 10.27000,
			11.91000, 13.64000, 15.47000, 17.38000, 19.38000, 21.46000,
			23.62000, 28.16000, 34.24000, 40.72000, 47.58000, 54.70000,
			62.33000, 70.16000, 78.28000, 86.65000};

Float_t kWEPLRatio_W[nPLEnergies] = // Tungsten
		  {4.601386, 4.674772, 4.725437, 4.762019, 4.797619, 4.832188,
			4.866607, 4.899836, 4.928949, 4.994805, 5.074945, 5.157056,
			5.237475, 5.315319, 5.389524, 5.462653, 5.532497, 5.598469,
			5.723503, 5.837424, 5.942037, 6.038318, 6.128404, 6.211735,
			6.361332, 6.493772, 6.605020, 6.698745, 6.797985, 6.899526,
			6.969639, 7.046255, 7.126619, 7.186265, 7.251098, 7.356054,
			7.477502, 7.578242, 7.676233, 7.758929, 7.843750, 7.917241,
			7.981595, 8.032967, 8.116592, 8.198502, 8.283439, 8.359890,
			8.409091, 8.472574, 8.562081, 8.641975, 8.709049, 8.769231,
			8.821998, 8.871701, 8.911441, 8.951669, 8.985036, 9.016775,
			9.045724, 9.097656, 9.150701, 9.198183, 9.238966, 9.274320,
			9.304829, 9.333381, 9.357690, 9.380842};

Float_t kPLFocal_Al[nPLEnergies] = // Aluminium
		   {0.015180, 0.017430, 0.019810, 0.022310, 0.024930, 0.027660,
			0.030500, 0.033460, 0.036530, 0.043000, 0.051720, 0.061130,
			0.071230, 0.082000, 0.093420, 0.105500, 0.118210, 0.131560,
			0.160100, 0.191080, 0.224460, 0.260190, 0.298220, 0.338530,
			0.425780, 0.521800, 0.626390, 0.739400, 0.860660, 0.990060,
			1.130000, 1.270000, 1.430000, 1.590000, 1.750000, 2.110000,
			2.600000, 3.140000, 3.720000, 4.340000, 5.000000, 5.710000,
			6.460000, 7.240000, 8.940000, 10.78000, 12.77000, 14.91000,
			17.20000, 19.62000, 24.85000, 30.61000, 36.86000, 43.58000,
			50.77000, 58.39000, 66.44000, 74.89000, 83.74000, 92.97000,
			102.5700, 122.8100, 149.9600, 179.0300, 209.8800, 242.3800,
			276.4000, 311.8400, 348.6100, 386.6000};

Float_t kWEPLRatio_Al[nPLEnergies] = // Aluminium
		   {1.749012, 1.764773, 1.772337, 1.775885, 1.778179, 1.780188,
			1.782295, 1.783622, 1.785108, 1.788837, 1.793697, 1.799444,
			1.805138, 1.811098, 1.817277, 1.823128, 1.829033, 1.834524,
			1.845034, 1.854668, 1.863406, 1.871440, 1.878881, 1.885741,
			1.897882, 1.908432, 1.915739, 1.920476, 1.928752, 1.939276,
			1.938053, 1.952756, 1.951049, 1.955975, 1.971429, 1.976303,
			1.984615, 1.987261, 1.994624, 2.002304, 2.008000, 2.010508,
			2.013932, 2.019337, 2.024609, 2.030612, 2.036805, 2.040912,
			2.043605, 2.046891, 2.053521, 2.058151, 2.062670, 2.066544,
			2.069529, 2.072444, 2.074955, 2.077447, 2.079412, 2.081317,
			2.083065, 2.086068, 2.089357, 2.092107, 2.094483, 2.096460,
			2.098300, 2.099891, 2.101259, 2.102561};
			
Float_t kPLFocal_PMMA[nPLEnergies] = // PMMA
		   {0.035100, 0.040440, 0.046040, 0.051900, 0.058020, 0.064410,
			0.071060, 0.077970, 0.085130, 0.100250, 0.120600, 0.142560,
			0.166120, 0.191230, 0.217890, 0.246060, 0.275730, 0.306870,
			0.373460, 0.445750, 0.523630, 0.606990, 0.695750, 0.789830,
			0.993440, 1.220000, 1.460000, 1.730000, 2.010000, 2.310000,
			2.630000, 2.970000, 3.330000, 3.700000, 4.100000, 4.930000,
			6.080000, 7.330000, 8.680000, 10.13000, 11.68000, 13.33000,
			15.08000, 16.91000, 20.86000, 25.17000, 29.82000, 34.82000,
			40.15000, 45.80000, 58.03000, 71.46000, 86.06000, 101.7700,
			118.5500, 136.3500, 155.1400, 174.8900, 195.5600, 217.1200,
			239.5400, 286.8100, 350.2400, 418.1700, 490.2500, 566.1800,
			645.6900, 728.5200, 814.4500, 903.2700};

Float_t kWEPLRatio_PMMA[nPLEnergies] = // PMMA
		   {0.756410, 0.760633, 0.762598, 0.763391, 0.764047, 0.764478,
			0.764987, 0.765423, 0.766005, 0.767282, 0.769237, 0.771605,
			0.774019, 0.776604, 0.779155, 0.781679, 0.784137, 0.786489,
			0.790955, 0.795042, 0.798770, 0.802204, 0.805347, 0.808250,
			0.813416, 0.816246, 0.821918, 0.820809, 0.825871, 0.831169,
			0.832700, 0.835017, 0.837838, 0.840541, 0.841463, 0.845842,
			0.848684, 0.851296, 0.854839, 0.857848, 0.859589, 0.861215,
			0.862732, 0.864577, 0.867689, 0.869686, 0.872233, 0.873923,
			0.875467, 0.876856, 0.879373, 0.881612, 0.883453, 0.884937,
			0.886293, 0.887495, 0.888617, 0.889588, 0.890417, 0.891212,
			0.891960, 0.893239, 0.894587, 0.895688, 0.896665, 0.897488,
			0.898217, 0.898850, 0.899405, 0.899897};

// THESE VALUES WERE CALCULATED USING THE OLD MATERIAL CONSTANTS

//const Int_t nPLEnergies = 29;
//Int_t kPLEnergies[nPLEnergies]
//	= {1, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130,
//	   140, 150, 160, 170, 180, 190, 200, 225, 250, 275, 300, 350, 400};

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

#endif