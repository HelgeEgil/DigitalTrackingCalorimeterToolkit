#ifndef Constants_h
#define Constants_h

#include <TObject.h>

const Int_t nx = 640;
const Int_t ny = 640;

const Int_t nLayers = 24;
const Int_t nTrackers = 4;

const Float_t dx = 0.03; // mm
const Float_t dy = 0.03; // mm
const Float_t dz = 3.84; // mm

const Float_t searchRadius = 50 * dx; // 20 if dimensionless
const Float_t secondSearchRadius = 50 * dx; // 20 if dimensionless

// Diffusion constants

const Double_t SpreadNumber = 0.1; // number of iterations in gaussian spread
const Double_t SpreadSigma  = 0.7; // sigma (in pixels) to spread

const Bool_t kCalorimeter = 0;
const Bool_t kTracker = 1;

const Bool_t kMC = 0;
const Bool_t kData = 1;

const Int_t kEventsPerRun = 150;

#endif
