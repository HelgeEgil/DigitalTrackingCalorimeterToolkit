#ifndef Analysis_h
#define Analysis_h

#include "Tracks.h"

Tracks * getTracks(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188);
void		drawBraggPeakFit(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188);
void 		drawBraggPeakGraphFit(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188);
void saveTracks(Tracks* tracks, Int_t dataType, Float_t energy);
Tracks* loadTracks(Int_t Runs, Int_t dataType, Float_t energy);
Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy);
void		DrawFrame2D(Int_t Runs, Int_t Layer);
//void		DrawData3D(Int_t RunFrom, Int_t RunTo);
//void		DrawRealData3D(Int_t RunFrom, Int_t RunTo);
void		drawDiffusionCheck(Int_t Runs, Int_t Layer);
void		drawTracks3D(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = false, Float_t energy = 188);
void		drawData3D(Int_t Runs);
//void		DrawRealTracks3D(Int_t Runs);
void		writeClusterFile(Int_t Runs, Int_t dataType, Float_t energy);
//void		MakeLayerPNGs(Int_t Runs, Int_t Diffusion);
void		drawClusterShapes(Int_t Runs, Bool_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188);
//void		CompareDataAndMC(Int_t Runs, Int_t DataRuns);
void 		Draw2DProjection(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy);
//void		GetTrackerStatistics(Int_t Events, Int_t Runs);
void		getTrackStatistics(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188);
// void   EdepHistogram();
Bool_t getCutTrackLength(Float_t energy, Track *track);
Bool_t getCutWEPL(Track *track);
Bool_t getCutChipNumber(Track *track);
Bool_t getCutBraggPeakInTrack(Track *track);
#endif
