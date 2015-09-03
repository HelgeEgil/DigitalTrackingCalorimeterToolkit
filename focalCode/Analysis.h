#ifndef Analysis_h
#define Analysis_h

#include "Tracks.h"

Tracks * getTracks(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Int_t energy = 190);
void		DrawBraggPeakFit(Int_t Runs, Int_t dataType = kMC, Int_t energy = 190);
TString getDataTypeString(Int_t dataType);
//void		DrawFrame2D(Int_t Runs, Int_t Layer);
//void		DrawRealFrame2D(Int_t Runs, Int_t Layer);
//void		DrawData3D(Int_t RunFrom, Int_t RunTo);
//void		DrawRealData3D(Int_t RunFrom, Int_t RunTo);
//void		DrawDiffusionCheck(Int_t Runs, Int_t Layer);
void		DrawTracks3D(Int_t Runs, Int_t dataType, Int_t frameType, Int_t energy);
//void		DrawRealTracks3D(Int_t Runs);
//void		WriteClusterFile(Int_t Runs);
//void		FindClusters(Int_t Runs, Int_t Layer);
//void		MakeLayerPNGs(Int_t Runs, Int_t Diffusion);
//void		DrawClusterShapes();
//void		CompareDataAndMC(Int_t Runs, Int_t DataRuns);
//void		Draw2DProjection(Int_t Runs, Int_t Events);
//void		GetTrackerStatistics(Int_t Events, Int_t Runs);
void		GetTrackStatistics(Int_t Runs, Int_t dataType = kMC, Int_t energy = 190);
// void   EdepHistogram();

#endif
