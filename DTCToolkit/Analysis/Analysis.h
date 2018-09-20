#ifndef Analysis_h
#define Analysis_h

class TCanvas;
class TH1F;
class TF1;

#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"

using namespace DTC;

void        drawTracksDeltaTheta(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness);
void        drawTracksDeltaThetaEachLayer(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness);
void			getTracksReconstructionEfficiency(Int_t dataType, Float_t energy = 250, Float_t degraderThickness = 0);
void			drawClusterShapes(Int_t Runs, Bool_t dataType = kMC, Float_t energy = 250, Float_t degraderThickness = 0);
void			drawFitScale(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy);
void   		drawTracksRangeHistogram(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 250, Float_t degraderThickness = 0, Int_t eventsPerRun = kEventsPerRun, Int_t outputFileIdx = -1, Bool_t drawFitResults = true, Bool_t doTracking = kDoTracking, Bool_t excludeNuclearInteractions = kFilterNuclearInteractions, Int_t skipTracks = 0);
void   		findTracksRangeAccuracy(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 250, Float_t degraderThickness = 0, Int_t eventsPerRun = kEventsPerRun, Int_t outputFileIdx = -1, Bool_t doTracking = kDoTracking, Bool_t excludeNuclearInteractions = kFilterNuclearInteractions);
void   		drawTracksDepthDose(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 250, Float_t degraderThickness = 0, Int_t eventsPerRun = kEventsPerRun, Bool_t doTracking = kDoTracking, Bool_t excludeNuclearInteractions = kFilterNuclearInteractions);
void			draw2DProjection(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy);
void			drawClusterSizeDistribution(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy);
void        compareChargeDiffusionModels(Int_t Runs, Bool_t recreate, Float_t energy);
void			drawTracks3D(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = false, Int_t switchLayer = 0, Float_t energy = 250, Float_t degraderThickness = 0, Int_t tracksperrun = kEventsPerRun, Bool_t doTracking = kDoTracking);
void			drawAlignmentCheck(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 250);
void			drawDiffusionCheck(Int_t Runs, Int_t Layer, Float_t energy);
void			drawFrame2D(Int_t dataType = kMC, Int_t Layer = 2, Float_t energy = 250, Float_t degraderThickness = 0);
void        drawDataProfile(Float_t energy);
// void        drawTrackingError(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = false, Float_t energy = 250, Float_t degraderThickness = 0);
//
// We want to make binary files each of the functions. First remove all unneccesary files. Then find a naming convention. I suggest:
// <verb><input><details>
// drawData2DFrame // Where Data is either Data or MC !!
// drawDataClusterShapes
// drawDataClusterSizes
// makeTracksFromData
// drawTracksDepthDose --input <file.root> --number N --dataType GATE (maybe standardize? or see from file)
// drawTracksRanges
// drawTracks3D
// savePreconFromTracks
// drawPreconRadiograph
//
#endif
