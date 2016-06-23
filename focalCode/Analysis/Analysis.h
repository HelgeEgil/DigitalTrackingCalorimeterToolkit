#ifndef Analysis_h
#define Analysis_h

class TCanvas;
class TH1F;
class TF1;

#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"

Tracks * getTracks(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188, Float_t *x = 0, Float_t *y = 0);
Float_t 		drawBraggPeakGraphFit(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188);
void	drawTungstenSpectrum(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188);
void saveTracks(Tracks* tracks, Int_t dataType, Float_t energy);
Tracks* loadTracks(Int_t Runs, Int_t dataType, Float_t energy);
Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy, Float_t *x = 0, Float_t *y = 0);
Hits * getEventIDs(Int_t Runs, Float_t energy);
void		drawFrame2D(Int_t Runs, Int_t Layer, Float_t energy);
void		drawDiffusionCheck(Int_t Runs, Int_t Layer, Float_t energy);
void		drawTracks3D(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = false, Float_t energy = 188);
void 		drawTrackAngleAtVaryingRunNumbers(Int_t dataType, Float_t energy);
void		drawData3D(Int_t Runs, Float_t energy);
void		writeClusterFile(Int_t Runs, Int_t dataType, Float_t energy);
void		drawClusterShapes(Int_t Runs, Bool_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188);
void 		draw2DProjection(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy);
void 		drawTrackRanges(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy);
void		getTrackStatistics(Int_t Runs, Int_t dataType = kMC, Bool_t recreate = 0, Float_t energy = 188, Int_t epr = 0);
void drawFitScale(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy);
Bool_t getCutTrackLength(Float_t energy, Track *track);
Bool_t getCutWEPL(Track *track);
Bool_t getCutChipNumber(Track *track);
Bool_t getCutBraggPeakInTrack(Track *track);
void drawIndividualGraphs(TCanvas *cGraph, TGraphErrors* outputGraph, Float_t fitEnergy, Float_t fitScale, Float_t fitError, Int_t fitIdx, Int_t skipIdx = 0, Float_t *x = 0, Float_t *y = 0);
Float_t doNGaussianFit( TH1F *h, Float_t *means, Float_t *sigmas);

#endif
