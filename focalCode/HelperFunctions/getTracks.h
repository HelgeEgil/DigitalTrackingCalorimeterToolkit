#ifndef getTracks_h
#define getTracks_h

#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"

Tracks * getTracks(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188, Float_t *x = 0, Float_t *y = 0);
void saveTracks(Tracks* tracks, Int_t dataType, Float_t energy);
Tracks* loadTracks(Int_t Runs, Int_t dataType, Float_t energy);
Tracks * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy, Float_t *x = 0, Float_t *y = 0);
#endif
