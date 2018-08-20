#ifndef getTracks_h
#define getTracks_h

#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"

// Tracks    * getTracksFOCAL(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188); // See pre-2018-08 version on GH for FOCAL reconstruction methods 
Tracks    * getTracksFromClusters(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188);
Clusters  * getClusters(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188);
void        saveTracks(Tracks* tracks, Int_t dataType, Float_t energy);
Tracks    * loadTracks(Int_t Runs, Int_t dataType, Float_t energy);
Tracks    * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy);
#endif
