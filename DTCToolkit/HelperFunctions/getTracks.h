#ifndef getTracks_h
#define getTracks_h

#include "Classes/Track/Tracks.h"
#include "Classes/Hit/Hits.h"

// Tracks    * getTracksFOCAL(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188); // See pre-2018-08 version on GH for FOCAL reconstruction methods 
Tracks    * getTracksFromClusters(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188, Float_t spotx = 0, Float_t spoty = 0, Clusters * c = nullptr, Int_t startRun = 0);
Tracks    * getTracksFromClustersVisualize(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188, Float_t spotx = 0, Float_t spoty = 0, Clusters * c = nullptr);
Tracks    * getTracksFromClustersMT(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188, Float_t spotx = 0, Float_t spoty = 0, Clusters * c = nullptr);
Clusters  * getClusters(Int_t Runs, Int_t dataType = kMC, Int_t frameType = kCalorimeter, Float_t energy = 188, Float_t spotx = 0, Float_t spoty = 0);
void        saveTracks(Tracks* tracks, Clusters * clusters, Int_t dataType, Float_t energy);
Tracks    * loadTracks(Int_t Runs, Int_t dataType, Float_t energy, Clusters *clusters);
Tracks    * loadOrCreateTracks(Bool_t recreate, Int_t Runs, Int_t dataType, Float_t energy, Float_t spotx = 0, Float_t spoty = 0, Clusters *c = nullptr, Int_t startRun = 0);
Hits      * diffuseHits(TRandom3 *gRandom, Hits * hits);
Hits      * diffuseHitsMT(Hits * hits);
#endif
