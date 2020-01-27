#ifndef DataInterface_h
#define DataInterface_h

#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TObjArray.h>


#include "GlobalConstants/Constants.h"
#include "GlobalConstants/RangeAndEnergyCalculations.h"
#include "Classes/Hit/Hit.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Track/Track.h"
#include "Classes/Track/Tracks.h"
// #include "Classes/Track/conversionFunctions.h"
#include "HelperFunctions/Tools.h"
#include "Classes/Layer/Layer.h"
#include "Classes/CalorimeterFrame/CalorimeterFrame.h"
#include "RootFiles/LinkDef.h"


class TH3F;
class TH2F;
class TRandom3;

using namespace std;

namespace DTC {
class Hits;

class DataInterface {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           PDGEncoding;
   Int_t           trackID;
   Int_t           parentID;
   Double_t        time;
   Float_t         edep;
   Float_t         stepLength;
   Float_t         posX;
   Float_t         posY;
   Float_t         posZ;
   Float_t         localPosX;
   Float_t         localPosY;
   Float_t         localPosZ;
   Int_t           baseID;
   Int_t           level1ID;
   Int_t           level2ID;
   Int_t           level3ID;
   Int_t           level4ID;
   Int_t           layerID;
   Int_t           photonID;
   Int_t           nPhantomCompton;
   Int_t           nCrystalCompton;
   Int_t           nPhantomRayleigh;
   Int_t           nCrystalRayleigh;
   Int_t           primaryID;
   Float_t         sourcePosX;
   Float_t         sourcePosY;
   Float_t         sourcePosZ;
   Float_t         spotPosX;
   Float_t         spotPosY;
   Int_t           sourceID;
   Int_t           eventID;
   Int_t           runID;
   Float_t         axialPos;
   Float_t         rotationAngle;
   Int_t           volumeID[10];
   Char_t          processName[15];
   Char_t          comptVolName[5];
   Char_t          RayleighVolName[5];

   // List of branches
   TBranch        *b_PDGEncoding;   //!
   TBranch        *b_trackID;   //!
   TBranch        *b_parentID;   //!
   TBranch        *b_time;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_stepLength;   //!
   TBranch        *b_posX;   //!
   TBranch        *b_posY;   //!
   TBranch        *b_posZ;   //!
   TBranch        *b_localPosX;   //!
   TBranch        *b_localPosY;   //!
   TBranch        *b_localPosZ;   //!
   TBranch        *b_baseID;   //!
   TBranch        *b_level1ID;   //!
   TBranch        *b_level2ID;   //!
   TBranch        *b_level3ID;   //!
   TBranch        *b_level4ID;   //!
   TBranch        *b_layerID;   //!
   TBranch        *b_photonID;   //!
   TBranch        *b_nPhantomCompton;   //!
   TBranch        *b_nCrystalCompton;   //!
   TBranch        *b_nPhantomRayleigh;   //!
   TBranch        *b_nCrystalRayleigh;   //!
   TBranch        *b_primaryID;   //!
   TBranch        *b_sourcePosX;   //!
   TBranch        *b_sourcePosY;   //!
   TBranch        *b_sourcePosZ;   //!
   TBranch        *b_sourceID;   //!
   TBranch        *b_eventID;   //!
   TBranch        *b_runID;   //!
   TBranch        *b_axialPos;   //!
   TBranch        *b_spotPosX;
   TBranch        *b_spotPosY;
   TBranch        *b_rotationAngle;   //!
   TBranch        *b_volumeID;   //!
   TBranch        *b_processName;   //!
   TBranch        *b_comptVolName;   //!
   TBranch        *b_RayleighVolName;   //!

   Long64_t primariesInSpot_;

   DataInterface(TTree *tree=0);
   virtual ~DataInterface();

   // Original functions made by GATE / G4
   virtual Int_t     Cut(Long64_t entry);
   virtual Int_t     GetEntry(Long64_t entry);
   virtual Long64_t  LoadTree(Long64_t entry);
   virtual void      Init(TTree *tree);
   virtual Bool_t    Notify();
   virtual void      Show(Long64_t entry = -1);

   // My own functions in DataInterface.C
   virtual Int_t     getMCFrame(Int_t runNo, CalorimeterFrame *cf);
   virtual Int_t     getMCFrame(Int_t runNo, Layer *l);
   virtual void      getMCClusters(Int_t runNo, Clusters *clusters = nullptr, Hits * hits = nullptr, Float_t spotPosX = 0, Float_t spotPosY = 0);
   virtual void      getMCClustersThreshold(Int_t runNo, Clusters *clusters = nullptr, Hits * hits = nullptr, Float_t spotPosX = 0, Float_t spotPosY = 0);
   virtual void      getDataFrame(Int_t runNo, CalorimeterFrame *cf, Int_t energy = 188);
   virtual void      getDataFrame(Int_t runNo, Layer *l, Int_t energy = 188);
   virtual void      getDataHits(Int_t runNo, Hits * hits, Int_t energy = 188);
//   virtual void      writeDataFrame(Int_t energy = 190); // Removed, se GH pre-2018-08
   virtual void      getMCData(Int_t runNo, TH3F* Frame3D);
   virtual void      getDataProfile(TH2F *hProfile, TH2F *hProjection, Int_t energy);
   virtual void      getEventIDs(Int_t runNo, Hits* hits);
   virtual Long64_t  findSpotIndex(Float_t useSpotX);
};
}
#endif

