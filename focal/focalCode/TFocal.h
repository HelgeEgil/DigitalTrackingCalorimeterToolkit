//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 23 09:22:08 2015 by ROOT version 5.34/24
// from TChain Hits/
//////////////////////////////////////////////////////////

#ifndef Focal_h
#define Focal_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <vector>
#include <TObjArray.h>
#include "Hit.h"
#include "Hits.h"
#include "Cluster.h"
#include "Clusters.h"
#include "Track.h"
#include "Tracks.h"
#include "CalorimeterFrame.h"
#include "TrackerFrame.h"
//#include "ClusterFinder.h"
//#include "TrackFinder.h"
#include "LinkDef.h"

class TH3F;
class TH2F;
class TRandom3;

using namespace std;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Focal {
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
   TBranch        *b_rotationAngle;   //!
   TBranch        *b_volumeID;   //!
   TBranch        *b_processName;   //!
   TBranch        *b_comptVolName;   //!
   TBranch        *b_RayleighVolName;   //!

   Focal(TTree *tree=0);
   virtual ~Focal();

   // Original functions made by GATE / G4
   virtual Int_t     Cut(Long64_t entry);
   virtual Int_t     GetEntry(Long64_t entry);
   virtual Long64_t  LoadTree(Long64_t entry);
   virtual void      Init(TTree *tree);
   virtual Bool_t    Notify();
   virtual void      Show(Long64_t entry = -1);

   // My own functions in TFocal.C
   virtual void getFrame3D(Int_t runNo, CalorimeterFrame *cf);
	//virtual Tracks * FindRealTracks(Int_t Runs, Clusters *restPoints, Int_t energy = -1);
	//virtual TrackerCollection * FindTrackerFrame(Int_t EventFrom, Int_t EventTo, Clusters *restPoints);
	//virtual void		GetRealFrame2D(Int_t Runs, Int_t Layer, TH2F *Frame2D);
   //virtual void      GetData3D(Int_t EventFrom, Int_t EventTo, TH3F *Frame3D);
	//virtual void		GetRealData3D(Int_t EventFrom, Int_t EventTo, TH3F *Frame3D);
	//virtual void		GetTrackerData(Int_t EventFrom, Int_t EventTo, vector<TH2F*> *Frame3D);
	//virtual void		GetClusterFrames(vector<Hits*> *clusterHitMap, Bool_t MC);
};

#endif

#ifdef Focal_cxx
Focal::Focal(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("Hits",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("Hits","");
      chain->Add("test.root/Hits");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Focal::~Focal()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Focal::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Focal::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Focal::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).


   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("PDGEncoding", &PDGEncoding, &b_PDGEncoding);
   fChain->SetBranchAddress("trackID", &trackID, &b_trackID);
   fChain->SetBranchAddress("parentID", &parentID, &b_parentID);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("edep", &edep, &b_edep);
   fChain->SetBranchAddress("stepLength", &stepLength, &b_stepLength);
   fChain->SetBranchAddress("posX", &posX, &b_posX);
   fChain->SetBranchAddress("posY", &posY, &b_posY);
   fChain->SetBranchAddress("posZ", &posZ, &b_posZ);
   fChain->SetBranchAddress("localPosX", &localPosX, &b_localPosX);
   fChain->SetBranchAddress("localPosY", &localPosY, &b_localPosY);
   fChain->SetBranchAddress("localPosZ", &localPosZ, &b_localPosZ);
   fChain->SetBranchAddress("baseID", &baseID, &b_baseID);
   fChain->SetBranchAddress("level1ID", &level1ID, &b_level1ID);
   fChain->SetBranchAddress("level2ID", &level2ID, &b_level2ID);
   fChain->SetBranchAddress("level3ID", &level3ID, &b_level3ID);
   fChain->SetBranchAddress("level4ID", &level4ID, &b_level4ID);
   fChain->SetBranchAddress("layerID", &layerID, &b_layerID);
   fChain->SetBranchAddress("photonID", &photonID, &b_photonID);
   fChain->SetBranchAddress("nPhantomCompton", &nPhantomCompton, &b_nPhantomCompton);
   fChain->SetBranchAddress("nCrystalCompton", &nCrystalCompton, &b_nCrystalCompton);
   fChain->SetBranchAddress("nPhantomRayleigh", &nPhantomRayleigh, &b_nPhantomRayleigh);
   fChain->SetBranchAddress("nCrystalRayleigh", &nCrystalRayleigh, &b_nCrystalRayleigh);
   fChain->SetBranchAddress("primaryID", &primaryID, &b_primaryID);
   fChain->SetBranchAddress("sourcePosX", &sourcePosX, &b_sourcePosX);
   fChain->SetBranchAddress("sourcePosY", &sourcePosY, &b_sourcePosY);
   fChain->SetBranchAddress("sourcePosZ", &sourcePosZ, &b_sourcePosZ);
   fChain->SetBranchAddress("sourceID", &sourceID, &b_sourceID);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("axialPos", &axialPos, &b_axialPos);
   fChain->SetBranchAddress("rotationAngle", &rotationAngle, &b_rotationAngle);
   fChain->SetBranchAddress("volumeID", volumeID, &b_volumeID);
   fChain->SetBranchAddress("processName", processName, &b_processName);
   fChain->SetBranchAddress("comptVolName", comptVolName, &b_comptVolName);
   fChain->SetBranchAddress("RayleighVolName", RayleighVolName, &b_RayleighVolName);
   Notify();
}

Bool_t Focal::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Focal::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Focal::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}



#endif // #ifdef Focal_cxx
