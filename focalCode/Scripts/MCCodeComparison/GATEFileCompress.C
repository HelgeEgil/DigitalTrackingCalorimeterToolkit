#include <TSystem.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>

using namespace std;

void Run()
{
   Float_t  lastZ = 0;
   Float_t  lastX = 0, lastY = 0;
   Int_t    lastEventID = -1, lastParentID = -1;
   Bool_t   lastProcessWasInelastic = false;
   Int_t    nominalEnergy;

   Int_t nMax = 1000;

//   TFile   *f = new TFile("Data/GATE/complex_geometry_emy_nolimit.root");
//   TFile   *fOut = new TFile("Data/GATE/compressed_complex_geometry_emy_nolimit.root", "recreate");

   for (Int_t j=0; j<19; j++) {
      nominalEnergy = (j+5) * 10;

      cout << "Compressing tree with " << nominalEnergy << " MeV -> ";
//      TFile *f = new TFile(Form("Data/GATE/ComplexGeometry/complex_%dMeV.root", nominalEnergy));
//      TFile *fOut = new TFile(Form("Data/GATE/ComplexGeometry/compressed_complex_%dMeV.root", nominalEnergy), "recreate");
      TFile *f = new TFile(Form("Data/GATE/Water/water_71_%dMeV.root", nominalEnergy));
      TFile *fOut = new TFile(Form("Data/GATE/Water/compressed_water_71_%dMeV.root", nominalEnergy), "recreate");
//      TFile *f = new TFile(Form("Data/GATE/Aluminium/aluminium_%dMeV.root", nominalEnergy));
//      TFile *fOut = new TFile(Form("Data/GATE/Aluminium/compressed_aluminium_%dMeV.root", nominalEnergy), "recreate");

      TTree   *tree = (TTree*) f->Get("Hits");
      TTree    treeOut("treeOut", "Compressed GATE tree");

      Float_t  allZ = 0;
      Int_t    fillN = 0;
      Float_t  x,y,z,edep, lastEdep;
      Int_t    parentID, eventID;
      Char_t   processName[17];

      tree->SetBranchAddress("posX",&x);
      tree->SetBranchAddress("posY",&y);
      tree->SetBranchAddress("posZ",&z);
      tree->SetBranchAddress("edep",&edep);
      tree->SetBranchAddress("eventID",&eventID);
      tree->SetBranchAddress("parentID",&parentID);
      tree->SetBranchAddress("processName",processName);

      treeOut.Branch("posX", &lastX, "posX/F");
      treeOut.Branch("posY", &lastY, "posY/F");
      treeOut.Branch("posZ", &lastZ, "posZ/F");
      treeOut.Branch("edep", &lastEdep, "edep/F");
      treeOut.Branch("eventID", &lastEventID, "eventID/I");
      treeOut.Branch("parentID", &lastParentID, "parentID/I");
      treeOut.Branch("isInelastic", &lastProcessWasInelastic, "isInelastic/O");

      for (Int_t i=0, N = tree->GetEntries(); i<N; ++i) {
         tree->GetEntry(i);

         if (parentID == 0) {
            if (lastEventID < 0) lastEventID = eventID;

            if (processName[0] == 'p') { // ProtonInelastic
               lastX = x;
               lastY = y;
               lastZ = z;
               lastEventID = eventID;
               lastParentID = parentID;
               lastEdep = edep;
               lastProcessWasInelastic = true;
               treeOut.Fill();
               lastEventID = -1; // We don't want these processes to be counted twice; here and next time when a new event ID is triggered
               lastProcessWasInelastic = false;
               continue;
            }

            if (eventID != lastEventID) { // new event ID - we know last interaction was the history's last
               allZ += lastZ; fillN++;
               treeOut.Fill();
            }

            lastX = x;
            lastY = y;
            lastZ = z;
            lastEventID = eventID;
            lastParentID = parentID;
            lastEdep = edep;
         }
      }

      // STORE LAST EVENT
      treeOut.Fill();
      treeOut.Write();

      delete f;
      delete fOut;

      cout << "range = " << allZ/fillN << " mm.\n";
   }
}
