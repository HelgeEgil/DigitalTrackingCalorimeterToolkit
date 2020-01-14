#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2C.h>
#include <list>
#include <vector>

using namespace std;

void combineSpots() {

   TFile         *fSimulationIn = nullptr;
   TFile         *fSimulationOut = new TFile("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Final_HeadPhantom_rotation90deg.root", "recreate");
   TTree         *treeSimulationIn = nullptr;
   TTree         *treeSimulationOut = new TTree("Hits", "Combined spots");

   Float_t        x,y,z,edep,outX, outY, outZ, outEdep, spotPosX, spotPosY;
   Int_t          parentID, eventID, outEventID, trackID, baseID, level1ID, outBaseID, outLevel1ID, outParentID, outTrackID, PDGEncoding, outPDGEncoding;
    
   Int_t    lastEventID = -1;
   Int_t    totalRunningTally = 0;
   
   treeSimulationOut->Branch("posX", &outX, "posX/F");
   treeSimulationOut->Branch("posY", &outY, "posY/F");
   treeSimulationOut->Branch("posZ", &outZ, "posZ/F");
   treeSimulationOut->Branch("edep",&outEdep, "edep/F");
   treeSimulationOut->Branch("eventID",&outEventID, "eventID/I");
   treeSimulationOut->Branch("parentID",&outParentID, "parentID/I");
   treeSimulationOut->Branch("PDGEncoding", &outPDGEncoding, "PDGEncoding/I");
   treeSimulationOut->Branch("baseID", &outBaseID, "level1ID/I");
   treeSimulationOut->Branch("level1ID", &outLevel1ID, "baseID/I");
   treeSimulationOut->Branch("spotPosX", &spotPosX, "spotPosX/F");
   treeSimulationOut->Branch("spotPosY", &spotPosY, "spotPosY/F");
  
   for (int spotX = -100; spotX <= 100; spotX += 100) {
      for (int spotY = -100; spotY <= 100; spotY += 10) {

         printf("Running @ spot (%03d, %03d)\n", spotX, spotY);
         fSimulationIn = new TFile(Form("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Final_HeadPhantom_spotx%04d_spoty%04d_rotation90deg.root", spotX, spotY), "READ");

         treeSimulationIn = (TTree*) fSimulationIn->Get("Hits");

         treeSimulationIn->SetBranchAddress("posX",&x);
         treeSimulationIn->SetBranchAddress("posY",&y);
         treeSimulationIn->SetBranchAddress("posZ",&z);
         treeSimulationIn->SetBranchAddress("edep",&edep);
         treeSimulationIn->SetBranchAddress("eventID",&eventID);
         treeSimulationIn->SetBranchAddress("parentID",&parentID);
         treeSimulationIn->SetBranchAddress("PDGEncoding", &PDGEncoding);
         treeSimulationIn->SetBranchAddress("baseID", &baseID);
         treeSimulationIn->SetBranchAddress("level1ID", &level1ID);
         
         Int_t    nIn = treeSimulationIn->GetEntries();
         printf("Found %d events in spot\n", nIn);
         lastEventID = -1;

         for (Int_t i=0; i<nIn; i++) {
            treeSimulationIn->GetEntry(i);

            if (lastEventID < 0) {
               lastEventID = eventID;
            }

            if (eventID != lastEventID) {
               totalRunningTally++;
            } 

            outX = x;
            outY = y;
            outZ = z;
            outEdep = edep;
            outEventID = totalRunningTally;
            outTrackID = trackID;
            outParentID = parentID;
            outLevel1ID = level1ID;
            outBaseID = baseID;
            outPDGEncoding = PDGEncoding;
            spotPosX = (float) spotX;
            spotPosY = (float) spotY;

            treeSimulationOut->Fill();

            lastEventID = eventID;
         }
      }
   }

   fSimulationOut->Write();
   printf("Merged %d primaries.\n", totalRunningTally);
}
