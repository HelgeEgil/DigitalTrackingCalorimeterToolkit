#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2C.h>
#include <list>
#include <vector>

using namespace std;

void makeCarbonHeliumBeamSpots() {

   TFile         *fSimulationInHelium = nullptr;
   TFile         *fSimulationInCarbon = nullptr;


   TFile         *fSimulationOut = new TFile("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Aluminium_CarbonHelium_Absorber3mm_Headphantom_760MeV.root", "recreate");
   TTree         *treeSimulationInHelium = nullptr;
   TTree         *treeSimulationInCarbon = nullptr;
   TTree         *treeSimulationOut = new TTree("Hits", "CarbonHelium beam");

   Float_t        x,y,z,edep,outX, outY, outZ, outEdep, spotPosX, spotPosY;
   Int_t          parentID, eventID, outEventID, trackID, baseID, level1ID, outBaseID, outLevel1ID, outParentID, outTrackID, PDGEncoding, outPDGEncoding;
         
   Int_t    lastEventIDCarbon = -1;
   Int_t    lastEventIDHelium = -1;
   Int_t    totalRunningTally = 0;
   
   treeSimulationOut->Branch("posX", &outX, "posX/F");
   treeSimulationOut->Branch("posY", &outY, "posY/F");
   treeSimulationOut->Branch("posZ", &outZ, "posZ/F");
   treeSimulationOut->Branch("edep",&outEdep, "edep/F");
   treeSimulationOut->Branch("eventID",&outEventID, "eventID/I");
   treeSimulationOut->Branch("trackID",&outTrackID, "trackID/I");
   treeSimulationOut->Branch("parentID",&outParentID, "parentID/I");
   treeSimulationOut->Branch("PDGEncoding", &outPDGEncoding, "PDGEncoding/I");
   treeSimulationOut->Branch("baseID", &outBaseID, "level1ID/I");
   treeSimulationOut->Branch("level1ID", &outLevel1ID, "baseID/I");
   treeSimulationOut->Branch("spotPosX", &spotPosX, "spotPosX/F");
   treeSimulationOut->Branch("spotPosY", &spotPosY, "spotPosY/F");
   
   for (int spotX = -72; spotX <= -30; spotX += 6) {
      for (int spotY = -72; spotY <= 72; spotY += 6) {

         lastEventIDCarbon = -1;
         lastEventIDHelium = -1;

         printf("Running @ spot (%03d, %03d)\n", spotX, spotY);
         fSimulationInHelium = new TFile(Form("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Aluminium_Helium_Absorber3mm_Headphantom_spotsxy_%03d_%03d_760MeV.root", spotX, spotY), "READ");
         fSimulationInCarbon = new TFile(Form("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Aluminium_Carbon_Absorber3mm_Headphantom_spotsxy_%03d_%03d_2280MeV.root", spotX, spotY), "READ");

         treeSimulationInHelium = (TTree*) fSimulationInHelium->Get("Hits");
         treeSimulationInCarbon = (TTree*) fSimulationInCarbon->Get("Hits");

         treeSimulationInHelium->SetBranchAddress("posX",&x);
         treeSimulationInHelium->SetBranchAddress("posY",&y);
         treeSimulationInHelium->SetBranchAddress("posZ",&z);
         treeSimulationInHelium->SetBranchAddress("edep",&edep);
         treeSimulationInHelium->SetBranchAddress("eventID",&eventID);
         treeSimulationInHelium->SetBranchAddress("trackID",&trackID);
         treeSimulationInHelium->SetBranchAddress("parentID",&parentID);
         treeSimulationInHelium->SetBranchAddress("PDGEncoding", &PDGEncoding);
         treeSimulationInHelium->SetBranchAddress("baseID", &baseID);
         treeSimulationInHelium->SetBranchAddress("level1ID", &level1ID);
         
         treeSimulationInCarbon->SetBranchAddress("posX",&x);
         treeSimulationInCarbon->SetBranchAddress("posY",&y);
         treeSimulationInCarbon->SetBranchAddress("posZ",&z);
         treeSimulationInCarbon->SetBranchAddress("edep",&edep);
         treeSimulationInCarbon->SetBranchAddress("eventID",&eventID);
         treeSimulationInCarbon->SetBranchAddress("trackID",&trackID);
         treeSimulationInCarbon->SetBranchAddress("parentID",&parentID);
         treeSimulationInCarbon->SetBranchAddress("PDGEncoding", &PDGEncoding);
         treeSimulationInCarbon->SetBranchAddress("baseID", &baseID);
         treeSimulationInCarbon->SetBranchAddress("level1ID", &level1ID);

         Int_t    nHelium = treeSimulationInHelium->GetEntries();
         Int_t    nCarbon = treeSimulationInCarbon->GetEntries();
         Int_t    idxHelium = 0;
         Int_t    idxCarbon = 0;
         Bool_t   isHelium = true;

         printf("Found %d carbon events and %d helium events in spot\n", nCarbon, nHelium);

         for (Int_t i=0; i<nHelium+nCarbon; i++) {
            if (idxHelium >= nHelium) break;
            if (idxCarbon >= nCarbon) break;

            if (isHelium) {
               treeSimulationInHelium->GetEntry(idxHelium);

               if (lastEventIDHelium < 0) {
                  lastEventIDHelium = eventID;
               }

               if (eventID != lastEventIDHelium) {
                  totalRunningTally++;
                  isHelium = false;
                  lastEventIDCarbon = -1;
               } 
               else {
                  idxHelium++;
                  lastEventIDHelium = eventID;
               }
            }

            if (!isHelium) {
               treeSimulationInCarbon->GetEntry(idxCarbon);

               if (lastEventIDCarbon < 0) {
                  lastEventIDCarbon = eventID;
               }

               if (lastEventIDCarbon != eventID) {
                  totalRunningTally++;
                  
                  if (eventID > lastEventIDHelium*5) { // was *10
                     isHelium = true;
                     treeSimulationInHelium->GetEntry(idxHelium++);
                     lastEventIDHelium = eventID;
                  }
                  else {
                     idxCarbon++;
                     lastEventIDCarbon = eventID;
                  }
               }
               else {
                  idxCarbon++;
                  lastEventIDCarbon = eventID;
               }
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
         }
      }
   }

   fSimulationOut->Write();
   printf("Merged %d Helium primaries and %d Carbon primaries (actual %d carbon events).\n", lastEventIDHelium, lastEventIDCarbon, totalRunningTally-lastEventIDHelium);
}
