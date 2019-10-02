#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2C.h>
#include <list>
#include <vector>

using namespace std;

void makeCarbonHeliumBeam() {
   TFile         *fSimulationInHelium = new TFile("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Aluminium_Helium_Absorber3mm_Degrader100mm_917MeV.root", "READ");
   TFile         *fSimulationInCarbon = new TFile("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Aluminium_Carbon_Absorber3mm_Degrader100mm_2751MeV_.root", "READ");
   TFile         *fSimulationOut = new TFile("/home/rttn/workspaceKDevelop/focal/focal/DTCToolkit/Data/MonteCarlo/DTC_Aluminium_CarbonHelium_Absorber3mm_Degrader100mm_917MeV.root", "recreate");
   TTree         *treeSimulationInHelium  = (TTree*) fSimulationInHelium->Get("Hits");
   TTree         *treeSimulationInCarbon  = (TTree*) fSimulationInCarbon->Get("Hits");
   TTree         *treeSimulationOut = new TTree("Hits", "CarbonHelium beam");

   Float_t        x,y,z,edep,outX, outY, outZ, outEdep;
   Int_t          parentID, eventID, outEventID, trackID, outParentID, outTrackID, PDGEncoding, outPDGEncoding;
   Int_t          idxHelium = 0;
   Int_t          idxCarbon = 0;
   
   treeSimulationInHelium->SetBranchAddress("posX",&x);
   treeSimulationInHelium->SetBranchAddress("posY",&y);
   treeSimulationInHelium->SetBranchAddress("posZ",&z);
   treeSimulationInHelium->SetBranchAddress("edep",&edep);
   treeSimulationInHelium->SetBranchAddress("eventID",&eventID);
   treeSimulationInHelium->SetBranchAddress("trackID",&trackID);
   treeSimulationInHelium->SetBranchAddress("parentID",&parentID);
   treeSimulationInHelium->SetBranchAddress("PDGEncoding", &PDGEncoding);
   
   treeSimulationInCarbon->SetBranchAddress("posX",&x);
   treeSimulationInCarbon->SetBranchAddress("posY",&y);
   treeSimulationInCarbon->SetBranchAddress("posZ",&z);
   treeSimulationInCarbon->SetBranchAddress("edep",&edep);
   treeSimulationInCarbon->SetBranchAddress("eventID",&eventID);
   treeSimulationInCarbon->SetBranchAddress("trackID",&trackID);
   treeSimulationInCarbon->SetBranchAddress("parentID",&parentID);
   treeSimulationInCarbon->SetBranchAddress("PDGEncoding", &PDGEncoding);

   treeSimulationOut->Branch("posX", &outX, "posX/F");
   treeSimulationOut->Branch("posY", &outY, "posY/F");
   treeSimulationOut->Branch("posZ", &outZ, "posZ/F");
   treeSimulationOut->Branch("edep",&outEdep, "edep/F");
   treeSimulationOut->Branch("eventID",&outEventID, "eventID/I");
   treeSimulationOut->Branch("trackID",&outTrackID, "trackID/I");
   treeSimulationOut->Branch("parentID",&outParentID, "parentID/I");
   treeSimulationOut->Branch("PDGEncoding", &outPDGEncoding, "PDGEncoding/I");

   Int_t    nHelium = treeSimulationInHelium->GetEntries();
   Int_t    nCarbon = treeSimulationInCarbon->GetEntries();
   Int_t    lastEventIDCarbon = -1;
   Int_t    lastEventIDHelium = -1;
   Bool_t   isHelium = true;

   for (Int_t i=0; i<nHelium+nCarbon; i++) {
      if (idxHelium >= nHelium) break;
      if (idxCarbon >= nCarbon) break;

      if (isHelium) {
         treeSimulationInHelium->GetEntry(idxHelium);

         if (lastEventIDHelium < 0) {
            lastEventIDHelium = eventID;
         }

         if (eventID != lastEventIDHelium) {
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
            
         if (eventID%10 == 0 && lastEventIDCarbon != eventID) {
            isHelium = true;
            treeSimulationInHelium->GetEntry(idxHelium++);
            lastEventIDHelium = eventID;
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
      outEventID = lastEventIDHelium + lastEventIDCarbon;
      outTrackID = trackID;
      outParentID = parentID;
      outPDGEncoding = PDGEncoding;
      treeSimulationOut->Fill();
   }

   treeSimulationOut->Write();
   printf("Merged %d Helium primaries and %d Carbon primaries.\n", lastEventIDHelium, lastEventIDCarbon);
}
