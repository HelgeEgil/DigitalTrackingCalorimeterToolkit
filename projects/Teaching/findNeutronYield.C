#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TStyle.h>

using namespace std;

void Run() {
   TFile  * f = new TFile("output/neutron_yield_W.root");
   TTree  * tree = (TTree*) f->Get("Hits");
   Float_t  z,edep,lastZ = -1, dE = 0;
   Int_t    eventID, parentID, pid, tid, nNeutrons = 0, lastEventID = -1, lastTrackID = -1;

   gStyle->SetTitleFont(22);
   
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("eventID", &eventID);
   tree->SetBranchAddress("parentID", &parentID);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("trackID", &tid);
   tree->SetBranchAddress("PDGEncoding", &pid);

   printf("Found %d entries.\n", tree->GetEntries());

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      
      if (parentID != 0) {
         if (pid == 2212) {
            if (lastTrackID != tid || lastEventID != eventID) {
               printf("Following neutron with event ID %d, trackID = %d and parentID %d\n", eventID, tid, parentID);
               nNeutrons++;
            }
         }
      }

      lastTrackID = tid;
      lastEventID = eventID;
   }

   printf("In total %d neutrons produced from %d primaries (yield = %.3e).\n", nNeutrons, lastEventID, float(nNeutrons) / float(lastEventID));
}
