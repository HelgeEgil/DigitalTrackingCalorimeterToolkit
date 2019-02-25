#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <list>
#include <vector>

using namespace std;

void convert_root_file() {
   TFile         *fCluster = new TFile("database_combined_clusters_filter.root", "READ");
   TFile         *fSimulationIn = new TFile("simulationToBeClustered.root", "READ");
   TFile         *fSimulationOut = new TFile("simulationToBeClustered_clustered.root", "recreate");
   TTree         *treeCluster = (TTree*) fCluster->Get("database");
   TTree         *treeSimulationIn  = (TTree*) fSimulationIn->Get("Hits");
   TTree         *treeSimulationOut = new TTree("Hits", "Clustered GATE tree");

   Double_t       x_mean, y_mean;
   Int_t          randomClusterIdx, binPos, idx_x;
   list<Int_t>   *hit_array = 0;
   TBranch       *b_hit_array;
   TRandom3      *gRand = new TRandom3(0);
   Float_t        pixel_size = 0.02929;
   
   Float_t        allZ = 0;
   Int_t          fillN = 0;
   Float_t        x,y,z,edep,outX, outY, outZ;
   Int_t          parentID, eventID, outEventID, outParentID, nClusters, timeInClockSpeed, PDGEncoding;
   Char_t         processName[17];

   TCanvas       *c = new TCanvas();
   TH2C          *hOrig = new TH2C("hOrig", "Original hitmap",1024,-15,15,512,-7.5,7.5);
   TH2C          *hNew  = new TH2C("hNew",  "New hitmap",     1024,-15,15,512,-7.5,7.5);
   Int_t          palette[2] = {kWhite, kBlack};
   
   c->Divide(2,1);
   gStyle->SetPalette(2, palette);

   treeSimulationIn->SetBranchAddress("posX",&x);
   treeSimulationIn->SetBranchAddress("posY",&y);
   treeSimulationIn->SetBranchAddress("posZ",&z);
   treeSimulationIn->SetBranchAddress("edep",&edep);
   treeSimulationIn->SetBranchAddress("eventID",&eventID);
   treeSimulationIn->SetBranchAddress("parentID",&parentID);
   treeSimulationIn->SetBranchAddress("PDGEncoding", &PDGEncoding);

   treeSimulationOut->Branch("posX", &outX, "posX/F");
   treeSimulationOut->Branch("posY", &outY, "posY/F");
   treeSimulationOut->Branch("posZ", &outZ, "posZ/F");
   treeSimulationOut->Branch("eventID", &outEventID, "eventID/I");
   treeSimulationOut->Branch("parentID", &outParentID, "parentID/I");
   treeSimulationOut->Branch("clockTime", &timeInClockSpeed, "clockTime/I");

   treeCluster->SetBranchAddress("x_mean", &x_mean);
   treeCluster->SetBranchAddress("y_mean", &y_mean);
   treeCluster->SetBranchAddress("hit_array", &hit_array, &b_hit_array);

   nClusters = treeCluster->GetEntriesFast();
   randomClusterIdx = gRand->Integer(nClusters);
   treeCluster->GetEntry(randomClusterIdx);

   for (Int_t i=0, N = treeSimulationIn->GetEntries(); i<N; ++i) {
      treeSimulationIn->GetEntry(i);
      
      timeInClockSpeed = 400 * eventID;
      outParentID = parentID;
      outEventID = eventID;
      outZ = z;

      if (z<2) hOrig->Fill(x,y);

      randomClusterIdx = gRand->Integer(nClusters);
      treeCluster->GetEntry(randomClusterIdx);

      idx_x = 0;
      if (abs(PDGEncoding) == 2212) {
         for (Int_t n : *hit_array) { // loop x
            for (int binPosPow = 0; binPosPow < 10; binPosPow++) { // loop y
               binPos = pow(2, binPosPow);
               if (binPos & n) {
                  outX = x + (idx_x - x_mean) * pixel_size;
                  outY = y + (binPosPow - y_mean) * pixel_size;
                  treeSimulationOut->Fill();
                  if (z<2) hNew->Fill(outX, outY);
               }
            }
            idx_x++;
         }
      }
   }

   treeSimulationOut->Write();

   c->cd(1);
   hOrig->Draw("COL");
   c->cd(2);
   hNew->Draw("COL");

   /*
   delete fSimulationIn;
   delete fSimulationOut;
   delete fCluster;
   */
}
