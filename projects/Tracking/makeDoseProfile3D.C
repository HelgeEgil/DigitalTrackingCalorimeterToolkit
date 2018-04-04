#include <TTree.h>
#include <TFile.h>
#include <TH3F.h>
#include <TView.h>
#include <TAxis3D.h>
#include <TPolyMarker3D.h>
#include <cmath>

using namespace std;

void makeDoseProfile3D() {
//   TFile   *f = new TFile("Data/DTC_Aluminium_Absorber3mm_Degrader250mm_250MeV.root");
   TFile   *f = new TFile("Data/RealisticPencilBeamThrough20cmWaterPhantom.root");
   TTree   *tree = (TTree*) f->Get("Hits");
   Float_t  x,y,z,edep;
   Int_t    iret, parentID;
   
   TCanvas *canvas = new TCanvas("canvas", "Dose profile", 1000, 1000);

   Float_t  fromx = -70, tox = 70;
   Float_t  fromy = -70, toy = 70;
   Float_t  fromz = 0, toz = 80;

   Float_t theta = 285, phi = 80;

   TView *view = TView::CreateView(1);
   view->SetRange(fromx, fromz, fromy, tox, toz, toy);
   view->SetView(theta, phi, 0, iret);
   
   Int_t pointStyle = 7, primaryPointIdx = 0, secondaryPointIdx = 0, maxPoints = 20000;
   TPolyMarker3D *primaryPoints = new TPolyMarker3D(maxPoints, pointStyle);
   TPolyMarker3D *secondaryPoints = new TPolyMarker3D(maxPoints, pointStyle);

   secondaryPoints->SetMarkerColor(kRed);

   tree->SetBranchAddress("posX", &x);
   tree->SetBranchAddress("posY", &y);
   tree->SetBranchAddress("posZ", &z);
   tree->SetBranchAddress("edep", &edep);
   tree->SetBranchAddress("parentID", &parentID);

   for (Int_t i=0; i<tree->GetEntries(); ++i) {
      tree->GetEntry(i);
      if (parentID == 0)   primaryPoints->SetPoint(primaryPointIdx++, x, z, y);
      else                 secondaryPoints->SetPoint(secondaryPointIdx++, x, z, y);

      if (primaryPointIdx >= maxPoints || secondaryPointIdx >= maxPoints) break;
   }

   primaryPoints->Draw();
   secondaryPoints->Draw();

   view->ShowAxis();
   TAxis3D *axis = TAxis3D::GetPadAxis();
   axis->SetTitle("3D view of tracks and clusters");
   axis->SetLabelColor(kBlack);
   axis->SetAxisColor(kBlack);
   axis->SetXTitle("X position [mm]");
   axis->SetYTitle("Depth [mm]");
   axis->SetZTitle("Y position [mm]");
   axis->SetTitleOffset(2);
}
