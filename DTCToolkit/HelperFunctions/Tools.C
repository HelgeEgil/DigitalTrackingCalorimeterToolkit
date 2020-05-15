#ifndef Tools_cxx
#define Tools_cxx

#include <vector>
#include <algorithm>

#include <TObject.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TLine.h>
#include <TCanvas.h>

#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "GlobalConstants/RangeAndEnergyCalculations.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Hit/Hit.h"
// #include "Classes/Track/conversionFunctions.h"


using namespace std;

Bool_t isItemInVector(Int_t i, vector<Int_t> *v) {
   return find(v->begin(), v->end(), i) != v->end();
}

void BinLogY(TH1 *h) {
   TAxis *axis = h->GetYaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i=0; i<=bins; i++) {
      new_bins[i] = pow(10, from + i*width);
   }
   axis->Set(bins, new_bins);
   delete[] new_bins;
}

Bool_t existsEnergyFile(Int_t energy) {
   Bool_t res = false;
   for (Int_t i=0; i<nEnergies; i++) {
      if (energies[i] == energy) res = true;
   }
   return res;
}

Float_t diffXY(Cluster *p, Hit *h) {
   if (!p || !h) return 1e5;

   Float_t diffx = p->getX() - h->getX();
   Float_t diffy = p->getY() - h->getY();

   return sqrt(pow(diffx, 2) + pow(diffy, 2));
}

Float_t diffXY(Cluster *p1, Cluster *p2) {
   if (!p1 || !p2)
      return -1;

   Double_t diffx = p2->getX() - p1->getX();
   Double_t diffy = p2->getY() - p1->getY();

   return sqrt(pow(diffx,2) + pow(diffy,2));
}

Float_t diffmmXY(Cluster *p1, Cluster *p2) {
   if (!p1 || !p2)
      return -1;

   Double_t diffx = p2->getXmm() - p1->getXmm();
   Double_t diffy = p2->getYmm() - p1->getYmm();

   return sqrt(pow(diffx,2) + pow(diffy,2));
}

Float_t diffmmXY(Hit *h1, Hit *h2) {
   if (!h1 || !h2)
      return -1;

   Double_t diffx = h2->getXmm() - h1->getXmm();
   Double_t diffy = h2->getYmm() - h1->getYmm();

   return sqrt(pow(diffx,2) + pow(diffy,2));
}

Float_t diffmmXYZ(Cluster *p1, Cluster *p2) {
   return sqrt(pow(p2->getXmm() - p1->getXmm(), 2) + 
            pow(p2->getYmm() - p1->getYmm(), 2) +
            pow(p2->getLayermm() - p1->getLayermm(), 2));
}

Double_t fitfunc_DBP(Double_t *v, Double_t *par) {
   // See PRC publication (H. Pettersen, 2017) for details about the accuracy of this approach

   Float_t depth = v[0];
   Float_t range = par[0];
   Float_t scale = par[1];

   Double_t fitval = 0;
   Float_t aa=alpha, pp=p;

//   if (depth > range) return 0;
   fitval = scale / pow(range - depth, 1-1/pp);

   if (std::isnan(fitval)) fitval = 0;

   return fitval;
}

Double_t double_landau(Double_t *v, Double_t *par) {
   Double_t x = v[0];
   Double_t sigma1 = par[1];
   Double_t sigma2 = par[2];
   Double_t mpv1 = par[0];
   Double_t mpv2 = mpv1 + dz; 
   if (kOutputUnit == kWEPL || kOutputUnit == kUnitEnergy)
      mpv2 = mpv1 + getWEPLFactorFromEnergy(run_energy) * dz;
   
   return TMath::Landau(x, mpv1, sigma1, 0) + TMath::Landau(x, mpv2, sigma2, 0);
}

char * getMaterialChar() {
   char *sMaterial;

   if (kMaterial == kTungsten)   sMaterial = (char*) "Tungsten";
   if (kMaterial == kAluminum)   sMaterial = (char*) "Aluminium";
   if (kMaterial == kPMMA)       sMaterial = (char*) "PMMA";
   if (kMaterial == kCarbon)     sMaterial = (char*) "Carbon";

   return sMaterial;
}

char * getMaterialAbbr() {
   char *sMaterial;

   if (kMaterial == kTungsten)   sMaterial = (char*) "W";
   if (kMaterial == kAluminum)   sMaterial = (char*) "Al";
   if (kMaterial == kPMMA)       sMaterial = (char*) "PMMA";
   if (kMaterial == kCarbon)     sMaterial = (char*) "C";

   return sMaterial;
}

char * getDataTypeChar(Int_t dataType) {
   char * sDataType = (char*) "";
   if       (dataType == kMC)    sDataType = (char*) "MC";
   else if  (dataType == kData)  sDataType = (char*) "Exp. data";

   return sDataType;
}

Int_t getFWxMInRange(TH1F* h, Float_t first, Float_t last, Int_t div) {
   // find the full with at maximum / div (2 = FWHM, 3 = FWTM, etc.)
   
   Int_t firstRangeBin = h->GetBin(first);
   Int_t lastRangeBin = h->GetBin(last);
   
// cout << "firstRangeBin = " << firstRangeBin << ", lastRangeBin = " << lastRangeBin << endl;
   
   Float_t maxBin = 0;
   Int_t maxIdx = 0;
   for (Int_t i=firstRangeBin; i<=lastRangeBin; i++) {
//       cout << h->GetBinContent(i) << endl;
      if (h->GetBinContent(i) > maxBin ) {
         maxBin = h->GetBinContent(i);
//          cout << "new max = " << maxBin << endl;
         maxIdx = i;
      }
   }
   
//    cout << "max value is " << maxBin << " at " << maxIdx << endl;
   
   Float_t firstCenter = 0, lastCenter = 0;
   for (Int_t i=firstRangeBin; i<=lastRangeBin; i++) {
      if (!firstCenter) {
         if (h->GetBinContent(i) > maxBin/div) {
            firstCenter = h->GetBinCenter(i);
//             cout << "firstCenter at " << firstCenter << endl;
         }
      }
      
      else {
         if (!lastCenter) {
            if (h->GetBinContent(i) < maxBin/div) {
               lastCenter = h->GetBinCenter(i-1);
//                cout << "lastCenter at " << firstCenter << endl;
            }
         }
      }
   }
   
   if (!firstCenter || !lastCenter) return 0;   
   Float_t width = lastCenter - firstCenter;

   return width * 2 * sqrt( 2 * log(div));
}

Int_t getMinimumTrackLength(Float_t energy) {
   Int_t minTL = 0;

   if (energy < 150) minTL = 5;
   else if (energy < 170) minTL = 10;
   else if (energy < 190) minTL = 15;
   else if (energy < 200) minTL = 20;
   else if (energy < 230) minTL = 23;

   return minTL;
}

Float_t quadratureAdd(Float_t a, Float_t b) {
   return sqrt(pow(a,2) + pow(b, 2));
}

void convertXYToWEPL(Float_t *x_energy, Float_t *y_energy, Int_t eventID) {
   Float_t WEPLFactor = getWEPLFactorFromEnergy(run_energy);
   Long64_t n=0;
   Long64_t j=0;

   for (Long64_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) {
      x_energy[i] *=  WEPLFactor;
   }
}

Float_t getEnergyFromXY(Float_t *x_energy, Float_t *y_energy, Int_t eventID) {
   Float_t WEPLFactor = getWEPLFactorFromEnergy(run_energy);
   Long64_t n=0;
   Long64_t j=0;

   for (Long64_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) {
      if (y_energy[i] == 0) {
         n = j;
         break;
      }
      j++;
   }

   return getEnergyFromWEPL(x_energy[eventID*sizeOfEventID + n-1]);
}

Double_t correctForEnergyParameterisation(Float_t energy) {
   // By using the range = alpha * pow ( energy, p ) parameterisation to find
   // the depth dose curve used in Bragg Peak fitting, a 0.2 % error is introduced
   // This may be fixed by calculating back the WEPL range using the r = a E^p method,
   // and then finding the corresponding energy using the cubic and more accurate method

   Double_t wepl = alpha_water * pow(energy, p_water);
   Double_t corrected_energy = getEnergyFromWEPL(wepl);

   return corrected_energy;
}

Float_t calculateMCS(Float_t energy, Float_t depth, Bool_t kUseDifferentX0) {
   Float_t gamma = (energy + proton_mass) / proton_mass;
   Float_t beta = sqrt(1 - pow(gamma, -2));
   Float_t momentum = gamma * beta * proton_mass;

   Float_t thisX0 = X0;
   if (kUseDifferentX0) {
      thisX0 = X0_firstlayer;
   }

   Float_t mcs = 13.6 / (beta * momentum) * sqrt(depth / thisX0) * (1 + 0.038 * log(depth / thisX0));

   return mcs; // radians
}

Float_t findMCSPixelRadiusForLayer(Float_t layer, Float_t E0) {
   Float_t currentEnergy = getEnergyAtTL(E0, getLayerPositionmm(layer));
   currentEnergy -= getAverageEnergyLoss(E0);
   
   Bool_t kUseDifferentX0 = false;

   Float_t depth = dz;
   if (layer == 0) {
      // Expected depth of scintillator + preabsorber (and no tungsten absorber before sensor)
      depth = 16.6;
      kUseDifferentX0 = true;
   }
   
   Float_t restRange = getTLFromEnergy(currentEnergy);
   if (restRange < depth) depth = restRange;

   Float_t mcs = calculateMCS(currentEnergy, depth, kUseDifferentX0);
   Float_t radius = depth * tan (mcs);

   return radius;
}

Float_t findMCSAtLayerRad(Int_t layer, Float_t E0) {
   Float_t energy = E0 - getAverageEnergyLoss(E0);
   Float_t restRange = getTLFromEnergy(energy);

   Float_t mcs = calculateMCS(energy, restRange, false);
   return mcs;
}

void fillMCSRadiusList(Float_t angleFactor) {
   Float_t pixelRadius;

   if (!run_energy) {
      cout << "\033[1mRUN_ENERGY NOT DEFINED! ASSUMING 180 MeV.\033[0m\n";
   }

   for (Int_t i=0; i<nLayers; i++) {
      pixelRadius = findMCSPixelRadiusForLayer(i, run_energy);
      if (!std::isnan(pixelRadius)) {
         mcs_radius_per_layer[i] = pixelRadius * angleFactor;
      }
   }
}

void multiplyRadiusFirstLayers(Float_t factor) {
   mcs_radius_per_layer[0] *= factor;
   mcs_radius_per_layer[1] *= factor;
}

Float_t getSearchRadiusForLayer(Int_t layer) {
   return mcs_radius_per_layer[layer];
}

Float_t getMCSAngleForLayer(Int_t layer) {
   Float_t radius = getSearchRadiusForLayer(layer);
   Float_t depth = dz;
   if (layer == 0) {
      depth = 16.6;
   }
   
   Float_t angle = atan( radius / depth );
   return angle * 180 / 3.14159265;
}

Float_t getEmpiricalMCSAngle(Int_t layer) {
   // defined as the dot product rule angle between the ingoing and outgoing vector in mm
   // This value is found through Analysis.C::findMCSAngles.

//   return kMCSFactor * 3 * mcs_radius_per_layer_empirical[layer];
   return kMCSFactor * 0.05;
}

Float_t getDotProductAngle(Cluster *a, Cluster *b, Cluster *c) {
   // If a == b, b is usually the first layer, but we want to parallel propagate the vector to calculate the change from ||
   // If not, the distance between the layers might be > dz, so calculate it.

//   Double_t in[3]  = {b->getXmm() - a->getXmm(), b->getYmm() - a->getYmm(), (a==b) ? dz : b->getLayermm() - a->getLayermm()};
//   Double_t out[3] = {c->getXmm() - b->getXmm(), c->getYmm() - b->getYmm(), c->getLayermm() - b->getLayermm()};

   Double_t in[3]  = {b->getXmm() - a->getXmm(), b->getYmm() - a->getYmm(), b->getLayermm() - a->getLayermm()};
   Double_t out[3] = {c->getXmm() - b->getXmm(), c->getYmm() - b->getYmm(), c->getLayermm() - b->getLayermm()};
   if (a == b) in[2] = out[2]; // NOT RIGHT WHEN USING TRACKER LAYERS
   if (b == c) out[2] = in[2];

   Double_t dot    = in[0] * out[0] + in[1] * out[1] + in[2] * out[2];
   Double_t scalar = sqrt(pow(out[0],2) + pow(out[1],2) + pow(out[2],2)) * sqrt(pow(in[0],2) + pow(in[1],2) + pow(in[2],2));
   
   if (std::isnan(acos(dot / scalar))) {
      if (!std::isnan(dot) && !std::isnan(scalar)) return 0; // floating point error on almost || track
   }

   return acos(dot / scalar);
}

Float_t getAverageEnergyLoss(Float_t energy) {
   Float_t dE_1 = getEnergyLossFromScintillators(energy, 1);
   Float_t dE_2 = getEnergyLossFromScintillators(energy, 2);
   Float_t dE_3 = getEnergyLossFromScintillators(energy, 3);
   Float_t dE_Al = getEnergyLossFromAluminumAbsorber(energy);

   Float_t probabilityOf1Scintillator = 0.5625;
   Float_t probabilityOf2Scintillators = 0.375;
   Float_t probabilityOf3Scintillators = 0.0675;

   Float_t avgEnergyLoss = dE_1 * probabilityOf1Scintillator +
                           dE_2 * probabilityOf2Scintillators +
                           dE_3 * probabilityOf3Scintillators +
                           dE_Al;

   return avgEnergyLoss;
}

Hit * sumHits(Hit * a, Hit * b) {
   Int_t x, y, layer, eventID;
   Float_t eDep;

   if (!a) cout << "!a\n";
   if (!b) cout << "!b\n";

   eDep = a->getEdep() + b->getEdep(); // sum 
   x = ( (int) a->getX() + (int) b->getX() + 0.5 ) / 2; // average
   y = ( (int) a->getY() + (int) b->getY() + 0.5 ) / 2; // average
   layer = a->getLayer();
   eventID = a->getEventID();

   if (a->getEventID() != b->getEventID() || a->getLayer() != b->getLayer()) {
      cout << "\033[1mASSERTION ERROR IN sumHits: THE EVENT IDS OR LAYER# ARE DIFFERENT!!\033[0m\n";
   }

   Hit * newHit = new Hit(x, y, layer, eventID, eDep);

   return newHit;
}

Cluster * getTrackExtrapolationToLayer(Track * track, Int_t layer) {
   if (track->GetEntriesFast() == 1) return nullptr;

   Int_t fromLayer = track->getLayer(track->GetEntriesFast() - 1);
   showDebug("Cluster::getTrackExtrapolationToLayer(" << track <<"," << layer << "): fromLayer = " << fromLayer << endl);

   return getTrackExtrapolationFromTo(track, fromLayer, layer);
}

Cluster * getRetrogradeTrackExtrapolationToLayer(Track * track, Int_t layer) {
   Int_t fromLayer = layer + 1;

   if (track->getIdxFromLayer(fromLayer) < 0) fromLayer++;
   if (track->getIdxFromLayer(fromLayer) < 0) return nullptr;
   
   return getTrackExtrapolationFromTo(track, fromLayer, layer);
}

Cluster  * getTrackExtrapolationCluster(Cluster *p1, Cluster *p2) {
   Int_t fromLayer = p2->getLayer();
   Int_t toLayer = p2->getLayer() + 1;
   Int_t diffLayer = toLayer - fromLayer;

   Cluster slope((p2->getX() - p1->getX()) / diffLayer, (p2->getY() - p1->getY()) / diffLayer);

   Float_t x = p2->getX() + diffLayer * slope.getX();
   Float_t y = p2->getY() + diffLayer * slope.getY();

   return new Cluster(x,y,toLayer);
}

Cluster * getTrackExtrapolationFromTo(Track * track, Int_t fromLayer, Int_t toLayer) {
   Int_t diffLayer = toLayer - fromLayer;
   Int_t sign = (0 < diffLayer) - (diffLayer < 0);
   Int_t vectorLayer = fromLayer - sign;

   Int_t fromIdx = track->getIdxFromLayer(fromLayer);
   Int_t vectorIdx = fromIdx - sign;
   if (!track->At(vectorIdx)) return nullptr;
   
   Int_t vectorDiff = abs(fromLayer - track->getLayer(vectorIdx));
   
   Cluster p1(track->getX(vectorIdx), track->getY(vectorIdx));
   Cluster p2(track->getX(fromIdx), track->getY(fromIdx));

   Cluster slope((p2.getX() - p1.getX()) / vectorDiff, (p2.getY() - p1.getY()) / vectorDiff);

   Float_t x = p2.getX() + diffLayer * slope.getX();
   Float_t y = p2.getY() + diffLayer * slope.getY();

   return new Cluster(x, y, toLayer);
}

Bool_t isPointOutOfBounds(Cluster *point, Float_t padding) {
   Bool_t isOutside;
   Float_t x = point->getX();
   Float_t y = point->getY();

   if (!point)
      isOutside = kTRUE;
   else {
      if (x < 0 - padding || x > nx + padding ||
          y < 0 - padding || y > ny + padding) {
         isOutside = kTRUE;
      }
      else
         isOutside = kFALSE;
   }

   return isOutside;
}

Float_t getEdepFromCS(Int_t cs) {
   Float_t edep;

   if (kUseExperimentalClustering)  edep = 0.10890 * pow(cs, 1.5384); // kev/um
   else                             edep = 0.13345 * pow(cs, 1.06576);

   return edep;
}

Int_t getCSFromEdep(Float_t edep) {
   Int_t cs;

   if (kUseExperimentalClustering)  cs = 4.2267 * pow(edep, 0.6500) + 0.5; // kev/um
   else                             cs = 6.6177 * pow(edep, 0.9383 + 0.5);

   return cs;
}

Bool_t isSameCluster(Cluster *a, Cluster *b) {
   if (!a || !b) return false;

   Bool_t x = (a->getX() == b->getX());
   Bool_t y = (a->getY() == b->getY());
   Bool_t z = (a->getLayer() == b->getLayer());

   return x*y*z;
}

void getPValues() {
   // create list with energies vs range
   TCanvas *cCustom = new TCanvas("cCustom", "Range fit for custom range-energy list", 1200, 900);
   TCanvas *cCustomInv = new TCanvas("cCustomInv", "Range fit for custom range-energy list", 1200, 900);
   
   TCanvas *cCustomW = new TCanvas("cCustomW", "Range fit for custom range-energy list tungsten", 1200, 900);
   TCanvas *cCustomInvW = new TCanvas("cCustomInvW", "Range fit for custom range-energy list tungsten", 1200, 900);

   TCanvas *cCustomAl = new TCanvas("cCustomAl", "Range fit for AL", 1200, 900);

   Float_t energies[19] = {10, 20, 30, 40, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400};
   Float_t ranges[13] = {1.047, 4.061, 8.639, 14.655, 22.022, 45.868, 76.789, 114.074, 157.148, 205.502, 258.75, 316.419, 378.225};
   Float_t rangesW[15] = {3.027, 5.989, 9.691, 13.983, 18.899, 24.586, 30.541, 37.219, 43.945, 51.257, 59.058, 67.172, 75.522, 84.338, 93.62};

   Float_t energiesAl[9] = {25, 50, 75, 100, 125, 150, 175, 200, 225};
   Float_t rangesAl[9] = {1.66, 9.48, 21.12, 36.13, 54.14, 74.96, 98.20, 123.85, 151.35};

   cCustom->cd();
   TGraph *graph_custom = new TGraph(13, energies, ranges);
   graph_custom->SetTitle("Range fit for custom range-energy list");
   graph_custom->Draw("A*");

   cCustomInv->cd();
   TGraph *graph_inv = new TGraph(13, ranges, energies);
   graph_inv->SetTitle("Energy fit for custom energy-range list");
   graph_inv->Draw("A*");

   cCustomW->cd();
   TGraph *graph_custom_w = new TGraph(15, energies, rangesW);
   graph_custom_w->SetTitle("Range fit for W");
   graph_custom_w->Draw("A*");

   cCustomInvW->cd();
   TGraph *graph_inv_w = new TGraph(15, rangesW, energies);
   graph_inv_w->SetTitle("Energy for for W");
   graph_inv_w->Draw("A*");

   cCustomAl->cd();
   TGraph *graphAl = new TGraph(9, energiesAl, rangesAl);
   graphAl->SetTitle("Range fit for Al;Energy [MeV];Range in Aluminuim DTC [mm]");
   graphAl->Draw("A*");

   TF1 *fCustom = new TF1("fCustom", "[0] * x * (1 + ([1] - [1]* exp(-[2] * x)) + ([3] - [3] * exp(-[4]*x)))", 10, 250);
   TF1 *fCustomInv = new TF1("fCustomInv", "x * ([0] * exp ( - [1] * x) + [2] * exp( - [3] * x) + [4] * exp( - [5] * x) + [6] * exp(-[7] * x) + [8] * exp(-[9] * x))", 1, 400);
   TF1 *fCustomW = new TF1("fCustomW", "[0] * x * (1 + ([1] - [1]* exp(-[2] * x)) + ([3] - [3] * exp(-[4]*x)))", 10, 250);
   TF1 *fCustomInvW = new TF1("fCustomInvW", "x * ([0] * exp ( - 1./[1] * x) + [2] * exp( - 1./[3] * x) + [4] * exp( - 1./[5] * x) + [6] * exp(-1./[7] * x) + [8] * exp(-1./[9] * x))", 1, 400);
   TF1 *fCustomAl = new TF1("fCustomAl", "[0] * x * (1 + ([1] - [1]* exp(-[2] * x)) + ([3] - [3] * exp(-[4]*x)))", 10, 250);
   TF1 *fCustomInvAl = new TF1("fCustomInv", "x * ([0] * exp ( - [1] * x) + [2] * exp( - [3] * x) + [4] * exp( - [5] * x) + [6] * exp(-[7] * x) + [8] * exp(-[9] * x))", 1, 400);
   fCustom->SetParameters(6.94656e-2, 15.14450027, 0.001260021, 29.84400076, 0.003260031);
   fCustomInv->SetParameters(9.663872, 1/0.975, 2.50472, 1/12.4999, 0.880745, 1/57.001, 0.419001, 1/106.501, 0.92732, 1/1067.2784);
   
   fCustomW->SetParameters(6.94656e-2/9, 15.14450027, 0.001260021, 29.84400076, 0.003260031);
   fCustomInvW->SetParameters(9.663872*9, 0.975/9, 2.50472*9, 12.4999/9, 0.880745*9, 57.001/9, 0.419001*9, 106.501/9, 0.92732*9, 1067.2784/9);
   
   fCustomAl->SetParameters(6.94656e-2/2, 15.14450027, 0.001260021, 29.84400076, 0.003260031);
   fCustomInvAl->SetParameters(9.663872*2, 1/0.975*2, 2.50472*2, 1/12.4999*2, 0.880745*2, 1/57.001*2, 0.419001*2, 1/106.501*2, 0.92732*2, 1/1067.2784*2);
   
   fCustom->SetParLimits(0, 0.02, 0.1);
   fCustom->SetParLimits(1, 5, 45);
   fCustom->SetParLimits(2, 0.0005, 0.0045);
   fCustom->SetParLimits(3, 10, 90);
   fCustom->SetParLimits(4, 0.001, 0.01);
   
   fCustomAl->SetParLimits(0, 0.01, 0.1);
   fCustomAl->SetParLimits(1, 5, 45);
   fCustomAl->SetParLimits(2, 0.0005, 0.0045);
   fCustomAl->SetParLimits(3, 10, 90);
   fCustomAl->SetParLimits(4, 0.001, 0.01);
   
   fCustomInv->SetParLimits(0, 0.1, 100);
   fCustomInv->SetParLimits(1, 0.001, 5);
   fCustomInv->SetParLimits(2, 0.1, 100);
   fCustomInv->SetParLimits(3, 0.001, 5);
   fCustomInv->SetParLimits(4, 0.1, 100);
   fCustomInv->SetParLimits(5, 0.001, 5);
   fCustomInv->SetParLimits(6, 0.1, 100);
   fCustomInv->SetParLimits(7, 0.001, 5);
   fCustomInv->SetParLimits(8, 0.1, 100);
   fCustomInv->SetParLimits(9, 0.001, 5);

   fCustomInvW->SetParLimits(0, 1, 50);
   fCustomInvW->SetParLimits(1, 0.01, 0.2);
   fCustomInvW->SetParLimits(2, 5, 30);
   fCustomInvW->SetParLimits(3, 1, 8);
   fCustomInvW->SetParLimits(4, 1, 15);
   fCustomInvW->SetParLimits(5, 5, 50);
   fCustomInvW->SetParLimits(6, 0.1, 10);
   fCustomInvW->SetParLimits(7, 0.5, 30);
   fCustomInvW->SetParLimits(8, 1, 20);
   fCustomInvW->SetParLimits(9, 30, 500);

   cCustom->cd();
   fCustom->SetNpx(500);
   graph_custom->Fit("fCustom", "M, B");
   cCustomInv->cd();
   fCustomInv->SetNpx(500);
   graph_inv->Fit("fCustomInv", "M, B");
   cCustomW->cd();
   fCustomW->SetNpx(500);
   graph_custom_w->Fit("fCustomW", "M, B");
   cCustomInvW->cd();
   fCustomInvW->SetNpx(500);
   graph_inv_w->Fit("fCustomInvW", "M, B");

   cCustomAl->cd();
   graphAl->Fit("fCustomAl", "M, B");

   Float_t sqrt1 = 0, sqrt2 = 0, sqrt3 = 0, sqrt4 = 0; 
   Float_t avgRangeSum = 0, avgRangeSumW = 0;
   for (Int_t i=0; i<9; i++) {
      sqrt1 += pow(ranges[i] - fCustom->Eval(energies[i]), 2);
      sqrt2 += pow(energies[i] - fCustomInv->Eval(ranges[i]), 2);
      sqrt3 += pow(rangesW[i] - fCustomW->Eval(energies[i]), 2);
      sqrt4 += pow(energies[i] - fCustomInvW->Eval(rangesW[i]), 2);
      avgRangeSum += fabs(ranges[i] - fCustom->Eval(energies[i]));
      avgRangeSumW += fabs(rangesW[i] - fCustomW->Eval(energies[i]));
   }

   avgRangeSum /= 9;
   avgRangeSumW /= 15;

   sqrt1 = sqrt(sqrt1);
   sqrt2 = sqrt(sqrt2);
   sqrt3 = sqrt(sqrt3);
   sqrt4 = sqrt(sqrt4);

   cout << "RSQ of Gompertz-function = " << sqrt1 << endl;
   cout << "RSQ of INV fitted Gompertz-function = " << sqrt2 << endl;
   cout << "RSQ of Gompertz-function in tungsten = " << sqrt3 << endl;
   cout << "RSQ of INV fitted Gompertz-function in tungsten = " << sqrt4 << endl;
   cout << endl << "Average range error for fitted Gompertz is " << avgRangeSum << " mm.\n";
   cout << "Average range error for fitted Gompertz in tungsten is " << avgRangeSumW << " mm.\n";

   Float_t E1 = fCustomInv->Eval(fCustom->Eval(10));
   Float_t E2 = fCustomInv->Eval(fCustom->Eval(30));
   Float_t E3 = fCustomInv->Eval(fCustom->Eval(50));

   Float_t E1w = fCustomInvW->Eval(fCustomW->Eval(10));
   Float_t E2w = fCustomInvW->Eval(fCustomW->Eval(30));
   Float_t E3w = fCustomInvW->Eval(fCustomW->Eval(50));

   cout << Form("From energy of 10, 30, 50 MeV, the calc-invcalc values are %.10f, %.10f, %.10f.\n", E1, E2, E3);
   cout << Form("TUNGSTEN - From energy of 10, 30, 50 MeV, the calc-invcalc values are %.10f, %.10f, %.10f.\n", E1w, E2w, E3w);

   cout << " --- \033[1m VALUES FOR WATER \033[0m --- " << endl;
   cout << "a1      = " << fCustom->GetParameter(0) << endl;
   cout << "b1      = " << fCustom->GetParameter(1) << endl;
   cout << "g1      = " << fCustom->GetParameter(2) << endl;
   cout << "b2      = " << fCustom->GetParameter(3) << endl;
   cout << "g2      = " << fCustom->GetParameter(4) << endl;
   cout << endl;
   cout << "c1      = " << fCustomInv->GetParameter(0) << endl;
   cout << "lambda1 = " << 1/fCustomInv->GetParameter(1) << endl;
   cout << "c2      = " << fCustomInv->GetParameter(2) << endl;
   cout << "lambda2 = " << 1/fCustomInv->GetParameter(3) << endl;
   cout << "c3      = " << fCustomInv->GetParameter(4) << endl;
   cout << "lambda3 = " << 1/fCustomInv->GetParameter(5) << endl;
   cout << "c4      = " << fCustomInv->GetParameter(6) << endl;
   cout << "lambda4 = " << 1/fCustomInv->GetParameter(7) << endl;
   cout << "c5      = " << fCustomInv->GetParameter(8) << endl;
   cout << "lambda5 = " << 1/fCustomInv->GetParameter(9) << endl;

   cout << endl << endl;

   cout << " --- \033[1m VALUES FOR TUNGSTEN \033[0m --- " << endl;
   cout << "a1      = " << fCustomW->GetParameter(0) << endl;
   cout << "b1      = " << fCustomW->GetParameter(1) << endl;
   cout << "g1      = " << fCustomW->GetParameter(2) << endl;
   cout << "b2      = " << fCustomW->GetParameter(3) << endl;
   cout << "g2      = " << fCustomW->GetParameter(4) << endl;
   cout << endl;
   cout << "c1      = " << fCustomInvW->GetParameter(0) << endl;
   cout << "lambda1 = " << 1/fCustomInvW->GetParameter(1) << endl;
   cout << "c2      = " << fCustomInvW->GetParameter(2) << endl;
   cout << "lambda2 = " << 1/fCustomInvW->GetParameter(3) << endl;
   cout << "c3      = " << fCustomInvW->GetParameter(4) << endl;
   cout << "lambda3 = " << 1/fCustomInvW->GetParameter(5) << endl;
   cout << "c4      = " << fCustomInvW->GetParameter(6) << endl;
   cout << "lambda4 = " << 1/fCustomInvW->GetParameter(7) << endl;
   cout << "c5      = " << fCustomInvW->GetParameter(8) << endl;
   cout << "lambda5 = " << 1/fCustomInvW->GetParameter(9) << endl;
}


Float_t max(Float_t a, Float_t b) {
   if (a > b) return a;
   else return b;
}

Float_t min(Float_t a, Float_t b) {
   if (a < b) return a;
   else return b;
}


void drawIndividualGraphs(TCanvas *cGraph, TGraphErrors* outputGraph, Float_t fitRange, Float_t fitScale, Float_t fitError, Int_t fitIdx, Int_t eventID, Float_t *x_energy, Float_t *y_energy) {
   cGraph->cd(fitIdx+1);
   Bool_t kDrawFit = true;
   Bool_t kDrawText = true;

   outputGraph->SetMinimum(0);
   if (kHelium) outputGraph->SetMaximum(30);
   else outputGraph->SetMaximum(20);
   outputGraph->SetTitle(Form("%.0f mm Al", kAbsorberThickness));
   
   outputGraph->GetXaxis()->SetTitle("Water Equivalent Thickness [mm]");
   
/*
   Float_t WEPLFactor = getWEPLFactorFromEnergy(run_energy);

   for (int i=0; i<outputGraph->GetN(); i++) {
      outputGraph->GetX()[i] *= WEPLFactor;
   }

   Float_t weplRange = WEPLFactor * fitRange;
   Float_t weplError = WEPLFactor * fitError;
*/
   outputGraph->GetYaxis()->SetTitle("Energy loss [keV/#mum]");
   outputGraph->GetYaxis()->SetTitleOffset(1);
   outputGraph->GetXaxis()->SetTitleSize(0.05);
   outputGraph->GetYaxis()->SetTitleSize(0.05);
   outputGraph->GetXaxis()->SetLabelSize(0.05);
   outputGraph->GetYaxis()->SetLabelSize(0.05);
   outputGraph->GetXaxis()->SetTitleFont(22);
   outputGraph->GetXaxis()->SetLabelFont(22);
   outputGraph->GetYaxis()->SetTitleFont(22);
   outputGraph->GetYaxis()->SetLabelFont(22);

//   Float_t low = getUnitFromEnergy(0);
//   Float_t high = getUnitFromEnergy(run_energy) * 1.7;

   outputGraph->GetXaxis()->SetLimits(0, 450);

   outputGraph->SetMarkerColor(4);
   outputGraph->SetMarkerStyle(21);
   outputGraph->Draw("PA");

   alpha = alpha_water;
   p = p_water;
   
   TF1 *func = new TF1("fit_BP", fitfunc_DBP, 0, 500, 2);
   func->SetParameters(fitRange, fitScale);
   func->SetNpx(2500);

   func->Draw("same");

   if (x_energy) {
      Long64_t n=0;
      Long64_t j=0;
      
      for (Long64_t i=eventID*sizeOfEventID; i<(eventID+1)*sizeOfEventID; i++) {
         if (y_energy[i] == 0) {
            n = j;
            break;
         }
         y_energy[i] *=  2;
         
         j++;
      }
      
      if (n>1) {
         TGraph *energyLoss = new TGraph(n, (x_energy + eventID*sizeOfEventID), (y_energy + eventID*sizeOfEventID));
         energyLoss->SetLineColor(kBlack); energyLoss->SetLineWidth(2);
         energyLoss->Draw("same");
      }
      
      Float_t lastX = x_energy[eventID*sizeOfEventID + n-1];
      TLine *l = new TLine(lastX, 0, lastX, 600);
      l->Draw();
      
      Float_t realEnergy = getEnergyFromWEPL(x_energy[eventID*sizeOfEventID + n-1]);
      TLatex *text2 = new TLatex(15, 28, Form("'Real' energy: %.1f #pm MeV", realEnergy));
      text2->SetTextSize(0.045);
      text2->SetTextFont(22);
      text2->Draw();
   }

   if (kDrawText) {
//      Float_t energy_error = getEnergyStragglingFromWEPLStraggling(fitRange, fitError); // remove this, assumes WEPL
//      Float_t energy_fit = getEnergyFromWEPL(fitRange); // remove this, assumes WEPL
//      TLatex *text = new TLatex(15, 18, Form("Estimated E_{0}: %.1f #pm %.1f mm", fitRange, fitError)); // Take back this
      TLatex *text = new TLatex(15, 25, Form("Fitted range: %.1f #pm %.1f mm WET", fitRange, fitError));
      text->SetTextSize(0.05);
      text->SetTextFont(22);
      text->Draw();
   }
   
   outputGraph->GetXaxis()->SetRangeUser(0, 450);

   cGraph->Update();
}

TF1 *  doSimpleGaussianFit (TH1F *h, Float_t *means, Float_t *sigmas, Int_t idx_txt) {
   Float_t estimated_energy = 0, estimated_range = 0;
   Float_t estimated_energy_error = 0;
   Float_t sum_constant = 0, sumSigma = 0;

   Float_t nominalMean = getUnitFromEnergy(run_energy);
   Float_t nominalSigma = getUnitStragglingFromEnergy(run_energy, 0);

   // nominalMean = getWETFromDegrader(run_degraderThickness);
   // WTF THIS IS A BAD BIAS TO HAVE

   nominalMean = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
   nominalSigma = 4;

   TAxis *axis = h->GetXaxis();
   TF1 *gauss = new TF1("gauss", "gaus");
   gauss->SetParameter(1, nominalMean);
   gauss->SetParameter(2, nominalSigma);
   if (idx_txt < 0) {
      h->Fit("gauss", "B,Q,W,M");
   }
   else {
      h->Fit("gauss", "B,Q,W,M,N");
   }

   Float_t mu, sigma;
   mu = gauss->GetParameter(1);
   sigma = fabs(gauss->GetParameter(2));

   // nominalSigma is too low at very low energies
   // fitted sigma may be too low at high energies, where fit can be ~= 1 layer
   Float_t sigmaToUse = fmax(sigma, nominalSigma);

   Float_t sigmaFrom = fmax(mu - 4 * sigmaToUse, 0);
   Float_t sigmaTo = fmin(mu + 4 * sigmaToUse, nominalMean * 1.4 + 10);

   Int_t binSigmaFrom = axis->FindBin(sigmaFrom);
   Int_t binSigmaTo = axis->FindBin(sigmaTo);
   
   if (binSigmaTo > h->GetNbinsX()) binSigmaTo = h->GetNbinsX();
   if (binSigmaTo < binSigmaFrom) binSigmaTo = h->GetNbinsX();

   printf("doSimpleGaussianFit: Searching from %.2f (bin %d)  to %.2f (bin %d).\n", sigmaFrom, binSigmaFrom, sigmaTo, binSigmaTo);

   Float_t squareMeanDifference = 0;
   Float_t empiricalMean = 0;
   Int_t N = 0;

   for (Int_t i=binSigmaFrom; i<=binSigmaTo; i++) {
      empiricalMean += h->GetBinContent(i) * axis->GetBinCenter(i);
      N += h->GetBinContent(i);
   } 
  
   empiricalMean /= N;

   for (Int_t i=binSigmaFrom; i<=binSigmaTo; i++) {
      squareMeanDifference += h->GetBinContent(i) * pow(axis->GetBinCenter(i) - empiricalMean, 2);
   }

   Float_t empiricalSigma = sqrt(abs(squareMeanDifference / N));
   cout << "The empirical estimated range is " << empiricalMean << " mm. (which is " << getEnergyFromUnit(empiricalMean) << " MeV).\n";
   cout << "The empirical standard deviation is " << empiricalSigma << " mm.\n";
   

   // CALCULATE FWTM
   Int_t bin1 = h->FindFirstBinAbove(h->GetMaximum()/10);
   Int_t bin2 = h->FindLastBinAbove(h->GetMaximum()/10);
   Float_t fwtm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);

   bin1 = h->FindFirstBinAbove(h->GetMaximum()/2);
   bin2 = h->FindLastBinAbove(h->GetMaximum()/2);
   Float_t fwhm = h->GetBinCenter(bin2) - h->GetBinCenter(bin1);

   Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;
 
   TString filename = "OutputFiles/result_makebraggpeakfit";
   if (idx_txt > 0) {
      filename.Append(Form("_idx%d.csv", idx_txt));
   }
   
   else {
      filename.Append(".csv");
   }

   printf("Calculated FWTM is %.2f mm. Actual integration yields %.2f mm.\n", empiricalSigma*4.2913, fwtm);

   ofstream file2(filename, ofstream::out | ofstream::app);

   // absorber thickness; energy; nominal range; estimated range; range sigma
   if (run_degraderThickness == 0) {
      file2 << readoutAbsorber << " " << run_energy << " " << getUnitFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
   }
   else {
      // file2 << readoutAbsorber << " " << run_degraderThickness << " " << getUnitFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
      file2 << readoutAbsorber << " " << run_degraderThickness << " " << getWETFromDegrader(run_degraderThickness) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << fwhm << " " << fwtm << endl;
   }

   file2.close();
   
   means[9] = empiricalMean;
   sigmas[9] = empiricalSigma;

   return gauss;

}

Float_t doNGaussianFit ( TH1F *h, Float_t *means, Float_t *sigmas) {
   TF1 *gauss;
   cout << "Energy " << run_energy << endl;
   cout << "Layer;\t Constant;\t Mean;\t\t Energy;\t Sigma;\t\t Fits in layer;\t Chi2/n\n";
   
   Float_t constant, mean, lEnergy, sigma;
   Float_t meanValueFromSumming = 0;
   
   Float_t array_mean[3] = {};
   Int_t array_layer[3] = {};
   Float_t array_constant[3] = {};
   Float_t array_sigma[3] = {};
   Float_t array_sum[3] = {};
   Float_t array_energy[3] = {};
   Float_t array_f[3] = {};
   Float_t array_chi2n[3] = {};
   
   TAxis *axis = h->GetXaxis();
   Float_t fullIntegral = h->Integral();
   Int_t binSearchFrom = axis->FindBin(getWEPLFromEnergy(run_energy * 0.9));
   Int_t binSearchTo   = axis->FindBin(getWEPLFromEnergy(run_energy * 1.1));
   Float_t partialIntegral = h->Integral(binSearchFrom, binSearchTo);

   for (Int_t i=binSearchFrom; i<binSearchTo; i++) {
      meanValueFromSumming += axis->GetBinCenter(i) * h->GetBinContent(i) / partialIntegral;
   }

   Float_t error = abs(100 * (meanValueFromSumming - getWEPLFromEnergy(run_energy)) / getWEPLFromEnergy(run_energy));
   cout << "\033[1mUsing the sum method, the weighted mean value of the distribution is " << meanValueFromSumming << " mm. (" << error << ")\033[0m\n";
   
   Float_t maxBinHeight = h->GetMaximum();

   Bool_t isLastLayer, wasLastLayer = false;;
   

   Int_t j=0;
   for (Int_t i=0; i<nLayers; i++) {
      if (getWEPLFromTL(getLayerPositionmm(i)) > getUnitFromEnergy(run_energy*1.2)) continue;
      if (getWEPLFromTL(getLayerPositionmm(i)) < getWEPLFromEnergy(run_energy * 0.8)) continue;

      Float_t searchFrom = getWEPLFromTL(getLayerPositionmm(i))+10;
      Float_t searchTo = getWEPLFromTL(getLayerPositionmm(i+1))+10;
      
      Int_t bmin = axis->FindBin(searchFrom);
      Int_t bmax = axis->FindBin(searchTo);
      
      Float_t integral = h->Integral(bmin,bmax);
      Float_t ratio = integral / fullIntegral;
      
      isLastLayer = ((getWEPLFromTL(getLayerPositionmm(i)) > getUnitFromEnergy(run_energy-10)) && !wasLastLayer && ratio > 0.01) ;
   
      if (i<=3) continue;
      if (ratio < 0.1 && !isLastLayer) continue;
      
      gauss = new TF1(Form("Gaus_%d", i), "gaus(0)", searchFrom, searchTo);
   
      gauss->SetLineWidth(2);

      sigma = 3;
      if (isLastLayer && ratio < 0.05) sigma = 0.2;
      
      gauss->SetParameters(10, (searchFrom+searchTo)/2, sigma);
      gauss->SetParLimits(0, 0, maxBinHeight);
      gauss->SetParLimits(1, searchFrom+5, searchTo);
      gauss->SetParLimits(2, 2, 9);
      
      h->Fit(gauss, "M, B, WW, Q, 0", "", searchFrom, searchTo);
      
      sigma = gauss->GetParameter(2);
      constant = gauss->GetParameter(0);
      mean = gauss->GetParameter(1);
      lEnergy = getEnergyFromUnit(mean);
      
      Float_t integralSigma = h->Integral(axis->FindBin(mean - sigma), axis->FindBin(mean + sigma));
      Float_t chi2 = gauss->GetChisquare();
      Float_t chi2n = chi2 / integral;
      Float_t chi2nSigma = chi2 / integralSigma;
      
      cout << Form("Searching from %.1f to %.1f, with midpoint at %.1f. Found best fit @ %.1f with chi2 = %.2f and chi2/n = %.2f and chi2/nSIGMA = %.2f, ratio = %.2f.\n", searchFrom, searchTo,(searchTo+searchFrom)/2 , mean, chi2, chi2n, chi2nSigma, ratio);

   /*
      if (run_energy > 182 && run_energy < 193) {
         if (chi2n > 200 || sigma < 2.5 || constant < 3) {
            delete gauss;
            continue;
         }
      }
   */

      if (ratio > 0.1) {
         gauss->SetLineColor(kRed);
         gauss->SetLineWidth(3);
//         gauss->Draw("same");
         cout << Form("%d;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f;\t %8.5f;\t %.2f\n", i, constant, mean, lEnergy, sigma, ratio, chi2n);
         
         if (j<3) {
            array_mean[j] = mean;
            array_layer[j] = i;
            array_constant[j] = constant;
            array_energy[j] = lEnergy;
            array_sigma[j] = sigma;
            array_f[j] = ratio;
            array_sum[j] = integralSigma;
         }
         
         sigmas[j] = sigma;
         means[j++] = mean;
      }
      wasLastLayer = isLastLayer;
   }

   Float_t estimated_energy = 0, estimated_range = 0;
   Float_t estimated_energy_error = 0;
   Float_t sum_constant = 0, sumSigma = 0;
   for (Int_t i=0 ; i<3; i++) {
      estimated_range += array_sum[i] * array_mean[i];
      sum_constant += array_sum[i];
      sumSigma += pow(array_sigma[i], 2);
   }

   sumSigma = sqrt(sumSigma);

   estimated_range /= sum_constant;

   estimated_energy = getEnergyFromUnit(estimated_range);
   estimated_energy_error = getEnergyFromWEPL(estimated_range + sumSigma/2) - getEnergyFromWEPL(estimated_range - sumSigma/2);
   cout << "ESTIMATED ENERGY FROM RUN IS " << estimated_energy << " +- " << estimated_energy_error << endl;
   cout << "Estimated range = " << estimated_range << " +- " << sumSigma << endl;
  
   Float_t lastMean = array_mean[1];
   Float_t lastSigma = array_sigma[1];
   if (lastMean == 0 || axis->FindBin(lastMean) > h->GetNbinsX()) {
      lastMean = array_mean[0];
      lastSigma = array_sigma[0];
   }

   printf("lastMean = %.2f. lastSigma = %.2f\n", lastMean, lastSigma);
   printf("sigmaFrom = %.2f\n", lastMean - 9 * lastSigma);

   Int_t binSigmaFrom = axis->FindBin(lastMean - 9*lastSigma);

   if (lastMean == 0) {
      binSigmaFrom = axis->FindBin(getWEPLFromEnergy(run_energy)*0.9);
      cout << "Setting binSigmaFrom to " << binSigmaFrom;
   }

   Float_t squareMeanDifference = 0;
   Float_t empiricalMean = 0;
   Int_t N = 0;

   printf("Looping from %d to %d.\n", binSigmaFrom, h->GetNbinsX());

   for (Int_t i=binSigmaFrom; i<=h->GetNbinsX(); i++) {
      empiricalMean += h->GetBinContent(i) * axis->GetBinCenter(i);
      N += h->GetBinContent(i);
   } 
  
   empiricalMean /= N;

   printf("empiricalMean = %.2f\n", empiricalMean);

   for (Int_t i=binSigmaFrom; i<=h->GetNbinsX(); i++) {
      squareMeanDifference += h->GetBinContent(i) * pow(axis->GetBinCenter(i) - empiricalMean, 2);
   }
   printf("squareMeanDifference = %.2f\n", squareMeanDifference);

   Float_t empiricalSigma = sqrt(abs(squareMeanDifference / N));
   cout << "The empirical estimated range is " << empiricalMean << " mm. (which is " << getEnergyFromWEPL(empiricalMean) << " MeV).\n";
   cout << "The empirical standard deviation is " << empiricalSigma << " mm.\n";
   
   ofstream file("OutputFiles/output_gauss.csv", ofstream::out | ofstream::app);
   // run_energy; layer[i], constant[i], mpv[i], energy[i], sigma[i], ratio[i];  i 1->3
   
   file << run_energy << "; ";
   for (Int_t j=0; j<3; j++) {
      file << Form("%d; %.5f; %.5f; %.5f; %.5f; %.5f; %.5f",
                array_layer[j], array_constant[j], array_mean[j],
                array_energy[j], array_sigma[j], array_f[j], array_chi2n[j]);
   }

   file << endl;
   file.close();
//    printf("Calculated FWTM is %.2f mm. Actual integration yields %.2f mm.\n", empiricalSigma*4.2);
   
   Float_t last_range = array_mean[0];
   if (array_constant[1] > 0.1) last_range = array_mean[1];
   if (array_constant[2] > 0.1) last_range = array_mean[2];

   Float_t nominalSigma = getWEPLStragglingFromEnergy(run_energy, 0);

   Float_t readoutAbsorber = (roundf(kAbsorberThickness) == kAbsorberThickness) ? kAbsorberThickness : kAbsorberThickness*10;

   ofstream file2("OutputFiles/result_makebraggpeakfit.csv", ofstream::out | ofstream::app);
   // absorber thickness; energy; nominal range; estimated range; range sigma
   if (run_degraderThickness == 0) {
      file2 << readoutAbsorber << " " << run_energy << " " << getWEPLFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
   }
   else {
      file2 << readoutAbsorber << " " << run_degraderThickness << " " << getWEPLFromEnergy(run_energy) << " " << empiricalMean << " " << nominalSigma << " " << empiricalSigma << " " << endl;
   }

   file2.close();

   means[9] = empiricalMean;
   sigmas[9] = empiricalSigma;

   return empiricalMean;
}

Bool_t getCutTrackLength(Float_t energy, Track *track) {
   Int_t minTL = getMinimumTrackLength(energy);
   Float_t TL = track->getTrackLengthmm();

   Bool_t cutTL = (TL > minTL) ? true : false;

   return cutTL;
}

Bool_t getCutWEPL(Track *track) {
   Float_t minTLWEPL = 150;
   Float_t WEPL = track->getWEPL();

   Bool_t cutTL = (WEPL > minTLWEPL) ? true : false;

   return cutTL;
}

Bool_t getCutChipNumber(Track *track) {
   Int_t x0 = track->getX(0);
   Int_t y0 = track->getY(0);
   Int_t chip = (x0 >= nx/2) + 2 * (y0 < ny/2);
   Bool_t cutChipNumber = (chip<2) ? true : false;

   return cutChipNumber;
}

Bool_t getCutBraggPeakInTrack(Track *track) {
   Float_t braggPeakRatio = 2.5;

   Int_t lastBin = track->GetEntriesFast() - 1;
   if (lastBin < 4) return false;

   Float_t lowStd = track->getStdSizeToIdx(lastBin-1);
   Int_t   nEmptyBins = track->getNMissingLayers();


   if (lowStd > 4) return false;
   if (nEmptyBins > 1) return false;

   Float_t rMean = track->getMeanSizeToIdx(lastBin - 1);
   Float_t rrMean = track->getMeanSizeToIdx(lastBin - 2);

   Float_t r = track->getSize(lastBin) / rMean;
   Float_t rr = track->getSize(lastBin-1) / rrMean;

   if (r > braggPeakRatio) return true;
   else {
      if (rr > braggPeakRatio) return true;
      else return false;
   }
}

Float_t getAngleAtSpot(Float_t spotPosInMM) {
   Float_t dIsocenter = 6600;
   return atan2(spotPosInMM, dIsocenter);
}

Float_t getAngleAtSpot(Float_t spotPosInMM_x, Float_t spotPosInMM_y) {
   // Returns in Rad
   Float_t dIsocenter = 6600;
   Float_t theta_x = atan2(spotPosInMM_x, dIsocenter);
   Float_t theta_y = atan2(spotPosInMM_y, dIsocenter);

   return sqrt(pow(theta_x, 2) + pow(theta_y, 2));
}

#endif
