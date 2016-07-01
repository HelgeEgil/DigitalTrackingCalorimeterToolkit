#ifndef FOCALCODE_TOOLS_H_
#define FOCALCODE_TOOLS_H_

#include <vector>
#include <TObject.h>

#include "../Classes/Hit/Hit.h"
#include "../Classes/Cluster/Cluster.h"

class TGraph;
class TH1F;

using namespace std;

Bool_t isItemInVector(Int_t i, vector<Int_t> *v);
Float_t diffXY(Cluster *p, Hit *h);
Float_t diffXY(Cluster *p1, Cluster *p2);
Float_t diffmmXY(Cluster *p1, Cluster *p2);
Float_t diffmmXY(Hit *h1, Hit *h2);
Float_t diffmmXYZ(Cluster *p1, Cluster *p2);
Bool_t existsEnergyFile(Int_t energy);
Double_t fitfunc_DBP(Double_t *v, Double_t *par);
Double_t double_landau(Double_t *v, Double_t *par);
char * getMaterialChar();
char * getDataTypeChar(Int_t dataType);
Int_t getMinimumTrackLength(Float_t energy);
Int_t getFWxMInRange(TH1F* h, Float_t first, Float_t last, Int_t div);
Float_t quadratureAdd(Float_t a, Float_t b);
Float_t getEnergyFromXY(Float_t *x_energy, Float_t *y_energy, Int_t eventID);
void convertXYToWEPL(Float_t *x_energy, Float_t *y_energy, Int_t eventID);
Double_t correctForEnergyParameterisation(Float_t energy);
Float_t getAverageEnergyLoss(Float_t energy);

Float_t calculateMCS(Float_t energy, Float_t depth);
Float_t findMCSPixelRadiusForLayer(Float_t layer, Float_t E0);
void fillMCSRadiusList(Float_t angleFactor);
void multiplyRadiusFirstLayers(Float_t factor);
Float_t getMCSAngleForLayer(Int_t layer);
Float_t findMCSAtLayerRad(Int_t layer, Float_t E0);
Float_t getSearchRadiusForLayer(Int_t layer);
Hit * sumHits(Hit * a, Hit * b);
Cluster * getTrackPropagationToLayer(Track * track, Int_t layer);
Cluster * getTrackPropagationFromTo(Track * track, Int_t fromLayer, Int_t toLayer);
Cluster * getRetrogradeTrackPropagationToLayer(Track * track, Int_t layer);
Bool_t isPointOutOfBounds(Cluster * point, Float_t padding = 0);
Bool_t isSameCluster(Cluster *a, Cluster *b);
Float_t max(Float_t a, Float_t b);

#endif /* FOCALCODE_TOOLS_H_ */
