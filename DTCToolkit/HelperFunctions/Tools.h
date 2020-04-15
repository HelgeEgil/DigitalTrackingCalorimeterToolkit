#ifndef FOCALCODE_TOOLS_H_
#define FOCALCODE_TOOLS_H_

#include <vector>
#include <TObject.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "../Classes/Hit/Hit.h"
#include "../Classes/Cluster/Cluster.h"

class TGraph;
class TH1F;
class TH1;

using namespace std;
using namespace DTC;

void        BinLogY(TH1 *h);
Bool_t      isItemInVector(Int_t i, vector<Int_t> *v);
Float_t     diffXY(Cluster *p, Hit *h);
Float_t     diffXY(Cluster *p1, Cluster *p2);
Float_t     diffmmXY(Cluster *p1, Cluster *p2);
Float_t     diffmmXY(Hit *h1, Hit *h2);
Float_t     diffmmXYZ(Cluster *p1, Cluster *p2);
Bool_t      existsEnergyFile(Int_t energy);
Double_t    fitfunc_DBP(Double_t *v, Double_t *par);
Double_t    double_landau(Double_t *v, Double_t *par);
char      * getMaterialChar();
char      * getMaterialAbbr();
char      * getDataTypeChar(Int_t dataType);
Int_t       getMinimumTrackLength(Float_t energy);
Int_t       getFWxMInRange(TH1F* h, Float_t first, Float_t last, Int_t div);
Float_t     quadratureAdd(Float_t a, Float_t b);
Float_t     getEnergyFromXY(Float_t *x_energy, Float_t *y_energy, Int_t eventID);
void        convertXYToWEPL(Float_t *x_energy, Float_t *y_energy, Int_t eventID);
Double_t    correctForEnergyParameterisation(Float_t energy);
Float_t     getAverageEnergyLoss(Float_t energy);
Float_t     calculateMCS(Float_t energy, Float_t depth);
Float_t     findMCSPixelRadiusForLayer(Float_t layer, Float_t E0);
void        fillMCSRadiusList(Float_t angleFactor = 1);
void        multiplyRadiusFirstLayers(Float_t factor);
Float_t     getMCSAngleForLayer(Int_t layer);
Float_t     findMCSAtLayerRad(Int_t layer, Float_t E0);
Float_t     getSearchRadiusForLayer(Int_t layer);
Float_t     getEmpiricalMCSAngle(Int_t layer);
Float_t     getDotProductAngle(Cluster *a, Cluster *b, Cluster *c);
Hit       * sumHits(Hit * a, Hit * b);
Cluster   * getTrackExtrapolationCluster(Cluster *p1, Cluster *p2);
Cluster   * getTrackExtrapolationToLayer(Track * track, Int_t layer);
Cluster   * getTrackExtrapolationFromTo(Track * track, Int_t fromLayer, Int_t toLayer);
Cluster   * getRetrogradeTrackExtrapolationToLayer(Track * track, Int_t layer);
Bool_t      isPointOutOfBounds(Cluster * point, Float_t padding = 0);
Float_t     getEdepFromCS(Int_t cs);
Int_t       getCSFromEdep(Float_t edep);
Bool_t      isSameCluster(Cluster *a, Cluster *b);
void        getPValues();
Float_t     max(Float_t a, Float_t b);
Float_t     min(Float_t a, Float_t b);
void			drawIndividualGraphs(TCanvas *cGraph, TGraphErrors* outputGraph, Float_t fitEnergy, Float_t fitScale, Float_t fitError, Int_t fitIdx, Int_t skipIdx = 0, Float_t *x = 0, Float_t *y = 0);
Float_t		doNGaussianFit( TH1F *h, Float_t *means, Float_t *sigmas);
TF1 *       doSimpleGaussianFit(TH1F *h, Float_t *means, Float_t *sigmas, Int_t idx_txt);
Bool_t		getCutTrackLength(Float_t energy, Track *track);
Bool_t		getCutWEPL(Track *track);
Bool_t		getCutChipNumber(Track *track);
Bool_t		getCutBraggPeakInTrack(Track *track);
Float_t     getAngleAtSpot(Float_t spotPosInMM);
Float_t     getAngleAtSpot(Float_t spotPosInMM_x, Float_t spotPosInMM_y);
#endif /* FOCALCODE_TOOLS_H_ */
