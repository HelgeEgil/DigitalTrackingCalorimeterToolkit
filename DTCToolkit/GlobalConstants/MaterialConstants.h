#ifndef MaterialConstants_h
#define MaterialConstants_h

#include <list>
#include <TTree.h>
#include <TTreeIndex.h>

#include <TObject.h>
#include <TSpline.h>
#include <TBranch.h>

namespace DTC {
class Cluster;
}

Int_t nLayers;

// range = alpha * energy ^ p
// range straggling is proportionate to alpha_prime
Double_t p, alpha, alpha_prime;
Double_t p_water, alpha_water, alpha_prime_water;
Double_t p_aluminum, alpha_aluminum, alpha_prime_aluminum;
Double_t p_pure_aluminum, alpha_pure_aluminum, alpha_prime_pure_aluminum;
Double_t p_tungsten, alpha_tungsten, alpha_prime_tungsten;
Double_t straggling_a, straggling_b; // straggling = straggling_a + straggling_b * range

Float_t p_material_high, p_material_low, alpha_material_high, alpha_material_low;

TSpline3 *splineDTC;
TSpline3 *splineDTCInv;
TSpline3 *splineWater;
TSpline3 *splineWaterInv;
TSpline3 *splineW;
TSpline3 *splineWInv;
TSpline3 *splineMaterial;
TSpline3 *splineMaterialInv;
TSpline3 *splinePureAl;
TSpline3 *splinePureAlInv;
TSpline3 *splineWET;
TSpline3 *splineWEPL;

Double_t X0, X0_tungsten, X0_aluminum, X0_pmma, X0_firstlayer;
Double_t proton_mass;

Float_t mcs_radius_per_layer[100]; // change if max(nLayers) > 100
Float_t mcs_radius_per_layer_empirical[100];

Double_t firstUpperLayerZ, firstLowerLayerZ;

TTree        * CDB_treeCluster;
//list<Int_t>  * CDB_hit_array = 0;
uint8_t        CDB_hit_array[10];
TBranch      * CDB_b_hit_array;
TTreeIndex   * CDB_index;
Int_t          CDB_clusterSize;
Float_t        CDB_x_mean, CDB_y_mean;
Int_t          CDB_sortIndex[50];
TFile        * CDB_fCluster; 

void     createSplines();
Double_t getLayerPositionmm(Int_t i);
Float_t  getSigmaEnergy(Int_t energy);
Bool_t   isChipLowResistivity(Int_t chipIdx);
Float_t  getChipCalibrationFactor(Int_t chipIdx);
Bool_t   isBadData(DTC::Cluster * estimatedPosition);

#endif
