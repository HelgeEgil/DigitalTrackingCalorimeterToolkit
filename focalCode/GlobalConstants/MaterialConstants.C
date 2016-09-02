#include <cmath>
#include <vector>
#include <iostream>

#include <TObject.h>

#include "MaterialConstants.h"
#include "Classes/Hit/Hit.h"
#include "Classes/Cluster/Cluster.h"
#include "Constants.h"

using namespace std;
	
// Full revision of range-energy and depth-dose relationship as of 2016-06
// Using Ulmer Rad Phys and Chem 76 (2007) 1089-1107
// In short:
// R = a1 E0 * (1 + sum_(k=1)^2 (bk - bk * exp(-gk * E0)))
// E = R * sum_(i=1)^5 (ck * exp(-lambdak * R))
//		 ( and  E(z) = (R-z) * sum_(i=1)^5 (ck * exp(-lambdak * (R-z))) )
// -dE/dz = E(z) / (R-z) - sum_(k=1)^5 (lambdak * ck * (R-z) * exp(-lambdak * (R-z)))
// The parameters a1, bk(1-2), gk(1-2), ck(1-5), lambdak(1-5) are found through range-energy fitting on different materials in GATE simulations, and defined below
// All calculations are performed in Classes/Track/conversionFunctions.C

void MaterialConstants() {

	a1_water = 0.0727504;
	b1_water = 44.7718;
	g1_water = 0.00108048;
	b2_water = 13.4883;
	g2_water = 0.00457689;

	c1_water = 10.7425;
	l1_water = 0.500803;
	c2_water = 1.7655;
	l2_water = 0.0734755;
	c3_water = 0.95865;
	l3_water = 0.0236093;
	c4_water = 0.723177;
	l4_water = 0.00662777;
	c5_water = 0.733621;
	l5_water = 0.000523239;

	a1_tungsten = 0.0128998;
	b1_tungsten = 39.5341;
	g1_tungsten = 0.00080173;
	b2_tungsten = 6.67826;
	g2_tungsten = 0.00687635;

	c1_tungsten = 1.01309;
	l1_tungsten = 22.4927;
	c2_tungsten = 12.6053;
	l2_tungsten = 0.472777;
	c3_tungsten = 4.31023;
	l3_tungsten = 0.0361097;
	c4_tungsten = 6.1832;
	l4_tungsten = 0.122736;
	c5_tungsten = 5.42993;
	l5_tungsten = 0.00293302;

	// USING VALUES FOR WATER FOR ALUMINIUM AND PMMA FIXME

	a1_aluminium = 0.0727504;
	b1_aluminium = 44.7718;
	g1_aluminium = 0.00108048;
	b2_aluminium = 13.4883;
	g2_aluminium = 0.00457689;
	c1_aluminium = 10.7425;
	l1_aluminium = 0.500803;
	c2_aluminium = 1.7655;
	l2_aluminium = 0.0734755;
	c3_aluminium = 0.95865;
	l3_aluminium = 0.0236093;
	c4_aluminium = 0.723177;
	l4_aluminium = 0.0662777;
	c5_aluminium = 0.733621;
	l5_aluminium = 0.000523239;

	a1_pmma = 0.0727504;
	b1_pmma = 44.7718;
	g1_pmma = 0.00108048;
	b2_pmma = 13.4883;
	g2_pmma = 0.00457689;
	c1_pmma = 10.7425;
	l1_pmma = 0.500803;
	c2_pmma = 1.7655;
	l2_pmma = 0.0734755;
	c3_pmma = 0.95865;
	l3_pmma = 0.0236093;
	c4_pmma = 0.723177;
	l4_pmma = 0.0662777;
	c5_pmma = 0.733621;
	l5_pmma = 0.000523239;

	firstUpperLayerZ = 0.3;
	firstLowerLayerZ = 0.6;

	// these values are still used for calculating the range and energy straggling!
	// The Bortfeld approximation will introduce small errors, but they are without importance (Ulmer 2007)

	p_water = 1.7547;
	alpha_water = 0.02387;
	alpha_prime_water = 0.0087;
	
	p_aluminum = 1.677;
	alpha_aluminum = 0.01738;
	alpha_prime_aluminum = 0.017102;
	
	p_tungsten = 1.6677;
	alpha_tungsten = 0.004461;
	alpha_prime_tungsten = 0.1086784;

	proton_mass = 938.27;
	X0_tungsten = 4.2; // 3.857;
	X0_firstlayer = 33.36;
	X0_aluminum = 5.88;
	X0_pmma = 16.52;

	for (Int_t i=0; i<100; i++) {
		mcs_radius_per_layer[i] = 0;
	}

	if (kMaterial == kTungsten) {
		nLayers = 41;
		p = p_tungsten;
		alpha = alpha_tungsten;
		alpha_prime = alpha_prime_tungsten;
		X0 = X0_tungsten;
		
		a1_material = a1_tungsten;
		b1_material = b1_tungsten;
		g1_material = g1_tungsten;
		b2_material = b2_tungsten;
		g2_material = g2_tungsten;
		c1_material = c1_tungsten;
		l1_material = l1_tungsten;
		c2_material = c2_tungsten;
		l2_material = l2_tungsten;
		c3_material = c3_tungsten;
		l3_material = l3_tungsten;
		c4_material = c4_tungsten;
		l4_material = l4_tungsten;
		c5_material = c5_tungsten;
		l5_material = l5_tungsten;
	}

	else if (kMaterial == kAluminum) {
		nLayers = 41;
		p = p_aluminum;
		alpha = alpha_aluminum;
		alpha_prime = alpha_prime_aluminum;
		X0 = X0_aluminum;

		a1_material = a1_aluminium;
		b1_material = b1_aluminium;
		g1_material = g1_aluminium;
		b2_material = b2_aluminium;
		g2_material = g2_aluminium;
		c1_material = c1_aluminium;
		l1_material = l1_aluminium;
		c2_material = c2_aluminium;
		l2_material = l2_aluminium;
		c3_material = c3_aluminium;
		l3_material = l3_aluminium;
		c4_material = c4_aluminium;
		l4_material = l4_aluminium;
		c5_material = c5_aluminium;
		l5_material = l5_aluminium;

	}

	else if (kMaterial == kPMMA) {
		nLayers = 65;
		alpha_prime = alpha_prime_water; // MUST VERIFY
		p = p_water; // MUST VERIFY
		alpha = alpha_water; // MUST VERIFY
		X0 = X0_pmma;
		
		a1_material = a1_pmma;
		b1_material = b1_pmma;
		g1_material = g1_pmma;
		b2_material = b2_pmma;
		g2_material = g2_pmma;
		c1_material = c1_pmma;
		l1_material = l1_pmma;
		c2_material = c2_pmma;
		l2_material = l2_pmma;
		c3_material = c3_pmma;
		l3_material = l3_pmma;
		c4_material = c4_pmma;
		l4_material = l4_pmma;
		c5_material = c5_pmma;
		l5_material = l5_pmma;
	}
}

Float_t getEnergyLossFromScintillators(Float_t energy, Int_t nScintillators) {
	Float_t energyLoss = 0;
	
 	if 	  (nScintillators == 0) energyLoss = 0;
 	else if (nScintillators == 1) energyLoss = 214.88444 * pow(energy, -0.7264956);
 	else if (nScintillators == 2) energyLoss = 341.79255 * pow(energy, -0.7357198);
 	else if (nScintillators == 3) energyLoss = 482.63531 * pow(energy, -0.7453624);
	
	else {
		cout << "ASSERTION ERROR ON NUMBER OF SCINTILLATORS!\n";
		energyLoss = 0;
	}
	
	return energyLoss;
}

Float_t getEnergyLossErrorFromScintillators(Int_t nScintillators) {
	Float_t energyLossError = 0;

	if 		(nScintillators == 0) energyLossError = 0;
	else if (nScintillators == 1) energyLossError = 0.30;
	else if (nScintillators == 2) energyLossError = 0.38;
	else if (nScintillators == 3) energyLossError = 0.44;
	
	return energyLossError;
}

Float_t getEnergyLossFromAluminumAbsorber(Float_t energy) {
	return 62.00441 * pow(energy, -0.7158403);
}

Float_t getEnergyLossErrorFromAluminumAbsorber() {
	return 0.135;
}

Double_t getLayerPositionmm(Int_t i) {
	Double_t z = 0;

	if (i>0) {
		z = ( firstUpperLayerZ + firstLowerLayerZ ) / 2 + dz * i;
	}

	return z;
}

Float_t getSigmaEnergy(Int_t energy) { 
	Float_t sigma_energy = 0;
	
	if 	  (energy == 188) sigma_energy = 0.3; // approx
	else if (energy == 180) sigma_energy = 2.3;
	else if (energy == 170) sigma_energy = 4.5; // was 4.5
	else if (energy == 160) sigma_energy = 8.3;
	else if (energy == 150) sigma_energy = 20;
	
	return sigma_energy;
}

Bool_t isChipLowResistivity(Int_t chipIdx) {
	Bool_t isHigh = false;

	if (chipIdx == 2  || chipIdx == 3  ||
		 chipIdx == 8  || chipIdx == 9  ||
		 chipIdx == 16 || chipIdx == 17 ||
		 chipIdx == 18 || chipIdx == 21 ||
		 chipIdx == 23 || chipIdx == 24 ||
		 chipIdx == 27) {
		isHigh = true;
	}

	return isHigh;
}
