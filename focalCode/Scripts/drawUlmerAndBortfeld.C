#include <vector>
#include <algorithm>

#include <TObject.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>

#include "Classes/Cluster/Cluster.h"
#include "Classes/Hit/Hit.h"
#include "GlobalConstants/Constants.h"

void drawUlmerAndBortfeld() {
	Float_t x[8000];
	Float_t dEdx[8000] = {0};
	Float_t dEdx2[8000] = {0};
	Float_t dEdx3[8000] = {0};
	Float_t E[8000] = {0};

	float aPhi1[8000] = {0};
	float aPhi2[8000] = {0};
	float aPhi3[8000] = {0};
	float aPhi4[8000] = {0};

	float factor = 1;

	float p = 1.6677;
	float a = 0.04461;

//	 float c1_water = 8.55221 / factor;
//	 float l1_water = 0.680604 * factor;
//	 float c2_water = 2.69391 / factor;
//	 float l2_water = 0.139472 * factor;
//	 float c3_water = 1.32594 / factor;
//	 float l3_water = 0.0383079 * factor;
//	 float c4_water = 0.894249 / factor;
//	 float l4_water = 0.0099829 * factor;
//	 float c5_water = 0.868646 / factor;
//	 float l5_water = 0.000805688 * factor;
/*
	 float c1_water = 76.9241;
	 float l1_water = 1/0.2;
	 float c2_water = 21.2753;
	 float l2_water = 1/0.96092;
	 float c3_water = 10.9277;
	 float l3_water = 1/2.97584;
	 float c4_water = 8.5464;
	 float l4_water = 1/10.3271;
	 float c5_water = 8.64322;
	 float l5_water = 1/125.208;
*/
	 
	 float c1_water = 96.63872;
	 float l1_water = 1/0.0975;
	 float c2_water = 25.0472;
	 float l2_water = 1/1.24999;
	 float c3_water = 8.80745;
	 float l3_water = 1/5.7001;
	 float c4_water = 4.19001;
	 float l4_water = 1/10.6501;
	 float c5_water = 9.2732;
	 float l5_water = 1/106.72784;

	 float a01 = 2.277463;
	 float a02 = 0.2431;
	 float a03 = 1.0295;
	 float a04 = 0.4053;
	 float a05 = 6.26751;
	 float a11 = -0.0018473;
	 float a12 = 0.0007;
	 float a13 = -0.00103;
	 float a14 = -0.0007;
	 float a15 = 0.00103;

//	float c1_water = 96.63872 / factor;
//	float l1_water = 1/0.0975 * factor;
//	float c2_water = 25.0472 / factor;
//	float l2_water = 1/1.24999 * factor;
//	float c3_water = 8.80745 / factor;
//	float l3_water = 1/5.7001 * factor;
//	float c4_water = 4.19001 / factor;
//	float l4_water = 1/10.6501 * factor;
//	float c5_water = 9.2732 / factor;
//	float l5_water = 1/106.72784 * factor;

	float energy = 250;
	float range = 37.8225;
	float sumLET = 0;

	double a1, b1, c1, d1, e1, a2, b2, c2, d2, e2;
	double Ez, Ek, dEdz, rz, rzCM;

	float cc1, cc2, cc3, cc4, cc5, pe, phi1, phi2, phi3, phi4, phi5, qp, theta, zmax, taurange, tau0;
	tau0 = 1e-5; // cm

	cout << "Energy is \033[1m " << energy << " MeV\033[0m, calculated range in water is \033[1m " << range << " mm\033[0m.\n";

	for (int i=0; i<8000; i++) {
		float z = i / 200.;
		x[i] = z;
	}

	for (int i=0; i<8000; i++) {
		float z = (float) i/200.;
		
		rz = range - z;

		if (rz<=0) {
			cout << "Break at \033[1m z = " << z << " cm\033[0m.\n";
			cout << "The total energy loss is \033[1m " << sumLET << " MeV\033[0m.\n";
			break;
		}

		a1 = rz * c1_water * exp(-l1_water * rz);
		b1 = rz * c2_water * exp(-l2_water * rz);
		c1 = rz * c3_water * exp(-l3_water * rz);
		d1 = rz * c4_water * exp(-l4_water * rz);
		e1 = rz * c5_water * exp(-l5_water * rz);

		Ez = a1 + b1 + c1 + d1 + e1;

		E[i] = Ez;

		a2 = l1_water * a1;
		b2 = l2_water * b1;
		c2 = l3_water * c1;
		d2 = l4_water * d1;
		e2 = l5_water * e1;

		Ek = a2 + b2 + c2 + d2 + e2;

		dEdz = (Ez / rz - Ek);

		dEdx[i] = dEdz / 200;
		sumLET += dEdz / 200;

		if (abs(rz) < 1) {
			cout << " z = " << z << " cm. -> E = " << Ez << endl;
		}

		// 2nd calculation method
		// Ulmer 2011 Rad Phys and Chem 80 378-389

		cc1 = a01 + a11 * energy;
		cc2 = a02 + a12 * energy;
		cc3 = a03 + a13 * energy;
		cc4 = a04 + a14 * energy;
		cc5 = a05 + a15 * energy;

		taurange = range * (2.117908559e-5 * energy + 0.919238854e-7 * pow(energy, 2)); // 0.41 cm
		zmax = range + taurange;
		pe = (a05 + a15 * energy);

		theta = 1;
		if (rz<=0) theta = 0;
		qp = 3.1415926536 * pe / zmax;

		phi1 = cc1 * exp(-pow(rz, 2) / pow(tau0, 2)) * theta;
		phi2 = 2 * cc2 * theta;
		phi3 = 2 * cc3 * exp(-3.1415*qp * rz * theta);
		phi4 = 2 * cc4 * pow(z / range, 2) * theta;

		Float_t redFactor = 40;

		aPhi1[i] = phi1 / redFactor;
		aPhi2[i] = phi2 / redFactor;
		aPhi3[i] = phi3 / redFactor;
		aPhi4[i] = phi4 / redFactor;
	
		dEdx2[i] = (phi1 + phi2 + phi3*5 + phi4) / redFactor;

		// 3d calculation method
		// Bortfeld

		dEdx3[i] = 1 / (50 * p * pow(a, 1./p) * pow(rz, 1-1./p));
	}

	cout << "zmax is \033[1m " << zmax << " cm \033[0m at max depth.\n";

	TCanvas *c = new TCanvas("c", "E(z) and dE(z)/dz as function of z", 1200, 800);
	TGraph *gEnergy = new TGraph(8000, x, E);
	TGraph *gLet = new TGraph(8000, x, dEdx);
	TGraph *gLet2 = new TGraph(8000, x, dEdx2);
	TGraph *gLet3 = new TGraph(8000, x, dEdx3);

	TGraph *gPhi1 = new TGraph(8000, x, aPhi1);
	TGraph *gPhi2 = new TGraph(8000, x, aPhi2);
	TGraph *gPhi3 = new TGraph(8000, x, aPhi3);
	TGraph *gPhi4 = new TGraph(8000, x, aPhi4);

	gEnergy->SetLineColor(kBlack);
	gLet->SetLineColor(kBlue);
	gLet2->SetLineColor(kRed);
	gLet3->SetLineColor(kGreen);

	gPhi1->SetLineColor(kBlack);
	gPhi2->SetLineColor(kBlack);
	gPhi3->SetLineColor(kBlack);
	gPhi4->SetLineColor(kBlack);

	gPad->SetLogy();

//	gEnergy->Draw("");
	gLet->Draw("");
	gLet2->Draw("same");
	gLet3->Draw("same");
//	gPhi1->Draw("same");
//	gPhi2->Draw("same");
//	gPhi3->Draw("same");
//	gPhi4->Draw("same");
}
