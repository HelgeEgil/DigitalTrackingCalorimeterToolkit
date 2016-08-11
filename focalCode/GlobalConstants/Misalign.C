#include <iostream>
#include <algorithm>

#include <TObject.h>

#include "GlobalConstants/Misalign.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"

Misalign::Misalign() {
	// Data retrieved from dropbox.com/sh/zv5cu0j1pjgb2ph/AAD2Yg72S9uek7mJN-P7TsE6a/backupgeometry.txt
	// I.e. the Utrecht experiment. Their unit is cm, I've multiplied every number with 10

	Float_t a[24] = {-0.086390, -0.001186, -0.006981, -0.120161, -0.349799, -0.086450,
					  	  -0.192518, -0.231557, -0.223723, -0.048009, -0.054521, -0.191975,
						  -0.246765, -0.207292, -0.414540, -0.311817, -0.388111, -0.345362,
						  -0.338518,  0		 , -0.491954,  0		  , -0.394087, -0.361302};

	Float_t b[24] = { 0.329784,  0.424458,  0.448127,  0.336280,  0.223954,  0.470487,
						   0.304410,  0.230047,  0.530765,  0.246614,  0.258031,  0.192016,
						   0.256506,  0.142987,  0.145595,  0.271740,  0.249492,  0.155917,
						   0.220785,  0.248612,  0.227604,  0       ,  0.221386,  0.202168};

	Float_t c[24] = {-0.657406, -0.400492, -0.585275, -0.362801, -0.384095, -0.406281,
					     -0.263581, -0.402140, -0.536376, -0.279049, -0.301387, -0.005010,
						  -0.198595,  0.082010, -0.308502, -0.220456, -0.211355, -0.337013,
						  -0.560832,  0       , -0.416552,  0       , -0.233913, -0.150391};

	Float_t d[24] = {-0.563372, -0.562861, -0.396891, -0.513463, -0.282952, -0.334688,
						  -0.440048, -0.380420, -0.227459, -0.450782, -0.219784, -0.235561,
						  -0.244287, -0.191445, -0.362194, -0.190465, -0.116786, -0.436257,
						  -0.222045, -0.080618, -0.301860,  0       , -0.130745, -0.074902};

	// Corrected values using my tests
	Float_t a_corr[7] = { 0.000,  0.101,  0.014,  0.004, -0.231,  0.162, -0.098}; 
	Float_t b_corr[7] = { 0.000,  0.101,  0.071, -0.548, -0.128,  0.095,  0.000};
	Float_t c_corr[7] = { 0.000,  0.064,  0.244, -0.019,  0.020, -0.110, -0.141};
	Float_t d_corr[7] = { 0.000,  0.156,  0.343,  0.191,  0.393,  0.356,  0.341};

	Float_t mean_1 = 0, mean_2 = 0;
	for (Int_t i=1; i<7; i++) {
		mean_1 += sqrt(pow(a_corr[i], 2) + pow(c_corr[i], 2)) / 6.;
		mean_2 += sqrt(pow(b_corr[i], 2) + pow(d_corr[i], 2)) / 6.;
	}

	std::copy(a, a+24, misalignX[0]);
	std::copy(b, b+24, misalignX[1]);
	std::copy(c, c+24, misalignY[0]);
	std::copy(d, d+24, misalignY[1]);
	
	for (int i=1; i<24; i++) {
		misalignX[0][i] = 0; misalignX[1][i] = 0;
		misalignY[0][i] = 0; misalignY[1][i] = 0;
	}

	std::copy(a_corr, a_corr+7, misalignX[0]);
	std::copy(b_corr, b_corr+7, misalignX[1]);
	std::copy(c_corr, c_corr+7, misalignY[0]);
	std::copy(d_corr, d_corr+7, misalignY[1]);
}

Misalign::~Misalign() {
}

void Misalign::correctClusters(Clusters * clusters) {
	Cluster *cluster = nullptr;
	Float_t	oldX, oldY, newX, newY;
	Int_t		layer;
	point	  align;

	for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
		cluster = clusters->At(i);
		if (!cluster) continue;

		oldX = cluster->getXmm();
		oldY = cluster->getYmm();
		align = getMisalign(cluster);

		newX = oldX + align.x;
		newY = oldY + align.y;

		cluster->setXmm(newX);
		cluster->setYmm(newY);
	}
}

point Misalign::getMisalign(Cluster * cluster) {
	Float_t	y = cluster->getYmm();
	Int_t		assembly = (y>0) ? 0 : 1;
	Int_t		layer = cluster->getLayer();

	return getMisalign(assembly, layer);
}

point Misalign::getMisalign(Int_t assembly, Int_t layer) {
	point p;
	p.x = misalignX[assembly][layer];
	p.y = misalignY[assembly][layer];
	return p;
}
