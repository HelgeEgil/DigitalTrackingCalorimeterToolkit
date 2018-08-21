#ifndef Misalign_cxx
#define Misalign_cxx

#include <iostream>
#include <fstream>
#include <algorithm>

#include <TObject.h>

#include "GlobalConstants/Misalign.h"
#include "Classes/Cluster/Cluster.h"
#include "Classes/Cluster/Clusters.h"

using namespace DTC;

Misalign::Misalign() {
   // Data retrieved from dropbox.com/sh/zv5cu0j1pjgb2ph/AAD2Yg72S9uek7mJN-P7TsE6a/backupgeometry.txt
   // I.e. the Utrecht experiment. Their unit is cm, I've multiplied every number with 10

   /*
   Float_t a[7] = { 0.000,  0.101,  0.014,  0.004, -0.231,  0.162, -0.098}; 
   Float_t b[7] = { 0.000,  0.101,  0.071, -0.548, -0.128,  0.095,  0.000};
   Float_t c[7] = { 0.000,  0.064,  0.244, -0.019,  0.020, -0.110, -0.141};
   Float_t d[7] = { 0.000,  0.156,  0.343,  0.191,  0.393,  0.356,  0.341};

   for (int i=1; i<24; i++) {
      misalignX[0][i] = 0; misalignX[1][i] = 0;
      misalignY[0][i] = 0; misalignY[1][i] = 0;
   }

   std::copy(a, a+7, misalignX[0]);
   std::copy(b, b+7, misalignX[1]);
   std::copy(c, c+7, misalignY[0]);
   std::copy(d, d+7, misalignY[1]);

   */
   
   chipAlignment chip;

   ifstream in;
   in.open("Data/ExperimentalData/Alignment.txt");
   Int_t nlines = 0;
   while (1) {
      
      in >> chip.idx >> chip.deltaX >> chip.deltaY >> chip.deltaTheta;

      if (!in.good()) break;

      chip.deltaX *= 10;
      chip.deltaY *= 10;

      chipAlignmentArray_[chip.idx] = chip;
   }
   in.close();
}

Misalign::~Misalign() {
}

void Misalign::correctClusters(Clusters * clusters) {
   Cluster      * cluster = nullptr;
   Float_t        oldX, oldY, newX, newY;
   Float_t        deltaX, deltaY, theta;
   Int_t          layer;
   chipAlignment  align;

   for (Int_t i=0; i<clusters->GetEntriesFast(); i++) {
      cluster = clusters->At(i);
      if (!cluster) continue;

      oldX  = cluster->getXmm();
      oldY  = cluster->getYmm();
      align = getMisalign(cluster);
      
      deltaX = align.deltaX;
      deltaY = align.deltaY;
      theta  = align.deltaTheta;

      newX = deltaX + oldX;
      newY = deltaY + oldY;

      newX = newX * cos(theta) - newY * sin(theta);
      newY = newX * sin(theta) + newY * cos(theta);

      cluster->setXmm(newX);
      cluster->setYmm(newY);
   }
}

chipAlignment Misalign::getMisalign(Cluster * cluster) {
    Int_t chipIdx = cluster->getChip();

    return chipAlignmentArray_[chipIdx];
}

#endif
