/*
 * Tools.h
 *
 *  Created on: Sep 2, 2015
 *      Author: local1
 */

#include <vector>
#include <TObject.h>
using namespace std;

#ifndef FOCALCODE_TOOLS_H_
#define FOCALCODE_TOOLS_H_

Bool_t isItemInVector(Int_t i, vector<Int_t> *v);
Float_t diffmm(Cluster *p1, Cluster *p2);
Float_t diffmm(Hit *h1, Hit *h2);

#endif /* FOCALCODE_TOOLS_H_ */
