#define UnitTests_cxx

#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TH2.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis3D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TStopwatch.h>
#include <TPaveStats.h>
#include <TView.h>
#include <TLeaf.h>
#include <TArrow.h>
#include <TF1.h>
#include <Math/ProbFunc.h>

#include "UnitTests/UnitTests.h"
#include "GlobalConstants/Constants.h"
#include "GlobalConstants/MaterialConstants.h"
#include "Classes/Hit/Hits.h"
#include "Classes/Cluster/Clusters.h"
#include "Classes/Track/Tracks.h"
#include "Classes/Track/conversionFunctions.h"
#include "Classes/DataInterface/DataInterface.h"
#include "HelperFunctions/Tools.h"

using namespace std;

void testHits() {
	Hits * someHits = new Hits();
	Hits * fullHits = new Hits();
	Hits * overloadHits = new Hits();
	Hits * errorHits = new Hits();

	int iSome = 0, iFull = 0, iOverload = 0, iError = 0;

	TRandom3 *gRandom = new TRandom3();

	for ( Int_t i=0; i<nLayers; i++) {

		cout << "Layer " << i << endl;

		if (i%2) {
			cout << "Adding Hit* to someHits\n";
			someHits->appendPoint(10, 10, i, 1, 5);
		}

		cout << "Adding Hit* to fullHits\n";
		fullHits->appendPoint(10, 10, i, 1, 5);
		cout << "Adding Hit* to overloadHits\n";
		overloadHits->appendPoint(10, 10, i, 1, 5);
		cout << "Adding Hit* to errorHits\n";
		errorHits->appendPoint(10, 10, i, 1, 5);
		
		if (i%3 != 0) {
			cout << "Adding 2x Hit* to overloadHits\n";
			overloadHits->appendPoint(20, 20, i, 1, 5);
			cout << "Adding 2x (wrong) Hit* to errorHits\n";
			errorHits->appendPoint(10, 10, i, 2, 5);
		}

	}


	cout << "Summing someHits " << endl;
	iSome = someHits->sumAllHitsInEachLayer();
	cout << "Summing fullHits " << endl;
	iFull = fullHits->sumAllHitsInEachLayer();
	cout << "Summing overloadHits " << endl;
	iOverload = overloadHits->sumAllHitsInEachLayer();
	cout << "Summing errorHits " << endl;
	iError = errorHits->sumAllHitsInEachLayer();

	cout << "someHits = " << iSome << endl;
	cout << "fullHits = " << iFull << endl;
	cout << "overloadHits = " << iOverload << endl;
	cout << "errorHits = " << iError << endl;

	cout << "Coordinate list for all some, full, overload, error: \n";
	for (Int_t i=0; i<overloadHits->GetEntriesFast(); i++) {
		if (someHits->At(i)) cout << *someHits->At(i) << ", ";
		if (fullHits->At(i)) cout << *fullHits->At(i) << ", ";
		if (overloadHits->At(i)) cout << *overloadHits->At(i) << ", ";
		if (errorHits->At(i)) cout << *errorHits->At(i) << ", ";
		cout << endl;
	}
	
	cout << "Edep list for all some, full, overload, error: \n";
	for (Int_t i=0; i<overloadHits->GetEntriesFast(); i++) {
		if (someHits->At(i)) cout << someHits->getEdep(i) << ", ";
		if (fullHits->At(i)) cout << fullHits->getEdep(i) << ", ";
		if (overloadHits->At(i)) cout << overloadHits->getEdep(i) << ", ";
		if (errorHits->At(i)) cout << errorHits->getEdep(i) << ", ";
		cout << endl;
	}
}
