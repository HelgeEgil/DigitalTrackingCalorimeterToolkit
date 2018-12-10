#include <TROOT.h>

using namespace std;

void drawTracksRangeHistogramScript(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness, Int_t eventsPerRun, Int_t outputFileIdx, Bool_t drawFitResults, Bool_t doTracking, Bool_t excludeNuclearInteractions, Int_t skipTracks);

void drawTracksRangeHistogramScript(Int_t Runs, Int_t dataType, Bool_t recreate, Float_t energy, Float_t degraderThickness, Int_t eventsPerRun, Int_t outputFileIdx, Bool_t drawFitResults, Bool_t doTracking, Bool_t excludeNuclearInteractions, Int_t skipTracks) {
	gROOT->ProcessLine(".x Load.C");
	gROOT->ProcessLine(Form("drawTracksRangeHistogram(%d, %d, %d, %.0f, %.0f, %d, %d, %d, %d, %d, %d)", Runs, dataType, recreate, energy, degraderThickness, eventsPerRun, outputFileIdx, drawFitResults, doTracking, excludeNuclearInteractions, skipTracks)); 
	gROOT->ProcessLine(".q");
}
