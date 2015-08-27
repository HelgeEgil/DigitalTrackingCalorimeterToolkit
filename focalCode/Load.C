void Load() {
/*
   gROOT->LoadMacro("Hit.C+");
   gROOT->LoadMacro("HitCollection.C+");
   gROOT->LoadMacro("Cluster.C+");
   gROOT->LoadMacro("ClusterCollection.C+");
   gROOT->LoadMacro("Track.C+");
   gROOT->LoadMacro("TrackCollection.C+");
   gROOT->LoadMacro("TFocal.C+");
*/
	gSystem->Load("libMathMore");
	gROOT->LoadMacro("wrapper.C+");
	gROOT->LoadMacro("Analysis.C+");
}
