void Load() {

   gROOT->LoadMacro("MaterialConstants.C+");
   gROOT->ProcessLine("MaterialConstants()");
   
//    gROOT->LoadMacro("Tools.C+");
//    gROOT->LoadMacro("Hit.C+");
//    gROOT->LoadMacro("Hits.C+");
//    gROOT->LoadMacro("Cluster.C+");
//    gROOT->LoadMacro("Clusters.C+");
//    gROOT->LoadMacro("Cluster_findTracks.C+");
//    gROOT->LoadMacro("Track.C+");
//    gROOT->LoadMacro("Tracks.C+");
//    gROOT->LoadMacro("TFocal.C+");
//    gROOT->LoadMacro("CalorimeterFrame.C+");
//    gROOT->LoadMacro("TrackerFrame.C");

//	gSystem->Load("libMathMore");
	gROOT->LoadMacro("wrapper.C+");
	gROOT->LoadMacro("Analysis.C+");
}
