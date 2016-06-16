void Load() {

   gROOT->LoadMacro("GlobalConstants/MaterialConstants.C+");
   gROOT->ProcessLine("MaterialConstants()");

	gROOT->LoadMacro("RootFiles/Wrapper.C+");
	gROOT->LoadMacro("Analysis/Analysis.C+");
	gROOT->LoadMacro("UnitTests/UnitTests.C+");
}
