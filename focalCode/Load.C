void Load() {
	
	// gDebug =3 ;

   gROOT->LoadMacro("GlobalConstants/MaterialConstants.C+");
   gROOT->ProcessLine("MaterialConstants()");

	gROOT->LoadMacro("RootFiles/Wrapper.C+");
	gROOT->LoadMacro("HelperFunctions/getTracks.C+");
	gROOT->LoadMacro("Analysis/Analysis.C+");
}
