void Load() {
	
	// gDebug =3 ;

   cout << "LOADING MATERIALCONSTANTS\n";
   gROOT->LoadMacro("GlobalConstants/MaterialConstants.C+");
   gROOT->ProcessLine("MaterialConstants()");

   cout << "LOADING WRAPPER\n";
	gROOT->LoadMacro("RootFiles/Wrapper.C+");
   cout << "LOADING GETTRACKS\n";
	gROOT->LoadMacro("HelperFunctions/getTracks.C+");
   cout << "LOADING ANALYSIS\n";
	gROOT->LoadMacro("Analysis/Analysis.C+");
}
