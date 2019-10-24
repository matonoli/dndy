//__________________________________________________________________________________________________
void runHelium3Analysis ()  {
    
    //Load Macro
    gROOT->LoadMacro("helium3Analysis.cxx");

   const char *inputFileName  = "../Trees/ReducedTreeNuclei_3He.root"; // path to the output of the reduced tree frame work for nuclei
    const char *inputListName  = "Results"; // Name of the list in the rootfile which containes the tree folder
    const char *outputFileName = "../Helium_LHC16q_FAST_woSDD.root";
    
    //Execute Analysis
    helium3Analysis analysis;
    analysis.SetFileNames (inputFileName,inputListName,outputFileName);
    analysis.GetInputTree();
    analysis.SetTreeBranches();
    analysis.CreateHistograms();
    analysis.ExecuteAnalysis();
    analysis.WriteOutputFile();
}
//__________________________________________________________________________________________________
