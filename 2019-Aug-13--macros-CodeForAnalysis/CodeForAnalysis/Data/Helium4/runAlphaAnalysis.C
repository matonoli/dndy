//__________________________________________________________________________________________________
void runAlphaAnalysis ()  {
    
    //Load Macro
    gROOT->LoadMacro("AlphaAnalysis.cxx");

   const char *inputFileName  = "../Trees/ReducedTreeNuclei_TOFlist_He4.root"; // path to the output of the reduced tree frame work for nuclei
    const char *inputListName  = "Results"; // Name of the list in the rootfile which containes the tree folder
    const char *outputFileName = "../Alpha_LHC16q_FAST_woSDD.root";
    
    //Execute Analysis
    AlphaAnalysis analysis;
    analysis.SetFileNames (inputFileName,inputListName,outputFileName);
    analysis.GetInputTree();
    analysis.SetTreeBranches();
    analysis.CreateHistograms();
    analysis.ExecuteAnalysis();
    analysis.WriteOutputFile();
}
//__________________________________________________________________________________________________
