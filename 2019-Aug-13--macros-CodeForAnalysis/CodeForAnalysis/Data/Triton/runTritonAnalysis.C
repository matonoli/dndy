//__________________________________________________________________________________________________
void runTritonAnalysis ()  {
    
    //Load Macro
    gROOT->LoadMacro("TritonAnalysis.cxx");

   const char *inputFileName  = "../Trees/ReducedTreeNuclei_Triton_TOFlist.root"; // path to the output of the reduced tree frame work for nuclei
    const char *inputListName  = "Results"; // Name of the list in the rootfile which containes the tree folder
    const char *outputFileName = "../Triton_LHC16q_FAST_woSDD.root";
    
    //Execute Analysis
    TritonAnalysis analysis;
    analysis.SetFileNames (inputFileName,inputListName,outputFileName);
    analysis.GetInputTree();
    analysis.SetTreeBranches();
    analysis.CreateHistograms();
    analysis.ExecuteAnalysis();
    analysis.WriteOutputFile();
}
//__________________________________________________________________________________________________
