#define ALIPHYSICS_VERSION "vAN-20190417_ROOT6-1"

#ifdef __CLING__
// ROOT6 specific includes
// Tell  ROOT where to find AliRoot headers
R__ADD_INCLUDE_PATH($ALICE_ROOT)
//#include <ANALYSIS/macros/AddTaskPIDResponse.C>

// Tell ROOT where to find AliPhysics headers
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "OADB/macros/AddTaskPhysicsSelection.C"
#include "OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"

// Load local macros
#include "AliAnalysisTaskTriggerEfficiency.h"

#endif

void RunTaskOnGrid()
{
  // Firstly, set some variables
  const char* analysisMode = "grid"; // "grid" or "proof"
  const char*  mode = "test"; // "test" "full" "terminate"
  const char*  user = "shornung";
  const char* work_dir = "TriggerEfficiency";
  const char* output_dir = "LHC18f3";
  Bool_t IsAOD = kTRUE;
  Bool_t IsMC = kTRUE;
  Bool_t JDLMerge = true;
  
  Int_t noffiles = 20; // Number of files per job (50)
  Int_t ttl = 15000; // maximum run time in seconds //10000
  Int_t cyclenumber = 1;
  
  // create and customize the alien handler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
  
  alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  alienHandler->SetAdditionalLibs("AliAnalysisTaskTriggerEfficiency.cxx AliAnalysisTaskTriggerEfficiency.h");
  alienHandler->SetAnalysisSource("AliAnalysisTaskTriggerEfficiency.cxx");
  // select the aliphysics version. all other packages
  // are LOADED AUTOMATICALLY!
  alienHandler->SetAliPhysicsVersion(ALIPHYSICS_VERSION);
  
  //MinBias (1:1 Stat): LHC18f3
  // FAST
  alienHandler->SetGridDataDir("/alice/sim/2018/LHC18f3_fast_1");
  //  alienHandler->SetGridDataDir("/alice/sim/2018/LHC18f3_fast_2");
  //CENT without SDD
  //         alienHandler->SetGridDataDir("/alice/sim/2018/LHC18f3_cent_woSDD_1");
  //   alienHandler->SetGridDataDir("/alice/sim/2018/LHC18f3_cent_woSDD_2");
  alienHandler->SetDataPattern("*/AOD202/*AOD.root"); // on-going refiltering
  
  
  if (!IsMC) {
    alienHandler->SetRunPrefix("000"); // IMPORTANT! TO BE USED FOR REAL DATA
    
    printf("!!! Task only for MC !!!\n");
    return -1;
  }
  
  Int_t runcycle[] = {0,32}; // {First run , last run (number of runs)}
  //LHC16q: 32 entries - RunList_LHC16q_pass1_CentralBarrelTracking_20170318_v1.txt: SSD SPD SDD V0 ZDC TPC
  Int_t runArray[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265335, 265334, 265332, 265309};
  
  //  Int_t runcycle[] = {0,31}; // {First run , last run (number of runs)} //LHC16q
  //  //LHC16q: 31 entries - RunList_LHC16q_pass1_CentralBarrelTracking_hadronPID_20171129_v2.txt: SSD SPD SDD V0 ZDC TPC TOF
  //  Int_t runArray[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};
  
  for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
  {
    if (i == sizeof(runArray) / sizeof(runArray[1])) break;
    alienHandler->AddRunNumber(runArray[i]);
  }
  
  alienHandler->SetSplitMaxInputFileNumber(noffiles);
  alienHandler->SetExecutable("TriggerEfficiency.sh");
  alienHandler->SetTTL(ttl);
  alienHandler->SetJDLName("TriggerEfficiency.jdl");
  alienHandler->SetOutputToRunNo(kTRUE);
  alienHandler->SetKeepLogs(kTRUE);
  
  // Merging: run with kTRUE to merge on grid
  // After re-running the jobs in SetRunMode("terminate")
  // set SetMergeViaJDL(kFALSE) to collect final results
  alienHandler->SetMaxMergeStages(1);
  alienHandler->SetMergeViaJDL(JDLMerge);
  // define the output folders
  alienHandler->SetGridWorkingDir(work_dir);
  alienHandler->SetGridOutputDir(output_dir);
  
  alienHandler->SetOverwriteMode();
  alienHandler->SetUser(user);
  alienHandler->SetAnalysisMacro("TriggerEffOnGRIDmacro.C");
  alienHandler->SetMasterResubmitThreshold(10); //resubmit if less than 10% successfull
  alienHandler->SetExecutableCommand("aliroot -b -q");
  alienHandler->SetInputFormat("xml-single");
  alienHandler->SetPrice(1);
  alienHandler->SetSplitMode("se");
  alienHandler->SetMergeExcludes("EventStat_temp.root");
  
  if (!strcmp(mode,"test")) { // check if mode is test mode via comparison of strings
    alienHandler->SetNtestFiles(1);
  }
  alienHandler->SetRunMode(mode);
  
  if (!alienHandler) return;
  
  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("TriggerEfficiency");
  mgr->SetGridHandler(alienHandler);
  
  if(IsAOD){
    AliAODInputHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
  } else{
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    if(IsMC){
      AliMCEventHandler* mcHandler = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
    }
  }
  
  // Physics selection task
  AliPhysicsSelectionTask *PhySel = AddTaskPhysicsSelection(IsMC,kTRUE);
  
  // Multiplicity task
  AliMultSelectionTask* multSelectionTask = AddTaskMultSelection();
  multSelectionTask->SelectCollisionCandidates(AliVEvent::kAnyINT);
  
  // since we will compile a class, tell root where to look for headers
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  
  // compile the class and load the add task macro
  gInterpreter->LoadMacro("AliAnalysisTaskTriggerEfficiency.cxx++g");
  AliAnalysisTaskTriggerEfficiency *task = reinterpret_cast<AliAnalysisTaskTriggerEfficiency*>(gInterpreter->ExecuteMacro("AddTaskTriggerEfficiency.C"));
  
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  //  mgr->SetDebugLevel(2);
  // Start analysis in grid.
  mgr->StartAnalysis(analysisMode);
  
}
