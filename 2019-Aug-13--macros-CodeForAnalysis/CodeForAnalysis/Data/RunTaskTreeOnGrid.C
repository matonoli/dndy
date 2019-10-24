#define ROOT_VERSION       "v5-34-30"
#define AliROOT_VERSION    "v5-09-09-1"
#define ALIPHYSICS_VERSION "vAN-20190109-1" // if AliPhysics version and local version differs change AddTask and AliAnalysisTask name


void RunTaskTreeOnGrid()
{
  // Firstly, set some variables
  const char* analysisMode = "grid"; // "grid" or "proof"
  const char*  mode = "terminate"; // "test" "full" "terminate"
  const char*  user = "shornung";
  char* work_dir = "NucleiTaskTree";
  char* output_dir = "pPb_3H_TOFlist_Cent";
  Bool_t enablePileupCuts=kTRUE;
  Bool_t IsAOD = kTRUE;
  Bool_t IsMC = false;
  Bool_t JDLMerge = false;


  Int_t noffiles = 50; // Number of files per job
  Int_t ttl = 12000; // maximum run time in seconds //10000
  Int_t cyclenumber = 1;

  // create and customize the alien handler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

  alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  alienHandler->SetAdditionalLibs("AliAnalysisTaskReducedTreeNuclei_Local.cxx AliAnalysisTaskReducedTreeNuclei_Local.h");
  alienHandler->SetAnalysisSource("AliAnalysisTaskReducedTreeNuclei_Local.cxx");
  // select the aliphysics version. all other packages
  // are LOADED AUTOMATICALLY!
  alienHandler->SetAliPhysicsVersion(ALIPHYSICS_VERSION);
  //   alienHandler->SetROOTVersion(ROOT_VERSION);
  //   alienHandler->SetAliROOTVersion(AliROOT_VERSION);

  // Data
  alienHandler->SetGridDataDir("/alice/data/2016/LHC16q");
  //CENT with SDD
  //      alienHandler->SetDataPattern("*/pass1_CENT_wSDD/AOD190/*AOD.root");
  //CENT without SDD
  alienHandler->SetDataPattern("*/pass1_CENT_woSDD/AOD190/*AOD.root");
  //FAST (without SDD)
  //  alienHandler->SetDataPattern("*/pass1_FAST/AOD190/*AOD.root");

  // Data
  //   alienHandler->SetGridDataDir("/alice/data/2016/LHC16t");
  //CENT with SDD
  //         alienHandler->SetDataPattern("*/pass1_CENT_wSDD/AOD190/*AOD.root");
  //CENT without SDD
  //   alienHandler->SetDataPattern("*/pass1_CENT_woSDD/AOD190/*AOD.root");
  //FAST (without SDD)
  //   alienHandler->SetDataPattern("*/pass1_FAST/AOD190/*AOD.root");

  // MC
  //Genereal purpose: LHC17f2a_fast_fix // LHC17f2a_cent_woSDD_fix // LHC17f2a_cent_fix
  //      alienHandler->SetGridDataDir("/alice/sim/2017/LHC17f2a_fast_fix");
  //LHC+Nuclei: LHC17d10_cent LHC17d10_fast
  //FAST (without SDD)
  //      alienHandler->SetGridDataDir("/alice/sim/2017/LHC17d10_fast");
  //   //CENT with SDD
  //   alienHandler->SetGridDataDir("/alice/sim/2017/LHC17d10_cent");
  //      alienHandler->SetDataPattern("*/AOD/*AOD.root");

  if (!IsMC) {
    alienHandler->SetRunPrefix("000"); // IMPORTANT! TO BE USED FOR REAL DATA
  }

  //  Int_t runcycle[] = {0,32}; // {First run , last run (number of runs)} //LHC16q
  //  //LHC16q: 32 entries - RunList_LHC16q_pass1_CentralBarrelTracking_20170318_v2.txt: SSD SPD SDD V0 ZDC TPC
  //  Int_t runArray[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265335, 265334, 265332, 265309};

  Int_t runcycle[] = {0,31}; // {First run , last run (number of runs)} //LHC16q
  //LHC16q: 31 entries - RunList_LHC16q_pass1_CentralBarrelTracking_hadronPID_20171129_v2.txt: SSD SPD SDD V0 ZDC TPC TOF
  Int_t runArray[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265334, 265332, 265309};

  //   Int_t runcycle[] = {0,4}; // {First run , last run (number of runs)} // LHC16t
  //   //   LHC16t: 4 entries - RunList_LHC16t_pass1_CentralBarrelTracking_20170202_v0.txt: SSD SPD SDD V0 ZDC TPC
  //   Int_t runArray[] = {267166, 267165, 267164, 267163};

  for (Int_t i =  runcycle[cyclenumber - 1]; i < runcycle[cyclenumber] ; i++)
  {
    if (i == sizeof(runArray) / sizeof(runArray[1])) break;
    alienHandler->AddRunNumber(runArray[i]);
  }

  alienHandler->SetSplitMaxInputFileNumber(noffiles);
  alienHandler->SetExecutable("NucleiTree.sh");
  alienHandler->SetTTL(ttl);
  alienHandler->SetJDLName("NucleiTree.jdl");
  alienHandler->SetOutputToRunNo(kTRUE);
  alienHandler->SetKeepLogs(kTRUE);

  // Merging: run with kTRUE to merge on grid
  // After re-running the jobs in SetRunMode("terminate")
  // set SetMergeViaJDL(kFALSE) to collect final results
  alienHandler->SetMergeViaJDL(JDLMerge);
  alienHandler->SetMaxMergeStages(1);
  //   alienHandler->SetMaxMergeFiles(50);
  //   alienHandler->SetMergeExcludes("AnalysisResults.root");
  // define the output folders
  alienHandler->SetGridWorkingDir(work_dir);
  alienHandler->SetGridOutputDir(output_dir);

  alienHandler->SetOverwriteMode();
  alienHandler->SetUser(user);
  alienHandler->SetAnalysisMacro("NucleiOnGRIDmacro.C");
  //   alienHandler->SetMasterResubmitThreshold(10); //resubmit if less than 10% successfull
  alienHandler->SetExecutableCommand("aliroot -b -q");
  alienHandler->SetInputFormat("xml-single");
  alienHandler->SetPrice(1);
  alienHandler->SetSplitMode("se");
  alienHandler->SetMergeExcludes("EventStat_temp.root");

  if (mode== "test") {
    alienHandler->SetNtestFiles(1);
  }
  alienHandler->SetRunMode(mode);

  if (!alienHandler) return;

  // header location
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("NucleiTree");
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
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *PhySel = AddTaskPhysicsSelection(IsMC,enablePileupCuts);

  // PID response task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  // IsMC, autoMCesd, tuneOnData, recopass, chachePID, detResponse, useTPCetaCorrection, useTPCMultiplicityCorrection, recoDataPass
  AddTaskPIDResponse(IsMC);

  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
//  AliMultSelectionTask *MultTask = AddTaskMultSelection();
//  MultTask->SelectCollisionCandidates(AliVEvent::kINT7);
  AddTaskMultSelection();

  // compile the class (locally)
  gROOT->LoadMacro("AliAnalysisTaskReducedTreeNuclei_Local.cxx++g");
  // load the addtask macro
  gROOT->LoadMacro("AddTaskReducedTreeNuclei_Local.C");
  // create an instance of your analysis task – adjust: data/MC , AOD/ESD
  AliAnalysisTaskReducedTreeNuclei_Local *task = AddTaskReducedTreeNuclei_Local("LHC16q");
  if (!mgr->InitAnalysis()) return;

  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis(analysisMode);

}
