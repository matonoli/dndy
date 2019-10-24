AliAnalysisTask *AddTaskTriggerEfficiency()
{
  // specify appendix: track cuts or differences in the tasks
  TString appendix("TriggerEff");
  
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  
  //create output containers
  TString containerName = mgr->GetCommonFileName();
  containerName += ":TriggerEff"; // create a subfolder in the file
  // Appendix to add the tack cuts to the folder name
  printf("container name: %s\n", containerName.Data());
  
  AliAnalysisTaskTriggerEfficiency *task = new AliAnalysisTaskTriggerEfficiency(Form("Task%s",appendix.Data()));
  //  configure task here
  // Trigger selection
  task->SelectCollisionCandidates(AliVEvent::kAnyINT);
  mgr->AddTask(task);
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  // create output containers for Results and QA
  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("Results_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("QA_%s", appendix.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, containerName.Data()));
  
  return task;
}
