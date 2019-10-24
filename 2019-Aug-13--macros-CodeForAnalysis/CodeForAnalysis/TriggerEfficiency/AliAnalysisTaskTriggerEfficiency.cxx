/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
// Task for the analysis of the trigger efficiency as a function of multiplicity
//
// Author:
//  Sebastian Hornung <Sebastian.Hornung@cern.ch>

// Root header
#include <TChain.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TH1D.h>

//AliRoot header
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"

//AliPhysics header
#include "AliAnalysisTaskTriggerEfficiency.h"

ClassImp(AliAnalysisTaskTriggerEfficiency) // classimp: necessary for root

//____________________________________________________________
AliAnalysisTaskTriggerEfficiency::AliAnalysisTaskTriggerEfficiency():
AliAnalysisTaskSE(),
MultiplicityPercentile(-1),
hGeneratedEvents(nullptr),
hTriggeredEvents(nullptr),
hSelectedEvents(nullptr),
hGeneratedSignal(nullptr),
hGeneratedSignalEventTriggered(nullptr),
hGeneratedSignalEventSelected(nullptr),
fOutput(nullptr),
fQAList(nullptr),
fAODMCHeader(nullptr),
fAODArrayMCParticles(nullptr),
fEventCuts()
{
}

//____________________________________________________________
AliAnalysisTaskTriggerEfficiency::AliAnalysisTaskTriggerEfficiency(const char *name):
AliAnalysisTaskSE(name),
MultiplicityPercentile(-1),
hGeneratedEvents(nullptr),
hTriggeredEvents(nullptr),
hSelectedEvents(nullptr),
hGeneratedSignal(nullptr),
hGeneratedSignalEventTriggered(nullptr),
hGeneratedSignalEventSelected(nullptr),
fOutput(nullptr),
fQAList(nullptr),
fAODMCHeader(nullptr),
fAODArrayMCParticles(nullptr),
fEventCuts()
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________
AliAnalysisTaskTriggerEfficiency::~AliAnalysisTaskTriggerEfficiency(){
  //
  // Destructor
  //
  
  // Delete output objects only if we are not running in PROOF mode because otherwise this produces a crash during merging
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(mgr && mgr->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
    if(fOutput) delete fOutput;
    if(fQAList) delete fQAList;
  }
}

//___________________________________________________
void AliAnalysisTaskTriggerEfficiency::UserCreateOutputObjects(){
  
  //Objects that are output are initialized in the function
  AliDebug(4, "Creating Output Objects");
  
  // Make lists for Output
  if(!fOutput) fOutput = new TList();
  fOutput->SetOwner();
  
  if(!fQAList) fQAList = new TList();
  fQAList->SetOwner();
  
  // QA histograms
  fEventCuts.AddQAplotsToList(fQAList); /// Add event selection QA plots
  
  // Trigger / event selection efficiency
  hGeneratedEvents = new TH1D ("hGeneratedEvents","Generated events vs multiplicity percentile",100,0.,100.);
  hTriggeredEvents = new TH1D ("hTriggeredEvents","Triggered events vs multiplicity percentile",100,0.,100.);
  hSelectedEvents = new TH1D ("hSelectedEvents","Selected events vs multiplicity percentile",100,0.,100.);

// Signal loss correction
  hGeneratedSignal = new TH2D ("hGeneratedSignal","Generated signal vs pt and multiplicity percentile",100,0.,100.,100,0,10);
  hGeneratedSignalEventTriggered = new TH2D ("hGeneratedSignalEventTriggered","Generated signal vs pt and multiplicity percentile (triggered events)",100,0.,100.,100,0,10);
  hGeneratedSignalEventSelected = new TH2D ("hGeneratedSignalEventSelected","Generated signal vs pt and multiplicity percentile (selected events)",100,0.,100.,100,0,10);
  
  fOutput->Add(hGeneratedEvents);
  fOutput->Add(hTriggeredEvents);
  fOutput->Add(hSelectedEvents);

  fOutput->Add(hGeneratedSignal);
  fOutput->Add(hGeneratedSignalEventTriggered);
  fOutput->Add(hGeneratedSignalEventSelected);
  
  PostData(1,fOutput);
  PostData(2,fQAList);
}

//___________________________________________________
void AliAnalysisTaskTriggerEfficiency::UserExec(Option_t *){ //called for each event
  
  AliDebug(4, "Starting Single Event Analysis");
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return;
  }
  
  // Get MC info (Task for MC studies only)
  fAODMCHeader = dynamic_cast<AliAODMCHeader *>(fInputEvent->FindListObject(AliAODMCHeader::StdBranchName()));
  if(!fAODMCHeader){
    AliError("No AliAODMCHeader");
    return;
  }
  fAODArrayMCParticles = dynamic_cast<TClonesArray *>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(!fAODArrayMCParticles){
    AliError("No AOD MC particles");
    return;
  }
  
  //   Multiplicity / centrality framework
  AliMultSelection *MultSelection = static_cast<AliMultSelection*>(fInputEvent->FindListObject("MultSelection"));
  if (MultSelection){
    // Percentile
    MultiplicityPercentile = (Double_t) MultSelection->GetMultiplicityPercentile("V0A"); // SPDTracklets V0A V0M (V0A for p-Pb)
  } else{
    //If this happens, re-check if AliMultSelectionTask ran before your task!
    AliInfo("AliMultSelection object not found!");
    return;
  }
  
  //  Fill generated events
  hGeneratedEvents->Fill(MultiplicityPercentile-0.01);
  
  //  Fill events with a given trigger
  // checking for event trigger
  Long64_t triggerMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t IsTrigINT7 = (triggerMask & AliVEvent::kINT7) ? true : false;
  
  if (IsTrigINT7) {
    hTriggeredEvents->Fill(MultiplicityPercentile-0.01);
  }
  
  //Event Cut (not needed but might be interessting to cross-check the event cut efficiency)
  Bool_t EventSelected = false;
  if (fEventCuts.AcceptEvent(fInputEvent)) {
    hSelectedEvents->Fill(MultiplicityPercentile-0.01);
    EventSelected = true;
  }
  
  for (Int_t iPartMC=0; iPartMC <= fAODArrayMCParticles->GetLast(); iPartMC++) {
    AliAODMCParticle* mcpart = (AliAODMCParticle*) fAODArrayMCParticles->At(iPartMC);
    if (!mcpart) {
      AliError("No MC particle found");
      continue;
    } else {
      // Check if nuclei
      if ( !(TMath::Abs(mcpart->GetPdgCode())==1000020030) ) continue; // He3 = 1000020030; 3H = 1000010030
      
      Double_t TruePt = mcpart->Pt();
      Double_t RapidityMC = mcpart->Y();
      // Correction to CMS rapidity
      RapidityMC = RapidityMC - 0.465;
      
      Double_t Source = -1;
      if (mcpart->IsPhysicalPrimary()){
        Source=0;
      } else if (mcpart->IsSecondaryFromMaterial()){
        Source=1;
      } else if (mcpart->IsSecondaryFromWeakDecay()){
        Source=2;
      }
      
      // Check if from weak decay of hypertriton
      if (mcpart->GetMother() > 0){
        AliAODMCParticle* mother = (AliAODMCParticle*) fAODArrayMCParticles->At(mcpart->GetMother());
        if (!mcpart) {
          AliError("No mother found");
        } else {
          while (TMath::Abs(mother->GetPdgCode()) == 1000020030) {
            if (mother->GetMother() > 0){
              mother = (AliAODMCParticle*) fAODArrayMCParticles->At(mother->GetMother());
              if (!mcpart) {
                AliError("No mother found");
                break;
              }
            } else{
              AliError("Mother with negative index found");
              break;
            }
          }
          if (TMath::Abs(mother->GetPdgCode()) == 1010010030){ // Hyper-triton
            Source = 2; // Correct Source
          }
        }
      }
      
      if (Source == 0) {
        hGeneratedSignal->Fill(MultiplicityPercentile-0.01,TruePt);
        if (IsTrigINT7) hGeneratedSignalEventTriggered->Fill(MultiplicityPercentile-0.01,TruePt);
        if (EventSelected) hGeneratedSignalEventSelected->Fill(MultiplicityPercentile-0.01,TruePt);
      }
    }
  }
  
  PostData(1,fOutput);
  PostData(2,fQAList);
}

//____________________________________________________________
void AliAnalysisTaskTriggerEfficiency::Terminate(Option_t *){
  // called at end of analysis
}
