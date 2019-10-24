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
// Task to calculate the trigger efficiency as a function of the multiplicity percentile
//
// Author:
//  Sebastian Hornung <Sebastian.Hornung@cern.ch>

#ifndef AliAnalysisTaskTriggerEfficiency_h
#define AliAnalysisTaskTriggerEfficiency_h

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TList;
class TClonesArray;
class TH1D;
class AliAODMCHeader;
class AliEventCuts;

class AliAnalysisTaskTriggerEfficiency: public AliAnalysisTaskSE {

public:
  // two class constructors
  AliAnalysisTaskTriggerEfficiency();
  AliAnalysisTaskTriggerEfficiency(const char *name);
  // class destructor
  virtual ~AliAnalysisTaskTriggerEfficiency();
  
  // called once at beginning or runtime
  virtual void UserCreateOutputObjects();
  // called for each event
  virtual void UserExec(Option_t* option);
  // called at end of analysis
  virtual void Terminate(Option_t* option);
  
  // Setter
  
private:
  // Variables
  Double_t MultiplicityPercentile;     // Multiplicity percentile
  
  // Results
  // Event number
  TH1D* hGeneratedEvents;              //!<! container for events as a function of multiplicity (without trigger)
  TH1D* hTriggeredEvents;              //!<! container for events as a function of multiplicity (with trigger)
  TH1D* hSelectedEvents;               //!<! container for events as a function of multiplicity (with trigger)

  // Generated signal (signal loss correction)
  TH2D* hGeneratedSignal;
  TH2D* hGeneratedSignalEventTriggered;
  TH2D* hGeneratedSignalEventSelected;
  
  // Objects
  TList *fOutput;                     //!<! Container for Task Output
  TList *fQAList;                     //!<! Container for QA Output
  AliAODMCHeader *fAODMCHeader;       //!<! MC info AOD
  TClonesArray *fAODArrayMCParticles; //!<! MC info particle AOD
  
  AliEventCuts fEventCuts;          // Event cuts
  
  AliAnalysisTaskTriggerEfficiency(const AliAnalysisTaskTriggerEfficiency&);
  AliAnalysisTaskTriggerEfficiency& operator=(const AliAnalysisTaskTriggerEfficiency&);

  ClassDef(AliAnalysisTaskTriggerEfficiency , 1);
};

#endif /* AliAnalysisTaskTriggerEfficiency_h */
