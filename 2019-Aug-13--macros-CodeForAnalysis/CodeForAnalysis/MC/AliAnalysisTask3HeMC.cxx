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
// Task for the analysis of the production of 3He as a function of multiplicity
//
// Author:
//  Sebastian Hornung <Sebastian.Hornung@cern.ch>

// Root header
#include <TChain.h>
#include <TList.h>

//AliRoot header
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"

#include "TH1F.h" // for dummy
#include "TH2D.h"
#include "TProfile2D.h"
#include "THnSparse.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"

//AliPhysics header
#include "AliAnalysisTask3HeMC.h"

//____________________________________________________________
AliAnalysisTask3HeMC::AliAnalysisTask3HeMC():
AliAnalysisTaskSE(),
kAODanalysis(kFALSE),
kHasMCdata(kFALSE),
kpp(kFALSE),
kpPb(kFALSE),
kPbPb(kFALSE),
fEtaMin(-0.9),
fEtaMax(0.9),
fMaxDCAxy(2.4),
fMaxDCAz(3.2),
fRejectKinkMother(kFALSE),
Multiplicity(-1),
MultiplicityPercentile(-1),
kTPCnSigmaCut(-5.0),
fHistMultEstZDep(NULL),
fHistTPCdEdxRigidity(NULL),
fHistTPCsigHe3(NULL),
fHistTPCsigHe4(NULL),
fHistTPCsigDeuteron(NULL),
fHistTPCsigTriton(NULL),
fHistTPCdEdxSigmaHe3(NULL),
fHistTPCdEdxSigmaHe4(NULL),
fHistTPCdEdxSigmaDeuteron(NULL),
fHistTPCdEdxSigmaTriton(NULL),
fHistPtTrueRecHe3(NULL),
fHistPtTrueRecHe4(NULL),
fHistPtTrueRecDeuteron(NULL),
fHistPtTrueRecTriton(NULL),
fHistTrueHe3(NULL),
fHistTrueHe4(NULL),
fHistTrueDeuteron(NULL),
fHistTrueTriton(NULL),
fHistTrueHyperTriton(NULL),
fEventCuts(),
fPIDResponse(NULL),
fUtils(NULL),
fOutput(NULL),
fQAList(NULL),
fAODMCHeader(NULL),
fAODArrayMCParticles(NULL)
{
}

//____________________________________________________________
AliAnalysisTask3HeMC::AliAnalysisTask3HeMC(const char *name):
AliAnalysisTaskSE(name),
kAODanalysis(kFALSE),
kHasMCdata(kFALSE),
kpp(kFALSE),
kpPb(kFALSE),
kPbPb(kFALSE),
fEtaMin(-0.9),
fEtaMax(0.9),
fMaxDCAxy(2.4),
fMaxDCAz(3.2),
fRejectKinkMother(kFALSE),
Multiplicity(-1),
MultiplicityPercentile(-1),
kTPCnSigmaCut(-5.0),
fHistMultEstZDep(NULL),
fHistMultValueVsPercientile(NULL),
fHistTPCdEdxRigidity(NULL),
fHistTPCsigHe3(NULL),
fHistTPCsigHe4(NULL),
fHistTPCsigDeuteron(NULL),
fHistTPCsigTriton(NULL),
fHistTPCdEdxSigmaHe3(NULL),
fHistTPCdEdxSigmaHe4(NULL),
fHistTPCdEdxSigmaDeuteron(NULL),
fHistTPCdEdxSigmaTriton(NULL),
fHistPtTrueRecHe3(NULL),
fHistPtTrueRecHe4(NULL),
fHistPtTrueRecDeuteron(NULL),
fHistPtTrueRecTriton(NULL),
fHistTrueHe3(NULL),
fHistTrueHe4(NULL),
fHistTrueDeuteron(NULL),
fHistTrueTriton(NULL),
fHistTrueHyperTriton(NULL),
fEventCuts(),
fPIDResponse(NULL),
fUtils(NULL),
fOutput(NULL),
fQAList(NULL),
fAODMCHeader(NULL),
fAODArrayMCParticles(NULL)
{

  fUtils = new AliAnalysisUtils();
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//____________________________________________________________
AliAnalysisTask3HeMC::~AliAnalysisTask3HeMC(){
  //
  // Destructor
  //

  // Delete output objects only if we are not running in PROOF mode because otherwise this produces a crash during merging
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(mgr && mgr->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
    if(fOutput) delete fOutput;
    if(fQAList) delete fQAList;
    if(fUtils) delete fUtils;
  }
}

//___________________________________________________
void AliAnalysisTask3HeMC::UserCreateOutputObjects(){ //Objects that are output are initialized in the function
  AliDebug(4, "Creating Output Objects");

  // Make lists for Output
  if(!fOutput) fOutput = new TList();
  fOutput->SetOwner();

  if(!fQAList) fQAList = new TList();
  fQAList->SetOwner();

  //   // Automatic determination of the analysis mode
  //   AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  //   if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
  //      SetAODAnalysis();
  //   } else {
  //      SetESDAnalysis();
  //      if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()) // possible for AOD too?
  //         SetHasMCData();
  //   }
  // Analysis mode
  AliInfo(Form("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD"));
  AliInfo(Form("MC Data available %s\n", HasMCData() ? "Yes" : "No"));

  fEventCuts.AddQAplotsToList(fQAList); /// Add event selection QA plots

  // QA histograms
  fHistMultEstZDep = new TH2D("fHistMultEstZDep","Multiplicity vs z_{vertex}", 100,0, 100, 201, -10, 10);
  fHistMultEstZDep->GetXaxis()->SetTitle("Multiplicity");
  fHistMultEstZDep->GetYaxis()->SetTitle("z_{Prim. vertex}");
  fHistMultEstZDep->Sumw2();
  fQAList->Add(fHistMultEstZDep);

  fHistMultValueVsPercientile = new TH2D("fHistMultValueVsPercientile","Multiplicity vs z_{vertex}",100, 0, 1000, 100, 0, 100);
  fHistMultValueVsPercientile->GetXaxis()->SetTitle("Multiplicity");
  fHistMultValueVsPercientile->GetYaxis()->SetTitle("Multiplicity Percentile");
  fHistMultValueVsPercientile->Sumw2();
  fQAList->Add(fHistMultValueVsPercientile);

  fHistTPCdEdxRigidity = new TH2D("fHistTPCdEdxRigidity","TPC dE/dx distribution",100,0.1,10,150, 0, 1500);
  fHistTPCdEdxRigidity->GetXaxis()->SetTitle("#it{p/z} (GeV/#it{c})");
  fHistTPCdEdxRigidity->GetYaxis()->SetTitle("TPC d#it{E}/d#it{x} (a.u.)");
  fHistTPCdEdxRigidity->Sumw2();
  fQAList->Add(fHistTPCdEdxRigidity);

  fHistTPCsigHe3 = new TProfile2D("fHistTPCsigHe3", "hist to check cut on TPC sigma for He3", 100, 0.1, 10.0, 150, 0.0, 1500);
  fHistTPCsigHe3->GetXaxis()->SetTitle("#it{p/z} (GeV/#it{c})");
  fHistTPCsigHe3->GetYaxis()->SetTitle("TPC d#it{E}/d#it{x} (a.u.)");
  fQAList->Add(fHistTPCsigHe3);

  fHistTPCsigDeuteron = new TProfile2D("fHistTPCsigDeuteron", "hist to check cut on TPC sigma for deuteron", 100, 0.1, 10.0, 150, 0.0, 1500);
  fHistTPCsigDeuteron->GetXaxis()->SetTitle("#it{p/z} (GeV/#it{c})");
  fHistTPCsigDeuteron->GetYaxis()->SetTitle("TPC d#it{E}/d#it{x} (a.u.)");
  fQAList->Add(fHistTPCsigDeuteron);

  fHistTPCsigTriton = new TProfile2D("fHistTPCsigTriton", "hist to check cut on TPC sigma for triton", 100, 0.1, 10.0, 150, 0.0, 1500);
  fHistTPCsigTriton->GetXaxis()->SetTitle("#it{p/z} (GeV/#it{c})");
  fHistTPCsigTriton->GetYaxis()->SetTitle("TPC d#it{E}/d#it{x} (a.u.)");
  fQAList->Add(fHistTPCsigTriton);

  //Result
  const Int_t nDim=11;
  const Int_t nBinMultiplicity=100;
  const Int_t nBinRigidity=100;
  const Int_t nBinTPCdEdx=80;
  const Int_t nBinPt=75;
  const Int_t nBinCharge=2;
  const Int_t nBinDCAxy=240;
  const Int_t nBinDCAz=240;
  const Int_t nBinRapidity=20;
  const Int_t nPtBins=15;
  const Int_t nITSBins=6;
  const Int_t nTPCclusterBins=20;

  //										mult  Rig	dE/dx		pT	Charge	DCAxy		DCAz	y	 source	ITShits	TPCcluster
  Double_t BinEdge_min[nDim] = {0,		0,		 -4,		0.5,	-2,	 -1.5,	-1.5,	-1.5,	-0.5,		0.5,		60};
  Double_t BinEdge_max[nDim] = {100,	10,		4,		8.0,	 2,		1.5,	 1.5,	 0.5,	 2.5,		6.5,		160};
  Int_t nBin[nDim] = {nBinMultiplicity,nBinRigidity,nBinTPCdEdx,nBinPt,nBinCharge, nBinDCAxy, nBinDCAz, nBinRapidity, 3,nITSBins,nTPCclusterBins};
  fHistTPCdEdxSigmaHe3 = new THnSparseD("fHistTPCdEdxSigmaHe3","TPC dE/dx distribution [sigma_{3He}]",nDim,nBin,BinEdge_min,BinEdge_max);
  fHistTPCdEdxSigmaHe3->GetAxis(0)->SetTitle("Multiplicity");
  fHistTPCdEdxSigmaHe3->GetAxis(1)->SetTitle("#it{p/z} (GeV/#it{c})");
  fHistTPCdEdxSigmaHe3->GetAxis(2)->SetTitle("TPC d#it{E}/d#it{x} (#sigma_{3He})");
  fHistTPCdEdxSigmaHe3->GetAxis(3)->SetTitle("#it{p}_{T}^{corr} (GeV/#it{c})");
  fHistTPCdEdxSigmaHe3->GetAxis(4)->SetTitle("(Anit-)Particle");
  fHistTPCdEdxSigmaHe3->GetAxis(5)->SetTitle("DCA in xy (cm)");
  fHistTPCdEdxSigmaHe3->GetAxis(6)->SetTitle("DCA in z (cm)");
  fHistTPCdEdxSigmaHe3->GetAxis(7)->SetTitle("Rapidity (lab frame)");
  fHistTPCdEdxSigmaHe3->GetAxis(8)->SetTitle("Source");
  fHistTPCdEdxSigmaHe3->GetAxis(9)->SetTitle("ITS hits");
  fHistTPCdEdxSigmaHe3->GetAxis(10)->SetTitle("TPC cluster");
  fHistTPCdEdxSigmaHe3->Sumw2();

  fHistTPCdEdxSigmaDeuteron = (THnSparseD*) fHistTPCdEdxSigmaHe3->Clone("fHistTPCdEdxSigmaDeuteron");
  fHistTPCdEdxSigmaDeuteron->SetTitle("TPC dE/dx distribution [sigma_{d}]");
  fHistTPCdEdxSigmaDeuteron->GetAxis(2)->SetTitle("TPC d#it{E}/d#it{x} (#sigma_{d})");

  //                             		 mult  Rig  	 dE/dx  	 pT    Charge   DCAxy   DCAz   y   source   ITShits   TPCcluster	TOFsigma
  Double_t BinEdgeTriton_min[nDim+1] = {0,      0,     -4,     0.5,   -2,    -1.5,   -1.5,   -1.5,   -0.5,      0.5, 		  60,			-5};
  Double_t BinEdgeTriton_max[nDim+1] = {100,   10,      4,      8.,   2,      1.5,   1.5,   0.5,      2.5,      6.5, 		  160,		 5};
  Int_t nBinTriton[nDim+1] = {nBinMultiplicity,nBinRigidity,nBinTPCdEdx,nBinPt,nBinCharge, nBinDCAxy, nBinDCAz, nBinRapidity, 3,nITSBins,nTPCclusterBins,40};

  fHistTPCdEdxSigmaTriton =new THnSparseD("fHistTPCdEdxSigmaTriton","TPC dE/dx distribution [sigma_{triton}]",nDim+1,nBinTriton,BinEdgeTriton_min,BinEdgeTriton_max);
  fHistTPCdEdxSigmaTriton->GetAxis(0)->SetTitle("Multiplicity");
  fHistTPCdEdxSigmaTriton->GetAxis(1)->SetTitle("#it{p/z} (GeV/#it{c})");
  fHistTPCdEdxSigmaTriton->GetAxis(2)->SetTitle("TPC d#it{E}/d#it{x} (#sigma_{triton})");
  fHistTPCdEdxSigmaTriton->GetAxis(3)->SetTitle("#it{p}_{T}^{corr} (GeV/#it{c})");
  fHistTPCdEdxSigmaTriton->GetAxis(4)->SetTitle("(Anit-)Particle");
  fHistTPCdEdxSigmaTriton->GetAxis(5)->SetTitle("DCA in xy (cm)");
  fHistTPCdEdxSigmaTriton->GetAxis(6)->SetTitle("DCA in z (cm)");
  fHistTPCdEdxSigmaTriton->GetAxis(7)->SetTitle("Rapidity (lab frame)");
  fHistTPCdEdxSigmaTriton->GetAxis(8)->SetTitle("Source");
  fHistTPCdEdxSigmaTriton->GetAxis(9)->SetTitle("ITS hits");
  fHistTPCdEdxSigmaTriton->GetAxis(10)->SetTitle("TPC cluster");
  fHistTPCdEdxSigmaTriton->GetAxis(11)->SetTitle("TOF t (#sigma_{triton})");
  fHistTPCdEdxSigmaTriton->Sumw2();

  fHistTPCdEdxSigmaHe4 = (THnSparseD*) fHistTPCdEdxSigmaTriton->Clone("fHistTPCdEdxSigmaHe4");
  fHistTPCdEdxSigmaHe4->SetTitle("TPC dE/dx distribution [sigma_{^{4}He}]");
  fHistTPCdEdxSigmaHe4->GetAxis(2)->SetTitle("TPC d#it{E}/d#it{x} (#sigma_{^{4}He})");
  fHistTPCdEdxSigmaHe4->GetAxis(11)->SetTitle("TOF t (#sigma_{^{4}He})");
  fHistTPCdEdxSigmaHe4->Sumw2();

  fOutput->Add(fHistTPCdEdxSigmaHe3);
  fOutput->Add(fHistTPCdEdxSigmaHe4);
  fOutput->Add(fHistTPCdEdxSigmaDeuteron);
  fOutput->Add(fHistTPCdEdxSigmaTriton);

  if(HasMCData()){
    fHistPtTrueRecHe3 = new TH2D("fHistPtTrueRecHe3", "(#it{p}_{T}^{true} - #it{p}_{T}^{rec}) vs #it{p}_{T}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
    fHistPtTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{T}^{rec} (GeV/#it{c})");
    fHistPtTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{T}^{true}- #it{p}_{T}^{rec}) (GeV/#it{c})");
    fOutput->Add(fHistPtTrueRecHe3);

    fHistPtTrueRecHe4 = (TH2D*) fHistPtTrueRecHe3->Clone("fHistPtTrueRecHe4");
    fHistPtTrueRecHe4->SetName("fHistPtTrueRecHe4");
    fOutput->Add(fHistPtTrueRecHe4);

    fHistPtTrueRecDeuteron = (TH2D*) fHistPtTrueRecHe3->Clone("fHistPtTrueRecDeuteron");
    fHistPtTrueRecDeuteron->SetName("fHistPtTrueRecDeuteron");
    fOutput->Add(fHistPtTrueRecDeuteron);

    fHistPtTrueRecTriton = (TH2D*) fHistPtTrueRecHe3->Clone("fHistPtTrueRecTriton");
    fHistPtTrueRecTriton->SetName("fHistPtTrueRecTriton");
    fOutput->Add(fHistPtTrueRecTriton);

    const Int_t nDimMC=5;
    //                           mult        pTtrue    Source        Charge      y
    Double_t BinEdgeMC_min[nDim] = {0, 			0.5,	  	-0.5, 	   	-2,	      -1.5};
    Double_t BinEdgeMC_max[nDim] = {100,		8.0,	  	2.5,   	  		2,		  0.5};
    Int_t nBinMC[nDimMC] = {nBinMultiplicity,nBinPt,3,nBinCharge, nBinRapidity};

    fHistTrueHe3 = new THnSparseD("fHistTrueHe3", "MC truth He3", nDimMC, nBinMC, BinEdgeMC_min, BinEdgeMC_max);
    fHistTrueHe3->GetAxis(0)->SetTitle("Multiplicity");
    fHistTrueHe3->GetAxis(1)->SetTitle("#it{p}_{T}^{corr} (GeV/#it{c})");
    fHistTrueHe3->GetAxis(2)->SetTitle("Source");
    fHistTrueHe3->GetAxis(3)->SetTitle("(Anit-)Particle");
    fHistTrueHe3->GetAxis(4)->SetTitle("Rapidity (lab frame)");

    fHistTrueHe4 = (THnSparseD*) fHistTrueHe3->Clone("fHistTrueHe4");
    fHistTrueHe4->SetTitle("MC truth He4");

    fHistTrueDeuteron = (THnSparseD*) fHistTrueHe3->Clone("fHistTrueDeuteron");
    fHistTrueDeuteron->SetTitle("MC truth deuteron");

    fHistTrueTriton = (THnSparseD*) fHistTrueHe3->Clone("fHistTrueTriton");
    fHistTrueTriton->SetTitle("MC truth triton");

    fHistTrueHyperTriton = (THnSparseD*) fHistTrueHe3->Clone("fHistTrueHyperTriton");
    fHistTrueHyperTriton->SetTitle("MC truth hyper-triton");

    fOutput->Add(fHistTrueHe3);
    fOutput->Add(fHistTrueHe4);
    fOutput->Add(fHistTrueDeuteron);
    fOutput->Add(fHistTrueTriton);
    fOutput->Add(fHistTrueHyperTriton);
  }

  PostData(1,fOutput);
  PostData(2,fQAList);
}

//___________________________________________________
void AliAnalysisTask3HeMC::UserExec(Option_t *){ //called for each event

  AliDebug(4, "Starting Single Event Analysis");
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return;
  }

  if(IsESDanalysis() && HasMCData()){
    AliDebug(3, Form("MC Event: %p", fMCEvent));
    if(!fMCEvent){
      AliError("No MC Event, but MC Data required");
      return;
    }
    // Protect against missing MC trees
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!mcH){
      AliError("No MC Event Handler available");
      return;
    }
    if(!mcH->InitOk()) return;
    if(!mcH->TreeK()) return;
    if(!mcH->TreeTR()) return;
  }

  if(IsAODanalysis() && HasMCData()){
    //       take MC info
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
  }

  //PID response
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse) {
    AliError("No PID Response found");
    return;
  }

  //Event Cut
  if (!fEventCuts.AcceptEvent(fInputEvent)) {
    PostData(2, fQAList);
    return;
  }

  // Multiplicity / centrality framework
  AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
  if(!MultSelection->IsEventSelected()) return;
  if (MultSelection){
    AliMultEstimator* MultiplicityEstimator = NULL;

    // Percentile
    MultiplicityPercentile = (Double_t) MultSelection->GetMultiplicityPercentile("V0A"); // SPDTracklets V0A V0M (V0A for p-Pb)
    // Value
    MultiplicityEstimator = MultSelection->GetEstimator("V0A"); // Multiplicity estimater based on SPD tracklets
    Multiplicity = (Double_t) MultiplicityEstimator->GetValue(); // returns multiplicity only for events used in calibration

    fHistMultValueVsPercientile->Fill(Multiplicity, MultiplicityPercentile);
  }else{
    //If this happens, re-check if AliMultSelectionTask ran before your task!
    AliInfo("AliMultSelection object not found!");
    return;
  }
  Double_t zVertex = fInputEvent->GetPrimaryVertex()->GetZ();
  fHistMultEstZDep->Fill(MultiplicityPercentile, zVertex);

  // Run analysis for AOD
  if (IsAODanalysis()) {
    ProcessAOD();
  }

  PostData(1,fOutput);
  PostData(2,fQAList);
}

//____________________________________________________________
void AliAnalysisTask3HeMC::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //
}

//___________________________________________________
void AliAnalysisTask3HeMC::ProcessAOD(){ // Loop over tracks and perform analysis
  AliDebug(3, "Processing AOD Event");

  AliAODEvent *fAOD = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!fAOD){
    AliError("AOD Event required for AOD Analysis");
    return;
  }

  Int_t fNumberOfVertices = fAOD->GetNumberOfVertices();
  Double_t fListOfmotherkink[fNumberOfVertices];
  Int_t fNumberOfMotherkink = 0;
  for(Int_t ivertex=0; ivertex < fNumberOfVertices; ivertex++){
    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()==AliAODVertex::kKink){
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      fListOfmotherkink[fNumberOfMotherkink] = idmother;
      fNumberOfMotherkink++;
    }
  }

  //loop over MC particles
  if (HasMCData()) {
    for (Int_t iPartMC=0; iPartMC <= fAODArrayMCParticles->GetLast(); iPartMC++) {
      AliAODMCParticle* mcpart = (AliAODMCParticle*) fAODArrayMCParticles->At(iPartMC);
      if (!mcpart) {
        AliError("No MC particle found");
        continue;
      } else {
        Int_t ParticleType = 0;
        if (TMath::Abs(mcpart->GetPdgCode())== 1000020030) ParticleType = 1; //He3
        if (TMath::Abs(mcpart->GetPdgCode())== 1000020040) ParticleType = 2; //He4
        if (TMath::Abs(mcpart->GetPdgCode())== 1000010020) ParticleType = 3; //deuteron
        if (TMath::Abs(mcpart->GetPdgCode())== 1000010030) ParticleType = 4; //triton
        if (TMath::Abs(mcpart->GetPdgCode())== 1010010030) ParticleType = 5; //hyper-triton
        if (ParticleType == 0) continue;

        Double_t TruePt = mcpart->Pt();
        Double_t Charge = mcpart->GetPdgCode()/TMath::Abs(mcpart->GetPdgCode()); // mcpart->Charge() sometimes gives back the wrong charge (±99 instead of ±3)
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

        // Check if 3He is from weak decay of hypertriton
        if (ParticleType == 1) {
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
                Source = 2; // Correct source
              }
            }
          }
        }

        // Check if triton is from weak decay of hypertriton
        if (ParticleType == 4) {
          if (mcpart->GetMother() > 0){
            AliAODMCParticle* mother = (AliAODMCParticle*) fAODArrayMCParticles->At(mcpart->GetMother());
            if (!mcpart) {
              AliError("No mother found");
            } else {
              while (TMath::Abs(mother->GetPdgCode()) == 1000010030) {
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
                Source = 2; // Correct source
              }
            }
          }
        }

        Double_t FillSparseMC[] = {MultiplicityPercentile,TruePt,Source,Charge/2.,RapidityMC};
        if (ParticleType == 1) {
          fHistTrueHe3->Fill(FillSparseMC);
        } else if (ParticleType == 2){
          fHistTrueHe4->Fill(FillSparseMC);
        } else if (ParticleType == 3){
          fHistTrueDeuteron->Fill(FillSparseMC);
        } else if (ParticleType == 4){
          fHistTrueTriton->Fill(FillSparseMC);
        } else if (ParticleType == 5){
          fHistTrueHyperTriton->Fill(FillSparseMC);
        }
      }
    }
  }

  Int_t nTracks(fAOD->GetNumberOfTracks());
  for(Int_t iTrack=0; iTrack < nTracks; iTrack++) { // loop over all the tracks

    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track) continue;

    Int_t ParticleType = 0; // needed since He3 and He4 are not separated in the TPC for high pT
    Double_t TruePt = -1;
    Double_t Source = -1;

    Bool_t HyperTritonFeedDown = 0;

    if (HasMCData()) {
      AliAODMCParticle* mcpart = (AliAODMCParticle*) fAODArrayMCParticles->At(TMath::Abs(track->GetLabel()));
      if (!mcpart) {
        AliError("No MC particle found");
        continue;
      } else {
        if (TMath::Abs(mcpart->GetPdgCode())== 1000020030) ParticleType = 1; //He3
        if (TMath::Abs(mcpart->GetPdgCode())== 1000020040) ParticleType = 2; //He4
        if (TMath::Abs(mcpart->GetPdgCode())== 1000010020) ParticleType = 3; //deuteron
        if (TMath::Abs(mcpart->GetPdgCode())== 1000010030) ParticleType = 4; //triton
        if (ParticleType == 0) continue;

        TruePt = mcpart->Pt();

        if (mcpart->IsPhysicalPrimary()){
          Source=0;
        } else if (mcpart->IsSecondaryFromMaterial()){
          Source=1;
        } else if (mcpart->IsSecondaryFromWeakDecay()){
          Source=2;
        }

        // Check if 3He is from weak decay of hypertriton
        if (ParticleType == 1) {
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
                HyperTritonFeedDown = kTRUE;
                Source = 2; //Correct source of 3He to weak decay -> Can be used to calculate the efficiency of having a 3He from hypertrition
              }
            }
          }
        }

        // Check if triton is from weak decay of hypertriton
        if (ParticleType == 4) {
          if (mcpart->GetMother() > 0){
            AliAODMCParticle* mother = (AliAODMCParticle*) fAODArrayMCParticles->At(mcpart->GetMother());
            if (!mcpart) {
              AliError("No mother found");
            } else {
              while (TMath::Abs(mother->GetPdgCode()) == 1000010030) {
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
                Source = 2; // Correct source
              }
            }
          }
        }

      }
    }

    // Geomatrical acceptance
    if((track->Eta() < fEtaMin) || (track->Eta() > fEtaMax)) continue;

    // check AOD filter bit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;

    // Reject kink mothers
    if(fRejectKinkMother){
      Bool_t kinkmotherpass = kTRUE;
      for(Int_t kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++){
        if(track->GetID() == fListOfmotherkink[kinkmother]){
          kinkmotherpass = kFALSE;
          continue;
        }
      }
      if(!kinkmotherpass) continue;
    }

    Double_t DCAxy = GetDCAxy(track);
    Double_t DCAz = GetDCAz(track);
    // additional track cuts
    if(TMath::Abs(DCAxy) > fMaxDCAxy) continue;
    if(TMath::Abs(DCAz) > fMaxDCAz) continue;

    if (track->GetTPCNcls() < 60) continue;
    if(track->GetITSNcls() < 1) continue;

    Double_t ITShits = (Double_t) track->GetITSNcls();
    Double_t TPCcluster = (Double_t) track->GetTPCNcls()+0.001; // cluster value shifted to have bins with >= x

    // Get particle properties
    Double_t TrackPt = track->Pt(); // particle with charge 1
    Double_t Charge = (Double_t) track->Charge(); // not the physical charge is always -1 or 1; more like sign of charge
    Double_t Rigidity = TMath::Abs(track->P());
    Double_t dEdx = track->GetTPCsignal();

    // Fill TPC QA histogram (performance)
    fHistTPCdEdxRigidity->Fill(Rigidity,dEdx);

    //TPC nSigma
    Double_t dEdxSigmaHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
    Double_t dEdxSigmaHe4 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kAlpha);
    Double_t dEdxSigmaDeuteron = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
    Double_t dEdxSigmaTriton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);

    Double_t TOFSigmaTriton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kTriton);
    Double_t TOFSigma4He =  fPIDResponse->NumberOfSigmasTOF(track, AliPID::kAlpha);

    // basic PID cut
    Double_t Rapidity = track->Y() - 0.465;
    Double_t FillSparse[] = {MultiplicityPercentile,Rigidity,dEdxSigmaDeuteron,TrackPt,Charge,DCAxy,DCAz,Rapidity,Source, ITShits, TPCcluster};

    //Fill deuteron candidates
    if (ParticleType==3){
      if (track->GetTPCNcls() >= 70 && track->GetITSNcls() > 1 && HasMCData() && TMath::Abs(dEdxSigmaDeuteron) < 3){
        if (Rapidity < 0 && Rapidity > -1. ) {
          fHistPtTrueRecDeuteron->Fill(TrackPt, (TruePt-TrackPt));
        }
      }
      fHistTPCsigHe3->Fill(Rigidity, dEdx, dEdxSigmaDeuteron);
      fHistTPCdEdxSigmaDeuteron->Fill(FillSparse);
    }

    //Fill triton candidates
    if (ParticleType==4){
      if (track->GetTPCNcls() >= 70 && track->GetITSNcls() > 1 && HasMCData() && TMath::Abs(dEdxSigmaTriton) < 3 ){ // && TMath::Abs(TOFSigmaTriton) < 3){
        if (Rapidity < 0 && Rapidity > -1. ) {
          fHistPtTrueRecTriton->Fill(TrackPt, (TruePt-TrackPt));
        }
      }
      fHistTPCsigHe3->Fill(Rigidity, dEdx, dEdxSigmaTriton);
      // Fill container for efficiency determination
      Double_t FillSparseTriton[] = {MultiplicityPercentile,Rigidity,dEdxSigmaTriton,TrackPt,Charge,DCAxy,DCAz,Rapidity,Source, ITShits, TPCcluster, TOFSigmaTriton};
      fHistTPCdEdxSigmaTriton->Fill(FillSparseTriton);
    }


    //below only particles with charge 2
    TrackPt = 2*track->Pt();
    FillSparse[3]= 2*track->Pt();

    if (dEdxSigmaTriton > kTPCnSigmaCut){
      Rapidity = CalculateRapidity(track->Px(), track->Py(), track->Pz(), "3He");
      //Fill 3He candidates
      if (ParticleType==1){
        if (track->GetTPCNcls() >= 70 && track->GetITSNcls() > 1 && HasMCData() && TMath::Abs(dEdxSigmaHe3) < 3){
          if (Rapidity < 0 && Rapidity > -1. ) {
            fHistPtTrueRecHe3->Fill(TrackPt, (TruePt-TrackPt));
          }
        }

        fHistTPCsigHe3->Fill(Rigidity, dEdx, dEdxSigmaHe3);
        FillSparse[2] = dEdxSigmaHe3;
        FillSparse[7] = CalculateRapidity(track->Px(), track->Py(), track->Pz(), "3He");
        FillSparse[3] = 2*track->Pt() + 0.00068 - 0.3641*TMath::Exp(-0.2386*TMath::Power(2*track->Pt(),3.0));
        // Fill container for efficiency determination
        fHistTPCdEdxSigmaHe3->Fill(FillSparse);
      }

      //Fill 4He candidates
      if (ParticleType==2){
        if (track->GetTPCNcls() >= 70 && track->GetITSNcls() > 1 && HasMCData() && TMath::Abs(dEdxSigmaHe4) < 3){
          if (Rapidity < 0 && Rapidity > -1. ) {
            fHistPtTrueRecHe4->Fill(TrackPt, (TruePt-TrackPt));
          }
        }

        TrackPt = 2*track->Pt() - 0.00295074 + 0.332899*TMath::Exp(-0.815359*2*track->Pt());
        Rapidity = CalculateRapidity(track->Px(), track->Py(), track->Pz(), "4He");
        // Fill container for efficiency determination
        Double_t FillSparse4He[] = {MultiplicityPercentile,Rigidity,dEdxSigmaHe4,TrackPt,Charge,DCAxy,DCAz,Rapidity,Source, ITShits, TPCcluster, TOFSigma4He};

        fHistTPCdEdxSigmaHe4->Fill(FillSparse4He);
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTask3HeMC::BinLogAxis(TAxis *axis) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];

  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
}
//________________________________________________________________________
Double_t AliAnalysisTask3HeMC::GetDCAxy (AliAODTrack *track)  {

  AliAODEvent *fAODevent = dynamic_cast<AliAODEvent *>(fInputEvent);

  Double_t impactParameter[2];
  Double_t covarianceMatrix[3];
  if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

  Double_t DCAxy = impactParameter[0];

  return DCAxy;
}
//________________________________________________________________________
Double_t AliAnalysisTask3HeMC::GetDCAz (AliAODTrack *track)  {

  AliAODEvent *fAODevent = dynamic_cast<AliAODEvent *>(fInputEvent);
  Double_t impactParameter[2];
  Double_t covarianceMatrix[3];
  if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

  Double_t DCAz = impactParameter[1];

  return DCAz;
}
//________________________________________________________________________
Double_t AliAnalysisTask3HeMC::CalculateRapidity(Double_t px, Double_t py, Double_t pz, TString type){
  // calculation of the rapidity
  Double_t p = TMath::Sqrt(px*px + py*py+ pz*pz);
  Double_t energy = 0;
  if (type == "4He") {
    energy = TMath::Sqrt(p*2*p*2 + AliPID::ParticleMass(AliPID::kAlpha)*AliPID::ParticleMass(AliPID::kAlpha));
  } else {
    energy = TMath::Sqrt(p*2*p*2 + AliPID::ParticleMass(AliPID::kHe3)*AliPID::ParticleMass(AliPID::kHe3));
  }
  Double_t rap = 0.5*TMath::Log((energy + pz*2)/(energy - pz*2));
  // Correction to CMS rapidity
  rap = rap - 0.465;

  return rap;
}
