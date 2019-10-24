#include "TritonAnalysis.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TList.h"
#include "stdio.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFractionFitter.h"

Double_t legsize = 0.04;
Double_t labelsize = 0.07;
TString pTLabel= "#it{p}_{T} (GeV/#it{c})";
TString SpectraLabel= "#frac{d^{2}N}{d#it{y} d#it{p}_{T}} (GeV/#it{c})";
TString RatioLabel = "^{3}#bar{H} / ^{3}H";
Double_t PlotMinPt = 1;
Double_t PlotMaxPt = 3.5;

Double_t DCAcutXY;
Double_t DCAcutZ;
Double_t minRapidity;
Double_t maxRapidity;
Int_t ITShits;
Int_t TPCcluster;
Double_t TPCselectMin; // TPC selection window lower boundary
Double_t TPCselectMax; // TPC selection window upper boundary
Double_t TOFselect; // TOF selection window

Bool_t SecondaryGaus = kTRUE;

// Multiplicity
Bool_t MultiDep = false; // Switch for multiplicity dependent analysis (true = multiplicity depenent , false = multiplicity integrated)
Double_t MultiplicityPercentile[2] = {0,100}; // Set the boundaries of the multiplicity percentile of interest

// Binning
const Int_t nPtBins=4;
Double_t Binning[]={1., 1.5, 2., 2.5, 3.0};

//___________________________________________________________________________________________
TritonAnalysis::TritonAnalysis():
fInputFileName(NULL),
fInputListName(NULL),
fOutputFileName(NULL),
tree(NULL),
SaveToFile(0),
magFieldSign(0),
multPercentile_V0M(0),
multPercentile_V0A(0),
multPercentile_V0C(0),
multPercentile_OnlineV0M(0),
multPercentile_OnlineV0A(0),
multPercentile_OnlineV0C(0),
multPercentile_ADM(0),
multPercentile_ADA(0),
multPercentile_ADC(0),
multPercentile_SPDClusters(0),
multPercentile_SPDTracklets(0),
multPercentile_RefMult05(0),
multPercentile_RefMult08(0),
multPercentile_CL1(0),
multPercentile_ZNA(0),
Ntrk_V0M(0),
Ntrk_V0A(0),
Ntrk_V0C(0),
Ntrk_OnlineV0M(0),
Ntrk_OnlineV0A(0),
Ntrk_OnlineV0C(0),
Ntrk_ADM(0),
Ntrk_ADA(0),
Ntrk_ADC(0),
Ntrk_SPDClusters(0),
Ntrk_SPDTracklets(0),
Ntrk_RefMult05(0),
Ntrk_RefMult08(0),
nVertexContributors(0),
xVertex(0),
yVertex(0),
zVertex(0),
px(0),
py(0),
pz(0),
TPCmomentum(0),
TRDmomentum(0),
integratedLength(0),
timeOfFlight(0),
beta(0),
gamma(0),
mass(0),
trackID(0),
eta(0),
phi(0),
theta(0),
y(0),
q(0),
dcaxy(0),
dcaz(0),
nTPC_Clusters(0),
nTRD_Clusters(0),
nITS_Clusters(0),
nTPC_FindableClusters(0),
nTPC_CrossedRows(0),
nTPC_Clusters_dEdx(0),
HasPointOnITSLayer0(0),
HasPointOnITSLayer1(0),
HasPointOnITSLayer2(0),
HasPointOnITSLayer3(0),
HasPointOnITSLayer4(0),
HasPointOnITSLayer5(0),
HasSharedPointOnITSLayer0(0),
HasSharedPointOnITSLayer1(0),
HasSharedPointOnITSLayer2(0),
HasSharedPointOnITSLayer3(0),
HasSharedPointOnITSLayer4(0),
HasSharedPointOnITSLayer5(0),
chi2_TPC(0),
chi2_NDF(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
TRDsignal(0),
HMPIDsignal(0),
nSigmaITS_He4(0),
nSigmaTPC_He4(0),
nSigmaTOF_He4(0),
nSigmaTRD_He4(0),
nSigmaHMPID_He4(0),
nSigmaITS_He3(0),
nSigmaTPC_He3(0),
nSigmaTOF_He3(0),
nSigmaTRD_He3(0),
nSigmaHMPID_He3(0),
nSigmaITS_Trit(0),
nSigmaTPC_Trit(0),
nSigmaTOF_Trit(0),
nSigmaTRD_Trit(0),
nSigmaHMPID_Trit(0),
nSigmaITS_Deut(0),
nSigmaTPC_Deut(0),
nSigmaTOF_Deut(0),
nSigmaTRD_Deut(0),
nSigmaHMPID_Deut(0),
nSigmaITS_Prot(0),
nSigmaTPC_Prot(0),
nSigmaTOF_Prot(0),
nSigmaTRD_Prot(0),
nSigmaHMPID_Prot(0),
nSigmaITS_Pion(0),
nSigmaTPC_Pion(0),
nSigmaTOF_Pion(0),
nSigmaTRD_Pion(0),
nSigmaHMPID_Pion(0),
nSigmaITS_Kaon(0),
nSigmaTPC_Kaon(0),
nSigmaTOF_Kaon(0),
nSigmaTRD_Kaon(0),
nSigmaHMPID_Kaon(0),
nSigmaITS_Elec(0),
nSigmaTPC_Elec(0),
nSigmaTOF_Elec(0),
nSigmaTRD_Elec(0),
nSigmaHMPID_Elec(0),
PtCorr(NULL),
NumberOfEvents(0),
fHistoRapidity(NULL),
fHistoMean(NULL),
fHistoWidth(NULL),
fHistoPtSpectrum(NULL),
fHistoAntiPtSpectrum(NULL),
fHistoRatioSpectra(NULL)
{
  gStyle->SetOptStat(false);         // switch off stat box
  gStyle->SetOptFit(kTRUE);         // switch on fit box
}

//___________________________________________________________________________________________
TritonAnalysis::~TritonAnalysis()
{
  delete tree;
  delete PtCorr;
  delete fHistoRapidity;
  delete fHistoMean;
  delete fHistoWidth;
  delete fHistoPtSpectrum;
  delete fHistoAntiPtSpectrum;
  delete fHistoRatioSpectra;
}

//___________________________________________________________________________________________
void TritonAnalysis::GetInputTree()  {

  TFile *inputFile = TFile::Open (Form("%s",fInputFileName));
  TList *inputList = (TList*) inputFile->Get(fInputListName);

  tree = (TTree*) inputFile->GetObjectChecked("Trees/reducedTree_Triton","TTree"); // tree is stored in a subfolder called Trees

  // Get Number of events
  TH1F* EventSelection = (TH1F*) inputList->FindObject("histoEventSelection");
  NumberOfEvents = EventSelection->GetBinContent(3);
  printf("Number of events: %d \n", NumberOfEvents);

  if (MultiDep) {
    // Get Number of events for multiplicity depenedent case
    TH1F* histoEventMultiplicity = (TH1F*) inputList->FindObject("histoEventMultiplicity");
    NumberOfEvents = histoEventMultiplicity->Integral(histoEventMultiplicity->FindBin(MultiplicityPercentile[0]+0.01),histoEventMultiplicity->FindBin(MultiplicityPercentile[1]-0.01));
    printf("Number of events: %d \n", NumberOfEvents);
  }
}

//___________________________________________________________________________________________
void TritonAnalysis::SetTreeBranches()  {

  tree -> SetBranchAddress("magFieldSign",&magFieldSign);
  tree -> SetBranchAddress("multPercentile_V0M",&multPercentile_V0M);
  tree -> SetBranchAddress("multPercentile_V0A",&multPercentile_V0A);
  tree -> SetBranchAddress("multPercentile_V0C",&multPercentile_V0C);
  //   tree -> SetBranchAddress("multPercentile_OnlineV0M",&multPercentile_OnlineV0M);
  //   tree -> SetBranchAddress("multPercentile_OnlineV0A",&multPercentile_OnlineV0A);
  //   tree -> SetBranchAddress("multPercentile_OnlineV0C",&multPercentile_OnlineV0C);
  //   tree -> SetBranchAddress("multPercentile_ADM",&multPercentile_ADM);
  //   tree -> SetBranchAddress("multPercentile_ADA",&multPercentile_ADA);
  //   tree -> SetBranchAddress("multPercentile_ADC",&multPercentile_ADC);
  //   tree -> SetBranchAddress("multPercentile_SPDClusters",&multPercentile_SPDClusters);
  tree -> SetBranchAddress("multPercentile_SPDTracklets",&multPercentile_SPDTracklets);
  //   tree -> SetBranchAddress("multPercentile_RefMult05",&multPercentile_RefMult05);
  //   tree -> SetBranchAddress("multPercentile_RefMult08",&multPercentile_RefMult08);
  tree -> SetBranchAddress("multPercentile_CL1",&multPercentile_CL1);
  tree -> SetBranchAddress("multPercentile_ZNA",&multPercentile_ZNA);
  //   tree -> SetBranchAddress("Ntrk_OnlineV0M",&Ntrk_OnlineV0M);
  //   tree -> SetBranchAddress("Ntrk_OnlineV0A",&Ntrk_OnlineV0A);
  //   tree -> SetBranchAddress("Ntrk_OnlineV0C",&Ntrk_OnlineV0C);
  tree -> SetBranchAddress("Ntrk_V0M",&Ntrk_V0M);
  tree -> SetBranchAddress("Ntrk_V0A",&Ntrk_V0A);
  tree -> SetBranchAddress("Ntrk_V0C",&Ntrk_V0C);
  //   tree -> SetBranchAddress("Ntrk_ADM",&Ntrk_ADM);
  //   tree -> SetBranchAddress("Ntrk_ADA",&Ntrk_ADA);
  //   tree -> SetBranchAddress("Ntrk_ADC",&Ntrk_ADC);
  //   tree -> SetBranchAddress("Ntrk_SPDClusters",&Ntrk_SPDClusters);
  tree -> SetBranchAddress("Ntrk_SPDTracklets",&Ntrk_SPDTracklets);
  //   tree -> SetBranchAddress("Ntrk_RefMult05",&Ntrk_RefMult05);
  //   tree -> SetBranchAddress("Ntrk_RefMult08",&Ntrk_RefMult08);
  //   tree -> SetBranchAddress("nVertexContributors",&nVertexContributors);
  tree -> SetBranchAddress("xVertex",&xVertex);
  tree -> SetBranchAddress("yVertex",&yVertex);
  tree -> SetBranchAddress("zVertex",&zVertex);
  tree -> SetBranchAddress("px",&px);
  tree -> SetBranchAddress("py",&py);
  tree -> SetBranchAddress("pz",&pz);
  tree -> SetBranchAddress("TPCmomentum",&TPCmomentum);
  tree -> SetBranchAddress("TRDmomentum",&TRDmomentum);
  tree -> SetBranchAddress("integratedLength",&integratedLength);
  tree -> SetBranchAddress("timeOfFlight",&timeOfFlight);
  tree -> SetBranchAddress("beta",&beta);
  tree -> SetBranchAddress("gamma",&gamma);
  tree -> SetBranchAddress("mass",&mass);
  tree -> SetBranchAddress("trackID",&trackID);
  tree -> SetBranchAddress("eta",&eta);
  tree -> SetBranchAddress("phi",&phi);
  tree -> SetBranchAddress("theta",&theta);
  tree -> SetBranchAddress("y",&y);
  tree -> SetBranchAddress("q",&q);
  tree -> SetBranchAddress("dcaxy",&dcaxy);
  tree -> SetBranchAddress("dcaz",&dcaz);
  tree -> SetBranchAddress("nTPC_Clusters",&nTPC_Clusters);
  tree -> SetBranchAddress("nTRD_Clusters",&nTRD_Clusters);
  tree -> SetBranchAddress("nITS_Clusters",&nITS_Clusters);
  tree -> SetBranchAddress("nTPC_FindableClusters",&nTPC_FindableClusters);
  tree -> SetBranchAddress("nTPC_CrossedRows",&nTPC_CrossedRows);
  tree -> SetBranchAddress("nTPC_Clusters_dEdx",&nTPC_Clusters_dEdx);
  tree -> SetBranchAddress("HasPointOnITSLayer0",&HasPointOnITSLayer0);
  tree -> SetBranchAddress("HasPointOnITSLayer1",&HasPointOnITSLayer1);
  tree -> SetBranchAddress("HasPointOnITSLayer2",&HasPointOnITSLayer2);
  tree -> SetBranchAddress("HasPointOnITSLayer3",&HasPointOnITSLayer3);
  tree -> SetBranchAddress("HasPointOnITSLayer4",&HasPointOnITSLayer4);
  tree -> SetBranchAddress("HasPointOnITSLayer5",&HasPointOnITSLayer5);
  tree -> SetBranchAddress("HasSharedPointOnITSLayer0",&HasSharedPointOnITSLayer0);
  tree -> SetBranchAddress("HasSharedPointOnITSLayer1",&HasSharedPointOnITSLayer1);
  tree -> SetBranchAddress("HasSharedPointOnITSLayer2",&HasSharedPointOnITSLayer2);
  tree -> SetBranchAddress("HasSharedPointOnITSLayer3",&HasSharedPointOnITSLayer3);
  tree -> SetBranchAddress("HasSharedPointOnITSLayer4",&HasSharedPointOnITSLayer4);
  tree -> SetBranchAddress("HasSharedPointOnITSLayer5",&HasSharedPointOnITSLayer5);
  tree -> SetBranchAddress("chi2_TPC",&chi2_TPC);
  tree -> SetBranchAddress("chi2_NDF",&chi2_NDF);
  tree -> SetBranchAddress("chi2_ITS",&chi2_ITS);
  tree -> SetBranchAddress("ITSsignal",&ITSsignal);
  tree -> SetBranchAddress("TPCsignal",&TPCsignal);
  tree -> SetBranchAddress("TOFsignal",&TOFsignal);
  tree -> SetBranchAddress("TRDsignal",&TRDsignal);
  tree -> SetBranchAddress("HMPIDsignal",&HMPIDsignal);
  //   tree -> SetBranchAddress("nSigmaITS_He4",&nSigmaITS_He4);
  //   tree -> SetBranchAddress("nSigmaTPC_He4",&nSigmaTPC_He4);
  //   tree -> SetBranchAddress("nSigmaTOF_He4",&nSigmaTOF_He4);
  //   tree -> SetBranchAddress("nSigmaTRD_He4",&nSigmaTRD_He4);
  //   tree -> SetBranchAddress("nSigmaHMPID_He4",&nSigmaHMPID_He4);
  tree -> SetBranchAddress("nSigmaITS_He3",&nSigmaITS_He3);
  tree -> SetBranchAddress("nSigmaTPC_He3",&nSigmaTPC_He3);
  tree -> SetBranchAddress("nSigmaTOF_He3",&nSigmaTOF_He3);
  tree -> SetBranchAddress("nSigmaTRD_He3",&nSigmaTRD_He3);
  tree -> SetBranchAddress("nSigmaHMPID_He3",&nSigmaHMPID_He3);
  tree -> SetBranchAddress("nSigmaITS_Trit",&nSigmaITS_Trit);
  tree -> SetBranchAddress("nSigmaTPC_Trit",&nSigmaTPC_Trit);
  tree -> SetBranchAddress("nSigmaTOF_Trit",&nSigmaTOF_Trit);
  tree -> SetBranchAddress("nSigmaTRD_Trit",&nSigmaTRD_Trit);
  tree -> SetBranchAddress("nSigmaHMPID_Trit",&nSigmaHMPID_Trit);
  tree -> SetBranchAddress("nSigmaITS_Deut",&nSigmaITS_Deut);
  tree -> SetBranchAddress("nSigmaTPC_Deut",&nSigmaTPC_Deut);
  tree -> SetBranchAddress("nSigmaTOF_Deut",&nSigmaTOF_Deut);
  tree -> SetBranchAddress("nSigmaTRD_Deut",&nSigmaTRD_Deut);
  tree -> SetBranchAddress("nSigmaHMPID_Deut",&nSigmaHMPID_Deut);
  //   tree -> SetBranchAddress("nSigmaITS_Prot",&nSigmaITS_Prot);
  //   tree -> SetBranchAddress("nSigmaTPC_Prot",&nSigmaTPC_Prot);
  //   tree -> SetBranchAddress("nSigmaTOF_Prot",&nSigmaTOF_Prot);
  //   tree -> SetBranchAddress("nSigmaTRD_Prot",&nSigmaTRD_Prot);
  //   tree -> SetBranchAddress("nSigmaHMPID_Prot",&nSigmaHMPID_Prot);
  //   tree -> SetBranchAddress("nSigmaITS_Pion",&nSigmaITS_Pion);
  //   tree -> SetBranchAddress("nSigmaTPC_Pion",&nSigmaTPC_Pion);
  //   tree -> SetBranchAddress("nSigmaTOF_Pion",&nSigmaTOF_Pion);
  //   tree -> SetBranchAddress("nSigmaTRD_Pion",&nSigmaTRD_Pion);
  //   tree -> SetBranchAddress("nSigmaHMPID_Pion",&nSigmaHMPID_Pion);
  //   tree -> SetBranchAddress("nSigmaITS_Kaon",&nSigmaITS_Kaon);
  //   tree -> SetBranchAddress("nSigmaTPC_Kaon",&nSigmaTPC_Kaon);
  //   tree -> SetBranchAddress("nSigmaTOF_Kaon",&nSigmaTOF_Kaon);
  //   tree -> SetBranchAddress("nSigmaTRD_Kaon",&nSigmaTRD_Kaon);
  //   tree -> SetBranchAddress("nSigmaHMPID_Kaon",&nSigmaHMPID_Kaon);
  //   tree -> SetBranchAddress("nSigmaITS_Elec",&nSigmaITS_Elec);
  //   tree -> SetBranchAddress("nSigmaTPC_Elec",&nSigmaTPC_Elec);
  //   tree -> SetBranchAddress("nSigmaTOF_Elec",&nSigmaTOF_Elec);
  //   tree -> SetBranchAddress("nSigmaTRD_Elec",&nSigmaTRD_Elec);
  //   tree -> SetBranchAddress("nSigmaHMPID_Elec",&nSigmaHMPID_Elec);

}
//___________________________________________________________________________________________
void TritonAnalysis::CreateHistograms()  {

  TString pTLabel= "#it{p}_{T} (GeV/#it{c})";

  fHistoPtSpectrum = new TH1D ("fHistoPtSpectrum","",nPtBins,Binning);
  fHistoPtSpectrum->GetXaxis()->SetTitle(pTLabel);
  fHistoPtSpectrum->GetXaxis()->SetTitleSize(labelsize);
  fHistoPtSpectrum->GetXaxis()->SetLabelSize(labelsize+0.01);
  fHistoPtSpectrum->GetYaxis()->SetTitle("#frac{d^{2}N}{d#it{y} d#it{p}_{T}} (GeV/#it{c})");
  fHistoPtSpectrum->GetYaxis()->SetTitleSize(labelsize+0.01);
  fHistoPtSpectrum->GetYaxis()->SetTitleOffset(0.8);
  fHistoPtSpectrum->GetYaxis()->SetLabelSize(labelsize);

  fHistoAntiPtSpectrum = new TH1D ("fHistoAntiPtSpectrum","",nPtBins,Binning);
  fHistoRatioSpectra = new TH1D ("fHistoRatioSpectra","",nPtBins,Binning);
  fHistoRatioSpectra->GetXaxis()->SetTitle(pTLabel);
  fHistoRatioSpectra->GetXaxis()->SetTitleSize(labelsize+0.01);
  fHistoRatioSpectra->GetXaxis()->SetLabelSize(labelsize);
  fHistoRatioSpectra->GetYaxis()->SetTitle("#bar{^{3}H} / ^{3}H");
  fHistoRatioSpectra->GetYaxis()->SetTitleSize(labelsize+0.01);
  fHistoRatioSpectra->GetYaxis()->SetTitleOffset(0.8);
  fHistoRatioSpectra->GetYaxis()->SetLabelSize(labelsize);

  fHistoPtSpectrum->Sumw2();
  fHistoAntiPtSpectrum->Sumw2();
  fHistoRatioSpectra->Sumw2();

  fHistoRapidity = new TH1D ("fHistoRapidity","Rapidity",30,-1.5,1.5);
  fHistoMean = new TH1D ("fHistoMean","Mean of the triton distribution",24,-3.,3.);
  fHistoWidth = new TH1D ("fHistoWidth","Width of the triton distribution",12,0.,3.);

}
//___________________________________________________________________________________________
void TritonAnalysis::ExecuteAnalysis()  {

  // defining size of pad ---------------------------
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);

  SaveToFile = false; // Save efficiency and pT slice fits to pdf
  Bool_t DrawTPCselection = false; // possibility to display the TPC selection in a plot
  Bool_t DrawDCADistribution = false; // Possibility to diplay the DCA distributions in xy and z for all candidates (before removal of triton)
  Bool_t DrawRapidity = false;

  Bool_t CorrectTOF = kTRUE;
  Bool_t EvaluatePrimaryFraction = kTRUE; // Estimate primary fraction via fit of the DCA disribution
  Bool_t UseFitter = kTRUE; //Primary fraction estimated using TFractionFitter, if false (1 - SecondaryFit/Total) is used
  Bool_t UseMeanAndSigmaFromData = kTRUE;
  Bool_t SubtractHypertritonFeedDown = kTRUE;


  // Cuts
  Double_t SelectTPCPlot = 6.; // Selection window for the TPC dE/dx plot used to fit and remove triton
  Double_t minFit = -6.; // Lower boundary for the triton fit interval
  Double_t maxFit = -2.; // Upper boundary for the triton fit interval
  Double_t maxEta = 0.9; // Pseudorapidity range

  // Parameters varied to evaluate the systematic uncertainties
  minRapidity = -1.;
  maxRapidity = 0.;
  ITShits = 2;
  TPCcluster = 120;

  DCAcutXY = 0.1;
  DCAcutZ = 1.;

  TPCselectMin = -3.0; // TPC selection window
  TPCselectMax = 3.0; // TPC selection window
  TOFselect = 3.; // TOF selection window

  // Evaluate the rapidity distributrion
  if (DrawRapidity) {
    TCanvas *cRapidity = new TCanvas("cRapidity","Rapidity");
    TLegend* LegALICERap = new TLegend(0.55,0.8,0.95,0.925);
    LegALICERap->SetTextSize(legsize);
    LegALICERap->SetFillColor(0);
    LegALICERap->SetFillStyle(0);
    LegALICERap->SetBorderSize(0);
    LegALICERap->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
    LegALICERap->AddEntry((TObject*)0,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");

    if (MultiDep){
      tree->Draw("CalculateRapidity(y)>>fHistoRapidity(30,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && multPercentile_V0A > %f && multPercentile_V0A < %f && TMath::Abs(nSigmaTOF_Trit) < %f", maxEta, ITShits, TPCcluster, TPCselectMin, TPCselectMax, DCAcutZ, DCAcutXY,MultiplicityPercentile[0], MultiplicityPercentile[1], TOFselect));
    } else {
      tree->Draw("CalculateRapidity(y)>>fHistoRapidity(30,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && TMath::Abs(nSigmaTOF_Trit) < %f", maxEta, ITShits, TPCcluster, TPCselectMin, TPCselectMax, DCAcutZ, DCAcutXY, TOFselect));
    }

    fHistoRapidity = (TH1D*) gDirectory->Get("fHistoRapidity");
    fHistoRapidity->SetTitle("");
    fHistoRapidity->GetXaxis()->SetTitle("#it{y}");
    fHistoRapidity->GetYaxis()->SetTitle("counts");
    fHistoRapidity->DrawCopy("E");
    LegALICERap->Draw();

    TLine UpperBoundary (0, 0, 0,fHistoRapidity->GetMaximum());
    TLine LowerBoundary (-1, 0, -1, fHistoRapidity->GetMaximum());
    UpperBoundary.SetLineColor(kGreen+2);
    LowerBoundary.SetLineColor(kGreen+2);
    UpperBoundary.SetLineStyle(2);
    LowerBoundary.SetLineStyle(2);
    UpperBoundary.DrawClone("same");
    LowerBoundary.DrawClone("same");
    cRapidity->SaveAs("Rapidity.pdf");
  }

  // Load efficiency
  // Default for the  multiplicity integrated result
  TFile *EffFile = TFile::Open ("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Triton/Efficiency.root");

  // For the evaluation of the systematic uncertainties due to the PID selection
  //   TFile *EffFile = TFile::Open (Form("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Triton/Systematics/PID/Efficiency_TPC_%.1f_%.1f.root",TPCselectMin, TPCselectMax));

  //   TFile *EffFile = TFile::Open (Form("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Triton/Systematics/PID/Efficiency_TOF_%.1f.root",TOFselect));

  // For the evaluation of the systematic uncertainties due tracking
  //   TFile *EffFile = TFile::Open (Form("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Triton/Systematics/Tracking/Efficiency_ITS%dTPC%d_%.2fy%.2f.root",ITShits, TPCcluster, minRapidity, maxRapidity));

  // For the evaluation of the systematic uncertainties due the secondary contamination via DCA variation
  //   TFile *EffFile = TFile::Open (Form("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Triton/Systematics/Secondaries/Efficiency_DCAxy%.2fz%.2f.root",DCAcutXY,DCAcutZ));

  TH1D* hEffFAST = (TH1D*) EffFile->GetObjectChecked("hEffTritonTOF", "TH1D");
  TH1D* hEffFASTAnti = (TH1D*) EffFile->GetObjectChecked("hEffAntiTritonTOF", "TH1D");

  // rebin efficiency
  for (Int_t i=1; i<=hEffFASTAnti->GetNbinsX(); i++) {
    hEffFAST->SetBinContent(i, hEffFAST->GetBinContent(i)*hEffFAST->GetBinWidth(i));
    hEffFAST->SetBinError(i, hEffFAST->GetBinError(i)*hEffFAST->GetBinWidth(i));

    hEffFASTAnti->SetBinContent(i, hEffFASTAnti->GetBinContent(i)*hEffFASTAnti->GetBinWidth(i));
    hEffFASTAnti->SetBinError(i, hEffFASTAnti->GetBinError(i)*hEffFASTAnti->GetBinWidth(i));
  }
  hEffFAST = (TH1D*) hEffFAST->Rebin(nPtBins,"hEffFAST",Binning);
  hEffFASTAnti = (TH1D*) hEffFASTAnti->Rebin(nPtBins,"hEffFASTAnti",Binning);

  hEffFAST->Scale(1,"width");
  hEffFASTAnti->Scale(1,"width");

  TCanvas* cEff = new TCanvas("cEff","Efficiencies");
  TLegend* LegALICE = new TLegend(0.12,0.75,0.5,0.85);
  LegALICE->SetTextSize(legsize);
  LegALICE->SetFillColor(0);
  LegALICE->SetFillStyle(0);
  LegALICE->SetBorderSize(0);
  LegALICE->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
  LegALICE->AddEntry((TObject*)0,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");

  gPad->SetTicks();
  TLegend* EffLeg = new TLegend(0.2,0.65,0.6,0.85,"ALICE #it{Work in progress}");
  EffLeg->SetTextSize(legsize);
  EffLeg->SetFillColor(0);
  EffLeg->SetFillStyle(0);
  EffLeg->SetBorderSize(0);
  EffLeg->AddEntry((TObject*)0,"#kern[-0.35]{p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}","");
  EffLeg->AddEntry(hEffFAST, "#varepsilon_{t}","l");
  EffLeg->AddEntry(hEffFASTAnti, "#varepsilon_{#bar{^{3}H}}","l");

  hEffFAST->SetTitle("");
  hEffFASTAnti->SetTitle("");
  hEffFAST->GetXaxis()->SetTitleSize(legsize+0.01);
  hEffFAST->GetXaxis()->SetLabelSize(legsize);
  hEffFAST->GetYaxis()->SetTitleSize(legsize+0.01);
  hEffFAST->GetYaxis()->SetLabelSize(legsize);
  hEffFAST->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hEffFASTAnti->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hEffFAST->GetYaxis()->SetTitle("#it{Acceptance x Efficieny}");
  hEffFAST->GetYaxis()->SetRangeUser(0, 1);
  hEffFAST->DrawCopy();
  hEffFASTAnti->DrawCopy("same");
  EffLeg->Draw();

  TH1D* hEffAnti = hEffFASTAnti->Clone("hEffAnti");
  TH1D* hEff = hEffFAST->Clone("hEff");
  if (SaveToFile) cEff->SaveAs("Efficiency.pdf");

  TCanvas* cEffRatio = new TCanvas("cEffRatio","Ratio of efficiencies");
  gPad->SetTicks();
  TH1D *hEffRatio = (TH1D*) hEffAnti->Clone("hEffRatio");
  hEffRatio->GetYaxis()->SetTitleOffset(0.8);
  hEffRatio->GetYaxis()->SetTitle("#varepsilon_{#bar{^{3}H}} / #varepsilon_{^{3}H}");
  hEffRatio->GetXaxis()->SetTitleSize(legsize+0.01);
  hEffRatio->GetXaxis()->SetLabelSize(legsize);
  hEffRatio->GetYaxis()->SetTitleSize(legsize+0.01);
  hEffRatio->GetYaxis()->SetLabelSize(legsize);
  hEffRatio->Divide(hEff);
  hEffRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hEffRatio->DrawCopy();
  LegALICE->Draw();
  if (SaveToFile) cEffRatio->SaveAs("EfficiencyRatio.pdf");

  // QA plot: TPC selection
  TLine UpperBoundary (0, 3, 10, 3);
  TLine LowerBoundary (0, -3, 10, -3);
  UpperBoundary.SetLineColor(kGreen+2);
  LowerBoundary.SetLineColor(kGreen+2);
  UpperBoundary.SetLineStyle(2);
  LowerBoundary.SetLineStyle(2);

  if (DrawTPCselection){ // possibility to display the TPC selection in a plot
    TH2D* TPCBuffer; // Buffer histogram for TPC dE/dx selection plot

    TCanvas* cTPC = new TCanvas("cTPC","TPC Selection triton");
    gPad->SetTicks();
    tree->Draw("nSigmaTPC_Trit:CorrectPt(px, py)>>TritonBuffer(100,0,10,100,-5,5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && CorrectPt(px, py) < 10 && TMath::Abs(nSigmaTOF_Trit) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, TOFselect));
    TPCBuffer = (TH2D*) gDirectory->Get("TritonBuffer");
    TPCBuffer->SetTitle("");
    TPCBuffer->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TPCBuffer->GetYaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}H}}{#sigma_{d#it{E} / d#it{x}}}");
    TPCBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetXaxis()->SetLabelSize(legsize);
    TPCBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetYaxis()->SetLabelSize(legsize);
    TPCBuffer->DrawCopy("");
    UpperBoundary.DrawClone("same");
    LowerBoundary.DrawClone("same");
    LegALICE->Draw();

    TCanvas* cTPCAnti = new TCanvas("cTPCAnti","TPC Selection^{3}#bar{H}");
    gPad->SetTicks();
    tree->Draw("nSigmaTPC_Trit:CorrectPt(px, py)>>AntiTritonBuffer(100,0,10,100,-5,5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && CorrectPt(px, py) < 10 && TMath::Abs(nSigmaTOF_Trit) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, TOFselect));
    TPCBuffer = (TH2D*) gDirectory->Get("AntiTritonBuffer");
    TPCBuffer->SetTitle("");
    TPCBuffer->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TPCBuffer->GetYaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}H}}{#sigma_{d#it{E} / d#it{x}}}");
    TPCBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetXaxis()->SetLabelSize(legsize);
    TPCBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetYaxis()->SetLabelSize(legsize);
    TPCBuffer->DrawCopy("");
    UpperBoundary.DrawClone("same");
    LowerBoundary.DrawClone("same");
    LegALICE->Draw();

    cTPC->SaveAs("TPCselection.pdf");
    cTPCAnti->SaveAs("TPCselection_Anti.pdf");
  }

  TF1* TOFGaus = new TF1("TOFGaus","gaus",-10,10); // To fit nsigma dist.

  Double_t Ptbin = 0;
  // TOF distribution
  Double_t meanTOF[nPtBins];
  Double_t sigmaTOF[nPtBins];
  for (Int_t jbins=0; jbins < nPtBins; jbins++) {
    meanTOF[jbins] = 0;
    sigmaTOF[jbins] = 1;
  }
  if(CorrectTOF){
    TCanvas* cTOF[nPtBins];
    TH1D *hTempTOF[nPtBins];
    for (Int_t i =0; i < nPtBins; i++) {
      Ptbin = (Binning[i]+Binning[i+1])/2;
      cTOF[i] = new TCanvas();
      cTOF[i]->SetTitle(Form("TOF %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
      cTOF[i]->Divide(2,1);
      cTOF[i]->cd(1);
      gPad->SetTicks();
      tree->Draw(Form("nSigmaTOF_Trit>>TOFtemp_%d(70,-7,7)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, TPCselectMin, TPCselectMax, Binning[i], Binning[i+1]));
      hTempTOF[i] = (TH1D*) gDirectory->Get(Form("TOFtemp_%d",i));
      hTempTOF[i]->SetTitle("Wide DCA for calibration");
      hTempTOF[i]->GetXaxis()->SetTitleSize(legsize+0.01);
      hTempTOF[i]->GetXaxis()->SetLabelSize(legsize);
      hTempTOF[i]->GetYaxis()->SetTitleSize(legsize+0.01);
      hTempTOF[i]->GetYaxis()->SetLabelSize(legsize);
      hTempTOF[i]->GetXaxis()->SetTitle("#frac{TOF #it{t}-#it{t}_{triton}}{#sigma_{triton}}");
      hTempTOF[i]->GetYaxis()->SetTitle("Counts");
      hTempTOF[i]->GetXaxis()->SetTitleOffset(1.1);
      hTempTOF[i]->GetYaxis()->SetTitleOffset(1.5);
      hTempTOF[i]->DrawCopy();

      TFitResultPtr Fit = hTempTOF[i]->Fit(TOFGaus, "NMLQ", "",-3.,3.);
      TOFGaus->SetRange(-10,10);
      TOFGaus->Draw("same");
      if (Ptbin < 2.5 ) {
        meanTOF[i] = TOFGaus->GetParameter(1);
        sigmaTOF[i] = TOFGaus->GetParameter(2);
      } else {
        meanTOF[i] = 0;
        sigmaTOF[i] = 1;
      }
      cTOF[i]->cd(2);
      gPad->SetTicks();
      tree->Draw(Form("(nSigmaTOF_Trit-%f)/%f>>TOFtemp_%d(70,-7,7)",meanTOF[i], sigmaTOF[i], i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f ", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, TPCselectMin, TPCselectMax, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY));

      hTempTOF[i] = (TH1D*) gDirectory->Get(Form("TOFtemp_%d",i));
      hTempTOF[i]->SetTitle("Calibrated sample with analysis selection");
      hTempTOF[i]->GetXaxis()->SetTitleSize(legsize+0.01);
      hTempTOF[i]->GetXaxis()->SetLabelSize(legsize);
      hTempTOF[i]->GetYaxis()->SetTitleSize(legsize+0.01);
      hTempTOF[i]->GetYaxis()->SetLabelSize(legsize);
      hTempTOF[i]->GetXaxis()->SetTitle("#frac{TOF #it{t}-#it{t}_{triton}}{#sigma_{triton}} (corrected)");
      hTempTOF[i]->GetYaxis()->SetTitle("Counts");
      hTempTOF[i]->GetXaxis()->SetTitleOffset(1.1);
      hTempTOF[i]->GetYaxis()->SetTitleOffset(1.5);
      hTempTOF[i]->DrawCopy();
      cTOF[i]->SaveAs(Form("TOF_%.2fGeV.pdf",Ptbin));
    }
  }

  //Check DCA_XY
  TCanvas* cDCAxy[nPtBins];
  TCanvas* cDCAz[nPtBins];
  Double_t PrimaryFraction[nPtBins];
  for (Int_t jbin=0; jbin < nPtBins; jbin++) {
    PrimaryFraction[jbin]=1;
  }
  Double_t PrimaryFractionUncertainty = 0;
  if (DrawDCADistribution || EvaluatePrimaryFraction) {
    TH1D* TritonBuffer;
    TH1D* AntiTritonBuffer;

    // Buffer for the DCA templates
    TH1D* DCATemplatePrimary;
    TH1D* DCATemplateSecondary;
    TObjArray *MCDCATemplates = new TObjArray(2);        // MC histograms are put in this array
    TFractionFitter* DCAFitter;
    // Buffer for template fit results
    Int_t FitStatus = -1;
    Double_t Secondary = 0., Primary = 1., ErrorPrimary = 0., ErrorSecondary = 0.;
    TF1* SecondaryTritionWoPeak = new TF1("SecondaryTritionWoPeak", BackgroundFitWithoutSignalRegion,-6,6,3);
    if (SecondaryGaus) {
      SecondaryTritionWoPeak->SetParameters(100,0,20);
      TF1* SecondaryTrition = new TF1("SecondaryTrition","gaus",-6,6);
    } else {
      TF1* SecondaryTrition = new TF1("SecondaryTrition","pol2",-6,6);
    }

    TFile *PrimaryTemplateFile = TFile::Open ("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Triton/DCATemplates_Triton_wTOF_New.root");
    //DCA output file
    TFile *DCAOutputFile = new TFile ("Triton_DCAxy.root","RECREATE");
    DCAOutputFile->cd();

    for (Int_t i =0; i < nPtBins; i++) {
      Ptbin = (Binning[i]+Binning[i+1])/2;
      if (Ptbin > 2.5) continue; // change to 3 GeV for DCA variation
      cDCAxy[i] = new TCanvas(Form("DCAxy_%d",i), Form("DCA %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
      gPad->SetTicks();
      gPad->SetLogy();
      Double_t LowerTPCforDCA = -3.0;

      if (MultiDep){
        tree->Draw("dcaxy>>TritonBuffer(60,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f && TMath::Abs(dcaz) <= %f && multPercentile_V0A > %f && multPercentile_V0A < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], meanTOF[i], sigmaTOF[i]*TOFselect, DCAcutZ, MultiplicityPercentile[0], MultiplicityPercentile[1]));
        // Anti-Particle
        tree->Draw("dcaxy>>AntiTritonBuffer(60,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f && TMath::Abs(dcaz) <= %f && multPercentile_V0A > %f && multPercentile_V0A < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], meanTOF[i], sigmaTOF[i]*TOFselect, DCAcutZ, MultiplicityPercentile[0], MultiplicityPercentile[1]));
      } else {
        // Particle
        tree->Draw("dcaxy>>TritonBuffer(40,-1.0,1.0)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f && TMath::Abs(dcaz) <= %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], meanTOF[i], sigmaTOF[i]*TOFselect, DCAcutZ));
        // Anti-Particle
        tree->Draw("dcaxy>>AntiTritonBuffer(40,-1.0,1.0)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f && TMath::Abs(dcaz) <= %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], meanTOF[i], sigmaTOF[i]*TOFselect, DCAcutZ));
      }

      TritonBuffer = (TH1D*) gDirectory->Get("TritonBuffer");
      AntiTritonBuffer = (TH1D*) gDirectory->Get("AntiTritonBuffer");
      TritonBuffer->SetTitle("");
      TritonBuffer->GetXaxis()->SetTitle("DCA_{xy} (cm)");
      TritonBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
      TritonBuffer->GetXaxis()->SetLabelSize(legsize);
      TritonBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
      TritonBuffer->GetYaxis()->SetLabelSize(legsize);
      TritonBuffer->GetYaxis()->SetTitle("Counts");

      AntiTritonBuffer->SetTitle("");
      AntiTritonBuffer->GetXaxis()->SetTitle("DCA_{xy} (cm)");
      AntiTritonBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
      AntiTritonBuffer->GetXaxis()->SetLabelSize(legsize);
      AntiTritonBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
      AntiTritonBuffer->GetYaxis()->SetLabelSize(legsize);
      AntiTritonBuffer->GetYaxis()->SetTitle("Counts");

      AntiTritonBuffer->SetMarkerColor(kRed);
      AntiTritonBuffer->SetLineColor(kRed);
      TritonBuffer->SetMarkerColor(kBlue);
      TritonBuffer->SetLineColor(kBlue);
      AntiTritonBuffer->SetMarkerStyle(20);
      AntiTritonBuffer->SetMarkerSize(0.5);
      TritonBuffer->SetMarkerStyle(24);
      TritonBuffer->SetMarkerSize(0.5);
      TritonBuffer->GetYaxis()->SetRangeUser(1,100);
      TritonBuffer->DrawCopy("");
      AntiTritonBuffer->DrawCopy("same");

      if (EvaluatePrimaryFraction){
        TFitResultPtr FitSecBkg = TritonBuffer->Fit(SecondaryTritionWoPeak, "NMLQ", "", -3, 3);
        SecondaryTritionWoPeak->SetLineColor(kBlack);
        SecondaryTritionWoPeak->Draw("same");
        SecondaryTrition->SetParameters(SecondaryTritionWoPeak->GetParameters());

//        //Load templates (data driven)
//        DCATemplatePrimary = (TH1D*) AntiTritonBuffer->Clone(Form("hDCAxyPrimary_%.1f", Ptbin)); // Use Anti-Triton as proxy for the primary Triton
//        
//        // Use fit to the DCAxy distribution as proxy for the secondary distribution (Functional form to be varied for syst.)
//        DCATemplateSecondary = (TH1D*) TritonBuffer->Clone(Form("hDCAxyMaterial_%.1f", Ptbin));
//        DCATemplateSecondary->TH1::Reset();
//        DCATemplatePrimary->Sumw2();
//        DCATemplateSecondary->Sumw2();
//        DCATemplateSecondary->Add(SecondaryTrition);
//        for (Int_t jbin=0; jbin<DCATemplateSecondary->GetNbinsX(); jbin++) {
//          DCATemplateSecondary->SetBinError(jbin,TMath::Sqrt(DCATemplateSecondary->GetBinContent(jbin)));
//        }

        //Load templates (MC)
        DCATemplatePrimary = (TH1D*) PrimaryTemplateFile->GetObjectChecked(Form("hDCAxy_Triton_Primary_%.2f", Ptbin), "TH1D");
        DCATemplatePrimary->Rebin(4); // Adjust binning
        DCATemplateSecondary = (TH1D*) PrimaryTemplateFile->GetObjectChecked(Form("hDCAxy_Triton_Material_%.2f", Ptbin), "TH1D");
        DCATemplateSecondary->Rebin(4); // Adjust binning

        //Empty TObjArray
        MCDCATemplates->Clear();
        //Fill TObjArray
        MCDCATemplates->Add(DCATemplatePrimary);
        MCDCATemplates->Add(DCATemplateSecondary);

        if (UseFitter) {
          //Fit templates
          DCAFitter = new TFractionFitter(TritonBuffer, MCDCATemplates,"Q");
          DCAFitter->Constrain(1,0.0,1.0);
          DCAFitter->Constrain(2,0.0,1.0);
          FitStatus = DCAFitter->Fit();
          printf("Status of the DCA template fit for %.2f GeV/c: %d \n",Ptbin, FitStatus);               // perform the fit and prints status
          if (FitStatus == 0) {                       // check on fit status
            TH1F* result = (TH1F*) DCAFitter->GetPlot();
            result->SetMarkerColor(kMagenta-3);
            result->SetLineColor(kMagenta-3);
            result->SetMarkerSize(0.5);
            result->Draw("same");
            // Get fractions
            DCAFitter->GetResult(0, Primary, ErrorPrimary);
            DCAFitter->GetResult(1, Secondary, ErrorSecondary);
            // Get adjusted DCA templates
            TH1D* hPrim = (TH1D*)DCAFitter->GetMCPrediction(0);
            TH1D* hSec = (TH1D*)DCAFitter->GetMCPrediction(1);
            Float_t DataIntegral = result->Integral();
            hPrim->Scale(DataIntegral * Primary / hPrim->Integral());
            hSec->Scale(DataIntegral * Secondary / hSec->Integral());
            hSec->SetLineColor(kGreen-3);
            hPrim->SetLineColor(kOrange-3);
            hSec->SetMarkerColor(kGreen-3);
            hPrim->SetMarkerColor(kOrange-3);
            hSec->DrawCopy("same");
            hPrim->DrawCopy("same");

            PrimaryFraction[i] = hPrim->Integral(hPrim->FindBin(-(DCAcutXY-0.001)), hPrim->FindBin((DCAcutXY-0.001))) / result->Integral(result->FindBin(-(DCAcutXY-0.001)), result->FindBin((DCAcutXY-0.001)));

            Double_t relativeFracionError = ErrorPrimary/Primary;
            Double_t IntegralError;
            Primary = hPrim->IntegralAndError(hPrim->FindBin(-(DCAcutXY-0.001)), hPrim->FindBin((DCAcutXY-0.001)),IntegralError);
            IntegralError = IntegralError/Primary;
            ErrorPrimary = Primary*sqrt(pow(relativeFracionError,2)+pow(IntegralError,2));

            Double_t relativeFracionError = ErrorSecondary/Secondary;
            IntegralError=0;
            Secondary = hSec->IntegralAndError(hSec->FindBin(-(DCAcutXY-0.001)), hSec->FindBin((DCAcutXY-0.001)),IntegralError);
            IntegralError = IntegralError/Secondary;
            ErrorSecondary = Secondary*sqrt(pow(relativeFracionError,2)+pow(IntegralError,2));

            PrimaryFractionUncertainty = abs(ErrorPrimary*Secondary-Primary*ErrorSecondary)/pow(Primary+Secondary,2);

            TritonBuffer->Write(Form("TritonDCAxy%.2fGeV",Ptbin));
            AntiTritonBuffer->Write(Form("AntiTritonDCAxy%.2fGeV",Ptbin));
            hPrim->Write(Form("PrimaryTemplate%.2fGeV",Ptbin));
            hSec->Write(Form("SecondaryTemplate%.2fGeV",Ptbin));
            result->Write(Form("TemplateFit%.2fGeV",Ptbin));
          }
        } else {
          PrimaryFraction[i] = 1.-(SecondaryTrition->Integral(-DCAcutXY,DCAcutXY) / TritonBuffer->Integral(TritonBuffer->FindBin(-(DCAcutXY-0.001)),TritonBuffer->FindBin((DCAcutXY-0.001)),"width"));
          PrimaryFractionUncertainty = SecondaryTrition->IntegralError(-DCAcutXY,DCAcutXY) / TritonBuffer->Integral(TritonBuffer->FindBin(-(DCAcutXY-0.01)),TritonBuffer->FindBin((DCAcutXY-0.01)),"width"); //Dummy, not complet
        }

        printf("Primary fraction:  %.3f ±  %.3f\n", PrimaryFraction[i],PrimaryFractionUncertainty);
        printf("Anti-Triton fraction: %.2f \n", AntiTritonBuffer->Integral(AntiTritonBuffer->FindBin(-DCAcutXY), AntiTritonBuffer->FindBin(DCAcutXY),"width") / TritonBuffer->Integral(TritonBuffer->FindBin(-DCAcutXY), TritonBuffer->FindBin(DCAcutXY),"width"));
      }
      LegALICE->Draw();
      cDCAxy[i]->SaveAs(Form("DCAxy%.2fGeV.pdf",Ptbin));

      if (DrawDCADistribution) {
        cDCAz[i] = new TCanvas(Form("DCAz_%d",i), Form("DCA %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
        gPad->SetTicks();
        gPad->SetLogy();
        // Particle
        tree->Draw("dcaz>>TritonBuffer(60,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaxy) <= %f && TMath::Abs(nSigmaTOF_Trit) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], DCAcutXY, TOFselect));
        // Anti-Particle
        tree->Draw("dcaz>>AntiTritonBuffer(60,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && nSigmaTPC_Trit > %f && nSigmaTPC_Trit < %f  &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaxy) <= %f && TMath::Abs(nSigmaTOF_Trit) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], DCAcutXY, TOFselect));

        TritonBuffer = (TH1D*) gDirectory->Get("TritonBuffer");
        AntiTritonBuffer = (TH1D*) gDirectory->Get("AntiTritonBuffer");
        TritonBuffer->SetTitle("");
        TritonBuffer->GetXaxis()->SetTitle("DCA_z (cm)");
        TritonBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
        TritonBuffer->GetXaxis()->SetLabelSize(legsize);
        TritonBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
        TritonBuffer->GetYaxis()->SetLabelSize(legsize);
        AntiTritonBuffer->SetMarkerColor(kRed);
        AntiTritonBuffer->SetLineColor(kRed);
        TritonBuffer->SetMarkerColor(kBlue);
        TritonBuffer->SetLineColor(kBlue);
        AntiTritonBuffer->SetMarkerStyle(20);
        AntiTritonBuffer->SetMarkerSize(0.5);
        TritonBuffer->SetMarkerStyle(24);
        TritonBuffer->SetMarkerSize(0.5);
        TritonBuffer->DrawCopy("");
        AntiTritonBuffer->DrawCopy("same");
        LegALICE->Draw();
        cDCAz[i]->SaveAs(Form("DCAz%.2fGeV.pdf",Ptbin));
      }
    }
    DCAOutputFile->Close();
  }

  TCanvas* cFit[nPtBins];
  TCanvas* cFitAnti[nPtBins];
  TH1D *hTempTPC[nPtBins];
  TH1D *hTempBackgroundSubtracted[nPtBins];

  TF1* Deuteron = new TF1("Deuteron","gaus",-10,10); // To fit the contamination
  Deuteron->SetParameter(0, 10); // constant
  Deuteron->SetParameter(1, -5); // mean
  Deuteron->SetParameter(2, 2); // width
  Deuteron->SetLineColor(kRed);

  // yield extraction
  TF1* Gaus = new TF1("Gaus","gaus",-10,10); // To fit nsigma dist.
  Gaus->SetParameter(1, 0);
  Gaus->SetParameter(2, 1);

  Int_t AntiTritonCounter = 0;
  for (Int_t i =0; i < nPtBins; i++) {
    Ptbin = (Binning[i]+Binning[i+1])/2;
    Double_t mean = 0;
    Double_t sigma = 1;
    Double_t nBackBin = 0;

    // evaluation of (primary) triton
    cFit[i] = new TCanvas();
    cFit[i]->SetTitle(Form("Fit %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
    cFit[i]->Divide(2,1);
    cFit[i]->cd(1);
    gPad->SetTicks();
    gPad->SetTopMargin(0.07);

    if (MultiDep){
      tree->Draw(Form("nSigmaTPC_Trit>>Bin%d(%d,%f,%f)",i,(int) SelectTPCPlot*4,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && TMath::Abs(nSigmaTPC_Trit) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && multPercentile_V0A > %f && multPercentile_V0A < %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY,MultiplicityPercentile[0], MultiplicityPercentile[1], meanTOF[i], sigmaTOF[i]*TOFselect));

    } else {
      tree->Draw(Form("nSigmaTPC_Trit>>Bin%d(%d,%f,%f)",i,(int) SelectTPCPlot*4,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && TMath::Abs(nSigmaTPC_Trit) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY, meanTOF[i], sigmaTOF[i]*TOFselect));
    }
    hTempTPC[i] = (TH1D*) gDirectory->Get(Form("Bin%d",i));
    hTempTPC[i]->GetXaxis()->SetRangeUser(-SelectTPCPlot, SelectTPCPlot);
    hTempTPC[i]->SetTitle(Form("Triton %.2f < #it{p}_{T} < %.2f",Binning[i], Binning[i+1]));
    hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}H}}{#sigma_{d#it{E} / d#it{x}}}");
    hTempTPC[i]->GetXaxis()->SetTitleOffset(1.1);
    hTempTPC[i]->GetYaxis()->SetTitle("Counts");
    hTempTPC[i]->GetYaxis()->SetTitleOffset(1.5);
    hTempTPC[i]->DrawCopy("E");
    LegALICE->Draw();

    if (Ptbin > 2.) {
      TFitResultPtr FitDeuteron = hTempTPC[i]->Fit(Deuteron, "NQL", "", minFit, maxFit);
      Deuteron->SetLineStyle(2);
      Deuteron->SetRange(-SelectTPCPlot, SelectTPCPlot);
      Deuteron->DrawCopy("same");
      Deuteron->SetRange(minFit, maxFit);
      Deuteron->SetLineStyle(1);
      Deuteron->DrawCopy("same");
    }
    cFit[i]->cd(2);
    gPad->SetTicks();
    gPad->SetTopMargin(0.07);

    Deuteron->SetRange(-5, 5);
    hTempBackgroundSubtracted[i] = (TH1D*) hTempTPC[i]->Clone("hTempBackgroundSubtracted");
    if (Ptbin > 2.) {
      if (Deuteron->GetParameter(0) > 0) { // Do not subtract if function is negative
        hTempBackgroundSubtracted[i]->Add(Deuteron, -1);
      }
    }
    hTempBackgroundSubtracted[i]->GetXaxis()->SetRangeUser(-3.5, 3.5);
    hTempBackgroundSubtracted[i]->DrawCopy("E");

    TFitResultPtr Fit = hTempBackgroundSubtracted[i]->Fit(Gaus, "NMQL", "", -3., 3.);
    Gaus->DrawCopy("same");
    fHistoMean->Fill(Gaus->GetParameter(1));
    fHistoWidth->Fill(Gaus->GetParameter(2));
    if (UseMeanAndSigmaFromData) {
      mean = Gaus->GetParameter(1);
      sigma = Gaus->GetParameter(2);
    }
    Double_t nTritonBin = hTempBackgroundSubtracted[i]->Integral(hTempBackgroundSubtracted[i]->FindBin(mean+TPCselectMin*sigma), hTempBackgroundSubtracted[i]->FindBin(mean+TPCselectMax*sigma));
    nTritonBin = nTritonBin*PrimaryFraction[i]; // correct for secondaries
    fHistoPtSpectrum->Fill(Ptbin, nTritonBin);
    fHistoPtSpectrum->SetBinError(fHistoPtSpectrum->FindBin(Ptbin), TMath::Sqrt(nTritonBin));
    if (SaveToFile) cFit[i]->SaveAs(Form("Triton%.2fGeV.pdf",Ptbin));

    if (Ptbin > 2.) {
      Double_t NallCand = hTempTPC[i]->Integral(hTempTPC[i]->FindBin(mean+TPCselectMin*sigma), hTempTPC[i]->FindBin(mean+TPCselectMax*sigma));
      Double_t NContamination = Deuteron->Integral(mean+TPCselectMin*sigma,mean+TPCselectMax*sigma);
      printf("Fraction of contamination: %.1f % \n", NContamination/NallCand *100);
    }

//    // Input for the PID plot in the paper
//    if (Ptbin == 2.25) {
//      tree->Draw(Form("(nSigmaTPC_Trit-%f)/%f>>Bin%d(%d,%f,%f)",mean, sigma, i,(int) SelectTPCPlot*4,-SelectTPCPlot/sigma,SelectTPCPlot/sigma),Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && TMath::Abs(nSigmaTPC_Trit) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY, meanTOF[i], sigmaTOF[i]*TOFselect));
//      hTempTPC[i] = (TH1D*) gDirectory->Get(Form("Bin%d",i));
////      hTempTPC[i]->GetXaxis()->SetRangeUser(-SelectTPCPlot/sigma, SelectTPCPlot/sigma);
//      hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
//      hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
//      hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
//      hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
//      hTempTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}H}}{#sigma_{d#it{E} / d#it{x}}}");
//      hTempTPC[i]->GetXaxis()->SetTitleOffset(1.1);
//      hTempTPC[i]->GetYaxis()->SetTitle("Counts");
//      hTempTPC[i]->GetYaxis()->SetTitleOffset(1.5);
//      hTempTPC[i]->SetTitle(Form("%.1f GeV/#it{c} < #it{p}_{T} < %.1f GeV/#it{c}",Binning[i], Binning[i+1]));
//      Deuteron->SetParameter(1,(Deuteron->GetParameter(1)-mean)/sigma);
//      Deuteron->SetParameter(2,Deuteron->GetParameter(2)/sigma);
//      hTempTPC[i]->GetListOfFunctions()->Add(Deuteron);
//      hTempTPC[i]->SaveAs(Form("Triton%.2fGeV_Sigma.root",Ptbin));
//    }

    // evaluation of anti-triton
    cFitAnti[i] = new TCanvas();
    cFitAnti[i]->SetTitle(Form("Fit %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
    cFitAnti[i]->Divide(2,1);
    cFitAnti[i]->cd(1);
    gPad->SetTicks();
    gPad->SetTopMargin(0.07);

    if (MultiDep){
      tree->Draw(Form("nSigmaTPC_Trit>>AntiBin%d(%d,%f,%f)",i,(int) SelectTPCPlot*4,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && TMath::Abs(nSigmaTPC_Trit) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && multPercentile_V0A > %f && multPercentile_V0A < %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY,MultiplicityPercentile[0], MultiplicityPercentile[1], meanTOF[i], sigmaTOF[i]*TOFselect));
    } else {
      tree->Draw(Form("nSigmaTPC_Trit>>AntiBin%d(%d,%f,%f)",i,(int) SelectTPCPlot*4,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && TMath::Abs(nSigmaTPC_Trit) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && TMath::Abs(nSigmaTOF_Trit-%f) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY, meanTOF[i], sigmaTOF[i]*TOFselect));
    }
    hTempTPC[i] = (TH1D*) gDirectory->Get(Form("AntiBin%d",i));
    hTempTPC[i]->GetXaxis()->SetRangeUser(-SelectTPCPlot, SelectTPCPlot);
    hTempTPC[i]->SetTitle(Form("Anti-triton %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
    hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}H}}{#sigma_{d#it{E} / d#it{x}}}");
    hTempTPC[i]->GetXaxis()->SetTitleOffset(1.1);
    hTempTPC[i]->GetYaxis()->SetTitle("Counts");
    hTempTPC[i]->GetYaxis()->SetTitleOffset(1.5);
    hTempTPC[i]->DrawCopy("E");
    LegALICE->Draw();

    if (Ptbin > 2.) {
      TFitResultPtr FitDeuteron = hTempTPC[i]->Fit(Deuteron, "NQL", "", minFit, maxFit);
      Deuteron->SetLineStyle(2);
      Deuteron->SetRange(-SelectTPCPlot, SelectTPCPlot);
      Deuteron->DrawCopy("same");
      Deuteron->SetRange(minFit, maxFit);
      Deuteron->SetLineStyle(1);
      Deuteron->DrawCopy("same");
    }

    cFitAnti[i]->cd(2);
    gPad->SetTicks();
    gPad->SetTopMargin(0.07);

    Deuteron->SetRange(-5, 5);
    hTempBackgroundSubtracted[i] = (TH1D*) hTempTPC[i]->Clone("hTempBackgroundSubtracted");
    if (Ptbin > 2.) {
      if (Deuteron->GetParameter(0) > 0) { // Do not subtract if function is negative
        hTempBackgroundSubtracted[i]->Add(Deuteron, -1);
      }
    }

    hTempBackgroundSubtracted[i]->GetXaxis()->SetRangeUser(-3.5, 3.5);
    hTempBackgroundSubtracted[i]->DrawCopy("E");

    TFitResultPtr Fit = hTempBackgroundSubtracted[i]->Fit(Gaus, "NQL", "", -3., 3.);
    Gaus->DrawCopy("same");
    fHistoMean->Fill(Gaus->GetParameter(1));
    fHistoWidth->Fill(Gaus->GetParameter(2));
    if (UseMeanAndSigmaFromData) {
      mean = Gaus->GetParameter(1);
      sigma = Gaus->GetParameter(2);
    }
    Double_t nTritonBin = hTempBackgroundSubtracted[i]->Integral(hTempBackgroundSubtracted[i]->FindBin(mean+TPCselectMin*sigma), hTempBackgroundSubtracted[i]->FindBin(mean+TPCselectMax*sigma));
    fHistoAntiPtSpectrum->Fill(Ptbin, nTritonBin);
    fHistoAntiPtSpectrum->SetBinError(fHistoPtSpectrum->FindBin(Ptbin), TMath::Sqrt(nTritonBin));
    if (SaveToFile) cFitAnti[i]->SaveAs(Form("Anti-Triton%.2fGeV.pdf",Ptbin));
    AntiTritonCounter += nTritonBin;


    if (Ptbin > 2.) {
      Double_t NallCand = hTempTPC[i]->Integral(hTempTPC[i]->FindBin(mean+TPCselectMin*sigma), hTempTPC[i]->FindBin(mean+TPCselectMax*sigma));
      Double_t NContamination = Deuteron->Integral(mean+TPCselectMin*sigma,mean+TPCselectMax*sigma);
      printf("Fraction of contamination (Anti): %.1f % \n", NContamination/NallCand *100);
    }
  }

  //    normalization (number of events, bin width)
  Double_t RapidityInterval = TMath::Abs(maxRapidity-minRapidity);
  printf("Rapidity interval: %.1f \n", RapidityInterval);
  fHistoPtSpectrum -> Scale(1./(NumberOfEvents*RapidityInterval),"width");
  fHistoAntiPtSpectrum -> Scale(1./(NumberOfEvents*RapidityInterval),"width");

  // correct for eff.
  fHistoPtSpectrum -> Divide(hEff);
  fHistoAntiPtSpectrum -> Divide(hEffAnti);

  if (SubtractHypertritonFeedDown) {
    // subtract hypertriton contamination
    TFile *FeedDownFile = TFile::Open ("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/HypertritonFeedDown.root");
    TH1D* hFeedDownFraction = (TH1D*) FeedDownFile->GetObjectChecked("hFeedDownFractionHelium", "TH1D");
    TH1D* hFeedDownFractionAnti = (TH1D*) FeedDownFile->GetObjectChecked("hFeedDownFractionAntiHelium", "TH1D");

    for (Int_t i=1; i<=hFeedDownFraction->GetNbinsX(); i++) {
      hFeedDownFraction->SetBinContent(i, hFeedDownFraction->GetBinContent(i)*hFeedDownFraction->GetBinWidth(i));
      hFeedDownFractionAnti->SetBinContent(i, hFeedDownFractionAnti->GetBinContent(i)*hFeedDownFractionAnti->GetBinWidth(i));

      hFeedDownFraction->SetBinError(i, hFeedDownFraction->GetBinError(i)*hFeedDownFraction->GetBinWidth(i));
      hFeedDownFractionAnti->SetBinError(i, hFeedDownFractionAnti->GetBinError(i)*hFeedDownFractionAnti->GetBinWidth(i));
    }
    hFeedDownFraction = (TH1D*) hFeedDownFraction->Rebin(nPtBins,"hFeedDownFraction",Binning);
    hFeedDownFractionAnti = (TH1D*) hFeedDownFractionAnti->Rebin(nPtBins,"hFeedDownFractionAnti",Binning);

    hFeedDownFraction->Scale(-0.5,"width"); // Factor of 0.5 because the BR of hypertriton to triton is half the one to helium
    hFeedDownFractionAnti->Scale(-0.5,"width"); // Factor of 0.5 because the BR of hypertriton to triton is half the one to helium

    hFeedDownFraction->Multiply(fHistoPtSpectrum);
    hFeedDownFractionAnti->Multiply(fHistoAntiPtSpectrum);

    TCanvas *cFeedDownSubtraction = new TCanvas("cFeedDownSubtraction","Feed down subtraction",800,600);
    gPad->SetLogy();
    gPad->SetTicks();

    fHistoAntiPtSpectrum->SetMarkerColor(kRed);
    fHistoAntiPtSpectrum->SetLineColor(kRed);
    fHistoPtSpectrum->SetMarkerColor(kBlue);
    fHistoPtSpectrum->SetLineColor(kBlue);

    fHistoAntiPtSpectrum->SetMarkerStyle(20);
    fHistoAntiPtSpectrum->SetMarkerSize(0.5);
    fHistoPtSpectrum->SetMarkerStyle(24);
    fHistoPtSpectrum->SetMarkerSize(0.5);

    fHistoPtSpectrum->DrawCopy("EX0");
    fHistoAntiPtSpectrum->DrawCopy("EX0same");

    hFeedDownFraction->SetMarkerColor(kMagenta+3);
    hFeedDownFraction->SetLineColor(kMagenta+3);
    hFeedDownFractionAnti->SetMarkerColor(kAzure+3);
    hFeedDownFractionAnti->SetLineColor(kAzure+3);
    hFeedDownFractionAnti->SetMarkerStyle(20);
    hFeedDownFractionAnti->SetMarkerSize(0.5);
    hFeedDownFraction->SetMarkerStyle(24);
    hFeedDownFraction->SetMarkerSize(0.5);

    hFeedDownFraction->DrawCopy("EX0same");
    hFeedDownFractionAnti->DrawCopy("EX0same");

    fHistoPtSpectrum->Add(hFeedDownFraction);
    fHistoAntiPtSpectrum->Add(hFeedDownFractionAnti);

    fHistoPtSpectrum->DrawCopy("EX0same");
    fHistoAntiPtSpectrum->DrawCopy("EX0same");
  }

  // Anti- /  ratio
  fHistoRatioSpectra -> Divide(fHistoAntiPtSpectrum,fHistoPtSpectrum);

  printf("\n %d #bar{t} found \n", AntiTritonCounter); // print the total number of anti- candidates in the terminal

}
//___________________________________________________________________________________________
void TritonAnalysis::WriteOutputFile()  {

  TCanvas *cTritonQA = new TCanvas("cTritonQA","Mean and width of the triton distribution");
  fHistoMean->SetMarkerColor(kRed);
  fHistoMean->SetLineColor(kRed);
  fHistoMean->SetMarkerStyle(20);
  fHistoMean->SetMarkerSize(0.5);

  fHistoWidth->SetMarkerColor(kBlue);
  fHistoWidth->SetLineColor(kBlue);
  fHistoWidth->SetMarkerStyle(20);
  fHistoWidth->SetMarkerSize(0.5);

  fHistoMean->DrawCopy("E");
  fHistoWidth->DrawCopy("Esame");

  TCanvas *cSpectra = new TCanvas("cSpectra","Spectra and Ratio");
  cSpectra->Divide(1,2,0,0); // 0 are margins
  cSpectra->cd(1);
  gStyle->SetPadBottomMargin(0);
  gPad->SetLogy();
  gPad->SetTicks();

  TH1F* hAxis = new TH1F("hAxis","",(PlotMaxPt-PlotMinPt)*10,PlotMinPt,PlotMaxPt);
  hAxis->GetYaxis()->SetRangeUser(8e-9,5e-6);
  hAxis->GetXaxis()->SetTitle(pTLabel);
  hAxis->GetYaxis()->SetTitle(SpectraLabel);
  hAxis->GetXaxis()->SetNdivisions(505);

  hAxis->GetXaxis()->SetTitleSize(labelsize+0.01);
  hAxis->GetXaxis()->SetLabelSize(labelsize);
  hAxis->GetYaxis()->SetTitleSize(labelsize+0.01);
  hAxis->GetYaxis()->SetTitleOffset(0.8);
  hAxis->GetYaxis()->SetLabelSize(labelsize);

  TLegend* leg = new TLegend(0.2,0.1,0.4,0.3);
  leg->SetTextSize(labelsize);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(fHistoAntiPtSpectrum, "^{3}#bar{H}","l");
  leg->AddEntry(fHistoPtSpectrum, "^{3}H","l");
  if (MultiDep) {
    leg->AddEntry((TObject*)0,Form("Multiplicity: (%.0f - %.0f) %%", MultiplicityPercentile[0], MultiplicityPercentile[1]),"");
  }

  fHistoAntiPtSpectrum->SetMarkerColor(kRed);
  fHistoAntiPtSpectrum->SetLineColor(kRed);
  fHistoPtSpectrum->SetMarkerColor(kBlue);
  fHistoPtSpectrum->SetLineColor(kBlue);

  fHistoAntiPtSpectrum->SetMarkerStyle(20);
  fHistoAntiPtSpectrum->SetMarkerSize(0.5);
  fHistoPtSpectrum->SetMarkerStyle(24);
  fHistoPtSpectrum->SetMarkerSize(0.5);
  fHistoPtSpectrum->GetXaxis()->SetRangeUser(1.5,3.);
  
  hAxis->GetYaxis()->SetRangeUser(9e-8,2e-6);
  hAxis->DrawCopy();
  fHistoPtSpectrum -> DrawCopy("EX0same");
  fHistoAntiPtSpectrum -> DrawCopy("EX0same");

  //Systematic uncertainties
  // Bin											  1.25	1.75	 2.25		2.75
  Double_t AntiTritonSyst[] = {0.18,	0.12,	 0.20,	0.32};
  Double_t TritonSyst[] =     {0.23,	0.09,	 0.14,	0.20};

  TH1D* fHistoPtSpectrumSyst = (TH1D*) fHistoPtSpectrum->Clone("fHistoPtSpectrumSyst");
  for (int i=1; i < = fHistoPtSpectrumSyst->GetNbinsX(); i++) {
    fHistoPtSpectrumSyst->SetBinError(i, fHistoPtSpectrumSyst->GetBinContent(i)*TritonSyst[i-1]);
  }
  TH1D* fHistoAntiPtSpectrumSyst = (TH1D*) fHistoAntiPtSpectrum->Clone("fHistoAntiPtSpectrumSyst");
  for (int i=1; i < = fHistoAntiPtSpectrumSyst->GetNbinsX(); i++) {
    fHistoAntiPtSpectrumSyst->SetBinError(i, fHistoAntiPtSpectrumSyst->GetBinContent(i)*AntiTritonSyst[i-1]);
  }

  TH1D* fHistoRatioSpectraSyst = (TH1D*) fHistoAntiPtSpectrumSyst->Clone("fHistoRatioSpectraSyst");
  fHistoRatioSpectraSyst -> Divide(fHistoAntiPtSpectrumSyst,fHistoPtSpectrumSyst);
  
  //  correct ratio systematic uncertainties (without correlated part)
  Double_t RatioSystTriton[8] = {0.26, 0.12, 0.23, 0.37};
  Double_t RatioError = 0;
  for (int i=1; i <= fHistoRatioSpectraSyst->GetNbinsX(); i++) {
    if (fHistoRatioSpectraSyst->GetBinCenter(i)<3) continue;
    RatioError = RatioSystTriton[i-1]*fHistoRatioSpectraSyst->GetBinContent(i);
    fHistoRatioSpectraSyst->SetBinError(i, RatioError);
  }

  fHistoPtSpectrumSyst -> SetFillStyle(0);
  fHistoAntiPtSpectrumSyst -> SetFillStyle(0);
  fHistoRatioSpectraSyst -> SetFillStyle(0);

  fHistoPtSpectrumSyst -> DrawCopy("E2same");
  fHistoAntiPtSpectrumSyst -> DrawCopy("E2same");

  TLegend* LegALICE = new TLegend(0.65,0.75,0.9,0.9);
  LegALICE->SetTextSize(labelsize);
  LegALICE->SetFillColor(0);
  LegALICE->SetFillStyle(0);
  LegALICE->SetBorderSize(0);
  LegALICE->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
  LegALICE->AddEntry((TObject*)0,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
  LegALICE->Draw();
  leg->Draw("same");
  gPad->Update();
  cSpectra->cd(2);
  gStyle->SetPadTopMargin(0.);
  gStyle->SetPadBottomMargin(0.25);
  gPad->SetTicks();
  fHistoRatioSpectra->SetMarkerColor(kRed);
  fHistoRatioSpectra->SetLineColor(kRed);
  fHistoRatioSpectra->SetMarkerStyle(20);
  fHistoRatioSpectra->SetMarkerSize(0.5);

  hAxis->GetYaxis()->SetTitle(RatioLabel);
  hAxis->GetYaxis()->SetNdivisions(505);
  hAxis->GetYaxis()->SetRangeUser(0,1.9);
  hAxis->GetYaxis()->CenterTitle();
  hAxis->GetYaxis()->SetLabelOffset(0.02);
  hAxis->DrawCopy();

  fHistoRatioSpectra -> DrawCopy("EX0same");
  fHistoRatioSpectraSyst -> DrawCopy("E2same");
  TLine *line = new TLine(PlotMinPt,1,PlotMaxPt,1);
  line->SetLineStyle(2);
  line->Draw("same");
  cSpectra->SaveAs("SpectraAndRatio.pdf");

  // default
  TFile *outputFile = new TFile (Form("%s",fOutputFileName),"RECREATE");

  // Save result for evaluation of systematics
  //   TFile *outputFile = new TFile (Form("../Results/Triton/TPCTOF/Systematics/Secondaries/Variations/Triton_LHC16q_FAST_woSDD_DCAxy%.2fz%.2f.root",DCAcutXY,DCAcutZ),"RECREATE");

  //   TFile *outputFile = new TFile (Form("../Results/Triton/TPCTOF/Systematics/Tracking/Variations/Triton_LHC16q_FAST_woSDD_ITS%dTPC%d_%.2fy%.2f.root",ITShits, TPCcluster, minRapidity, maxRapidity),"RECREATE");

  // TPC PID
  //   TFile *outputFile = new TFile (Form("../Results/Triton/TPCTOF/Systematics/PID/Variations/Triton_LHC16q_FAST_woSDD_TPC_%.1f_%.1f.root", TPCselectMin, TPCselectMax),"RECREATE");
  // TOF PID
  //   TFile *outputFile = new TFile (Form("../Results/Triton/TPCTOF/Systematics/PID/Variations/Triton_LHC16q_FAST_woSDD_TOF_%.1f.root", TOFselect),"RECREATE");

  outputFile -> cd();
  fHistoPtSpectrum -> Write();
  fHistoAntiPtSpectrum -> Write();
  fHistoRatioSpectra -> Write();
  //systematic uncertainties
  fHistoPtSpectrumSyst -> Write();
  fHistoAntiPtSpectrumSyst -> Write();
  fHistoRatioSpectraSyst -> Write();
  leg->Write();
  outputFile ->Close();
}
//___________________________________________________________________________________________
Double_t CalculateRapidity(Double_t y){
  return  y - 0.465;
}

//___________________________________________________________________________________________
Double_t CorrectPt(Double_t px, Double_t py){
  return TMath::Sqrt(px*px + py*py);
}
//___________________________________________________________________________________________
Double_t BackgroundFitWithoutSignalRegion(Double_t *x, Double_t *par){
  if (TMath::Abs(x[0]) < 0.1) {
    TF1::RejectPoint();
    return 0;
  }
  if (SecondaryGaus) {
    return par[0]*TMath::Gaus(x[0], par[1], par[2]); // gaus
  } else{
    return par[0] + par[1]*x[0]+ par[2]*x[0]*x[0]; // pol2
  }
}
