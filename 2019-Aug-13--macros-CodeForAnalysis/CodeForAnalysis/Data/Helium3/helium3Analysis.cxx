#include "helium3Analysis.h"
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
TString RatioLabel = "^{3}#bar{He} / ^{3}He";
Double_t PlotMinPt = 1;
Double_t PlotMaxPt = 5.5;

Double_t DCAcutXY;
Double_t DCAcutZ;
Double_t minRapidity;
Double_t maxRapidity;
Int_t ITShits;
Int_t TPCcluster;
Double_t TPCselectMin;
Double_t TPCselectMax;

Bool_t SecondaryGaus = kTRUE;

// Multiplicity
Bool_t MultiDep = false; // Switch for multiplicity dependent analysis (true = multiplicity depenent , false = multiplicity integrated)
Double_t MultiplicityPercentile[2] = {40,100}; // Set the boundaries of the multiplicity percentile of interest

// Binning for the multiplicity integrated result
const Int_t nPtBins=8;
Double_t Binning[]={1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6.,7.};

//// Binning for 0-10 % multiplicity
//const Int_t nPtBins = 6;
//Double_t Binning[]={1., 1.5, 2., 2.5, 3., 3.5, 5.};

//// Binning for 10-20 % and 20-40 % multiplicity
//const Int_t nPtBins = 5;
//Double_t Binning[]={1., 1.5, 2., 2.5, 3., 5.};

//// Binning for 40-100 % multiplicity
//const Int_t nPtBins = 4;
//Double_t Binning[]={1., 1.5, 2., 3., 5.};

//___________________________________________________________________________________________
Double_t BackgroundFitWithoutSignalRegion(Double_t *x, Double_t *par)
{
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

//___________________________________________________________________________________________
helium3Analysis::helium3Analysis():
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
fHistoHe3Mean(NULL),
fHistoHe3Width(NULL),
fHistoHe3PtSpectrum(NULL),
fHistoAntiHe3PtSpectrum(NULL),
fHistoRatioSpectra(NULL)
{
  gStyle->SetOptStat(0);         // switch off stat box
  gStyle->SetOptFit(1111);         // switch on fit box
}

//___________________________________________________________________________________________
helium3Analysis::~helium3Analysis()
{
  delete tree;
  delete PtCorr;
  delete fHistoRapidity;
  delete fHistoHe3Mean;
  delete fHistoHe3Width;
  delete fHistoHe3PtSpectrum;
  delete fHistoAntiHe3PtSpectrum;
  delete fHistoRatioSpectra;
}

//___________________________________________________________________________________________
void helium3Analysis::GetInputTree()  {

  TFile *inputFile = TFile::Open (Form("%s",fInputFileName));
  TList *inputList = (TList*) inputFile->Get(fInputListName);

  tree = (TTree*) inputFile->GetObjectChecked("Trees/reducedTree_Helium","TTree"); // tree is stored in a subfolder called Trees

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
void helium3Analysis::SetTreeBranches()  {

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
  tree -> SetBranchAddress("nSigmaITS_He4",&nSigmaITS_He4);
  tree -> SetBranchAddress("nSigmaTPC_He4",&nSigmaTPC_He4);
  tree -> SetBranchAddress("nSigmaTOF_He4",&nSigmaTOF_He4);
  tree -> SetBranchAddress("nSigmaTRD_He4",&nSigmaTRD_He4);
  tree -> SetBranchAddress("nSigmaHMPID_He4",&nSigmaHMPID_He4);
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
  //   tree -> SetBranchAddress("nSigmaITS_Deut",&nSigmaITS_Deut);
  //   tree -> SetBranchAddress("nSigmaTPC_Deut",&nSigmaTPC_Deut);
  //   tree -> SetBranchAddress("nSigmaTOF_Deut",&nSigmaTOF_Deut);
  //   tree -> SetBranchAddress("nSigmaTRD_Deut",&nSigmaTRD_Deut);
  //   tree -> SetBranchAddress("nSigmaHMPID_Deut",&nSigmaHMPID_Deut);
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
void helium3Analysis::CreateHistograms()  {

  fHistoHe3PtSpectrum = new TH1D ("fHistoHe3PtSpectrum","",nPtBins,Binning);
  fHistoHe3PtSpectrum->GetXaxis()->SetTitle(pTLabel);
  fHistoHe3PtSpectrum->GetXaxis()->SetTitleSize(labelsize+0.01);
  fHistoHe3PtSpectrum->GetXaxis()->SetLabelSize(labelsize);
  fHistoHe3PtSpectrum->GetYaxis()->SetTitle(SpectraLabel);
  fHistoHe3PtSpectrum->GetYaxis()->SetTitleSize(labelsize+0.01);
  fHistoHe3PtSpectrum->GetYaxis()->SetTitleOffset(0.8);
  fHistoHe3PtSpectrum->GetYaxis()->SetLabelSize(labelsize);

  fHistoAntiHe3PtSpectrum = new TH1D ("fHistoAntiHe3PtSpectrum","",nPtBins,Binning);

  fHistoRatioSpectra = new TH1D ("fHistoRatioSpectra","",nPtBins,Binning);
  fHistoRatioSpectra->GetXaxis()->SetTitle(pTLabel);
  fHistoRatioSpectra->GetXaxis()->SetTitleSize(labelsize+0.01);
  fHistoRatioSpectra->GetXaxis()->SetLabelSize(labelsize);
  fHistoRatioSpectra->GetYaxis()->SetTitle(RatioLabel);
  fHistoRatioSpectra->GetYaxis()->SetTitleSize(labelsize+0.01);
  fHistoRatioSpectra->GetYaxis()->SetTitleOffset(0.8);
  fHistoRatioSpectra->GetYaxis()->SetLabelSize(labelsize);

  fHistoHe3PtSpectrum->Sumw2();
  fHistoAntiHe3PtSpectrum->Sumw2();
  fHistoRatioSpectra->Sumw2();

  fHistoRapidity = new TH1D ("fHistoRapidity","Rapidity",30,-1.5,1.5);
  fHistoHe3Mean = new TH1D ("fHistoHe3Mean","Mean of the helium distribution",31,-3.,3.);
  fHistoHe3Width = new TH1D ("fHistoHe3Width","Width of the helium distribution",15,0.,3.);

}
//___________________________________________________________________________________________
void helium3Analysis::ExecuteAnalysis()  {

  // defining size of pad ---------------------------
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);

  SaveToFile = false; // Save efficiency and pT slice fits to pdf
  Bool_t DrawRapidity = false;
  Bool_t DrawTPCselection = false; // possibility to display the TPC selection in a plot
  Bool_t DrawDCADistribution = false; // Draw DCA distributions (DCAxy is also drawn to estimate the primary fraction)
  Bool_t DrawTOF = false;

  Bool_t EvaluatePrimaryFraction = false; // Estimate primary fraction via fit of the DCA disribution
  if (MultiDep) EvaluatePrimaryFraction = false; // Turn off estimation for multiplicity intervals because of too low statitics (use values from integrated sample)
  Bool_t UseFitter = kTRUE; //Primary fraction estimated using TFractionFitter, if false (1 - SecondaryFit/Total) is used
  Bool_t GausFit = kTRUE; // Triton fit using a Gaussian if true, Exponential function if false
  Bool_t UseMeanAndSigmaFromData = kTRUE;
  Bool_t SubtractHypertritonFeedDown = kTRUE;

  // Cuts
  Double_t SelectTPCPlot = 5; // Selection window for the TPC dE/dx plot used to fit and remove triton
  Double_t minFit = -5.; // Lower boundary for the triton fit interval
  Double_t maxFit = -3.; // Upper boundary for the triton fit interval
  Double_t maxEta = 0.9; // Pseudorapidity range

  // Parameters varied to evaluate the systematic uncertainties
  minRapidity = -1.;
  maxRapidity = 0.;
  ITShits = 2;
  TPCcluster = 70;

  DCAcutXY = 0.1;
  DCAcutZ = 1.; // DCAz cut in filter bit 3.2cm

  Double_t TPCselect = 3.0; // TPC selection window
  TPCselectMin = -3.; // TPC PID (for systematics)
  TPCselectMax = 3.; // TPC PID (for systematics)

//  // Effect of the DCA cuts
//  TCanvas* cDCA2D[nPtBins];
//  TH2D* hDCA2D;
//  Double_t Ptbin = 0;
//  for (Int_t i =0; i < nPtBins; i++) {
//    Ptbin = (Binning[i]+Binning[i+1])/2;
//    cDCA2D[i] = new TCanvas(Form("cDCA2D_%d",i), Form("DCA 2D %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
//    gPad->SetTicks();
//
//    tree->Draw("dcaz:dcaxy>>DCA2D(48,-2.4,2.4,64,-3.2,3.2)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && nSigmaTPC_He3 > %f && nSigmaTPC_He3 < %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f ", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, TPCselectMin, TPCselectMax, Binning[i], Binning[i+1]));
//
//    hDCA2D = (TH2D*) gDirectory->Get("DCA2D");
//
//    hDCA2D->SetTitle("");
//    hDCA2D->GetXaxis()->SetTitle("DCA_{xy} (cm)");
//    hDCA2D->GetXaxis()->SetTitleSize(legsize+0.01);
//    hDCA2D->GetXaxis()->SetLabelSize(legsize);
//    hDCA2D->GetYaxis()->SetTitleSize(legsize+0.01);
//    hDCA2D->GetYaxis()->SetLabelSize(legsize);
//    hDCA2D->GetYaxis()->SetTitle("DCA_{z} (cm)");
//
//    hDCA2D->DrawCopy("");
//
//    cout << "Cut effect: " << 1. - hDCA2D->Integral(hDCA2D->GetXaxis()->FindBin(-0.1), hDCA2D->GetXaxis()->FindBin(0.1), hDCA2D->GetYaxis()->FindBin(-1.), hDCA2D->GetYaxis()->FindBin(1.)) / hDCA2D->Integral(hDCA2D->GetXaxis()->FindBin(-2.4), hDCA2D->GetXaxis()->FindBin(2.4), hDCA2D->GetYaxis()->FindBin(-3.2), hDCA2D->GetYaxis()->FindBin(3.2)) << endl;
//  }

  if (DrawRapidity) {
    // Evaluate the rapidity distributrion
    TCanvas *cRapidity = new TCanvas("cRapidity","Rapidity", 800,600);
    TLegend* LegALICERap = new TLegend(0.55,0.8,0.95,0.925);
    LegALICERap->SetTextSize(legsize);
    LegALICERap->SetFillColor(0);
    LegALICERap->SetFillStyle(0);
    LegALICERap->SetBorderSize(0);
    LegALICERap->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
    LegALICERap->AddEntry((TObject*)0,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");

    if (MultiDep){
      tree->Draw("CalculateRapidity(px,py,pz)>>fHistoRapidity(30,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(nSigmaTPC_He3) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && multPercentile_V0A >= %f && multPercentile_V0A < %f", maxEta, ITShits, TPCcluster, TPCselect, DCAcutZ, DCAcutXY,MultiplicityPercentile[0], MultiplicityPercentile[1]));
    } else {
      tree->Draw("CalculateRapidity(px,py,pz)>>fHistoRapidity(30,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(nSigmaTPC_He3) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f", maxEta, ITShits, TPCcluster, TPCselect, DCAcutZ, DCAcutXY));
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
  TFile *EffFile = TFile::Open ("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Efficiency.root");

  // For the evaluation of the systematic uncertainties due to the PID selection
  //      TFile *EffFile = TFile::Open (Form("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Systematics/PID/Efficiency_TPC_%.1f_%.1f.root",TPCselectMin,TPCselectMax));

  // For the evaluation of the systematic uncertainties due tracking
//  TFile *EffFile = TFile::Open (Form("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Systematics/Tracking/Efficiency_ITS%dTPC%d_%.2fy%.2f.root",ITShits, TPCcluster, minRapidity, maxRapidity));

  //   For the evaluation of the systematic uncertainties due the secondary contamination via DCA variation
  //      TFile *EffFile = TFile::Open (Form("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Systematics/DCA(Secondaries)/Efficiency_DCAxy%.2fz%.2f.root",DCAcutXY,DCAcutZ));

  TH1D* hEffFAST = (TH1D*) EffFile->GetObjectChecked("hEffHelium", "TH1D");
  TH1D* hEffFASTAnti = (TH1D*) EffFile->GetObjectChecked("hEffAntiHelium", "TH1D");

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
  TLegend* LegALICE = new TLegend(0.55,0.8,0.95,0.925);
  LegALICE->SetTextSize(legsize);
  LegALICE->SetFillColor(0);
  LegALICE->SetFillStyle(0);
  LegALICE->SetBorderSize(0);
  LegALICE->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
  LegALICE->AddEntry((TObject*)0,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");

  gPad->SetTicks();
  TLegend* EffLeg = new TLegend(0.5,0.2,0.9,0.4,"ALICE #it{Work in progress}");
  EffLeg->SetTextSize(legsize);
  EffLeg->SetFillColor(0);
  EffLeg->SetFillStyle(0);
  EffLeg->SetBorderSize(0);
  EffLeg->AddEntry((TObject*)0,"#kern[-0.35]{p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV}","");
  EffLeg->AddEntry(hEffFAST, "#varepsilon_{^{3}He}","l");
  EffLeg->AddEntry(hEffFASTAnti, "#varepsilon_{Anti-^{3}He}","l");

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

  TH1D* hEffAnti = (TH1D*) hEffFASTAnti->Clone("hEffAnti");
  TH1D* hEff = (TH1D*) hEffFAST->Clone("hEff");
  if (SaveToFile) cEff->SaveAs("Efficiency.pdf");

  TCanvas* cEffRatio = new TCanvas("cEffRatio","Ratio of efficiencies");
  gPad->SetTicks();
  TH1D *hEffRatio = (TH1D*) hEffAnti->Clone("hEffRatio");
  hEffRatio->GetYaxis()->SetTitleOffset(0.8);
  hEffRatio->GetYaxis()->SetTitle("#varepsilon_{Anti-^{3}He} / #varepsilon_{^{3}He}");
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

    TCanvas* cTPC3He = new TCanvas("cTPC3He","TPC Selection ^{3}He");
    gPad->SetTicks();
    tree->Draw("nSigmaTPC_He3:CorrectPt(px, py)>>HeBuffer(100,0,10,100,-5,5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && CorrectPt(px, py) < 10 && TMath::Abs(nSigmaTPC_He3) < 5", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster));
    TPCBuffer = (TH2D*) gDirectory->Get("HeBuffer");
    TPCBuffer->SetTitle("");
    TPCBuffer->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TPCBuffer->GetYaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}He}}{#sigma_{d#it{E} / d#it{x}}}");
    TPCBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetXaxis()->SetLabelSize(legsize);
    TPCBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetYaxis()->SetLabelSize(legsize);
    TPCBuffer->DrawCopy("");
    UpperBoundary.DrawClone("same");
    LowerBoundary.DrawClone("same");
    LegALICE->Draw();

    TCanvas* cTPCAnti3He = new TCanvas("cTPCAnti3He","TPC Selection ^{3}#bar{He}");
    gPad->SetTicks();
    tree->Draw("nSigmaTPC_He3:CorrectPt(px, py)>>AntiHeBuffer(100,0,10,100,-5,5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && CorrectPt(px, py) < 10 && TMath::Abs(nSigmaTPC_He3) < 5", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster));
    TPCBuffer = (TH2D*) gDirectory->Get("AntiHeBuffer");
    TPCBuffer->SetTitle("");
    TPCBuffer->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TPCBuffer->GetYaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}He}}{#sigma_{d#it{E} / d#it{x}}}");
    TPCBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetXaxis()->SetLabelSize(legsize);
    TPCBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
    TPCBuffer->GetYaxis()->SetLabelSize(legsize);
    TPCBuffer->DrawCopy("");
    UpperBoundary.DrawClone("same");
    LowerBoundary.DrawClone("same");
    LegALICE->Draw();

    cTPC3He->SaveAs("TPCselection.pdf");
    cTPCAnti3He->SaveAs("TPCselection_Anti.pdf");
  }

  //Check DCA_XY
  TCanvas* cDCAxy[nPtBins];
  TCanvas* cDCAz[nPtBins];
  Double_t PrimaryFraction[nPtBins];
  for (Int_t j=0; j < nPtBins; j++){
    PrimaryFraction[j] = 1;
    if(j==0) PrimaryFraction[j] = 0.43;
    if(j==1) PrimaryFraction[j] = 0.73;
    if(j==2) PrimaryFraction[j] = 0.945;
  }
  Double_t PrimaryFractionUncertainty;

  Double_t Ptbin = 0;
  if (DrawDCADistribution || EvaluatePrimaryFraction) {
    TH1D* HeBuffer;
    TH1D* AntiHeBuffer;

    // File with templates from MC
    TFile *PrimaryTemplateFile = TFile::Open ("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/DCATemplates_All.root");
    //DCA output file
    TFile *DCAOutputFile = new TFile ("Helium_DCAxy.root","RECREATE");
    DCAOutputFile->cd();

    // Buffer for the DCA templates
    TH1D* DCATemplatePrimary;
    TH1D* DCATemplateSecondary;
    TObjArray *MCDCATemplates = new TObjArray(2);        // MC histograms are put in this array
    TFractionFitter* DCAFitter;
    // Buffer for template fit results
    Int_t FitStatus = -1;
    Double_t Secondary = 0., Primary = 1., ErrorPrimary = 0., ErrorSecondary = 0.;

    TF1* SecondaryWoPeak = new TF1("SecondaryWoPeak", BackgroundFitWithoutSignalRegion,-6,6,3);
    
    TF1* SecondaryHe;
    if (SecondaryGaus) {
      SecondaryWoPeak->SetParameters(100,0,20);
      SecondaryHe = new TF1("SecondaryHe","gaus",-6,6);
    } else {
      SecondaryHe = new TF1("SecondaryHe","pol2",-6,6);
    }

    for (Int_t i =0; i < nPtBins; i++) {
      Ptbin = (Binning[i]+Binning[i+1])/2;
      if (Ptbin > 2.5) continue; // change to 3 GeV for DCA variation
      cDCAxy[i] = new TCanvas(Form("DCAxy_%d",i), Form("DCA %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
      gPad->SetTicks();
//      gPad->SetLogy();
      Double_t LowerTPCforDCA = -2.0;
      if (Ptbin > 2.0) LowerTPCforDCA = -2.5;

      if (MultiDep) {
        tree->Draw("dcaxy>>HeBuffer(40,-1.0,1.0)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && nSigmaTPC_He3 > %f && nSigmaTPC_He3 < %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && multPercentile_V0A >= %f && multPercentile_V0A < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], DCAcutZ, MultiplicityPercentile[0], MultiplicityPercentile[1]));
        // Anti-Particle
        tree->Draw("dcaxy>>AntiHeBuffer(40,-1.0,1.0)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 &&  nSigmaTPC_He3 > %f && nSigmaTPC_He3 < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && multPercentile_V0A >= %f && multPercentile_V0A < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], DCAcutZ, MultiplicityPercentile[0], MultiplicityPercentile[1]));
      }else{
        // Particle
        tree->Draw("dcaxy>>HeBuffer(40,-1.0,1.0)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && nSigmaTPC_He3 > %f && nSigmaTPC_He3 < %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], DCAcutZ));
        // Anti-Particle
        tree->Draw("dcaxy>>AntiHeBuffer(40,-1.0,1.0)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 &&  nSigmaTPC_He3 > %f && nSigmaTPC_He3 < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1], DCAcutZ));
      }
      HeBuffer = (TH1D*) gDirectory->Get("HeBuffer");
      AntiHeBuffer = (TH1D*) gDirectory->Get("AntiHeBuffer");

      HeBuffer->SetTitle("");
      HeBuffer->GetXaxis()->SetTitle("DCA_{xy} (cm)");
      HeBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
      HeBuffer->GetXaxis()->SetLabelSize(legsize);
      HeBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
      HeBuffer->GetYaxis()->SetLabelSize(legsize);
      HeBuffer->GetYaxis()->SetTitle("Counts");

      AntiHeBuffer->SetTitle("");
      AntiHeBuffer->GetXaxis()->SetTitle("DCA_{xy} (cm)");
      AntiHeBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
      AntiHeBuffer->GetXaxis()->SetLabelSize(legsize);
      AntiHeBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
      AntiHeBuffer->GetYaxis()->SetLabelSize(legsize);
      AntiHeBuffer->GetYaxis()->SetTitle("Counts");

      AntiHeBuffer->SetMarkerColor(kRed);
      AntiHeBuffer->SetLineColor(kRed);
      HeBuffer->SetMarkerColor(kBlue);
      HeBuffer->SetLineColor(kBlue);
      AntiHeBuffer->SetMarkerStyle(20);
      AntiHeBuffer->SetMarkerSize(0.5);
      HeBuffer->SetMarkerStyle(24);
      HeBuffer->SetMarkerSize(0.5);
      HeBuffer->GetYaxis()->SetRangeUser(0, 80);
      HeBuffer->DrawCopy("");
      AntiHeBuffer->DrawCopy("same");

      if (EvaluatePrimaryFraction) {
        TFitResultPtr FitSecBkg = HeBuffer->Fit(SecondaryWoPeak, "QNL", "", -3, 3);
        SecondaryWoPeak->SetLineColor(kBlack);
        SecondaryWoPeak->DrawCopy("same");
        SecondaryHe->SetParameters(SecondaryWoPeak->GetParameters());
        SecondaryHe->DrawCopy("same");

        //Load templates (MC) -> d or t used -> pT_d = 2* pT_He
        DCATemplatePrimary = (TH1D*) PrimaryTemplateFile->GetObjectChecked(Form("hDCAxy_Primary_%.2f", Ptbin), "TH1D");
        DCATemplatePrimary->Rebin(4); // Adjust binning
//        DCATemplateSecondary = (TH1D*) PrimaryTemplateFile->GetObjectChecked(Form("hDCAxy_Triton_Material_%.3f", Ptbin/2.), "TH1D");
        DCATemplateSecondary = (TH1D*) PrimaryTemplateFile->GetObjectChecked(Form("hDCAxy_Deuteron_Material_%.3f", Ptbin/2.), "TH1D");
        DCATemplateSecondary->Rebin(4); // Adjust binning

//        // Data driven approach
//        // Use Anti-He as proxy for the primary He
//        DCATemplatePrimary = (TH1D*) AntiHeBuffer->Clone(Form("hDCAxyPrimary_%.1f", Ptbin));
//
//        // Use fit to the DCAxy distribution as proxy for the secondary distribution (Functional form to be varied for syst.)
//        DCATemplateSecondary = (TH1D*) HeBuffer->Clone(Form("hDCAxyMaterial_%.1f", Ptbin));
//        DCATemplateSecondary->TH1::Reset();
//        DCATemplateSecondary->Add(SecondaryHe);
//        for (Int_t jbin=0; jbin<DCATemplateSecondary->GetNbinsX(); jbin++) {
//          DCATemplateSecondary->SetBinError(jbin,TMath::Sqrt(DCATemplateSecondary->GetBinContent(jbin)));
//        }
//        DCATemplatePrimary->Sumw2();
//        DCATemplateSecondary->Sumw2();

        if (UseFitter) {
          //Empty TObjArray
          MCDCATemplates->Clear();
          //Fill TObjArray
          MCDCATemplates->Add(DCATemplatePrimary);
          MCDCATemplates->Add(DCATemplateSecondary);
          //Fit templates
          DCAFitter = new TFractionFitter(HeBuffer, MCDCATemplates,"Q");
          DCAFitter->Constrain(1,0.0,1.0);
          DCAFitter->Constrain(2,0.0,1.0);
          TVirtualFitter::SetMaxIterations(10000); // Restrict the number of iterations
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
            Float_t DataIntegral = result->Integral(); // Integral of the total template fit
            hPrim->Scale(DataIntegral * Primary / hPrim->Integral());
            hSec->Scale(DataIntegral * Secondary / hSec->Integral());
            hSec->SetLineColor(kGreen-3);
            hPrim->SetLineColor(kOrange-3);
            hSec->SetMarkerColor(kGreen-3);
            hPrim->SetMarkerColor(kOrange-3);
            hSec->DrawCopy("EX0same");
            hPrim->DrawCopy("EX0same");
            PrimaryFraction[i] = hPrim->Integral(hPrim->FindBin(-(DCAcutXY-0.001)), hPrim->FindBin((DCAcutXY-0.001))) / result->Integral(result->FindBin(-(DCAcutXY-0.001)), result->FindBin((DCAcutXY-0.001)));

            Double_t relativeFracionError = ErrorPrimary/Primary;
            Double_t IntegralError;
            Primary = hPrim->IntegralAndError(hPrim->FindBin(-(DCAcutXY-0.001)), hPrim->FindBin((DCAcutXY-0.001)),IntegralError);
            IntegralError = IntegralError/Primary;
            ErrorPrimary = Primary*sqrt(pow(relativeFracionError,2)+pow(IntegralError,2));

            Double_t relativeFracionErrorSec = ErrorSecondary/Secondary;
            IntegralError=0;
            Secondary = hSec->IntegralAndError(hSec->FindBin(-(DCAcutXY-0.001)), hSec->FindBin((DCAcutXY-0.001)),IntegralError);
            IntegralError = IntegralError/Secondary;
            ErrorSecondary = Secondary*sqrt(pow(relativeFracionErrorSec,2)+pow(IntegralError,2));

            PrimaryFractionUncertainty = abs(ErrorPrimary*Secondary-Primary*ErrorSecondary)/pow(Primary+Secondary,2);

            HeBuffer->GetYaxis()->UnZoom();
            HeBuffer->Write(Form("HeliumDCAxy%.2fGeV",Ptbin));
            AntiHeBuffer->Write(Form("AntiHeliumDCAxy%.2fGeV",Ptbin));
            hPrim->Write(Form("PrimaryTemplate%.2fGeV",Ptbin));
            hSec->Write(Form("SecondaryTemplate%.2fGeV",Ptbin));
            result->Write(Form("TemplateFit%.2fGeV",Ptbin));
          }
        } else {
          PrimaryFraction[i] = 1.-(SecondaryHe->Integral(-DCAcutXY,DCAcutXY)/HeBuffer->Integral(HeBuffer->FindBin(-(DCAcutXY-0.001)),HeBuffer->FindBin((DCAcutXY-0.001)),"width"));
          PrimaryFractionUncertainty = SecondaryHe->IntegralError(-DCAcutXY,DCAcutXY) / HeBuffer->Integral(HeBuffer->FindBin(-(DCAcutXY-0.01)),HeBuffer->FindBin((DCAcutXY-0.01)),"width"); //Dummy, not complet
        }


        printf("Primary fraction:  %.3f ±  %.3f\n", PrimaryFraction[i],PrimaryFractionUncertainty);
        printf("Anti-He fraction: %.2f \n", AntiHeBuffer->Integral(AntiHeBuffer->FindBin(-DCAcutXY), AntiHeBuffer->FindBin(DCAcutXY)) / HeBuffer->Integral(HeBuffer->FindBin(-DCAcutXY), HeBuffer->FindBin(DCAcutXY)));
      }
      LegALICE->Draw();
      if (DrawDCADistribution) cDCAxy[i]->SaveAs(Form("DCAxy%.2fGeV.pdf",Ptbin));

      if (DrawDCADistribution){
        cDCAz[i] = new TCanvas(Form("DCAz_%d",i), Form("DCA %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
        gPad->SetTicks();
//        gPad->SetLogy();
        // Particle
        tree->Draw("dcaz>>HeBuffer(60,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && nSigmaTPC_He3 > %f && nSigmaTPC_He3 < %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1]));
        // Anti-Particle
        tree->Draw("dcaz>>AntiHeBuffer(60,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && nSigmaTPC_He3 > %f && nSigmaTPC_He3 < %f  &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, LowerTPCforDCA, TPCselectMax, Binning[i], Binning[i+1]));

        HeBuffer = (TH1D*) gDirectory->Get("HeBuffer");
        AntiHeBuffer = (TH1D*) gDirectory->Get("AntiHeBuffer");
        HeBuffer->SetTitle("");
        HeBuffer->GetXaxis()->SetTitle("DCA_{z} (cm)");
        HeBuffer->GetXaxis()->SetTitleSize(legsize+0.01);
        HeBuffer->GetXaxis()->SetLabelSize(legsize);
        HeBuffer->GetYaxis()->SetTitleSize(legsize+0.01);
        HeBuffer->GetYaxis()->SetLabelSize(legsize);
        AntiHeBuffer->SetMarkerColor(kRed);
        AntiHeBuffer->SetLineColor(kRed);
        HeBuffer->SetMarkerColor(kBlue);
        HeBuffer->SetLineColor(kBlue);
        AntiHeBuffer->SetMarkerStyle(20);
        AntiHeBuffer->SetMarkerSize(0.5);
        HeBuffer->SetMarkerStyle(24);
        HeBuffer->SetMarkerSize(0.5);
        HeBuffer->GetYaxis()->SetRangeUser(0, 80);
        HeBuffer->DrawCopy("");
        AntiHeBuffer->DrawCopy("same");
        LegALICE->Draw();
        cDCAz[i]->SaveAs(Form("DCAz%.2fGeV.pdf",Ptbin));
      }
    }
    DCAOutputFile->Close();
  }

  TF1* Triton;
  if (GausFit) {
    //Gaussian function
    Triton = new TF1("Triton","gaus",-6,6);
    Triton->SetParameter(0, 10); // constant
    Triton->SetParameter(1, -5); // mean
    Triton->SetParameter(1, 2); // width
    Triton->SetLineColor(kRed);
  } else {
    //Exponential function
    Triton = new TF1("Triton","[0]*TMath::Exp(x*[1])",-10,10);
    Triton->SetParameter(0, -1); // constant
    Triton->SetParameter(1, -1.5); // slope
    Triton->SetParLimits(1, -1e3, -0.2); // Limits for the slope parameter
    Triton->SetLineColor(kRed);
    Triton->SetParNames("const","sigma");
    minFit = -4.4; // Lower boundary for the triton fit interval reduced to avoid flattening at the mean of the distribution
  }

  TF1* Gaus = new TF1("Gaus","gaus",-3,3);
  TCanvas* cFit[nPtBins];
  TCanvas* cFitAnti[nPtBins];
  TCanvas* cTOF[nPtBins];
  TH1D *hTempTPC[nPtBins];
  TH1D *hTempTritonSubtracted[nPtBins];
  Int_t AntiHeCounter = 0;

  // Helium
  for (Int_t i =0; i < nPtBins; i++) {
    Ptbin = (Binning[i]+Binning[i+1])/2;
    Double_t mean = 0;
    Double_t sigma = 1;

    if(DrawTOF){
      cTOF[i] = new TCanvas();
      cTOF[i]->SetTitle(Form("TOF %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
      cTOF[i]->Divide(2,1);
      cTOF[i]->cd(1);
      gPad->SetTicks();
      tree->Draw("nSigmaTOF_He3", Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && TMath::Abs(nSigmaTPC_He3) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f ", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, TPCselect, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY));
      cTOF[i]->cd(2);
      gPad->SetTicks();
      tree->Draw("nSigmaTOF_He3", Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && TMath::Abs(nSigmaTPC_He3) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f ", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, TPCselect, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY));
    }

    // evaluation of (primary) 3He
    cFit[i] = new TCanvas();
    cFit[i]->SetTitle(Form("Fit %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
    cFit[i]->Divide(2,1);
    cFit[i]->cd(1);
    gPad->SetTicks();
    gPad->SetTopMargin(0.07);

    if (MultiDep){
      tree->Draw(Form("nSigmaTPC_He3>>Bin%d(%d,%f,%f)",i,(int) SelectTPCPlot*10,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && TMath::Abs(nSigmaTPC_He3) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && multPercentile_V0A >= %f && multPercentile_V0A < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY,MultiplicityPercentile[0], MultiplicityPercentile[1]));
    } else {
      tree->Draw(Form("nSigmaTPC_He3>>Bin%d(%d,%f,%f)",i,(int) SelectTPCPlot*10,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && TMath::Abs(nSigmaTPC_He3) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY));
    }
    hTempTPC[i] = (TH1D*) gDirectory->Get(Form("Bin%d",i));
    hTempTPC[i]->GetXaxis()->SetRangeUser(-SelectTPCPlot, SelectTPCPlot);
    hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}He}}{#sigma_{d#it{E} / d#it{x}}}");
    hTempTPC[i]->GetXaxis()->SetTitleOffset(1.1);
    hTempTPC[i]->GetYaxis()->SetTitle("Counts");
    hTempTPC[i]->GetYaxis()->SetTitleOffset(1.5);
    hTempTPC[i]->SetTitle(Form("%.2f < #it{p}_{T} < %.2f",Binning[i], Binning[i+1]));
    hTempTPC[i]->DrawCopy("E");

    TFitResultPtr FitTriton = hTempTPC[i]->Fit(Triton, "QNL", "", minFit -0.01, maxFit + 0.01);
    Triton->SetLineStyle(2);
    Triton->SetRange(minFit, 4);
    Triton->DrawCopy("same");
    Triton->SetRange(minFit, maxFit);
    Triton->SetLineStyle(1);
    Triton->DrawCopy("same");

    TLegend* LegALICE2 = new TLegend(0.3,0.8,0.9,0.9);
    LegALICE2->SetTextSize(legsize);
    LegALICE2->SetFillColor(0);
    LegALICE2->SetFillStyle(0);
    LegALICE2->SetBorderSize(0);
    LegALICE2->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
    LegALICE2->AddEntry((TObject*)0,"^{3}He, p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    LegALICE2->Draw();

    cFit[i]->cd(2);
    gPad->SetTicks();
    gPad->SetTopMargin(0.07);
    Triton->SetRange(-10, 10);
    hTempTritonSubtracted[i] = (TH1D*) hTempTPC[i]->Clone("hTempTritonSubtracted");

    if (Triton->GetParameter(0) > 0) { // Do not subtract if function is negative
      hTempTritonSubtracted[i]->Add(Triton, -1);
    }
    hTempTritonSubtracted[i]->GetXaxis()->SetRangeUser(-3.5, 3.5);
    hTempTritonSubtracted[i]->DrawCopy("E");
    TFitResultPtr FitHe3 = hTempTritonSubtracted[i]->Fit(Gaus, "QNL", "", -2, 3);
    Gaus->DrawCopy("same");
    fHistoHe3Mean->Fill(Gaus->GetParameter(1));
    fHistoHe3Width->Fill(Gaus->GetParameter(2));
    if (UseMeanAndSigmaFromData) {
      mean = Gaus->GetParameter(1);
      sigma = Gaus->GetParameter(2);
    }
    Double_t nHeliumBin = hTempTritonSubtracted[i]->Integral(hTempTritonSubtracted[i]->FindBin(mean+TPCselectMin*sigma), hTempTritonSubtracted[i]->FindBin(mean+TPCselectMax*sigma));
    nHeliumBin = nHeliumBin*PrimaryFraction[i];
    fHistoHe3PtSpectrum->Fill(Ptbin, nHeliumBin);
    fHistoHe3PtSpectrum->SetBinError(fHistoHe3PtSpectrum->FindBin(Ptbin), TMath::Sqrt( nHeliumBin ) );
    if (SaveToFile) cFit[i]->SaveAs(Form("Helium%.2fGeV.pdf",Ptbin));

    Double_t NallCand = hTempTPC[i]->Integral(hTempTPC[i]->FindBin(mean+TPCselectMin*sigma), hTempTPC[i]->FindBin(mean+TPCselectMax*sigma));
    Double_t NContamination = Triton->Integral(mean+TPCselectMin*sigma,mean+TPCselectMax*sigma);
    printf("Fraction of contamination: %.1f %%\n", NContamination/NallCand *100);

//    // Input for the PID plot in the paper
//    if (Ptbin == 1.75) {
//      tree->Draw(Form("(nSigmaTPC_He3-%f)/%f>>Bin%d(%d,%f,%f)", mean, sigma, i,(int) SelectTPCPlot*10,-SelectTPCPlot/sigma,SelectTPCPlot/sigma),Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q > 0 && TMath::Abs(nSigmaTPC_He3) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY));
//      hTempTPC[i] = (TH1D*) gDirectory->Get(Form("Bin%d",i));
//      hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
//      hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
//      hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
//      hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
//      hTempTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}He}}{#sigma_{d#it{E} / d#it{x}}}");
//      hTempTPC[i]->GetXaxis()->SetTitleOffset(1.1);
//      hTempTPC[i]->GetYaxis()->SetTitle("Counts");
//      hTempTPC[i]->GetYaxis()->SetTitleOffset(1.5);
//      hTempTPC[i]->SetTitle(Form("%.1f GeV/#it{c} < #it{p}_{T} < %.1f GeV/#it{c}",Binning[i], Binning[i+1]));
//      Triton->SetParameter(1,(Triton->GetParameter(1)-mean)/sigma);
//      Triton->SetParameter(2,Triton->GetParameter(2)/sigma);
//      hTempTPC[i]->GetListOfFunctions()->Add(Triton);
//      hTempTPC[i]->SaveAs(Form("Helium%.2fGeV_Sigma.root",Ptbin));
//    }

    // evaluation of anti-3He
    cFitAnti[i] = new TCanvas();
    cFitAnti[i]->SetTitle(Form("Fit %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
    gPad->SetMargin(0.1, 0.05, 0.2, 0.1);
    gPad->SetTicks();

    if (MultiDep){
      tree->Draw(Form("nSigmaTPC_He3>>AntiBin%d(%d,%f,%f)",i,(int) SelectTPCPlot*10,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && TMath::Abs(nSigmaTPC_He3) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && multPercentile_V0A >= %f && multPercentile_V0A < %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY,MultiplicityPercentile[0], MultiplicityPercentile[1]));

    } else {
      tree->Draw(Form("nSigmaTPC_He3>>AntiBin%d(%d,%f,%f)",i,(int) SelectTPCPlot*10,-SelectTPCPlot,SelectTPCPlot),Form("TMath::Abs(eta) <= %f && CalculateRapidity(px,py,pz) >= %f && CalculateRapidity(px,py,pz) < %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && q < 0 && TMath::Abs(nSigmaTPC_He3) < %f &&  %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, SelectTPCPlot, Binning[i], Binning[i+1], DCAcutZ, DCAcutXY));
    }
    hTempTPC[i] = (TH1D*) gDirectory->Get(Form("AntiBin%d",i));
    hTempTPC[i]->GetXaxis()->SetRangeUser(-SelectTPCPlot, SelectTPCPlot);
    hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{3}He}}{#sigma_{d#it{E} / d#it{x}}}");
    hTempTPC[i]->GetXaxis()->SetTitleOffset(1.3);
    hTempTPC[i]->GetYaxis()->SetTitle("Counts");
    hTempTPC[i]->GetYaxis()->SetTitleOffset(1.);
    hTempTPC[i]->SetTitle(Form("%.2f < #it{p}_{T} < %.2f",Binning[i], Binning[i+1]));
    hTempTPC[i]->DrawCopy("E");

    TLegend* LegALICE3 = new TLegend(0.55,0.7,0.95,0.8);
    LegALICE3->SetTextSize(legsize);
    LegALICE3->SetFillColor(0);
    LegALICE3->SetFillStyle(0);
    LegALICE3->SetBorderSize(0);
    LegALICE3->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
    LegALICE3->AddEntry((TObject*)0,"^{3}#bar{He},p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    LegALICE3->Draw();
    TFitResultPtr FitAntiHe3 = hTempTPC[i]->Fit(Gaus, "QNL", "", -2, 3);
    Gaus->DrawCopy("same");
    fHistoHe3Mean->Fill(Gaus->GetParameter(1));
    fHistoHe3Width->Fill(Gaus->GetParameter(2));
    if (UseMeanAndSigmaFromData) {
      mean = Gaus->GetParameter(1);
      sigma = Gaus->GetParameter(2);
    }
    Double_t nAntiHeliumBin = hTempTPC[i]->Integral(hTempTPC[i]->FindBin(mean+TPCselectMin*sigma), hTempTPC[i]->FindBin(mean+TPCselectMax*sigma));
    fHistoAntiHe3PtSpectrum->Fill(Ptbin, nAntiHeliumBin);
    fHistoAntiHe3PtSpectrum->SetBinError(fHistoHe3PtSpectrum->FindBin(Ptbin), TMath::Sqrt( nAntiHeliumBin ) );
    if (SaveToFile) cFitAnti[i]->SaveAs(Form("AntiHelium%.2fGeV.pdf",Ptbin));
    AntiHeCounter += nAntiHeliumBin;
  }

  // normalization (number of events, bin width)
  Double_t RapidityInterval = TMath::Abs(maxRapidity-minRapidity);
  printf("Rapidity interval: %.1f \n", RapidityInterval);
  fHistoHe3PtSpectrum -> Scale(1./(NumberOfEvents*RapidityInterval),"width");
  fHistoAntiHe3PtSpectrum -> Scale(1./(NumberOfEvents*RapidityInterval),"width");

  // correct for eff.
  fHistoHe3PtSpectrum -> Divide(hEff);
  fHistoAntiHe3PtSpectrum -> Divide(hEffAnti);

  if (SubtractHypertritonFeedDown) {
    // subtract hypertriton contamination
    TFile *FeedDownFile = TFile::Open ("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/HypertritonFeedDown.root");
    TH1D* hFeedDownFractionHelium = (TH1D*) FeedDownFile->GetObjectChecked("hFeedDownFractionHelium", "TH1D");
    TH1D* hFeedDownFractionAntiHelium = (TH1D*) FeedDownFile->GetObjectChecked("hFeedDownFractionAntiHelium", "TH1D");

    for (Int_t i=1; i<=hFeedDownFractionHelium->GetNbinsX(); i++) {
      hFeedDownFractionHelium->SetBinContent(i, hFeedDownFractionHelium->GetBinContent(i)*hFeedDownFractionHelium->GetBinWidth(i));
      hFeedDownFractionAntiHelium->SetBinContent(i, hFeedDownFractionAntiHelium->GetBinContent(i)*hFeedDownFractionAntiHelium->GetBinWidth(i));

      hFeedDownFractionHelium->SetBinError(i, hFeedDownFractionHelium->GetBinError(i)*hFeedDownFractionHelium->GetBinWidth(i));
      hFeedDownFractionAntiHelium->SetBinError(i, hFeedDownFractionAntiHelium->GetBinError(i)*hFeedDownFractionAntiHelium->GetBinWidth(i));
    }
    hFeedDownFractionHelium = (TH1D*) hFeedDownFractionHelium->Rebin(nPtBins,"hFeedDownFractionHelium",Binning);
    hFeedDownFractionAntiHelium = (TH1D*) hFeedDownFractionAntiHelium->Rebin(nPtBins,"hFeedDownFractionAntiHelium",Binning);

    hFeedDownFractionHelium->Scale(-1,"width");
    hFeedDownFractionAntiHelium->Scale(-1,"width");

    hFeedDownFractionHelium->Multiply(fHistoHe3PtSpectrum);
    hFeedDownFractionAntiHelium->Multiply(fHistoAntiHe3PtSpectrum);

    TCanvas *cFeedDownSubtraction = new TCanvas("cFeedDownSubtraction","Feed down subtraction",800,600);
    gPad->SetLogy();
    gPad->SetTicks();

    fHistoAntiHe3PtSpectrum->SetMarkerColor(kRed);
    fHistoAntiHe3PtSpectrum->SetLineColor(kRed);
    fHistoHe3PtSpectrum->SetMarkerColor(kBlue);
    fHistoHe3PtSpectrum->SetLineColor(kBlue);

    fHistoAntiHe3PtSpectrum->SetMarkerStyle(20);
    fHistoAntiHe3PtSpectrum->SetMarkerSize(0.5);
    fHistoHe3PtSpectrum->SetMarkerStyle(24);
    fHistoHe3PtSpectrum->SetMarkerSize(0.5);

    fHistoHe3PtSpectrum->DrawCopy("EX0");
    fHistoAntiHe3PtSpectrum->DrawCopy("EX0same");

    hFeedDownFractionHelium->SetMarkerColor(kMagenta+3);
    hFeedDownFractionHelium->SetLineColor(kMagenta+3);
    hFeedDownFractionAntiHelium->SetMarkerColor(kAzure+3);
    hFeedDownFractionAntiHelium->SetLineColor(kAzure+3);
    hFeedDownFractionAntiHelium->SetMarkerStyle(20);
    hFeedDownFractionAntiHelium->SetMarkerSize(0.5);
    hFeedDownFractionHelium->SetMarkerStyle(24);
    hFeedDownFractionHelium->SetMarkerSize(0.5);

    hFeedDownFractionHelium->DrawCopy("EX0same");
    hFeedDownFractionAntiHelium->DrawCopy("EX0same");

    fHistoHe3PtSpectrum->Add(hFeedDownFractionHelium);
    fHistoAntiHe3PtSpectrum->Add(hFeedDownFractionAntiHelium);

    fHistoHe3PtSpectrum->DrawCopy("EX0same");
    fHistoAntiHe3PtSpectrum->DrawCopy("EX0same");
  }


  // Anti-3He / 3He ratio
  fHistoRatioSpectra -> Divide(fHistoAntiHe3PtSpectrum,fHistoHe3PtSpectrum);

  printf("\n %d Anti-^{3}He found \n", AntiHeCounter); // print the total number of anti-3He candidates in the terminal

}
//___________________________________________________________________________________________
void helium3Analysis::WriteOutputFile()  {

  TCanvas *cHeliumQA = new TCanvas("cHeliumQA","Mean and width of the helium distribution");
  fHistoHe3Mean->SetMarkerColor(kRed);
  fHistoHe3Mean->SetLineColor(kRed);
  fHistoHe3Mean->SetMarkerStyle(20);
  fHistoHe3Mean->SetMarkerSize(0.5);

  fHistoHe3Width->SetMarkerColor(kBlue);
  fHistoHe3Width->SetLineColor(kBlue);
  fHistoHe3Width->SetMarkerStyle(20);
  fHistoHe3Width->SetMarkerSize(0.5);

  fHistoHe3Mean->DrawCopy("E");
  fHistoHe3Width->DrawCopy("Esame");

  TCanvas *cSpectra = new TCanvas("cSpectra","Spectra and Ratio",800,600);
  cSpectra->Divide(1,2,0,0); // 0 are margins
  cSpectra->cd(1);
  gStyle->SetPadBottomMargin(0);
  gPad->SetLogy();
  gPad->SetTicks();

  TH1F* hAxis = new TH1F("hAxis","",(PlotMaxPt-PlotMinPt)*10,PlotMinPt,PlotMaxPt);
  hAxis->GetYaxis()->SetRangeUser(8e-9,5e-6);
  if (MultiDep && MultiplicityPercentile[0]== 40) hAxis->GetYaxis()->SetRangeUser(2e-9,5e-6);
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
  leg->AddEntry(fHistoAntiHe3PtSpectrum, "^{3}#bar{He}","l");
  leg->AddEntry(fHistoHe3PtSpectrum, "^{3}He","l");
  if (MultiDep) {
    leg->AddEntry((TObject*)0,Form("Multiplicity: (%.0f - %.0f) %%", MultiplicityPercentile[0], MultiplicityPercentile[1]),"");
  }
  fHistoAntiHe3PtSpectrum->SetMarkerColor(kRed);
  fHistoAntiHe3PtSpectrum->SetLineColor(kRed);
  fHistoHe3PtSpectrum->SetMarkerColor(kBlue);
  fHistoHe3PtSpectrum->SetLineColor(kBlue);

  fHistoAntiHe3PtSpectrum->SetMarkerStyle(20);
  fHistoAntiHe3PtSpectrum->SetMarkerSize(0.5);
  fHistoHe3PtSpectrum->SetMarkerStyle(24);
  fHistoHe3PtSpectrum->SetMarkerSize(0.5);
  fHistoHe3PtSpectrum->GetXaxis()->SetRangeUser(1.5,5.);
  hAxis->DrawCopy();
  fHistoHe3PtSpectrum -> DrawCopy("EX0same");
  fHistoAntiHe3PtSpectrum -> DrawCopy("EX0same");

  //Systematic uncertainties
  // Bin                   1.25  1.75  2.25  2.75  3.25  3.75  4.25  4.75
  Double_t AntiHeSyst[] = {0.25, 0.07, 0.07, 0.06, 0.06, 0.06, 0.06, 0.06};
  Double_t HeSyst[] =     {0.28, 0.13, 0.10, 0.10, 0.05, 0.05, 0.05, 0.05};
  
  //  correct ratio systematic uncertainties (without correlated part)
  Double_t RatioSyst[] = {0.36, 0.14, 0.10, 0.10, 0.03, 0.03, 0.03, 0.03};

  if (MultiDep) {
    if (MultiplicityPercentile[0]==0) {
      HeSyst[5] = (HeSyst[5] + HeSyst[6] + HeSyst[7])/3.;
      AntiHeSyst[5] = (AntiHeSyst[5] + AntiHeSyst[6] + AntiHeSyst[7])/3.;
      
      RatioSyst[5] = (RatioSyst[5] + RatioSyst[6] + RatioSyst[7])/3.;
    } else if (MultiplicityPercentile[0]==10 || MultiplicityPercentile[0]==20){
      HeSyst[4] = (HeSyst[4] + HeSyst[5] + HeSyst[6] + HeSyst[7])/4.;
      AntiHeSyst[4] = (AntiHeSyst[4] + AntiHeSyst[5] + AntiHeSyst[6] + AntiHeSyst[7])/4.;
      
      RatioSyst[4] = (RatioSyst[4] + RatioSyst[5] + RatioSyst[6] + RatioSyst[7])/4.;
    } else if (MultiplicityPercentile[0]==40){
      HeSyst[2] = (HeSyst[2] + HeSyst[3])/2.;
      HeSyst[3] = (HeSyst[4] + HeSyst[5] + HeSyst[6] + HeSyst[7])/4.;
      AntiHeSyst[2] = (AntiHeSyst[2] + AntiHeSyst[3])/2.;
      AntiHeSyst[3] = (AntiHeSyst[4] + AntiHeSyst[5] + AntiHeSyst[6] + AntiHeSyst[7])/4.;

      RatioSyst[2] = (RatioSyst[2] + RatioSyst[3])/2.;
      RatioSyst[3] = (RatioSyst[4] + RatioSyst[5] + RatioSyst[6] + RatioSyst[7])/4.;
    }
  }

  TH1D* fHistoHe3PtSpectrumSyst = (TH1D*) fHistoHe3PtSpectrum->Clone("fHistoHe3PtSpectrumSyst");
  for (int i=1; i <= fHistoHe3PtSpectrumSyst->GetNbinsX(); i++) {
    fHistoHe3PtSpectrumSyst->SetBinError(i, fHistoHe3PtSpectrumSyst->GetBinContent(i)*HeSyst[i-1]);
  }
  TH1D* fHistoAntiHe3PtSpectrumSyst = (TH1D*) fHistoAntiHe3PtSpectrum->Clone("fHistoAntiHe3PtSpectrumSyst");
  for (int i=1; i <= fHistoAntiHe3PtSpectrumSyst->GetNbinsX(); i++) {
    fHistoAntiHe3PtSpectrumSyst->SetBinError(i, fHistoAntiHe3PtSpectrumSyst->GetBinContent(i)*AntiHeSyst[i-1]);
  }

  TH1D* fHistoRatioSpectraSyst = (TH1D*) fHistoAntiHe3PtSpectrumSyst->Clone("fHistoRatioSpectraSyst");
  fHistoRatioSpectraSyst -> Divide(fHistoAntiHe3PtSpectrumSyst,fHistoHe3PtSpectrumSyst);

  // Treat syst uncertainties as 100% correlated above 3 GeV/c
  Double_t RatioError;
  for (int i=1; i <= fHistoRatioSpectraSyst->GetNbinsX(); i++) {
    if (fHistoAntiHe3PtSpectrumSyst->GetBinCenter(i)<3) continue;
    RatioError = RatioSyst[i-1]*fHistoRatioSpectraSyst->GetBinContent(i);
    fHistoRatioSpectraSyst->SetBinError(i, RatioError);
  }

  fHistoHe3PtSpectrumSyst -> SetFillStyle(0);
  fHistoAntiHe3PtSpectrumSyst -> SetFillStyle(0);
  fHistoRatioSpectraSyst -> SetFillStyle(0);

  fHistoHe3PtSpectrum->GetXaxis()->SetRangeUser(1.5,5.);

  fHistoHe3PtSpectrumSyst -> DrawCopy("E2same");
  fHistoAntiHe3PtSpectrumSyst -> DrawCopy("E2same");

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
  gPad->SetTicks();
  gStyle->SetPadTopMargin(0.);
  gStyle->SetPadBottomMargin(0.35);
  fHistoRatioSpectra->SetMarkerColor(kRed);
  fHistoRatioSpectra->SetLineColor(kRed);
  fHistoRatioSpectra->SetMarkerStyle(20);
  fHistoRatioSpectra->SetMarkerSize(0.5);
  fHistoRatioSpectra->GetXaxis()->SetRangeUser(1.5,5.);
  fHistoRatioSpectraSyst->GetXaxis()->SetRangeUser(1.5,5.);

  hAxis->GetYaxis()->SetTitle(RatioLabel);
  hAxis->GetYaxis()->SetNdivisions(505);
  hAxis->GetYaxis()->CenterTitle();
  hAxis->GetYaxis()->SetLabelOffset(0.02);

  if (!MultiDep) {
    hAxis->GetYaxis()->SetRangeUser(0.5,3.5);
  } else{
    if (MultiplicityPercentile[0]==0) {
      hAxis->GetYaxis()->SetRangeUser(0.5,2.2);
    } else if (MultiplicityPercentile[0]==10){
      hAxis->GetYaxis()->SetRangeUser(0.5,4.7);
    } else if (MultiplicityPercentile[0]==20){
      hAxis->GetYaxis()->SetRangeUser(0.5,1.5);
    } else{
      hAxis->GetYaxis()->SetRangeUser(0.3,4.5);
    }
  }

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
  //   TFile *outputFile = new TFile (Form("../Results/Systematic/DCA/Variations/Helium_LHC16q_FAST_woSDD_DCAxy%.2fz%.2f.root",DCAcutXY,DCAcutZ),"RECREATE");

//  TFile *outputFile = new TFile (Form("../Results/Systematic/Tracking/Variations/Helium_LHC16q_FAST_woSDD_ITS%dTPC%d_%.2fy%.2f.root",ITShits, TPCcluster, minRapidity, maxRapidity),"RECREATE");

  //      TFile *outputFile = new TFile (Form("../Results/Systematic/PID/Variations/Helium_LHC16q_FAST_woSDD_TPC%.1f_%.1f.root",TPCselectMin,TPCselectMax),"RECREATE");

  outputFile -> cd();
  fHistoHe3PtSpectrum -> Write();
  fHistoAntiHe3PtSpectrum -> Write();
  fHistoRatioSpectra -> Write();
  //systematic uncertainties
  fHistoHe3PtSpectrumSyst -> Write();
  fHistoAntiHe3PtSpectrumSyst -> Write();
  fHistoRatioSpectraSyst -> Write();
  leg->Write();
  outputFile ->Close();
}
//___________________________________________________________________________________________
Double_t CalculateRapidity(Double_t px, Double_t py, Double_t pz){
  // calculation of the rapidity for 3He
  Double_t p = TMath::Sqrt(px*px + py*py+ pz*pz);
  Double_t energyHe3 = TMath::Sqrt(p*2*p*2 + AliPID::ParticleMass(AliPID::kHe3)*AliPID::ParticleMass(AliPID::kHe3));
  Double_t rapHe3 = 0.5*TMath::Log((energyHe3 + pz*2)/(energyHe3 - pz*2));
  // Correction with system rapidity
  rapHe3 = rapHe3 - 0.465;

  return rapHe3;
}

//___________________________________________________________________________________________
Double_t CorrectPt(Double_t px, Double_t py){

  //      TF1* PtCorr = new TF1("PtCorr","[0]+[1]*exp(x^3*[2])",0,10);
  //   // new extracted from (pt_true - pT_rec) CENT fit to profile in between 1.3 and 9.5 GeV/c
  //            PtCorr->SetParameter(0, 6.55239e-04);   //   +/-   4.07690e-05
  //            PtCorr->SetParameter(1, -3.69300e-01);      //   +/-   9.36303e-04
  //            PtCorr->SetParameter(2, -2.40396e-01);         //  +/-   5.34414e-04
  //
  //      // new extracted from (pt_true - pT_rec) FAST fit to profile in between 1.3 and 9.5 GeV/c
  //         PtCorr->SetParameter(0, 7.05001e-04);   //   +/-   4.02415e-05
  //         PtCorr->SetParameter(1, -3.59028e-01);      //   +/-   9.35793e-04
  //         PtCorr->SetParameter(2, -2.36844e-01);         //  +/-   5.37830e-04
  //   =>
  //      PtCorr->SetParameter(0, 0.00068);   //   +/-   0.00004
  //      PtCorr->SetParameter(1, -0.3641);         //  +/-   0.0009
  //      PtCorr->SetParameter(2, -0.2386);      //   +/-   0.0005

  Double_t pT = TMath::Sqrt(px*px + py*py);
  Double_t correctedpT = 2*pT + 0.00068 - 0.3641 * TMath::Exp(-0.2386* TMath::Power(2*pT,3.0));

  return correctedpT;
}
