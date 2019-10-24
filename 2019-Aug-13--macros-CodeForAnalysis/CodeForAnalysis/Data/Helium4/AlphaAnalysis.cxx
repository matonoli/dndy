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

// needed for limit calculation
#include "TFeldmanCousins2.h" // modified root class to allow to take systematic uncertainties into account
#include "TFeldmanCousins2.cxx" // modified root class to allow to take systematic uncertainties into account
#include "/Users/tidus/alice/ali-master/AliPhysics/PWG/Tools/AliPWGFunc.h"

#include "AlphaAnalysis.h"

Double_t legsize = 0.05;
Double_t labelsize = 0.07;

Double_t DCAcutXY;
Double_t DCAcutZ;
Double_t minRapidity;
Double_t maxRapidity;
Int_t ITShits;
Int_t TPCcluster;

// Binning
const Int_t nPtBins=1;
Double_t Binning[]={2., 10.};

ClassImp(AlphaAnalysis);

//___________________________________________________________________________________________
AlphaAnalysis::AlphaAnalysis():
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
AlphaAnalysis::~AlphaAnalysis()
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
void AlphaAnalysis::GetInputTree()  {
  
  TFile *inputFile = TFile::Open (Form("%s",fInputFileName));
  TList *inputList = (TList*) inputFile->Get(fInputListName);
  
  tree = (TTree*) inputFile->GetObjectChecked("Trees/reducedTree_Helium","TTree"); // tree is stored in a subfolder called Trees
  
  // Get Number of events
  TH1F* EventSelection = (TH1F*) inputList->FindObject("histoEventSelection");
  NumberOfEvents = EventSelection->GetBinContent(3);
  printf("Number of events: %d \n", NumberOfEvents);
}
//___________________________________________________________________________________________
void AlphaAnalysis::SetTreeBranches()  {
  
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
  //   tree -> SetBranchAddress("nSigmaITS_Trit",&nSigmaITS_Trit);
  //   tree -> SetBranchAddress("nSigmaTPC_Trit",&nSigmaTPC_Trit);
  //   tree -> SetBranchAddress("nSigmaTOF_Trit",&nSigmaTOF_Trit);
  //   tree -> SetBranchAddress("nSigmaTRD_Trit",&nSigmaTRD_Trit);
  //   tree -> SetBranchAddress("nSigmaHMPID_Trit",&nSigmaHMPID_He4);
  //      tree -> SetBranchAddress("nSigmaITS_Deut",&nSigmaITS_Deut);
  //      tree -> SetBranchAddress("nSigmaTPC_Deut",&nSigmaTPC_Deut);
  //      tree -> SetBranchAddress("nSigmaTOF_Deut",&nSigmaTOF_Deut);
  //      tree -> SetBranchAddress("nSigmaTRD_Deut",&nSigmaTRD_Deut);
  //      tree -> SetBranchAddress("nSigmaHMPID_Deut",&nSigmaHMPID_Deut);
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
void AlphaAnalysis::CreateHistograms()  {
  
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
  fHistoRatioSpectra->GetYaxis()->SetTitle("#bar{^{4}He} / ^{4}He");
  fHistoRatioSpectra->GetYaxis()->SetTitleSize(labelsize+0.01);
  fHistoRatioSpectra->GetYaxis()->SetTitleOffset(0.8);
  fHistoRatioSpectra->GetYaxis()->SetLabelSize(labelsize);
  
  fHistoPtSpectrum->Sumw2();
  fHistoAntiPtSpectrum->Sumw2();
  fHistoRatioSpectra->Sumw2();
  
  fHistoRapidity = new TH1D ("fHistoRapidity","Rapidity",30,-1.5,1.5);
  fHistoMean = new TH1D ("fHistoMean","Mean of the Alpha distribution",24,-3.,3.);
  fHistoWidth = new TH1D ("fHistoWidth","Width of the Alpha distribution",12,0.,3.);
  
}
//___________________________________________________________________________________________
void AlphaAnalysis::ExecuteAnalysis()  {
  
  Double_t SystematicUncertainy = 0.2; // Taken from 4He analysis in Pb--Pb ()
  
  // defining size of pad ---------------------------
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  
  SaveToFile = false; // Save efficiency and pT slice fits to pdf
  Bool_t DrawRapidity = false;
  
  // Track selection
  Double_t maxEta = 0.9; // Pseudorapidity range
  minRapidity = -1.;
  maxRapidity = 0.;
  ITShits = 2;
  TPCcluster = 70;
  DCAcutXY = 2.4;
  DCAcutZ = 3.2;
  
  Bool_t DoTPConly = false; // not reasonable (efficiency too steep below 2.5 GeV/c and above about 3 GeV/c 3He starts to merge in)
  Double_t TOFselect = 5;
  Double_t TPCselect = 3;
  
  /// Blast wave parameterisation for 4He
  AliPWGFunc* Functions = new AliPWGFunc();
  Functions->SetVarType(AliPWGFunc::kdNdpt);
  
  //Use mT-exponential instead of BW
  TF1* pTshape = Functions->GetMTExp(3.73, 0.597, 1, "#it{m}_{T} exponential");
  pTshape->FixParameter(1, 0.597); // restrict beta to positive values and below 1
  pTshape->FixParameter(2, 1.); // fix temperatur
  pTshape->SetLineStyle(2);
  pTshape->SetLineWidth(2);
  pTshape->SetLineColor(kBlue);
  
  //  TF1* pTshape = Functions->GetBGBW(3.73, 0.75, 0.161, 2.95, 100., "BlastWave");
  //  pTshape->FixParameter(1, 0.75); // restrict beta to positive values and below 1
  //  pTshape->FixParameter(2, 0.161); // fix temperatur
  //  pTshape->FixParameter(3, 2.95); // fix n
  
  /// Load efficiency file and get anti-4He efficiency
  TFile *EffFile = TFile::Open ("/Users/tidus/Documents/Studium/PhD/Nuclei/Results/MC_FAST/Efficiency_Alpha.root");
  TH1D* hEffAnti4He = (TH1D*) EffFile->GetObjectChecked("hEffAnti4He", "TH1D");
  hEffAnti4He->SetTitle("");
  hEffAnti4He->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hEffAnti4He->GetYaxis()->SetTitle("Acceptance #times efficiency");
  hEffAnti4He->GetYaxis()->SetLabelSize(legsize);
  hEffAnti4He->GetXaxis()->SetLabelSize(legsize);
  hEffAnti4He->GetYaxis()->SetTitleSize(legsize+0.01);
  hEffAnti4He->GetXaxis()->SetTitleSize(legsize+0.01);
  hEffAnti4He->GetYaxis()->SetTitleOffset(1);
  hEffAnti4He->GetYaxis()->SetNdivisions(505);
  hEffAnti4He->GetXaxis()->SetNdivisions(505);
  
  // fit efficiency
  TF1* EffParam = new TF1("EffParam","[0]+[1]*exp([2]*x)", 0,10);
  const Double_t StartParam[3] = {0.6, -14, -2};
  EffParam->SetParameters(StartParam);
  TCanvas* cEff = new TCanvas("cEff","Efficiencies", 800, 800);
  cEff->SetTicks();
  TLegend* LegALICE = new TLegend(0.30,0.6,0.75,0.85);
  LegALICE->SetTextSize(legsize);
  LegALICE->SetFillColor(0);
  LegALICE->SetFillStyle(0);
  LegALICE->SetBorderSize(0);
  LegALICE->AddEntry((TObject*)0,"ALICE","");
  LegALICE->AddEntry((TObject*)0,"p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
  LegALICE->AddEntry((TObject*)0,"#minus1 #scale[1.2]{#leq} #it{y}_{cms} < 0","");
  LegALICE->AddEntry((TObject*)0,"^{4}#bar{He}","");
  
  hEffAnti4He->Fit(EffParam,"","",2,10);
  LegALICE->Draw();
  cEff->SaveAs("4He_Efficiency.eps");
  
  TF1* weightedEff = new TF1("weightedEff", WeightedEfficiency,0,10,8);
  weightedEff->FixParameter(0, pTshape->GetParameter(0));
  weightedEff->FixParameter(1, pTshape->GetParameter(1));
  weightedEff->FixParameter(2, pTshape->GetParameter(2));
  //  weightedEff->FixParameter(3, pTshape->GetParameter(3)); // needed for BW
  //  weightedEff->SetParameter(4, pTshape->GetParameter(4));
  
  weightedEff->FixParameter(5, EffParam->GetParameter(0)); // fix eff. param.
  weightedEff->FixParameter(6, EffParam->GetParameter(1));
  weightedEff->FixParameter(7, EffParam->GetParameter(2));
  Double_t Efficiency = weightedEff->Integral(Binning[0], Binning[1]) / pTshape->Integral(Binning[0], Binning[1]);
  printf("Average efficiency for 4He in %f GeV/c < #it{p}_{T} < %f GeV/c: %f \n", Binning[0], Binning[1], Efficiency);
  // Extrapolate to full pT
  Double_t IntegralFraction = pTshape->Integral(Binning[0], Binning[1])/pTshape->Integral(0,10); // Fraction of the total yield inside the used pT range (Not yet implemented)
  printf("Fraction of yield in %f GeV/c < #it{p}_{T} < %f GeV/c: %f \n \n", Binning[0], Binning[1], IntegralFraction);
  
  // Evaluate the rapidity distributrion
  if (DrawRapidity) {
    TCanvas *cRapidity = new TCanvas("cRapidity","Rapidity");
    TLegend* LegALICERap = new TLegend(0.55,0.8,0.95,0.925);
    LegALICERap->SetTextSize(legsize);
    LegALICERap->SetFillColor(0);
    LegALICERap->SetFillStyle(0);
    LegALICERap->SetBorderSize(0);
    LegALICERap->AddEntry((TObject*)0,"ALICE #it{Work in progress}","");
    LegALICERap->AddEntry((TObject*)0,"p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","");
    
    tree->Draw("CalculateRapidity(y)>>fHistoRapidity(30,-1.5,1.5)",Form("TMath::Abs(eta) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaz) <= %f && TMath::Abs(dcaxy) <= %f && mass > %f && mass < %f", maxEta, ITShits, TPCcluster, DCAcutZ, DCAcutXY, LowerMass, UpperMass));
    
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
  
  Double_t Ptbin = 0;
  Int_t AntiAlphaCounter = 0;
  Double_t BackgroundEstimate = 0;
  Double_t BackgroundEstimate2 = 0;
  
  TF1* Background = new TF1("Background", "gaus", -5.,5.);
  TF1* Background2 = new TF1("Background", "expo", -5.,5.);
  
  TCanvas* cTOF[nPtBins];
  TCanvas* cTPC[nPtBins];
  TCanvas* cTOFTPC[nPtBins];
  
  TH1D *hTempTOF[nPtBins];
  TH1D *hTempTPC[nPtBins];
  TH2D* hTempTOFTPC[nPtBins];
  
  for (Int_t i =0; i < nPtBins; i++) {
    Ptbin = (Binning[i]+Binning[i+1])/2;
    
    //    // TOF distribution
    //    cTOF[i] = new TCanvas();
    //    cTOF[i]->SetTitle(Form("TOF %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
    //    cTOF[i]->Divide(2,1);
    //    cTOF[i]->cd(1);
    //    gPad->SetTicks();
    //    tree->Draw(Form("nSigmaTOF_He4>>TOFtemp_%d(100,-5.,5.)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q < 0 && abs(nSigmaTPC_He4) < 5", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1]));
    //    hTempTOF[i] = (TH1D*) gDirectory->Get(Form("TOFtemp_%d",i));
    //    hTempTOF[i]->SetTitle("Anti-alpha");
    //    hTempTOF[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOF[i]->GetXaxis()->SetLabelSize(legsize);
    //    hTempTOF[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOF[i]->GetYaxis()->SetLabelSize(legsize);
    //    hTempTOF[i]->GetXaxis()->SetTitle("#it{t}_{TOF}-#it{t}_{^{4}He} (#sigma_{TOF})");
    //    hTempTOF[i]->GetYaxis()->SetTitle("Counts");
    //    hTempTOF[i]->GetXaxis()->SetTitleOffset(1.1);
    //    hTempTOF[i]->GetYaxis()->SetTitleOffset(1.5);
    //    hTempTOF[i]->DrawCopy();
    //
    //    cTOF[i]->cd(2);
    //    gPad->SetTicks();
    //    tree->Draw(Form("nSigmaTOF_He4>>TOFtemp_%d(100,-5.,5..)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q > 0 && abs(nSigmaTPC_He4) < 5", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1]));
    //    hTempTOF[i] = (TH1D*) gDirectory->Get(Form("TOFtemp_%d",i));
    //    hTempTOF[i]->SetTitle("Alpha");
    //    hTempTOF[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOF[i]->GetXaxis()->SetLabelSize(legsize);
    //    hTempTOF[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOF[i]->GetYaxis()->SetLabelSize(legsize);
    //    hTempTOF[i]->GetXaxis()->SetTitle("#it{t}_{TOF}-#it{t}_{^{4}He} (#sigma_{TOF})");
    //    hTempTOF[i]->GetYaxis()->SetTitle("Counts");
    //    hTempTOF[i]->GetXaxis()->SetTitleOffset(1.1);
    //    hTempTOF[i]->GetYaxis()->SetTitleOffset(1.5);
    //    hTempTOF[i]->DrawCopy();
    
    // TPC distribution
    cTPC[i] = new TCanvas("cTPC",Form("TPC %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]),800,800);
    //    cTPC[i]->Divide(2,1);
    //    cTPC[i]->cd(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetTicks();
    
    if (DoTPConly) {
      tree->Draw(Form("nSigmaTPC_He4>>TPCtemp_%d(100,-5.,5.)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q < 0", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1]));
    }else {
      tree->Draw(Form("nSigmaTPC_He4>>TPCtemp_%d(100,-5.,5.)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q < 0 && abs(nSigmaTOF_He4)<5", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1], TOFselect));
    }
    hTempTPC[i] = (TH1D*) gDirectory->Get(Form("TPCtemp_%d",i));
    hTempTPC[i]->SetTitle("");
    hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
    hTempTPC[i]->GetYaxis()->SetNdivisions(505);
    hTempTPC[i]->GetXaxis()->SetNdivisions(505);
    hTempTPC[i]->GetXaxis()->SetTitle("#it{n}_{#sigma}^{TPC}(^{4}He)"); //"#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{4}He}}{#sigma_{d#it{E} / d#it{x}}}"
    hTempTPC[i]->GetYaxis()->SetTitle("Counts");
    hTempTPC[i]->GetYaxis()->SetTitleOffset(1.2);
    hTempTPC[i]->GetXaxis()->SetTitleOffset(1.1);
    hTempTPC[i]->SetMarkerStyle(20);
    hTempTPC[i]->SetMarkerSize(1.5);
    hTempTPC[i]->SetMarkerColor(1);
    hTempTPC[i]->SetLineColor(1);
    hTempTPC[i]->DrawCopy("EX0");
    LegALICE->Draw();
    
    hTempTPC[i]->Fit(Background,"0","",-5,-3);
    Background->SetLineStyle(1);
    Background->SetLineColor(kRed);
    Background->DrawCopy("same");
    BackgroundEstimate = Background->Integral(-3., 3.);
    hTempTPC[i]->Fit(Background2,"N","",-5,-3);
    Background2->SetLineStyle(2);
    Background2->SetLineColor(kBlue);
    //    Background2->DrawCopy("same");
    BackgroundEstimate2 = Background2->Integral(-3., 3.);
    
    TLegend* LegTPC = new TLegend(0.4,0.4,0.7,0.55);
    LegTPC->SetTextSize(legsize);
    LegTPC->SetFillColor(0);
    LegTPC->SetFillStyle(0);
    LegTPC->SetBorderSize(0);
    LegTPC->AddEntry(hTempTPC[i],"Data","p");
    LegTPC->AddEntry(Background,"Background","l");
    LegTPC->Draw();
    
    AntiAlphaCounter += hTempTPC[i]->Integral(hTempTPC[i]->FindBin(-TPCselect), hTempTPC[i]->FindBin(TPCselect));
    cTPC[i]->SaveAs("4He_TPC_selection_wTOF.eps");
    
    //    cTPC[i]->cd(2);
    //    gPad->SetTicks();
    //    if (DoTPConly) {
    //      tree->Draw(Form("nSigmaTPC_He4>>TPCtemp_%d(100,-5.,5.)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q > 0", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1]));
    //    } else {
    //      tree->Draw(Form("nSigmaTPC_He4>>TPCtemp_%d(100,-5.,5.)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q > 0 && abs(nSigmaTOF_He4)<5 ", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1], TOFselect));
    //    }
    //    hTempTPC[i] = (TH1D*) gDirectory->Get(Form("TPCtemp_%d",i));
    //    hTempTPC[i]->SetTitle("Alpha");
    //    hTempTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    //    hTempTPC[i]->GetXaxis()->SetLabelSize(legsize);
    //    hTempTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    //    hTempTPC[i]->GetYaxis()->SetLabelSize(legsize);
    //    hTempTPC[i]->GetYaxis()->SetNdivisions(505);
    //    hTempTPC[i]->GetXaxis()->SetNdivisions(505);
    //    hTempTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{4}He}}{#sigma_{d#it{E} / d#it{x}}}");
    //    hTempTPC[i]->GetYaxis()->SetTitle("Counts");
    //    hTempTPC[i]->GetYaxis()->SetTitleOffset(1.2);
    //    hTempTPC[i]->DrawCopy();
    //
    //    hTempTPC[i]->Fit(Background,"QN","",-5,-3);
    //    Background->SetLineStyle(2);
    //    Background->SetLineColor(kRed);
    //    Background->DrawCopy("same");
    //    hTempTPC[i]->Fit(Background2,"QN","",-5,-3);
    //    Background2->SetLineStyle(2);
    //    Background2->SetLineColor(kBlue);
    //    Background2->DrawCopy("same");
    
    
    //    // TOF vs TPC
    //    cTOFTPC[i] = new TCanvas();
    //    cTOFTPC[i]->SetTitle(Form("TOF vs TPC %f < #it{p}_{T} < %f",Binning[i], Binning[i+1]));
    //    cTOFTPC[i]->Divide(2,1);
    //    cTOFTPC[i]->cd(1);
    //    gPad->SetTicks();
    //    tree->Draw(Form("nSigmaTOF_He4:nSigmaTPC_He4>>TOFTPCtemp_%d(100,-5.,5., 100, -5., 5.)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q < 0", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1]));
    //    hTempTOFTPC[i] = (TH2D*) gDirectory->Get(Form("TOFTPCtemp_%d",i));
    //    hTempTOFTPC[i]->SetTitle("Anti-alpha");
    //    hTempTOFTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOFTPC[i]->GetXaxis()->SetLabelSize(legsize);
    //    hTempTOFTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOFTPC[i]->GetYaxis()->SetLabelSize(legsize);
    //    hTempTOFTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{4}He}}{#sigma_{d#it{E} / d#it{x}}}");
    //    hTempTOFTPC[i]->GetYaxis()->SetTitle("#it{t}_{TOF}-#it{t}_{^{4}He} (#sigma_{TOF})");
    //    hTempTPC[i]->GetZaxis()->SetTitle("Counts");
    //    hTempTOFTPC[i]->GetXaxis()->SetTitleOffset(1.1);
    //    hTempTOFTPC[i]->GetYaxis()->SetTitleOffset(1.5);
    //    hTempTOFTPC[i]->DrawCopy("COLZ");
    //
    //    cTOFTPC[i]->cd(2);
    //    gPad->SetTicks();
    //    tree->Draw(Form("nSigmaTOF_He4:nSigmaTPC_He4>>TOFTPCtemp_%d(100,-5.,5., 100, -5., 5.)",i), Form("TMath::Abs(eta) <= %f && CalculateRapidity(y) >= %f && CalculateRapidity(y) <= %f && nITS_Clusters >= %d && nTPC_Clusters >= %d && TMath::Abs(dcaxy)<= %f && TMath::Abs(dcaz)<= %f && %f <= CorrectPt(px, py) && CorrectPt(px, py) < %f && q > 0", maxEta, minRapidity, maxRapidity, ITShits, TPCcluster, DCAcutXY, DCAcutZ, Binning[i], Binning[i+1]));
    //    hTempTOFTPC[i] = (TH2D*) gDirectory->Get(Form("TOFTPCtemp_%d",i));
    //    hTempTOFTPC[i]->SetTitle("Alpha");
    //    hTempTOFTPC[i]->GetXaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOFTPC[i]->GetXaxis()->SetLabelSize(legsize);
    //    hTempTOFTPC[i]->GetYaxis()->SetTitleSize(legsize+0.01);
    //    hTempTOFTPC[i]->GetYaxis()->SetLabelSize(legsize);
    //    hTempTOFTPC[i]->GetXaxis()->SetTitle("#frac{d#it{E} / d#it{x} - #LTd#it{E} / d#it{x}#GT_{^{4}He}}{#sigma_{d#it{E} / d#it{x}}}");
    //    hTempTOFTPC[i]->GetYaxis()->SetTitle("#it{t}_{TOF}-#it{t}_{^{4}He} (#sigma_{TOF})");
    //    hTempTPC[i]->GetZaxis()->SetTitle("Counts");
    //    hTempTOFTPC[i]->GetXaxis()->SetTitleOffset(1.1);
    //    hTempTOFTPC[i]->GetYaxis()->SetTitleOffset(1.5);
    //    hTempTOFTPC[i]->DrawCopy("COLZ");
  }
  
  //    normalization (number of events, bin width)
  Double_t RapidityInterval = TMath::Abs(maxRapidity-minRapidity);
  printf("Rapidity interval: %.1f \n", RapidityInterval);
  printf("\n %d #bar{^4He} found, expected background: %.2e (gaus), %.2e (expo) \n", AntiAlphaCounter, BackgroundEstimate, BackgroundEstimate2); // print the total number of anti- candidates in the terminal
  
  TFeldmanCousins2* Limit = new TFeldmanCousins2();
  Double_t UpperLimit4He = Limit->CalculateUpperLimit(AntiAlphaCounter, BackgroundEstimate,SystematicUncertainy);
  
  // Normalise by number of events, rapidity
  UpperLimit4He = UpperLimit4He/(NumberOfEvents*RapidityInterval*Efficiency);
  // Extrapolate to full pT
  UpperLimit4He = UpperLimit4He/IntegralFraction;
  printf("Upper limit on the integrated helium-4 yield (bkg=gaus): %.2e \n", UpperLimit4He);
  
  Double_t UpperLimit4HeExp = Limit->CalculateUpperLimit(AntiAlphaCounter, BackgroundEstimate2,SystematicUncertainy);
  // Normalise by number of events, rapidity
  UpperLimit4HeExp = UpperLimit4HeExp/(NumberOfEvents*RapidityInterval*Efficiency);
  // Extrapolate to full pT
  UpperLimit4HeExp = UpperLimit4HeExp/IntegralFraction;
  printf("Upper limit on the integrated helium-4 yield (bkg=exp): %.2e \n", UpperLimit4HeExp);
}
//___________________________________________________________________________________________
void AlphaAnalysis::WriteOutputFile()  {
}
//___________________________________________________________________________________________
Double_t CalculateRapidity(Double_t y){
  return  y - 0.465;
}

//___________________________________________________________________________________________
Double_t CorrectPt(Double_t px, Double_t py){
  Double_t pT = TMath::Sqrt(px*px + py*py);
  Double_t correctedpT = 2*pT- 0.00295074 + 0.332899*TMath::Exp(-0.815359*2*pT);
  
  return TMath::Sqrt(px*px + py*py);
}
//_____________________________________________________________________________________
Double_t WeightedEfficiency(Double_t *x, Double_t *par){
  
  AliPWGFunc* Functions = new AliPWGFunc();
  Functions->SetVarType(AliPWGFunc::kdNdpt);
  
  TF1* pTshape = Functions->GetMTExp(3.73, 0.597, 1, "#it{m}_{T} exponential");
  
  //  TF1* pTshape = Functions->GetBGBW(3.73, 0.7, 0.16, 2, 100, "BlastWave");
  pTshape->SetParameter(0, par[0]);
  pTshape->SetParameter(1, par[1]);
  pTshape->SetParameter(2, par[2]);
  //  pTshape->SetParameter(3, par[3]);
  //  pTshape->SetParameter(4, par[4]);
  
  TF1* EffParam = new TF1("EffParam","[0]+[1]*exp([2]*x)", 0,10);
  EffParam->SetParameter(0, par[5]);
  EffParam->SetParameter(1, par[6]);
  EffParam->SetParameter(2, par[7]);
  
  return pTshape(x[0])*EffParam(x[0]); // exponential
}
