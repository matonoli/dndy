#include <stdio.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFitResultPtr.h>
#include "/Users/tidus/alice/ali-master/AliPhysics/PWG/Tools/AliPWGFunc.h"
#include "YieldMean.C"

using namespace std;

Double_t minpT = 0.;
Double_t maxpT = 10.;

Double_t minPtMeasured = 1.5;
Double_t maxPtMeasured = 5.;

Double_t GetMaximum(Double_t input[], Int_t ArraySize){
  Double_t max = input[0];
  for (Int_t i=1; i < ArraySize; i++) {
    if(input[i] > max) max = input[i];
  }
  return max;
}

Double_t GetMinimum(Double_t input[], Int_t ArraySize){
  Double_t min = input[0];
  for (Int_t i=1; i < ArraySize; i++) {
    if(input[i] < min) min = input[i];
  }
  return min;
}

void MakeFits_YieldMean(){

  TString file1="Helium_LHC16q_FAST_woSDD_0_10.root",
  file2="Helium_LHC16q_FAST_woSDD_Mult_Int.root",
  file3="Helium_LHC16q_FAST_woSDD_10_20.root",
  file4="Helium_LHC16q_FAST_woSDD_20_40.root",
  file5="Helium_LHC16q_FAST_woSDD_40_100.root";

  TString Objectname = "fHistoHe3PtSpectrum"; // 3He
  //  TString Objectname = "fHistoAntiHe3PtSpectrum"; // Anti-3He (Tsallis does not converge for 0-10%)

  TString Objectname2 = "fHistoAntiHe3PtSpectrum"; // needed if average spectrum used

  Bool_t average = true;
  Bool_t CorrelatedUncertainties = true; // between particle and anti-particle: up to 3 GeV primary fraction dominate -> uncorrelated, above tracking dominate -> correlated

  Bool_t RemoveBinByBinCorrelatedSyst = false; // for mean pT calculation

  TString LegendTitle = "#bf{ALICE Work in progress}"; // ALICE Work in progress / ALICE Preliminary
  TString AnalysisName = "(^{3}He + ^{3}#bar{He})/2, p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV";
  //  TString AnalysisName = "^{3}He, p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV";
  //  TString AnalysisName = "^{3}#bar{He}, p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV";

  TString RatioLabel = "(^{3}He + ^{3}#bar{He}) / (p + #bar{p})";
  //  TString RatioLabel = "2^{3}He / (p + #bar{p})";
  //  TString RatioLabel = "2^{3}#bar{He} / (p + #bar{p})";

  TString MultiplicityLabel = "#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{#lbar#it{#eta}_{lab}#lbar< 0.5}";

  Int_t nFiles = 5;

  //style area
  Double_t legsize = 0.04;
  Double_t textsize = 0.05;
  gSystem->Load("mystyle_C.so");
  mystyle();

  TString pTLabel= "#it{p}_{T} (GeV/#it{c})";
  TString yLable = "1/#it{N}_{evt} d^{2}#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-1}"; // for spectra

  // Fit functions
  // start parameters
  Double_t mass = 2.81; // Helium-3 mass
  Double_t beta = 0.8;
  Double_t T = 0.55;
  Double_t TTsallis = 0.156;
  Double_t n = 2.5;
  Double_t norm = 5e-05;
  Double_t normTsallis = 2.e-6;
  // Set function variable
  AliPWGFunc* Functions = new AliPWGFunc();
  Functions->SetVarType(AliPWGFunc::kdNdpt);
  TF1* fFitFunctions[7];
  //pT exponential
  fFitFunctions[0] = Functions->GetPTExp(T, norm, "#it{p}_{T} exponential");
  fFitFunctions[0]->SetLineStyle(2);
  fFitFunctions[0]->SetLineWidth(2);
  fFitFunctions[0]->SetLineColor(kCyan-2);
  //mT exponential
  fFitFunctions[1] = Functions->GetMTExp(mass, T, norm, "#it{m}_{T} exponential");
  fFitFunctions[1]->SetLineStyle(2);
  fFitFunctions[1]->SetLineWidth(2);
  fFitFunctions[1]->SetLineColor(kBlue);
  // Boltzmann
  fFitFunctions[2] = Functions->GetBoltzmann(mass, T, norm, "Boltzmann");
  fFitFunctions[2]->SetLineStyle(2);
  fFitFunctions[2]->SetLineWidth(2);
  fFitFunctions[2]->SetLineColor(kBlack);
  // Bose-Einstein
  fFitFunctions[3] = Functions->GetBoseEinstein(mass, T, norm, "Bose-Einstein");
  fFitFunctions[3]->SetLineStyle(2);
  fFitFunctions[3]->SetLineWidth(2);
  fFitFunctions[3]->SetLineColor(kOrange-3);
  // Fermi-Dirac
  fFitFunctions[4] = Functions->GetFermiDirac(mass, T, norm, "Fermi-Dirac");
  fFitFunctions[4]->SetLineStyle(2);
  fFitFunctions[4]->SetLineWidth(2);
  fFitFunctions[4]->SetLineColor(kGreen-3);

  // blast wave (fit only converging if n and T are fixed)
  fFitFunctions[5] = Functions->GetBGBW(mass, beta, TTsallis, n, normTsallis, "Boltzmann-Gibbs blast wave");
  fFitFunctions[5]->SetLineStyle(2);
  fFitFunctions[5]->SetLineWidth(2);
  fFitFunctions[5]->SetParLimits(1,0,1.); // restrict beta to positive values and below 1
  fFitFunctions[5]->SetParLimits(2, 0.1, 0.18); // set limit for T
  fFitFunctions[5]->FixParameter(2, 0.161); // fix temperatur
  fFitFunctions[5]->SetParLimits(3, 1., 3.5); // set limit for n
  fFitFunctions[5]->SetParameter(3, 3.); //  n
  fFitFunctions[5]->SetLineColor(kCyan+1);

  // Lévy-Tsallis
  fFitFunctions[6] = Functions->GetLevi(mass, TTsallis, n, normTsallis, "Lévy-Tsallis");
  fFitFunctions[6]->SetLineStyle(2);
  fFitFunctions[6]->SetLineWidth(2);
  fFitFunctions[6]->SetLineColor(kMagenta-3);
  fFitFunctions[6]->SetParLimits(1, 0, 20.); // set limit for n
  fFitFunctions[6]->SetParameter(1, 10);
  fFitFunctions[6]->SetParLimits(2, 0.1, 0.5); // set limit for T
  fFitFunctions[6]->SetParameter(2, 0.47); //
  fFitFunctions[6]->FixParameter(3,mass); // fix mass to He mass

  const Int_t NumberOfFitFunctions = 5; // 7 with Tsallis; 6 is without Tsallis; 5 without Tsallis and BW

  Double_t MeanMultiplicity[4] = {40.6, 30.5, 23.2, 10.1}; // taken from old analysis to be updated
  Double_t MultiplicityUncertainty[4] = {0.9, 0.7, 0.5 , 0.2}; // taken from old analysis to be updated

  // Proton integrated yield				0-5	5-10	10-20	20-40	40-60	60-80	80-100
  Double_t ProtonIntegratedYield[7] = {2.280446e+00, 1.853576e+00, 1.577375e+00, 1.221875e+00, 8.621290e-01, 5.341445e-01, 2.307767e-01};
  Double_t ProtonIntegratedYieldStatisticalUncertainties[7] = {4.394273e-03, 3.901768e-03, 2.739986e-03, 1.735134e-03, 1.381181e-03, 1.034190e-03, 6.969543e-04};
  Double_t ProtonIntegratedYieldSystematicUncertainties[7] = {1.585968e-01, 1.284426e-01, 1.084793e-01, 8.190972e-02, 5.608663e-02, 3.408551e-02, 1.450569e-02};
  Double_t ProtonIntegratedYieldUncorrelatedUncertainties[7] = {7.001060e-02, 5.548104e-02, 4.467687e-02, 3.347605e-02, 2.265973e-02, 1.377331e-02, 6.911692e-03};

  Double_t ProtonBuffer=0;
  Double_t ProtonBufferStat=0;
  Double_t ProtonBufferSyst=0;
  Double_t ProtonBufferWeight=0;
  for (Int_t j=0; j <7; j++) {
    // Add uncorrelated uncertainty to the systematic uncertainties
    //      ProtonIntegratedYieldSystematicUncertainties[j] = TMath::Sqrt(ProtonIntegratedYieldSystematicUncertainties[j]*ProtonIntegratedYieldSystematicUncertainties[j] +  ProtonIntegratedYieldUncorrelatedUncertainties[j]*ProtonIntegratedYieldUncorrelatedUncertainties[j]);
    if (j==2) {
      ProtonBufferWeight=2;
    } else if (j>2){
      ProtonBufferWeight=4;
    } else{
      ProtonBufferWeight=1;
    }
    ProtonBuffer += ProtonIntegratedYield[j]*ProtonBufferWeight;
    ProtonBufferSyst += ProtonIntegratedYieldSystematicUncertainties[j]*ProtonBufferWeight;
    ProtonBufferStat+= ProtonIntegratedYieldStatisticalUncertainties[j]*ProtonIntegratedYieldStatisticalUncertainties[j]*ProtonBufferWeight*ProtonBufferWeight;
  }

  //change to analysis multiplicity intervals
  ProtonIntegratedYield[0] = (ProtonIntegratedYield[0] + ProtonIntegratedYield[1])/2.;
  ProtonIntegratedYield[1] = ProtonIntegratedYield[2];
  ProtonIntegratedYield[2] = ProtonIntegratedYield[3];
  ProtonIntegratedYield[3] = (ProtonIntegratedYield[4] + ProtonIntegratedYield[5] + ProtonIntegratedYield[6])/3.;
  //stat. unc.
  ProtonIntegratedYieldStatisticalUncertainties[0] = TMath::Sqrt(ProtonIntegratedYieldStatisticalUncertainties[0]*ProtonIntegratedYieldStatisticalUncertainties[0] + ProtonIntegratedYieldStatisticalUncertainties[1]*ProtonIntegratedYieldStatisticalUncertainties[1])/2.;
  ProtonIntegratedYieldStatisticalUncertainties[1] = ProtonIntegratedYieldStatisticalUncertainties[2];
  ProtonIntegratedYieldStatisticalUncertainties[2] = ProtonIntegratedYieldStatisticalUncertainties[3];
  ProtonIntegratedYieldStatisticalUncertainties[3] = TMath::Sqrt(ProtonIntegratedYieldStatisticalUncertainties[4]*ProtonIntegratedYieldStatisticalUncertainties[4] + ProtonIntegratedYieldStatisticalUncertainties[5]*ProtonIntegratedYieldStatisticalUncertainties[5] + ProtonIntegratedYieldStatisticalUncertainties[6]*ProtonIntegratedYieldStatisticalUncertainties[6])/3.;
  //syst. unc.
  ProtonIntegratedYieldSystematicUncertainties[0] = (ProtonIntegratedYieldSystematicUncertainties[0] + ProtonIntegratedYieldSystematicUncertainties[1])/2.;
  ProtonIntegratedYieldSystematicUncertainties[1] = ProtonIntegratedYieldSystematicUncertainties[2];
  ProtonIntegratedYieldSystematicUncertainties[2] = ProtonIntegratedYieldSystematicUncertainties[3];
  ProtonIntegratedYieldSystematicUncertainties[3] = (ProtonIntegratedYieldSystematicUncertainties[4] + ProtonIntegratedYieldSystematicUncertainties[5] + ProtonIntegratedYieldSystematicUncertainties[6])/3.;

  ProtonIntegratedYield[4] = ProtonBuffer/20.;
  ProtonIntegratedYieldStatisticalUncertainties[4] = TMath::Sqrt(ProtonBufferStat)/20.;
  ProtonIntegratedYieldSystematicUncertainties[4] = ProtonBufferSyst/20.;

  for (Int_t j=0; j <5; j++) {
    printf("Proton Yield %10e ± %10e (stat.) ± %10e (syst.) \n", ProtonIntegratedYield[j], ProtonIntegratedYieldStatisticalUncertainties[j], ProtonIntegratedYieldSystematicUncertainties[j]);
    // convert to relative uncertainty
    ProtonIntegratedYieldStatisticalUncertainties[j] = ProtonIntegratedYieldStatisticalUncertainties[j]/ProtonIntegratedYield[j];
    ProtonIntegratedYieldSystematicUncertainties[j] = ProtonIntegratedYieldSystematicUncertainties[j]/ProtonIntegratedYield[j];
  }
  printf("\n");

  // declare empty histograms
  TH1D* hMultInterval2;
  TH1D* hMultInterval2Syst;
  TH1D* hMultInterval3;
  TH1D* hMultInterval3Syst;
  TH1D* hMultInterval4;
  TH1D* hMultInterval4Syst;

  // Load data
  TFile MultIntervalfile (file1.Data());
  TH1D* hMultInterval = (TH1D*) MultIntervalfile.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hMultIntervalB = (TH1D*) MultIntervalfile.GetObjectChecked(Objectname2,"TH1D");
    hMultInterval->Add(hMultIntervalB);
    hMultInterval->Scale(0.5);
  }
  setOPT_hists(hMultInterval,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval->SetTitle("");
  hMultInterval->SetName("MultInterval");
  TH1D* hMultIntervalSyst = (TH1D*) MultIntervalfile.GetObjectChecked(Objectname+"Syst","TH1D");

  // Create histogram with only bin-by-bin uncorrelated uncertainties
  // Change to bin-by-bin uncorrelated syst (0-10%)
  Double_t BinByBinUncorSyst[] = {0.12,0.10,0.09,0.09,0.05,0.05};
  Double_t BinByBinUncorSystAnti[] = {0.05,0.04,0.04,0.04,0.05,0.05};
  TH1D* hMultIntervalSystUncor = (TH1D*) hMultIntervalSyst->Clone("hMultIntervalSystUncor");

  for (int i = 1; i <=hMultIntervalSyst->GetNbinsX(); i++) {
    hMultIntervalSystUncor->SetBinError(i, hMultIntervalSyst->GetBinContent(i)*BinByBinUncorSyst[i-1]); // correct syst histogram
  }

  if (average){
    TH1D* hMultIntervalSystA = (TH1D*) hMultIntervalSyst->Clone("hMultIntervalSystA");
    TH1D* hMultIntervalSystB = (TH1D*) MultIntervalfile.GetObjectChecked(Objectname2+"Syst","TH1D");

    TH1D* hMultIntervalSystUncorA = (TH1D*) hMultIntervalSystUncor->Clone("hMultIntervalSystUncorA");
    TH1D* hMultIntervalSystUncorB = (TH1D*) hMultIntervalSystB->Clone("hMultIntervalSystUncorB");

    for (int i = 1; i <=hMultIntervalSystB->GetNbinsX(); i++) {
      hMultIntervalSystUncorB->SetBinError(i, hMultIntervalSystB->GetBinContent(i)*BinByBinUncorSystAnti[i-1]);
    }

    hMultIntervalSyst->Add(hMultIntervalSystB);
    hMultIntervalSyst->Scale(0.5);

    Double_t SystUncA = 0;
    Double_t SystUncB = 0;
    if (CorrelatedUncertainties) {
      for (int i =0; i <= hMultIntervalSyst->GetNbinsX(); i++) {
        if (hMultIntervalSyst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultIntervalSystA->GetBinError(i);
        SystUncB = hMultIntervalSystB->GetBinError(i);
        hMultIntervalSyst->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }

    hMultIntervalSystUncor->Add(hMultIntervalSystUncorB);
    hMultIntervalSystUncor->Scale(0.5);
    if (CorrelatedUncertainties) {
      for (int i=0; i <= hMultIntervalSystUncor->GetNbinsX(); i++) {
        if (hMultIntervalSystUncor->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultIntervalSystUncorA->GetBinError(i);
        SystUncB = hMultIntervalSystUncorB->GetBinError(i);
        hMultIntervalSystUncor->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }
  }
  setOPT_hists(hMultIntervalSyst,pTLabel, yLable,510,20,0.5,kRed);
  hMultIntervalSyst->SetFillStyle(0);

  setOPT_hists(hMultIntervalSystUncor,pTLabel, yLable,510,20,0.5,kRed);
  hMultIntervalSystUncor->SetFillStyle(0);

  TFile Denominatorfile (file2.Data());
  TH1D* hIntegrated = (TH1D*) Denominatorfile.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hIntegratedB = (TH1D*) Denominatorfile.GetObjectChecked(Objectname2,"TH1D");
    hIntegrated->Add(hIntegratedB);
    hIntegrated->Scale(0.5);
  }
  setOPT_hists(hIntegrated,pTLabel, yLable,510,20,0.5,kRed);
  hIntegrated->SetTitle("");
  hIntegrated->SetName("Denominator");
  TH1D* hIntegratedSyst = (TH1D*) Denominatorfile.GetObjectChecked(Objectname+"Syst","TH1D");

  // Create histogram with only bin-by-bin uncorrelated uncertainties
  // Change to bin-by-bin uncorrelated syst (0-100%)
  Double_t BinByBinUncorSystMinBias[] = {0.12,0.10,0.09,0.09,0.05,0.05,0.05,0.05};
  Double_t BinByBinUncorSystMinBiasAnti[] = {0.05,0.04,0.04,0.04,0.05,0.05,0.05,0.05};
  TH1D* hIntegratedSystUncor = (TH1D*) hIntegratedSyst->Clone("hIntegratedSystUncor");

  for (int i = 1; i<=hIntegratedSyst->GetNbinsX(); i++) {
    hIntegratedSystUncor->SetBinError(i, hIntegratedSyst->GetBinContent(i)*BinByBinUncorSystMinBias[i-1]);
  }

  if (average){
    TH1D* hIntegratedSystA = (TH1D*) hIntegratedSyst->Clone("hIntegratedSystA");
    TH1D* hIntegratedSystB = (TH1D*) Denominatorfile.GetObjectChecked(Objectname2+"Syst","TH1D");

    TH1D* hIntegratedSystUncorA = (TH1D*) hIntegratedSystUncor->Clone("hIntegratedSystUncorA");
    TH1D* hIntegratedSystUncorB = (TH1D*) hIntegratedSystB->Clone("hIntegratedSystUncorB");

    for (int i = 1; i <=hIntegratedSystB->GetNbinsX(); i++) {
      hIntegratedSystUncorB->SetBinError(i, hIntegratedSystB->GetBinContent(i)*BinByBinUncorSystMinBiasAnti[i-1]);
    }

    hIntegratedSyst->Add(hIntegratedSystB);
    hIntegratedSyst->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hIntegratedSyst->GetNbinsX(); i++) {
        if (hIntegratedSyst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hIntegratedSystA->GetBinError(i);
        SystUncB = hIntegratedSystB->GetBinError(i);
        hIntegratedSyst->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }

    hIntegratedSystUncor->Add(hIntegratedSystUncorB);
    hIntegratedSystUncor->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hIntegratedSystUncor->GetNbinsX(); i++) {
        if (hIntegratedSystUncor->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hIntegratedSystUncorA->GetBinError(i);
        SystUncB = hIntegratedSystUncorB->GetBinError(i);
        hIntegratedSystUncor->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }

  }
  setOPT_hists(hIntegratedSyst,pTLabel, yLable,510,20,0.5,kRed);
  hIntegratedSyst->SetFillStyle(0);

  setOPT_hists(hIntegratedSystUncor,pTLabel, yLable,510,20,0.5,kRed);
  hIntegratedSystUncor->SetFillStyle(0);

  TFile MultInterval2file (file3.Data());
  hMultInterval2 = (TH1D*) MultInterval2file.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hMultInterval2B = (TH1D*) MultInterval2file.GetObjectChecked(Objectname2,"TH1D");
    hMultInterval2->Add(hMultInterval2B);
    hMultInterval2->Scale(0.5);
  }
  setOPT_hists(hMultInterval2,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval2->SetTitle("");
  hMultInterval2->SetName("MultInterval");
  hMultInterval2Syst = (TH1D*) MultInterval2file.GetObjectChecked(Objectname+"Syst","TH1D");

  // Change to bin-by-bin uncorrelated syst (10-20%)
  Double_t BinByBinUncorSyst2[] = {0.12,0.10,0.09,0.09,0.05};
  Double_t BinByBinUncorSyst2Anti[] = {0.05,0.04,0.04,0.04,0.05};

    TH1D* hMultInterval2SystUncor = (TH1D*) hMultInterval2Syst->Clone("hMultInterval2SystUncor");
  for (int i = 1; i <=hMultInterval2SystUncor->GetNbinsX(); i++) {
    hMultInterval2SystUncor->SetBinError(i, hMultInterval2Syst->GetBinContent(i)*BinByBinUncorSyst2[i-1]); // correct syst histogram
  }

  if (average){
    TH1D* hMultInterval2SystA = (TH1D*) hMultInterval2Syst->Clone("hMultInterval2SystA");
    TH1D* hMultInterval2SystB = (TH1D*) MultInterval2file.GetObjectChecked(Objectname2+"Syst","TH1D");

    TH1D* hMultInterval2SystUncorA = (TH1D*) hMultInterval2SystUncor->Clone("hMultInterval2SystUncorA");
    TH1D* hMultInterval2SystUncorB = (TH1D*) hMultInterval2SystB->Clone("hMultInterval2SystUncorB");

    for (int i = 1; i <=hMultInterval2SystUncorB->GetNbinsX(); i++) {
      hMultInterval2SystUncorB->SetBinError(i, hMultInterval2SystB->GetBinContent(i)*BinByBinUncorSyst2Anti[i-1]);
    }

    hMultInterval2Syst->Add(hMultInterval2SystB);
    hMultInterval2Syst->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hMultInterval2Syst->GetNbinsX(); i++) {
        if (hMultInterval2Syst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultInterval2SystA->GetBinError(i);
        SystUncB = hMultInterval2SystB->GetBinError(i);
        hMultInterval2Syst->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }

    hMultInterval2SystUncor->Add(hMultInterval2SystUncorB);
    hMultInterval2SystUncor->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hMultInterval2SystUncor->GetNbinsX(); i++) {
        if (hMultInterval2SystUncor->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultInterval2SystUncorA->GetBinError(i);
        SystUncB = hMultInterval2SystUncorB->GetBinError(i);
        hMultInterval2SystUncor->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }
  }
  setOPT_hists(hMultInterval2Syst,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval2Syst->SetFillStyle(0);

  setOPT_hists(hMultInterval2SystUncor,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval2SystUncor->SetFillStyle(0);

  TFile MultInterval3file (file4.Data());
  hMultInterval3 = (TH1D*) MultInterval3file.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hMultInterval3B = (TH1D*) MultInterval3file.GetObjectChecked(Objectname2,"TH1D");
    hMultInterval3->Add(hMultInterval3B);
    hMultInterval3->Scale(0.5);
  }
  setOPT_hists(hMultInterval3,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval3->SetTitle("");
  hMultInterval3->SetName("MultInterval");
  hMultInterval3Syst = (TH1D*) MultInterval3file.GetObjectChecked(Objectname+"Syst","TH1D");

  // Change to bin-by-bin uncorrelated syst (20-40%)
  Double_t BinByBinUncorSyst3[] = {0.12,0.10,0.09,0.09,0.05};
  Double_t BinByBinUncorSyst3Anti[] = {0.05,0.04,0.04,0.04,0.05};
  TH1D* hMultInterval3SystUncor = (TH1D*) hMultInterval3Syst->Clone("hMultInterval3SystUncor");
  for (int i = 1; i <=hMultInterval3SystUncor->GetNbinsX(); i++) {
    hMultInterval3SystUncor->SetBinError(i, hMultInterval3Syst->GetBinContent(i)*BinByBinUncorSyst3[i-1]); // correct syst histogram
  }

  if (average){
    TH1D* hMultInterval3SystA = (TH1D*) hMultInterval3Syst->Clone("hMultInterval3SystA");
    TH1D* hMultInterval3SystB = (TH1D*) MultInterval3file.GetObjectChecked(Objectname2+"Syst","TH1D");

    TH1D* hMultInterval3SystUncorA = (TH1D*) hMultInterval3SystUncor->Clone("hMultInterval3SystUncorA");
    TH1D* hMultInterval3SystUncorB = (TH1D*) hMultInterval3SystB->Clone("hMultInterval3SystUncorB");
    for (int i = 1; i <=hMultInterval3SystUncorB->GetNbinsX(); i++) {
      hMultInterval3SystUncorB->SetBinError(i, hMultInterval3SystB->GetBinContent(i)*BinByBinUncorSyst3Anti[i-1]);
    }

    hMultInterval3Syst->Add(hMultInterval3SystB);
    hMultInterval3Syst->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hMultInterval3Syst->GetNbinsX(); i++) {
        if (hMultInterval3Syst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultInterval3SystA->GetBinError(i);
        SystUncB = hMultInterval3SystB->GetBinError(i);
        hMultInterval3Syst->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }

    hMultInterval3SystUncor->Add(hMultInterval3SystUncorB);
    hMultInterval3SystUncor->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hMultInterval3SystUncor->GetNbinsX(); i++) {
        if (hMultInterval3SystUncor->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultInterval3SystUncorA->GetBinError(i);
        SystUncB = hMultInterval3SystUncorB->GetBinError(i);
        hMultInterval3SystUncor->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }
  }
  setOPT_hists(hMultInterval3Syst,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval3Syst->SetFillStyle(0);

  setOPT_hists(hMultInterval3SystUncor,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval3SystUncor->SetFillStyle(0);

  TFile MultInterval4file (file5.Data());
  hMultInterval4 = (TH1D*) MultInterval4file.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hMultInterval4B = (TH1D*) MultInterval4file.GetObjectChecked(Objectname2,"TH1D");
    hMultInterval4->Add(hMultInterval4B);
    hMultInterval4->Scale(0.5);
  }
  setOPT_hists(hMultInterval4,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval4->SetTitle("");
  hMultInterval4->SetName("MultInterval");
  hMultInterval4Syst = (TH1D*) MultInterval4file.GetObjectChecked(Objectname+"Syst","TH1D");

  // Change to bin-by-bin uncorrelated syst (20-40%)
  Double_t BinByBinUncorSyst4[] = {0.12,0.10,0.09,0.05};
  Double_t BinByBinUncorSyst4Anti[] = {0.05,0.04,0.04,0.05};

  TH1D* hMultInterval4SystUncor = (TH1D*) hMultInterval4Syst->Clone("hMultInterval4SystUncor");
  for (int i = 1; i <=hMultInterval4SystUncor->GetNbinsX(); i++) {
    hMultInterval4SystUncor->SetBinError(i, hMultInterval4Syst->GetBinContent(i)*BinByBinUncorSyst4[i-1]); // correct syst histogram
  }

  if (average){
    TH1D* hMultInterval4SystA = (TH1D*) hMultInterval4Syst->Clone("hMultInterval4SystA");
    TH1D* hMultInterval4SystB = (TH1D*) MultInterval4file.GetObjectChecked(Objectname2+"Syst","TH1D");

    TH1D* hMultInterval4SystUncorA = (TH1D*) hMultInterval4SystUncor->Clone("hMultInterval4SystUncorA");
    TH1D* hMultInterval4SystUncorB = (TH1D*) hMultInterval4SystB->Clone("hMultInterval4SystUncorB");

    for (int i = 1; i <=hMultInterval4SystUncorB->GetNbinsX(); i++) {
      hMultInterval4SystUncorB->SetBinError(i, hMultInterval4SystB->GetBinContent(i)*BinByBinUncorSyst4Anti[i-1]);
    }

    hMultInterval4Syst->Add(hMultInterval4SystB);
    hMultInterval4Syst->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hMultInterval4Syst->GetNbinsX(); i++) {
        if (hMultInterval4Syst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultInterval4SystA->GetBinError(i);
        SystUncB = hMultInterval4SystB->GetBinError(i);
        hMultInterval4Syst->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }

    hMultInterval4SystUncor->Add(hMultInterval4SystUncorB);
    hMultInterval4SystUncor->Scale(0.5);
    if (CorrelatedUncertainties) {
      Double_t SystUncA = 0;
      Double_t SystUncB = 0;
      for (int i =0; i <= hMultInterval4SystUncor->GetNbinsX(); i++) {
        if (hMultInterval4SystUncor->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hMultInterval4SystUncorA->GetBinError(i);
        SystUncB = hMultInterval4SystUncorB->GetBinError(i);
        hMultInterval4SystUncor->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }
  }

  setOPT_hists(hMultInterval4Syst,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval4Syst->SetFillStyle(0);

  setOPT_hists(hMultInterval4SystUncor,pTLabel, yLable,510,20,0.5,kRed);
  hMultInterval4SystUncor->SetFillStyle(0);

  // Empty bins outside the valid range (Important because YieldMean macro integrates full histogram even if the fit range is different)
  for (int jBin=1; jBin <= hIntegrated->GetNbinsX(); jBin++) {
    if (hIntegrated->GetBinCenter(jBin) > minPtMeasured && hIntegrated->GetBinCenter(jBin) < maxPtMeasured) continue;
    hIntegrated->SetBinContent(jBin,0);
    hIntegrated->SetBinError(jBin,0);
    hIntegratedSyst->SetBinContent(jBin,0);
    hIntegratedSyst->SetBinError(jBin,0);
    hIntegratedSystUncor->SetBinContent(jBin,0);
    hIntegratedSystUncor->SetBinError(jBin,0);
  }
  for (int jBin=1; jBin <= hMultInterval->GetNbinsX(); jBin++) {
    if (hMultInterval->GetBinCenter(jBin) > minPtMeasured && hMultInterval->GetBinCenter(jBin) < maxPtMeasured) continue;
    hMultInterval->SetBinContent(jBin,0);
    hMultInterval->SetBinError(jBin,0);
    hMultIntervalSyst->SetBinContent(jBin,0);
    hMultIntervalSyst->SetBinError(jBin,0);
    hMultIntervalSystUncor->SetBinContent(jBin,0);
    hMultIntervalSystUncor->SetBinError(jBin,0);
  }
  for (int jBin=1; jBin <= hMultInterval2->GetNbinsX(); jBin++) {
    if (hMultInterval2->GetBinCenter(jBin) > minPtMeasured && hMultInterval2->GetBinCenter(jBin) < maxPtMeasured) continue;
    hMultInterval2->SetBinContent(jBin,0);
    hMultInterval2->SetBinError(jBin,0);
    hMultInterval2Syst->SetBinContent(jBin,0);
    hMultInterval2Syst->SetBinError(jBin,0);
    hMultInterval2SystUncor->SetBinContent(jBin,0);
    hMultInterval2SystUncor->SetBinError(jBin,0);
  }
  for (int jBin=1; jBin <= hMultInterval3->GetNbinsX(); jBin++) {
    if (hMultInterval3->GetBinCenter(jBin) > minPtMeasured && hMultInterval3->GetBinCenter(jBin) < maxPtMeasured) continue;
    hMultInterval3->SetBinContent(jBin,0);
    hMultInterval3->SetBinError(jBin,0);
    hMultInterval3Syst->SetBinContent(jBin,0);
    hMultInterval3Syst->SetBinError(jBin,0);
    hMultInterval3SystUncor->SetBinContent(jBin,0);
    hMultInterval3SystUncor->SetBinError(jBin,0);
  }
  for (int jBin=1; jBin <= hMultInterval4->GetNbinsX(); jBin++) {
    if (hMultInterval4->GetBinCenter(jBin) > minPtMeasured && hMultInterval4->GetBinCenter(jBin) < maxPtMeasured) continue;
    hMultInterval4->SetBinContent(jBin,0);
    hMultInterval4->SetBinError(jBin,0);
    hMultInterval4Syst->SetBinContent(jBin,0);
    hMultInterval4Syst->SetBinError(jBin,0);
    hMultInterval4SystUncor->SetBinContent(jBin,0);
    hMultInterval4SystUncor->SetBinError(jBin,0);
  }

  hIntegrated->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hIntegratedSyst->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hIntegratedSystUncor->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);

  hMultInterval->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultIntervalSyst->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultIntervalSystUncor->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);

  hMultInterval2->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultInterval2Syst->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultInterval2SystUncor->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);

  hMultInterval3->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultInterval3Syst->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultInterval3SystUncor->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);

  hMultInterval4->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultInterval4Syst->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);
  hMultInterval4SystUncor->GetXaxis()->SetRangeUser(minPtMeasured,maxPtMeasured);

  TH1D* Axis = new TH1D("", "", 100,minpT,maxpT);
  setOPT_hists(Axis,pTLabel, yLable,510,20,0.5);
  Axis->GetYaxis()->SetRangeUser(5e-9,5e-4);
  Axis->GetXaxis()->SetLabelSize(0.035);
  Axis->GetYaxis()->SetLabelSize(0.035);

  TH1* hFitResults[NumberOfFitFunctions];
  Double_t IntegratedYield[NumberOfFitFunctions];
  Double_t IntegratedYieldStat[NumberOfFitFunctions];
  Double_t IntegratedYieldSystCor[NumberOfFitFunctions];
  Double_t IntegratedYieldSystUnCor[NumberOfFitFunctions];

  Double_t AverageIntegratedYield = 0;
  Double_t AverageIntegratedYieldStat = 0;
  Double_t AverageIntegratedYieldSyst = 0;
  Double_t AverageIntegratedYieldSpreadSyst = 0;

  Double_t MeanPt[NumberOfFitFunctions];
  Double_t MeanPtStat[NumberOfFitFunctions];
  Double_t MeanPtSystCor[NumberOfFitFunctions];
  Double_t MeanPtSystUncor[NumberOfFitFunctions];

  TGraphErrors* gMeanPt = new TGraphErrors();
  TGraphErrors* gMeanPtSyst = new TGraphErrors();

  Double_t AverageMeanPt = 0;
  Double_t AverageMeanPtStat = 0;
  Double_t AverageMeanPtSyst = 0;
  Double_t AverageMeanPtSpreadSyst = 0;

  for (int i=0; i<NumberOfFitFunctions; i++) {
    hFitResults[i]=YieldMean(hMultInterval, hMultIntervalSyst, hMultIntervalSystUncor, fFitFunctions[i], minpT, maxpT, 0.01, 0.1, "0QI","log.root",minPtMeasured,maxPtMeasured);

    // Integrated yield
    IntegratedYield[i] = hFitResults[i]->GetBinContent(1);
    IntegratedYieldStat[i] = hFitResults[i]->GetBinContent(2);
    IntegratedYieldSystUnCor[i] = hFitResults[i]->GetBinContent(3);
    IntegratedYieldSystCor[i] = hFitResults[i]->GetBinContent(4);
    // Mean pT
    MeanPt[i] = hFitResults[i]->GetBinContent(5);
    MeanPtStat[i] = hFitResults[i]->GetBinContent(6);
    MeanPtSystUncor[i] = hFitResults[i]->GetBinContent(7);
    MeanPtSystCor[i] = hFitResults[i]->GetBinContent(8);

  }

  AverageMeanPt = 0;
  AverageMeanPtStat = 0;
  AverageMeanPtSyst = 0;
  AverageIntegratedYield = 0;
  AverageIntegratedYieldStat = 0;
  AverageIntegratedYieldSyst = 0;
  for (int i=0; i<NumberOfFitFunctions; i++) {
    if (i == 5) continue;
    if (i == 6) continue;

    AverageIntegratedYield += IntegratedYield[i];
    AverageIntegratedYieldStat += IntegratedYieldStat[i];
    AverageIntegratedYieldSyst += TMath::Sqrt(IntegratedYieldSystUnCor[i]*IntegratedYieldSystUnCor[i]+IntegratedYieldSystCor[i]*IntegratedYieldSystCor[i]);

    AverageMeanPt += MeanPt[i];
    AverageMeanPtStat += MeanPtStat[i];
    AverageMeanPtSyst += TMath::Sqrt(MeanPtSystUncor[i]*MeanPtSystUncor[i]+MeanPtSystCor[i]*MeanPtSystCor[i]);
  }

  AverageIntegratedYield = AverageIntegratedYield/5.;
  AverageIntegratedYieldStat = AverageIntegratedYieldStat/5.;
  AverageIntegratedYieldSyst = AverageIntegratedYieldSyst/5.;
  AverageIntegratedYieldSpreadSyst = (GetMaximum(IntegratedYield,NumberOfFitFunctions)-GetMinimum(IntegratedYield,NumberOfFitFunctions))/TMath::Sqrt(12.);

  AverageMeanPt = AverageMeanPt/5.;
  AverageMeanPtStat = AverageMeanPtStat/5.;
  AverageMeanPtSyst = AverageMeanPtSyst/5.;
  AverageMeanPtSpreadSyst = (GetMaximum(MeanPt,NumberOfFitFunctions)-GetMinimum(MeanPt,NumberOfFitFunctions))/TMath::Sqrt(12.);

  printf("\n");
  printf("Multiplicity 0 - 10 %%\n");
  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("<p_{T}> = %.2f ± %.2f ± %.2f ± %.2f(%s) \n", MeanPt[i], MeanPtStat[i], MeanPtSystUncor[i], MeanPtSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (syst.) ± %.2f (spread) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, AverageMeanPtSyst, AverageMeanPtSpreadSyst);
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (total syst.) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, sqrt(pow(AverageMeanPtSyst,2)+pow(AverageMeanPtSpreadSyst,2)));
  printf("\n");

  // Fill Graph for mean pT vs multiplicity
  gMeanPt->SetPoint(gMeanPt->GetN(), MeanMultiplicity[0],AverageMeanPt);
  gMeanPt->SetPointError(gMeanPt->GetN()-1,MultiplicityUncertainty[0],AverageMeanPtStat);

  AverageMeanPtSyst = TMath::Sqrt(AverageMeanPtSyst*AverageMeanPtSyst + AverageMeanPtSpreadSyst*AverageMeanPtSpreadSyst);

  gMeanPtSyst->SetPoint(gMeanPtSyst->GetN(), MeanMultiplicity[0],AverageMeanPt);
  gMeanPtSyst->SetPointError(gMeanPtSyst->GetN()-1,MultiplicityUncertainty[0],AverageMeanPtSyst);

  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("dN/dy = %10e ± %10e ± %10e ± %10e (%s)\n",  IntegratedYield[i], IntegratedYieldStat[i], IntegratedYieldSystUnCor[i], IntegratedYieldSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (syst.) ± %10e (spread) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, AverageIntegratedYieldSyst, AverageIntegratedYieldSpreadSyst);
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (total syst.) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst));


  AverageIntegratedYieldStat = AverageIntegratedYieldStat/AverageIntegratedYield;
  AverageIntegratedYieldSyst = TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst)/AverageIntegratedYield;

  Double_t RatioToProtonValue = 0;
  Double_t RatioToProtonStat = 0;
  Double_t RatioToProtonSyst = 0;

  TGraphErrors* gRatioToProton = new TGraphErrors();
  TGraphErrors* gRatioToProtonSyst = new TGraphErrors();

  RatioToProtonValue = 2*AverageIntegratedYield/ProtonIntegratedYield[0];
  RatioToProtonStat = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldStat*AverageIntegratedYieldStat + ProtonIntegratedYieldStatisticalUncertainties[0]*ProtonIntegratedYieldStatisticalUncertainties[0]);

  RatioToProtonSyst = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + ProtonIntegratedYieldSystematicUncertainties[0]*ProtonIntegratedYieldSystematicUncertainties[0]);

  gRatioToProton->SetPoint(gRatioToProton->GetN(), MeanMultiplicity[0],RatioToProtonValue);
  gRatioToProton->SetPointError(gRatioToProton->GetN()-1,MultiplicityUncertainty[0],RatioToProtonStat);

  gRatioToProtonSyst->SetPoint(gRatioToProtonSyst->GetN(), MeanMultiplicity[0],RatioToProtonValue);
  gRatioToProtonSyst->SetPointError(gRatioToProtonSyst->GetN()-1,MultiplicityUncertainty[0],RatioToProtonSyst);

  printf("\n");
  printf("%10e ± %10e (stat.) ± %10e (syst) (Ratio to proton)\n", RatioToProtonValue, RatioToProtonStat,RatioToProtonSyst);

  for (int i=0; i<NumberOfFitFunctions; i++) {
    hFitResults[i]=YieldMean(hMultInterval2, hMultInterval2Syst, hMultInterval2SystUncor, fFitFunctions[i], minpT, maxpT, 0.01, 0.1, "0QI","log.root",minPtMeasured,maxPtMeasured);
    // Integrated yield
    IntegratedYield[i] = hFitResults[i]->GetBinContent(1);
    IntegratedYieldStat[i] = hFitResults[i]->GetBinContent(2);
    IntegratedYieldSystUnCor[i] = hFitResults[i]->GetBinContent(3);
    IntegratedYieldSystCor[i] = hFitResults[i]->GetBinContent(4);
    // Mean pT
    MeanPt[i] = hFitResults[i]->GetBinContent(5);
    MeanPtStat[i] = hFitResults[i]->GetBinContent(6);
    MeanPtSystUncor[i] = hFitResults[i]->GetBinContent(7);
    MeanPtSystCor[i] = hFitResults[i]->GetBinContent(8);
  }

  AverageMeanPt = 0;
  AverageMeanPtStat = 0;
  AverageMeanPtSyst = 0;
  AverageIntegratedYield = 0;
  AverageIntegratedYieldStat = 0;
  AverageIntegratedYieldSyst = 0;
  for (int i=0; i<NumberOfFitFunctions; i++) {
    if (i == 5) continue;
    if (i == 6) continue;
    AverageIntegratedYield += IntegratedYield[i];
    AverageIntegratedYieldStat += IntegratedYieldStat[i];
    AverageIntegratedYieldSyst += TMath::Sqrt(IntegratedYieldSystUnCor[i]*IntegratedYieldSystUnCor[i]+IntegratedYieldSystCor[i]*IntegratedYieldSystCor[i]);

    AverageMeanPt += MeanPt[i];
    AverageMeanPtStat += MeanPtStat[i];
    AverageMeanPtSyst += TMath::Sqrt( MeanPtSystUncor[i]*MeanPtSystUncor[i]+MeanPtSystCor[i]*MeanPtSystCor[i] );
  }
  AverageIntegratedYield = AverageIntegratedYield/5.;
  AverageIntegratedYieldStat = AverageIntegratedYieldStat/5.;
  AverageIntegratedYieldSyst = AverageIntegratedYieldSyst/5.;
  AverageIntegratedYieldSpreadSyst = (GetMaximum(IntegratedYield,NumberOfFitFunctions)-GetMinimum(IntegratedYield,NumberOfFitFunctions))/TMath::Sqrt(12.);

  AverageMeanPt = AverageMeanPt/5.;
  AverageMeanPtStat = AverageMeanPtStat/5.;
  AverageMeanPtSyst = AverageMeanPtSyst/5.;
  AverageMeanPtSpreadSyst = (GetMaximum(MeanPt,NumberOfFitFunctions)-GetMinimum(MeanPt,NumberOfFitFunctions))/TMath::Sqrt(12.);

  printf("\n");
  printf("Multiplicity 10 - 20 %%\n");
  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("<p_{T}> = %.2f ± %.2f ± %.2f ± %.2f(%s) \n", MeanPt[i], MeanPtStat[i], MeanPtSystUncor[i], MeanPtSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (syst.) ± %.2f (spread) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, AverageMeanPtSyst, AverageMeanPtSpreadSyst);
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (total syst.) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, sqrt(pow(AverageMeanPtSyst,2)+pow(AverageMeanPtSpreadSyst,2)));
  printf("\n");

  // Fill Graph for mean pT vs multiplicity
  gMeanPt->SetPoint(gMeanPt->GetN(), MeanMultiplicity[1],AverageMeanPt);
  gMeanPt->SetPointError(gMeanPt->GetN()-1,MultiplicityUncertainty[1],AverageMeanPtStat);

  AverageMeanPtSyst = TMath::Sqrt(AverageMeanPtSyst*AverageMeanPtSyst + AverageMeanPtSpreadSyst*AverageMeanPtSpreadSyst);

  gMeanPtSyst->SetPoint(gMeanPtSyst->GetN(), MeanMultiplicity[1],AverageMeanPt);
  gMeanPtSyst->SetPointError(gMeanPtSyst->GetN()-1,MultiplicityUncertainty[1],AverageMeanPtSyst);

  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("dN/dy = %10e ± %10e ± %10e ± %10e (%s)\n",  IntegratedYield[i], IntegratedYieldStat[i], IntegratedYieldSystUnCor[i], IntegratedYieldSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (syst.) ± %10e (spread) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, AverageIntegratedYieldSyst, AverageIntegratedYieldSpreadSyst);
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (total syst.) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst));

  AverageIntegratedYieldStat = AverageIntegratedYieldStat/AverageIntegratedYield;
  AverageIntegratedYieldSyst = TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst)/AverageIntegratedYield;

  RatioToProtonValue = 0;
  RatioToProtonStat = 0;
  RatioToProtonSyst = 0;

  RatioToProtonValue = 2*AverageIntegratedYield/ProtonIntegratedYield[1];
  RatioToProtonStat = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldStat*AverageIntegratedYieldStat + ProtonIntegratedYieldStatisticalUncertainties[1]*ProtonIntegratedYieldStatisticalUncertainties[1]);

  RatioToProtonSyst = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + ProtonIntegratedYieldSystematicUncertainties[1]*ProtonIntegratedYieldSystematicUncertainties[1]);

  gRatioToProton->SetPoint(gRatioToProton->GetN(), MeanMultiplicity[1],RatioToProtonValue);
  gRatioToProton->SetPointError(gRatioToProton->GetN()-1,MultiplicityUncertainty[1],RatioToProtonStat);

  gRatioToProtonSyst->SetPoint(gRatioToProtonSyst->GetN(), MeanMultiplicity[1],RatioToProtonValue);
  gRatioToProtonSyst->SetPointError(gRatioToProtonSyst->GetN()-1,MultiplicityUncertainty[1],RatioToProtonSyst);

  printf("\n");
  printf("%10e ± %10e (stat.) ± %10e (syst) (Ratio to proton)\n", RatioToProtonValue, RatioToProtonStat,RatioToProtonSyst);

  for (int i=0; i<NumberOfFitFunctions; i++) {
    hFitResults[i]=YieldMean(hMultInterval3, hMultInterval3Syst, hMultInterval3SystUncor, fFitFunctions[i], minpT, maxpT, 0.01, 0.1, "0QI","log.root",minPtMeasured,maxPtMeasured);
    // Integrated yield
    IntegratedYield[i] = hFitResults[i]->GetBinContent(1);
    IntegratedYieldStat[i] = hFitResults[i]->GetBinContent(2);
    IntegratedYieldSystUnCor[i] = hFitResults[i]->GetBinContent(3);
    IntegratedYieldSystCor[i] = hFitResults[i]->GetBinContent(4);
    // Mean pT
    MeanPt[i] = hFitResults[i]->GetBinContent(5);
    MeanPtStat[i] = hFitResults[i]->GetBinContent(6);
    MeanPtSystUncor[i] = hFitResults[i]->GetBinContent(7);
    MeanPtSystCor[i] = hFitResults[i]->GetBinContent(8);
  }

  AverageMeanPt = 0;
  AverageMeanPtStat = 0;
  AverageMeanPtSyst = 0;
  AverageIntegratedYield = 0;
  AverageIntegratedYieldStat = 0;
  AverageIntegratedYieldSyst = 0;
  for (int i=0; i<NumberOfFitFunctions; i++) {
    if (i == 5) continue;
    if (i == 6) continue;
    AverageIntegratedYield += IntegratedYield[i];
    AverageIntegratedYieldStat += IntegratedYieldStat[i];
    AverageIntegratedYieldSyst += TMath::Sqrt(IntegratedYieldSystUnCor[i]*IntegratedYieldSystUnCor[i]+IntegratedYieldSystCor[i]*IntegratedYieldSystCor[i]);

    AverageMeanPt += MeanPt[i];
    AverageMeanPtStat += MeanPtStat[i];
    AverageMeanPtSyst += TMath::Sqrt(MeanPtSystUncor[i]*MeanPtSystUncor[i]+MeanPtSystCor[i]*MeanPtSystCor[i]);
  }
  AverageIntegratedYield = AverageIntegratedYield/5.;
  AverageIntegratedYieldStat = AverageIntegratedYieldStat/5.;
  AverageIntegratedYieldSyst = AverageIntegratedYieldSyst/5.;
  AverageIntegratedYieldSpreadSyst = (GetMaximum(IntegratedYield,NumberOfFitFunctions)-GetMinimum(IntegratedYield,NumberOfFitFunctions))/TMath::Sqrt(12.);

  AverageMeanPt = AverageMeanPt/5.;
  AverageMeanPtStat = AverageMeanPtStat/5.;
  AverageMeanPtSyst = AverageMeanPtSyst/5.;
  AverageMeanPtSpreadSyst = (GetMaximum(MeanPt,NumberOfFitFunctions)-GetMinimum(MeanPt,NumberOfFitFunctions))/TMath::Sqrt(12.);

  printf("\n");
  printf("Multiplicity 20 - 40 %%\n");
  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("<p_{T}> = %.2f ± %.2f ± %.2f ± %.2f(%s) \n", MeanPt[i], MeanPtStat[i], MeanPtSystUncor[i], MeanPtSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (syst.) ± %.2f (spread) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, AverageMeanPtSyst, AverageMeanPtSpreadSyst);
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (total syst.) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, sqrt(pow(AverageMeanPtSyst,2)+pow(AverageMeanPtSpreadSyst,2)));
  printf("\n");

  // Fill Graph for mean pT vs multiplicity
  gMeanPt->SetPoint(gMeanPt->GetN(), MeanMultiplicity[2],AverageMeanPt);
  gMeanPt->SetPointError(gMeanPt->GetN()-1,MultiplicityUncertainty[2],AverageMeanPtStat);

  AverageMeanPtSyst = TMath::Sqrt(AverageMeanPtSyst*AverageMeanPtSyst + AverageMeanPtSpreadSyst*AverageMeanPtSpreadSyst);

  gMeanPtSyst->SetPoint(gMeanPtSyst->GetN(), MeanMultiplicity[2],AverageMeanPt);
  gMeanPtSyst->SetPointError(gMeanPtSyst->GetN()-1,MultiplicityUncertainty[2],AverageMeanPtSyst);

  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("dN/dy = %10e ± %10e ± %10e ± %10e (%s)\n",  IntegratedYield[i], IntegratedYieldStat[i], IntegratedYieldSystUnCor[i], IntegratedYieldSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (syst.) ± %10e (spread) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, AverageIntegratedYieldSyst, AverageIntegratedYieldSpreadSyst);
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (total syst.) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst));


  AverageIntegratedYieldStat = AverageIntegratedYieldStat/AverageIntegratedYield;

  AverageIntegratedYieldSyst = TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst)/AverageIntegratedYield; // 0.114 - relative uncertainty evaluated from 20-100%

  RatioToProtonValue = 0;
  RatioToProtonStat = 0;
  RatioToProtonSyst = 0;

  RatioToProtonValue = 2*AverageIntegratedYield/ProtonIntegratedYield[2];
  RatioToProtonStat = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldStat*AverageIntegratedYieldStat + ProtonIntegratedYieldStatisticalUncertainties[2]*ProtonIntegratedYieldStatisticalUncertainties[2]);

  RatioToProtonSyst = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + ProtonIntegratedYieldSystematicUncertainties[2]*ProtonIntegratedYieldSystematicUncertainties[2]);

  gRatioToProton->SetPoint(gRatioToProton->GetN(), MeanMultiplicity[2],RatioToProtonValue);
  gRatioToProton->SetPointError(gRatioToProton->GetN()-1,MultiplicityUncertainty[2],RatioToProtonStat);

  gRatioToProtonSyst->SetPoint(gRatioToProtonSyst->GetN(), MeanMultiplicity[2],RatioToProtonValue);
  gRatioToProtonSyst->SetPointError(gRatioToProtonSyst->GetN()-1,MultiplicityUncertainty[2],RatioToProtonSyst);

  printf("\n");
  printf("%10e ± %10e (stat.) ± %10e (syst) (Ratio to proton)\n", RatioToProtonValue, RatioToProtonStat,RatioToProtonSyst);

  printf("Parameter: #beta=%f, n=%f, norm=%f\n",fFitFunctions[5]->GetParameter(1),fFitFunctions[5]->GetParameter(3), fFitFunctions[5]->GetParameter(4));

  Double_t StatUnc, SystUnc;
  // Fix to fit to 20-40%
  TH1D* hMultInterval3TotalUnc = (TH1D*) hMultInterval3Syst->Clone("hMultIntervalTotalUnc");
  for (int i=0; i <= hMultInterval3TotalUnc->GetNbinsX(); i++) {
    StatUnc = hMultInterval3->GetBinError(i);
    SystUnc = hMultInterval3Syst->GetBinError(i);
    hMultInterval3TotalUnc->SetBinError(i,TMath::Sqrt(StatUnc*StatUnc+SystUnc*SystUnc));
  }

  hMultInterval3TotalUnc->Fit(fFitFunctions[5], "0QI","",minPtMeasured,maxPtMeasured);
  hMultInterval3TotalUnc->Fit(fFitFunctions[6], "0QI","",minPtMeasured,maxPtMeasured);
  fFitFunctions[5]->FixParameter(3, fFitFunctions[5]->GetParameter(3)); // BW
  fFitFunctions[6]->FixParameter(1, fFitFunctions[6]->GetParameter(1)); // Tsallis

  for (int i=0; i<NumberOfFitFunctions; i++) {
    hFitResults[i]=YieldMean(hMultInterval4, hMultInterval4Syst, hMultInterval4SystUncor, fFitFunctions[i], minpT, maxpT, 0.01, 0.1, "0QI","log.root",minPtMeasured,maxPtMeasured);
    // Integrated yield
    IntegratedYield[i] = hFitResults[i]->GetBinContent(1);
    IntegratedYieldStat[i] = hFitResults[i]->GetBinContent(2);
    IntegratedYieldSystUnCor[i] = hFitResults[i]->GetBinContent(3);
    IntegratedYieldSystCor[i] = hFitResults[i]->GetBinContent(4);
    // Mean pT
    MeanPt[i] = hFitResults[i]->GetBinContent(5);
    MeanPtStat[i] = hFitResults[i]->GetBinContent(6);
    MeanPtSystUncor[i] = hFitResults[i]->GetBinContent(7);
    MeanPtSystCor[i] = hFitResults[i]->GetBinContent(8);
  }

  AverageMeanPt = 0;
  AverageMeanPtStat = 0;
  AverageMeanPtSyst = 0;
  AverageIntegratedYield = 0;
  AverageIntegratedYieldStat = 0;
  AverageIntegratedYieldSyst = 0;
  for (int i=0; i<NumberOfFitFunctions; i++) {
    if (i == 5) continue;
    if (i == 6) continue;
    AverageIntegratedYield += IntegratedYield[i];
    AverageIntegratedYieldStat += IntegratedYieldStat[i];
    AverageIntegratedYieldSyst += TMath::Sqrt(IntegratedYieldSystUnCor[i]*IntegratedYieldSystUnCor[i]+IntegratedYieldSystCor[i]*IntegratedYieldSystCor[i]);

    AverageMeanPt += MeanPt[i];
    AverageMeanPtStat += MeanPtStat[i];
    AverageMeanPtSyst += TMath::Sqrt(MeanPtSystUncor[i]*MeanPtSystUncor[i]+MeanPtSystCor[i]*MeanPtSystCor[i]);
  }
  AverageIntegratedYield = AverageIntegratedYield/5.;
  AverageIntegratedYieldStat = AverageIntegratedYieldStat/5.;
  AverageIntegratedYieldSyst = AverageIntegratedYieldSyst/5.;
  if (NumberOfFitFunctions > 6) IntegratedYield[6] = IntegratedYield[4]; // Override Tsallis result to remove it from the spread uncertainty
  AverageIntegratedYieldSpreadSyst = (GetMaximum(IntegratedYield,NumberOfFitFunctions)-GetMinimum(IntegratedYield,NumberOfFitFunctions))/TMath::Sqrt(12.);

  AverageMeanPt = AverageMeanPt/5.;
  AverageMeanPtStat = AverageMeanPtStat/5.;
  AverageMeanPtSyst = AverageMeanPtSyst/5.;
  if (NumberOfFitFunctions > 6) MeanPt[6] = MeanPt[4]; // Override Tsallis result to remove it from the spread uncertainty
  AverageMeanPtSpreadSyst = (GetMaximum(MeanPt,NumberOfFitFunctions)-GetMinimum(MeanPt,NumberOfFitFunctions))/TMath::Sqrt(12.);

  printf("\n");
  printf("Multiplicity 40 - 100 %%\n");
  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("<p_{T}> = %.2f ± %.2f ± %.2f ± %.2f(%s) \n", MeanPt[i], MeanPtStat[i], MeanPtSystUncor[i], MeanPtSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (syst.) ± %.2f (spread) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, AverageMeanPtSyst, AverageMeanPtSpreadSyst);
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (total syst.) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, sqrt(pow(AverageMeanPtSyst,2)+pow(AverageMeanPtSpreadSyst,2)));
  printf("\n");

  // Fill Graph for mean pT vs multiplicity
  gMeanPt->SetPoint(gMeanPt->GetN(), MeanMultiplicity[3],AverageMeanPt);
  gMeanPt->SetPointError(gMeanPt->GetN()-1,MultiplicityUncertainty[3],AverageMeanPtStat);

  AverageMeanPtSyst = TMath::Sqrt(AverageMeanPtSyst*AverageMeanPtSyst + AverageMeanPtSpreadSyst*AverageMeanPtSpreadSyst);

  gMeanPtSyst->SetPoint(gMeanPtSyst->GetN(), MeanMultiplicity[3],AverageMeanPt);
  gMeanPtSyst->SetPointError(gMeanPtSyst->GetN()-1,MultiplicityUncertainty[3],AverageMeanPtSyst);

  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("dN/dy = %10e ± %10e ± %10e ± %10e (%s)\n",  IntegratedYield[i], IntegratedYieldStat[i], IntegratedYieldSystUnCor[i], IntegratedYieldSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (syst.) ± %10e (spread) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, AverageIntegratedYieldSyst, AverageIntegratedYieldSpreadSyst);

  printf("dN/dy = %10e ± %10e (stat.) ± %10e (total syst.) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst));

  AverageIntegratedYieldStat = AverageIntegratedYieldStat/AverageIntegratedYield;
  AverageIntegratedYieldSyst = TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst)/AverageIntegratedYield;

  RatioToProtonValue = 0;
  RatioToProtonStat = 0;
  RatioToProtonSyst = 0;

  RatioToProtonValue = 2*AverageIntegratedYield/ProtonIntegratedYield[3];
  RatioToProtonStat = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldStat*AverageIntegratedYieldStat + ProtonIntegratedYieldStatisticalUncertainties[3]*ProtonIntegratedYieldStatisticalUncertainties[3]);

  RatioToProtonSyst = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + ProtonIntegratedYieldSystematicUncertainties[3]*ProtonIntegratedYieldSystematicUncertainties[3]);

  gRatioToProton->SetPoint(gRatioToProton->GetN(), MeanMultiplicity[3],RatioToProtonValue);
  gRatioToProton->SetPointError(gRatioToProton->GetN()-1,MultiplicityUncertainty[3],RatioToProtonStat);

  gRatioToProtonSyst->SetPoint(gRatioToProtonSyst->GetN(), MeanMultiplicity[3],RatioToProtonValue);
  gRatioToProtonSyst->SetPointError(gRatioToProtonSyst->GetN()-1,MultiplicityUncertainty[3],RatioToProtonSyst);

  printf("\n");
  printf("%10e ± %10e (stat.) ± %10e (syst) (Ratio to proton)\n", RatioToProtonValue, RatioToProtonStat,RatioToProtonSyst);

  // BW
  fFitFunctions[5]->ReleaseParameter(3); // unfix n
  fFitFunctions[5]->SetParLimits(3, 1, 5); // set limit for n

  //Tsallis
  fFitFunctions[6]->ReleaseParameter(1); // unfix n
  fFitFunctions[6]->SetParLimits(1, 0., 20); // set limit for n

  for (int i=0; i<NumberOfFitFunctions; i++) {
    hFitResults[i]=YieldMean(hIntegrated, hIntegratedSyst, hIntegratedSystUncor, fFitFunctions[i], minpT, maxpT, 0.01, 0.1,"0QI","log.root",minPtMeasured,maxPtMeasured);
    // Integrated yield
    IntegratedYield[i] = hFitResults[i]->GetBinContent(1);
    IntegratedYieldStat[i] = hFitResults[i]->GetBinContent(2);
    IntegratedYieldSystUnCor[i] = hFitResults[i]->GetBinContent(3);
    IntegratedYieldSystCor[i] = hFitResults[i]->GetBinContent(4);
    // Mean pT
    MeanPt[i] = hFitResults[i]->GetBinContent(5);
    MeanPtStat[i] = hFitResults[i]->GetBinContent(6);
    MeanPtSystUncor[i] = hFitResults[i]->GetBinContent(7);
    MeanPtSystCor[i] = hFitResults[i]->GetBinContent(8);
  }

  AverageMeanPt = 0;
  AverageMeanPtStat = 0;
  AverageMeanPtSyst = 0;
  AverageIntegratedYield = 0;
  AverageIntegratedYieldStat = 0;
  AverageIntegratedYieldSyst = 0;
  for (int i=0; i<NumberOfFitFunctions; i++) {
    if (i == 5) continue;
    if (i == 6) continue;
    AverageIntegratedYield += IntegratedYield[i];
    AverageIntegratedYieldStat += IntegratedYieldStat[i];
    AverageIntegratedYieldSyst += TMath::Sqrt(IntegratedYieldSystUnCor[i]*IntegratedYieldSystUnCor[i]+IntegratedYieldSystCor[i]*IntegratedYieldSystCor[i]);

    AverageMeanPt += MeanPt[i];
    AverageMeanPtStat += MeanPtStat[i];
    AverageMeanPtSyst += TMath::Sqrt(MeanPtSystUncor[i]*MeanPtSystUncor[i]+MeanPtSystCor[i]*MeanPtSystCor[i]);
  }
  AverageIntegratedYield = AverageIntegratedYield/5.;
  AverageIntegratedYieldStat = AverageIntegratedYieldStat/5.;
  AverageIntegratedYieldSyst = AverageIntegratedYieldSyst/5.;
  AverageIntegratedYieldSpreadSyst = (GetMaximum(IntegratedYield,NumberOfFitFunctions)-GetMinimum(IntegratedYield,NumberOfFitFunctions))/TMath::Sqrt(12.);

  AverageMeanPt = AverageMeanPt/5.;
  AverageMeanPtStat = AverageMeanPtStat/5.;
  AverageMeanPtSyst = AverageMeanPtSyst/5.;
  AverageMeanPtSpreadSyst = (GetMaximum(MeanPt,NumberOfFitFunctions)-GetMinimum(MeanPt,NumberOfFitFunctions))/TMath::Sqrt(12.);

  printf("\n");
  printf("Multiplicity integrated\n");
  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("<p_{T}> = %.2f ± %.2f ± %.2f ± %.2f(%s) \n", MeanPt[i], MeanPtStat[i], MeanPtSystUncor[i], MeanPtSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (syst.) ± %.2f (spread) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, AverageMeanPtSyst, AverageMeanPtSpreadSyst);
  printf("<p_{T}> = %.2f ± %.2f (stat.) ± %.2f (total syst.) (Average mean pT)\n", AverageMeanPt, AverageMeanPtStat, sqrt(pow(AverageMeanPtSyst,2)+pow(AverageMeanPtSpreadSyst,2)));

  printf("\n");
  for (int i=0; i<NumberOfFitFunctions; i++) {
    printf("dN/dy = %10e ± %10e ± %10e ± %10e (%s)\n",  IntegratedYield[i], IntegratedYieldStat[i], IntegratedYieldSystUnCor[i], IntegratedYieldSystCor[i], fFitFunctions[i]->GetName());
  }
  printf("dN/dy = %10e ± %10e (stat.) ± %10e (syst.) ± %10e (spread) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, AverageIntegratedYieldSyst, AverageIntegratedYieldSpreadSyst);

  printf("dN/dy = %10e ± %10e (stat.) ± %10e (total syst.) (Average integrated yield)\n", AverageIntegratedYield, AverageIntegratedYieldStat, TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst));

  RatioToProtonValue = 0;
  RatioToProtonStat = 0;
  RatioToProtonSyst = 0;

  AverageIntegratedYieldStat = AverageIntegratedYieldStat/AverageIntegratedYield;
  AverageIntegratedYieldSyst = TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + AverageIntegratedYieldSpreadSyst*AverageIntegratedYieldSpreadSyst)/AverageIntegratedYield;

  RatioToProtonValue = 2*AverageIntegratedYield/ProtonIntegratedYield[4];
  RatioToProtonStat = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldStat*AverageIntegratedYieldStat + ProtonIntegratedYieldStatisticalUncertainties[4]*ProtonIntegratedYieldStatisticalUncertainties[4]);

  RatioToProtonSyst = RatioToProtonValue*TMath::Sqrt(AverageIntegratedYieldSyst*AverageIntegratedYieldSyst + ProtonIntegratedYieldSystematicUncertainties[4]*ProtonIntegratedYieldSystematicUncertainties[4]);

  printf("\n");
  printf("%10e ± %10e (stat.) ± %10e (syst) (Ratio to proton)\n", RatioToProtonValue, RatioToProtonStat,RatioToProtonSyst);

  //Ratio to proton
  TCanvas* cRatioToProton = new TCanvas("cRatioToProton","Ratio to proton",800,600);
  cRatioToProton->SetTicks();
  TH1D* Axis2 = new TH1D("", "", 50,0,50);
  setOPT_hists(Axis2,MultiplicityLabel, RatioLabel,510,20,0.5);
  Axis2->GetYaxis()->SetRangeUser(5e-7,5.5e-6);
  Axis2->GetXaxis()->SetLabelSize(0.035);
  Axis2->GetYaxis()->SetLabelSize(0.035);
  setOPT_graph(gRatioToProton,MultiplicityLabel, RatioLabel,510,20,1);
  setOPT_graph(gRatioToProtonSyst,MultiplicityLabel, RatioLabel,510,20,1);
  gRatioToProtonSyst->SetFillStyle(0);

  TLegend* legRatioToProton = plotLegend("left_top",LegendTitle, 1, 0.3, 0.03,-0.03);
  legRatioToProton->SetFillStyle(0);
  legRatioToProton->SetBorderSize(0);
  legRatioToProton->SetTextSize(legsize);
  legRatioToProton->AddEntry(gRatioToProton,"p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV","lp");
  Axis2->DrawCopy();
  gRatioToProtonSyst->DrawClone("2same");
  gRatioToProton->DrawClone("PEsame");
  legRatioToProton->Draw();
  cRatioToProton->SaveAs("RatioToProton.pdf");

  TFile *outputFile = new TFile ("RatioToProton.root","RECREATE");
  outputFile -> cd();
  gRatioToProtonSyst->Write("RatioToProtonSyst");
  gRatioToProton->Write("RatioToProton");
  outputFile->Close();

  //<pT> vs multiplicity
  TCanvas* cMeanPt = new TCanvas("cMeanPt","<#it{p}_{T}> vs multiplicity",800,600);
  cMeanPt->SetTicks();
  TH1D* Axis3 = new TH1D("", "", 50,0,50);
  setOPT_hists(Axis3,MultiplicityLabel, "<#it{p}_{T}>",510,20,0.5);
  Axis3->GetYaxis()->SetRangeUser(0,5);
  Axis3->GetXaxis()->SetLabelSize(0.035);
  Axis3->GetYaxis()->SetLabelSize(0.035);
  gMeanPtSyst->SetFillStyle(0);

  TLegend* legMeanPt = plotLegend("left_top",LegendTitle, 1, 0.3, 0.03,-0.03);
  legMeanPt->SetFillStyle(0);
  legMeanPt->SetBorderSize(0);
  legMeanPt->SetTextSize(legsize);
  legMeanPt->AddEntry(gMeanPt,AnalysisName,"lp");
  Axis3->DrawCopy();
  gMeanPtSyst->DrawClone("2same");
  gMeanPt->DrawClone("PEsame");
  legMeanPt->Draw();
  cMeanPt->SaveAs("MeanPt.pdf");

  TFile *outputFile2 = new TFile ("MeanPt.root","Update");
  outputFile2 -> cd();
  gMeanPtSyst->Write("MeanPtSyst_Average"); // MeanPtSyst_AntiHelium MeanPtSyst_Helium MeanPtSyst_Average
  gMeanPt->Write("MeanPt_Average"); // MeanPt_AntiHelium MeanPt_Helium MeanPt_Average
  outputFile2->Close();
}
