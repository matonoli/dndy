#include <stdio.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "../mystyle.C"

using namespace std;

Double_t minpT = 0.4;
Double_t maxpT = 1.75;

TH1D* ChangeToPtOverA(TH1D* spectrum){
  Int_t nPtBins = spectrum->GetNbinsX();
  Double_t* Binning = (Double_t*) spectrum->GetXaxis()->GetXbins()->GetArray();
  for (int i=0; i < nPtBins+1; i++) {
    Binning[i] = Binning[i]/3.;
  }
  TString name = spectrum->GetName();
  name += "PtOverA";
  TH1D* HistoPtOverA = new TH1D(name.Data(),"", nPtBins, Binning);

  for (int i=0; i<nPtBins+1; i++) {
    HistoPtOverA->SetBinContent(i,spectrum->GetBinContent(i));
    HistoPtOverA->SetBinError(i,spectrum->GetBinError(i));
  }
  return HistoPtOverA;
}

TH1D* InterpolateProton(TH1D* Proton, TH1D* Helium){

  TH1D* Interpolated = (TH1D*) Helium->Clone("InterpolatedProton");
  int ProtonBinLow, ProtonBinUp;
  double pTLow, pTUp, pT, y, estat;
  for (int i = 0; i < Helium->GetNbinsX()+1; i++) {
    pTLow = Helium->GetBinLowEdge(i);
    pTUp = Helium->GetBinLowEdge(i+1);
    pT = Helium->GetBinCenter(i);
    ProtonBinLow = Proton->FindBin(pTLow);
    ProtonBinUp = Proton->FindBin(pTUp);
    y = Proton->GetBinContent(ProtonBinLow) + (Proton->GetBinContent(ProtonBinUp) - Proton->GetBinContent(ProtonBinLow))/(Proton->GetBinCenter(ProtonBinUp) - Proton->GetBinCenter(ProtonBinLow))*(pT - Proton->GetBinCenter(ProtonBinLow));
    estat = Proton->GetBinError(ProtonBinLow) + (Proton->GetBinError(ProtonBinUp) - Proton->GetBinError(ProtonBinLow))/(Proton->GetBinCenter(ProtonBinUp) - Proton->GetBinCenter(ProtonBinLow))*(pT - Proton->GetBinCenter(ProtonBinLow));

    Interpolated->SetBinContent(i,y);
    Interpolated->SetBinError(i,estat);
  }
  return Interpolated;
}

void CalculateB3(){
  
  TString file1 = "../Helium_LHC16q_FAST_woSDD_0_10.root",
          file2="../Helium_LHC16q_FAST_woSDD_Mult_Int.root",
          file3="../Helium_LHC16q_FAST_woSDD_10_20.root",
          file4="../Helium_LHC16q_FAST_woSDD_20_40.root",
          file5="../Helium_LHC16q_FAST_woSDD_40_100.root";

  bool DrawCoal = true;

  TString Objectname = "fHistoHe3PtSpectrum"; // 3He
  TString ObjectnameTriton = "fHistoPtSpectrum"; // triton

//  TString Objectname = "fHistoAntiHe3PtSpectrum"; // Anti-3He
//  TString ObjectnameTriton = "fHistoAntiPtSpectrum"; // Anti-triton

  Bool_t average= true; // triton and low pT He systematic uncertainties uncorrelated
  Double_t SystUncA, SystUncB;
  TString Objectname2, ObjectnameTriton2;
  if (average){
    Objectname2 = "fHistoAntiHe3PtSpectrum";
    ObjectnameTriton2 = "fHistoAntiPtSpectrum"; // Anti-triton
  }

  TString LegendTitle = "ALICE Work in progress"; // ALICE Work in progress ALICE Preliminary
  TString AnalysisName = "#kern[-0.19]{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}";
  TString RapidityInfo = "#kern[-0.3]{#minus1 #scale[1.2]{#leq} #it{y}_{cms} < 0}";
  TString NucleiType = "#kern[-0.3]{(^{3}He + ^{3}#bar{He})/2}"; // ^{3}#bar{He}; ^{3}He; (^{3}He + ^{3}#bar{He})/2
  Int_t nFiles = 5;

  Double_t textsize = 0.05;
  Double_t legsize = 0.045;

  //style area
//  gSystem->Load("../mystyle_C.so"); //root 5
  mystyle();

  TString pTLabel= "#it{p}_{T} (GeV/#it{c})";
  TString PtOverALabel= "#it{p}_{T} / A (GeV/#it{c})";
  TString yLable = "1/(2#pi#it{p}_{T}#it{N}_{evt}) d^{2}#it{N}/(d#it{y}d#it{p}_{T}) (GeV/#it{c})^{-2}";
  TString B3Lable = "#it{B}_{3} (GeV^{4}/#it{c}^{6})";
  TString RatioLabel = "^{3}H / ^{3}He"; // "^{3}#bar{H} / ^{3}#bar{He}"

  // load proton
  TFile ProtonFile ("Proton_inveriantyield_pPb_0-40.root");
  TDirectory *d = ProtonFile.GetDirectory("Table 7");

  TH1D* hProton_0_5 = (TH1D*) d->GetObjectChecked("Hist1D_y1","TH1F");
  TH1D* hProton_0_5_stat = (TH1D*) d->GetObjectChecked("Hist1D_y1_e1","TH1F");
  TH1D* hProton_0_5_syst = (TH1D*) d->GetObjectChecked("Hist1D_y1_e2","TH1F");
  hProton_0_5_stat->SetTitle("");
  hProton_0_5_syst->SetTitle("");
  hProton_0_5_stat->SetName("hProton_0_5_stat");
  hProton_0_5_syst->SetName("hProton_0_5_syst");

  Double_t Value, StatUnc, SystUnc;
  for (int i = 0; i < hProton_0_5->GetNbinsX(); i++) {
    Value = hProton_0_5->GetBinContent(i);
    StatUnc = hProton_0_5_stat->GetBinContent(i);
    SystUnc = hProton_0_5_syst->GetBinContent(i);

    hProton_0_5_stat->SetBinContent(i,Value);
    hProton_0_5_stat->SetBinError(i,StatUnc);

    hProton_0_5_syst->SetBinContent(i,Value);
    hProton_0_5_syst->SetBinError(i,SystUnc);
  }

  setOPT_hists(hProton_0_5_stat,pTLabel, yLable,510,20,1);
  setOPT_hists(hProton_0_5_syst,pTLabel, yLable,510,20,1);
  hProton_0_5_syst->SetFillStyle(0);

  TH1D* hProton_5_10 = (TH1D*) d->GetObjectChecked("Hist1D_y2","TH1F");
  TH1D* hProton_5_10_stat = (TH1D*) d->GetObjectChecked("Hist1D_y2_e1","TH1F");
  TH1D* hProton_5_10_syst = (TH1D*) d->GetObjectChecked("Hist1D_y2_e2","TH1F");
  hProton_5_10_stat->SetTitle("");
  hProton_5_10_syst->SetTitle("");
  hProton_5_10_stat->SetName("hProton_5_10_stat");
  hProton_5_10_syst->SetName("hProton_5_10_syst");

  for (int i = 0; i < hProton_5_10->GetNbinsX(); i++) {
    Value = hProton_5_10->GetBinContent(i);
    StatUnc = hProton_5_10_stat->GetBinContent(i);
    SystUnc = hProton_5_10_syst->GetBinContent(i);

    hProton_5_10_stat->SetBinContent(i,Value);
    hProton_5_10_stat->SetBinError(i,StatUnc);

    hProton_5_10_syst->SetBinContent(i,Value);
    hProton_5_10_syst->SetBinError(i,SystUnc);
  }

  setOPT_hists(hProton_5_10_stat,pTLabel, yLable,510,20,1);
  setOPT_hists(hProton_5_10_syst,pTLabel, yLable,510,20,1);
  hProton_5_10_syst->SetFillStyle(0);

  TH1D* hProton_0_10_stat = (TH1D*) hProton_0_5_stat->Clone("hProton_0_10_stat");
  hProton_0_10_stat->Add(hProton_5_10_stat);
  hProton_0_10_stat->Scale(0.5);

  TH1D* hProton_0_10_syst = (TH1D*) hProton_0_5_syst->Clone("hProton_0_10_syst");
  hProton_0_10_syst->Add(hProton_5_10_syst);
  hProton_0_10_syst->Scale(0.5);

  TH1D* hProton_10_20 = (TH1D*) d->GetObjectChecked("Hist1D_y3","TH1F");
  TH1D* hProton_10_20_stat = (TH1D*) d->GetObjectChecked("Hist1D_y3_e1","TH1F");
  TH1D* hProton_10_20_syst = (TH1D*) d->GetObjectChecked("Hist1D_y3_e2","TH1F");
  hProton_10_20_stat->SetTitle("");
  hProton_10_20_syst->SetTitle("");
  hProton_10_20_stat->SetName("hProton_10_20_stat");
  hProton_10_20_syst->SetName("hProton_10_20_syst");

  for (int i = 0; i < hProton_10_20->GetNbinsX(); i++) {
    Value = hProton_10_20->GetBinContent(i);
    StatUnc = hProton_10_20_stat->GetBinContent(i);
    SystUnc = hProton_10_20_syst->GetBinContent(i);

    hProton_10_20_stat->SetBinContent(i,Value);
    hProton_10_20_stat->SetBinError(i,StatUnc);

    hProton_10_20_syst->SetBinContent(i,Value);
    hProton_10_20_syst->SetBinError(i,SystUnc);
  }

  setOPT_hists(hProton_10_20_stat,pTLabel, yLable,510,20,1);
  setOPT_hists(hProton_10_20_syst,pTLabel, yLable,510,20,1);
  hProton_10_20_syst->SetFillStyle(0);

  TH1D* hProton_20_40 = (TH1D*) d->GetObjectChecked("Hist1D_y4","TH1F");
  TH1D* hProton_20_40_stat = (TH1D*) d->GetObjectChecked("Hist1D_y4_e1","TH1F");
  TH1D* hProton_20_40_syst = (TH1D*) d->GetObjectChecked("Hist1D_y4_e2","TH1F");
  hProton_20_40_stat->SetTitle("");
  hProton_20_40_syst->SetTitle("");
  hProton_20_40_stat->SetName("hProton_20_40_stat");
  hProton_20_40_syst->SetName("hProton_20_40_syst");

  for (int i = 0; i < hProton_20_40->GetNbinsX(); i++) {
    Value = hProton_20_40->GetBinContent(i);
    StatUnc = hProton_20_40_stat->GetBinContent(i);
    SystUnc = hProton_20_40_syst->GetBinContent(i);

    hProton_20_40_stat->SetBinContent(i,Value);
    hProton_20_40_stat->SetBinError(i,StatUnc);

    hProton_20_40_syst->SetBinContent(i,Value);
    hProton_20_40_syst->SetBinError(i,SystUnc);
  }

  setOPT_hists(hProton_20_40_stat,pTLabel, yLable,510,20,1);
  setOPT_hists(hProton_20_40_syst,pTLabel, yLable,510,20,1);
  hProton_20_40_syst->SetFillStyle(0);

  TFile ProtonFile2 ("Proton_inveriantyield_pPb_40-100.root");
  TDirectory *d2 = ProtonFile2.GetDirectory("Table 8");

  TH1D* hProton_40_60 = (TH1D*) d2->GetObjectChecked("Hist1D_y1","TH1F");
  TH1D* hProton_40_60_stat = (TH1D*) d2->GetObjectChecked("Hist1D_y1_e1","TH1F");
  TH1D* hProton_40_60_syst = (TH1D*) d2->GetObjectChecked("Hist1D_y1_e2","TH1F");
  hProton_40_60_stat->SetTitle("");
  hProton_40_60_syst->SetTitle("");
  hProton_40_60_stat->SetName("hProton_40_60_stat");
  hProton_40_60_syst->SetName("hProton_40_60_syst");

  for (int i = 0; i < hProton_40_60->GetNbinsX(); i++) {
    Value = hProton_40_60->GetBinContent(i);
    StatUnc = hProton_40_60_stat->GetBinContent(i);
    SystUnc = hProton_40_60_syst->GetBinContent(i);

    hProton_40_60_stat->SetBinContent(i,Value);
    hProton_40_60_stat->SetBinError(i,StatUnc);

    hProton_40_60_syst->SetBinContent(i,Value);
    hProton_40_60_syst->SetBinError(i,SystUnc);
  }

  setOPT_hists(hProton_40_60_stat,pTLabel, yLable,510,20,1);
  setOPT_hists(hProton_40_60_syst,pTLabel, yLable,510,20,1);
  hProton_40_60_syst->SetFillStyle(0);

  TH1D* hProton_60_80 = (TH1D*) d2->GetObjectChecked("Hist1D_y2","TH1F");
  TH1D* hProton_60_80_stat = (TH1D*) d2->GetObjectChecked("Hist1D_y2_e1","TH1F");
  TH1D* hProton_60_80_syst = (TH1D*) d2->GetObjectChecked("Hist1D_y2_e2","TH1F");
  hProton_60_80_stat->SetTitle("");
  hProton_60_80_syst->SetTitle("");
  hProton_60_80_stat->SetName("hProton_60_80_stat");
  hProton_60_80_syst->SetName("hProton_60_80_syst");

  for (int i = 0; i < hProton_60_80->GetNbinsX(); i++) {
    Value = hProton_60_80->GetBinContent(i);
    StatUnc = hProton_60_80_stat->GetBinContent(i);
    SystUnc = hProton_60_80_syst->GetBinContent(i);

    hProton_60_80_stat->SetBinContent(i,Value);
    hProton_60_80_stat->SetBinError(i,StatUnc);

    hProton_60_80_syst->SetBinContent(i,Value);
    hProton_60_80_syst->SetBinError(i,SystUnc);
  }

  setOPT_hists(hProton_60_80_stat,pTLabel, yLable,510,20,1);
  setOPT_hists(hProton_60_80_syst,pTLabel, yLable,510,20,1);
  hProton_60_80_syst->SetFillStyle(0);

  TH1D* hProton_80_100 = (TH1D*) d2->GetObjectChecked("Hist1D_y3","TH1F");
  TH1D* hProton_80_100_stat = (TH1D*) d2->GetObjectChecked("Hist1D_y3_e1","TH1F");
  TH1D* hProton_80_100_syst = (TH1D*) d2->GetObjectChecked("Hist1D_y3_e2","TH1F");
  hProton_80_100_stat->SetTitle("");
  hProton_80_100_syst->SetTitle("");
  hProton_80_100_stat->SetName("hProton_80_100_stat");
  hProton_80_100_syst->SetName("hProton_80_100_syst");

  for (int i = 0; i < hProton_80_100->GetNbinsX(); i++) {
    Value = hProton_80_100->GetBinContent(i);
    StatUnc = hProton_80_100_stat->GetBinContent(i);
    SystUnc = hProton_80_100_syst->GetBinContent(i);

    hProton_80_100_stat->SetBinContent(i,Value);
    hProton_80_100_stat->SetBinError(i,StatUnc);

    hProton_80_100_syst->SetBinContent(i,Value);
    hProton_80_100_syst->SetBinError(i,SystUnc);
  }

  setOPT_hists(hProton_80_100_stat,pTLabel, yLable,510,20,1);
  setOPT_hists(hProton_80_100_syst,pTLabel, yLable,510,20,1);
  hProton_80_100_syst->SetFillStyle(0);

  TH1D* hProton_40_100_stat = (TH1D*) hProton_40_60_stat->Clone("hProton_40_100_stat");
  hProton_40_100_stat->Add(hProton_60_80_stat);
  hProton_40_100_stat->Add(hProton_80_100_stat);
  hProton_40_100_stat->Scale(1./3.);

  TH1D* hProton_40_100_syst = (TH1D*) hProton_40_60_syst->Clone("hProton_40_100_syst");
  hProton_40_100_syst->Add(hProton_60_80_syst);
  hProton_40_100_syst->Add(hProton_80_100_syst);
  hProton_40_100_syst->Scale(1./3.);

  TCanvas* ctest = new TCanvas("ctest","Interpolation",800,600);

  hProton_0_10_stat->Scale(0.5);
  hProton_10_20_stat->Scale(0.5);
  hProton_20_40_stat->Scale(0.5);
  hProton_40_100_stat->Scale(0.5);

  hProton_0_10_syst->Scale(0.5);
  hProton_10_20_syst->Scale(0.5);
  hProton_20_40_syst->Scale(0.5);
  hProton_40_100_syst->Scale(0.5);

  hProton_0_10_syst->DrawCopy("E2");
  hProton_10_20_syst->DrawCopy("E2same");
  hProton_20_40_syst->DrawCopy("E2same");
  hProton_40_100_syst->DrawCopy("E2same");

  hProton_0_10_stat->DrawCopy("Esame");
  hProton_10_20_stat->DrawCopy("Esame");
  hProton_20_40_stat->DrawCopy("Esame");
  hProton_40_100_stat->DrawCopy("Esame");

  // weighted average for multplicity intergrated
  hProton_20_40_stat->Scale(2.);
  hProton_40_100_stat->Scale(6.);

  hProton_20_40_syst->Scale(2.);
  hProton_40_100_syst->Scale(6.);

  TH1D* hProton_MultInt_stat = (TH1D*) hProton_0_10_stat->Clone("hProton_MultInt_stat");
  hProton_MultInt_stat->Add(hProton_10_20_stat);
  hProton_MultInt_stat->Add(hProton_20_40_stat);
  hProton_MultInt_stat->Add(hProton_40_100_stat);
  hProton_MultInt_stat->Scale(1./10.);

  TH1D* hProton_MultInt_syst = (TH1D*) hProton_0_10_syst->Clone("hProton_MultInt_syst");
  hProton_MultInt_syst->Add(hProton_10_20_syst);
  hProton_MultInt_syst->Add(hProton_20_40_syst);
  hProton_MultInt_syst->Add(hProton_40_100_syst);
  hProton_MultInt_syst->Scale(1./10.);

  // weighted average for multplicity intergrated (undo)
  hProton_20_40_stat->Scale(1./2.);
  hProton_40_100_stat->Scale(1./6.);

  hProton_20_40_syst->Scale(1./2.);
  hProton_40_100_syst->Scale(1./6.);
  
  // load Helium
  TFile Numeratorfile (file1.Data());
  TH1D* hNumerator = (TH1D*) Numeratorfile.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hNumeratorB = (TH1D*) Numeratorfile.GetObjectChecked(Objectname2,"TH1D");
    hNumerator->Add(hNumeratorB);
    hNumerator->Scale(0.5);
  }
  setOPT_hists(hNumerator,pTLabel, yLable,510,20,1,kRed);
  hNumerator->SetTitle("");
  hNumerator->SetName("Numerator");
  TH1D* hNumeratorSyst = (TH1D*) Numeratorfile.GetObjectChecked(Objectname+"Syst","TH1D");
  if (average){
    TH1D* hNumeratorSystA = (TH1D*) Numeratorfile.GetObjectChecked(Objectname+"Syst","TH1D");
    TH1D* hNumeratorSystB = (TH1D*) Numeratorfile.GetObjectChecked(Objectname2+"Syst","TH1D");
    hNumeratorSyst->Add(hNumeratorSystB);
    hNumeratorSyst->Scale(0.5);
    SystUncA = 0;
    SystUncB = 0;
    for (int i =0; i <= hNumeratorSyst->GetNbinsX(); i++) {
      if (hNumeratorSyst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
      SystUncA = hNumeratorSystA->GetBinError(i);
      SystUncB = hNumeratorSystB->GetBinError(i);
      hNumeratorSyst->SetBinError(i,0.5*(SystUncA+SystUncB));
    }
  }
  setOPT_hists(hNumeratorSyst,pTLabel, yLable,510,20,1,kRed);
  hNumeratorSyst->SetFillStyle(0);

  // convert to invariant yield
  for (int i=0; i <= hNumerator->GetNbinsX(); i++) {
    hNumerator->SetBinContent(i, hNumerator->GetBinContent(i)/hNumerator->GetBinCenter(i));
    hNumerator->SetBinError(i, hNumerator->GetBinError(i)/hNumerator->GetBinCenter(i));

    hNumeratorSyst->SetBinContent(i, hNumeratorSyst->GetBinContent(i)/hNumeratorSyst->GetBinCenter(i));
    hNumeratorSyst->SetBinError(i, hNumeratorSyst->GetBinError(i)/hNumeratorSyst->GetBinCenter(i));

  }
  hNumerator->Scale(1./(2.*TMath::Pi()));
  hNumeratorSyst->Scale(1./(2.*TMath::Pi()));

  TH1D* hNumeratorRig = (TH1D*) ChangeToPtOverA(hNumerator);
  setOPT_hists(hNumeratorRig,PtOverALabel, B3Lable,510,21,1,kRed);
  TH1D* hNumeratorSystRig = (TH1D*) ChangeToPtOverA(hNumeratorSyst);
  setOPT_hists(hNumeratorSystRig,PtOverALabel, B3Lable,510,21,1,kRed);

  TFile Denominatorfile (file2.Data());
  TH1D* hDenominator = (TH1D*) Denominatorfile.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hDenominatorB = (TH1D*) Denominatorfile.GetObjectChecked(Objectname2,"TH1D");
    hDenominator->Add(hDenominatorB);
    hDenominator->Scale(0.5);
  }
  setOPT_hists(hDenominator,pTLabel, yLable,510,20,1);
  hDenominator->SetTitle("");
  hDenominator->SetName("Denominator");
  TH1D* hDenominatorSyst = (TH1D*) Denominatorfile.GetObjectChecked(Objectname+"Syst","TH1D");
  if (average){
    TH1D* hDenominatorSystA = (TH1D*) Denominatorfile.GetObjectChecked(Objectname+"Syst","TH1D");
    TH1D* hDenominatorSystB = (TH1D*) Denominatorfile.GetObjectChecked(Objectname2+"Syst","TH1D");
    hDenominatorSyst->Add(hDenominatorSystB);
    hDenominatorSyst->Scale(0.5);
    SystUncA = 0;
    SystUncB = 0;
    for (int i =0; i <= hDenominatorSyst->GetNbinsX(); i++) {
      if (hDenominatorSyst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
      SystUncA = hDenominatorSystA->GetBinError(i);
      SystUncB = hDenominatorSystB->GetBinError(i);
      hDenominatorSyst->SetBinError(i,0.5*(SystUncA+SystUncB));
    }
  }
  setOPT_hists(hDenominatorSyst,pTLabel, yLable,510,20,1);
  hDenominatorSyst->SetFillStyle(0);

  // convert to invariant yield
  for (int i=0; i <= hDenominator->GetNbinsX(); i++) {
    hDenominator->SetBinContent(i, hDenominator->GetBinContent(i)/hDenominator->GetBinCenter(i));
    hDenominator->SetBinError(i, hDenominator->GetBinError(i)/hDenominator->GetBinCenter(i));

    hDenominatorSyst->SetBinContent(i, hDenominatorSyst->GetBinContent(i)/hDenominatorSyst->GetBinCenter(i));
    hDenominatorSyst->SetBinError(i, hDenominatorSyst->GetBinError(i)/hDenominatorSyst->GetBinCenter(i));

  }
  hDenominator->Scale(1./(2.*TMath::Pi()));
  hDenominatorSyst->Scale(1./(2.*TMath::Pi()));
  TH1D* hDenominatorRig = (TH1D*) ChangeToPtOverA(hDenominator);
  setOPT_hists(hDenominatorRig,PtOverALabel, B3Lable,510,20,1);
  TH1D* hDenominatorSystRig = (TH1D*) ChangeToPtOverA(hDenominatorSyst);
  setOPT_hists(hDenominatorSystRig,PtOverALabel, B3Lable,510,20,1);

  TFile Numerator2file (file3.Data());
  TH1D* hNumerator2 = (TH1D*) Numerator2file.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hNumerator2B = (TH1D*) Numerator2file.GetObjectChecked(Objectname2,"TH1D");
    hNumerator2->Add(hNumerator2B);
    hNumerator2->Scale(0.5);
  }
  setOPT_hists(hNumerator2,pTLabel, yLable,510,20,1,kBlue);
  hNumerator2->SetTitle("");
  hNumerator2->SetName("Numerator");
  
  TH1D* hNumerator2Syst = (TH1D*) Numerator2file.GetObjectChecked(Objectname+"Syst","TH1D");
  if (average){
    TH1D* hNumerator2SystA = (TH1D*) Numerator2file.GetObjectChecked(Objectname+"Syst","TH1D");
    TH1D* hNumerator2SystB = (TH1D*) Numerator2file.GetObjectChecked(Objectname2+"Syst","TH1D");
    hNumerator2Syst->Add(hNumerator2SystB);
    hNumerator2Syst->Scale(0.5);
    SystUncA = 0;
    SystUncB = 0;
    for (int i =0; i <= hNumerator2Syst->GetNbinsX(); i++) {
      if (hNumerator2Syst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
      SystUncA = hNumerator2SystA->GetBinError(i);
      SystUncB = hNumerator2SystB->GetBinError(i);
      hNumerator2Syst->SetBinError(i,0.5*(SystUncA+SystUncB));
    }
  }
  setOPT_hists(hNumerator2Syst,pTLabel, yLable,510,20,1,kBlue);
  hNumerator2Syst->SetFillStyle(0);
  
  // convert to invariant yield
  for (int i=0; i <= hNumerator2->GetNbinsX(); i++) {
    hNumerator2->SetBinContent(i, hNumerator2->GetBinContent(i)/hNumerator2->GetBinCenter(i));
    hNumerator2->SetBinError(i, hNumerator2->GetBinError(i)/hNumerator2->GetBinCenter(i));
    
    hNumerator2Syst->SetBinContent(i, hNumerator2Syst->GetBinContent(i)/hNumerator2Syst->GetBinCenter(i));
    hNumerator2Syst->SetBinError(i, hNumerator2Syst->GetBinError(i)/hNumerator2Syst->GetBinCenter(i));
    
  }
  hNumerator2->Scale(1./(2.*TMath::Pi()));
  hNumerator2Syst->Scale(1./(2.*TMath::Pi()));
  
  TH1D* hNumerator2Rig = (TH1D*) ChangeToPtOverA(hNumerator2);
  setOPT_hists(hNumerator2Rig,PtOverALabel, B3Lable,510,33,1.5,kBlue);
  TH1D* hNumerator2SystRig = (TH1D*) ChangeToPtOverA(hNumerator2Syst);
  setOPT_hists(hNumerator2SystRig,PtOverALabel, B3Lable,510,33,1.5,kBlue);

  TFile Numerator3file (file4.Data());
  TH1D* hNumerator3 = (TH1D*) Numerator3file.GetObjectChecked(Objectname,"TH1D");
  if (average){
    TH1D* hNumerator3B = (TH1D*) Numerator3file.GetObjectChecked(Objectname2,"TH1D");
    hNumerator3->Add(hNumerator3B);
    hNumerator3->Scale(0.5);
  }
  setOPT_hists(hNumerator3,pTLabel, yLable,510,20,1,kOrange-3);
  hNumerator3->SetTitle("");
  hNumerator3->SetName("Numerator");
  TH1D* hNumerator3Syst = (TH1D*) Numerator3file.GetObjectChecked(Objectname+"Syst","TH1D");
  if (average){
    TH1D* hNumerator3SystA = (TH1D*) Numerator3file.GetObjectChecked(Objectname+"Syst","TH1D");
    TH1D* hNumerator3SystB = (TH1D*) Numerator3file.GetObjectChecked(Objectname2+"Syst","TH1D");
    hNumerator3Syst->Add(hNumerator3SystB);
    hNumerator3Syst->Scale(0.5);
    SystUncA = 0;
    SystUncB = 0;
    for (int i =0; i <= hNumerator3Syst->GetNbinsX(); i++) {
      if (hNumerator3Syst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
      SystUncA = hNumerator3SystA->GetBinError(i);
      SystUncB = hNumerator3SystB->GetBinError(i);
      hNumerator3Syst->SetBinError(i,0.5*(SystUncA+SystUncB));
    }
  }
  setOPT_hists(hNumerator3Syst,pTLabel, yLable,510,20,1,kOrange-3);
  hNumerator3Syst->SetFillStyle(0);
  // convert to invariant yield
  for (int i=0; i <= hNumerator3->GetNbinsX(); i++) {
    hNumerator3->SetBinContent(i, hNumerator3->GetBinContent(i)/hNumerator3->GetBinCenter(i));
    hNumerator3->SetBinError(i, hNumerator3->GetBinError(i)/hNumerator3->GetBinCenter(i));
    
    hNumerator3Syst->SetBinContent(i, hNumerator3Syst->GetBinContent(i)/hNumerator3Syst->GetBinCenter(i));
    hNumerator3Syst->SetBinError(i, hNumerator3Syst->GetBinError(i)/hNumerator3Syst->GetBinCenter(i));
    
  }
  hNumerator3->Scale(1./(2.*TMath::Pi()));
  hNumerator3Syst->Scale(1./(2.*TMath::Pi()));
  TH1D* hNumerator3Rig = (TH1D*) ChangeToPtOverA(hNumerator3);
  setOPT_hists(hNumerator3Rig,PtOverALabel, B3Lable,510,34,1.25,kOrange-3);
  TH1D* hNumerator3SystRig = (TH1D*) ChangeToPtOverA(hNumerator3Syst);
  setOPT_hists(hNumerator3SystRig,PtOverALabel, B3Lable,510,34,1.25,kOrange-3);

    TFile Numerator4file (file5.Data());
    TH1D* hNumerator4 = (TH1D*) Numerator4file.GetObjectChecked(Objectname,"TH1D");
    if (average){
      TH1D* hNumerator4B = (TH1D*) Numerator4file.GetObjectChecked(Objectname2,"TH1D");
      hNumerator4->Add(hNumerator4B);
      hNumerator4->Scale(0.5);
    }
    setOPT_hists(hNumerator4,pTLabel, yLable,510,20,1,kMagenta-3);
    hNumerator4->SetTitle("");
    hNumerator4->SetName("Numerator");
    TH1D* hNumerator4Syst = (TH1D*) Numerator4file.GetObjectChecked(Objectname+"Syst","TH1D");
    if (average){
      TH1D* hNumerator4SystA = (TH1D*) Numerator4file.GetObjectChecked(Objectname+"Syst","TH1D");
      TH1D* hNumerator4SystB = (TH1D*) Numerator4file.GetObjectChecked(Objectname2+"Syst","TH1D");
      hNumerator4Syst->Add(hNumerator4SystB);
      hNumerator4Syst->Scale(0.5);
      SystUncA = 0;
      SystUncB = 0;
      for (int i =0; i <= hNumerator4Syst->GetNbinsX(); i++) {
        if (hNumerator4Syst->GetBinCenter(i) < 3) continue; // systematic uncertianty dominated by primary fraction below pT = 3 GeV/c => not correlated
        SystUncA = hNumerator4SystA->GetBinError(i);
        SystUncB = hNumerator4SystB->GetBinError(i);
        hNumerator4Syst->SetBinError(i,0.5*(SystUncA+SystUncB));
      }
    }
    setOPT_hists(hNumerator4Syst,pTLabel, yLable,510,20,1,kMagenta-3);
    hNumerator4Syst->SetFillStyle(0);
    // convert to invariant yield
    for (int i=0; i <= hNumerator4->GetNbinsX(); i++) {
      hNumerator4->SetBinContent(i, hNumerator4->GetBinContent(i)/hNumerator4->GetBinCenter(i));
      hNumerator4->SetBinError(i, hNumerator4->GetBinError(i)/hNumerator4->GetBinCenter(i));

      hNumerator4Syst->SetBinContent(i, hNumerator4Syst->GetBinContent(i)/hNumerator4Syst->GetBinCenter(i));
      hNumerator4Syst->SetBinError(i, hNumerator4Syst->GetBinError(i)/hNumerator4Syst->GetBinCenter(i));

    }
    hNumerator4->Scale(1./(2.*TMath::Pi()));
    hNumerator4Syst->Scale(1./(2.*TMath::Pi()));

    TH1D* hNumerator4Rig = (TH1D*) ChangeToPtOverA(hNumerator4);
    setOPT_hists(hNumerator4Rig,PtOverALabel, B3Lable,510,43,2,kMagenta-3);
    TH1D* hNumerator4SystRig = (TH1D*) ChangeToPtOverA(hNumerator4Syst);
    setOPT_hists(hNumerator4SystRig,PtOverALabel, B3Lable,510,43,2,kMagenta-3);
  
  // triton
  TFile Denominatorfile2 ("../Triton/Triton_LHC16q_FAST_woSDD.root");
  TH1D* hDenominator2 = (TH1D*) Denominatorfile2.GetObjectChecked(ObjectnameTriton,"TH1D");
  if (average){
    TH1D* hDenominator2B = (TH1D*) Denominatorfile2.GetObjectChecked(ObjectnameTriton2,"TH1D");
    hDenominator2->Add(hDenominator2B);
    hDenominator2->Scale(0.5);
  }
  setOPT_hists(hDenominator2,pTLabel, yLable,510,21,1, kRed);
  hDenominator2->SetTitle("");
  hDenominator2->SetName("Denominator2");
  TH1D* hDenominatorSyst2 = (TH1D*) Denominatorfile2.GetObjectChecked(ObjectnameTriton+"Syst","TH1D");
  if (average){
    TH1D* hDenominatorSyst2A = (TH1D*) Denominatorfile2.GetObjectChecked(ObjectnameTriton+"Syst","TH1D");
    TH1D* hDenominatorSyst2B = (TH1D*) Denominatorfile2.GetObjectChecked(ObjectnameTriton2+"Syst","TH1D");
    hDenominatorSyst2->Add(hDenominatorSyst2B);
    hDenominatorSyst2->Scale(0.5);
  }
  setOPT_hists(hDenominatorSyst2,pTLabel, yLable,510,21,1);
  hDenominatorSyst2->SetFillStyle(0);

  // convert to invariant yield
  for (int i=0; i <= hDenominator->GetNbinsX(); i++) {
    hDenominator2->SetBinContent(i, hDenominator2->GetBinContent(i)/hDenominator2->GetBinCenter(i));
    hDenominator2->SetBinError(i, hDenominator2->GetBinError(i)/hDenominator2->GetBinCenter(i));

    hDenominatorSyst2->SetBinContent(i, hDenominatorSyst2->GetBinContent(i)/hDenominatorSyst2->GetBinCenter(i));
    hDenominatorSyst2->SetBinError(i, hDenominatorSyst2->GetBinError(i)/hDenominatorSyst2->GetBinCenter(i));

  }
  hDenominator2->Scale(1./(2.*TMath::Pi()));
  hDenominatorSyst2->Scale(1./(2.*TMath::Pi()));
  TH1D* hDenominatorRig2 = (TH1D*) ChangeToPtOverA(hDenominator2);
  setOPT_hists(hDenominatorRig2,PtOverALabel, B3Lable,510,21,1,kRed);
  TH1D* hDenominatorSystRig2 = (TH1D*) ChangeToPtOverA(hDenominatorSyst2);
  setOPT_hists(hDenominatorSystRig2,PtOverALabel, B3Lable,510,21,1,kRed);

  hNumeratorRig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);
  hNumeratorSystRig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);

  hDenominatorRig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);
  hDenominatorSystRig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);

  hNumerator2Rig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);
  hNumerator2SystRig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);

  hNumerator3Rig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);
  hNumerator3SystRig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);

  hNumerator4Rig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);
  hNumerator4SystRig->GetXaxis()->SetRangeUser(1.5/3., 5./3.);

  hDenominatorRig2->GetXaxis()->SetRangeUser(1.5/3., 5./3.);
  hDenominatorSystRig2->GetXaxis()->SetRangeUser(1.5/3., 5./3.);


  // Triton-over-3Helium ratio
  const Int_t nPtBins=4;
  Double_t Binning[]={1./3., 1.5/3., 2./3.,2.5/3.,3.0/3.};

  TH1D* hHeliumRebinned = (TH1D*) hDenominatorRig->Rebin(nPtBins,"hHelium",Binning);
  TH1D* hHeliumSystRebinned = (TH1D*) hDenominatorSystRig->Rebin(nPtBins,"hHeliumSyst",Binning);
  TH1D* hHeliumTotalRebinned = (TH1D*)hHeliumSystRebinned->Clone("hHeliumTotalRebinned");
  Double_t TotalUncertainty = 0;

  for (Int_t iBin=0; iBin<hHeliumTotalRebinned->GetNbinsX(); iBin++) {
    TotalUncertainty = sqrt(hDenominatorRig->GetBinError(iBin)*hDenominatorRig->GetBinError(iBin)+ hDenominatorSystRig->GetBinError(iBin)*hDenominatorSystRig->GetBinError(iBin));
    hHeliumTotalRebinned->SetBinError(iBin,TotalUncertainty);
  }

  TH1D* hRatio = (TH1D*) hDenominatorRig2->Clone("hRatio");
  TH1D* hRatioSyst = (TH1D*) hDenominatorSystRig2->Clone("hRatioSyst");
  TH1D* hRatioTotal = (TH1D*) hDenominatorSystRig2->Clone("hRatioTotal");
  for (Int_t iBin=0; iBin<hRatioTotal->GetNbinsX(); iBin++) {
    TotalUncertainty = sqrt(hDenominatorSystRig2->GetBinError(iBin)*hDenominatorSystRig2->GetBinError(iBin)+ hDenominatorRig2->GetBinError(iBin)*hDenominatorRig2->GetBinError(iBin));
    hRatioTotal->SetBinError(iBin,TotalUncertainty);
  }
  setOPT_hists(hRatio,pTLabel, RatioLabel,510,20,1);
  setOPT_hists(hRatioSyst,pTLabel, RatioLabel,510,20,1);
  setOPT_hists(hRatioTotal,pTLabel, RatioLabel,510,21,1);

  hRatio->Divide(hHeliumRebinned);
  hRatioSyst->Divide(hHeliumSystRebinned);
  hRatioTotal->Divide(hHeliumTotalRebinned);

  hRatio->SetFillStyle(0);
  hRatioSyst->SetFillStyle(0);
  
  // Adjust Protons
  //stat
  TH1D* hProton_0_10_stat_Rebinned = (TH1D*) InterpolateProton(hProton_0_10_stat, hNumeratorRig);
  hProton_0_10_stat_Rebinned->DrawCopy("same EX0");
  TH1D* hProton_0_10_stat_3 = (TH1D*) hProton_0_10_stat_Rebinned->Clone("hProton_0_10_stat_3");
  hProton_0_10_stat_3->Multiply(hProton_0_10_stat_Rebinned);
  hProton_0_10_stat_3->Multiply(hProton_0_10_stat_Rebinned);

  TH1D* hProton_10_20_stat_Rebinned = (TH1D*) InterpolateProton(hProton_10_20_stat, hNumerator2Rig);
  hProton_10_20_stat_Rebinned->DrawCopy("same EX0");
  TH1D* hProton_10_20_stat_3 = (TH1D*) hProton_10_20_stat_Rebinned->Clone("hProton_10_20_stat_3");
  hProton_10_20_stat_3->Multiply(hProton_10_20_stat_Rebinned);
  hProton_10_20_stat_3->Multiply(hProton_10_20_stat_Rebinned);

  TH1D* hProton_20_40_stat_Rebinned = (TH1D*) InterpolateProton(hProton_20_40_stat, hNumerator3Rig);
  hProton_20_40_stat_Rebinned->DrawCopy("same EX0");
  TH1D* hProton_20_40_stat_3 = (TH1D*) hProton_20_40_stat_Rebinned->Clone("hProton_20_40_stat_3");
  hProton_20_40_stat_3->Multiply(hProton_20_40_stat_Rebinned);
  hProton_20_40_stat_3->Multiply(hProton_20_40_stat_Rebinned);

  TH1D* hProton_40_100_stat_Rebinned = (TH1D*) InterpolateProton(hProton_40_100_stat, hNumerator4Rig);
  hProton_40_100_stat_Rebinned->DrawCopy("same EX0");
  TH1D* hProton_40_100_stat_3 = (TH1D*) hProton_40_100_stat_Rebinned->Clone("hProton_40_100_stat_3");
  hProton_40_100_stat_3->Multiply(hProton_40_100_stat_Rebinned);
  hProton_40_100_stat_3->Multiply(hProton_40_100_stat_Rebinned);

  //syst
  TH1D* hProton_0_10_syst_Rebinned = (TH1D*) InterpolateProton(hProton_0_10_syst, hNumeratorRig);
  hProton_0_10_syst_Rebinned->SetFillStyle(0);
  hProton_0_10_syst_Rebinned->DrawCopy("E2same");
  TH1D* hProton_0_10_syst_3 = (TH1D*) hProton_0_10_syst_Rebinned->Clone("hProton_0_10_syst_3");
  hProton_0_10_syst_3->Multiply(hProton_0_10_syst_Rebinned);
  hProton_0_10_syst_3->Multiply(hProton_0_10_syst_Rebinned);
  hProton_0_10_syst_3->SetFillStyle(0);

  TH1D* hProton_10_20_syst_Rebinned = (TH1D*) InterpolateProton(hProton_10_20_syst, hNumerator2Rig);
  hProton_10_20_syst_Rebinned->SetFillStyle(0);
  hProton_10_20_syst_Rebinned->DrawCopy("E2same");
  TH1D* hProton_10_20_syst_3 = (TH1D*) hProton_10_20_syst_Rebinned->Clone("hProton_10_20_syst_3");
  hProton_10_20_syst_3->Multiply(hProton_10_20_syst_Rebinned);
  hProton_10_20_syst_3->Multiply(hProton_10_20_syst_Rebinned);
  hProton_10_20_syst_3->SetFillStyle(0);

  TH1D* hProton_20_40_syst_Rebinned = (TH1D*) InterpolateProton(hProton_20_40_syst, hNumerator3Rig);
  hProton_20_40_syst_Rebinned->SetFillStyle(0);
  hProton_20_40_syst_Rebinned->DrawCopy("E2same");
  TH1D* hProton_20_40_syst_3 = (TH1D*) hProton_20_40_syst_Rebinned->Clone("hProton_20_40_syst_3");
  hProton_20_40_syst_3->Multiply(hProton_20_40_syst_Rebinned);
  hProton_20_40_syst_3->Multiply(hProton_20_40_syst_Rebinned);
  hProton_20_40_syst_3->SetFillStyle(0);

  TH1D* hProton_40_100_syst_Rebinned = (TH1D*) InterpolateProton(hProton_40_100_syst, hNumerator4Rig);
  hProton_40_100_syst_Rebinned->SetFillStyle(0);
  hProton_40_100_syst_Rebinned->DrawCopy("E2same");
  TH1D* hProton_40_100_syst_3 = (TH1D*) hProton_40_100_syst_Rebinned->Clone("hProton_40_100_syst_3");
  hProton_40_100_syst_3->Multiply(hProton_40_100_syst_Rebinned);
  hProton_40_100_syst_3->Multiply(hProton_40_100_syst_Rebinned);
  hProton_40_100_syst_3->SetFillStyle(0);

  // Multiplicity integrated
  TH1D* hProton_MultInt_syst_Rebinned = (TH1D*) InterpolateProton(hProton_MultInt_syst, hDenominatorRig);
  hProton_MultInt_syst_Rebinned->SetFillStyle(0);
  hProton_MultInt_syst_Rebinned->DrawCopy("E2same");
  TH1D* hProton_MultInt_syst_3 = (TH1D*) hProton_MultInt_syst_Rebinned->Clone("hProton_MultInt_syst_3");
  hProton_MultInt_syst_3->Multiply(hProton_MultInt_syst_Rebinned);
  hProton_MultInt_syst_3->Multiply(hProton_MultInt_syst_Rebinned);
  hProton_MultInt_syst_3->SetFillStyle(0);

  TH1D* hProton_MultInt_stat_Rebinned = (TH1D*) InterpolateProton(hProton_MultInt_stat, hDenominatorRig);
  hProton_MultInt_stat_Rebinned->SetFillStyle(0);
  hProton_MultInt_stat_Rebinned->DrawCopy("E2same");
  TH1D* hProton_MultInt_stat_3 = (TH1D*) hProton_MultInt_stat_Rebinned->Clone("hProton_MultInt_stat_3");
  hProton_MultInt_stat_3->Multiply(hProton_MultInt_stat_Rebinned);
  hProton_MultInt_stat_3->Multiply(hProton_MultInt_stat_Rebinned);
  hProton_MultInt_stat_3->SetFillStyle(0);

  hNumeratorRig->Divide(hProton_0_10_stat_3);
  hNumerator2Rig->Divide(hProton_10_20_stat_3);
  hNumerator3Rig->Divide(hProton_20_40_stat_3);
  hNumerator4Rig->Divide(hProton_40_100_stat_3);

  // syst unc. propagated with quadratic sum but correlated should be summed linear => sqrt(3) applied should be 3 => scale with sqrt(3)
  for (int i =0; i <= hProton_0_10_syst_3->GetNbinsX(); i++) {
    hProton_0_10_syst_3->SetBinError(i,hProton_0_10_syst_3->GetBinError(i)*TMath::Sqrt(3.));
  }

  for (int i =0; i <= hProton_10_20_syst_3->GetNbinsX(); i++) {
    hProton_10_20_syst_3->SetBinError(i,hProton_10_20_syst_3->GetBinError(i)*TMath::Sqrt(3.));
  }

  for (int i =0; i <= hProton_20_40_syst_3->GetNbinsX(); i++) {
    hProton_20_40_syst_3->SetBinError(i,hProton_20_40_syst_3->GetBinError(i)*TMath::Sqrt(3.));
  }

  for (int i =0; i <= hProton_40_100_syst_3->GetNbinsX(); i++) {
    hProton_40_100_syst_3->SetBinError(i,hProton_40_100_syst_3->GetBinError(i)*TMath::Sqrt(3.));
  }

  for (int i =0; i <= hProton_MultInt_syst_3->GetNbinsX(); i++) {
    hProton_MultInt_syst_3->SetBinError(i,hProton_MultInt_syst_3->GetBinError(i)*TMath::Sqrt(3.));
  }


  hNumeratorSystRig->Divide(hProton_0_10_syst_3);
  hNumerator2SystRig->Divide(hProton_10_20_syst_3);
  hNumerator3SystRig->Divide(hProton_20_40_syst_3);
  hNumerator4SystRig->Divide(hProton_40_100_syst_3);

  hDenominatorRig->Divide(hProton_MultInt_stat_3);
  hDenominatorSystRig->Divide(hProton_MultInt_syst_3);

  // Multiplicity integrated (triton)
  TH1D* hProton_MultInt_syst_Rebinned2 = (TH1D*) InterpolateProton(hProton_MultInt_syst, hDenominatorRig2);
  hProton_MultInt_syst_Rebinned2->SetFillStyle(0);
  hProton_MultInt_syst_Rebinned2->DrawCopy("E2same");
  TH1D* hProton_MultInt_syst2_3 = (TH1D*) hProton_MultInt_syst_Rebinned2->Clone("hProton_MultInt_syst2_3");
  hProton_MultInt_syst2_3->Multiply(hProton_MultInt_syst_Rebinned2);
  hProton_MultInt_syst2_3->Multiply(hProton_MultInt_syst_Rebinned2);
  hProton_MultInt_syst2_3->SetFillStyle(0);

  TH1D* hProton_MultInt_stat_Rebinned2 = (TH1D*) InterpolateProton(hProton_MultInt_stat, hDenominatorRig2);
  hProton_MultInt_stat_Rebinned2->SetFillStyle(0);
  hProton_MultInt_stat_Rebinned2->DrawCopy("E2same");
  TH1D* hProton_MultInt_stat2_3 = (TH1D*) hProton_MultInt_stat_Rebinned2->Clone("hProton_MultInt_stat2_3");
  hProton_MultInt_stat2_3->Multiply(hProton_MultInt_stat_Rebinned2);
  hProton_MultInt_stat2_3->Multiply(hProton_MultInt_stat_Rebinned2);
  hProton_MultInt_stat2_3->SetFillStyle(0);

  // syst unc. propagated with quadratic sum but correlated should be summed linear => sqrt(3) applied should be 3 => scale with sqrt(3)
  for (int i =0; i <= hProton_MultInt_syst2_3->GetNbinsX(); i++) {
    hProton_MultInt_syst2_3->SetBinError(i,hProton_MultInt_syst2_3->GetBinError(i)*TMath::Sqrt(3.));
  }

  hDenominatorRig2->Divide(hProton_MultInt_stat2_3);
  hDenominatorSystRig2->Divide(hProton_MultInt_syst2_3);

  hNumeratorSystRig->SetFillStyle(0);
  hNumerator2SystRig->SetFillStyle(0);
  hNumerator3SystRig->SetFillStyle(0);
  hNumerator4SystRig->SetFillStyle(0);
  hDenominatorSystRig->SetFillStyle(0);
  hDenominatorSystRig2->SetFillStyle(0);

  TH1D* Axis = new TH1D("", "", 100,minpT,maxpT);
  setOPT_hists(Axis,PtOverALabel, B3Lable,505,20,1);
  Axis->GetXaxis()->SetLabelSize(textsize);
  Axis->GetYaxis()->SetLabelSize(textsize);
  Axis->GetYaxis()->SetTitleOffset(1.1);
  Axis->GetYaxis()->SetRangeUser(8e-6,1e-3);

  // 0-10 %
  TCanvas* cNumerator = new TCanvas("cNumerator","0-10 %",800,600);
  gPad->SetTicks();
  gPad->SetLogy();

  TLegend* legSingle = plotLegend("left_top",LegendTitle, 0.5, 0.6, 0.05,-0.05);
  legSingle->SetFillStyle(0);
  legSingle->SetBorderSize(0);
  legSingle->SetTextSize(legsize);
  legSingle->AddEntry((TObject*)0,AnalysisName.Data(),"");
  legSingle->AddEntry((TObject*)0,NucleiType.Data(),"");
  legSingle->AddEntry(hNumeratorRig,"0#minus10%");

  Axis->DrawCopy();
  hNumeratorSystRig->DrawCopy("same E2");
  hNumeratorRig->DrawCopy("same EX0");
  legSingle->Draw("same");
  cNumerator->Update();
//  cNumerator->SaveAs("Single_0_10.eps");

  // 10-20 %
  TCanvas* cNumerator2 = new TCanvas("cNumerator2","10-20 %",800,600);
  gPad->SetTicks();
  gPad->SetLogy();

  legSingle->Clear();
  legSingle->SetHeader(LegendTitle);
  legSingle->AddEntry((TObject*)0,AnalysisName.Data(),"");
  legSingle->AddEntry((TObject*)0,NucleiType.Data(),"");
  legSingle->AddEntry(hNumerator2Rig,"10#minus20%");

  Axis->DrawCopy();
  hNumerator2SystRig->DrawCopy("same E2");
  hNumerator2Rig->DrawCopy("same EX0");
  legSingle->Draw("same");
  cNumerator2->Update();
//  cNumerator2->SaveAs("Single_10_20.eps");

  // 20-40 %
  TCanvas* cNumerator3 = new TCanvas("cNumerator3","20-40 %",800,600);
  gPad->SetTicks();
  gPad->SetLogy();

  legSingle->Clear();
  legSingle->SetHeader(LegendTitle);
  legSingle->AddEntry((TObject*)0,AnalysisName.Data(),"");
  legSingle->AddEntry(hNumerator3Rig,"20#minus40%");

  Axis->DrawCopy();
  hNumerator3SystRig->DrawCopy("same E2");
  hNumerator3Rig->DrawCopy("same EX0");
  legSingle->Draw("same");
  cNumerator3->Update();
//  cNumerator3->SaveAs("Single_20_40.eps");

  // 40-100 %
  TCanvas* cNumerator4 = new TCanvas("cNumerator4","40-100 %",800,600);
  gPad->SetTicks();
  gPad->SetLogy();

  legSingle->Clear();
  legSingle->SetHeader(LegendTitle);
  legSingle->AddEntry((TObject*)0,AnalysisName.Data(),"");
  legSingle->AddEntry((TObject*)0,NucleiType.Data(),"");
  legSingle->AddEntry(hNumerator4Rig,"40#minus100%");

  Axis->DrawCopy();
  hNumerator4SystRig->DrawCopy("same E2");
  hNumerator4Rig->DrawCopy("same EX0");
  legSingle->Draw("same");
  cNumerator4->Update();
//  cNumerator4->SaveAs("Single_40_100.eps");

  // Integrated
  TCanvas* cDenominator = new TCanvas("cDenominator","Multiplicity integrated",800,600);
  gPad->SetTicks();
  gPad->SetLogy();

  legSingle->Clear();
  legSingle->SetHeader(LegendTitle);
  legSingle->AddEntry((TObject*)0,AnalysisName.Data(),"");
  legSingle->AddEntry((TObject*)0,NucleiType.Data(),"");
  legSingle->AddEntry(hDenominator,"Multiplicity integrated");

  Axis->DrawCopy();
  hDenominatorSystRig->DrawCopy("E2same");
  hDenominatorRig->DrawCopy("same EX0");
  legSingle->Draw("same");
  cDenominator->Update();
//  cDenominator->SaveAs("Single_MultInt.eps");

  // Integrated triton
  TCanvas* cDenominator2 = new TCanvas("cDenominator2","Multiplicity integrated (triton)",800,600);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.5, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1, 0.5);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();               // pad1 becomes the current pad
  gPad->SetTicks();
  gPad->SetLogy();

  Double_t legsize2 = 0.1;

  TH1D* Axis2 = new TH1D("", "", 100,minpT,maxpT);
  setOPT_hists(Axis2,PtOverALabel, B3Lable,505,20,1);
  Axis2->GetYaxis()->SetLabelSize(legsize2);
  Axis2->GetYaxis()->SetTitleSize(legsize2+0.025);
  Axis2->GetYaxis()->SetTitleOffset(0.55);
  Axis2->GetYaxis()->SetRangeUser(8e-6,1e-3);
  Axis2->GetYaxis()->SetRangeUser(4e-5,9e-4);

  TLegend* legSingle2 = plotLegend("left_top",LegendTitle, 0.5, 0.7, 0.03,-0.05);
  legSingle2->SetFillStyle(0);
  legSingle2->SetBorderSize(0);
  legSingle2->SetTextSize(legsize2-0.01);
  legSingle2->SetHeader("ALICE");
  legSingle2->AddEntry((TObject*)0,"#kern[-0.105]{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}","");
  legSingle2->AddEntry((TObject*)0,"#kern[-0.175]{#minus1 #scale[1.2]{#leq} #it{y}_{cms} < 0}","");

  TLegend* legSingle2b = plotLegend("mid_top","", 0.5, 0.5, 0.05, -0.1);
  legSingle2b->SetFillStyle(0);
  legSingle2b->SetBorderSize(0);
  legSingle2b->SetTextSize(legsize2-0.01);
  legSingle2b->AddEntry(hDenominatorRig,"(^{3}He + ^{3}#bar{He})/2","p");
  legSingle2b->AddEntry(hDenominatorRig2,"(^{3}H + ^{3}#bar{H})/2","p");

  Axis2->DrawCopy();
  hDenominatorSystRig->DrawCopy("E2same");
  hDenominatorSystRig2->DrawCopy("E2same");
  hDenominatorRig->DrawCopy("same EX0");
  hDenominatorRig2->DrawCopy("same EX0");
  legSingle2->Draw("same");
  legSingle2b->Draw("same");

  pad2->cd();       // pad2 becomes the current pad
  gPad->SetTicks();

  TH1D* RatioAxis = new TH1D("", "", 100,minpT,maxpT);
  setOPT_hists(RatioAxis,PtOverALabel, RatioLabel,505,20,1);
  RatioAxis->GetXaxis()->SetTitleOffset(1.);
  RatioAxis->GetXaxis()->SetTitleSize(legsize2+0.025);
  RatioAxis->GetXaxis()->SetLabelSize(legsize2);
  RatioAxis->GetYaxis()->SetTitleOffset(0.55);
  RatioAxis->GetYaxis()->SetLabelOffset(0.02);
  RatioAxis->GetYaxis()->SetTitleSize(legsize2+0.025);
  RatioAxis->GetYaxis()->SetLabelSize(legsize2);
  RatioAxis->GetYaxis()->SetNdivisions(503);
//  RatioAxis->GetXaxis()->CenterTitle();
  RatioAxis->GetYaxis()->CenterTitle();
  RatioAxis->GetYaxis()->SetRangeUser(0.5,2.5);

  TLine *line = new TLine(minpT,1.,maxpT,1.);
  line->SetLineStyle(2);

  Double_t minpTCoal = 0.475;
  Double_t maxpTCoal = 1.025;
  //Coalescence expectation (Bellini, Kalweit)
  Double_t xCol[2] = {minpTCoal, maxpTCoal};
  Double_t yCol[2] = {1.52, 1.52}; //1.46  1.58
  Double_t exCol[2] = {0, 0};
  Double_t eyCol[2] = {0.06, 0.06};
  TGraphErrors* grCol = new TGraphErrors(2, xCol, yCol, exCol, eyCol);
  setOPT_graph(grCol,pTLabel, RatioLabel,505,20,0.5,kGreen-2);
  grCol->SetFillColor(kGreen-9);
  grCol->SetLineWidth(2);
  grCol->SetLineStyle(7);

  //Coalescence expectation (Sun, Ko, Dönigus)
  Double_t yCol2Body[2] = {1.06097, 1.06097};
  Double_t eyCol2Body[2] = {0.00142385, 0.00142385};
  TGraphErrors* TwoBodyCoal = new TGraphErrors(2, xCol, yCol2Body, exCol, eyCol2Body);
  setOPT_graph(TwoBodyCoal,pTLabel, RatioLabel,505,20,0.5,kViolet-2);
  TwoBodyCoal->SetFillColor(kMagenta-10);
  TwoBodyCoal->SetLineWidth(2);
  TwoBodyCoal->SetLineStyle(5);


  Double_t yCol3Body[2] = {1.14308, 1.14308};
  Double_t eyCol3Body[2] = {0.00335148, 0.00335148};
  TGraphErrors* ThreeBodyCoal = new TGraphErrors(2, xCol, yCol3Body, exCol, eyCol3Body);
  setOPT_graph(ThreeBodyCoal,pTLabel, RatioLabel,505,20,0.5,kOrange+7);
  ThreeBodyCoal->SetFillColor(kOrange-9);
  ThreeBodyCoal->SetLineWidth(2);
  ThreeBodyCoal->SetLineStyle(3);

  TLegend* legRatio = plotLegend("right_top","Coalescence expectations", 0.8, 1.2, -0.11, -0.05);
  legRatio->SetFillStyle(0);
  legRatio->SetBorderSize(0);
  legRatio->SetTextSize(legsize2-0.01);
  legRatio->AddEntry(grCol,"Bellini, Kalweit","lf"); // , #scale[0.7]{arXiv:1807.05894}
  legRatio->AddEntry(TwoBodyCoal,"Sun et al. (Two-body)","fl");
  legRatio->AddEntry(ThreeBodyCoal,"Sun et al. (Three-body) ","fl");
//  legRatio->AddEntry((TObject*)0,"Work in progress","");

  RatioAxis->DrawCopy();
  if (DrawCoal) {
    grCol->Draw("C3same");
    grCol->Draw("LXsame");
    TwoBodyCoal->Draw("C3same");
    TwoBodyCoal->Draw("LXsame");
    ThreeBodyCoal->Draw("C3same");
    ThreeBodyCoal->Draw("LXsame");
  }
  hRatioSyst->DrawCopy("E2same");
  hRatio->DrawCopy("EX0same");
  line->Draw("same");
  if (DrawCoal) legRatio->Draw();
  hRatioTotal->Fit("pol0","N");

  cDenominator2->SaveAs("Single_MultInt_Triton.eps");

  TFile B3File("B3_pPb_5TeV.root","RECREATE");
  B3File.cd();
  hNumeratorRig->Write("B3_3He_0_10");
  hNumeratorSystRig->Write("B3_3He_0_10_Syst");

  hNumerator2Rig->Write("B3_3He_10_20");
  hNumerator2SystRig->Write("B3_3He_10_20_Syst");

  hNumerator3Rig->Write("B3_3He_20_40");
  hNumerator3SystRig->Write("B3_3He_20_40_Syst");

  hNumerator4Rig->Write("B3_3He_40_100");
  hNumerator4SystRig->Write("B3_3He_40_100_Syst");

  hDenominatorRig->Write("B3_3He_MultiplicityIntegrated");
  hDenominatorSystRig->Write("B3_3He_MultiplicityIntegrated_Syst");

  hDenominatorRig2->Write("B3_Triton_MultiplicityIntegrated");
  hDenominatorSystRig2->Write("B3_Triton_MultiplicityIntegrated_Syst");
  B3File.Close();

  // Shift spectra by scaling
  hDenominatorRig->Scale(16);
  hDenominatorSystRig->Scale(16);

  hNumerator4Rig->Scale(8);
  hNumerator4SystRig->Scale(8);

  hNumerator3Rig->Scale(4);
  hNumerator3SystRig->Scale(4);

  hNumerator2Rig->Scale(2);
  hNumerator2SystRig->Scale(2);

  // here you can change the name for the legend
  TString Numeratorname = "0#minus10%";
  TString Denominatorname = "Minimum-bias (x16)";
  TString Numeratorname2 ="10#minus20% (x2)";
  TString Numeratorname3 = "20#minus40% (x4)";
  TString Numeratorname4 = "40#minus100% (x8)";
  if (file3!="") {
    hNumerator2->SetFillStyle(0);
  }
  if (file4!="") {
    hNumerator3->SetFillStyle(0);
  }
  if (file5!="") {
    hNumerator4->SetFillStyle(0);
  }

  TCanvas* c1 = new TCanvas("c1","Summary",800,600);
  gPad->SetRightMargin(0.01);
  gPad->SetTicks();
  gPad->SetLogy();

  hNumerator->GetXaxis()->SetLabelSize(textsize);
  hNumerator->GetYaxis()->SetLabelSize(textsize);
  Axis->GetYaxis()->SetRangeUser(1.2e-5,0.2);

  Axis->DrawCopy();
  hNumeratorSystRig->DrawCopy("same E2");
  hNumeratorRig->DrawCopy("same EX0");

  hDenominatorSystRig->DrawCopy("same E2");
  hDenominatorRig->DrawCopy("same EX0");

  if (file3!="") {
    hNumerator2SystRig->DrawCopy("same E2");
    hNumerator2Rig->DrawCopy("same EX0");
  }
  if (file4!="") {
    hNumerator3SystRig->DrawCopy("same E2");
    hNumerator3Rig->DrawCopy("same EX0");
  }
  if (file5!="") {
    hNumerator4SystRig->DrawCopy("same E2");
    hNumerator4Rig->DrawCopy("same EX0");
  }

  TLegend* legALICE = plotLegend("left_top","ALICE", 1., 0.7,0.01,-0.03);
  legALICE->SetFillStyle(0);
  legALICE->SetBorderSize(0);
  legALICE->SetTextSize(legsize);
  legALICE->AddEntry((TObject*)0,AnalysisName.Data(),"");
  legALICE->AddEntry((TObject*)0,RapidityInfo.Data(),"");
  legALICE->AddEntry((TObject*)0,NucleiType.Data(),"");
  legALICE->Draw();

  TLegend* leg = plotLegend("right_top","", 1., 0.1*(nFiles+1),-0.09,-0.03);
  leg->SetNColumns(2);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(legsize);
  leg->AddEntry((TObject*)0,"#kern[-0.08]{V0A Multiplicity Classes}","");
  leg->AddEntry((TObject*)0," ","");
  leg->AddEntry(hNumeratorSystRig,Numeratorname,"p");
  if (file3!="") {
    leg->AddEntry(hNumerator2SystRig,Numeratorname2,"p");
  }
  if (file4!="") {
    leg->AddEntry(hNumerator3SystRig,Numeratorname3,"p");
  }
  if (file5!="") {
    leg->AddEntry(hNumerator4SystRig,Numeratorname4,"p");
  }
  leg->AddEntry(hDenominatorSystRig,Denominatorname,"p");
  leg->Draw();
  c1->SaveAs("Summary_B3.eps");
}
