#include <stdio.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "../mystyle.C"

using namespace std;

Double_t minpT = 0;
Double_t maxpT = 6;

void MakeMeanPtPlots(){
  
  bool plotTriton = false;
  bool DrawLineFit = false;
  bool CheckMassScaling = false;
  
  TString file1="MeanPt.root";
  
  TString LegendTitle = "ALICE"; // "ALICE Work in progress";
  TString AnalysisName = "#kern[-0.03]{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}";
  Int_t nFiles = 5;
  
  //style area
  Double_t legsize = 0.04;
  Double_t textsize = 0.05;
  Double_t MarkerSize = 1.75;
  mystyle();
  
  TString MultiplicityLabel = "#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{#lbar#it{#eta}_{lab}#lbar < 0.5}";
  TString yLable = "#LT#it{p}_{T}#GT (GeV/#it{c})";
  
  TFile InputFile (file1.Data());
  
  //Helium
  TGraphErrors* gHelium = (TGraphErrors*) InputFile.GetObjectChecked("MeanPt_Helium","TGraphErrors");
  setOPT_graph(gHelium,MultiplicityLabel, yLable,510,20,MarkerSize,kBlue);
  TGraphErrors* gHeliumSyst = (TGraphErrors*) InputFile.GetObjectChecked("MeanPtSyst_Helium","TGraphErrors");
  setOPT_graph(gHeliumSyst,MultiplicityLabel, yLable,510,20,MarkerSize,kBlue);
  gHeliumSyst->SetFillStyle(0);
  
  //Anti-Helium
  TGraphErrors* gAntiHelium = (TGraphErrors*) InputFile.GetObjectChecked("MeanPt_AntiHelium","TGraphErrors");
  setOPT_graph(gAntiHelium,MultiplicityLabel, yLable,510,21,MarkerSize,kGreen-3);
  TGraphErrors* gAntiHeliumSyst = (TGraphErrors*) InputFile.GetObjectChecked("MeanPtSyst_AntiHelium","TGraphErrors");
  setOPT_graph(gAntiHeliumSyst,MultiplicityLabel, yLable,510,21,MarkerSize,kGreen-3);
  gAntiHeliumSyst->SetFillStyle(0);
  
  //Average
  TGraphErrors* gAverage = (TGraphErrors*) InputFile.GetObjectChecked("MeanPt_Average","TGraphErrors");
  setOPT_graph(gAverage,MultiplicityLabel, yLable,510,20,MarkerSize,kRed);
  TGraphErrors* gAverageSyst = (TGraphErrors*) InputFile.GetObjectChecked("MeanPtSyst_Average","TGraphErrors");
  setOPT_graph(gAverageSyst,MultiplicityLabel, yLable,510,20,MarkerSize,kRed);
  gAverageSyst->SetFillStyle(0);
  
  for (Int_t iPoint=0; iPoint<gHelium->GetN(); iPoint++) {
    gHelium->SetPointError(iPoint,0,gHelium->GetErrorY(iPoint));
  }
  for (Int_t iPoint=0; iPoint<gAntiHelium->GetN(); iPoint++) {
    gAntiHelium->SetPointError(iPoint,0,gAntiHelium->GetErrorY(iPoint));
  }
  for (Int_t iPoint=0; iPoint<gAverage->GetN(); iPoint++) {
    gAverage->SetPointError(iPoint,0,gAverage->GetErrorY(iPoint));
  }
  
  TCanvas* cMeanPt = new TCanvas("cMeanPt","Mean p_{T} vs multiplicity",800, 600);
  gPad->SetTicks();
  
  TH1D* Axis = new TH1D("", "", 50,0,50);
  setOPT_hists(Axis,MultiplicityLabel, yLable,510,20,0.5);
  Axis->GetYaxis()->SetRangeUser(1,3);
  Axis->GetXaxis()->SetLabelSize(0.05);
  Axis->GetYaxis()->SetLabelSize(0.05);
  Axis->GetXaxis()->SetTitleOffset(1.05);
  
  TLegend* leg = plotLegend("left_top",LegendTitle, 0.35, 0.6,0.05,-0.03);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->SetTextSize(legsize);
  leg->AddEntry((TObject*)0,AnalysisName.Data(),"");
  leg->AddEntry((TObject*)0,"","");
  leg->AddEntry(gHelium,"^{3}He","p");
  leg->AddEntry(gAverage,"(^{3}He + ^{3}#bar{He})/2","p");
  leg->AddEntry(gAntiHelium,"^{3}#bar{He}","p");
  
  Axis->DrawCopy();
  gHeliumSyst->DrawClone("2same");
  gAntiHeliumSyst->DrawClone("2same");
  gAverageSyst->DrawClone("2same");
  
  gStyle->SetErrorX(0);
  gHelium->DrawClone("Zpsame");
  gAntiHelium->DrawClone("Zpsame");
  gAverage->DrawClone("Zpsame");
  leg->DrawClone();
  cMeanPt->SaveAs("Comparison_MeanPt.eps");
  
  // All measurements
  TCanvas* cOverview = new TCanvas("cOverview","Mean p_{T} vs multiplicity (all measurements)",800, 800);
  gPad->SetTicks();
  gPad->SetLogx();
  gPad->SetBottomMargin(0.175);
  
  TH1D* Axis2 = new TH1D("", "", 4999,1,5000);
  setOPT_hists(Axis2,MultiplicityLabel, yLable,510,20,0.5);
  Axis2->GetYaxis()->SetNdivisions(505);
  Axis2->GetYaxis()->SetRangeUser(1,4.5);
  Axis2->GetXaxis()->SetLabelSize(0.05);
  Axis2->GetYaxis()->SetLabelSize(0.05);
  Axis2->GetXaxis()->SetTitleOffset(1.3);
  Axis2->GetYaxis()->SetTitleOffset(1);
  
  Axis2->GetXaxis()->SetRangeUser(4,2000);
  
  TLegend* leg2 = plotLegend("left_top",LegendTitle, 0.8, 0.95,0.03,-0.03);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(legsize);
  // pp result
  Double_t ppMult = 6;
  Double_t ppMultUnc = 0.2;
  Double_t ppValue = 1.6;
  Double_t ppUnc = 0.4;
  TGraphErrors *gpp = new TGraphErrors(1,&ppMult,&ppValue,0,&ppUnc);
  gpp->SetName("ppHe");
  gpp->SetTitle("ppHe");
  gpp->SetFillColor(kGreen+2);
  gpp->SetLineColor(kGreen+2);
  gpp->SetMarkerColor(kGreen+2);
  gpp->SetMarkerStyle(33);
  gpp->SetMarkerSize(MarkerSize+0.75);
  
  Double_t ppUncSyst = 0.04;
  TGraphErrors *gppSyst = new TGraphErrors(1,&ppMult,&ppValue,&ppMultUnc,&ppUncSyst);
  gppSyst->SetName("ppHeSyst");
  gppSyst->SetTitle("ppHeSyst");
  gppSyst->SetFillColor(kGreen+2);
  gppSyst->SetLineColor(kGreen+2);
  gppSyst->SetMarkerColor(kGreen+2);
  gppSyst->SetMarkerStyle(33);
  gppSyst->SetMarkerSize(MarkerSize+0.75);
  gppSyst->SetFillStyle(0);
  
  // PbPb result 2.76 TeV
  Double_t Mult[2] = {1206.7, 266};
  Double_t MultUnc[2] = {46.67, 9.83};
  Double_t PbPbValue[2] = {2.83, 2.65};
  Double_t PbPbUnc[2] = {0.05, 0.06};
  TGraphErrors *gPbPb2 = new TGraphErrors(2,Mult, PbPbValue, 0, PbPbUnc);
  gPbPb2->SetName("PbPb");
  gPbPb2->SetTitle("PbPb");
  gPbPb2->SetFillColor(kBlue);
  gPbPb2->SetLineColor(kBlue);
  gPbPb2->SetMarkerColor(kBlue);
  gPbPb2->SetMarkerStyle(21);
  gPbPb2->SetMarkerSize(MarkerSize);
  
  Double_t PbPbUncSyst[2] = {0.45, 0.45};
  TGraphErrors *gPbPb2Syst = new TGraphErrors(2,Mult, PbPbValue, MultUnc, PbPbUncSyst);
  gPbPb2Syst->SetName("PbPbSyst");
  gPbPb2Syst->SetTitle("PbPbSyst");
  gPbPb2Syst->SetFillColor(kBlue);
  gPbPb2Syst->SetLineColor(kBlue);
  gPbPb2Syst->SetMarkerColor(kBlue);
  gPbPb2Syst->SetMarkerStyle(21);
  gPbPb2Syst->SetMarkerSize(MarkerSize);
  gPbPb2Syst->SetFillStyle(0);
  leg2->AddEntry( (TObject*) NULL,"#kern[-0.135]{#minus1 #scale[1.2]{#leq} #it{y}_{cms} < 0}","");
  leg2->AddEntry(gAverage,"(^{3}He + ^{3}#bar{He})/2, p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","p");
  leg2->AddEntry( (TObject*) NULL,"#kern[-0.2]{|#it{y_{cms}}| < 0.5}","");
  leg2->AddEntry(gPbPb2,"^{3}He, Pb#minusPb #sqrt{#it{s}_{NN}} = 2.76 TeV","p"); // , #scale[0.7]{Phys.Rev. C93 (2016) 2, 024917}
  leg2->AddEntry(gpp,"^{3}#bar{He}, pp #sqrt{#it{s}} = 7 TeV","p"); // , #scale[0.7]{Phys.Rev. C97 (2018) 2, 024615}


  
  //  // PbPb result 5.02 TeV
  //  Double_t Mult5[3] ={1764, 826, 131.944};
  //  Double_t MultUnc5[3] ={49.75, 22, 6.208};
  //  Double_t PbPbValue5[3] = {3.06, 2.84, 1.99};
  //  Double_t PbPbUnc5[3] = {0.09, 0.09, 0.12};
  //  TGraphErrors *gPbPb5 = new TGraphErrors(3,Mult5, PbPbValue5, 0, PbPbUnc5);
  //  gPbPb5->SetName("PbPb");
  //  gPbPb5->SetTitle("PbPb");
  //  gPbPb5->SetFillColor(1);
  //  gPbPb5->SetLineColor(1);
  //  gPbPb5->SetMarkerColor(1);
  //  gPbPb5->SetMarkerStyle(20);
  
  //  Double_t PbPbUncSyst5[3] = {0.34, 0.28, 0.11};
  //  TGraphErrors *gPbPb5Syst = new TGraphErrors(3,Mult5, PbPbValue5, MultUnc5, PbPbUncSyst5);
  //  gPbPb5Syst->SetName("PbPbSyst");
  //  gPbPb5Syst->SetTitle("PbPbSyst");
  //  gPbPb5Syst->SetFillColor(1);
  //  gPbPb5Syst->SetLineColor(1);
  //  gPbPb5Syst->SetMarkerColor(1);
  //  gPbPb5Syst->SetMarkerStyle(20);
  //  gPbPb5Syst->SetFillStyle(0);
  //  leg2->AddEntry(gPbPb5,"^{3}He, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, #scale[0.7]{Preliminary}","p");
  
  
  // all data points in one graph with total uncertainty
  TGraphErrors *gAll = (TGraphErrors*)gAverageSyst->Clone("gAll");
  for (Int_t iPoint=0; iPoint<gAverage->GetN(); iPoint++) {
    gAll->SetPointError(iPoint, gAverage->GetErrorX(iPoint), sqrt(pow(gAverage->GetErrorY(iPoint),2)+pow(gAverageSyst->GetErrorY(iPoint),2)));
  }
  
  // add pp
  gAll->SetPoint(gAll->GetN(), ppMult, ppValue);
  gAll->SetPointError(gAll->GetN()-1, ppMultUnc, sqrt(pow(ppUnc,2)+pow(ppUncSyst,2)));
  
  // add PbPb result 2.76 TeV
  gAll->SetPoint(gAll->GetN(), Mult[0], PbPbValue[0]);
  gAll->SetPointError(gAll->GetN()-1, MultUnc[0], sqrt(pow(PbPbUnc[0],2)+pow(PbPbUncSyst[0],2)));
  
  gAll->SetPoint(gAll->GetN(), Mult[1], PbPbValue[1]);
  gAll->SetPointError(gAll->GetN()-1, MultUnc[1], sqrt(pow(PbPbUnc[1],2)+pow(PbPbUncSyst[1],2)));
  
  //  // add PbPb result 5 TeV
  //  gAll->SetPoint(gAll->GetN(), Mult5[0], PbPbValue5[0]);
  //  gAll->SetPointError(gAll->GetN()-1, MultUnc5[0], sqrt(pow(PbPbUnc5[0],2)+pow(PbPbUncSyst5[0],2)));
  //
  //  gAll->SetPoint(gAll->GetN(), Mult5[1], PbPbValue5[1]);
  //  gAll->SetPointError(gAll->GetN()-1, MultUnc5[1], sqrt(pow(PbPbUnc5[1],2)+pow(PbPbUncSyst5[1],2)));
  //
  //  gAll->SetPoint(gAll->GetN(), Mult5[2], PbPbValue5[2]);
  //  gAll->SetPointError(gAll->GetN()-1, MultUnc5[2], sqrt(pow(PbPbUnc5[2],2)+pow(PbPbUncSyst5[2],2)));
  
  Axis2->DrawCopy();
  //  gAll->DrawClone("2same");
  gAverageSyst->DrawClone("2same");
  gppSyst->DrawClone("2same");
  gPbPb2Syst->DrawClone("2same");
  //  gPbPb5Syst->DrawClone("2same");
  
  if (plotTriton) {
    
    //3H pPb result (24.05.2019
    Double_t TritonMult = 17.8;
    Double_t TritonMultUnc = 0.4;
    Double_t TritonValue = 2.81;
    Double_t TritonUnc = 0.43;
    TGraphErrors *Triton = new TGraphErrors(1,&TritonMult,&TritonValue,0,&TritonUnc);
    Triton->SetName("Triton");
    Triton->SetTitle("Triton");
    Triton->SetFillColor(kRed+2);
    Triton->SetLineColor(kRed+2);
    Triton->SetMarkerColor(kRed+2);
    Triton->SetMarkerStyle(24);
    
    TritonUnc = 0.67; // Syst. unc.
    TGraphErrors *TritonSyst = new TGraphErrors(1,&TritonMult,&TritonValue,&TritonMultUnc,&TritonUnc);
    TritonSyst->SetName("TritonSyst");
    TritonSyst->SetTitle("TritonSyst");
    TritonSyst->SetFillColor(kRed+2);
    TritonSyst->SetLineColor(kRed+2);
    TritonSyst->SetMarkerColor(kRed+2);
    TritonSyst->SetMarkerStyle(24);
    TritonSyst->SetFillStyle(0);
    
    leg2->AddEntry(Triton,"(^{3}H + ^{3}#bar{H})/2, p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV","p");
    
    TritonSyst->Draw("2same");
    Triton->Draw("Zpsame");
  }
  
  gAverage->DrawClone("Zpsame");
  gpp->DrawClone("Zpsame");
  gPbPb2->DrawClone("Zpsame");
  //  gPbPb5->DrawClone("Zpsame");
  //
  
  leg2->DrawClone();
  cOverview->SaveAs("Summary_MeanPt.eps");
  
  if (DrawLineFit) {
    TF1* Param = new TF1("Param","[0]+[1]*TMath::Log10(x)", 0,10000);
    Param->SetLineColor(1);
    Param->SetLineStyle(2);
    const Double_t StartParam[3] = {1.2, 0.5};
    Param->SetParameters(StartParam);
    //  Param->FixParameter(0,0);
    gAll->Fit(Param,"0EM");
    Param->Draw("same");
    printf("chi2/ndf: %.2f / %d \n", Param->GetChisquare(), Param->GetNDF() );
  }
  
  
  // <pt>/m vs multiplicity
  TCanvas* cComparisonParticles = new TCanvas("cComparisonParticles","<p_{T}>/m vs multiplicity for different particles ",800, 800);
  cComparisonParticles->SetTicks();
  cComparisonParticles->SetLogx();
  
  // Merge all collision systems in one graph
  TGraphErrors *gAll_Stat = (TGraphErrors*)gAverage->Clone("gAll_Stat");
  TGraphErrors *gAll_Syst = (TGraphErrors*)gAverageSyst->Clone("gAll_Syst");

//  // add pp
//  gAll_Stat->SetPoint(gAll_Stat->GetN(), ppMult, ppValue);
//  gAll_Stat->SetPointError(gAll_Stat->GetN()-1, ppMultUnc, ppUnc);
//  gAll_Syst->SetPoint(gAll_Syst->GetN(), ppMult, ppValue);
//  gAll_Syst->SetPointError(gAll_Syst->GetN()-1, ppMultUnc, ppUncSyst);
//
//  // add PbPb result 2.76 TeV
//  gAll_Stat->SetPoint(gAll_Stat->GetN(), Mult[0], PbPbValue[0]);
//  gAll_Stat->SetPointError(gAll_Stat->GetN()-1, MultUnc[0], PbPbUnc[0]);
//
//  gAll_Syst->SetPoint(gAll_Syst->GetN(), Mult[0], PbPbValue[0]);
//  gAll_Syst->SetPointError(gAll_Syst->GetN()-1, MultUnc[0], PbPbUncSyst[0]);
//
//  gAll_Stat->SetPoint(gAll_Stat->GetN(), Mult[1], PbPbValue[1]);
//  gAll_Stat->SetPointError(gAll_Stat->GetN()-1, MultUnc[1], PbPbUnc[1]);
//
//  gAll_Syst->SetPoint(gAll_Syst->GetN(), Mult[1], PbPbValue[1]);
//  gAll_Syst->SetPointError(gAll_Syst->GetN()-1, MultUnc[1], PbPbUncSyst[1]);
  
  setOPT_graph(gAll_Syst, MultiplicityLabel, yLable,510,20,MarkerSize,kRed);
  setOPT_graph(gAll_Stat, MultiplicityLabel, yLable,510,20,MarkerSize,kRed);
  
  if (CheckMassScaling) {
    // Divide by the 3He mass (2.810 GeV/c^2)
    Double_t He3Mass = 2.810;
    Double_t x,y;
    for (int ipoint3He=0; ipoint3He<gAll_Stat->GetN() ; ipoint3He++) {
      gAll_Stat->GetPoint(ipoint3He,x,y);
      gAll_Stat->SetPoint(ipoint3He, x, y/He3Mass);
      gAll_Stat->SetPointError(ipoint3He, gAll_Stat->GetErrorX(ipoint3He), gAll_Stat->GetErrorY(ipoint3He)/He3Mass);
      
      gAll_Syst->SetPoint(ipoint3He, x, y/He3Mass);
      gAll_Syst->SetPointError(ipoint3He, gAll_Syst->GetErrorX(ipoint3He), gAll_Syst->GetErrorY(ipoint3He)/He3Mass);
    }
  }
  // deuteron <pT>
  constexpr int kDeutBin = 5;
  double mult[kDeutBin]{40.6,30.5,23.2,16.1,7.1};
  double mult_err[kDeutBin]{0.9,0.7,0.5,0.4,0.2};
  double deuteron_mean_pt[kDeutBin]{1.55,1.44,1.38,1.23,1.07};
  double deuteron_mean_pt_stat[kDeutBin]{0.0122,0.0134,0.0103,0.0046,0.0113};
  double deuteron_mean_pt_syst[kDeutBin]{0.0856,0.0931,0.0789,0.0668,0.0712};
  
  if (CheckMassScaling) {
    // Divide by the deuteron mass (1.876 GeV/c^2)
    Double_t DeuteronMass = 1.876;
    for (int ipoint=0; ipoint < kDeutBin; ipoint++) {
      deuteron_mean_pt[ipoint] = deuteron_mean_pt[ipoint]/DeuteronMass;
      deuteron_mean_pt_stat[ipoint] = deuteron_mean_pt_stat[ipoint]/DeuteronMass;
      deuteron_mean_pt_syst[ipoint] = deuteron_mean_pt_syst[ipoint]/DeuteronMass;
    }
  }
  
  TGraphErrors *gDeuteron_stat = new TGraphErrors(kDeutBin, mult, deuteron_mean_pt, 0, deuteron_mean_pt_stat);
  setOPT_graph(gDeuteron_stat, MultiplicityLabel, yLable,510,21,MarkerSize,kBlue);
  TGraphErrors *gDeuteron_syst = new TGraphErrors(kDeutBin, mult, deuteron_mean_pt, mult_err, deuteron_mean_pt_syst);
  setOPT_graph(gDeuteron_syst, MultiplicityLabel, yLable,510,21,MarkerSize,kBlue);
  gDeuteron_syst->SetFillStyle(0);
  
  // Data taken from arXiv:1307.6796
  // proton <pT>
  constexpr int kProtBin = 7;
  double mult2[kProtBin]{45,36.2,30.5,23.2,16.1,9.8,4.4};
  double mult2_err[kProtBin]{1,0.8,0.7,0.5,0.4,0.2,0.1};
  double proton_mean_pt[kProtBin]{1.248,1.223,1.186,1.1132,1.053,0.9607,0.8208};
  double proton_mean_pt_stat[kProtBin]{0.004,0.004,0.004,0.003,0.002,0.0024,0.0029};
  double proton_mean_pt_syst[kProtBin]{0.043,0.049,0.054,0.044,0.036,0.0304,0.0238};
  
  if (CheckMassScaling) {
    // Divide by the proton mass (0.938 GeV/c^2)
    Double_t ProtonMass = 0.938;
    for (int ipoint=0; ipoint < kProtBin; ipoint++) {
      proton_mean_pt[ipoint] = proton_mean_pt[ipoint]/ProtonMass;
      proton_mean_pt_stat[ipoint] = proton_mean_pt_stat[ipoint]/ProtonMass;
      proton_mean_pt_syst[ipoint] = proton_mean_pt_syst[ipoint]/ProtonMass;
    }
  }
  
  TGraphErrors *gProton_stat = new TGraphErrors(kProtBin, mult2, proton_mean_pt, 0, proton_mean_pt_stat);
  setOPT_graph(gProton_stat, MultiplicityLabel, yLable,510,33,MarkerSize+0.75,kGreen-3);
  TGraphErrors *gProton_syst = new TGraphErrors(kProtBin, mult2, proton_mean_pt, mult2_err, proton_mean_pt_syst);
  setOPT_graph(gProton_syst, MultiplicityLabel, yLable,510,33,MarkerSize+0.75,kGreen-3);
  gProton_syst->SetFillStyle(0);
  
  // kaon <pT>
  double kaon_mean_pt[kProtBin]{0.9366,0.9177,0.901,0.8764,0.8317,0.7722,0.6809};
  double kaon_mean_pt_stat[kProtBin]{0.0048,0.0043,0.003,0.0024,0.0023,0.0022,0.0026};
  double kaon_mean_pt_syst[kProtBin]{0.0962,0.0789,0.0743,0.0702,0.063,0.0521,0.0352};
  
  if (CheckMassScaling) {
    // Divide by the kaon mass (0.498 GeV/c^2)
    Double_t KaonMass = 0.498;
    for (int ipoint=0; ipoint < kProtBin; ipoint++) {
      kaon_mean_pt[ipoint] = kaon_mean_pt[ipoint]/KaonMass;
      kaon_mean_pt_stat[ipoint] = kaon_mean_pt_stat[ipoint]/KaonMass;
      kaon_mean_pt_syst[ipoint] = kaon_mean_pt_syst[ipoint]/KaonMass;
    }
  }
  
  TGraphErrors *gKaon_stat = new TGraphErrors(kProtBin, mult2, kaon_mean_pt, 0, kaon_mean_pt_stat);
  setOPT_graph(gKaon_stat, MultiplicityLabel, yLable,510,34,MarkerSize+0.25,kMagenta-3);
  TGraphErrors *gKaon_syst = new TGraphErrors(kProtBin, mult2, kaon_mean_pt, mult2_err, kaon_mean_pt_syst);
  setOPT_graph(gKaon_syst, MultiplicityLabel, yLable,510,34,MarkerSize+0.25,kMagenta-3);
  gKaon_syst->SetFillStyle(0);
  
  // pion <pT>
  double pion_mean_pt[kProtBin]{0.5453,0.5375,0.5293,0.5142,0.4944,0.4705,0.4336};
  double pion_mean_pt_stat[kProtBin]{0.0004,0.0005,0.0003,0.0003,0.0003,0.0004,0.0005};
  double pion_mean_pt_syst[kProtBin]{0.0154,0.0147,0.0145,0.0134,0.013,0.012,0.0098};
  
  if (CheckMassScaling) {
    // Divide by the pion mass (0.140 GeV/c^2)
    Double_t pionMass = 0.140;
    Double_t PlottingFactor = 0.6;
    for (int ipoint=0; ipoint < kProtBin; ipoint++) {
      pion_mean_pt[ipoint] = pion_mean_pt[ipoint]/pionMass*PlottingFactor;
      pion_mean_pt_stat[ipoint] = pion_mean_pt_stat[ipoint]/pionMass*PlottingFactor;
      pion_mean_pt_syst[ipoint] = pion_mean_pt_syst[ipoint]/pionMass*PlottingFactor;
    }
  }
  
  TGraphErrors *gpion_stat = new TGraphErrors(kProtBin, mult2, pion_mean_pt, 0, pion_mean_pt_stat);
  setOPT_graph(gpion_stat, MultiplicityLabel, yLable,510,43,MarkerSize+1,kOrange-3);
  TGraphErrors *gpion_syst = new TGraphErrors(kProtBin, mult2, pion_mean_pt, mult2_err, pion_mean_pt_syst);
  setOPT_graph(gpion_syst, MultiplicityLabel, yLable,510,43,MarkerSize+1,kOrange-3);
  gpion_syst->SetFillStyle(0);
  
  // lambda <pT>
  double lambda_mean_pt[kProtBin]{1.378,1.254,1.324,1.271,1.197,1.098,0.9711};
  double lambda_mean_pt_stat[kProtBin]{0.005,0.005,0.004,0.004,0.004,0.004,0.0056};
  double lambda_mean_pt_syst[kProtBin]{0.04,0.037,0.039,0.032,0.029,0.036,0.0463};
  
  if (CheckMassScaling) {
    // Divide by the lambda mass (1.116 GeV/c^2)
    Double_t lambdaMass = 1.116;
    for (int ipoint=0; ipoint < kProtBin; ipoint++) {
      lambda_mean_pt[ipoint] = lambda_mean_pt[ipoint]/lambdaMass;
      lambda_mean_pt_stat[ipoint] = lambda_mean_pt_stat[ipoint]/lambdaMass;
      lambda_mean_pt_syst[ipoint] = lambda_mean_pt_syst[ipoint]/lambdaMass;
    }
  }
  
  TGraphErrors *glambda_stat = new TGraphErrors(kProtBin, mult2, lambda_mean_pt, 0, lambda_mean_pt_stat);
  setOPT_graph(glambda_stat, MultiplicityLabel, yLable,510,29,MarkerSize+0.75,kGreen+3);
  TGraphErrors *glambda_syst = new TGraphErrors(kProtBin, mult2, lambda_mean_pt, mult2_err, lambda_mean_pt_syst);
  setOPT_graph(glambda_syst, MultiplicityLabel, yLable,510,29,MarkerSize+0.75,kGreen+3);
  glambda_syst->SetFillStyle(0);
  
  TH1D* Axis3 = new TH1D("", "", 4999,1,5000);
  setOPT_hists(Axis3,MultiplicityLabel,yLable,510,20,0.5);
  if (CheckMassScaling) {
    setOPT_hists(Axis3,MultiplicityLabel,"#LT#it{p}_{T}#GT/#it{m} (#it{c})",510,20,0.5);
  }
  Axis3->GetYaxis()->SetNdivisions(505);
  Axis3->GetXaxis()->SetLabelSize(0.05);
  Axis3->GetYaxis()->SetLabelSize(0.05);
  Axis3->GetXaxis()->SetTitleOffset(1.3);
  Axis3->GetYaxis()->SetTitleOffset(1.);
  Axis3->GetYaxis()->SetRangeUser(0,3);
  if (CheckMassScaling) Axis3->GetYaxis()->SetRangeUser(0,3.5);
  Axis3->GetXaxis()->SetRangeUser(3,60); // only pp and p–Pb
  
  cComparisonParticles->SetBottomMargin(0.175);
  
  Axis3->DrawCopy();
  gDeuteron_syst->DrawClone("2same");
  gProton_syst->DrawClone("2same");
  gKaon_syst->DrawClone("2same");
  gpion_syst->DrawClone("2same");
  glambda_syst->DrawClone("2same");
  gAll_Syst->DrawClone("2same");
  
  gDeuteron_stat->DrawClone("Zpsame");
  gProton_stat->DrawClone("Zpsame");
  gKaon_stat->DrawClone("Zpsame");
  gpion_stat->DrawClone("Zpsame");
  glambda_stat->DrawClone("Zpsame");
  gAll_Stat->DrawClone("Zpsame");
  
  TLegend* leg3 = plotLegend("left_top",LegendTitle, 0.9, 0.7,0.03,-0.03);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->SetTextSize(legsize);
  leg3->SetNColumns(2);
  leg3->AddEntry((TObject*)nullptr,AnalysisName,"");
  leg3->AddEntry((TObject*)nullptr,"","");

  if (CheckMassScaling) {
    leg3->AddEntry(gpion_stat,"#pi^{+} + #pi^{#minus} (x0.6)","p");
  } else {
    leg3->AddEntry(gpion_stat,"#pi^{+} + #pi^{#minus}","p");
  }
  leg3->AddEntry(gKaon_stat,"K^{+} + K^{#minus}","p");
  leg3->AddEntry(gProton_stat,"p + #bar{p}","p");
  leg3->AddEntry(glambda_stat,"#Lambda + #bar{#Lambda}","p");
  leg3->AddEntry(gDeuteron_stat,"d + #bar{d}","p");
  leg3->AddEntry(gAll_Stat,"^{3}He + ^{3}#bar{He}","p");
  leg3->Draw();
  
  if (CheckMassScaling) {
    cComparisonParticles->SaveAs("MeanPtOverMass.eps");
  } else {
    cComparisonParticles->SaveAs("MeanPtDifferentSpecies.eps");
  }
}
