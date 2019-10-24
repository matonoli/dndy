#include <stdio.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "../mystyle.C"

using namespace std;

void MakeMeanPtvsMassPlots(){
  
  bool DrawLineFit = true;
  
  TString file1="MeanPt.root";
  
  TString LegendTitle = "ALICE"; // "ALICE Work in progress";
  TString AnalysisName = "#kern[-0.15]{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}";
  TString Rapidity = "#kern[-0.25]{#minus1 #scale[1.2]{#leq} #it{y}_{cms} < 0}";
  Int_t nFiles = 5;
  
  //style area
  Double_t legsize = 0.035;
  Double_t textsize = 0.05;
  Double_t MarkerSize = 1.5;
  mystyle();
  
  TString MultiplicityLabel = "#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{#lbar#it{#eta}_{lab}#lbar < 0.5}";
  
  TString MassLabel = "#it{m} (GeV/#it{c}^{2})";
  TString yLable = "#LT#it{p}_{T}#GT (GeV/#it{c})";
  
  TFile InputFile (file1.Data());
  
  //Average
  TGraphErrors* gAverage = (TGraphErrors*) InputFile.GetObjectChecked("MeanPt_Average","TGraphErrors");
  setOPT_graph(gAverage,MultiplicityLabel, yLable,510,20,MarkerSize,kRed);
  TGraphErrors* gAverageSyst = (TGraphErrors*) InputFile.GetObjectChecked("MeanPtSyst_Average","TGraphErrors");
  setOPT_graph(gAverageSyst,MultiplicityLabel, yLable,510,20,MarkerSize,kRed);
  gAverageSyst->SetFillStyle(0);
  
  // Create TGraphs wit <pT> vs mass for different multiplicities
  Double_t pionMass = 0.140;
  Double_t KaonMass = 0.498;
  Double_t ProtonMass = 0.938;
  Double_t lambdaMass = 1.116;
  Double_t DeuteronMass = 1.876;
  Double_t He3Mass = 2.810;
  
  // deuteron <pT>
  constexpr int kDeutBin = 5;
  double mult[kDeutBin]{40.6,30.5,23.2,16.1,7.1};
  double mult_err[kDeutBin]{0.9,0.7,0.5,0.4,0.2};
  double deuteron_mean_pt[kDeutBin]{1.55,1.44,1.38,1.23,1.07};
  double deuteron_mean_pt_stat[kDeutBin]{0.0122,0.0134,0.0103,0.0046,0.0113};
  double deuteron_mean_pt_syst[kDeutBin]{0.0856,0.0931,0.0789,0.0668,0.0712};

  // Data taken from arXiv:1307.6796
  // proton <pT>
  constexpr int kProtBin = 7;
  double mult2[kProtBin]{45,36.2,30.5,23.2,16.1,9.8,4.4};
  double mult2_err[kProtBin]{1,0.8,0.7,0.5,0.4,0.2,0.1};
  double proton_mean_pt[kProtBin]{1.248,1.223,1.186,1.1132,1.053,0.9607,0.8208};
  double proton_mean_pt_stat[kProtBin]{0.004,0.004,0.004,0.003,0.002,0.0024,0.0029};
  double proton_mean_pt_syst[kProtBin]{0.043,0.049,0.054,0.044,0.036,0.0304,0.0238};
  // kaon <pT>
  double kaon_mean_pt[kProtBin]{0.9366,0.9177,0.901,0.8764,0.8317,0.7722,0.6809};
  double kaon_mean_pt_stat[kProtBin]{0.0048,0.0043,0.003,0.0024,0.0023,0.0022,0.0026};
  double kaon_mean_pt_syst[kProtBin]{0.0962,0.0789,0.0743,0.0702,0.063,0.0521,0.0352};
  // pion <pT>
  double pion_mean_pt[kProtBin]{0.5453,0.5375,0.5293,0.5142,0.4944,0.4705,0.4336};
  double pion_mean_pt_stat[kProtBin]{0.0004,0.0005,0.0003,0.0003,0.0003,0.0004,0.0005};
  double pion_mean_pt_syst[kProtBin]{0.0154,0.0147,0.0145,0.0134,0.013,0.012,0.0098};
  // lambda <pT>
  double lambda_mean_pt[kProtBin]{1.378,1.254,1.324,1.271,1.197,1.098,0.9711};
  double lambda_mean_pt_stat[kProtBin]{0.005,0.005,0.004,0.004,0.004,0.004,0.0056};
  double lambda_mean_pt_syst[kProtBin]{0.04,0.037,0.039,0.032,0.029,0.036,0.0463};
  
  
  const int NParticles = 6;
  double masses[NParticles] = {pionMass,KaonMass,ProtonMass, lambdaMass, DeuteronMass, He3Mass};
  double massesError[NParticles] = {0.1*pionMass,0.1*KaonMass,0.1*ProtonMass, 0.1*lambdaMass, 0.1*DeuteronMass, 0.1*He3Mass}; // dummy errors on the masses to have a width for the syst. unc. box

  // mult = 9.8±0.2 / 10.1 (3He + deuteron)
  // no deuterons, multiplicity of 3He only approximatly the same
  double x,y,stat,syst;
  gAverage->GetPoint(3,x,y);
  stat = gAverage->GetErrorY(3);
  syst = gAverage->GetErrorY(3);
  double Mult10[NParticles] = {pion_mean_pt[5], kaon_mean_pt[5], proton_mean_pt[5], lambda_mean_pt[5], (deuteron_mean_pt[3]+2*deuteron_mean_pt[4])/3., y};
  double Mult10Stat[NParticles] = {pion_mean_pt_stat[5], kaon_mean_pt_stat[5], proton_mean_pt_stat[5], lambda_mean_pt_stat[5], (deuteron_mean_pt_stat[3]+2*deuteron_mean_pt_stat[4])/3., stat};
  double Mult10Syst[NParticles] = {pion_mean_pt_syst[5], kaon_mean_pt_syst[5], proton_mean_pt_syst[5], lambda_mean_pt_syst[5], (deuteron_mean_pt_syst[3]+2*deuteron_mean_pt_syst[4])/3., syst};
  
  // total uncertainty
  double Mult10Tot[NParticles] = { sqrt(pow(pion_mean_pt_stat[5],2) + pow(pion_mean_pt_syst[5],2)),
    sqrt(pow(kaon_mean_pt_stat[5],2)+pow(kaon_mean_pt_syst[5],2)),
    sqrt(pow(proton_mean_pt_stat[5],2)+pow(proton_mean_pt_syst[5],2)),
    sqrt(pow(lambda_mean_pt_stat[5],2)+pow(lambda_mean_pt_syst[5],2)),
    (sqrt(pow(deuteron_mean_pt_stat[3],2)+pow(deuteron_mean_pt_syst[3],2))+2*sqrt(pow(deuteron_mean_pt_stat[4],2)+pow(deuteron_mean_pt_syst[4],2)))/3.,
     sqrt(pow(stat,2)+pow(syst,2))};

  TGraphErrors* gMult10 = new TGraphErrors(NParticles, masses, Mult10, nullptr, Mult10Stat);
  TGraphErrors* gMult10Syst = new TGraphErrors(NParticles,masses, Mult10, massesError, Mult10Syst);
  setOPT_graph(gMult10,MassLabel, yLable,510,20,MarkerSize);
  setOPT_graph(gMult10Syst,MassLabel, yLable,510,20,MarkerSize);
  gMult10Syst->SetFillStyle(0);
  
  // mult = 23.2±0.5
  gAverage->GetPoint(2,x,y);
  stat = gAverage->GetErrorY(2);
  syst = gAverage->GetErrorY(2);
  double Mult23[NParticles] = {pion_mean_pt[3], kaon_mean_pt[3], proton_mean_pt[3], lambda_mean_pt[3], deuteron_mean_pt[2], y};
  double Mult23Stat[NParticles] = {pion_mean_pt_stat[3], kaon_mean_pt_stat[3], proton_mean_pt_stat[3], lambda_mean_pt_stat[3], deuteron_mean_pt_stat[2], stat};
  double Mult23Syst[NParticles] = {pion_mean_pt_syst[3], kaon_mean_pt_syst[3], proton_mean_pt_syst[3], lambda_mean_pt_syst[3], deuteron_mean_pt_syst[2], syst};

  double Mult23Tot[NParticles] = { sqrt(pow(pion_mean_pt_stat[3],2) + pow(pion_mean_pt_syst[3],2)),
    sqrt(pow(kaon_mean_pt_stat[3],2)+pow(kaon_mean_pt_syst[3],2)),
    sqrt(pow(proton_mean_pt_stat[3],2)+pow(proton_mean_pt_syst[3],2)),
    sqrt(pow(lambda_mean_pt_stat[3],2)+pow(lambda_mean_pt_syst[3],2)),
    sqrt(pow(deuteron_mean_pt_stat[2],2)+pow(deuteron_mean_pt_syst[2],2)),
    sqrt(pow(stat,2)+pow(syst,2))};
  
  TGraphErrors* gMult23 = new TGraphErrors(NParticles, masses, Mult23, nullptr, Mult23Stat);
  TGraphErrors* gMult23Syst = new TGraphErrors(NParticles,masses, Mult23, massesError, Mult23Syst);
  setOPT_graph(gMult23,MassLabel, yLable,510,21,MarkerSize,kRed);
  setOPT_graph(gMult23Syst,MassLabel, yLable,510,21,MarkerSize,kRed);
  gMult23Syst->SetFillStyle(0);
  
  // mult = 30.5±0.7
  gAverage->GetPoint(1,x,y);
  stat = gAverage->GetErrorY(1);
  syst = gAverage->GetErrorY(1);
  double Mult30[NParticles] = {pion_mean_pt[2], kaon_mean_pt[2], proton_mean_pt[2], lambda_mean_pt[2], deuteron_mean_pt[1], y};
  double Mult30Stat[NParticles] = {pion_mean_pt_stat[2], kaon_mean_pt_stat[2], proton_mean_pt_stat[2], lambda_mean_pt_stat[2], deuteron_mean_pt_stat[1], stat};
  double Mult30Syst[NParticles] = {pion_mean_pt_syst[2], kaon_mean_pt_syst[2], proton_mean_pt_syst[2], lambda_mean_pt_syst[2], deuteron_mean_pt_syst[1], syst};
  
  double Mult30Tot[NParticles] = { sqrt(pow(pion_mean_pt_stat[2],2) + pow(pion_mean_pt_syst[2],2)),
    sqrt(pow(kaon_mean_pt_stat[2],2)+pow(kaon_mean_pt_syst[2],2)),
    sqrt(pow(proton_mean_pt_stat[2],2)+pow(proton_mean_pt_syst[2],2)),
    sqrt(pow(lambda_mean_pt_stat[2],2)+pow(lambda_mean_pt_syst[2],2)),
    sqrt(pow(deuteron_mean_pt_stat[1],2)+pow(deuteron_mean_pt_syst[1],2)),
    sqrt(pow(stat,2)+pow(syst,2))};
  
  TGraphErrors* gMult30 = new TGraphErrors(NParticles, masses, Mult30, nullptr, Mult30Stat);
  TGraphErrors* gMult30Syst = new TGraphErrors(NParticles,masses, Mult30, massesError, Mult30Syst);
  setOPT_graph(gMult30,MassLabel, yLable,510,33,MarkerSize+1.,kBlue);
  setOPT_graph(gMult30Syst,MassLabel, yLable,510,33,MarkerSize+1.,kBlue);
  gMult30Syst->SetFillStyle(0);
  
  // mult = 40.6±0.9
  gAverage->GetPoint(0,x,y);
  stat = gAverage->GetErrorY(0);
  syst = gAverage->GetErrorY(0);
  double Mult41[NParticles] = {(pion_mean_pt[1]+pion_mean_pt[0])/2., (kaon_mean_pt[1]+kaon_mean_pt[0])/2., (proton_mean_pt[1]+proton_mean_pt[0])/2., (lambda_mean_pt[1]+lambda_mean_pt[0])/2, deuteron_mean_pt[0], y};
  double Mult41Stat[NParticles] = {(pion_mean_pt_stat[1]+pion_mean_pt_stat[0])/2., (kaon_mean_pt_stat[1]+kaon_mean_pt_stat[0])/2., (proton_mean_pt_stat[1]+proton_mean_pt_stat[0])/2., (lambda_mean_pt_stat[1]+lambda_mean_pt_stat[0])/2, deuteron_mean_pt_stat[0], stat};
  double Mult41Syst[NParticles] = {(pion_mean_pt_syst[1]+pion_mean_pt_syst[0])/2., (kaon_mean_pt_syst[1]+kaon_mean_pt_syst[0])/2., (proton_mean_pt_syst[1]+proton_mean_pt_syst[0])/2., (lambda_mean_pt_syst[1]+lambda_mean_pt_syst[0])/2, deuteron_mean_pt_syst[0], syst};
  
  double Mult41Tot[NParticles] = { (sqrt(pow(pion_mean_pt_stat[1],2) + pow(pion_mean_pt_syst[1],2))+sqrt(pow(pion_mean_pt_stat[0],2) + pow(pion_mean_pt_syst[0],2)))/2.,
    (sqrt(pow(kaon_mean_pt_stat[1],2) + pow(kaon_mean_pt_syst[1],2))+sqrt(pow(kaon_mean_pt_stat[0],2) + pow(kaon_mean_pt_syst[0],2)))/2.,
    (sqrt(pow(proton_mean_pt_stat[1],2) + pow(proton_mean_pt_syst[1],2))+sqrt(pow(proton_mean_pt_stat[0],2) + pow(proton_mean_pt_syst[0],2)))/2.,
    (sqrt(pow(lambda_mean_pt_stat[1],2) + pow(lambda_mean_pt_syst[1],2))+sqrt(pow(lambda_mean_pt_stat[0],2) + pow(lambda_mean_pt_syst[0],2)))/2.,
    sqrt(pow(deuteron_mean_pt_stat[0],2)+pow(deuteron_mean_pt_syst[0],2)),
    sqrt(pow(stat,2)+pow(syst,2))};
  
  TGraphErrors* gMult41 = new TGraphErrors(NParticles, masses, Mult41, nullptr, Mult41Stat);
  TGraphErrors* gMult41Syst = new TGraphErrors(NParticles,masses, Mult41, massesError, Mult41Syst);
  setOPT_graph(gMult41,MassLabel, yLable,510,34,MarkerSize,kMagenta-3);
  setOPT_graph(gMult41Syst,MassLabel, yLable,510,34,MarkerSize,kMagenta-3);
  gMult41Syst->SetFillStyle(0);
  
  // <pt> vs m
  TCanvas* cMeanPtvsMass = new TCanvas("cMeanPtvsMass","<p_{T}> vs mass",800, 800);
  cMeanPtvsMass->SetTicks();

  TH1D* Axis3 = new TH1D("", "", 300,0,3);
  setOPT_hists(Axis3,MassLabel,yLable,510,20,0.5);
  Axis3->GetYaxis()->SetNdivisions(505);
  Axis3->GetXaxis()->SetLabelSize(0.05);
  Axis3->GetYaxis()->SetLabelSize(0.05);
  Axis3->GetXaxis()->SetTitleOffset(1.2);
  Axis3->GetYaxis()->SetTitleOffset(1.);
  Axis3->GetYaxis()->SetRangeUser(0,4);

  Axis3->DrawCopy();
  gMult10Syst->DrawClone("2same");
  gMult10->DrawClone("Zpsame");

  gMult23Syst->DrawClone("2same");
  gMult23->DrawClone("Zpsame");

  gMult30Syst->DrawClone("2same");
  gMult30->DrawClone("Zpsame");

  gMult41Syst->DrawClone("2same");
  gMult41->DrawClone("Zpsame");
  
  TLegend* leg = plotLegend("left_top",LegendTitle, 0.9, 0.9, 0.02, -0.02);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(legsize);
  leg->AddEntry((TObject*)nullptr,AnalysisName,"");
  leg->AddEntry((TObject*)nullptr,Rapidity,"");
  leg->AddEntry(gMult10,"#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{#lbar#it{#eta}_{lab}#lbar < 0.5} = 9.8#pm0.2 (10.1#pm0.2 for d, ^{3}He)","p");
  leg->AddEntry(gMult23,"#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{#lbar#it{#eta}_{lab}#lbar < 0.5} = 23.2#pm0.5","p");
  leg->AddEntry(gMult30,"#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{#lbar#it{#eta}_{lab}#lbar < 0.5} = 30.5#pm0.7","p");
  leg->AddEntry(gMult41,"#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{#lbar#it{#eta}_{lab}#lbar < 0.5} = 40.6#pm0.9","p");
  leg->Draw();

  if (DrawLineFit) {
    
    // Graphs with total uncertainty - sqrt(stat^2 + syst^2)
    TGraphErrors* gMult10Total = new TGraphErrors(NParticles,masses, Mult10, massesError, Mult10Tot);
    TGraphErrors* gMult23Total = new TGraphErrors(NParticles,masses, Mult23, massesError, Mult23Tot);
    TGraphErrors* gMult30Total = new TGraphErrors(NParticles,masses, Mult30, massesError, Mult30Tot);
    TGraphErrors* gMult41Total = new TGraphErrors(NParticles,masses, Mult41, massesError, Mult41Tot);
    
    TF1* Param = new TF1("Param","pol1", 0,3);
    Param->SetLineStyle(2);
    Param->SetLineColor(1);
    gMult10Total->Fit(Param,"0EM","",0,1.5);
    Param->DrawClone("same");
    printf("10: chi2/ndf: %.2f / %d = %.2f \n", Param->GetChisquare(), Param->GetNDF(), Param->GetChisquare()/Param->GetNDF());
    Param->SetLineColor(kRed);
    gMult23Total->Fit(Param,"0EM","",0,1.5);
    Param->DrawClone("same");
    printf("23: chi2/ndf: %.2f / %d = %.2f \n", Param->GetChisquare(), Param->GetNDF(), Param->GetChisquare()/Param->GetNDF());
    Param->SetLineColor(kBlue);
    gMult30Total->Fit(Param,"0EM","",0,1.5);
    Param->DrawClone("same");
    printf("30: chi2/ndf: %.2f / %d = %.2f \n", Param->GetChisquare(), Param->GetNDF(), Param->GetChisquare()/Param->GetNDF());
    Param->SetLineColor(kMagenta-3);
    gMult41Total->Fit(Param,"0EM","",0,1.5);
    Param->DrawClone("same");
    printf("41: chi2/ndf: %.2f / %d = %.2f \n", Param->GetChisquare(), Param->GetNDF(), Param->GetChisquare()/Param->GetNDF());
  }
  
//  Particle names
  TPaveText* PionName = new TPaveText(0.16,0.22,0.17,0.23);
  PionName->SetFillStyle(0);
  PionName->SetBorderSize(0);
  PionName->SetTextFont(42);
  PionName->SetTextSize(legsize);
  PionName->AddText("#pi^{#pm}");
  PionName->Draw();

  TPaveText* KaonName = new TPaveText(0.52,0.22,0.54,0.23);
  KaonName->SetFillStyle(0);
  KaonName->SetBorderSize(0);
  KaonName->SetTextFont(42);
  KaonName->SetTextSize(legsize);
  KaonName->AddText("K^{#pm}");
  KaonName->Draw();

  TPaveText* ProtonName = new TPaveText(0.88,0.22,1.03,0.23);
  ProtonName->SetFillStyle(0);
  ProtonName->SetBorderSize(0);
  ProtonName->SetTextFont(42);
  ProtonName->SetTextSize(legsize);
  ProtonName->AddText("p/#bar{p}");
  ProtonName->Draw();

  TPaveText* LambdaName = new TPaveText(1.06,0.25,1.25,0.26);
  LambdaName->SetFillStyle(0);
  LambdaName->SetBorderSize(0);
  LambdaName->SetTextFont(42);
  LambdaName->SetTextSize(legsize);
  LambdaName->AddText("#Lambda/#bar{#Lambda}");
  LambdaName->Draw();

  TPaveText* DeuteronName = new TPaveText(1.84,0.25,1.87,0.26);
  DeuteronName->SetFillStyle(0);
  DeuteronName->SetBorderSize(0);
  DeuteronName->SetTextFont(42);
  DeuteronName->SetTextSize(legsize);
  DeuteronName->AddText("d/#bar{d}");
  DeuteronName->Draw();

  TPaveText* He3Name = new TPaveText(2.55,0.25,2.9,0.26);
  He3Name->SetFillStyle(0);
  He3Name->SetBorderSize(0);
  He3Name->SetTextFont(42);
  He3Name->SetTextSize(legsize);
  He3Name->AddText("^{3}He/^{3}#bar{He}");
  He3Name->Draw();

  
  cMeanPtvsMass->SaveAs("MeanPtvsMass.eps");
}
