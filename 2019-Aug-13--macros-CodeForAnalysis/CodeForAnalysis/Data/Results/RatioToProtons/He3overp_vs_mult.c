//#include "mystyle.C" // root6

void He3overp_vs_mult(){

  bool plotXeXe = false;
  bool plotpp13TeV = false;
  bool plotTriton = true;
  bool DrawTheorie = true;

  bool Transparent = false;

  gSystem->Load("../mystyle_C.so"); //root5 only
  mystyle();
  Double_t legsize = 0.05;
  Double_t MarkerSize = 1.25;

  gStyle->SetOptStat(0);         // switch off stat box
  gStyle->SetOptFit(1);

  //=========Macro generated from canvas: ratio_cv/
  //=========  (Tue Apr 24 17:06:25 2018) by ROOT version6.10/08
  TCanvas *ratio_cv = new TCanvas("ratio_cv", "",10,758,800,600);
  ratio_cv->Range(-289.1566,-2.158531e-06,2120.482,1.582923e-05);
  gPad->SetTicks();
  gPad->SetLogx();
  gPad->SetLogy();

  TH1F *Axis = new TH1F("hframe__1","",3000,0,3000);
  Axis->GetXaxis()->SetRangeUser(4,25000);
  Axis->SetMinimum(2e-8);
  Axis->SetMaximum(3e-05);
  gPad->SetBottomMargin(0.17);
  Axis->GetXaxis()->SetLabelSize(legsize);
  Axis->GetYaxis()->SetLabelSize(legsize);
  setOPT_hists(Axis,"#it{p}_{T} (GeV/#it{c})", "Acceptance x efficiency",505,20,0.5);
  Axis->GetYaxis()->SetNdivisions(505);
  Axis->GetXaxis()->SetTitle("#LTd#it{N}_{ch}/d#it{#eta}_{lab}#GT_{|#it{#eta}_{lab}| < 0.5}");
  Axis->GetYaxis()->SetTitle("Ratio of integrated yields");
  Axis->GetXaxis()->SetTitleOffset(1.25);
  Axis->GetYaxis()->SetTitleOffset(1.1);
  Axis->Draw(" ");

  if (DrawTheorie) {
    TFile FISTfile ("He3-to-p_ThermalFIST.root");
    TCanvas* cFIST = (TCanvas*) FISTfile.GetObjectChecked("ratio_cv","TCanvas");
    TGraph* gFISThighT = (TGraph*) cFIST->GetListOfPrimitives()->FindObject("CE-HRG-highT");
    TGraph* gFIST = (TGraph*) cFIST->GetListOfPrimitives()->FindObject("CE-HRG-1");
    TGraph* gFIST2 = (TGraph*) cFIST->GetListOfPrimitives()->FindObject("CE-HRG-2");
    TGraphErrors* gFISTband = (TGraphErrors*) cFIST->GetListOfPrimitives()->FindObject("CE-HRG-band");

    TLegend* legFIST = new TLegend(0.55,0.45,0.73,0.71, "CSM (Thermal-FIST)");
    legFIST->SetBorderSize(0);
    legFIST->SetFillStyle(0);
    legFIST->SetTextSize(0.0375);

    legFIST->AddEntry(gFIST,"#it{T} = 155 MeV, #it{V}_{c} = d#it{V} / d#it{y}","l");
    legFIST->AddEntry(gFIST2,"#it{T} = 155 MeV, #it{V}_{c} = 3 d#it{V} / d#it{y}","l");

    TFile CoalFile ("He3_to_p_SunKoDoenigusCoal_published.root");
    TCanvas* cCoal = (TCanvas*) CoalFile.GetObjectChecked("horst","TCanvas");
    TGraph* TwoBody = (TGraph*) cCoal->GetListOfPrimitives()->FindObject("two_body_coal");
    TGraph* TwoBodyUp = (TGraph*) cCoal->GetListOfPrimitives()->FindObject("two_body_coal_uncert1");
    TGraph* TwoBodyDown = (TGraph*) cCoal->GetListOfPrimitives()->FindObject("two_body_coal_uncert2");
    Double_t x1, x2, x3, y1, y2;
    Double_t UncUp, UncDown;

    TGraphAsymmErrors* CoalsecenceTwoBody2 = new TGraphAsymmErrors(TwoBody->GetN());
    for (Int_t iPoint = 0; iPoint < TwoBody->GetN(); iPoint++) {
      TwoBody->GetPoint(iPoint,x1,y1);
      TwoBodyUp->GetPoint(iPoint,x2,UncUp);
      TwoBodyDown->GetPoint(iPoint,x3,UncDown);
      UncUp = UncUp-y1;
      UncDown = y1-UncDown;
      //      printf("X values: %.1f, %.1f, %.1f Uncertainties: %.2e, %.2e \n", x1, x2, x3, UncDown, UncUp);
      CoalsecenceTwoBody2->SetPoint(iPoint,x1, y1);
      CoalsecenceTwoBody2->SetPointError(iPoint,0,0, UncDown, UncUp);
    }

    TGraph* ThreeBody = (TGraph*) cCoal->GetListOfPrimitives()->FindObject("three_body_coal");
    TGraph* ThreeBodyUp = (TGraph*) cCoal->GetListOfPrimitives()->FindObject("three_body_coal_uncert1");
    TGraph* ThreeBodyDown = (TGraph*) cCoal->GetListOfPrimitives()->FindObject("three_body_coal_uncert2");

    TGraphAsymmErrors* CoalsecenceThreeBody2 = new TGraphAsymmErrors(ThreeBody->GetN());
    for (Int_t iPoint = 0; iPoint < ThreeBody->GetN(); iPoint++) {
      ThreeBody->GetPoint(iPoint,x1,y1);
      ThreeBodyUp->GetPoint(iPoint,x2,UncUp);
      ThreeBodyDown->GetPoint(iPoint,x3,UncDown);
      UncUp = UncUp-y1;
      UncDown = y1-UncDown;
      //      printf("X values: %.1f, %.1f, %.1f Uncertainties: %.2e, %.2e \n", x1, x2, x3, UncDown, UncUp);
      CoalsecenceThreeBody2->SetPoint(iPoint,x1, y1);
      CoalsecenceThreeBody2->SetPointError(iPoint,0,0, UncDown, UncUp);
    }

    if(Transparent){
      gFISTband->SetFillColorAlpha(kGray+2,0.25);
    }else{
      gFISTband->SetFillColor(kGray);
//      gFISTband->SetFillStyle(3002);
    }

    gFISTband->DrawClone("3same");
    gFIST->DrawClone("lsame");
    gFIST2->DrawClone("lsame");
    //    gFISThighT->DrawClone("lsame");

    //    TF1* CoalsecenceThreeBody = new TF1("Coalsecence","7.1e-6/pow(1+pow(1.24/(0.83*pow(x,1/3)),2),3)",0,3000); //Sun Ko Döngus
    //    TF1* CoalsecenceTwoBody = new TF1("Coalsecence","7.1e-6/(pow(1+pow(1.15/(0.83*pow(x,1/3)),2),3/2)*pow(1+pow(1.6/(0.83*pow(x,1/3)),2),3/2))",0,3000); //Sun Ko Döngus
    //    CoalsecenceThreeBody->DrawClone("same");
    //    CoalsecenceTwoBody->DrawClone("same");

    if(Transparent){
      CoalsecenceTwoBody2->SetFillColorAlpha(kViolet-2, 0.25);
    }else{
      CoalsecenceTwoBody2->SetFillColor(kMagenta-10);
//      CoalsecenceTwoBody2->SetFillStyle(3002);
    }

    CoalsecenceTwoBody2->SetLineStyle(5);
    CoalsecenceTwoBody2->SetLineWidth(2);
    CoalsecenceTwoBody2->SetLineColor(kViolet-2);
    CoalsecenceTwoBody2->DrawClone("C3same");
    if (!Transparent){
      for (int i=0; i <CoalsecenceTwoBody2->GetN(); i++) {
        CoalsecenceTwoBody2->SetPointError(i,0,0,0,0);
      }
      CoalsecenceTwoBody2->DrawClone("Csame");
    }

    if(Transparent){
      CoalsecenceThreeBody2->SetFillColorAlpha(kOrange-3, 0.25);
    }else{
      CoalsecenceThreeBody2->SetFillColor(kOrange-9);
//      CoalsecenceThreeBody2->SetFillStyle(3002);
    }

    CoalsecenceThreeBody2->SetLineStyle(7);
    CoalsecenceThreeBody2->SetLineWidth(2);
    CoalsecenceThreeBody2->SetLineColor(kOrange+7);
    CoalsecenceThreeBody2->DrawClone("C3same");
    if (!Transparent){
      for (int i=0; i <CoalsecenceThreeBody2->GetN(); i++) {
        CoalsecenceThreeBody2->SetPointError(i,0,0,0,0);
      }
      CoalsecenceThreeBody2->DrawClone("Csame");
    }

    legFIST->AddEntry((TObject*)0,"#kern[-0.3]{Coalescence}","");
    legFIST->AddEntry(CoalsecenceThreeBody2,"Three-body","lf");
    legFIST->AddEntry(CoalsecenceTwoBody2,"Two-body","lf");
    legFIST->DrawClone("same");
  }

  Double_t Graph0_fx1001[3] = {
    1764,
    826,
    131.944};
  Double_t Graph0_fy1001[3] = {
    6.649351e-06,
    8.464292e-06,
    9.353635e-06};
  Double_t Graph0_fex1001[3] = {
    0,
    0,
    0};
  Double_t Graph0_fey1001[3] = {
    2.47524e-07,
    4.500016e-07,
    4.549451e-07};
  TGraphErrors *gre = new TGraphErrors(3,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
  gre->SetName("Graph0");
  gre->SetTitle("Graph");
  gre->SetFillColor(1);
  gre->SetMarkerStyle(20);

  TH1F *Graph_Graph1001 = new TH1F("Graph_Graph1001","Graph",100,0,1927.206);
  Graph_Graph1001->SetMinimum(6.061152e-06);
  Graph_Graph1001->SetMaximum(1.014926e-05);
  Graph_Graph1001->SetDirectory(0);
  Graph_Graph1001->SetStats(0);

  Graph_Graph1001->SetLineColor(1);
  Graph_Graph1001->GetXaxis()->SetLabelFont(42);
  Graph_Graph1001->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph1001->GetXaxis()->SetTitleOffset(1.05);
  Graph_Graph1001->GetXaxis()->SetTitleFont(42);
  Graph_Graph1001->GetYaxis()->SetLabelFont(42);
  Graph_Graph1001->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph1001->GetYaxis()->SetTitleOffset(1.1);
  Graph_Graph1001->GetYaxis()->SetTitleFont(42);
  Graph_Graph1001->GetZaxis()->SetLabelFont(42);
  Graph_Graph1001->GetZaxis()->SetTitleSize(0.05);
  Graph_Graph1001->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph_Graph1001);

  //  gre->Draw("pz");

  Double_t Graph1_fx1002[3] = {
    1764,
    826,
    131.944};
  Double_t Graph1_fy1002[3] = {
    6.649351e-06,
    8.464292e-06,
    9.353635e-06};
  Double_t Graph1_fex1002[3] = {
    49.75,
    22,
    6.208};
  Double_t Graph1_fey1002[3] = {
    1.392313e-06,
    1.700516e-06,
    1.517712e-06};
  gre = new TGraphErrors(3,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
  gre->SetName("Graph1");
  gre->SetTitle("Graph");
  gre->SetFillColor(1);
  gre->SetFillStyle(0);
  gre->SetMarkerStyle(20);

  TH1F *Graph_Graph1002 = new TH1F("Graph_Graph1002","Graph",100,0,1982.551);
  Graph_Graph1002->SetMinimum(4.695607e-06);
  Graph_Graph1002->SetMaximum(1.143278e-05);
  Graph_Graph1002->SetDirectory(0);
  Graph_Graph1002->SetStats(0);

  Graph_Graph1002->SetLineColor(1);
  Graph_Graph1002->GetXaxis()->SetLabelFont(42);
  Graph_Graph1002->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph1002->GetXaxis()->SetTitleOffset(1.05);
  Graph_Graph1002->GetXaxis()->SetTitleFont(42);
  Graph_Graph1002->GetYaxis()->SetLabelFont(42);
  Graph_Graph1002->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph1002->GetYaxis()->SetTitleOffset(1.1);
  Graph_Graph1002->GetYaxis()->SetTitleFont(42);
  Graph_Graph1002->GetZaxis()->SetLabelFont(42);
  Graph_Graph1002->GetZaxis()->SetTitleSize(0.05);
  Graph_Graph1002->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph_Graph1002);

  //  gre->Draw("p2");

  Double_t Graph2_fx1003[2] = {
    1206.7,
    266};
  Double_t Graph2_fy1003[2] = {
    1.0611e-05,
    8.341e-06};
  Double_t Graph2_fex1003[2] = {
    0,
    0};
  Double_t Graph2_fey1003[2] = {
    3.59e-07,
    4.03e-07};
  gre = new TGraphErrors(2,Graph2_fx1003,Graph2_fy1003,Graph2_fex1003,Graph2_fey1003);
  gre->SetName("Graph2");
  gre->SetTitle("Graph");
  gre->SetFillColor(1);
  gre->SetMarkerSize(MarkerSize);

  gre->SetLineColor(kBlue);
  gre->SetMarkerColor(kBlue);
  gre->SetMarkerStyle(21);

  TH1F *Graph_Graph1003 = new TH1F("Graph_Graph1003","Graph",100,171.93,1300.77);
  Graph_Graph1003->SetMinimum(7.6348e-06);
  Graph_Graph1003->SetMaximum(1.12732e-05);
  Graph_Graph1003->SetDirectory(0);
  Graph_Graph1003->SetStats(0);

  Graph_Graph1003->SetLineColor(kBlue);
  Graph_Graph1003->GetXaxis()->SetLabelFont(42);
  Graph_Graph1003->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph1003->GetXaxis()->SetTitleOffset(1.05);
  Graph_Graph1003->GetXaxis()->SetTitleFont(42);
  Graph_Graph1003->GetYaxis()->SetLabelFont(42);
  Graph_Graph1003->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph1003->GetYaxis()->SetTitleOffset(1.1);
  Graph_Graph1003->GetYaxis()->SetTitleFont(42);
  Graph_Graph1003->GetZaxis()->SetLabelFont(42);
  Graph_Graph1003->GetZaxis()->SetTitleSize(0.05);
  Graph_Graph1003->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph_Graph1003);

  gre->Draw("pz");

  Double_t Graph3_fx1004[2] = {
    1206.7,
    266};
  Double_t Graph3_fy1004[2] = {
    1.0611e-05,
    8.341e-06};
  Double_t Graph3_fex1004[2] = {
    46.67,
    9.83};
  Double_t Graph3_fey1004[2] = {
    2.594953e-06,
    2.377406e-06};
  gre = new TGraphErrors(2,Graph3_fx1004,Graph3_fy1004,Graph3_fex1004,Graph3_fey1004);
  gre->SetName("Graph3");
  gre->SetTitle("Graph");
  gre->SetFillColor(kBlue);
  gre->SetFillStyle(0);
  gre->SetLineColor(kBlue);
  gre->SetMarkerColor(kBlue);
  gre->SetMarkerStyle(21);
  gre->SetMarkerSize(MarkerSize);

  TH1F *Graph_Graph1004 = new TH1F("Graph_Graph1004","Graph",100,156.45,1353.09);
  Graph_Graph1004->SetMinimum(5.239358e-06);
  Graph_Graph1004->SetMaximum(1.393019e-05);
  Graph_Graph1004->SetDirectory(0);
  Graph_Graph1004->SetStats(0);

  Graph_Graph1004->SetLineColor(kBlue);
  Graph_Graph1004->GetXaxis()->SetLabelFont(42);
  Graph_Graph1004->GetXaxis()->SetTitleSize(0.05);
  Graph_Graph1004->GetXaxis()->SetTitleOffset(1.05);
  Graph_Graph1004->GetXaxis()->SetTitleFont(42);
  Graph_Graph1004->GetYaxis()->SetLabelFont(42);
  Graph_Graph1004->GetYaxis()->SetTitleSize(0.05);
  Graph_Graph1004->GetYaxis()->SetTitleOffset(1.1);
  Graph_Graph1004->GetYaxis()->SetTitleFont(42);
  Graph_Graph1004->GetZaxis()->SetLabelFont(42);
  Graph_Graph1004->GetZaxis()->SetTitleSize(0.05);
  Graph_Graph1004->GetZaxis()->SetTitleFont(42);
  gre->SetHistogram(Graph_Graph1004);

  gre->Draw("p2");

  TFile HeliumFile ("RatioToProton_YieldMean_Average.root");
  TGraphErrors* RatioToProtonSyst = (TGraphErrors*) HeliumFile.GetObjectChecked("RatioToProtonSyst","TGraphErrors");
  RatioToProtonSyst->SetFillStyle(0);
  TGraphErrors* RatioToProton = (TGraphErrors*) HeliumFile.GetObjectChecked("RatioToProton","TGraphErrors");

  RatioToProtonSyst->SetLineColor(kRed);
  RatioToProton->SetMarkerColor(kRed);
  RatioToProton->SetMarkerSize(MarkerSize);
  RatioToProton->SetLineColor(kRed);
  RatioToProtonSyst->Draw("2same");
  RatioToProton->Draw("ZPsame");

  Double_t xValue=0;
  Double_t yValue=0;
  for (int i = 0; i < RatioToProtonSyst->GetN(); i++) {
    RatioToProton->GetPoint (i, xValue, yValue);
    RatioToProton->SetPointError(i,0,RatioToProton->GetErrorY(i));
    printf("x= %.2f: %.2e Statistical uncertainty: %.2e; systematic uncertainty: %.2e; total uncertainty %.2e \n", xValue, yValue, RatioToProton->GetErrorY(i), RatioToProtonSyst->GetErrorY(i), TMath::Sqrt(RatioToProton->GetErrorY(i)*RatioToProton->GetErrorY(i)+ RatioToProtonSyst->GetErrorY(i)*RatioToProtonSyst->GetErrorY(i)));
  }

  Double_t ppMult = 6;
  Double_t ppMultUnc = 0.2;
  Double_t ppValue = 2*0.11*1e-6/(0.123+0.124);
  Double_t ppUnc = TMath::Sqrt(TMath::Power(TMath::Sqrt(2)*0.0002/(0.124+0.123),2) + (0.6/1.1)*(0.6/1.1))*ppValue;
  TGraphErrors *ppHe = new TGraphErrors(1,&ppMult,&ppValue,0,&ppUnc);
  ppHe->SetName("ppHe");
  ppHe->SetTitle("ppHe");
  ppHe->SetFillColor(kGreen+2);
  ppHe->SetLineColor(kGreen+2);
  ppHe->SetMarkerColor(kGreen+2);
  ppHe->SetMarkerStyle(33);
  ppHe->SetMarkerSize(MarkerSize+0.5);

  ppUnc = TMath::Sqrt((0.009*0.009+0.01*0.01)/(0.124+0.123)*(0.124+0.123) + (0.2/1.1)*(0.2/1.1))*ppValue;
  TGraphErrors *ppHeSyst = new TGraphErrors(1,&ppMult,&ppValue,&ppMultUnc,&ppUnc);
  ppHeSyst->SetName("ppHeSyst");
  ppHeSyst->SetTitle("ppHeSyst");
  ppHeSyst->SetFillColor(kGreen+2);
  ppHeSyst->SetLineColor(kGreen+2);
  ppHeSyst->SetMarkerColor(kGreen+2);
  ppHeSyst->SetMarkerStyle(33);
  ppHeSyst->SetFillStyle(0);
  ppHeSyst->SetMarkerSize(MarkerSize+0.5);

  ppHeSyst->Draw("2same");
  ppHe->Draw("ZPsame");

  TLegend *leg = new TLegend(0.2,0.2,0.4,0.45,"ALICE");
  leg->SetBorderSize(0);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.0375);

  TLegendEntry* entry;
  leg->AddEntry(RatioToProton,"(^{3}He + ^{3}#bar{He}) / (p + #bar{p}), p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV (#minus1 #scale[1.2]{#leq} #it{y}_{cms} < 0)","p");
  //  entry=leg->AddEntry("Graph1","2 #scale[0.9]{#bullet} ^{3}He / (p + #bar{p}), Pb#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV, #scale[0.7]{Preliminary}","p");
  //  entry->SetLineColor(1);
  //  entry->SetLineStyle(1);
  //  entry->SetLineWidth(1);
  //  entry->SetMarkerColor(1);
  //  entry->SetMarkerStyle(20);
  //  entry->SetMarkerSize(MarkerSize);
  //  entry->SetTextFont(42);
  entry=leg->AddEntry("Graph3","2 #scale[0.9]{#bullet} ^{3}He / (p + #bar{p}), Pb#minusPb #sqrt{#it{s}_{NN}} = 2.76 TeV (|#it{y_{cms}}| < 0.5)", "p"); //, #scale[0.7]{Phys.Rev. C93 (2016) 2, 024917}","p"); // reference for presentations
  leg->AddEntry(ppHe,"2 #scale[0.9]{#bullet} ^{3}#bar{He} / (p + #bar{p}), pp #sqrt{#it{s}} = 7 TeV (|#it{y_{cms}}| < 0.5)", "p"); //, #scale[0.7]{Phys.Rev. C97 (2018) 2, 024615}","p"); // reference for presentations


  if (plotXeXe) {
    // First Xe–Xe result
    Double_t XeXeMult = 331;
    Double_t XeXeMultUnc = 3;
    Double_t XeXeValue = 14.77e-6;
    Double_t XeXeUnc = 1.31e-6;
    TGraphErrors *XeXeHe = new TGraphErrors(1,&XeXeMult,&XeXeValue,0,&XeXeUnc);
    XeXeHe->SetName("XeXeHe");
    XeXeHe->SetTitle("XeXeHe");
    XeXeHe->SetFillColor(kRed-2);
    XeXeHe->SetLineColor(kRed-2);
    XeXeHe->SetMarkerColor(kRed-2);
    XeXeHe->SetMarkerStyle(24);

    XeXeUnc = 4.46e-6; // Syst. unc.
    TGraphErrors *XeXeHeSyst = new TGraphErrors(1,&XeXeMult,&XeXeValue,&XeXeMultUnc,&XeXeUnc);
    XeXeHeSyst->SetName("XeXeHeSyst");
    XeXeHeSyst->SetTitle("XeXeHeSyst");
    XeXeHeSyst->SetFillColor(kRed-2);
    XeXeHeSyst->SetLineColor(kRed-2);
    XeXeHeSyst->SetMarkerColor(kRed-2);
    XeXeHeSyst->SetMarkerStyle(24);
    XeXeHeSyst->SetFillStyle(0);

    XeXeHeSyst->Draw("2same");
    XeXeHe->Draw("ZPsame");

    leg->AddEntry(XeXeHe,"(^{3}He + ^{3}#bar{He}) / (p + #bar{p}), Xe-Xe #sqrt{#it{s}} = 5.44 TeV (|#it{y_{cms}}| < 0.5)","p");
  }

  if (plotTriton) {
    Double_t TritonMult = 17.8;
    Double_t TritonMultUnc = 0.4;
    Double_t TritonValue = 3.31e-6;
    Double_t TritonUnc = 6.47e-07;
    TGraphErrors *Triton = new TGraphErrors(1,&TritonMult,&TritonValue,0,&TritonUnc);
    Triton->SetName("Triton");
    Triton->SetTitle("Triton");
    Triton->SetFillColor(kRed+2);
    Triton->SetLineColor(kRed+2);
    Triton->SetMarkerColor(kRed+2);
    Triton->SetMarkerStyle(24);
    Triton->SetMarkerSize(MarkerSize);

    TritonUnc = 7.25e-07; // Syst. unc.
    TGraphErrors *TritonSyst = new TGraphErrors(1,&TritonMult,&TritonValue,&TritonMultUnc,&TritonUnc);
    TritonSyst->SetName("TritonSyst");
    TritonSyst->SetTitle("TritonSyst");
    TritonSyst->SetFillColor(kRed+2);
    TritonSyst->SetLineColor(kRed+2);
    TritonSyst->SetMarkerColor(kRed+2);
    TritonSyst->SetMarkerStyle(24);
    TritonSyst->SetMarkerSize(MarkerSize);
    TritonSyst->SetFillStyle(0);

    TritonSyst->Draw("2same");
    Triton->Draw("ZPsame");

    leg->AddEntry(Triton,"(^{3}H + ^{3}#bar{H}) / (p + #bar{p}), p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV (#minus1 #scale[1.2]{#leq} #it{y}_{cms} < 0)","p");
  }

  if(plotpp13TeV){
    Double_t pp13Mult[2] = {18.685, 5.98278};
    Double_t pp13MultUnc[2] = {0.25, 0.09};
    Double_t pp13NullUnc[2] = {0,0};
    Double_t pp13Value[2] = {2.14377e-06/1.00157, 3.60773e-07/0.309989};
    Double_t pp13Unc[2] = {pp13Value[0]*sqrt((0.000452349/1.00157)**2+(1.87504e-07/2.14377e-06)**2), pp13Value[1]*sqrt((9.3491e-05/0.309989)**2+(5.11704e-08/3.60773e-07)**2)};
    TGraphErrors *pp13He = new TGraphErrors(2,pp13Mult,pp13Value,pp13NullUnc,pp13Unc);
    pp13He->SetName("pp13He");
    pp13He->SetTitle("pp13He");
    pp13He->SetFillColor(kAzure+7);
    pp13He->SetLineColor(kAzure+7);
    pp13He->SetMarkerColor(kAzure+7);
    pp13He->SetMarkerStyle(24);

    Double_t pp13SystUnc[2] = {pp13Value[0]*sqrt((0.038984/1.00157)**2+(4.23529e-07/2.14377e-06)**2), pp13Value[1]*sqrt((0.00744822/0.309989)**2+(8.60749e-08/3.60773e-07)**2)};
    TGraphErrors *pp13HeSyst = new TGraphErrors(2,pp13Mult,pp13Value,pp13MultUnc,pp13SystUnc);
    pp13HeSyst->SetName("pp13HeSyst");
    pp13HeSyst->SetTitle("pp13HeSyst");
    pp13HeSyst->SetFillColor(kAzure+7);
    pp13HeSyst->SetLineColor(kAzure+7);
    pp13HeSyst->SetMarkerColor(kAzure+7);
    pp13HeSyst->SetMarkerStyle(24);
    pp13HeSyst->SetFillStyle(0);

    pp13HeSyst->Draw("2same");
    pp13He->Draw("ZPsame");

    leg->AddEntry(pp13He,"(^{3}He + ^{3}#bar{He}) / (p + #bar{p}), pp #sqrt{#it{s}} = 13 TeV (|#it{y_{cms}}| < 0.5)","p");
  }

  leg->Draw();
  ratio_cv->Modified();
  ratio_cv->cd();
  ratio_cv->SetSelected(ratio_cv);
  ratio_cv->SaveAs("HeliumOverProton.eps");
  if(Transparent) ratio_cv->SaveAs("HeliumOverProton.pdf");

}
