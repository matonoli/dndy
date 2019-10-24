#include <TString.h>
#include <TFile.h>
#include <TList.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TKey.h>

#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>

#include "mystyle.C"

TList* GetQAListGRID(const TString str,const TString dir);
TList* GetResultListGRID(const TString str,const TString dir);

void PerformOfflineAnalysis(){
  
  Double_t HyperTritonBR = 0.25;
  Double_t HyperTritonFraction = 0.25; // evaluated by extrapolation of measurements in Pb-Pb collisions at 2.76 TeV (and 5.02 TeV)
  
  //style area
  mystyle();
  Double_t legsize = 0.05;
  gStyle->SetOptStat(0);         // switch off stat box
  Double_t ylow = -1.;
  Double_t yup = 0.;
  
  TString LegendTitle = "ALICE";
  TString AnalysisName = "#kern[-0.15]{p#minusPb #sqrt{#it{s}_{NN}} = 5.02 TeV}";
  TString RapidityInfo = "#kern[-0.25]{#minus1 < #it{y}_{cms} < 0}";
  
  Int_t nITS = 2;
  Int_t nTPC = 70; //> 80; >75; >70; >65; >60
  
  Double_t DCAxyCut = 0.1;
  Double_t DCAzCut = 1.;
  
  Double_t minTPC = -3.0;
  Double_t maxTPC =  3.0;
  Double_t TOFselect = 3.;   // triton
  
  Bool_t CalculateAphaEfficiency = false;
  Int_t nTPC4He = 70;
  Double_t TPCselect4He = 3.0;
  Double_t TOFselect4He = 5.;
  
  Bool_t saveEff = false;
  //      TString EfficiencyFileName = Form("Efficiency_DCAxy%.2fz%.2f.root",DCAxyCut,DCAzCut); // DCA
  //      TString EfficiencyFileName = Form("Efficiency_ITS%dTPC%d_%.2fy%.2f.root",nITS,nTPC,ylow,yup); // ITS TPC cluster Variation
  //      TString EfficiencyFileName = Form("Efficiency_TPC_%.1f_%.1f.root",minTPC,maxTPC); // TPC PID
  //   TString EfficiencyFileName = Form("Efficiency_TOF_%.1f.root",TOFselect); // TOF PID
  TString EfficiencyFileName = "Efficiency.root";
  
  
  Bool_t WriteToFile = false;
  Bool_t DrawDCA = false;
  Bool_t HyperTriton = true;
  
  const TString dirnameData = "He3MCTest";
  
  //  TString mcsource = "../Results/AnalysisResults_Fast_TOFlist_4He.root";
  TString mcsource = "../Results/AnalysisResults_Fast_CorrectRapidity.root";
  //   TString mcsource = "../Results/MinBias/AnalysisResults_LHC17f2ab_New.root"; // For secondary templates
  //   TString mcsource = "../Results/MinBias/AnalysisResults_LHC18f3.root"; // For secondary templates
  
  
  //   TString mcsource = "../Results/peripheral_PbPb/MaterialBudget/AnalysisResults_LHC16h7c_New.root"; // For material budget studies  Pb--Pb
  //   TString mcsource = "../Results/peripheral_PbPb/MaterialBudget/AnalysisResults_LHC17d5a_New.root"; // For material budget studies  Pb--Pb
  //      TString mcsource = "../Results/peripheral_PbPb/MaterialBudget/AnalysisResults_LHC17d5b_New.root"; // For material budget studies Pb--Pb
  
  //   TString mcsource = "../Results/peripheral_PbPb/HadronicCrossSection/AnalysisResults_LHC16h7c_G4.root"; // GEANT4 Pb--Pb
  
  TList *ResultListMC = GetResultListGRID(mcsource,dirnameData);
  
  // p_T correction
  //_______________________________________________________________________________
  TH2D* fHistPtTrueRecHe3 = (TH2D*) ResultListMC->FindObject("fHistPtTrueRecHe3");
  TCanvas* cPtCorrection = new TCanvas("cPtCorrection","Correction of the trannsverse momentum");
  gPad->SetTicks();
  gPad->SetLogz();
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.2);
  gPad->SetBottomMargin(0.15);
  fHistPtTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{T}^{ rec} (GeV/#it{c})");
  fHistPtTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{T}^{ true}- #it{p}_{T}^{ rec}) (GeV/#it{c})");
  fHistPtTrueRecHe3->GetZaxis()->SetTitle("Counts");
  fHistPtTrueRecHe3->SetTitle("");
  
  fHistPtTrueRecHe3->GetXaxis()->SetLabelSize(0.05);
  fHistPtTrueRecHe3->GetYaxis()->SetLabelSize(0.05);
  fHistPtTrueRecHe3->GetXaxis()->SetLabelOffset(0.006);
  fHistPtTrueRecHe3->GetYaxis()->SetLabelOffset(0.006);
  
  fHistPtTrueRecHe3->GetXaxis()->SetTitleSize(0.06);
  fHistPtTrueRecHe3->GetYaxis()->SetTitleSize(0.06);
  fHistPtTrueRecHe3->GetXaxis()->SetTitleOffset(1.1);
  
  fHistPtTrueRecHe3->GetZaxis()->SetLabelSize(0.05);
  fHistPtTrueRecHe3->GetZaxis()->SetLabelOffset(0.006);
  fHistPtTrueRecHe3->GetZaxis()->SetTitleSize(0.06);
  
  fHistPtTrueRecHe3->DrawCopy("COLZ");
  
  TF1* PtCorr = new TF1("PtCorr","[0]+[1]*exp([2]*TMath::Power(x,3))",0,10);
  PtCorr->SetParameter(0, 0.00068);
  PtCorr->SetParameter(1, -0.3641);
  PtCorr->SetParameter(2, -0.2386);
  PtCorr->SetLineColor(kBlack);
  
  TProfile* profile = (TProfile*) fHistPtTrueRecHe3->ProfileX();
  profile->SetMarkerStyle(20);
  profile->SetMarkerSize(0.7);
  profile->SetMarkerColor(kRed);
  profile->DrawCopy("same");
  
  TFitResultPtr FitTsalisResult = fHistPtTrueRecHe3->Fit(PtCorr, "NM", "",1.3,9.5);
  PtCorr->DrawClone("same");
  if(WriteToFile) cPtCorrection->SaveAs("PtCorrection.eps");
  
  // Helium
  THnSparseD *fPID = (THnSparseD*) ResultListMC->FindObject("fHistTPCdEdxSigmaHe3");
  THnSparseD *fPIDTrue = (THnSparseD*) ResultListMC->FindObject("fHistTrueHe3");
  // Deuteron
  THnSparseD *fPIDDeuteron = (THnSparseD*) ResultListMC->FindObject("fHistTPCdEdxSigmaDeuteron");
  // Triton
  THnSparseD *fPIDTriton = (THnSparseD*) ResultListMC->FindObject("fHistTPCdEdxSigmaTriton");
  THnSparseD *fPIDTrueTriton = (THnSparseD*) ResultListMC->FindObject("fHistTrueTriton");
  
  //      // Select multiplicity
  //      Double_t Multiplicity[2]={0,100};
  //      fPID->GetAxis(0)->SetRangeUser(Multiplicity[0] +0.01, Multiplicity[1] -0.001);
  //      fPIDTrue->GetAxis(0)->SetRangeUser(Multiplicity[0] +0.01, Multiplicity[1] -0.001);
  
  //      fPIDDeuteron->GetAxis(0)->SetRangeUser(Multiplicity[0] +0.01, Multiplicity[1] -0.001);
  //      fPIDTriton->GetAxis(0)->SetRangeUser(Multiplicity[0] +0.01, Multiplicity[1] -0.001);
  //      fPIDTrueTriton->GetAxis(0)->SetRangeUser(Multiplicity[0] +0.01, Multiplicity[1] -0.001);
  
  fPID->GetAxis(7)->SetRangeUser(ylow, yup); // select rapidity
  fPIDTrue->GetAxis(4)->SetRangeUser(ylow, yup); // select rapidity
  
  fPID->GetAxis(9)->SetRangeUser(nITS - 0.49, 7); // ITS hits
  fPID->GetAxis(10)->SetRangeUser(nTPC,170); // TPC cluster (bins filled with >= X)
  fPID->GetAxis(2)->SetRangeUser(minTPC+0.001, maxTPC-0.001); // select TPC PID cut
  
  // Triton
  fPIDTriton->GetAxis(7)->SetRangeUser(ylow, yup); // select rapidity
  fPIDTrueTriton->GetAxis(4)->SetRangeUser(ylow, yup); // select rapidity
  
  fPIDTriton->GetAxis(9)->SetRangeUser(nITS - 0.49, 7); // ITS hits
  fPIDTriton->GetAxis(10)->SetRangeUser(nTPC,170); // TPC cluster (bins filled with >= X)
  fPIDTriton->GetAxis(2)->SetRangeUser(minTPC+0.001, maxTPC-0.001); // TPC PID cut
  
  // DCA MC studies
  if (DrawDCA) {
    Int_t nPtBins = fPIDTriton->GetAxis(3)->GetNbins();
    TCanvas* cDCA[nPtBins];
    TH1D* hDCAxyPrim;
    TH1D* hDCAzPrim;
    TH1D* hDCAxyMat;
    TH1D* hDCAzMat;
    TH1D* hDCAxyDec;
    TH1D* hDCAzDec;
    
    fPIDTriton->GetAxis(2)->SetRangeUser(minTPC+0.001, maxTPC-0.001); // select TPC PID cut
    fPIDTriton->GetAxis(4)->SetRangeUser(0.01, 1.99); // select particle
    fPIDTriton->GetAxis(7)->SetRangeUser(ylow, yup); // select rapidity
    fPIDTriton->GetAxis(9)->SetRangeUser(nITS - 0.49, 7); // ITS hits
    fPIDTriton->GetAxis(10)->SetRangeUser(nTPC,170); // TPC cluster (bins filled with >= X)
    fPIDTriton->GetAxis(11)->SetRangeUser(-TOFselect+0.001, TOFselect-0.001); // TOF PID cut
    
    TFile *DCAoutputFile = new TFile ("DCATemplates_Triton.root","UPDATE");
    DCAoutputFile -> cd();
    
    Double_t BinCentralPt = 0;
    for (Int_t i =1; i <= nPtBins; i=i+2) {
      fPIDTriton->GetAxis(6)->SetRangeUser(-DCAzCut, DCAzCut); // select DCAz
      fPIDTriton->GetAxis(5)->SetRangeUser(-1., 1); // select DCAxy
      
      fPIDTriton->GetAxis(3)->SetRange(i, i+1); //pT
      BinCentralPt = (fPIDTriton->GetAxis(3)->GetBinCenter(i)+fPIDTriton->GetAxis(3)->GetBinCenter(i+1))/2.;
      fPIDTriton->GetAxis(8)->SetRangeUser(-0.49, 0.49); // select primary
      hDCAxyPrim = fPIDTriton->Projection(5);
      hDCAzPrim = fPIDTriton->Projection(6);
      hDCAxyPrim->SetName("hDCAxyPrim");
      hDCAzPrim->SetName("hDCAzPrim");
      hDCAxyPrim->SetLineColor(kRed);
      hDCAzPrim->SetLineColor(kRed);
      
      fPIDTriton->GetAxis(8)->SetRangeUser(0.51, 1.49); // select secondary - material
      hDCAxyMat = fPIDTriton->Projection(5);
      hDCAzMat = fPIDTriton->Projection(6);
      hDCAxyMat->SetName("hDCAxyMat");
      hDCAzMat->SetName("hDCAzMat");
      hDCAxyMat->SetLineColor(kBlue);
      hDCAzMat->SetLineColor(kBlue);
      
      fPIDTriton->GetAxis(8)->SetRangeUser(1.51, 2.49); // select secondary - weak decays
      hDCAxyDec = fPIDTriton->Projection(5);
      hDCAzDec = fPIDTriton->Projection(6);
      hDCAxyDec->SetName("hDCAxyDec");
      hDCAzDec->SetName("hDCAzDec");
      
      cDCA[i] = new TCanvas(Form("cDCA%d",i),Form("DCA %.2f",BinCentralPt));
      gPad->SetTicks();
      gPad->SetLogy();
      hDCAxyPrim->DrawCopy();
      hDCAxyMat->DrawCopy("same");
      hDCAxyDec->DrawCopy("same");
      
      //         if(WriteToFile){
      //            cDCA[i]->SaveAs(Form("DCA_%.2f.pdf", BinCentralPt));
      hDCAxyPrim->Write(Form("hDCAxy_Triton_Primary_%.2f", BinCentralPt));
      //            hDCAxyMat->Write(Form("hDCAxy_Triton_Material_%.2f", BinCentralPt));
      //            hDCAxyDec->Write(Form("hDCAxyDecay_%.2f", BinCentralPt));
      //         }
    }
    DCAoutputFile->Close();
    
    fPIDTriton->GetAxis(8)->UnZoom();
    fPIDTriton->GetAxis(3)->UnZoom();
  }
  
  // Efficiency
  // (anti-)3He
  fPIDTrue->GetAxis(2)->SetRangeUser(-0.49, 0.49); // select primary
  fPID->GetAxis(8)->SetRangeUser(-0.49, 0.49); // select primary
  // 3He
  fPIDTrue->GetAxis(3)->SetRangeUser(0.01, 1.99); // select particle
  TH1D* hAll = fPIDTrue->Projection(1);
  hAll->SetName("hAll");
  hAll->GetXaxis()->SetRangeUser(0, 7);
  hAll->SetMarkerColor(kBlack);
  hAll->SetLineColor(kBlack);
  
  // Anti - 3He
  fPIDTrue->GetAxis(3)->SetRangeUser(-1.99, -0.01); // select anti - particle
  TH1D* hAllAnti = fPIDTrue->Projection(1);
  hAllAnti->SetName("hAllAnti");
  hAllAnti->GetXaxis()->SetRangeUser(0, 7);
  hAllAnti->SetMarkerColor(kBlack);
  hAllAnti->SetLineColor(kBlack);
  
  // additional cuts
  fPID->GetAxis(5)->SetRangeUser(-DCAxyCut, DCAxyCut); // select DCAxy
  fPID->GetAxis(6)->SetRangeUser(-DCAzCut, DCAzCut); // select DCAz
  
  // 3He
  fPID->GetAxis(4)->SetRangeUser(0.01, 1.99); // select particle
  TH1D* hSelected = fPID->Projection(3);
  hSelected->SetName("hSelected");
  hSelected->GetXaxis()->SetRangeUser(0, 7);
  hSelected->SetMarkerColor(kBlue);
  hSelected->SetLineColor(kBlue);
  
  // Anti-3He
  fPID->GetAxis(4)->SetRangeUser(-1.99, -0.01); // select anti - particle
  TH1D* hSelectedAnti = fPID->Projection(3);
  hSelectedAnti->SetName("hSelectedAnti");
  hSelectedAnti->GetXaxis()->SetRangeUser(0, 7);
  hSelectedAnti->SetMarkerColor(kRed);
  hSelectedAnti->SetLineColor(kRed);
  
  TCanvas* cEff = new TCanvas("cEff","Efficiency of ^{3}He");
  cEff->Divide(2,1);
  cEff->cd(1);
  gPad->SetTicks();
  hAll->GetYaxis()->SetRangeUser(0, 6e4);
  hAll->DrawCopy("E");
  hSelected->DrawCopy("E same");
  
  cEff->cd(2);
  gPad->SetTicks();
  
  // Correct binning to match the one chosen for the analysis
  // Binning for the multiplicity integrated result
  const Int_t nbins=10;
  Double_t Binning[]={1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 7.};
  
  //   Binning for the comparison with the old result at 5 TeV (analysis performed by Natasha Sharma)
  //   Int_t nbins = 7;
  //   Double_t Binning[]={0.5, 1., 1.5, 2., 2.5, 3., 3.5, 5., 5.5, 6., 6.5};
  
  //   // Binning for 0-10 % multiplicity
  //      const Int_t nbins = 6;
  //      Double_t Binning[]={1., 1.5, 2., 2.5, 3., 3.5, 5.};
  
  //   // Binning for 10-20 % and 20-40 % multiplicity
  //      const Int_t nbins = 5;
  //      Double_t Binning[]={1., 1.5, 2., 2.5, 3., 5.};
  
  
  //   // Binning for 40-100 % multiplicity
  //      const Int_t nbins = 4;
  //      Double_t Binning[]={1., 1.5, 2., 3., 5.};
  
  hSelected = (TH1D*) hSelected->Rebin(nbins, hSelected->GetName() , Binning);
  hAll = (TH1D*) hAll->Rebin(nbins, hAll->GetName() , Binning);
  
  TH1D* hEffHelium = (TH1D*)hSelected->Clone("hEffHelium");
  hEffHelium->Divide(hSelected,hAll,1,1,"B");
  hEffHelium->GetYaxis()->SetRangeUser(0, 1);
  hEffHelium->DrawCopy("E");
  
  if(WriteToFile) cEff->SaveAs("Efficiency3He.pdf");
  
  TCanvas* cEffAnti = new TCanvas("cEffAnti","Efficiency of Anti-^{3}He");
  cEffAnti->Divide(2,1);
  cEffAnti->cd(1);
  gPad->SetTicks();
  hAllAnti->GetYaxis()->SetRangeUser(0, 6e4);
  hAllAnti->DrawCopy("E");
  hSelectedAnti->DrawCopy("E same");
  cEffAnti->cd(2);
  gPad->SetTicks();
  
  // Rebin also selected particles spectra
  hSelectedAnti = (TH1D*) hSelectedAnti->Rebin(nbins, hSelectedAnti->GetName() , Binning);
  hAllAnti = (TH1D*) hAllAnti->Rebin(nbins, hAllAnti->GetName() , Binning);
  
  TH1D* hEffAntiHelium = (TH1D*)hSelectedAnti->Clone("hEffAntiHelium");
  hEffAntiHelium->Divide(hSelectedAnti,hAllAnti,1,1,"B");
  hEffAntiHelium->GetYaxis()->SetRangeUser(0, 1);
  hEffAntiHelium->DrawCopy("E");
  if(WriteToFile) cEffAnti->SaveAs("EfficiencyAnti3He.pdf");
  
  if(HyperTriton){
    TFile *ContaminationOutputFile = new TFile ("HypertritonFeedDown.root","RECREATE");
    ContaminationOutputFile -> cd();
    
    // 3He from Hypertriton efficicency
    // Hyper-Triton
    //      THnSparseD *fPIDHyperTriton = (THnSparseD*) ResultListMC->FindObject("fHistTrueHyperTriton");
    //      fPIDHyperTriton->GetAxis(4)->SetRangeUser(ylow, yup); // select rapidity
    //      fPIDHyperTriton->GetAxis(2)->SetRangeUser(-0.49, 0.49); // select primary
    //      fPIDHyperTriton->GetAxis(3)->SetRangeUser(0.01, 1.99); // select particle
    
    fPIDTrue->GetAxis(2)->SetRangeUser(1.51, 2.49); // select secondary weak decay
    fPIDTrue->GetAxis(3)->SetRangeUser(0.01, 1.99); // select particle
    TH1D* hAllHyperTriton = fPIDTrue->Projection(1);
    hAllHyperTriton->SetName("hAllHyperTriton");
    hAllHyperTriton->GetXaxis()->SetRangeUser(0, 7);
    hAllHyperTriton->SetMarkerColor(kBlack);
    hAllHyperTriton->SetLineColor(kBlack);
    
    fPIDTrue->GetAxis(3)->SetRangeUser(-1.99, -0.01); // select anti - particle
    TH1D* hAllAntiHyperTriton = fPIDTrue->Projection(1);
    hAllAntiHyperTriton->SetName("hAllAnti");
    hAllAntiHyperTriton->GetXaxis()->SetRangeUser(0, 7);
    hAllAntiHyperTriton->SetMarkerColor(kBlack);
    hAllAntiHyperTriton->SetLineColor(kBlack);
    
    fPID->GetAxis(8)->SetRangeUser(1.51, 2.49); // select weak decays
    
    fPID->GetAxis(4)->SetRangeUser(0.01, 1.99); // select particle
    TH1D* hSelectedHyperTriton = fPID->Projection(3);
    hSelectedHyperTriton->SetName("hSelectedHyperTriton");
    hSelectedHyperTriton->GetXaxis()->SetRangeUser(0, 7);
    hSelectedHyperTriton->SetMarkerColor(kBlue);
    hSelectedHyperTriton->SetLineColor(kBlue);
    
    fPID->GetAxis(4)->SetRangeUser(-1.99, -0.01); // select anti - particle
    TH1D* hSelectedAntiHyperTriton = fPID->Projection(3);
    hSelectedAntiHyperTriton->SetName("hSelectedAntiHyperTriton");
    hSelectedAntiHyperTriton->GetXaxis()->SetRangeUser(0, 7);
    hSelectedAntiHyperTriton->SetMarkerColor(kRed);
    hSelectedAntiHyperTriton->SetLineColor(kRed);
    
    TCanvas* cEffAntiHyperTriton = new TCanvas("cEffAntiHyperTriton","Efficiency of Anti-^{3}He from Anti-Hypertriton");
    cEffAntiHyperTriton->Divide(2,1);
    cEffAntiHyperTriton->cd(1);
    gPad->SetTicks();
    hAllAntiHyperTriton->GetYaxis()->SetRangeUser(0, 6e4);
    hAllAntiHyperTriton->DrawCopy("E");
    hSelectedAntiHyperTriton->DrawCopy("E same");
    cEffAntiHyperTriton->cd(2);
    gPad->SetTicks();
    
    // Rebin also selected particles spectra
    hSelectedAntiHyperTriton = (TH1D*) hSelectedAntiHyperTriton->Rebin(nbins, hSelectedAntiHyperTriton->GetName() , Binning);
    hAllAntiHyperTriton = (TH1D*) hAllAntiHyperTriton->Rebin(nbins, hAllAntiHyperTriton->GetName() , Binning);
    
    TH1D* hEffAntiHyperTriton = (TH1D*)hSelectedAntiHyperTriton->Clone("hEffAntiHyperTriton");
    hEffAntiHyperTriton->Divide(hSelectedAntiHyperTriton,hAllAntiHyperTriton,1,1,"B");
    hEffAntiHyperTriton->GetYaxis()->SetRangeUser(0, 1);
    hEffAntiHyperTriton->DrawCopy("E");
    cEffAntiHyperTriton->SaveAs("EfficiencyAntiHyperTriton.pdf");
    
    TCanvas* cEffHyperTriton = new TCanvas("cEffHyperTriton","Efficiency of ^{3}He from Hypertriton");
    cEffHyperTriton->Divide(2,1);
    cEffHyperTriton->cd(1);
    gPad->SetTicks();
    hAllHyperTriton->GetYaxis()->SetRangeUser(0, 6e4);
    hAllHyperTriton->DrawCopy("E");
    hSelectedHyperTriton->DrawCopy("E same");
    cEffHyperTriton->cd(2);
    gPad->SetTicks();
    
    // Rebin also selected particles spectra
    hSelectedHyperTriton = (TH1D*) hSelectedHyperTriton->Rebin(nbins, hSelectedHyperTriton->GetName() , Binning);
    hAllHyperTriton = (TH1D*) hAllHyperTriton->Rebin(nbins, hAllHyperTriton->GetName() , Binning);
    
    TH1D* hEffHyperTriton = (TH1D*)hSelectedHyperTriton->Clone("hEffHyperTriton");
    hEffHyperTriton->Divide(hSelectedHyperTriton,hAllHyperTriton,1,1,"B");
    hEffHyperTriton->GetYaxis()->SetRangeUser(0, 1);
    hEffHyperTriton->DrawCopy("E");
    cEffHyperTriton->SaveAs("EfficiencyHyperTriton.pdf");
    
    TH1I* Axis = new TH1I("Axis","",55,0,5.5);
    Axis->GetXaxis()->SetRangeUser(1.,5.5);
    Axis->GetYaxis()->SetRangeUser(0.55, 0.65);
    setOPT_hists(Axis,"#it{p}_{T} (GeV/#it{c})", "#varepsilon_{feed-down} / #varepsilon_{^{3}He}",505,20,0.5);
    Axis->GetXaxis()->SetLabelSize(legsize);
    Axis->GetYaxis()->SetLabelSize(legsize);
    Axis->GetXaxis()->SetTitleSize(legsize+0.01);
    Axis->GetYaxis()->SetTitleSize(legsize+0.01);
    Axis->GetYaxis()->SetTitleOffset(1.2);
    Axis->GetYaxis()->SetNdivisions(502);
    
    TCanvas* cEffRatio = new TCanvas("cEffRatio","Secondary over primary eff",800,600);
    TH1D* hSecOverPrim = (TH1D*)hEffHyperTriton->Clone("hSecOverPrim");
    hSecOverPrim->Divide(hEffHyperTriton,hEffHelium,1,1);
    hSecOverPrim->GetXaxis()->SetRangeUser(1.5,5.);
    Axis->DrawCopy();
    setOPT_hists(hSecOverPrim,"#it{p}_{T} (GeV/#it{c})", "#varepsilon_{feed-down}/#varepsilon_{^{3}He}",505,20,1,kBlue);
    hSecOverPrim->DrawCopy("EX0same");
    
    TH1D* hSecOverPrimAnti = (TH1D*)hEffAntiHyperTriton->Clone("hSecOverPrimAnti");
    hSecOverPrimAnti->Divide(hEffAntiHyperTriton,hEffAntiHelium,1,1);
    hSecOverPrimAnti->GetXaxis()->SetRangeUser(1.5,5.);
    setOPT_hists(hSecOverPrimAnti,"#it{p}_{T} (GeV/#it{c})", "#varepsilon_{feed-down}/#varepsilon_{^{3}He}",505,20,1,kRed);
    hSecOverPrimAnti->DrawCopy("EX0same");
  
    TLegend* leg = plotLegend("left_bottom","", 1.3, 0.8,0.02,0.03);
    leg->SetTextSize(legsize);
    leg->AddEntry((TObject*)0,AnalysisName,"");
    leg->AddEntry((TObject*)0,RapidityInfo,"");
    leg->AddEntry(hSecOverPrim, "^{3}He","p");
    leg->AddEntry(hSecOverPrimAnti, "^{3}#bar{He}","p");
    leg->Draw();
    gPad->Update();
    
    cEffRatio->SaveAs("FeedDownEfficiencyRatio.pdf");
    
    //    calculated feed-down
    hSecOverPrim->Scale(HyperTritonBR*HyperTritonFraction);
    hSecOverPrimAnti->Scale(HyperTritonBR*HyperTritonFraction);
    
    hSecOverPrim->Write("hFeedDownFractionHelium");
    hSecOverPrimAnti->Write("hFeedDownFractionAntiHelium");
    ContaminationOutputFile->Close();
  }
  
  // (anti-)triton
  //MC truth
  fPIDTrueTriton->GetAxis(2)->SetRangeUser(-0.49, 0.49); // select primary
  fPIDTriton->GetAxis(8)->SetRangeUser(-0.49, 0.49); // select primary
  // triton
  fPIDTrueTriton->GetAxis(3)->SetRangeUser(0.01, 1.99); // select particle
  TH1D* hAllTriton = fPIDTrueTriton->Projection(1);
  hAllTriton->SetName("hAllTriton");
  hAllTriton->GetXaxis()->SetRangeUser(0, 7);
  hAllTriton->SetMarkerColor(kBlack);
  hAllTriton->SetLineColor(kBlack);
  
  // Anti - triton
  fPIDTrueTriton->GetAxis(3)->SetRangeUser(-1.99, -0.01); // select anti - particle
  TH1D* hAllAntiTriton = fPIDTrueTriton->Projection(1);
  hAllAntiTriton->SetName("hAllAntiTriton");
  hAllAntiTriton->GetXaxis()->SetRangeUser(0, 7);
  hAllAntiTriton->SetMarkerColor(kBlack);
  hAllAntiTriton->SetLineColor(kBlack);
  
  // MC after track selection
  // additional cuts
  fPIDTriton->GetAxis(5)->SetRangeUser(-DCAxyCut, DCAxyCut); // select DCAxy
  fPIDTriton->GetAxis(6)->SetRangeUser(-DCAzCut, DCAzCut); // select DCAz
  
  // triton
  fPIDTriton->GetAxis(4)->SetRangeUser(0.01, 1.99); // select particle
  TH1D* hSelectedTriton = fPIDTriton->Projection(3);
  hSelectedTriton->SetName("hSelectedTriton");
  hSelectedTriton->GetXaxis()->SetRangeUser(0, 7);
  hSelectedTriton->SetMarkerColor(kBlue);
  hSelectedTriton->SetLineColor(kBlue);
  
  // anti-triton
  fPIDTriton->GetAxis(4)->SetRangeUser(-1.99, -0.01); // select anti - particle
  TH1D* hSelectedAntiTriton = fPIDTriton->Projection(3);
  hSelectedAntiTriton->SetName("hSelectedAntiTriton");
  hSelectedAntiTriton->GetXaxis()->SetRangeUser(0, 7);
  hSelectedAntiTriton->SetMarkerColor(kRed);
  hSelectedAntiTriton->SetLineColor(kRed);
  
  TCanvas* cEffTriton = new TCanvas("cEffTriton","Efficiency of triton");
  cEffTriton->Divide(2,1);
  cEffTriton->cd(1);
  gPad->SetTicks();
  hAllTriton->GetYaxis()->SetRangeUser(0, 6e4);
  hAllTriton->DrawCopy("E");
  hSelectedTriton->DrawCopy("E same");
  
  cEffTriton->cd(2);
  gPad->SetTicks();
  
  // Correct binning to match the one chosen for the analysis
  //TPC only
  const Int_t nPtBins=5;
  Double_t BinningTriton[]={0.5, 0.75, 1., 1.25, 1.5, 2.};
  
  hSelectedTriton = (TH1D*) hSelectedTriton->Rebin(nPtBins, hSelectedTriton->GetName() , BinningTriton);
  hAllTriton = (TH1D*) hAllTriton->Rebin(nPtBins, hAllTriton->GetName() , BinningTriton);
  
  TH1D* hEffTriton = (TH1D*)hSelectedTriton->Clone("hEffTriton");
  hEffTriton->Divide(hSelectedTriton,hAllTriton,1,1,"B");
  hEffTriton->GetYaxis()->SetRangeUser(0, 1);
  hEffTriton->DrawCopy("E");
  
  if(WriteToFile) cEffTriton->SaveAs("EfficiencyTriton.pdf");
  
  TCanvas* cEffAntiTriton = new TCanvas("cEffAntiTriton","Efficiency of Anti-triton");
  cEffAntiTriton->Divide(2,1);
  cEffAntiTriton->cd(1);
  gPad->SetTicks();
  hAllAntiTriton->GetYaxis()->SetRangeUser(0, 6e4);
  hAllAntiTriton->DrawCopy("E");
  hSelectedAntiTriton->DrawCopy("E same");
  cEffAntiTriton->cd(2);
  gPad->SetTicks();
  
  // Rebin also selected particles spectra
  hSelectedAntiTriton = (TH1D*) hSelectedAntiTriton->Rebin(nPtBins, hSelectedAntiTriton->GetName() , BinningTriton);
  hAllAntiTriton = (TH1D*) hAllAntiTriton->Rebin(nPtBins, hAllAntiTriton->GetName() , BinningTriton);
  
  TH1D* hEffAntiTriton = (TH1D*)hSelectedAntiTriton->Clone("hEffAntiTriton");
  hEffAntiTriton->Divide(hSelectedAntiTriton,hAllAntiTriton,1,1,"B");
  hEffAntiTriton->GetYaxis()->SetRangeUser(0, 1);
  hEffAntiTriton->DrawCopy("E");
  if(WriteToFile) cEffAntiTriton->SaveAs("EfficiencyAntiTriton.pdf");
  
  fPIDTriton->GetAxis(11)->SetRangeUser(-TOFselect+0.001, TOFselect-0.001); // TOF PID cut
  // (anti-)triton
  //MC truth
  fPIDTrueTriton->GetAxis(2)->SetRangeUser(-0.49, 0.49); // select primary
  fPIDTriton->GetAxis(8)->SetRangeUser(-0.49, 0.49); // select primary
  // triton
  fPIDTrueTriton->GetAxis(3)->SetRangeUser(0.01, 1.99); // select particle
  TH1D* hAllTritonTOF = fPIDTrueTriton->Projection(1);
  hAllTritonTOF->SetName("hAllTritonTOF");
  hAllTritonTOF->GetXaxis()->SetRangeUser(0, 7);
  hAllTritonTOF->SetMarkerColor(kBlack);
  hAllTritonTOF->SetLineColor(kBlack);
  
  // Anti - triton
  fPIDTrueTriton->GetAxis(3)->SetRangeUser(-1.99, -0.01); // select anti - particle
  TH1D* hAllAntiTritonTOF = fPIDTrueTriton->Projection(1);
  hAllAntiTritonTOF->SetName("hAllAntiTritonTOF");
  hAllAntiTritonTOF->GetXaxis()->SetRangeUser(0, 7);
  hAllAntiTritonTOF->SetMarkerColor(kBlack);
  hAllAntiTritonTOF->SetLineColor(kBlack);
  
  // MC after track selection
  // additional cuts
  fPIDTriton->GetAxis(5)->SetRangeUser(-DCAxyCut, DCAxyCut); // select DCAxy
  fPIDTriton->GetAxis(6)->SetRangeUser(-DCAzCut, DCAzCut); // select DCAz
  
  // triton
  fPIDTriton->GetAxis(4)->SetRangeUser(0.01, 1.99); // select particle
  TH1D* hSelectedTritonTOF = fPIDTriton->Projection(3);
  hSelectedTritonTOF->SetName("hSelectedTritonTOF");
  hSelectedTritonTOF->GetXaxis()->SetRangeUser(0, 7);
  hSelectedTritonTOF->SetMarkerColor(kBlue);
  hSelectedTritonTOF->SetLineColor(kBlue);
  
  // anti-triton
  fPIDTriton->GetAxis(4)->SetRangeUser(-1.99, -0.01); // select anti - particle
  TH1D* hSelectedAntiTritonTOF = fPIDTriton->Projection(3);
  hSelectedAntiTritonTOF->SetName("hSelectedAntiTritonTOF");
  hSelectedAntiTritonTOF->GetXaxis()->SetRangeUser(0, 7);
  hSelectedAntiTritonTOF->SetMarkerColor(kRed);
  hSelectedAntiTritonTOF->SetLineColor(kRed);
  
  TCanvas* cEffTritonTOF = new TCanvas("cEffTritonTOF","Efficiency of triton (TOF-TPC)");
  cEffTritonTOF->Divide(2,1);
  cEffTritonTOF->cd(1);
  gPad->SetTicks();
  hAllTritonTOF->GetYaxis()->SetRangeUser(0, 6e4);
  hAllTritonTOF->DrawCopy("E");
  hSelectedTritonTOF->DrawCopy("E same");
  
  cEffTritonTOF->cd(2);
  gPad->SetTicks();
  
  // Correct binning to match the one chosen for the analysis
  // Binning for TPC-TOF
  const Int_t nPtBinsTOF=8.;
  Double_t BinningTritonTOF[]={1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3.0};
  
  hSelectedTritonTOF = (TH1D*) hSelectedTritonTOF->Rebin(nPtBinsTOF, hSelectedTritonTOF->GetName() , BinningTritonTOF);
  hAllTritonTOF = (TH1D*) hAllTritonTOF->Rebin(nPtBinsTOF, hAllTritonTOF->GetName() , BinningTritonTOF);
  
  TH1D* hEffTritonTOF = (TH1D*)hSelectedTritonTOF->Clone("hEffTritonTOF");
  hEffTritonTOF->Divide(hSelectedTritonTOF,hAllTritonTOF,1,1,"B");
  hEffTritonTOF->GetYaxis()->SetRangeUser(0, 1);
  hEffTritonTOF->DrawCopy("E");
  
  if(WriteToFile) cEffTritonTOF->SaveAs("EfficiencyTritonTOF.pdf");
  
  TCanvas* cEffAntiTritonTOF = new TCanvas("cEffAntiTritonTOF","Efficiency of Anti-triton (TOF-TPC)");
  cEffAntiTritonTOF->Divide(2,1);
  cEffAntiTritonTOF->cd(1);
  gPad->SetTicks();
  hAllAntiTritonTOF->GetYaxis()->SetRangeUser(0, 6e4);
  hAllAntiTritonTOF->DrawCopy("E");
  hSelectedAntiTritonTOF->DrawCopy("E same");
  cEffAntiTritonTOF->cd(2);
  gPad->SetTicks();
  
  // Rebin also selected particles spectra
  hSelectedAntiTritonTOF = (TH1D*) hSelectedAntiTritonTOF->Rebin(nPtBinsTOF, hSelectedAntiTritonTOF->GetName() , BinningTritonTOF);
  hAllAntiTritonTOF = (TH1D*) hAllAntiTritonTOF->Rebin(nPtBinsTOF, hAllAntiTritonTOF->GetName() , BinningTritonTOF);
  
  TH1D* hEffAntiTritonTOF = (TH1D*)hSelectedAntiTritonTOF->Clone("hEffAntiTritonTOF");
  hEffAntiTritonTOF->Divide(hSelectedAntiTritonTOF,hAllAntiTritonTOF,1,1,"B");
  hEffAntiTritonTOF->GetYaxis()->SetRangeUser(0, 1);
  hEffAntiTritonTOF->DrawCopy("E");
  if(WriteToFile) cEffAntiTritonTOF->SaveAs("EfficiencyAntiTritonTOF.pdf");
  
  if (saveEff) {
    TFile *outputFile = new TFile (EfficiencyFileName.Data(),"RECREATE");
    outputFile -> cd();
    hEffHelium -> Write();
    hEffAntiHelium -> Write();
    hEffTriton -> Write();
    hEffAntiTriton -> Write();
    hEffTritonTOF -> Write();
    hEffAntiTritonTOF -> Write();
    outputFile->Close();
  }
  
  TCanvas* cCombareEff = new TCanvas("cCombareEff","Efficiency of (anti-)^{3}He");
  gPad->SetTicks();
  hEffHelium->DrawCopy("E");
  hEffAntiHelium->DrawCopy("Esame");
  if(WriteToFile) cCombareEff->SaveAs("Efficiencies.pdf");
  
  TCanvas* cCombareEffTriton = new TCanvas("cCombareEffTriton","Efficiency of (anti-)triton");
  gPad->SetTicks();
  hEffTriton->DrawCopy("E");
  hEffAntiTriton->DrawCopy("Esame");
  if(WriteToFile) cCombareEffTriton->SaveAs("EfficienciesTriton.pdf");
  
  TCanvas* cCombareEffTritonTOF = new TCanvas("cCombareEffTritonTOF","Efficiency of (anti-)triton (TOF-TPC)");
  gPad->SetTicks();
  hEffTritonTOF->DrawCopy("E");
  hEffAntiTritonTOF->DrawCopy("Esame");
  if(WriteToFile) cCombareEffTritonTOF->SaveAs("EfficienciesTritonTOF.pdf");
  
  if (CalculateAphaEfficiency) {
    
    // p_T correction
    //_______________________________________________________________________________
    TH2D* fHistPtTrueRecHe4 = (TH2D*) ResultListMC->FindObject("fHistPtTrueRecHe4");
    TCanvas* cPtCorrection2 = new TCanvas("cPtCorrection2","Correction of the trannsverse momentum (^{4}He)");
    gPad->SetTicks();
    gPad->SetLogz();
    fHistPtTrueRecHe4->DrawCopy("COLZ");
    
    TF1* PtCorr2 = new TF1("PtCorr2","[0]+[1]*exp([2]*x)",0,10);
    PtCorr2->SetParameter(0, 0.00068);
    PtCorr2->SetParameter(1, -0.3641);
    PtCorr2->SetParameter(2, -0.2386);
    PtCorr2->SetLineColor(kBlack);
    
    TProfile* profile2 = (TProfile*) fHistPtTrueRecHe4->ProfileX();
    profile2->SetMarkerStyle(20);
    profile2->SetMarkerSize(0.7);
    profile2->SetMarkerColor(kRed);
    profile2->DrawCopy("same");
    
    TFitResultPtr FitTsalisResult = fHistPtTrueRecHe4->Fit(PtCorr2, "NM", "",1.3,9.5);
    PtCorr2->DrawClone("same");
    if(WriteToFile) cPtCorrection2->SaveAs("PtCorrection_Alpha.pdf");
    
    // Triton
    THnSparseD *fPID4He = (THnSparseD*) ResultListMC->FindObject("fHistTPCdEdxSigmaHe4");
    THnSparseD *fPIDTrue4He = (THnSparseD*) ResultListMC->FindObject("fHistTrueHe4");
    
    // Triton
    fPID4He->GetAxis(7)->SetRangeUser(ylow, yup); // select rapidity
    fPIDTrue4He->GetAxis(4)->SetRangeUser(ylow, yup); // select rapidity
    
    fPID4He->GetAxis(9)->SetRangeUser(nITS - 0.49, 7); // ITS hits
    fPID4He->GetAxis(10)->SetRangeUser(nTPC4He ,170); // TPC cluster (bins filled with >= X)
    fPID4He->GetAxis(2)->SetRangeUser(-(TPCselect4He-0.001), TPCselect4He-0.001); // TPC PID cut
    fPID4He->GetAxis(11)->SetRangeUser(-TOFselect4He+0.001, TOFselect4He-0.001); // TOF PID cut
    
    // (anti-)triton
    //MC truth
    fPIDTrue4He->GetAxis(2)->SetRangeUser(-0.49, 0.49); // select primary
    fPID4He->GetAxis(8)->SetRangeUser(-0.49, 0.49); // select primary
    // triton
    fPIDTrue4He->GetAxis(3)->SetRangeUser(0.01, 1.99); // select particle
    TH1D* hAll4He = fPIDTrue4He->Projection(1);
    hAll4He->SetName("hAll4He");
    hAll4He->GetXaxis()->SetRangeUser(0, 10);
    hAll4He->SetMarkerColor(kBlack);
    hAll4He->SetLineColor(kBlack);
    
    // Anti - triton
    fPIDTrue4He->GetAxis(3)->SetRangeUser(-1.99, -0.01); // select anti - particle
    TH1D* hAllAnti4He = fPIDTrue4He->Projection(1);
    hAllAnti4He->SetName("hAllAnti4He");
    hAllAnti4He->GetXaxis()->SetRangeUser(0, 10);
    hAllAnti4He->SetMarkerColor(kOrange-3);
    hAllAnti4He->SetLineColor(kOrange-3);
    
    // triton
    fPID4He->GetAxis(4)->SetRangeUser(0.01, 1.99); // select particle
    TH1D* hSelected4He = fPID4He->Projection(3);
    hSelected4He->SetName("hSelected4He");
    hSelected4He->GetXaxis()->SetRangeUser(0, 10);
    hSelected4He->SetMarkerColor(kBlue);
    hSelected4He->SetLineColor(kBlue);
    
    // anti-triton
    fPID4He->GetAxis(4)->SetRangeUser(-1.99, -0.01); // select anti - particle
    TH1D* hSelectedAnti4He = fPID4He->Projection(3);
    hSelectedAnti4He->SetName("hSelectedAnti4He");
    hSelectedAnti4He->GetXaxis()->SetRangeUser(0, 10);
    hSelectedAnti4He->SetMarkerColor(kRed);
    hSelectedAnti4He->SetLineColor(kRed);
    
    TCanvas* cEff4He = new TCanvas("cEff4He","Efficiency of ^{4}He (TOF-TPC)");
    cEff4He->Divide(2,1);
    cEff4He->cd(1);
    gPad->SetTicks();
    hAll4He->GetYaxis()->SetRangeUser(0, 8e3);
    hAll4He->DrawCopy("E");
    hAllAnti4He->DrawCopy("Esame");
    hSelected4He->DrawCopy("E same");
    hSelectedAnti4He->DrawCopy("E same");
    
    cEff4He->cd(2);
    gPad->SetTicks();
    
    //    const Int_t nPtBins4He = 1;
    //    Double_t BinningTriton4He[]={2.,10.0};
    
    //    hSelected4He = (TH1D*) hSelected4He->Rebin(nPtBins4He, hSelected4He->GetName() , BinningTriton4He);
    //    hAll4He = (TH1D*) hAll4He->Rebin(nPtBins4He, hAll4He->GetName() , BinningTriton4He);
    
    TH1D* hEff4He = (TH1D*)hSelected4He->Clone("hEff4He");
    hEff4He->Divide(hSelected4He,hAll4He,1,1,"B");
    hEff4He->GetYaxis()->SetRangeUser(0, 1);
    hEff4He->DrawCopy("E");
    
    //    hSelectedAnti4He = (TH1D*) hSelectedAnti4He->Rebin(nPtBins4He, hSelectedAnti4He->GetName() , BinningTriton4He);
    //    hAllAnti4He = (TH1D*) hAllAnti4He->Rebin(nPtBins4He, hAllAnti4He->GetName() , BinningTriton4He);
    
    TH1D* hEffAnti4He = (TH1D*)hSelectedAnti4He->Clone("hEffAnti4He");
    hEffAnti4He->Divide(hSelectedAnti4He,hAllAnti4He,1,1,"B");
    hEffAnti4He->GetYaxis()->SetRangeUser(0, 1);
    hEffAnti4He->DrawCopy("Esame");
    if(WriteToFile) cEff4He->SaveAs("Efficiency_4He.pdf");
    
    TFile *outputFile2 = new TFile ("Efficiency_Alpha.root","RECREATE");
    outputFile2 -> cd();
    hEff4He->Write("hEff4He");
    hEffAnti4He->Write("hEffAnti4He");
    outputFile2->Close();
    
  }
  
}

//_________________________________________________________________________
TList* GetResultListGRID(const TString str,const TString dir){
  
  TFile *f = TFile::Open(str.Data());
  if(!f || f->IsZombie()){
    printf("Could not read file %s\n",str.Data());
    return NULL ;
  }
  if(f->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",str.Data());
  }
  TDirectory *d = f->GetDirectory(dir.Data());
  if(!d || d->IsZombie()){
    printf("Could not read file %s\n",dir.Data());
    return NULL ;
  }
  if(d->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",dir.Data());
  }
  
  TKey *k;
  TIter next(d->GetListOfKeys());
  while ((k = dynamic_cast<TKey *>(next()))){TString s(k->GetName()); if(s.Contains("Results")) break;}
  if(!k){
    printf("Output container not found\n");
    f->Close(); delete f;
    return NULL;
  }
  TList *returnlist = dynamic_cast<TList *>(k->ReadObj());
  f->Close(); delete f;
  return returnlist;
}

//_________________________________________________________________________
TList* GetQAListGRID(const TString str,const TString dir){
  
  TFile *f = TFile::Open(str.Data());
  if(!f || f->IsZombie()){
    printf("Could not read file %s\n",str.Data());
    return NULL ;
  }
  if(f->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",str.Data());
  }
  TDirectory *d = f->GetDirectory(dir.Data());
  if(!d || d->IsZombie()){
    printf("Could not read file %s\n",dir.Data());
    return NULL ;
  }
  if(d->TestBit(TFile::kRecovered)){
    printf("File \"%s\" is corrupt!\n",dir.Data());
  }
  
  TKey *k;
  TIter next(d->GetListOfKeys());
  while ((k = dynamic_cast<TKey *>(next()))){
    TString s(k->GetName());
    if(s.Contains("QA")) break;
  }
  if(!k){
    printf("Output container not found\n");
    f->Close(); delete f;
    return NULL;
  }
  TList *returnlist = dynamic_cast<TList *>(k->ReadObj());
  f->Close(); delete f;
  return returnlist;
}
