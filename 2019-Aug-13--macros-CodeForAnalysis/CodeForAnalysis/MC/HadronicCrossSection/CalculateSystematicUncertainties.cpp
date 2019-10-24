#include <stdio.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

using namespace std;

void CalculateSystematicUncertainties(){

   gROOT->SetStyle("Plain");
   //style area
   gSystem->Load("/Users/tidus/Documents/Studium/PhD/Nuclei/GITRepository/AliAnalysisHe/reducedTreeNuclei/OfflineFramework/Results/mystyle_C.so");
   mystyle();

   TString pTLabel= "#it{p}_{T} (GeV/#it{c})";
   TString yLable = "Acceptance x Efficiency";

   TFile DefaultFile ("Efficiency_GEANT3_TPC70cl.root","READ");
   if (!DefaultFile || DefaultFile.IsZombie()) {
      printf("Cannot open GEANT3 file \n");
      return;
   }
   TH1D* hDefaultTriton = (TH1D*) DefaultFile.GetObjectChecked("hEffTriton","TH1D");
   TH1D* hDefaultAntiTriton = (TH1D*) DefaultFile.GetObjectChecked("hEffAntiTriton","TH1D");

   const Int_t nPtBins=3;
   Double_t Binning[]={0.5, 1., 1.5, 2.}; // triton (TPConly)

   // rebin efficiency
   for (Int_t i=1; i<=hDefaultTriton->GetNbinsX(); i++) {
      hDefaultTriton->SetBinContent(i, hDefaultTriton->GetBinContent(i)*hDefaultTriton->GetBinWidth(i));
      hDefaultAntiTriton->SetBinContent(i, hDefaultAntiTriton->GetBinContent(i)*hDefaultAntiTriton->GetBinWidth(i));
   }
   hDefaultTriton = (TH1D*) hDefaultTriton->Rebin(nPtBins,"hEffTriton",Binning);
   hDefaultAntiTriton = (TH1D*) hDefaultAntiTriton->Rebin(nPtBins,"hEffAntiTriton",Binning);

   hDefaultTriton->Scale(1,"width");
   hDefaultAntiTriton->Scale(1,"width");

   TH1D* hDefaultTritonTOF = (TH1D*) DefaultFile.GetObjectChecked("hEffTritonTOF","TH1D");
   TH1D* hDefaultAntiTritonTOF = (TH1D*) DefaultFile.GetObjectChecked("hEffAntiTritonTOF","TH1D");

   const Int_t nPtBinsTOF=4;
   Double_t BinningTOF[]={1., 1.5, 2. ,2.5 ,3.0}; // triton (TPC-TOF)

   // rebin efficiency
   for (Int_t i=1; i<=hDefaultTritonTOF->GetNbinsX(); i++) {
      hDefaultTritonTOF->SetBinContent(i, hDefaultTritonTOF->GetBinContent(i)*hDefaultTritonTOF->GetBinWidth(i));
      hDefaultAntiTritonTOF->SetBinContent(i, hDefaultAntiTritonTOF->GetBinContent(i)*hDefaultAntiTritonTOF->GetBinWidth(i));
   }
   hDefaultTritonTOF = (TH1D*) hDefaultTritonTOF->Rebin(nPtBinsTOF,"hEffTritonTOF",BinningTOF);
   hDefaultAntiTritonTOF = (TH1D*) hDefaultAntiTritonTOF->Rebin(nPtBinsTOF,"hEffAntiTritonTOF",BinningTOF);

   hDefaultTritonTOF->Scale(1,"width");
   hDefaultAntiTritonTOF->Scale(1,"width");

   TH1D* hDefaultHelium = (TH1D*) DefaultFile.GetObjectChecked("hEffHelium","TH1D");
   TH1D* hDefaultAntiHelium = (TH1D*) DefaultFile.GetObjectChecked("hEffAntiHelium","TH1D");

   TFile GEANT4file ("Efficiency_GEANT4_TPC70cl.root","READ");
   if (!GEANT4file || GEANT4file.IsZombie()) {
      printf("Cannot open GEANT4 file \n");
      return;
   }
   TH1D* hTritonGEANT4 = (TH1D*) GEANT4file.GetObjectChecked("hEffTriton","TH1D");
   TH1D* hAntiTritonGEANT4 = (TH1D*) GEANT4file.GetObjectChecked("hEffAntiTriton","TH1D");;

   // rebin efficiency
   for (Int_t i=1; i<=hTritonGEANT4->GetNbinsX(); i++) {
      hTritonGEANT4->SetBinContent(i, hTritonGEANT4->GetBinContent(i)*hTritonGEANT4->GetBinWidth(i));
      hAntiTritonGEANT4->SetBinContent(i, hAntiTritonGEANT4->GetBinContent(i)*hAntiTritonGEANT4->GetBinWidth(i));
   }
   hTritonGEANT4 = (TH1D*) hTritonGEANT4->Rebin(nPtBins,"hEffTriton",Binning);
   hAntiTritonGEANT4 = (TH1D*) hAntiTritonGEANT4->Rebin(nPtBins,"hEffAntiTriton",Binning);

   hTritonGEANT4->Scale(1,"width");
   hAntiTritonGEANT4->Scale(1,"width");

   TH1D* hTritonTOFGEANT4 = (TH1D*) GEANT4file.GetObjectChecked("hEffTritonTOF","TH1D");
   TH1D* hAntiTritonTOFGEANT4 = (TH1D*) GEANT4file.GetObjectChecked("hEffAntiTritonTOF","TH1D");

   // rebin efficiency
   for (Int_t i=1; i<=hTritonTOFGEANT4->GetNbinsX(); i++) {
      hTritonTOFGEANT4->SetBinContent(i, hTritonTOFGEANT4->GetBinContent(i)*hTritonTOFGEANT4->GetBinWidth(i));
      hAntiTritonTOFGEANT4->SetBinContent(i, hAntiTritonTOFGEANT4->GetBinContent(i)*hAntiTritonTOFGEANT4->GetBinWidth(i));
   }
   hTritonTOFGEANT4 = (TH1D*) hTritonTOFGEANT4->Rebin(nPtBinsTOF,"hEffTritonTOF",BinningTOF);
   hAntiTritonTOFGEANT4 = (TH1D*) hAntiTritonTOFGEANT4->Rebin(nPtBinsTOF,"hEffAntiTritonTOF",BinningTOF);

   hTritonTOFGEANT4->Scale(1,"width");
   hAntiTritonTOFGEANT4->Scale(1,"width");


   TH1D* hHeliumGEANT4 = (TH1D*) GEANT4file.GetObjectChecked("hEffHelium","TH1D");
   TH1D* hAntiHeliumGEANT4 = (TH1D*) GEANT4file.GetObjectChecked("hEffAntiHelium","TH1D");

   //Get number of bins
   Double_t minpT = 1.;
   Double_t maxpT = 7.;
   Int_t Firstbin = hDefaultHelium->FindBin(minpT+0.01);
   Int_t Lastbin = hDefaultHelium->FindBin(maxpT-0.01);
   Int_t nbins = Lastbin-Firstbin+1;

   Double_t Uncertainty = 0;
   Double_t UncertaintyAnti = 0;

   for (Int_t ibin=Firstbin; ibin<=nbins; ibin++) {
      Uncertainty = TMath::Abs(hHeliumGEANT4->GetBinContent(ibin) - hDefaultHelium->GetBinContent(ibin))/(hDefaultHelium->GetBinContent(ibin))*0.5*100;
      UncertaintyAnti = TMath::Abs(hAntiHeliumGEANT4->GetBinContent(ibin) - hDefaultAntiHelium->GetBinContent(ibin))/(hDefaultAntiHelium->GetBinContent(ibin))*0.5*100;

      printf(" %.2f GeV/c: %.2f %% (Helium); %.2f %% (Anti-Helium) \n", hDefaultAntiHelium->GetXaxis()->GetBinCenter(ibin), Uncertainty, UncertaintyAnti);
   }
   printf("\n");

   minpT = 1.;
   maxpT= 3.;
   Int_t Firstbin = hDefaultTritonTOF->FindBin(minpT+0.01);
   Int_t Lastbin = hDefaultTritonTOF->FindBin(maxpT-0.01);
   Int_t nbins = Lastbin-Firstbin+1;

   for (Int_t ibin=Firstbin; ibin<=nbins; ibin++) {
      Uncertainty = TMath::Abs(hTritonTOFGEANT4->GetBinContent(ibin) - hDefaultTritonTOF->GetBinContent(ibin))/(hDefaultTritonTOF->GetBinContent(ibin))*0.5*100;
      UncertaintyAnti = TMath::Abs(hAntiTritonTOFGEANT4->GetBinContent(ibin) - hDefaultAntiTritonTOF->GetBinContent(ibin))/(hDefaultAntiTritonTOF->GetBinContent(ibin))*0.5*100;

      printf(" %.2f GeV/c: %.2f %% (Triton, TOF); %.2f %% (Anti-Triton, TOF) \n", hDefaultTritonTOF->GetXaxis()->GetBinCenter(ibin), Uncertainty, UncertaintyAnti);
   }
   printf("\n");

   minpT = 0.5;
   maxpT= 2.;
   //Get number of bins
   Int_t Firstbin = hDefaultTriton->FindBin(minpT+0.01);
   Int_t Lastbin = hDefaultTriton->FindBin(maxpT-0.01);
   Int_t nbins = Lastbin-Firstbin+1;

   for (Int_t ibin=Firstbin; ibin<=nbins; ibin++) {
      Uncertainty = TMath::Abs(hTritonGEANT4->GetBinContent(ibin) - hDefaultAntiTriton->GetBinContent(ibin))/(hDefaultAntiTriton->GetBinContent(ibin))*0.5*100;
      UncertaintyAnti = TMath::Abs(hAntiTritonGEANT4->GetBinContent(ibin) - hDefaultAntiTriton->GetBinContent(ibin))/(hDefaultAntiTriton->GetBinContent(ibin))*0.5*100;
      
      printf(" %.2f GeV/c:  %.2f %% (Triton); %.2f %% (Anti-Triton) \n", hDefaultTriton->GetXaxis()->GetBinCenter(ibin), Uncertainty, UncertaintyAnti);
   }

   TH1D* Axis = new TH1D("Axis","",75,0,7.5);
   setOPT_hists(Axis,pTLabel, yLable,510,20,0.5);
   Axis->GetYaxis()->SetRangeUser(0,1);

   TCanvas* cHelium = new TCanvas("cHelium","Efficiencies for different material budgets (helium)",800, 600);
   gPad->SetTicks();
	Axis->DrawCopy();

   setOPT_hists(hDefaultHelium,pTLabel, yLable,510,20,0.75,kViolet-3);
   setOPT_hists(hHeliumGEANT4,pTLabel, yLable,510,21,0.75,kBlue);
	hDefaultHelium->DrawCopy("same");
   hHeliumGEANT4->DrawCopy("same");

   setOPT_hists(hDefaultAntiHelium,pTLabel, yLable,510,20,0.75,kOrange-3);
   setOPT_hists(hAntiHeliumGEANT4,pTLabel, yLable,510,21,0.75,kRed);
   hDefaultAntiHelium->DrawCopy("same");
   hAntiHeliumGEANT4->DrawCopy("same");

   TLegend* leg = plotLegend("right_bottom","(Anti-)helium efficiencies in Pb-Pb @ 5.02 TeV", 1.2, 1., -0.05, 0.02);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->AddEntry((TObject*)0,"#kern[-0.32]{(50 - 90) % centrality}","");
   leg->AddEntry(hDefaultHelium,"^{3}He - GEANT3 (LHC16h7c)");
   leg->AddEntry(hHeliumGEANT4,"^{3}He - GEANT4 (LHC16h7c_g4_2)");
   leg->AddEntry(hDefaultAntiHelium,"^{3}#bar{He} - GEANT3 (LHC16h7c)");
   leg->AddEntry(hAntiHeliumGEANT4,"^{3}#bar{He} - GEANT4 (LHC16h7c_g4_2)");
   leg->AddEntry((TObject*)0,"Work in progress","");
   leg->DrawClone();

   TCanvas* cTritonTOF = new TCanvas("cTritonTOF","Efficiencies for different material budgets (triton, TOF-TPC)",800, 600);
   gPad->SetTicks();
   Axis->GetXaxis()->SetRangeUser(0,3.5);
   Axis->DrawCopy();

   setOPT_hists(hDefaultTritonTOF,pTLabel, yLable,510,20,0.75,kViolet-3);
   setOPT_hists(hTritonTOFGEANT4,pTLabel, yLable,510,21,0.75,kBlue);
   hDefaultTritonTOF->DrawCopy("same");
   hTritonTOFGEANT4->DrawCopy("same");

   setOPT_hists(hDefaultAntiTritonTOF,pTLabel, yLable,510,20,0.75,kOrange-3);
   setOPT_hists(hAntiTritonTOFGEANT4,pTLabel, yLable,510,21,0.75,kRed);
   hDefaultAntiTritonTOF->DrawCopy("same");
   hAntiTritonTOFGEANT4->DrawCopy("same");

   TLegend* leg = plotLegend("left_top","(Anti-)triton efficiencies in Pb-Pb @ 5.02 TeV (TOF-TPC)", 1.2, 1., 0.05, -0.05);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->AddEntry((TObject*)0,"#kern[-0.32]{(50 - 90) % centrality}","");
   leg->AddEntry(hDefaultTritonTOF,"t - GEANT3 (LHC16h7c)");
   leg->AddEntry(hTritonTOFGEANT4,"t - GEANT4 (LHC16h7c_g4_2)");
   leg->AddEntry(hDefaultAntiTritonTOF,"#bar{t} - GEANT3 (LHC16h7c)");
   leg->AddEntry(hAntiTritonTOFGEANT4,"#bar{t} - GEANT4 (LHC16h7c_g4_2)");
   leg->AddEntry((TObject*)0,"Work in progress","");
   leg->DrawClone();

   TCanvas* cTriton = new TCanvas("cTriton","Efficiencies for different material budgets (triton)",800, 600);
   Axis->GetXaxis()->SetRangeUser(0,2.5);
   Axis->DrawCopy();
   
   setOPT_hists(hDefaultTriton,pTLabel, yLable,510,20,0.75,kViolet-3);
   setOPT_hists(hTritonGEANT4,pTLabel, yLable,510,21,0.75,kBlue);
   hDefaultTriton->DrawCopy("same");
   hTritonGEANT4->DrawCopy("same");
   
   setOPT_hists(hDefaultAntiTriton,pTLabel, yLable,510,20,0.75,kOrange-3);
   setOPT_hists(hAntiTritonGEANT4,pTLabel, yLable,510,21,0.75,kRed);
   hDefaultAntiTriton->DrawCopy("same");
   hAntiTritonGEANT4->DrawCopy("same");

   TLegend* leg2 = plotLegend("left_top","(Anti-)triton efficiencies in Pb-Pb @ 5.02 TeV (TPC)", 1.2, 0.2, 0.05, -0.05);
   leg2->SetFillStyle(0);
   leg2->SetBorderSize(0);
   leg2->SetTextSize(0.03);
   leg2->AddEntry((TObject*)0,"#kern[-0.34]{(50 - 90) % centrality}","");
   leg2->DrawClone();

   TLegend* leg = plotLegend("right_bottom","", 1., 0.8, -0.05, 0.02);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->AddEntry(hDefaultTriton,"t - GEANT3 (LHC16h7c)");
   leg->AddEntry(hTritonGEANT4,"t - GEANT4 (LHC16h7c_g4_2)");
   leg->AddEntry(hDefaultAntiTriton,"#bar{t} - GEANT3 (LHC16h7c)");
   leg->AddEntry(hAntiTritonGEANT4,"#bar{t} - GEANT4 (LHC16h7c_g4_2)");
   leg->AddEntry((TObject*)0,"Work in progress","");
   leg->DrawClone();


//// Scaling GEANT3 to GEANT4
//  TCanvas* cAbsorptionCorrection = new TCanvas("cAbsorptionCorrection","Absorption correction", 800,600);
//  gPad->SetTicks();
//
//  TH1D* hEffRatioHe = (TH1D*) hAntiHeliumGEANT4->Clone("hEffRatioHe");
//  hEffRatioHe->Divide(hDefaultAntiHelium);
//  TH1D* hEffRatioAntiHe = (TH1D*) hHeliumGEANT4->Clone("hEffRatioAntiHe");
//  hEffRatioAntiHe->Divide(hDefaultHelium);
//
//  TH1D* hCorrect= (TH1D*) hEffRatioAntiHe->Clone("hCorrect");
//  hCorrect->Divide(hEffRatioHe);
//  setOPT_hists(hCorrect,pTLabel, "Correction factor",510,21,0.75,kBlack);
//
//  hEffRatioHe->DrawCopy("E");
//  hEffRatioAntiHe->DrawCopy("Esame");
//  hCorrect->DrawCopy("Esame");
//
//  TFile DataFile ("/Users/tidus/Documents/Studium/PhD/Nuclei/GITRepository/AliAnalysisHe/reducedTreeNuclei/OfflineFramework/Results/Helium_LHC16q_FAST_woSDD_Mult_Int.root");
//  TString Objectname = "fHistoHe3PtSpectrum"; // 3He
//  TString Objectname2 = "fHistoAntiHe3PtSpectrum"; // Anti-3He
//
//  TH1D* hHe = (TH1D*) DataFile.GetObjectChecked(Objectname,"TH1D");
//  setOPT_hists(hHe,pTLabel, yLable,510,20,0.5,kBlue);
//  hHe->SetTitle("");
//  hHe->SetName("hHe");
//  TH1D* hHeSyst = (TH1D*) DataFile.GetObjectChecked(Objectname+"Syst","TH1D");
//  setOPT_hists(hHeSyst,pTLabel, yLable,510,20,0.5,kBlue);
//  hHeSyst->SetFillStyle(0);
//
//  TH1D* hAntiHe = (TH1D*) DataFile.GetObjectChecked(Objectname2,"TH1D");
//  setOPT_hists(hAntiHe,pTLabel, yLable,510,20,0.5,kRed);
//  hAntiHe->SetTitle("");
//  hAntiHe->SetName("hAntiHe");
//  TH1D* hAntiHeSyst = (TH1D*) DataFile.GetObjectChecked(Objectname2+"Syst","TH1D");
//  setOPT_hists(hAntiHeSyst,pTLabel, yLable,510,20,0.5,kRed);
//  hAntiHeSyst->SetFillStyle(0);
//
//  const Int_t nPtBins=8;
//  Double_t BinningData[]={1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6.,7.};
//
//  // rebin efficiency
//  for (Int_t i=1; i<=hEffRatioHe->GetNbinsX(); i++) {
//    hEffRatioHe->SetBinContent(i, hEffRatioHe->GetBinContent(i)*hEffRatioHe->GetBinWidth(i));
//    hEffRatioHe->SetBinError(i, hEffRatioHe->GetBinError(i)*hEffRatioHe->GetBinWidth(i));
//    hEffRatioAntiHe->SetBinContent(i, hEffRatioAntiHe->GetBinContent(i)*hEffRatioAntiHe->GetBinWidth(i));
//    hEffRatioAntiHe->SetBinError(i, hEffRatioAntiHe->GetBinError(i)*hEffRatioAntiHe->GetBinWidth(i));
//  }
//  hEffRatioHe = (TH1D*) hEffRatioHe->Rebin(nPtBins,"hEffRatioHe",BinningData);
//  hEffRatioAntiHe = (TH1D*) hEffRatioAntiHe->Rebin(nPtBins,"hEffRatioAntiHe",BinningData);
//
//  hEffRatioHe->Scale(1,"width");
//  hEffRatioAntiHe->Scale(1,"width");
//
//  TH1D* hHeSystUncor = (TH1D*) hHeSyst->Clone("hHeSystUncor");
//  TH1D* hAntiHeSystUncor = (TH1D*) hAntiHeSyst->Clone("hAntiHeSystUncor");
//
//  hHeSyst->Divide(hEffRatioHe);
//  hHe->Divide(hEffRatioHe);
//
//  hAntiHeSyst->Divide(hEffRatioAntiHe);
//  hAntiHe->Divide(hEffRatioAntiHe);
//
//  TCanvas* cCorrected = new TCanvas("cCorrected","Corrected Spectra",800,600);
//  gPad->SetTicks();
//  gPad->SetLogy();
//  hHeSyst->DrawCopy("E2");
//  hAntiHeSyst->DrawCopy("E2same");
//  hHe->DrawCopy("EX0same");
//  hAntiHe->DrawCopy("EX0same");
//
//  TCanvas* cChange = new TCanvas("cChange","Change to spectra",800,600);
//  gPad->SetTicks();
//  hHeSyst->Divide(hHeSystUncor);
//  hAntiHeSyst->Divide(hAntiHeSystUncor);
//  hHeSyst->DrawCopy("E2");
//  hAntiHeSyst->DrawCopy("E2same");

}
