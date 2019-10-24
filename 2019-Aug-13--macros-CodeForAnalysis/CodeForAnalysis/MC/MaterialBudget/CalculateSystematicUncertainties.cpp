//
//  MakeOverlayHFE.cpp - 3 spectra
//
//
//  Created by tidus on 30.03.16.
//  Please move data to the same folder as the macro
//

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

   TFile DefaultFile ("Efficiency_LHC16h7c.root","READ");
   if (!DefaultFile || DefaultFile.IsZombie()) {
      printf("Cannot open Efficiency_LHC16h7c.root \n");
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

   TH1D* hVariationTriton[2];
   TH1D* hVariationAntiTriton[2];

   TH1D* hVariationTritonTOF[2];
   TH1D* hVariationAntiTritonTOF[2];

   TH1D* hVariationHelium[2];
   TH1D* hVariationAntiHelium[2];

   TFile VariationFile1 ("Efficiency_LHC17d5a.root","READ");
   if (!VariationFile1 || VariationFile1.IsZombie()) {
      printf("Cannot open Efficiency_LHC17d5a.root \n");
      return;
   }
   hVariationTriton[0] = (TH1D*) VariationFile1.GetObjectChecked("hEffTriton","TH1D");
   hVariationAntiTriton[0] = (TH1D*) VariationFile1.GetObjectChecked("hEffAntiTriton","TH1D");

   // rebin efficiency
   for (Int_t i=1; i<=hVariationTriton[0]->GetNbinsX(); i++) {
      hVariationTriton[0]->SetBinContent(i, hVariationTriton[0]->GetBinContent(i)*hVariationTriton[0]->GetBinWidth(i));
      hVariationAntiTriton[0]->SetBinContent(i, hVariationAntiTriton[0]->GetBinContent(i)*hVariationAntiTriton[0]->GetBinWidth(i));
   }
   hVariationTriton[0] = (TH1D*) hVariationTriton[0]->Rebin(nPtBins,"hEffTriton",Binning);
   hVariationAntiTriton[0] = (TH1D*) hVariationAntiTriton[0]->Rebin(nPtBins,"hEffAntiTriton",Binning);

   hVariationTriton[0]->Scale(1,"width");
   hVariationAntiTriton[0]->Scale(1,"width");

   hVariationTritonTOF[0] = (TH1D*) VariationFile1.GetObjectChecked("hEffTritonTOF","TH1D");
   hVariationAntiTritonTOF[0] = (TH1D*) VariationFile1.GetObjectChecked("hEffAntiTritonTOF","TH1D");

   // rebin efficiency
   for (Int_t i=1; i<=hVariationTritonTOF[0]->GetNbinsX(); i++) {
      hVariationTritonTOF[0]->SetBinContent(i, hVariationTritonTOF[0]->GetBinContent(i)*hVariationTritonTOF[0]->GetBinWidth(i));
      hVariationAntiTritonTOF[0]->SetBinContent(i, hVariationAntiTritonTOF[0]->GetBinContent(i)*hVariationAntiTritonTOF[0]->GetBinWidth(i));
   }
   hVariationTritonTOF[0] = (TH1D*) hVariationTritonTOF[0]->Rebin(nPtBinsTOF,"hEffTritonTOF",BinningTOF);
   hVariationAntiTritonTOF[0] = (TH1D*) hVariationAntiTritonTOF[0]->Rebin(nPtBinsTOF,"hEffAntiTritonTOF",BinningTOF);

   hVariationTritonTOF[0]->Scale(1,"width");
   hVariationAntiTritonTOF[0]->Scale(1,"width");

   hVariationHelium[0] = (TH1D*) VariationFile1.GetObjectChecked("hEffHelium","TH1D");
   hVariationAntiHelium[0] = (TH1D*) VariationFile1.GetObjectChecked("hEffAntiHelium","TH1D");
   
   TFile VariationFile2 ("Efficiency_LHC17d5b.root","READ");
   if (!VariationFile2 || VariationFile2.IsZombie()) {
      printf("Cannot open Efficiency_LHC17d5b.root \n");
      return;
   }
   hVariationTriton[1] = (TH1D*) VariationFile2.GetObjectChecked("hEffTriton","TH1D");
   hVariationAntiTriton[1] = (TH1D*) VariationFile2.GetObjectChecked("hEffAntiTriton","TH1D");

   // rebin efficiency
   for (Int_t i=1; i<=hVariationTriton[1]->GetNbinsX(); i++) {
      hVariationTriton[1]->SetBinContent(i, hVariationTriton[1]->GetBinContent(i)*hVariationTriton[1]->GetBinWidth(i));
      hVariationAntiTriton[1]->SetBinContent(i, hVariationAntiTriton[1]->GetBinContent(i)*hVariationAntiTriton[1]->GetBinWidth(i));
   }
   hVariationTriton[1] = (TH1D*) hVariationTriton[1]->Rebin(nPtBins,"hEffTriton",Binning);
   hVariationAntiTriton[1] = (TH1D*) hVariationAntiTriton[1]->Rebin(nPtBins,"hEffAntiTriton",Binning);

   hVariationTriton[1]->Scale(1,"width");
   hVariationAntiTriton[1]->Scale(1,"width");

   hVariationTritonTOF[1] = (TH1D*) VariationFile2.GetObjectChecked("hEffTritonTOF","TH1D");
   hVariationAntiTritonTOF[1] = (TH1D*) VariationFile2.GetObjectChecked("hEffAntiTritonTOF","TH1D");

   // rebin efficiency
   for (Int_t i=1; i<=hVariationTritonTOF[1]->GetNbinsX(); i++) {
      hVariationTritonTOF[1]->SetBinContent(i, hVariationTritonTOF[1]->GetBinContent(i)*hVariationTritonTOF[1]->GetBinWidth(i));
      hVariationAntiTritonTOF[1]->SetBinContent(i, hVariationAntiTritonTOF[1]->GetBinContent(i)*hVariationAntiTritonTOF[1]->GetBinWidth(i));
   }
   hVariationTritonTOF[1] = (TH1D*) hVariationTritonTOF[1]->Rebin(nPtBinsTOF,"hEffTritonTOF",BinningTOF);
   hVariationAntiTritonTOF[1] = (TH1D*) hVariationAntiTritonTOF[1]->Rebin(nPtBinsTOF,"hEffAntiTritonTOF",BinningTOF);

   hVariationTritonTOF[1]->Scale(1,"width");
   hVariationAntiTritonTOF[1]->Scale(1,"width");

   hVariationHelium[1] = (TH1D*) VariationFile2.GetObjectChecked("hEffHelium","TH1D");
   hVariationAntiHelium[1] = (TH1D*) VariationFile2.GetObjectChecked("hEffAntiHelium","TH1D");

   //Get number of bins
   Double_t minpT = 1.;
   Double_t maxpT = 7.;
   Int_t Firstbin = hDefaultHelium->FindBin(minpT+0.01);
   Int_t Lastbin = hDefaultHelium->FindBin(maxpT-0.01);
   Int_t nbins = Lastbin-Firstbin+1;

   Double_t Minimum = 0;
   Double_t Maximum = 0;
   Double_t Uncertainty = 0;
   Double_t UncertaintyAnti = 0;

   const Int_t index[3] = {1,2,3};
   Double_t Array[3];
   for (Int_t ibin=Firstbin; ibin<=nbins; ibin++) {
      Double_t Array[3] = {hVariationHelium[0]->GetBinContent(ibin),hDefaultHelium->GetBinContent(ibin),hVariationHelium[1]->GetBinContent(ibin)};
      Array= TMath::Sort(3,Array,index,0);
      Uncertainty = (Array[index[2]]-Array[index[0]])/(2.*hDefaultHelium->GetBinContent(ibin))*100;

      Double_t Array[3] = {hVariationAntiHelium[0]->GetBinContent(ibin),hDefaultAntiHelium->GetBinContent(ibin),hVariationAntiHelium[1]->GetBinContent(ibin)};
      Array= TMath::Sort(3,Array,index,0);
      UncertaintyAnti = (Array[index[2]]-Array[index[0]])/(2.*hDefaultAntiHelium->GetBinContent(ibin))*100;

      printf(" %.2f GeV/c: %.2f %% (Helium); %.2f %% (Anti-Helium) \n", hDefaultAntiHelium->GetXaxis()->GetBinCenter(ibin), Uncertainty, UncertaintyAnti);
   }
   printf("\n");

   minpT = 1.;
   maxpT= 3.;
   Int_t Firstbin = hDefaultTritonTOF->FindBin(minpT+0.01);
   Int_t Lastbin = hDefaultTritonTOF->FindBin(maxpT-0.01);
   Int_t nbins = Lastbin-Firstbin+1;
   for (Int_t ibin=Firstbin; ibin<=nbins; ibin++) {
      Double_t Array[3] = {hVariationTritonTOF[0]->GetBinContent(ibin),hDefaultTritonTOF->GetBinContent(ibin),hVariationTritonTOF[1]->GetBinContent(ibin)};
      Array= TMath::Sort(3,Array,index,0);
      Uncertainty = (Array[index[2]]-Array[index[0]])/(2.*hDefaultTritonTOF->GetBinContent(ibin))*100;

      Double_t Array[3] = {hVariationAntiTritonTOF[0]->GetBinContent(ibin),hDefaultAntiTritonTOF->GetBinContent(ibin),hVariationAntiTritonTOF[1]->GetBinContent(ibin)};
      Array= TMath::Sort(3,Array,index,0);
      UncertaintyAnti = (Array[index[2]]-Array[index[0]])/(2.*hDefaultAntiTritonTOF->GetBinContent(ibin))*100;

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
      Double_t Array[3] = {hVariationTriton[0]->GetBinContent(ibin),hDefaultTriton->GetBinContent(ibin),hVariationTriton[1]->GetBinContent(ibin)};
      Array= TMath::Sort(3,Array,index,0);
      Uncertainty = (Array[index[2]]-Array[index[0]])/(2.*hDefaultTriton->GetBinContent(ibin))*100;
      
      Double_t Array[3] = {hVariationAntiTriton[0]->GetBinContent(ibin),hDefaultAntiTriton->GetBinContent(ibin),hVariationAntiTriton[1]->GetBinContent(ibin)};
      Array= TMath::Sort(3,Array,index,0);
      UncertaintyAnti = (Array[index[2]]-Array[index[0]])/(2.*hDefaultAntiTriton->GetBinContent(ibin))*100;
      
      printf(" %.2f GeV/c:  %.2f %% (Triton); %.2f %% (Anti-Triton) \n", hDefaultTriton->GetXaxis()->GetBinCenter(ibin), Uncertainty, UncertaintyAnti);
   }

   TH1D* Axis = new TH1D("Axis","",75,0,7.5);
   setOPT_hists(Axis,pTLabel, yLable,510,20,0.5);
   Axis->GetYaxis()->SetRangeUser(0,1);

   TCanvas* cHelium = new TCanvas("cHelium","Efficiencies for different material budgets (helium)",800, 600);
   gPad->SetTicks();
	Axis->DrawCopy();

   setOPT_hists(hDefaultHelium,pTLabel, yLable,510,20,0.75,kViolet-3);
   setOPT_hists(hVariationHelium[0],pTLabel, yLable,510,21,0.75,kBlue);
   setOPT_hists(hVariationHelium[1],pTLabel, yLable,510,22,0.75,kAzure-2);
	hDefaultHelium->DrawCopy("same");
   hVariationHelium[0]->DrawCopy("same");
   hVariationHelium[1]->DrawCopy("same");

   setOPT_hists(hDefaultAntiHelium,pTLabel, yLable,510,20,0.75,kOrange-3);
   setOPT_hists(hVariationAntiHelium[0],pTLabel, yLable,510,21,0.75,kRed);
   setOPT_hists(hVariationAntiHelium[1],pTLabel, yLable,510,22,0.75,kPink-2);
   hDefaultAntiHelium->DrawCopy("same");
   hVariationAntiHelium[0]->DrawCopy("same");
   hVariationAntiHelium[1]->DrawCopy("same");

   TLegend* leg = plotLegend("right_bottom","(Anti-)helium efficiencies in Pb-Pb @ 5.02 TeV", 1.2, 1., -0.05, 0.02);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->AddEntry((TObject*)0,"#kern[-0.32]{(50 - 90) % centrality}","");
   leg->AddEntry(hDefaultHelium,"^{3}He - Normal matrial budget (LHC16h7c)");
   leg->AddEntry(hVariationHelium[0],"^{3}He - Matrial budget + 4.5 % (LHC17d5a)");
   leg->AddEntry(hVariationHelium[1],"^{3}He - Matrial budget - 4.5 % (LHC17d5b)");
   leg->AddEntry(hDefaultAntiHelium,"^{3}#bar{He} - Normal matrial budget (LHC16h7c)");
   leg->AddEntry(hVariationAntiHelium[0],"^{3}#bar{He} - Matrial budget + 4.5 % (LHC17d5a)");
   leg->AddEntry(hVariationAntiHelium[1],"^{3}#bar{He} - Matrial budget - 4.5 % (LHC17d5b)");
   leg->AddEntry((TObject*)0,"Work in progress","");
   leg->DrawClone();

   TCanvas* cTritonTOF = new TCanvas("cTritonTOF","Efficiencies for different material budgets (triton, TOF-TPC)",800, 600);
   gPad->SetTicks();
   Axis->GetXaxis()->SetRangeUser(0,3.5);
   Axis->DrawCopy();

   setOPT_hists(hDefaultTritonTOF,pTLabel, yLable,510,20,0.75,kViolet-3);
   setOPT_hists(hVariationTritonTOF[0],pTLabel, yLable,510,21,0.75,kBlue);
   setOPT_hists(hVariationTritonTOF[1],pTLabel, yLable,510,22,0.75,kAzure-2);
   hDefaultTritonTOF->DrawCopy("same");
   hVariationTritonTOF[0]->DrawCopy("same");
   hVariationTritonTOF[1]->DrawCopy("same");

   setOPT_hists(hDefaultAntiTritonTOF,pTLabel, yLable,510,20,0.75,kOrange-3);
   setOPT_hists(hVariationAntiTritonTOF[0],pTLabel, yLable,510,21,0.75,kRed);
   setOPT_hists(hVariationAntiTritonTOF[1],pTLabel, yLable,510,22,0.75,kPink-2);
   hDefaultAntiTritonTOF->DrawCopy("same");
   hVariationAntiTritonTOF[0]->DrawCopy("same");
   hVariationAntiTritonTOF[1]->DrawCopy("same");

   TLegend* leg = plotLegend("left_top","(Anti-)triton efficiencies in Pb-Pb @ 5.02 TeV (TOF-TPC)", 1.2, 1., 0.05, -0.05);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->AddEntry((TObject*)0,"#kern[-0.32]{(50 - 90) % centrality}","");
   leg->AddEntry(hDefaultTritonTOF,"t - Normal matrial budget (LHC16h7c)");
   leg->AddEntry(hVariationTritonTOF[0],"t - Matrial budget + 4.5 % (LHC17d5a)");
   leg->AddEntry(hVariationTritonTOF[1],"t - Matrial budget - 4.5 % (LHC17d5b)");
   leg->AddEntry(hDefaultAntiTritonTOF,"#bar{t} - Normal matrial budget (LHC16h7c)");
   leg->AddEntry(hVariationAntiTritonTOF[0],"#bar{t} - Matrial budget + 4.5 % (LHC17d5a)");
   leg->AddEntry(hVariationAntiTritonTOF[1],"#bar{t} - Matrial budget - 4.5 % (LHC17d5b)");
   leg->AddEntry((TObject*)0,"Work in progress","");
   leg->DrawClone();

   TCanvas* cTriton = new TCanvas("cTriton","Efficiencies for different material budgets (triton)",800, 600);
   Axis->GetXaxis()->SetRangeUser(0,2.5);
   Axis->DrawCopy();
   
   setOPT_hists(hDefaultTriton,pTLabel, yLable,510,20,0.75,kViolet-3);
   setOPT_hists(hVariationTriton[0],pTLabel, yLable,510,21,0.75,kBlue);
   setOPT_hists(hVariationTriton[1],pTLabel, yLable,510,22,0.75,kAzure-2);
   hDefaultTriton->DrawCopy("same");
   hVariationTriton[0]->DrawCopy("same");
   hVariationTriton[1]->DrawCopy("same");
   
   setOPT_hists(hDefaultAntiTriton,pTLabel, yLable,510,20,0.75,kOrange-3);
   setOPT_hists(hVariationAntiTriton[0],pTLabel, yLable,510,21,0.75,kRed);
   setOPT_hists(hVariationAntiTriton[1],pTLabel, yLable,510,22,0.75,kPink-2);
   hDefaultAntiTriton->DrawCopy("same");
   hVariationAntiTriton[0]->DrawCopy("same");
   hVariationAntiTriton[1]->DrawCopy("same");

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
   leg->AddEntry(hDefaultTriton,"t - Normal matrial budget (LHC16h7c)");
   leg->AddEntry(hVariationTriton[0],"t - Matrial budget + 4.5 % (LHC17d5a)");
   leg->AddEntry(hVariationTriton[1],"t - Matrial budget - 4.5 % (LHC17d5b)");
   leg->AddEntry(hDefaultAntiTriton,"#bar{t} - Normal matrial budget (LHC16h7c)");
   leg->AddEntry(hVariationAntiTriton[0],"#bar{t} - Matrial budget + 4.5 % (LHC17d5a)");
   leg->AddEntry(hVariationAntiTriton[1],"#bar{t} - Matrial budget - 4.5 % (LHC17d5b)");
   leg->AddEntry((TObject*)0,"Work in progress","");
   leg->DrawClone();
}
