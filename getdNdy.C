#include <iostream>

using namespace std;

Double_t
LevyTsallis_Func(const Double_t *x, const Double_t *p)
{
  /* dN/dpt */
  
  Double_t pt = x[0];
  Double_t mass = p[0];
  Double_t mt = TMath::Sqrt(pt * pt + mass * mass);
  Double_t n = p[1];
  Double_t C = p[2];
  Double_t norm = p[3];

  Double_t part1 = (n - 1.) * (n - 2.);
  Double_t part2 = n * C * (n * C + mass * (n - 2.));
  Double_t part3 = part1 / part2;
  Double_t part4 = 1. + (mt - mass) / n / C;
  Double_t part5 = TMath::Power(part4, -n);
  return pt * norm * part3 * part5;
}

TF1 *
LevyTsallis(const Char_t *name, Double_t mass, Double_t n = 5., Double_t C = 0.1, Double_t norm = 1.)
{
  
  TF1 *fLevyTsallis = new TF1(name, LevyTsallis_Func, 0., 10., 4);
  fLevyTsallis->SetParameters(mass, n, C, norm);
  fLevyTsallis->SetParNames("mass", "n", "C", "norm");
  fLevyTsallis->FixParameter(0, mass);
  fLevyTsallis->SetParLimits(1, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(2, 1.e-3, 1.e3);
  fLevyTsallis->SetParLimits(3, 1.e-6, 1.e8);
  return fLevyTsallis;
}

void getdNdy() {

	gROOT->ProcessLine(".L YieldMean.C");

	TDirectoryFile* spectra[6] = { Jetty_V0M, Reference_V0M, Isotropic_V0M, 
		Jetty_Trks, Reference_Trks, Isotropic_Trks };

	Double_t dndy[6]; Double_t dndy_stat[6]; Double_t dndy_sys[6];
	Double_t pt[6]; Double_t pt_stat[6]; Double_t pt_sys[6];

	TF1* func = LevyTsallis("LT",0.13957);
	TH1* hout;

	for (int iS = 0; iS < 6; iS++) {
		spectra[iS]->cd();
		hout = YieldMean(hPion_stat,hPion_tot_sys,hPion_unc_sys,func,0., 20.);
		hout->Draw();
		printf("yield %f +- %f , <pt> %f +- %f \n", hout->GetBinContent(kYield), hout->GetBinContent(kYieldStat),
		hout->GetBinContent(kMean), hout->GetBinContent(kMeanStat) );


		dndy[iS] 		= hout->GetBinContent(kYield);
		dndy_stat[iS] 	= hout->GetBinContent(kYieldStat);
		dndy_sys[iS] 	= hout->GetBinContent(kYieldSysHi);

		pt[iS] 			= hout->GetBinContent(kMean);
		pt_stat[iS] 	= hout->GetBinContent(kMeanStat);
		pt_sys[iS] 		= hout->GetBinContent(kMeanSysHi);
	}

	TCanvas* c1 = new TCanvas("c1","",800,800);
	TGraphErrors* gr[6]; TGraphErrors* gr_sys[6]; 
	
	for (int iS = 0; iS < 6; ++iS)	{
		
		printf("dndy: %f +- %f +- %f ; pt: %f +- %f +- %f \n", 
			dndy[iS], dndy_stat[iS], dndy_sys[iS],
			pt[iS], pt_stat[iS], pt_sys[iS]);

		
		gr[iS] = new TGraphErrors(1, &dndy[iS], &pt[iS], &dndy_stat[iS], &pt_stat[iS]);
		gr_sys[iS] = new TGraphErrors(1, &dndy[iS], &pt[iS], &dndy_sys[iS], &pt_sys[iS]);
		(iS < 3) ? gr[iS]->SetMarkerStyle(20) : gr[iS]->SetMarkerStyle(33);
		(iS < 3) ? gr[iS]->SetMarkerSize(1.8) : gr[iS]->SetMarkerSize(2.8);
		if (iS == 0 || iS == 3) gr[iS]->SetMarkerColor(kRed);
		if (iS == 1 || iS == 4) gr[iS]->SetMarkerColor(kBlack);
		if (iS == 2 || iS == 5) gr[iS]->SetMarkerColor(kBlue);
		if (iS == 0 || iS == 3) gr_sys[iS]->SetMarkerColor(kRed);
		if (iS == 1 || iS == 4) gr_sys[iS]->SetMarkerColor(kBlack);
		if (iS == 2 || iS == 5) gr_sys[iS]->SetMarkerColor(kBlue);
		
	}

	gr[0]->SetTitle("V0M 0-10%, Jetty"); gr[1]->SetTitle("V0M 0-10%"); gr[2]->SetTitle("V0M 0-10%, Iso");
	gr[3]->SetTitle("CL1 0-10%, Jetty"); gr[4]->SetTitle("CL1 0-10%"); gr[5]->SetTitle("CL1 0-10%, Iso");

	//TGraphErrors* gr_tot = new TGraphErrors(6, dndy, pt, dndy_stat, pt_stat);
	TMultiGraph* gr_tot = new TMultiGraph();
	for (int iS = 0; iS < 6; iS++) { gr_tot->Add(gr_sys[iS],"2"); gr_tot->Add(gr[iS],"");}
	gr_tot->SetTitle("; <dN_{#pi}/dy>; <p_{T}>");
	gr_tot->SetMinimum(0.55); gr_tot->SetMaximum(0.65);
	gr_tot->Draw("APZ"); 
	gr_tot->GetXaxis()->SetLimits(10., 22.);
	gPad->Modified();
	gPad->Update();

	TLegend* leg = c1->BuildLegend();
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	leg->SetMargin(0.25);
	leg->SetEntrySeparation(0.5);
	leg->SetNColumns(2);
	leg->Clear();
	leg->AddEntry(gr[0],gr[0]->GetTitle(),"p");
	leg->AddEntry(gr[3],gr[3]->GetTitle(),"p");
	leg->AddEntry(gr[1],gr[1]->GetTitle(),"p");
	leg->AddEntry(gr[4],gr[4]->GetTitle(),"p");
	leg->AddEntry(gr[2],gr[2]->GetTitle(),"p");
	leg->AddEntry(gr[5],gr[5]->GetTitle(),"p");
	leg->Draw();

}

