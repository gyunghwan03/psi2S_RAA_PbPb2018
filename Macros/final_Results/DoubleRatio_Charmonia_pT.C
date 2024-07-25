#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;

void DoubleRatio_Charmonia_pT()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid_2S = new TFile(Form("roots/RAA_psi2S_midRap_pT_%.f.root",40.));
	TFile *f_fwd_2S = new TFile(Form("roots/RAA_psi2S_forRap_pT_%.f.root",40.));

	TFile *f_mid_Jpsi = new TFile("./roots/RAA_JPsi_midRap_pT.root");
	TFile *f_fwd_Jpsi = new TFile("./roots/RAA_JPsi_forRap_pT.root");

	TH1D *h_midPR_2S = (TH1D*) f_mid_2S->Get("hRAA_PR");
	TH1D *h_fwdPR_2S = (TH1D*) f_fwd_2S->Get("hRAA_PR");
	TH1D *h_midPR_Jpsi = (TH1D*) f_mid_Jpsi->Get("hRAA_PR");
	TH1D *h_fwdPR_Jpsi = (TH1D*) f_fwd_Jpsi->Get("hRAA_PR");

	const int nPtBins_mid = 6;
	const int nPtBins_fwd = 4;
	double ptBin_mid[nPtBins_mid+1] = {6.5,9,12,15,20,25,40};
	double ptBin_fwd[nPtBins_fwd+1] = {3.5,6.5,9,12,40};

	double x_mid[nPtBins_mid]; double binWidth_mid[nPtBins_mid];
	double x_fwd[nPtBins_fwd]; double binWidth_fwd[nPtBins_fwd];

	double midPR_2S[nPtBins_mid];
	double midPR_2S_Err[nPtBins_mid];
	double midPR_2S_Sys[nPtBins_mid];
	double fwdPR_2S[nPtBins_fwd];
	double fwdPR_2S_Err[nPtBins_fwd];
	double fwdPR_2S_Sys[nPtBins_fwd];

	double midPR_Jpsi[nPtBins_mid];
	double midPR_Jpsi_Err[nPtBins_mid];
	double midPR_Jpsi_Sys[nPtBins_mid];
	double fwdPR_Jpsi[nPtBins_fwd];
	double fwdPR_Jpsi_Err[nPtBins_fwd];
	double fwdPR_Jpsi_Sys[nPtBins_fwd];

	double DR_midPR[nPtBins_mid];
	double DR_fwdPR[nPtBins_fwd];
	double DR_midPR_Err[nPtBins_mid];
	double DR_fwdPR_Err[nPtBins_fwd];

	for (int i=0; i<nPtBins_mid; i++){
		midPR_2S[i] = h_midPR_2S->GetBinContent(i+1);
		midPR_2S_Err[i] = h_midPR_2S->GetBinError(i+1);
		midPR_Jpsi[i] = h_midPR_Jpsi->GetBinContent(i+1);
		midPR_Jpsi_Err[i] = h_midPR_Jpsi->GetBinError(i+1);
		x_mid[i] = (ptBin_mid[i+1]+ptBin_mid[i])/2;
		binWidth_mid[i] = (ptBin_mid[i+1]-ptBin_mid[i])/2;
		DR_midPR[i] = midPR_2S[i]/midPR_Jpsi[i];
		DR_midPR_Err[i] = TMath::Sqrt(TMath::Power(midPR_2S_Err[i],2)+TMath::Power(midPR_Jpsi_Err[i],2));

		//SysPR[i] = midPR_2S[i]*(hSys_PR->GetBinContent(i+1));
	}

	for (int i=0; i<nPtBins_fwd; i++){
		fwdPR_2S[i] = h_fwdPR_2S->GetBinContent(i+1);
		fwdPR_2S_Err[i] = h_fwdPR_2S->GetBinError(i+1);
		fwdPR_Jpsi[i] = h_fwdPR_Jpsi->GetBinContent(i+1);
		fwdPR_Jpsi_Err[i] = h_fwdPR_Jpsi->GetBinError(i+1);
		x_fwd[i] = (ptBin_fwd[i+1]+ptBin_fwd[i])/2;
		binWidth_fwd[i] = (ptBin_fwd[i+1]-ptBin_fwd[i])/2;
		DR_fwdPR[i] = fwdPR_2S[i]/fwdPR_Jpsi[i];
		DR_fwdPR_Err[i] = TMath::Sqrt(TMath::Power(fwdPR_2S_Err[i],2)+TMath::Power(fwdPR_Jpsi_Err[i],2));
		//SysPR_fwd[i] = fwdPR_2S[i]*(hSys_fwd_PR
	}

	TGraphErrors *g_midPR;
	TGraphErrors *g_fwdPR;

	g_midPR = new TGraphErrors(nPtBins_mid, x_mid, DR_midPR, binWidth_mid, DR_midPR_Err);
	g_fwdPR = new TGraphErrors(nPtBins_fwd, x_fwd, DR_fwdPR, binWidth_fwd, DR_fwdPR_Err);

	int iPeriod = 101;
	int iPos = 33;

	double text_x = 0.18;
	double text_y = 0.8;
	double y_diff = 0.08;
	float pos_x = 0.21;
	float pos_y = 0.87;
	float pos_y_diff = 0.051;
	int text_color = 1;
	float text_size = 25;

	TCanvas *c1 = new TCanvas("c1","",900,800);
	c1->cd();

	g_midPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midPR->GetXaxis()->CenterTitle();
	g_midPR->GetYaxis()->SetTitle("Double Ratio");
	g_midPR->GetYaxis()->CenterTitle();
	g_midPR->SetTitle();
	g_midPR->GetXaxis()->SetLimits(0.,40);
	g_midPR->SetMinimum(0.);
	g_midPR->SetMaximum(1.44);

	g_midPR->SetMarkerColor(kBlue+2);
	g_midPR->SetLineColor(kBlue+2);
	g_midPR->SetMarkerStyle(20);
	g_midPR->SetMarkerSize(1.4);

	g_fwdPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_fwdPR->GetXaxis()->CenterTitle();
	g_fwdPR->GetYaxis()->SetTitle("Double Ratio");
	g_fwdPR->GetYaxis()->CenterTitle();
	g_fwdPR->SetTitle();
	g_fwdPR->GetXaxis()->SetLimits(0.,40);
	g_fwdPR->SetMinimum(0.);
	g_fwdPR->SetMaximum(1.44);

	g_fwdPR->SetMarkerColor(kRed+2);
	g_fwdPR->SetLineColor(kRed+2);
	g_fwdPR->SetMarkerStyle(20);
	g_fwdPR->SetMarkerSize(1.4);

	g_midPR->Draw("AP");
	g_fwdPR->Draw("P");

	TLegend *leg1 = new TLegend(0.42,0.84,0.69,0.71);
    leg1->SetTextSize(text_size);
    leg1->SetTextFont(43);
    leg1->SetBorderSize(0);
	leg1->AddEntry(g_midPR,"#bf{Prompt}, |y| < 1.6,  Cent. 0-90%");
	leg1->AddEntry(g_fwdPR,"#bf{Prompt}, 1.6 < |y| < 2.4, Cent. 0-90%");
	leg1->Draw("SAME");
	jumSun(0,1,40,1);
	CMS_lumi_v2mass(c1, iPeriod, iPos);
	c1->SaveAs("figs/DoubleRatio_Prompt_pT.pdf");
    c1->SaveAs("figs/DoubleRatio_Prompt_pT.png");
}
