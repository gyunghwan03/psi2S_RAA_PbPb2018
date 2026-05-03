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

void DoubleRatio_Charmonia_pT_midRap(bool isSys=true)
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid_2S = new TFile(Form("roots/RAA_psi2S_midRap_pT_%.f.root",40.));
	TFile *f_ALICE = new TFile("roots/DoubleRatio_ALICE_pT.root");
	TFile *f_Old = new TFile("roots/DoubleRatio_pT_midRap_HIN_16_004.root");

	TFile *f_Sys1S = new TFile("../syst_summary_Jpsi/syst_roots/total_syst.root");
    TFile *f_Sys2S = new TFile("../syst_summary_psi2S/syst_roots/total_syst.root");

	TFile *f_mid_Jpsi = new TFile("./roots/RAA_JPsi_midRap_pT.root");

	TH1D *h_midPR_2S = (TH1D*) f_mid_2S->Get("hRAA_PR");
	TH1D *h_midPR_Jpsi = (TH1D*) f_mid_Jpsi->Get("hRAA_PR");
	TH1D *h_midNP_2S = (TH1D*) f_mid_2S->Get("hRAA_NP");
	TH1D *h_midNP_Jpsi = (TH1D*) f_mid_Jpsi->Get("hRAA_NP");
	TH1D *h_midPR_Sys1S = (TH1D*) f_Sys1S->Get("mid_pt_PR");
    TH1D *h_midPR_Sys2S = (TH1D*) f_Sys2S->Get("mid_pt_PR");
    TH1D *h_midNP_Sys1S = (TH1D*) f_Sys1S->Get("mid_pt_PR");
    TH1D *h_midNP_Sys2S = (TH1D*) f_Sys2S->Get("mid_pt_PR");

	const int nPtBins_mid = 6;
	double ptBin_mid[nPtBins_mid+1] = {6.5,9,12,15,20,25,40};

	double x_mid[nPtBins_mid]; double binWidth_mid[nPtBins_mid];

	double midPR_2S[nPtBins_mid];
	double midPR_2S_Err[nPtBins_mid];
	double midPR_2S_Sys[nPtBins_mid];

	double midPR_Jpsi[nPtBins_mid];
	double midPR_Jpsi_Err[nPtBins_mid];
	double midPR_Jpsi_Sys[nPtBins_mid];

	double DR_midPR[nPtBins_mid];
	double DR_midPR_Err[nPtBins_mid];
	double DR_midPR_Sys[nPtBins_mid];

	double midNP_2S[nPtBins_mid];
	double midNP_2S_Err[nPtBins_mid];
	double midNP_2S_Sys[nPtBins_mid];

	double midNP_Jpsi[nPtBins_mid];
	double midNP_Jpsi_Err[nPtBins_mid];
	double midNP_Jpsi_Sys[nPtBins_mid];

	double DR_midNP[nPtBins_mid];
	double DR_midNP_Err[nPtBins_mid];
	double DR_midNP_Sys[nPtBins_mid];

	for (int i=0; i<nPtBins_mid; i++){
		midPR_2S[i] = h_midPR_2S->GetBinContent(i+1);
		midPR_2S_Err[i] = h_midPR_2S->GetBinError(i+1);
		midPR_Jpsi[i] = h_midPR_Jpsi->GetBinContent(i+1);
		midPR_Jpsi_Err[i] = h_midPR_Jpsi->GetBinError(i+1);
		midPR_Jpsi_Sys[i] = midPR_Jpsi[i]*(h_midPR_Sys1S->GetBinContent(i+1));
        midPR_2S_Sys[i] = midPR_2S[i]*(h_midPR_Sys2S->GetBinContent(i+1));

		midNP_2S[i] = h_midNP_2S->GetBinContent(i+1);
		midNP_2S_Err[i] = h_midNP_2S->GetBinError(i+1);
		midNP_Jpsi[i] = h_midNP_Jpsi->GetBinContent(i+1);
		midNP_Jpsi_Err[i] = h_midNP_Jpsi->GetBinError(i+1);
		midNP_Jpsi_Sys[i] = midNP_Jpsi[i]*(h_midNP_Sys1S->GetBinContent(i+1));
        midNP_2S_Sys[i] = midNP_2S[i]*(h_midNP_Sys2S->GetBinContent(i+1));

		x_mid[i] = (ptBin_mid[i+1]+ptBin_mid[i])/2;
		binWidth_mid[i] = (ptBin_mid[i+1]-ptBin_mid[i])/2;

		DR_midPR[i] = midPR_2S[i]/midPR_Jpsi[i];
		DR_midPR_Err[i] = DR_midPR[i]*(TMath::Sqrt(TMath::Power(midPR_2S_Err[i]/midPR_2S[i],2)+TMath::Power(midPR_Jpsi_Err[i]/midPR_Jpsi[i],2)));
		DR_midPR_Sys[i] = DR_midPR[i]*(TMath::Sqrt(TMath::Power(midPR_Jpsi_Sys[i]/midPR_Jpsi[i],2)+TMath::Power(midPR_2S_Sys[i]/midPR_2S[i],2)));
		DR_midNP[i] = midNP_2S[i]/midNP_Jpsi[i];
		DR_midNP_Err[i] = DR_midNP[i]*(TMath::Sqrt(TMath::Power(midNP_2S_Err[i]/midNP_2S[i],2)+TMath::Power(midNP_Jpsi_Err[i]/midNP_Jpsi[i],2)));
		DR_midNP_Sys[i] = DR_midNP[i]*(TMath::Sqrt(TMath::Power(midNP_Jpsi_Sys[i]/midNP_Jpsi[i],2)+TMath::Power(midNP_2S_Sys[i]/midNP_2S[i],2)));

		//SysPR[i] = midPR_2S[i]*(hSys_PR->GetBinContent(i+1));
	}


	TGraphErrors *g_midPR; TGraphErrors *g_midPR_Sys;
	TGraphErrors *g_midNP; TGraphErrors *g_midNP_Sys;
	TGraphAsymmErrors *g_alice = (TGraphAsymmErrors*)f_ALICE->Get("Table 6/Graph1D_y1");
	TGraphAsymmErrors *g_alice_Sys = (TGraphAsymmErrors*)f_ALICE->Get("Table 6/Graph1D_y1");
	TGraphAsymmErrors *g_old = (TGraphAsymmErrors*)f_Old->Get("Table 1/Graph1D_y1");
	TGraphAsymmErrors *g_old_Sys = (TGraphAsymmErrors*)f_Old->Get("Table 1/Graph1D_y1");

	TH1D *h_old_Err= (TH1D*)f_Old->Get("Table 1/Hist1D_y1_e1");
	TH1D *h_old_Sys= (TH1D*)f_Old->Get("Table 1/Hist1D_y1_e2");
	double old_sys[5]; double binWidth_old[5]; double old_err[5];
	for(int i=0; i<5; i++){
		old_err[i] = h_old_Err->GetBinContent(i+1);
		old_sys[i] = h_old_Sys->GetBinContent(i+1);
		binWidth_old[i] = (h_old_Sys->GetBinWidth(i+1))/2;
		if(isSys==1){
			g_old_Sys->SetPointError(i,binWidth_old[i],binWidth_old[i],old_sys[i],old_sys[i]);
			g_old->SetPointError(i,0,0,old_err[i],old_err[i]); 
		}
	}
	TH1D *h_alice_Err = (TH1D*)f_ALICE->Get("Table 6/Hist1D_y1_e1");
	TH1D *h_alice_Sys = (TH1D*)f_ALICE->Get("Table 6/Hist1D_y1_e2");
	double alice_sys[4]; double binWidth_alice[4]; double alice_err[4];
	for(int i=0; i<4; i++){
		alice_err[i] = h_alice_Err->GetBinContent(i+1);
		alice_sys[i] = h_alice_Sys->GetBinContent(i+1);
		binWidth_alice[i] = (h_alice_Sys->GetBinWidth(i+1))/2;
		if(isSys==1){
			g_alice_Sys->SetPointError(i,binWidth_alice[i],binWidth_alice[i],alice_sys[i],alice_sys[i]);
			g_alice->SetPointError(i,0,0,alice_err[i],alice_err[i]);
		}
	}

	g_midPR = new TGraphErrors(nPtBins_mid, x_mid, DR_midPR, binWidth_mid, DR_midPR_Err);
	g_midNP = new TGraphErrors(nPtBins_mid, x_mid, DR_midNP, binWidth_mid, DR_midNP_Err);
	g_midPR_Sys = new TGraphErrors(nPtBins_mid, x_mid, DR_midPR, binWidth_mid, DR_midPR_Sys);
	g_midNP_Sys = new TGraphErrors(nPtBins_mid, x_mid, DR_midNP, binWidth_mid, DR_midNP_Sys);

	if(isSys==1){
		g_midPR = new TGraphErrors(nPtBins_mid, x_mid, DR_midPR, 0, DR_midPR_Err);
		g_midNP = new TGraphErrors(nPtBins_mid, x_mid, DR_midNP, 0, DR_midNP_Err);
	}



	int iPeriod = 101;
	int iPos = 33;

	double text_x = 0.18;
	double text_y = 0.8;
	double y_diff = 0.08;
	float pos_x = 0.21;
	float pos_y = 0.87;
	float pos_y_diff = 0.051;
	int text_color = 1;
	float text_size = 20;

	TCanvas *c1 = new TCanvas("c1","",900,800);
	c1->cd();

	g_midPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midPR->GetXaxis()->CenterTitle();
	g_midPR->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_midPR->GetYaxis()->CenterTitle();
	g_midPR->SetTitle();
	g_midPR->GetXaxis()->SetLimits(0.,40);
	g_midPR->SetMinimum(0.);
	g_midPR->SetMaximum(1.44);

	g_midPR->SetMarkerColor(kBlue+2);
	g_midPR->SetLineColor(kBlue+2);
	g_midPR->SetLineWidth(2);
	g_midPR->SetMarkerStyle(20);
	g_midPR->SetMarkerSize(1.6);
	g_midPR_Sys->SetLineColor(kBlue-4);
    g_midPR_Sys->SetFillColorAlpha(kBlue-9,0.40);

	g_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    g_old->GetXaxis()->CenterTitle();
    g_old->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
    g_old->GetYaxis()->CenterTitle();
    g_old->SetTitle();
    g_old->GetXaxis()->SetLimits(0.,40);
    g_old->SetMinimum(0.);
    g_old->SetMaximum(1.44);

    g_old->SetMarkerColor(15);
    g_old->SetLineColor(15);
    g_old->SetMarkerStyle(30);
    g_old->SetMarkerSize(2.0);
	g_old_Sys->SetLineColor(15);
    g_old_Sys->SetFillColorAlpha(15,0.40);

	g_alice->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_alice->GetXaxis()->CenterTitle();
	g_alice->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_alice->GetYaxis()->CenterTitle();
	g_alice->SetTitle();
	g_alice->GetXaxis()->SetLimits(0.,40);
	g_alice->SetMinimum(0.);
	g_alice->SetMaximum(1.44);

	g_alice->SetMarkerColor(45);
    g_alice->SetLineColor(45);
    g_alice->SetMarkerStyle(33);
    g_alice->SetMarkerSize(2.0);
	g_alice_Sys->SetLineColor(45);
	g_alice_Sys->SetFillColorAlpha(45,0.40);

	g_midPR->Draw("AP");
	g_old->Draw("P");
	g_alice->Draw("P");
	if(isSys==1) {
		g_midPR_Sys->Draw("5");
		g_old_Sys->Draw("5");
		g_alice_Sys->Draw("5");
	}

	TLegend *leg1 = new TLegend(0.17,0.87,0.44,0.73);
    leg1->SetTextSize(text_size);
    leg1->SetTextFont(43);
    leg1->SetBorderSize(0);
	leg1->AddEntry(g_midPR,"#bf{Prompt}, |y| < 1.6,  Cent. 0-90%", "pe");
	leg1->AddEntry(g_old,"#bf{Prompt} (HIN-16-004), |y| < 1.6, Cent. 0-100%", "pe");
	leg1->AddEntry(g_alice,"#bf{Inclusive} (ALICE), 2.5 < y < 4, Cent. 0-90%", "pe");
	leg1->Draw("SAME");
	jumSun(0,1,40,1);
	CMS_lumi_v2mass(c1, iPeriod, iPos);
	c1->SaveAs("figs/DoubleRatio_Prompt_midRap_pT.pdf");
    c1->SaveAs("figs/DoubleRatio_Prompt_midRap_pT.png");


	TCanvas *c2 = new TCanvas("c2","",900,800);
	c2->cd();

	g_midNP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midNP->GetXaxis()->CenterTitle();
	g_midNP->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_midNP->GetYaxis()->CenterTitle();
	g_midNP->SetTitle();
	g_midNP->GetXaxis()->SetLimits(0.,40);
	g_midNP->SetMinimum(0.);
	g_midNP->SetMaximum(1.44);

	g_midNP->SetMarkerColor(kRed+3);
	g_midNP->SetLineColor(kRed+3);
	g_midNP->SetLineWidth(2);
	g_midNP->SetMarkerStyle(21);
	g_midNP->SetMarkerSize(1.6);
	g_midNP_Sys->SetLineColor(kRed-2);
	g_midNP_Sys->SetFillColorAlpha(kRed-9,0.40);

	g_midNP->Draw("AP");
	g_alice->Draw("P");
	if(isSys) {
		g_midNP_Sys->Draw("5");
		g_alice_Sys->Draw("5");
	}

	TLegend *leg2 = new TLegend(0.17,0.87,0.44,0.73);
    leg2->SetTextSize(text_size);
    leg2->SetTextFont(43);
    leg2->SetBorderSize(0);
	leg2->AddEntry(g_midNP,"#bf{NonPrompt}, |y| < 1.6,  Cent. 0-90%", "pe");
	leg2->AddEntry(g_alice,"#bf{Inclusive} (ALICE), 2.5 < y < 4, Cent. 0-90%", "pe");
	leg2->Draw("SAME");
	jumSun(0,1,40,1);
	CMS_lumi_v2mass(c2, iPeriod, iPos);
	c2->SaveAs("figs/DoubleRatio_NonPrompt_midRap_pT.pdf");
    c2->SaveAs("figs/DoubleRatio_NonPrompt_midRap_pT.png");
}
