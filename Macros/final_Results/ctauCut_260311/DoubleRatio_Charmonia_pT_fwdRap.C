#include "../../../rootFitHeaders.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"
#include "../../../cutsAndBin.h"
#include "../../../CMS_lumi_v2mass.C"
#include "../../../tdrstyle.C"
#include "../../../Style.h"
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

void DoubleRatio_Charmonia_pT_fwdRap(bool isSys=true)
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_fwd_2S = new TFile(Form("roots/RAA_psi2S_forRap_pT_%.f.root",40.));
	TFile *f_ALICE = new TFile("../roots/DoubleRatio_ALICE_pT.root");
	TFile *f_Old = new TFile("../roots/DoubleRatio_pT_fwdRap_HIN_16_004.root");

	TFile *f_Sys1S = new TFile("../../syst_summary_Jpsi/syst_roots/DoubleRatio_syst.root");
    TFile *f_Sys2S = new TFile("../../syst_summary_psi2S/syst_roots/DoubleRatio_syst.root");

	TFile *f_fwd_Jpsi = new TFile("roots/RAA_JPsi_forRap_pT.root");

	TH1D *h_fwdPR_2S = (TH1D*) f_fwd_2S->Get("hRAA_PR");
	TH1D *h_fwdPR_Jpsi = (TH1D*) f_fwd_Jpsi->Get("hRAA_PR");
	TH1D *h_old = (TH1D*) f_Old->Get("Table 3/Hist1D_y1");
	TH1D *h_oldErr = (TH1D*) f_Old->Get("Table 3/Hist1D_y1_e1");
	TH1D *h_oldSys = (TH1D*) f_Old->Get("Table 3/Hist1D_y1_e2");
	TH1D *h_fwdNP_2S = (TH1D*) f_fwd_2S->Get("hRAA_NP");
	TH1D *h_fwdNP_Jpsi = (TH1D*) f_fwd_Jpsi->Get("hRAA_NP");
	TH1D *h_fwdPR_Sys1S = (TH1D*) f_Sys1S->Get("fwd_pt_PR");
    TH1D *h_fwdPR_Sys2S = (TH1D*) f_Sys2S->Get("fwd_pt_PR");
    TH1D *h_fwdNP_Sys1S = (TH1D*) f_Sys1S->Get("fwd_pt_PR");
    TH1D *h_fwdNP_Sys2S = (TH1D*) f_Sys2S->Get("fwd_pt_PR");

	const int nPtBins_fwd = 4;
	double ptBin_fwd[nPtBins_fwd+1] = {3.5,6.5,9,12,40};

	double x_fwd[nPtBins_fwd]; double binWidth_fwd[nPtBins_fwd];

	double fwdPR_2S[nPtBins_fwd];
	double fwdPR_2S_Err[nPtBins_fwd];
	double fwdPR_2S_Sys[nPtBins_fwd];

	double fwdPR_Jpsi[nPtBins_fwd];
	double fwdPR_Jpsi_Err[nPtBins_fwd];
	double fwdPR_Jpsi_Sys[nPtBins_fwd];

	double DR_fwdPR[nPtBins_fwd];
	double DR_fwdPR_Err[nPtBins_fwd];
	double DR_fwdPR_Sys[nPtBins_fwd];

	double fwdNP_2S[nPtBins_fwd];
	double fwdNP_2S_Err[nPtBins_fwd];
	double fwdNP_2S_Sys[nPtBins_fwd];

	double fwdNP_Jpsi[nPtBins_fwd];
	double fwdNP_Jpsi_Err[nPtBins_fwd];
	double fwdNP_Jpsi_Sys[nPtBins_fwd];

	double DR_fwdNP[nPtBins_fwd];
	double DR_fwdNP_Err[nPtBins_fwd];
	double DR_fwdNP_Sys[nPtBins_fwd];

	cout << "/////////////// Forward Rapdiity ///////////////" << endl; 
	for (int i=0; i<nPtBins_fwd; i++){
		fwdPR_2S[i] = h_fwdPR_2S->GetBinContent(i+1);
		fwdPR_2S_Err[i] = h_fwdPR_2S->GetBinError(i+1);
		fwdPR_Jpsi[i] = h_fwdPR_Jpsi->GetBinContent(i+1);
		fwdPR_Jpsi_Err[i] = h_fwdPR_Jpsi->GetBinError(i+1);
		fwdPR_Jpsi_Sys[i] = fwdPR_Jpsi[i]*(h_fwdPR_Sys1S->GetBinContent(i+1));
        fwdPR_2S_Sys[i] = fwdPR_2S[i]*(h_fwdPR_Sys2S->GetBinContent(i+1));

		fwdNP_2S[i] = h_fwdNP_2S->GetBinContent(i+1);
		fwdNP_2S_Err[i] = h_fwdNP_2S->GetBinError(i+1);
		fwdNP_Jpsi[i] = h_fwdNP_Jpsi->GetBinContent(i+1);
		fwdNP_Jpsi_Err[i] = h_fwdNP_Jpsi->GetBinError(i+1);
		fwdNP_Jpsi_Sys[i] = fwdNP_Jpsi[i]*(h_fwdNP_Sys1S->GetBinContent(i+1));
		fwdNP_2S_Sys[i] = fwdNP_2S[i]*(h_fwdNP_Sys2S->GetBinContent(i+1));

		x_fwd[i] = (ptBin_fwd[i+1]+ptBin_fwd[i])/2;
		binWidth_fwd[i] = (ptBin_fwd[i+1]-ptBin_fwd[i])/2;
		DR_fwdPR[i] = fwdPR_2S[i]/fwdPR_Jpsi[i];
		DR_fwdPR_Err[i] = DR_fwdPR[i]*(TMath::Sqrt(TMath::Power(fwdPR_2S_Err[i]/fwdPR_2S[i],2)+TMath::Power(fwdPR_Jpsi_Err[i]/fwdPR_Jpsi[i],2)));
		DR_fwdPR_Sys[i] = DR_fwdPR[i]*(TMath::Sqrt(TMath::Power(fwdPR_Jpsi_Sys[i]/fwdPR_Jpsi[i],2)+TMath::Power(fwdPR_2S_Sys[i]/fwdPR_2S[i],2)));
		DR_fwdNP[i] = fwdNP_2S[i]/fwdNP_Jpsi[i];
		DR_fwdNP_Err[i] = DR_fwdNP[i]*(TMath::Sqrt(TMath::Power(fwdNP_2S_Err[i]/fwdNP_2S[i],2)+TMath::Power(fwdNP_Jpsi_Err[i]/fwdNP_Jpsi[i],2)));
		DR_fwdNP_Sys[i] = DR_fwdNP[i]*(TMath::Sqrt(TMath::Power(fwdNP_Jpsi_Sys[i]/fwdNP_Jpsi[i],2)+TMath::Power(fwdNP_2S_Sys[i]/fwdNP_2S[i],2)));
		cout << "ptBin	" << ptBin_fwd[i] << " - " << ptBin_fwd[i+1] << "	psi2S RAA :	" << fwdPR_2S[i] << "	Jpsi RAA :	" << fwdPR_Jpsi[i] << endl; 
		//SysPR_fwd[i] = fwdPR_2S[i]*(hSys_fwd_PR
	}

	TGraphErrors *g_fwdPR; TGraphErrors *g_fwdPR_Sys;
	TGraphErrors *g_fwdNP; TGraphErrors *g_fwdNP_Sys;
	TGraphAsymmErrors *g_alice = (TGraphAsymmErrors*)f_ALICE->Get("Table 6/Graph1D_y1");
	TGraphAsymmErrors *g_alice_Sys = (TGraphAsymmErrors*)f_ALICE->Get("Table 6/Graph1D_y1");
	TGraphErrors *g_old; TGraphErrors *g_old_Sys;

	g_fwdPR = new TGraphErrors(nPtBins_fwd, x_fwd, DR_fwdPR, binWidth_fwd, DR_fwdPR_Err);
	g_fwdNP = new TGraphErrors(nPtBins_fwd, x_fwd, DR_fwdNP, binWidth_fwd, DR_fwdNP_Err);
	g_fwdPR_Sys = new TGraphErrors(nPtBins_fwd, x_fwd, DR_fwdPR, binWidth_fwd, DR_fwdPR_Sys);
	g_fwdNP_Sys = new TGraphErrors(nPtBins_fwd, x_fwd, DR_fwdNP, binWidth_fwd, DR_fwdNP_Sys);
	if(isSys){
		g_fwdPR = new TGraphErrors(nPtBins_fwd, x_fwd, DR_fwdPR, 0, DR_fwdPR_Err);
		g_fwdNP = new TGraphErrors(nPtBins_fwd, x_fwd, DR_fwdNP, 0, DR_fwdNP_Err);
	}
	g_old = new TGraphErrors(); g_old_Sys = new TGraphErrors();
	
	g_old->SetPoint(0,h_old->GetBinCenter(2),h_old->GetBinContent(2));
	g_old->SetPoint(1,h_old->GetBinCenter(3),h_old->GetBinContent(3));
	g_old->SetPointError(0,h_old->GetBinWidth(2)/2,h_oldErr->GetBinContent(2));
	g_old->SetPointError(1,h_old->GetBinWidth(3)/2,h_oldErr->GetBinContent(3));
	if(isSys) {
		g_old->SetPointError(0,0,h_oldErr->GetBinContent(2));
		g_old->SetPointError(1,0,h_oldErr->GetBinContent(3));
	}
	g_old_Sys->SetPoint(0,h_old->GetBinCenter(2),h_old->GetBinContent(2));
	g_old_Sys->SetPoint(1,h_old->GetBinCenter(3),h_old->GetBinContent(3));
	g_old_Sys->SetPointError(0,h_old->GetBinWidth(2)/2,h_oldSys->GetBinContent(2));
	g_old_Sys->SetPointError(1,h_old->GetBinWidth(3)/2,h_oldSys->GetBinContent(3));

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

	TArrow *old_upper = new TArrow(4.75, 0, 4.75, 0.84, 0.027,"<-|");

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

	g_fwdPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    g_fwdPR->GetXaxis()->CenterTitle();
    g_fwdPR->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
    g_fwdPR->GetYaxis()->CenterTitle();
    g_fwdPR->SetTitle();
    g_fwdPR->GetXaxis()->SetLimits(0.,40);
    g_fwdPR->SetMinimum(0.);
    g_fwdPR->SetMaximum(1.44);
	g_fwdPR->SetMarkerColor(kBlue+2);
	g_fwdPR->SetLineColor(kBlue+2);
	g_fwdPR->SetLineWidth(2);
	g_fwdPR->SetMarkerStyle(20);
	g_fwdPR->SetMarkerSize(1.6);
	cout << "HERE" << endl;
	g_fwdPR_Sys->SetLineColor(kBlue-4);
	g_fwdPR_Sys->SetFillColorAlpha(kBlue-9,0.40);

	g_fwdNP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    g_fwdNP->GetXaxis()->CenterTitle();
    g_fwdNP->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
    g_fwdNP->GetYaxis()->CenterTitle();
    g_fwdNP->SetTitle();
    g_fwdNP->GetXaxis()->SetLimits(0.,40);
    g_fwdNP->SetMinimum(0.);
    g_fwdNP->SetMaximum(1.44);
	g_fwdNP->SetMarkerColor(kRed+3);
	g_fwdNP->SetLineColor(kRed+3);
	g_fwdNP->SetLineWidth(2);
	g_fwdNP->SetMarkerStyle(21);
	g_fwdNP->SetMarkerSize(1.6);
	g_fwdNP_Sys->SetLineColor(kRed-4);
	g_fwdNP_Sys->SetFillColorAlpha(kRed-9,0.40);

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
	g_alice_Sys->SetFillColorAlpha(45,0.4);

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

	TCanvas *c1 = new TCanvas("c1","",900,800);
	c1->cd();

	g_fwdPR->Draw("AP");
	g_alice->Draw("P");
	g_old->Draw("P");
	if(isSys){
		g_fwdPR_Sys->Draw("5");
		g_alice_Sys->Draw("5");
		g_old_Sys->Draw("5");
	}
	old_upper->Draw("");

	TLegend *leg1 = new TLegend(0.17,0.87,0.44,0.73);
    leg1->SetTextSize(text_size);
    leg1->SetTextFont(43);
    leg1->SetBorderSize(0);
	leg1->AddEntry(g_fwdPR,"#bf{Prompt}, 1.6 < |y| < 2.4, Cent. 0-90%", "pe");
	leg1->AddEntry(g_old, "#bf{Prompt} (HIN-16-004), 1.6 < |y| < 2.4, Cent. 0-100%", "pe");
	leg1->AddEntry(g_alice,"#bf{Inclusive} (ALICE), 2.5 < y < 4, Cent. 0-90%", "pe");
	leg1->Draw("SAME");
	jumSun(0,1,40,1);
	CMS_lumi_v2mass(c1, iPeriod, iPos);
	c1->SaveAs("figs/DoubleRatio_Prompt_fwdRap_pT.pdf");
    c1->SaveAs("figs/DoubleRatio_Prompt_fwdRap_pT.png");

	TCanvas *c2 = new TCanvas("c2","",900,800);
	c2->cd();

	g_fwdNP->Draw("AP");
	g_alice->Draw("P");
	if(isSys){
		g_fwdNP_Sys->Draw("5");
		g_alice_Sys->Draw("5");
	}

	TLegend *leg2 = new TLegend(0.17,0.87,0.44,0.73);
    leg2->SetTextSize(text_size);
    leg2->SetTextFont(43);
    leg2->SetBorderSize(0);
	leg2->AddEntry(g_fwdNP,"#bf{NonPrompt}, 1.6 < |y| < 2.4, Cent. 0-90%", "pe");
	leg2->AddEntry(g_alice,"#bf{Inclusive} (ALICE), 2.5 < y < 4, Cent. 0-90%", "pe");
	leg2->Draw("SAME");
	jumSun(0,1,40,1);
	CMS_lumi_v2mass(c2, iPeriod, iPos);
	c2->SaveAs("figs/DoubleRatio_NonPrompt_fwdRap_pT.pdf");
    c2->SaveAs("figs/DoubleRatio_NonPrompt_fwdRap_pT.png");
}
