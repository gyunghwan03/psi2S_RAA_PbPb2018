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

void DoubleRatio_Charmonia_Npart()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid_2S = new TFile("roots/RAA_psi2S_midRap_Npart_pt6p5_40_8Bins.root");
	TFile *f_fwd_2S = new TFile("roots/RAA_psi2S_forRap_Npart_pt3p5_40_4Bins.root");

	TFile *f_mid_Jpsi = new TFile("./roots/RAA_JPsi_midRap_Npart_8Bins.root");
	TFile *f_fwd_Jpsi = new TFile("./roots/RAA_JPsi_forRap_Npart_4Bins.root");

	TFile *f_Old = new TFile("./roots/DoubleRatio_midRap_HIN_16_004.root");

	TH1D *h_midPR_2S = (TH1D*) f_mid_2S->Get("hRAA_PR");
	TH1D *h_fwdPR_2S = (TH1D*) f_fwd_2S->Get("hRAA_PR");
	TH1D *h_midPR_Jpsi = (TH1D*) f_mid_Jpsi->Get("hRAA_PR");
	TH1D *h_fwdPR_Jpsi = (TH1D*) f_fwd_Jpsi->Get("hRAA_PR");

	TH1D *h_Old = (TH1D*) f_Old->Get("Table 5/Hist1D_y1");
	TH1D *h_Old_Err = (TH1D*) f_Old->Get("Table 5/Hist1D_y1_e1");
	TH1D *h_OldNpart = (TH1D*) f_Old->Get("Table 5/Hist1D_y2");

	const int nCentBins_mid = 8;
	const int nCentBins_fwd = 4;
	const int nCentBins_old = 6;
	double NpartBin_mid[nCentBins_mid+1] = {27.12,87.19,131.0,188.2,262.3,283.6,331.5,382.3};
	double NpartBin_fwd[nCentBins_fwd+1] = {27.12,109.1,225.2,356.9};
	double NpartBin_old[nCentBins_old];

	double x_mid[nCentBins_mid];
	double x_fwd[nCentBins_fwd];
	double binWidth_mid[nCentBins_mid]={4.3,4.3,4.3,4.3,4.3,4.3,4.3,4.3};
	double binWidth_old[nCentBins_old]={4.3,4.3,4.3,4.3,4.3,4.3};
	double binWidth_fwd[nCentBins_fwd]={4.3,4.3,4.3,4.3};

	double midPR_2S[nCentBins_mid];
	double midPR_2S_Err[nCentBins_mid];
	double midPR_2S_Sys[nCentBins_mid];
	double fwdPR_2S[nCentBins_fwd];
	double fwdPR_2S_Err[nCentBins_fwd];
	double fwdPR_2S_Sys[nCentBins_fwd];

	double midPR_Jpsi[nCentBins_mid];
	double midPR_Jpsi_Err[nCentBins_mid];
	double midPR_Jpsi_Sys[nCentBins_mid];
	double fwdPR_Jpsi[nCentBins_fwd];
	double fwdPR_Jpsi_Err[nCentBins_fwd];
	double fwdPR_Jpsi_Sys[nCentBins_fwd];

	double DR_midPR[nCentBins_mid];
	double DR_fwdPR[nCentBins_fwd];
	double DR_midPR_Err[nCentBins_mid];
	double DR_fwdPR_Err[nCentBins_fwd];

	double DR_old[nCentBins_old];
	double DR_old_Err[nCentBins_old];

	for (int i=0; i<nCentBins_mid; i++){
		midPR_2S[i] = h_midPR_2S->GetBinContent(i+1);
		midPR_2S_Err[i] = h_midPR_2S->GetBinError(i+1);
		midPR_Jpsi[i] = h_midPR_Jpsi->GetBinContent(i+1);
		midPR_Jpsi_Err[i] = h_midPR_Jpsi->GetBinError(i+1);
		DR_midPR[i] = midPR_2S[i]/midPR_Jpsi[i];
		DR_midPR_Err[i] = TMath::Sqrt(TMath::Power(midPR_2S_Err[i],2)+TMath::Power(midPR_Jpsi_Err[i],2));

		//SysPR[i] = midPR_2S[i]*(hSys_PR->GetBinContent(i+1));
	}

	for (int i=0; i<nCentBins_fwd; i++){
		fwdPR_2S[i] = h_fwdPR_2S->GetBinContent(i+1);
		fwdPR_2S_Err[i] = h_fwdPR_2S->GetBinError(i+1);
		fwdPR_Jpsi[i] = h_fwdPR_Jpsi->GetBinContent(i+1);
		fwdPR_Jpsi_Err[i] = h_fwdPR_Jpsi->GetBinError(i+1);
		DR_fwdPR[i] = fwdPR_2S[i]/fwdPR_Jpsi[i];
		DR_fwdPR_Err[i] = TMath::Sqrt(TMath::Power(fwdPR_2S_Err[i],2)+TMath::Power(fwdPR_Jpsi_Err[i],2));
		//SysPR_fwd[i] = fwdPR_2S[i]*(hSys_fwd_PR
	}

	for (int i=nCentBins_old-1; i>=0;i--){
		DR_old[nCentBins_old-i-1] = h_Old->GetBinContent(i+1);
		DR_old_Err[nCentBins_old-i-1] = h_Old_Err->GetBinContent(i+1);
		NpartBin_old[nCentBins_old-i-1] = h_OldNpart->GetBinContent(i+1);
		cout << "i : "<< nCentBins_old-i-1 << ",	DR :	" << DR_old[nCentBins_old-i-1] << ",	NPart :	" << NpartBin_old[nCentBins_old-i-1] << endl; 
	}

	TGraphErrors *g_midPR;
	TGraphErrors *g_fwdPR;
	TGraphErrors *g_old;

	g_midPR = new TGraphErrors(nCentBins_mid, NpartBin_mid, DR_midPR, 0, DR_midPR_Err);
	g_fwdPR = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, DR_fwdPR, 0, DR_fwdPR_Err);
	g_old = new TGraphErrors(nCentBins_old-1, NpartBin_old, DR_old, 0, DR_old_Err);

	TArrow *old_upper = new TArrow(NpartBin_old[5],0,NpartBin_old[5],0.55,0.027,"<-|");

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

	g_midPR->GetXaxis()->SetTitle("<N_{Part}>");
	g_midPR->GetXaxis()->CenterTitle();
	g_midPR->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_midPR->GetYaxis()->CenterTitle();
	g_midPR->SetTitle();
	g_midPR->GetXaxis()->SetLimits(0.,400);
	g_midPR->SetMinimum(0.);
	g_midPR->SetMaximum(1.44);

	g_midPR->SetMarkerColor(kBlue+2);
	g_midPR->SetLineColor(kBlue+2);
	g_midPR->SetMarkerStyle(20);
	g_midPR->SetMarkerSize(1.4);

	g_fwdPR->GetXaxis()->SetTitle("<N_{Part}>");
	g_fwdPR->GetXaxis()->CenterTitle();
	g_fwdPR->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_fwdPR->GetYaxis()->CenterTitle();
	g_fwdPR->SetTitle();
	g_fwdPR->GetXaxis()->SetLimits(0.,400);
	g_fwdPR->SetMinimum(0.);
	g_fwdPR->SetMaximum(1.44);

	g_fwdPR->SetMarkerColor(kRed+2);
	g_fwdPR->SetLineColor(kRed+2);
	g_fwdPR->SetMarkerStyle(21);
	g_fwdPR->SetMarkerSize(1.4);

	g_old->GetXaxis()->SetTitle("<N_{Part}>");
	g_old->GetXaxis()->CenterTitle();
	g_old->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_old->GetYaxis()->CenterTitle();
	g_old->SetTitle();
	g_old->GetXaxis()->SetLimits(0.,400);
	g_old->SetMinimum(0.);
	g_old->SetMaximum(1.44);

	g_old->SetMarkerColor(15);
	g_old->SetLineColor(15);
	g_old->SetMarkerStyle(30);
	g_old->SetMarkerSize(1.8);

	g_midPR->Draw("AP");
	g_fwdPR->Draw("P");
	g_old->Draw("P");
	old_upper->Draw("");

	TLegend *leg1 = new TLegend(0.17,0.87,0.44,0.77);
    leg1->SetTextSize(text_size);
    leg1->SetTextFont(43);
    leg1->SetBorderSize(0);
	leg1->AddEntry(g_midPR,"#bf{Prompt}, |y| < 1.6,  Cent. 0-90%");
	leg1->AddEntry(g_fwdPR,"#bf{Prompt}, 1.6 < |y| < 2.4, Cent. 0-90%");
	leg1->AddEntry(g_old,"#bf{Prompt(Old)}, 6.5 < p_{T} < 30 GeV/c, |y| <1.6, Cent. 0-100%");

	leg1->Draw("SAME");
	jumSun(0,1,400,1);
	CMS_lumi_v2mass(c1, iPeriod, iPos);
	c1->SaveAs("figs/DoubleRatio_Prompt_Npart.pdf");
	c1->SaveAs("figs/DoubleRatio_Prompt_Npart.png");
}
