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

void DoubleRatio_Charmonia_Npart_midRap_240807()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid_2S = new TFile("roots/RAA_psi2S_midRap_Npart_pt6p5_40.root");

	TFile *f_Old = new TFile("../roots/DoubleRatio_midRap_HIN_16_004.root");
	TFile *f_ATLAS = new TFile("../roots/DoubleRatio_ATLAS.root");
	TFile *f_ATLAS_NP = new TFile("../roots/DoubleRatio_ATLAS_NonPrompt.root");

	TFile *f_mid_Jpsi = new TFile("roots/RAA_JPsi_midRap_Npart.root");

	TFile *f_Sys1S = new TFile("../../syst_summary_Jpsi/syst_roots/DoubleRatio_syst.root");
	TFile *f_Sys2S = new TFile("../../syst_summary_psi2S/syst_roots/DoubleRatio_syst.root");

	TH1D *h_midPR_2S = (TH1D*) f_mid_2S->Get("hRAA_PR");
	TH1D *h_midPR_Jpsi = (TH1D*) f_mid_Jpsi->Get("hRAA_PR");
	TH1D *h_midNP_2S = (TH1D*) f_mid_2S->Get("hRAA_NP");
	TH1D *h_midNP_Jpsi = (TH1D*) f_mid_Jpsi->Get("hRAA_NP");
	TH1D *h_midPR_Sys1S = (TH1D*) f_Sys1S->Get("mid_cent_PR");
	TH1D *h_midPR_Sys2S = (TH1D*) f_Sys2S->Get("mid_cent_PR");
	TH1D *h_midNP_Sys1S = (TH1D*) f_Sys1S->Get("mid_cent_PR");
	TH1D *h_midNP_Sys2S = (TH1D*) f_Sys2S->Get("mid_cent_PR");

	TH1D *h_Old = (TH1D*) f_Old->Get("Table 5/Hist1D_y1");
    TH1D *h_Old_Err = (TH1D*) f_Old->Get("Table 5/Hist1D_y1_e1");
    TH1D *h_Old_Sys = (TH1D*) f_Old->Get("Table 5/Hist1D_y1_e2");
    TH1D *h_OldNpart = (TH1D*) f_Old->Get("Table 5/Hist1D_y2");

	const int nCentBins_mid = 6;
	const int nCentBins_fwd = 4;
	const int nCentBins_old = 6;
	double NpartBin_mid[nCentBins_mid+1] = {27.12,87.19,131.0,188.2,262.3,356.9};
	double NpartBin_fwd[nCentBins_fwd+1] = {27.12,109.1,225.2,356.9};
	double NpartBin_old[nCentBins_old];

	double x_mid[nCentBins_mid];
	double x_fwd[nCentBins_fwd];
	double binWidth_mid[nCentBins_mid]={4.3,4.3,4.3,4.3,4.3,4.3};
	double binWidth_fwd[nCentBins_fwd]={4.3,4.3,4.3,4.3};
	double binWidth_old[nCentBins_old]={4.3,4.3,4.3,4.3,4.3,4.3};

	double midPR_2S[nCentBins_mid];
	double midPR_2S_Err[nCentBins_mid];
	double midPR_2S_Sys[nCentBins_mid];

	double midPR_Jpsi[nCentBins_mid];
	double midPR_Jpsi_Err[nCentBins_mid];
	double midPR_Jpsi_Sys[nCentBins_mid];

	double DR_midPR[nCentBins_mid];
	double DR_midPR_Err[nCentBins_mid];

	double midNP_2S[nCentBins_mid];
	double midNP_2S_Err[nCentBins_mid];
	double midNP_2S_Sys[nCentBins_mid];

	double midNP_Jpsi[nCentBins_mid];
	double midNP_Jpsi_Err[nCentBins_mid];
	double midNP_Jpsi_Sys[nCentBins_mid];

	double DR_midNP[nCentBins_mid];
	double DR_midNP_Err[nCentBins_mid];

	double DR_old[nCentBins_old];
    double DR_old_Err[nCentBins_old];
    double DR_old_Sys[nCentBins_old];

	double DR_midPR_Sys[nCentBins_mid];
	double DR_midNP_Sys[nCentBins_mid];


	for (int i=0; i<nCentBins_mid; i++){
		midPR_2S[i] = h_midPR_2S->GetBinContent(i+1);
		midPR_2S_Err[i] = h_midPR_2S->GetBinError(i+1);
		midPR_Jpsi[i] = h_midPR_Jpsi->GetBinContent(i+1);
		midPR_Jpsi_Err[i] = h_midPR_Jpsi->GetBinError(i+1);
		midPR_Jpsi_Sys[i] = midPR_Jpsi[i]*(h_midPR_Sys1S->GetBinContent(i+1));
		midPR_2S_Sys[i] = midPR_2S[i]*(h_midPR_Sys2S->GetBinContent(i+1));
		DR_midPR[i] = midPR_2S[i]/midPR_Jpsi[i];
		DR_midPR_Err[i] = DR_midPR[i]*(TMath::Sqrt(TMath::Power(midPR_2S_Err[i]/midPR_2S[i],2)+TMath::Power(midPR_Jpsi_Err[i]/midPR_Jpsi[i],2)));
		DR_midPR_Sys[i] = DR_midPR[i]*(TMath::Sqrt(TMath::Power(midPR_Jpsi_Sys[i]/midPR_Jpsi[i],2)+TMath::Power(midPR_2S_Sys[i]/midPR_2S[i],2)));
		
		//cout << "Prompt DR	:	" << DR_midPR[i] << "	,	Err :	" << DR_midPR_Err[i] << endl;

		midNP_2S[i] = h_midNP_2S->GetBinContent(i+1);
		midNP_2S_Err[i] = h_midNP_2S->GetBinError(i+1);
		midNP_Jpsi[i] = h_midNP_Jpsi->GetBinContent(i+1);
		midNP_Jpsi_Err[i] = h_midNP_Jpsi->GetBinError(i+1);
		midNP_Jpsi_Sys[i] = midNP_Jpsi[i]*(h_midNP_Sys1S->GetBinContent(i+1));
		midNP_2S_Sys[i] = midNP_2S[i]*(h_midNP_Sys2S->GetBinContent(i+1));
		DR_midNP[i] = midNP_2S[i]/midNP_Jpsi[i];
		DR_midNP_Err[i] = DR_midNP[i]*(TMath::Sqrt(TMath::Power(midNP_2S_Err[i]/midNP_2S[i],2)+TMath::Power(midNP_Jpsi_Err[i]/midNP_Jpsi[i],2)));
		DR_midNP_Sys[i] = DR_midNP[i]*(TMath::Sqrt(TMath::Power(midNP_Jpsi_Sys[i]/midNP_Jpsi[i],2)+TMath::Power(midNP_2S_Sys[i]/midNP_2S[i],2)));

		//cout << "Jpsi Err :	" << midNP_Jpsi_Err[i] << "	,	psi2S Err :	" << midNP_2S_Err[i] << endl;
		//cout << "Prompt DR	:	" << DR_midPR[i] << "	,	Err :	" << DR_midPR_Err[i] << endl;

		cout << "Prompt Err :	" << DR_midPR_Err[i] << " ,	NonPrompt Err :	" << DR_midNP_Err[i] << endl;
		//SysPR[i] = midPR_2S[i]*(hSys_PR->GetBinContent(i+1));
	}

	for (int i=nCentBins_old-1; i>=0;i--){
        DR_old[nCentBins_old-i-1] = h_Old->GetBinContent(i+1);
        DR_old_Err[nCentBins_old-i-1] = h_Old_Err->GetBinContent(i+1);
        DR_old_Sys[nCentBins_old-i-1] = h_Old_Sys->GetBinContent(i+1);
        NpartBin_old[nCentBins_old-i-1] = h_OldNpart->GetBinContent(i+1);
        cout << "i : "<< nCentBins_old-i-1 << ",    DR Old :    " << DR_old[nCentBins_old-i-1] << ",    NPart : " << NpartBin_old[nCentBins_old-i-1] << endl;
    }

	TGraphErrors *g_midPR;
	TGraphErrors *g_midNP;
	TGraphErrors *g_midPR_Sys;
	TGraphErrors *g_midNP_Sys;
	TGraphErrors *g_old;
	TGraphErrors *g_old_Sys;
	TGraphErrors *g_ATLAS = (TGraphErrors*)f_ATLAS->Get("Table 16/Graph1D_y1");
	TGraphAsymmErrors *g_ATLAS_Sys = (TGraphAsymmErrors*)f_ATLAS->Get("Table 16/Graph1D_y1");
	TGraphErrors *g_ATLAS_NP = (TGraphErrors*)f_ATLAS_NP->Get("Table 17/Graph1D_y1");
	TGraphAsymmErrors *g_ATLAS_NP_Sys = (TGraphAsymmErrors*)f_ATLAS_NP->Get("Table 17/Graph1D_y1");
	
	TH1D *h_ATLAS_NP_Sys = (TH1D*)f_ATLAS_NP->Get("Table 17/Hist1D_y1_e2plus");
	TH1D *h_ATLAS_Sys = (TH1D*)f_ATLAS->Get("Table 16/Hist1D_y1_e2plus");
	double ATLAS_NP_Sys[3]; double ATLAS_Sys[5];

	
	for(int i=0; i<3; i++) {
		ATLAS_NP_Sys[i]=h_ATLAS_NP_Sys->GetBinContent(2*i+1);;
		//cout << ATLAS_NP_Sys[i] << endl;
		g_ATLAS_NP_Sys->SetPointError(i,4.3,4.3,ATLAS_NP_Sys[i],ATLAS_NP_Sys[i]);
	}
	
	for(int i=0; i<5; i++){
		ATLAS_Sys[i]=h_ATLAS_Sys->GetBinContent(2*i+1);
		cout << ATLAS_Sys[i] << endl;
		g_ATLAS_Sys->SetPointError(i,4.3,4.3,ATLAS_Sys[i],ATLAS_Sys[i]);
	}

	g_midPR = new TGraphErrors(nCentBins_mid, NpartBin_mid, DR_midPR, 0, DR_midPR_Err);
	g_midNP = new TGraphErrors(nCentBins_mid, NpartBin_mid, DR_midNP, 0, DR_midNP_Err);
    g_old = new TGraphErrors(nCentBins_old-1, NpartBin_old, DR_old, 0, DR_old_Err);
    g_old_Sys = new TGraphErrors(nCentBins_old-1, NpartBin_old, DR_old, binWidth_old, DR_old_Sys);
	g_midPR_Sys = new TGraphErrors(nCentBins_mid, NpartBin_mid, DR_midPR, binWidth_mid, DR_midPR_Sys);
	g_midNP_Sys = new TGraphErrors(nCentBins_mid, NpartBin_mid, DR_midNP, binWidth_mid, DR_midNP_Sys);


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
	float text_size = 20;

	TCanvas *c1 = new TCanvas("c1","",1000,700);
    TPad *pad1_1, *pad1_2;
    pad1_1 = new TPad("pad1_1", "", 0.00,0,0.83,1);
    pad1_1 -> SetRightMargin(0);
    pad1_2 = new TPad("pad1_2", "", 0.83,0,1,1);
    pad1_2 -> SetLeftMargin(0);
    pad1_1->SetTicks();
    pad1_2->SetTicks();
    c1->cd();
    pad1_1->Draw();
    pad1_2->Draw();

	pad1_1->cd();

	g_midPR->GetXaxis()->SetTitle("<N_{Part}>");
	g_midPR->GetXaxis()->CenterTitle();
	g_midPR->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_midPR->GetYaxis()->CenterTitle();
	g_midPR->SetTitle();
	g_midPR->GetXaxis()->SetLimits(0.,400);
	g_midPR->SetMinimum(0.);
	g_midPR->SetMaximum(1.80);

	g_midPR->SetMarkerColor(kBlue+2);
	g_midPR->SetLineColor(kBlue+2);
	g_midPR_Sys->SetLineColor(kBlue-4);
	g_midPR_Sys->SetFillColorAlpha(kBlue-9,0.40);
	//g_midPR->SetLineWidth(2);
	g_midPR->SetMarkerStyle(20);
	g_midPR->SetMarkerSize(1.6);

	g_old->GetXaxis()->SetTitle("<N_{Part}>");
    g_old->GetXaxis()->CenterTitle();
    g_old->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
    g_old->GetYaxis()->CenterTitle();
    g_old->SetTitle();
    g_old->GetXaxis()->SetLimits(0.,400);
    g_old->SetMinimum(0.);
    g_old->SetMaximum(1.80);

    g_old->SetMarkerColor(15);
    g_old->SetLineColor(15);
    g_old->SetMarkerStyle(30);
    g_old->SetMarkerSize(2.0);
	g_old_Sys->SetLineColor(15);
	g_old_Sys->SetFillColorAlpha(15,0.40);

	old_upper->SetLineColor(15);

	g_ATLAS->GetXaxis()->SetTitle("<N_{Part}>");
	g_ATLAS->GetXaxis()->CenterTitle();
	g_ATLAS->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_ATLAS->GetYaxis()->CenterTitle();
	g_ATLAS->SetTitle();
	g_ATLAS->GetXaxis()->SetLimits(0.,400);
	g_ATLAS->SetMinimum(0.);
	g_ATLAS->SetMaximum(1.80);
	g_ATLAS->SetMarkerColor(45);
	g_ATLAS->SetLineColor(45);
	g_ATLAS->SetMarkerStyle(33);
	g_ATLAS->SetMarkerSize(2.0);
	g_ATLAS_Sys->SetLineColor(45);
	g_ATLAS_Sys->SetFillColorAlpha(45,0.40);

	g_midPR->Draw("AP");
	g_midPR_Sys->Draw("5");
	//g_fwdPR->Draw("P");
	g_old->Draw("P");
	g_old_Sys->Draw("5");
	g_ATLAS->Draw("P");
	g_ATLAS_Sys->Draw("5");
    old_upper->Draw("");

	TLegend *leg1 = new TLegend(0.17,0.87,0.44,0.71);
	//SetLegendStyle(leg);
    leg1->SetFillColor(0);
    //leg->SetFillStyle(4000);
    leg1->SetBorderSize(0);
    leg1->SetMargin(0.2);
    leg1->SetTextSize(0.029);
    leg1->SetTextFont(42);
	leg1->AddEntry(g_midPR,"#bf{Prompt}, 6.5 < p_{T} < 40 GeV/c, |y| < 1.6,  Cent. 0-90%", "pe");
	//leg1->AddEntry(g_fwdPR,"#bf{Prompt}, 1.6 < |y| < 2.4, Cent. 0-90%");
	leg1->AddEntry(g_old,"#bf{Prompt (Old)}, 6.5 < p_{T} < 30 GeV/c, |y| < 1.6, Cent. 0-100%", "pe");
	leg1->AddEntry(g_ATLAS,"#bf{Prompt ATLAS}, 9 < p_{T} < 40 GeV/c, |y| < 2", "pe");
	leg1->Draw("SAME");
	jumSun(0,1,400,1);
	CMS_lumi_v2mass(pad1_1, iPeriod, iPos);

	pad1_2->cd();

	//double midPR2S_int = 0.205; double midPR2S_intErr = 0.0112;
	double midPR2S_int; double midPR2S_intErr;
	double midPRJpsi_int; double midPRJpsi_intErr;

	TH1D *h_midPR_2S_int = (TH1D*) f_mid_2S->Get("hmidPR_int");
	TH1D *h_midPR_Jpsi_int = (TH1D*) f_mid_Jpsi->Get("hmidPR_int");
	midPR2S_int = h_midPR_2S_int->GetBinContent(2);
	midPR2S_intErr = h_midPR_2S_int->GetBinError(2);
	midPRJpsi_int = h_midPR_Jpsi_int->GetBinContent(2);
	midPRJpsi_intErr = h_midPR_Jpsi_int->GetBinError(2);

	double DR_midPR_int = midPR2S_int/midPRJpsi_int;
	double DR_midPR_intErr = DR_midPR_int*( TMath::Sqrt( TMath::Power(midPR2S_intErr/midPR2S_int,2) + TMath::Power(midPRJpsi_intErr/midPRJpsi_int,2) ) );
	cout << "Double ratio : " << DR_midPR_int << endl;
	cout << "DR Error : " << DR_midPR_intErr << endl;

	TH1D *hmidPR_int = new TH1D("hmidPR_int", "", 3,0,1);
	hmidPR_int->SetBinContent(2, DR_midPR_int);
    hmidPR_int->SetBinError(2, DR_midPR_intErr);
    hmidPR_int->SetMarkerStyle(20);
    hmidPR_int->SetMarkerSize(1.5);
    hmidPR_int->SetLineColor(kBlue+2);
    hmidPR_int->SetMarkerColor(kBlue+2);

	hmidPR_int->GetYaxis()->SetRangeUser(0,1.80);
    hmidPR_int->GetYaxis()->SetLabelOffset(0);
    hmidPR_int->GetXaxis()->SetLabelOffset(0);
    hmidPR_int->GetXaxis()->SetLabelSize(0);
    hmidPR_int->GetXaxis()->SetTickLength(0);
    hmidPR_int->GetYaxis()->SetLabelSize(0.12);

	hmidPR_int->Draw("PE");
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff, text_color, text_size);
	jumSun(0,1,1,1);

	c1->SaveAs("figs/DoubleRatio_Prompt_Npart_midRap_240807.pdf");
	c1->SaveAs("figs/DoubleRatio_Prompt_Npart_midRap_240807.png");


	/////////////// Non Prompt //////////////////////
	TCanvas *c2 = new TCanvas("c2","",1000,700);
    TPad *pad2_1, *pad2_2;
    pad2_1 = new TPad("pad2_1", "", 0.00,0,0.83,1);
    pad2_1 -> SetRightMargin(0);
    pad2_2 = new TPad("pad2_2", "", 0.83,0,1,1);
    pad2_2 -> SetLeftMargin(0);
    pad2_1->SetTicks();
    pad2_2->SetTicks();
    c2->cd();
    pad2_1->Draw();
    pad2_2->Draw();

	pad2_1->cd();

	g_midNP->GetXaxis()->SetTitle("<N_{Part}>");
	g_midNP->GetXaxis()->CenterTitle();
	g_midNP->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_midNP->GetYaxis()->CenterTitle();
	g_midNP->SetTitle();
	g_midNP->GetXaxis()->SetLimits(0.,400);
	g_midNP->SetMinimum(0.);
	g_midNP->SetMaximum(1.80);

	g_midNP->SetMarkerColor(kRed+3);
	g_midNP->SetLineColor(kRed+3);
	//g_midNP->SetLineWidth(2);
	g_midNP->SetMarkerStyle(21);
	g_midNP->SetMarkerSize(1.6);
	g_midNP_Sys->SetLineColor(kRed-4);
	g_midNP_Sys->SetFillColorAlpha(kRed-9,0.40);

	g_ATLAS_NP->GetXaxis()->SetTitle("<N_{Part}>");
	g_ATLAS_NP->GetXaxis()->CenterTitle();
	g_ATLAS_NP->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_ATLAS_NP->GetYaxis()->CenterTitle();
	g_ATLAS_NP->SetTitle();
	g_ATLAS_NP->GetXaxis()->SetLimits(0.,400);
	g_ATLAS_NP->SetMinimum(0.);
	g_ATLAS_NP->SetMaximum(1.80);
	g_ATLAS_NP->SetMarkerColor(45);
	g_ATLAS_NP->SetLineColor(45);
	g_ATLAS_NP->SetMarkerStyle(33);
	g_ATLAS_NP->SetMarkerSize(2.0);
	g_ATLAS_NP_Sys->SetLineColor(44);
	g_ATLAS_NP_Sys->SetFillColorAlpha(44,0.40);

	g_midNP->Draw("AP");
	g_ATLAS_NP->Draw("P");
	g_midNP_Sys->Draw("5");
	g_ATLAS_NP_Sys->Draw("5");

	TLegend *leg2 = new TLegend(0.17,0.87,0.44,0.72);
	//SetLegendStyle(leg);
    leg2->SetFillColor(0);
    //leg->SetFillStyle(4000);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.2);
    leg2->SetTextSize(0.029);
    leg2->SetTextFont(42);
	leg2->AddEntry(g_midNP,"#bf{NonPrompt}, 6.5 < p_{T} < 40 GeV/c, |y| < 1.6,  Cent. 0-90%", "pe");
	leg2->AddEntry(g_ATLAS_NP,"#bf{NonPrompt ATLAS}, 9 < p_{T} < 40 GeV/c, |y| < 2", "pe");
	leg2->Draw("SAME");
	jumSun(0,1,400,1);
	CMS_lumi_v2mass(pad2_1, iPeriod, iPos);

	pad2_2->cd();

	//double midNP2S_int = 0.200; double midNP2S_intErr = 0.0109;
	double midNP2S_int = 0.294; double midNP2S_intErr = 0.0161;
	double midNPJpsi_int = 0.460; double midNPJpsi_intErr = 0.0210;

	double DR_midNP_int = midNP2S_int/midNPJpsi_int;
	double DR_midNP_intErr = DR_midNP_int*( TMath::Sqrt( TMath::Power(midNP2S_intErr/midNP2S_int,2) + TMath::Power(midNPJpsi_intErr/midNPJpsi_int,2) ) );

	TH1D *hmidNP_int = new TH1D("hmidNP_int", "", 3,0,1);
	hmidNP_int->SetBinContent(2, DR_midNP_int);
    hmidNP_int->SetBinError(2, DR_midNP_intErr);
    hmidNP_int->SetMarkerStyle(21);
    hmidNP_int->SetMarkerSize(1.5);
    hmidNP_int->SetLineColor(kRed+3);
    hmidNP_int->SetMarkerColor(kRed+3);

	hmidNP_int->GetYaxis()->SetRangeUser(0,1.80);
    hmidNP_int->GetYaxis()->SetLabelOffset(0);
    hmidNP_int->GetXaxis()->SetLabelOffset(0);
    hmidNP_int->GetXaxis()->SetLabelSize(0);
    hmidNP_int->GetXaxis()->SetTickLength(0);
    hmidNP_int->GetYaxis()->SetLabelSize(0.12);

	hmidNP_int->Draw("PE");
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff, text_color, text_size);
	jumSun(0,1,1,1);

	c2->SaveAs("figs/DoubleRatio_NonPrompt_Npart_midRap_240807.pdf");
	c2->SaveAs("figs/DoubleRatio_NonPrompt_Npart_midRap_240807.png");

	TCanvas *c3 = new TCanvas("c3","",1000,700);
    TPad *pad3_1, *pad3_2;
    pad3_1 = new TPad("pad3_1", "", 0.00,0,0.83,1);
    pad3_1 -> SetRightMargin(0);
    pad3_2 = new TPad("pad3_2", "", 0.83,0,1,1);
    pad3_2 -> SetLeftMargin(0);
    pad3_1->SetTicks();
    pad3_2->SetTicks();
    c3->cd();
    pad3_1->Draw();
    pad3_2->Draw();

	pad3_1->cd();

	g_midPR->GetXaxis()->SetTitle("<N_{Part}>");
	g_midPR->GetXaxis()->CenterTitle();
	g_midPR->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_midPR->GetYaxis()->CenterTitle();
	g_midPR->SetTitle();
	g_midPR->GetXaxis()->SetLimits(0.,400);
	g_midPR->SetMinimum(0.);
	g_midPR->SetMaximum(1.80);

	g_midPR->SetMarkerColor(kBlue+2);
	g_midPR->SetLineColor(kBlue+2);
	g_midPR_Sys->SetLineColor(kBlue-4);
	g_midPR_Sys->SetFillColorAlpha(kBlue-9,0.40);
	//g_midPR->SetLineWidth(2);
	g_midPR->SetMarkerStyle(20);
	g_midPR->SetMarkerSize(1.6);

	g_midNP->GetXaxis()->SetTitle("<N_{Part}>");
	g_midNP->GetXaxis()->CenterTitle();
	g_midNP->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_midNP->GetYaxis()->CenterTitle();
	g_midNP->SetTitle();
	g_midNP->GetXaxis()->SetLimits(0.,400);
	g_midNP->SetMinimum(0.);
	g_midNP->SetMaximum(1.80);

	g_midNP->SetMarkerColor(kRed+3);
	g_midNP->SetLineColor(kRed+3);
	//g_midNP->SetLineWidth(1);
	g_midNP->SetMarkerStyle(21);
	g_midNP->SetMarkerSize(1.6);
	g_midNP_Sys->SetLineColor(kRed-4);
	g_midNP_Sys->SetFillColorAlpha(kRed-9,0.40);

	g_midPR->Draw("AP");
	g_midPR_Sys->Draw("5");
	g_midNP->Draw("P");
	g_midNP_Sys->Draw("5");

	TLegend *leg3 = new TLegend(0.17,0.87,0.44,0.73);
	//SetLegendStyle(leg);
    leg3->SetFillColor(0);
    //leg->SetFillStyle(4000);
    leg3->SetBorderSize(0);
    leg3->SetMargin(0.2);
    leg3->SetTextSize(0.029);
    leg3->SetTextFont(42);
	leg3->AddEntry(g_midPR,"#bf{Prompt}", "pe");
	leg3->AddEntry(g_midNP,"#bf{Nonrompt}", "pe");
	leg3->Draw("SAME");
	jumSun(0,1,400,1);
	CMS_lumi_v2mass(pad3_1, iPeriod, iPos);

	pad3_2->cd();

	hmidPR_int->SetBinContent(2, DR_midPR_int);
    hmidPR_int->SetBinError(2, DR_midPR_intErr);
    hmidPR_int->SetMarkerStyle(20);
    hmidPR_int->SetMarkerSize(1.5);
    hmidPR_int->SetLineColor(kBlue+2);
    hmidPR_int->SetMarkerColor(kBlue+2);

	hmidPR_int->GetYaxis()->SetRangeUser(0,1.80);
    hmidPR_int->GetYaxis()->SetLabelOffset(0);
    hmidPR_int->GetXaxis()->SetLabelOffset(0);
    hmidPR_int->GetXaxis()->SetLabelSize(0);
    hmidPR_int->GetXaxis()->SetTickLength(0);
    hmidPR_int->GetYaxis()->SetLabelSize(0.12);


	hmidNP_int->SetBinContent(2, DR_midNP_int);
    hmidNP_int->SetBinError(2, DR_midNP_intErr);
    hmidNP_int->SetMarkerStyle(21);
    hmidNP_int->SetMarkerSize(1.5);
    hmidNP_int->SetLineColor(kRed+3);
    hmidNP_int->SetMarkerColor(kRed+3);

	hmidNP_int->GetYaxis()->SetRangeUser(0,1.80);
    hmidNP_int->GetYaxis()->SetLabelOffset(0);
    hmidNP_int->GetXaxis()->SetLabelOffset(0);
    hmidNP_int->GetXaxis()->SetLabelSize(0);
    hmidNP_int->GetXaxis()->SetTickLength(0);
    hmidNP_int->GetYaxis()->SetLabelSize(0.12);

	hmidNP_int->Draw("PE");
	hmidPR_int->Draw("PE SAME");
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff, text_color, text_size);
	jumSun(0,1,1,1);

	c3->SaveAs("figs/DoubleRatio_Npart_midRap_240807.pdf");
	c3->SaveAs("figs/DoubleRatio_Npart_midRap_240807.png");

}
