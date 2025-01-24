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

void DoubleRatio_Charmonia_Npart_fwdRap_241223(bool isSys=false)
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_fwd_2S = new TFile("roots/RAA_psi2S_forRap_Npart_pt3p5_40_241223.root");

	TFile *f_Old = new TFile("./roots/DoubleRatio_midRap_HIN_16_004.root");
	TFile *f_ATLAS = new TFile("./roots/DoubleRatio_ATLAS.root");
	TFile *f_ATLAS_NP = new TFile("./roots/DoubleRatio_ATLAS_NonPrompt.root");

	TFile *f_fwd_Jpsi = new TFile("./roots/RAA_JPsi_forRap_Npart_241223.root");

	TFile *f_Sys1S = new TFile("../syst_summary_Jpsi/syst_roots/total_syst.root");
    TFile *f_Sys2S = new TFile("../syst_summary_psi2S/syst_roots/total_syst.root");

	TH1D *h_fwdPR_2S = (TH1D*) f_fwd_2S->Get("hRAA_PR");
	TH1D *h_fwdPR_Jpsi = (TH1D*) f_fwd_Jpsi->Get("hRAA_PR");

	TH1D *h_fwdNP_2S = (TH1D*) f_fwd_2S->Get("hRAA_NP");
	TH1D *h_fwdNP_Jpsi = (TH1D*) f_fwd_Jpsi->Get("hRAA_NP");
	TH1D *h_fwdPR_Sys1S = (TH1D*) f_Sys1S->Get("fwd_cent_PR");
    TH1D *h_fwdPR_Sys2S = (TH1D*) f_Sys2S->Get("fwd_cent_PR");
    TH1D *h_fwdNP_Sys1S = (TH1D*) f_Sys1S->Get("fwd_cent_PR");
    TH1D *h_fwdNP_Sys2S = (TH1D*) f_Sys2S->Get("fwd_cent_PR");

	TH1D *h_Old = (TH1D*) f_Old->Get("Table 5/Hist1D_y1");
    TH1D *h_Old_Err = (TH1D*) f_Old->Get("Table 5/Hist1D_y1_e1");
    TH1D *h_OldNpart = (TH1D*) f_Old->Get("Table 5/Hist1D_y2");

	const int nCentBins_fwd = 5;
	const int nCentBins_old = 6;
	double centBin[nCentBins_fwd+1] = {0,10,20,30,50,90};
	double NpartBin_fwd[nCentBins_fwd+1] = {27.12,109.1,188.2,262.3,356.9};
	double NpartBin_old[nCentBins_old];

	double x_fwd[nCentBins_fwd];
	double binWidth_fwd[nCentBins_fwd]={4.3,4.3,4.3,4.3,4.3};
	double binWidth_old[nCentBins_old]={4.3,4.3,4.3,4.3,4.3,4.3};

	double fwdPR_2S[nCentBins_fwd];
	double fwdPR_2S_Err[nCentBins_fwd];
	double fwdPR_2S_Sys[nCentBins_fwd];

	double fwdPR_Jpsi[nCentBins_fwd];
	double fwdPR_Jpsi_Err[nCentBins_fwd];
	double fwdPR_Jpsi_Sys[nCentBins_fwd];

	double DR_fwdPR[nCentBins_fwd];
	double DR_fwdPR_Err[nCentBins_fwd];
	double DR_fwdPR_Sys[nCentBins_fwd];

	double fwdNP_2S[nCentBins_fwd];
	double fwdNP_2S_Err[nCentBins_fwd];
	double fwdNP_2S_Sys[nCentBins_fwd];

	double fwdNP_Jpsi[nCentBins_fwd];
	double fwdNP_Jpsi_Err[nCentBins_fwd];
	double fwdNP_Jpsi_Sys[nCentBins_fwd];

	double DR_fwdNP[nCentBins_fwd];
	double DR_fwdNP_Err[nCentBins_fwd];
	double DR_fwdNP_Sys[nCentBins_fwd];

	double DR_old[nCentBins_old];
    double DR_old_Err[nCentBins_old];

	for (int i=0; i<nCentBins_fwd; i++){
		fwdPR_2S[i] = h_fwdPR_2S->GetBinContent(i+1);
		fwdPR_2S_Err[i] = h_fwdPR_2S->GetBinError(i+1);
		fwdPR_Jpsi[i] = h_fwdPR_Jpsi->GetBinContent(i+1);
		fwdPR_Jpsi_Err[i] = h_fwdPR_Jpsi->GetBinError(i+1);
		fwdPR_Jpsi_Sys[i] = fwdPR_Jpsi[i]*(h_fwdPR_Sys1S->GetBinContent(i+1));
        fwdPR_2S_Sys[i] = fwdPR_2S[i]*(h_fwdPR_Sys2S->GetBinContent(i+1));
		DR_fwdPR[i] = fwdPR_2S[i]/fwdPR_Jpsi[i];
		DR_fwdPR_Err[i] = DR_fwdPR[i]*(TMath::Sqrt(TMath::Power(fwdPR_2S_Err[i]/fwdPR_2S[i],2)+TMath::Power(fwdPR_Jpsi_Err[i]/fwdPR_Jpsi[i],2)));
		DR_fwdPR_Sys[i] = DR_fwdPR[i]*(TMath::Sqrt(TMath::Power(fwdPR_Jpsi_Sys[i]/fwdPR_Jpsi[i],2)+TMath::Power(fwdPR_2S_Sys[i]/fwdPR_2S[i],2)));
		//SysPR_fwd[i] = fwdPR_2S[i]*(hSys_fwd_PR
		fwdNP_2S[i] = h_fwdNP_2S->GetBinContent(i+1);
		fwdNP_2S_Err[i] = h_fwdNP_2S->GetBinError(i+1);
		fwdNP_Jpsi[i] = h_fwdNP_Jpsi->GetBinContent(i+1);
		fwdNP_Jpsi_Err[i] = h_fwdNP_Jpsi->GetBinError(i+1);
		fwdNP_Jpsi_Sys[i] = fwdNP_Jpsi[i]*(h_fwdNP_Sys1S->GetBinContent(i+1));
        fwdNP_2S_Sys[i] = fwdNP_2S[i]*(h_fwdNP_Sys2S->GetBinContent(i+1));
		DR_fwdNP[i] = fwdNP_2S[i]/fwdNP_Jpsi[i];
		DR_fwdNP_Err[i] = DR_fwdNP[i]*(TMath::Sqrt(TMath::Power(fwdNP_2S_Err[i]/fwdNP_2S[i],2)+TMath::Power(fwdNP_Jpsi_Err[i]/fwdNP_Jpsi[i],2)));
		DR_fwdNP_Sys[i] = DR_fwdNP[i]*(TMath::Sqrt(TMath::Power(fwdNP_Jpsi_Sys[i]/fwdNP_Jpsi[i],2)+TMath::Power(fwdNP_2S_Sys[i]/fwdNP_2S[i],2)));
		cout << "cent." << centBin[nCentBins_fwd-i-1] << "-" << centBin[nCentBins_fwd-i] << "\tDouble Ratio : " << DR_fwdNP[i] << endl; 
	}

	for (int i=nCentBins_old-1; i>=0;i--){
        DR_old[nCentBins_old-i-1] = h_Old->GetBinContent(i+1);
        DR_old_Err[nCentBins_old-i-1] = h_Old_Err->GetBinContent(i+1);
        NpartBin_old[nCentBins_old-i-1] = h_OldNpart->GetBinContent(i+1);
    }

	TGraphErrors *g_fwdPR;
	TGraphErrors *g_fwdNP;
	TGraphErrors *g_fwdPR_Sys;
	TGraphErrors *g_fwdNP_Sys;
	TGraphErrors *g_old;
	TGraphErrors *g_old_Sys;
	TGraphErrors *g_ATLAS = (TGraphErrors*)f_ATLAS->Get("Table 16/Graph1D_y1");
	TGraphErrors *g_ATLAS_NP = (TGraphErrors*)f_ATLAS_NP->Get("Table 17/Graph1D_y1");
	TGraphAsymmErrors *g_ATLAS_NP_Sys = (TGraphAsymmErrors*)f_ATLAS_NP->Get("Table 17/Graph1D_y1");
	double ATLAS_NP_Sys[3];
	TH1D *h_ATLAS_NP_Sys = (TH1D*)f_ATLAS_NP->Get("Table 17/Hist1D_y1_e2plus");

	for(int i=0; i<3; i++) {
        ATLAS_NP_Sys[i]=h_ATLAS_NP_Sys->GetBinContent(2*i+1);;
        //cout << ATLAS_NP_Sys[i] << endl;
        g_ATLAS_NP_Sys->SetPointError(i,4.3,4.3,ATLAS_NP_Sys[i],ATLAS_NP_Sys[i]);
    }

	g_fwdPR = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, DR_fwdPR, 0, DR_fwdPR_Err);
	g_fwdNP = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, DR_fwdNP, 0, DR_fwdNP_Err);
	g_fwdPR_Sys = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, DR_fwdPR, binWidth_fwd, DR_fwdPR_Sys);
	g_fwdNP_Sys = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, DR_fwdNP, binWidth_fwd, DR_fwdNP_Sys);
    g_old = new TGraphErrors();
	g_old_Sys = new TGraphErrors();

	g_old->SetPoint(0,312,0.52);
	g_old->SetPointError(0,0,0.39);
	g_old_Sys->SetPoint(0,312,0.52);
	g_old_Sys->SetPointError(0,4.3,0.15);

	TArrow *old_upper1 = new TArrow(33, 0, 33, 0.72, 0.027,"<-|");
	TArrow *old_upper2 = new TArrow(160, 0, 160, 0.52, 0.027,"<-|");

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

	g_fwdPR->GetXaxis()->SetTitle("<N_{Part}>");
	g_fwdPR->GetXaxis()->CenterTitle();
	g_fwdPR->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_fwdPR->GetYaxis()->CenterTitle();
	g_fwdPR->SetTitle();
	g_fwdPR->GetXaxis()->SetLimits(0.,400);
	g_fwdPR->SetMinimum(0.);
	g_fwdPR->SetMaximum(1.44);

	g_fwdPR->SetMarkerColor(kBlue+2);
	g_fwdPR->SetLineColor(kBlue+2);
	//g_fwdPR->SetLineWidth(2);
	g_fwdPR->SetMarkerStyle(20);
	g_fwdPR->SetMarkerSize(1.6);
	//g_fwdPR_Sys->SetMarkerColor(kBlue-4);
	g_fwdPR_Sys->SetLineColor(kBlue-4);
    g_fwdPR_Sys->SetFillColorAlpha(kBlue-9,0.40);

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
    g_old->SetMarkerSize(2.0);
    g_old_Sys->SetFillColorAlpha(15,0.40);
    g_old_Sys->SetLineColor(15);

	old_upper1->SetLineColor(15);
	old_upper2->SetLineColor(15);

	g_ATLAS->GetXaxis()->SetTitle("<N_{Part}>");
	g_ATLAS->GetXaxis()->CenterTitle();
	g_ATLAS->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_ATLAS->GetYaxis()->CenterTitle();
	g_ATLAS->SetTitle();
	g_ATLAS->GetXaxis()->SetLimits(0.,400);
	g_ATLAS->SetMinimum(0.);
	g_ATLAS->SetMaximum(1.44);
	g_ATLAS->SetMarkerColor(45);
	g_ATLAS->SetLineColor(45);
	g_ATLAS->SetMarkerStyle(33);
	g_ATLAS->SetMarkerSize(2.0);

	g_fwdPR->Draw("AP");
	g_old->Draw("P");
	if(isSys){
		g_old_Sys->Draw("5");
		g_fwdPR_Sys->Draw("5"); }
	//g_ATLAS->Draw("P");
    old_upper1->Draw("");
    old_upper2->Draw("");

	TLegend *leg1 = new TLegend(0.17,0.87,0.44,0.77);
	leg1->SetFillColor(0);
    //leg->SetFillStyle(4000);
    leg1->SetBorderSize(0);
    leg1->SetMargin(0.2);
    leg1->SetTextSize(0.029);
    leg1->SetTextFont(42);
	leg1->AddEntry(g_fwdPR,"#bf{Prompt}, 3.5 < p_{T} < 40 GeV/c, 1.6 < |y| < 2.4, Cent. 0-90%");
	leg1->AddEntry(g_old,"#bf{Prompt (Old)}, 3 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4, Cent. 0-100%");
	//leg1->AddEntry(g_ATLAS,"#bf{Prompt ATLAS}, 9 < p_{T} < 40 GeV/c, |y| < 2");
	leg1->Draw("SAME");
	jumSun(0,1,400,1);
	CMS_lumi_v2mass(pad1_1, iPeriod, iPos);

	pad1_2->cd();

	double fwdPR2S_int = 0.276; double fwdPR2S_intErr = 0.0294;
    double fwdPRJpsi_int = 0.414; double fwdPRJpsi_intErr = 0.0026;

    double DR_fwdPR_int = fwdPR2S_int/fwdPRJpsi_int;
    double DR_fwdPR_intErr = DR_fwdPR_int*( TMath::Sqrt( TMath::Power(fwdPR2S_intErr/fwdPR2S_int,2) + TMath::Power(fwdPRJpsi_intErr/fwdPRJpsi_int,2) ) );

    TH1D *hfwdPR_int = new TH1D("hfwdPR_int", "", 3,0,1);
    hfwdPR_int->SetBinContent(2, DR_fwdPR_int);
    hfwdPR_int->SetBinError(2, DR_fwdPR_intErr);
    hfwdPR_int->SetMarkerStyle(20);
    hfwdPR_int->SetMarkerSize(1.5);
    hfwdPR_int->SetLineColor(kBlue+2);
    hfwdPR_int->SetMarkerColor(kBlue+2);

    hfwdPR_int->GetYaxis()->SetRangeUser(0,1.44);
    hfwdPR_int->GetYaxis()->SetLabelOffset(0);
    hfwdPR_int->GetXaxis()->SetLabelOffset(0);
    hfwdPR_int->GetXaxis()->SetLabelSize(0);
    hfwdPR_int->GetXaxis()->SetTickLength(0);
    hfwdPR_int->GetYaxis()->SetLabelSize(0.12);

    hfwdPR_int->Draw("PE");
    drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff, text_color, text_size);
    jumSun(0,1,1,1);

	c1->SaveAs("figs/DoubleRatio_Prompt_Npart_fwdRap_241223.pdf");
	c1->SaveAs("figs/DoubleRatio_Prompt_Npart_fwdRap_241223.png");

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

	g_fwdNP->GetXaxis()->SetTitle("<N_{Part}>");
	g_fwdNP->GetXaxis()->CenterTitle();
	g_fwdNP->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
	g_fwdNP->GetYaxis()->CenterTitle();
	g_fwdNP->SetTitle();
	g_fwdNP->GetXaxis()->SetLimits(0.,400);
	g_fwdNP->SetMinimum(0.);
	g_fwdNP->SetMaximum(1.44);

	g_fwdNP->SetMarkerColor(kRed+3);
	g_fwdNP->SetLineColor(kRed+3);
	//g_fwdNP->SetLineWidth(2);
	g_fwdNP->SetMarkerStyle(21);
	g_fwdNP->SetMarkerSize(1.6);
	g_fwdNP_Sys->SetMarkerColor(kRed-4);
	g_fwdNP_Sys->SetLineColor(kRed-4);
    g_fwdNP_Sys->SetFillColorAlpha(kRed-9,0.40);

	g_ATLAS_NP->GetXaxis()->SetTitle("<N_{Part}>");
    g_ATLAS_NP->GetXaxis()->CenterTitle();
    g_ATLAS_NP->GetYaxis()->SetTitle("(#psi(2S)/J/#psi)_{PbPb} / (#psi(2S)/J/#psi)_{pp}");
    g_ATLAS_NP->GetYaxis()->CenterTitle();
    g_ATLAS_NP->SetTitle();
    g_ATLAS_NP->GetXaxis()->SetLimits(0.,400);
    g_ATLAS_NP->SetMinimum(0.);
    g_ATLAS_NP->SetMaximum(1.44);
    g_ATLAS_NP->SetMarkerColor(45);
    g_ATLAS_NP->SetLineColor(45);
    g_ATLAS_NP->SetMarkerStyle(33);
    g_ATLAS_NP->SetMarkerSize(2.0);
	g_ATLAS_NP_Sys->SetLineColor(44);
    g_ATLAS_NP_Sys->SetFillColorAlpha(44,0.40);

	g_fwdNP->Draw("AP");
	g_ATLAS_NP->Draw("P");
	if(isSys) {
		g_fwdNP_Sys->Draw("5");
		g_ATLAS_NP_Sys->Draw("5"); }

	TLegend *leg2 = new TLegend(0.17,0.87,0.44,0.77);
	leg2->SetFillColor(0);
    //leg->SetFillStyle(4000);
    leg2->SetBorderSize(0);
    leg2->SetMargin(0.2);
    leg2->SetTextSize(0.029);
    leg2->SetTextFont(42);
	leg2->AddEntry(g_fwdNP,"#bf{NonPrompt}, 3.5 < p_{T} < 40 GeV/c, 1.6 < |y| < 2.4, Cent. 0-90%");
	//leg2->AddEntry(g_old,"#bf{NonPrompt (Old)}, 3 < p_{T} < 30 GeV/c, 1.6 < |y| < 2.4, Cent. 0-100%");
	leg2->AddEntry(g_ATLAS_NP,"#bf{NonPrompt ATLAS}, 9 < p_{T} < 40 GeV/c, |y| < 2");
	leg2->Draw("SAME");
	jumSun(0,1,400,1);
	CMS_lumi_v2mass(pad2_1, iPeriod, iPos);

	pad2_2->cd();

	double fwdNP2S_int = 0.238; double fwdNP2S_intErr = 0.0254;
    double fwdNPJpsi_int = 0.411; double fwdNPJpsi_intErr = 0.0026;

    double DR_fwdNP_int = fwdNP2S_int/fwdNPJpsi_int;
    double DR_fwdNP_intErr = DR_fwdNP_int*( TMath::Sqrt( TMath::Power(fwdNP2S_intErr/fwdNP2S_int,2) + TMath::Power(fwdNPJpsi_intErr/fwdNPJpsi_int,2) ) );

    TH1D *hfwdNP_int = new TH1D("hfwdNP_int", "", 3,0,1);
    hfwdNP_int->SetBinContent(2, DR_fwdNP_int);
    hfwdNP_int->SetBinError(2, DR_fwdNP_intErr);
    hfwdNP_int->SetMarkerStyle(21);
    hfwdNP_int->SetMarkerSize(1.5);
    hfwdNP_int->SetLineColor(kRed+3);
    hfwdNP_int->SetMarkerColor(kRed+3);

    hfwdNP_int->GetYaxis()->SetRangeUser(0,1.44);
    hfwdNP_int->GetYaxis()->SetLabelOffset(0);
    hfwdNP_int->GetXaxis()->SetLabelOffset(0);
    hfwdNP_int->GetXaxis()->SetLabelSize(0);
    hfwdNP_int->GetXaxis()->SetTickLength(0);
    hfwdNP_int->GetYaxis()->SetLabelSize(0.12);

    hfwdNP_int->Draw("PE");
    drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff, text_color, text_size);
    jumSun(0,1,1,1);

	c2->SaveAs("figs/DoubleRatio_NonPrompt_Npart_fwdRap_241223.pdf");
	c2->SaveAs("figs/DoubleRatio_NonPrompt_Npart_fwdRap_241223.png");

}
