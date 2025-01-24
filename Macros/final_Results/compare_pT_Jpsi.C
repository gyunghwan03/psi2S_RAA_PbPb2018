#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"

void compare_pT_Jpsi(bool isSys = false)
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid = new TFile("roots/RAA_JPsi_midRap_pT.root"); 
	TFile *f_fwd = new TFile("roots/RAA_JPsi_forRap_pT.root"); 
	TFile *fSys = new TFile("../syst_summary_Jpsi/syst_roots/total_syst.root");
	TFile *fPR_mid_old = new TFile("roots/RAA_PR_Jpsi_HIN_16_025_mid_pT.root");
	TFile *fPR_fwd_old = new TFile("roots/RAA_PR_Jpsi_HIN_16_025_fwd_pT.root");
	TFile *fNP_mid_old = new TFile("roots/RAA_NP_Jpsi_HIN_16_025_mid_pT.root");
	TFile *fNP_fwd_old = new TFile("roots/RAA_NP_Jpsi_HIN_16_025_fwd_pT.root");

	TH1D *h_midPR = (TH1D*) f_mid->Get("hRAA_PR");
	TH1D *h_fwdPR = (TH1D*) f_fwd->Get("hRAA_PR");

	TH1D *hSys_mid_PR = (TH1D*) fSys->Get("mid_pt_PR");
	TH1D *hSys_fwd_PR = (TH1D*) fSys->Get("fwd_pt_PR");

	TH1F *h_PR_Jpsi_mid_old    = (TH1F*) fPR_mid_old->Get("Table 18/Hist1D_y1");
    TH1F *h_PR_Jpsi_mid_oldErr = (TH1F*) fPR_mid_old->Get("Table 18/Hist1D_y1_e1");
    TH1F *h_PR_Jpsi_mid_oldSys = (TH1F*) fPR_mid_old->Get("Table 18/Hist1D_y1_e2");

    TH1F *h_PR_Jpsi_fwd_old    = (TH1F*) fPR_fwd_old->Get("Table 22/Hist1D_y1");
    TH1F *h_PR_Jpsi_fwd_oldErr = (TH1F*) fPR_fwd_old->Get("Table 22/Hist1D_y1_e1");
    TH1F *h_PR_Jpsi_fwd_oldSys = (TH1F*) fPR_fwd_old->Get("Table 22/Hist1D_y1_e2");

	TH1D *h_midNP = (TH1D*) f_mid->Get("hRAA_NP");
	TH1D *h_fwdNP = (TH1D*) f_fwd->Get("hRAA_NP");

	TH1D *hSys_mid_NP = (TH1D*) fSys->Get("mid_pt_NP");
	TH1D *hSys_fwd_NP = (TH1D*) fSys->Get("fwd_pt_NP");

	TH1F *h_NP_Jpsi_mid_old    = (TH1F*) fNP_mid_old->Get("Table 28/Hist1D_y1");
    TH1F *h_NP_Jpsi_mid_oldErr = (TH1F*) fNP_mid_old->Get("Table 28/Hist1D_y1_e1");
    TH1F *h_NP_Jpsi_mid_oldSys = (TH1F*) fNP_mid_old->Get("Table 28/Hist1D_y1_e2");

    TH1F *h_NP_Jpsi_fwd_old    = (TH1F*) fNP_fwd_old->Get("Table 29/Hist1D_y1");
    TH1F *h_NP_Jpsi_fwd_oldErr = (TH1F*) fNP_fwd_old->Get("Table 29/Hist1D_y1_e1");
    TH1F *h_NP_Jpsi_fwd_oldSys = (TH1F*) fNP_fwd_old->Get("Table 29/Hist1D_y1_e2");


	const int nPtBins=6;
	const int nPtBins_fwd=4;
	const int nBins_mid_old = 6;
	const int nBins_fwd_old = 10;
	double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,40};
	double ptBin_fwd[nPtBins_fwd+1] = {3.5,6.5,9,12,40};
	double ptBin_old[nPtBins] = {6.5,9,12,15,20,30};
	double ptBin_fwd_old[nPtBins_fwd] = {3.5,6.5,12,30};
	double x[nPtBins]; double binWidth[nPtBins]; 
	double x_fwd[nPtBins_fwd]; double binWidth_fwd[nPtBins_fwd];
	double x_old[nPtBins]; double binWidth_old[nPtBins-1];
	double x_fwd_old[nPtBins_fwd-1]; double binWidth_fwd_old[nPtBins_fwd];
	double x_NP_mid[nBins_mid_old]; double binWidth_NP_mid[nBins_mid_old];
	double x_NP_fwd[nBins_fwd_old]; double binWidth_NP_fwd[nBins_fwd_old];

	double midPR_new[nPtBins]; double midPR_new_Err[nPtBins]; double midPR_new_Sys[nPtBins]; 
	double fwdPR_new[nPtBins_fwd]; double fwdPR_new_Err[nPtBins_fwd]; double fwdPR_new_Sys[nPtBins_fwd]; 

	double midPR_Sys[nPtBins]; 
	double fwdPR_Sys[nPtBins_fwd];

	double SysPR_mid_old[nPtBins]; 
	double SysPR_fwd_old[nPtBins_fwd];

 	double midPR_old[nPtBins-1]; double midPR_old_Err[nPtBins-1]; double midPR_old_Sys[nPtBins-1];
	double fwdPR_old[nPtBins_fwd-1]; double fwdPR_old_Err[nPtBins_fwd-1]; double fwdPR_old_Sys[nPtBins_fwd-1];


	double midNP_new[nPtBins]; double midNP_new_Err[nPtBins]; double midNP_new_Sys[nPtBins]; 
	double fwdNP_new[nPtBins_fwd]; double fwdNP_new_Err[nPtBins_fwd]; double fwdNP_new_Sys[nPtBins_fwd]; 

	double midNP_Sys[nPtBins]; 
	double fwdNP_Sys[nPtBins_fwd];

	double SysNP_mid_old[nBins_mid_old]; 
	double SysNP_fwd_old[nBins_fwd_old];

	double midNP_old[nBins_mid_old]; double midNP_old_Err[nBins_mid_old]; double midNP_old_Sys[nBins_mid_old];
	double fwdNP_old[nBins_fwd_old]; double fwdNP_old_Err[nBins_fwd_old]; double fwdNP_old_Sys[nBins_fwd_old];

	for (int i=0; i<nPtBins; i++){
		midPR_new[i] = h_midPR->GetBinContent(i+1);
		midPR_new_Err[i] = h_midPR->GetBinError(i+1);
		midNP_new[i] = h_midNP->GetBinContent(i+1);
		midNP_new_Err[i] = h_midNP->GetBinError(i+1);
		x[i] = (ptBin[i+1]+ptBin[i])/2;
		binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		midPR_Sys[i] = midPR_new[i]*(hSys_mid_PR->GetBinContent(i+1));
		midNP_Sys[i] = midNP_new[i]*(hSys_mid_NP->GetBinContent(i+1));
		cout << "pt" << ptBin[i] << " - " << ptBin[i+1]  << "\tmid PR RAA : " << midPR_new[i] << " +/- " << midPR_new_Err[i] << "\t(stat.) +/- " << midPR_Sys[i] << " \t(syst.)" << endl;
	}
	for(int i=0; i<nPtBins-1; i++){
		x_old[i] = h_PR_Jpsi_mid_old->GetXaxis()->GetBinCenter(i+1);
		binWidth_old[i] = h_PR_Jpsi_mid_old->GetXaxis()->GetBinWidth(i+1)/2; 
		midPR_old[i] = h_PR_Jpsi_mid_old->GetBinContent(i+1);
		midPR_old_Err[i] = h_PR_Jpsi_mid_oldErr->GetBinContent(i+1);
		midPR_old_Sys[i] = h_PR_Jpsi_mid_oldSys->GetBinContent(i+1);
	}
	cout << " " << endl;
	for (int i=0; i<nPtBins_fwd; i++){
		fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
		fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
		fwdNP_new[i] = h_fwdNP->GetBinContent(i+1);
		fwdNP_new_Err[i] = h_fwdNP->GetBinError(i+1);
		x_fwd[i] = (ptBin_fwd[i+1]+ptBin_fwd[i])/2;
		binWidth_fwd[i] = (ptBin_fwd[i+1]-ptBin_fwd[i])/2;
		fwdPR_Sys[i] = fwdPR_new[i]*(hSys_fwd_PR->GetBinContent(i+1));
		fwdNP_Sys[i] = fwdNP_new[i]*(hSys_fwd_NP->GetBinContent(i+1));
		cout << "pt" << ptBin_fwd[i] << " - " << ptBin_fwd[i+1]  << "\tfwd PR RAA : " << fwdPR_new[i] << " +/- " << fwdPR_new_Err[i] << "\t(stat.) +/- " << fwdPR_Sys[i] << " \t(syst.)" << endl;
	}
	for(int i=0; i<nPtBins_fwd-1; i++){
		x_fwd_old[i] = h_PR_Jpsi_fwd_old->GetXaxis()->GetBinCenter(i+1);
		binWidth_fwd_old[i] = h_PR_Jpsi_fwd_old->GetXaxis()->GetBinWidth(i+1)/2; 
		fwdPR_old[i] = h_PR_Jpsi_fwd_old->GetBinContent(i+1);
		fwdPR_old_Err[i] = h_PR_Jpsi_fwd_oldErr->GetBinContent(i+1);
		fwdPR_old_Sys[i] = h_PR_Jpsi_fwd_oldSys->GetBinContent(i+1);
		cout << "fwdPR_old : " << fwdPR_old[i] << "\tx_fwd_old : " << x_fwd_old[i] << "\tbinWdith_fwd_old : " << binWidth_fwd_old[i] << endl;
	}
	for(int i=0; i<nBins_mid_old; i++){
		x_NP_mid[i] = h_NP_Jpsi_mid_old->GetXaxis()->GetBinCenter(i+1);
		binWidth_NP_mid[i] = h_NP_Jpsi_mid_old->GetXaxis()->GetBinWidth(i+1)/2; 
		midNP_old[i] = h_NP_Jpsi_mid_old->GetBinContent(i+1);
		midNP_old_Err[i] = h_NP_Jpsi_mid_oldErr->GetBinContent(i+1);
		midNP_old_Sys[i] = h_NP_Jpsi_mid_oldSys->GetBinContent(i+1);
	}
	for(int i=0; i<nBins_fwd_old; i++){
		x_NP_fwd[i] = h_NP_Jpsi_fwd_old->GetXaxis()->GetBinCenter(i+1);
		binWidth_NP_fwd[i] = h_NP_Jpsi_fwd_old->GetXaxis()->GetBinWidth(i+1)/2; 
		fwdNP_old[i] = h_NP_Jpsi_fwd_old->GetBinContent(i+1);
		fwdNP_old_Err[i] = h_NP_Jpsi_fwd_oldErr->GetBinContent(i+1);
		fwdNP_old_Sys[i] = h_NP_Jpsi_fwd_oldSys->GetBinContent(i+1);
	}

	TGraphErrors *g_midPR;
	TGraphErrors *g_midPR_old;
	TGraphErrors *g_fwdPR;
	TGraphErrors *g_fwdPR_old;
	TGraphErrors *g_midNP;
	TGraphErrors *g_midNP_old;
	TGraphErrors *g_fwdNP;
	TGraphErrors *g_fwdNP_old;

	if(!isSys) {
		g_midPR     = new TGraphErrors(nPtBins,   x,     midPR_new, binWidth,     midPR_new_Err);
		g_midPR_old = new TGraphErrors(nPtBins-1, x_old, midPR_old, binWidth_old, midPR_old_Err);

		g_fwdPR     = new TGraphErrors(nPtBins_fwd,    x_fwd,     fwdPR_new, binWidth_fwd,     fwdPR_new_Err);
		g_fwdPR_old = new TGraphErrors(nPtBins_fwd-1,  x_fwd_old, fwdPR_old, binWidth_fwd_old, fwdPR_old_Err);

		g_midNP     = new TGraphErrors(nPtBins,   x,     midNP_new, binWidth,     midNP_new_Err);
		g_fwdNP     = new TGraphErrors(nPtBins_fwd,    x_fwd,     fwdNP_new, binWidth_fwd,     fwdNP_new_Err);
		g_midNP_old = new TGraphErrors(nBins_mid_old, x_NP_mid, midNP_old, binWidth_NP_mid, midNP_old_Err);
		g_fwdNP_old = new TGraphErrors(nBins_fwd_old, x_NP_fwd, fwdNP_old, binWidth_NP_fwd, fwdNP_old_Err);
	}
	else {
		g_midPR 	  = new TGraphErrors(nPtBins,   x,      midPR_new, 0, midPR_new_Err);
		g_midPR_old = new TGraphErrors(nPtBins-1, x_old,  midPR_old, 0, midPR_old_Err);

		g_fwdPR     = new TGraphErrors(nPtBins_fwd,     x_fwd,     fwdPR_new, 0, fwdPR_new_Err);
		g_fwdPR_old = new TGraphErrors(nPtBins_fwd-1  , x_fwd_old, fwdPR_old, 0, fwdPR_old_Err);

		g_midNP     = new TGraphErrors(nPtBins,   x,     midNP_new, 0,     midNP_new_Err);
		g_fwdNP     = new TGraphErrors(nPtBins_fwd,    x_fwd,     fwdNP_new, 0,     fwdNP_new_Err);
		g_midNP_old = new TGraphErrors(nBins_mid_old , x_NP_mid, midNP_old, 0, midNP_old_Err);
		g_fwdNP_old = new TGraphErrors(nBins_fwd_old , x_NP_fwd, fwdNP_old, 0, fwdNP_old_Err);
	}
	TGraphErrors *g_midPRSys = new TGraphErrors(nPtBins,     x,     midPR_new, binWidth,     midPR_Sys);
	TGraphErrors *g_fwdPRSys = new TGraphErrors(nPtBins_fwd, x_fwd, fwdPR_new, binWidth_fwd, fwdPR_Sys);
	TGraphErrors *g_midPR_oldSys = new TGraphErrors(nPtBins-1,     x_old,     midPR_old, binWidth_old,     midPR_old_Sys);
	TGraphErrors *g_fwdPR_oldSys = new TGraphErrors(nPtBins_fwd-1, x_fwd_old, fwdPR_old, binWidth_fwd_old, fwdPR_old_Sys);
	TGraphErrors *g_midNPSys     = new TGraphErrors(nPtBins,     x,        midNP_new, binWidth,        midNP_Sys);
	TGraphErrors *g_fwdNPSys     = new TGraphErrors(nPtBins_fwd, x_fwd,    fwdNP_new, binWidth_fwd,    fwdNP_Sys);
	TGraphErrors *g_midNP_oldSys = new TGraphErrors(nBins_mid_old ,   x_NP_mid, midNP_old, binWidth_NP_mid, midNP_old_Sys);
	TGraphErrors *g_fwdNP_oldSys = new TGraphErrors(nBins_fwd_old ,   x_NP_fwd, fwdNP_old, binWidth_NP_fwd, fwdNP_old_Sys);

	int iPeriod = 101;
	int iPos = 33;

	double text_x = 0.18;
	double text_y = 0.8;
	double y_diff = 0.08;
	float pos_x = 0.20;
	float pos_y = 0.87;
	float pos_y_diff = 0.051;
	int text_color = 1;
	float text_size = 25;

	TCanvas *c1 = new TCanvas("c1","",900,800);
	c1->cd();

	g_midPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midPR->GetXaxis()->CenterTitle();
	g_midPR->GetYaxis()->SetTitle("R_{AA}");
	g_midPR->GetYaxis()->CenterTitle();
	g_midPR->SetTitle();
	g_midPR->GetXaxis()->SetLimits(0.,40.);
	g_midPR->SetMinimum(0.);
	g_midPR->SetMaximum(1.44);

	g_midPR->SetMarkerColor(kBlue+2);
	g_midPR->SetLineColor(kBlue+2);
	g_midPR->SetMarkerStyle(20);
	g_midPR->SetMarkerSize(1.4);
	g_midPRSys->SetLineColor(kBlue-4);
	g_midPRSys->SetFillColorAlpha(kBlue-9,0.40);
	g_midPR_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midPR_old->GetXaxis()->CenterTitle();
	g_midPR_old->GetYaxis()->SetTitle("R_{AA}");
	g_midPR_old->GetYaxis()->CenterTitle();
	g_midPR_old->SetTitle();
	g_midPR_old->GetXaxis()->SetLimits(0.,40.);
	g_midPR_old->SetMinimum(0.);
	g_midPR_old->SetMaximum(1.44);
	g_midPR_old->SetMarkerStyle(25);
	g_midPR_old->SetMarkerSize(1.4);
	g_midPR_old->SetMarkerColor(15);
	g_midPR_old->SetLineColor(15);

	g_midPR_oldSys->SetLineColor(15);
	g_midPR_oldSys->SetFillColorAlpha(15,0.4);


	g_midPR->Draw("AP");
	g_midPR_old->Draw("P");
	if(isSys) {
	g_midPRSys->Draw("5");
	g_midPR_oldSys->Draw("5");
	}

	TLegend *leg1 = new TLegend(0.52,0.61,0.76,0.69);
	leg1->SetTextSize(20);
	leg1->SetTextFont(43);
	leg1->SetBorderSize(0);
	leg1->AddEntry(g_midPR,"#bf{Prompt J/#psi}, New Data", "pe");
	leg1->AddEntry(g_midPR_old,"#bf{Prompt J/#psi}, HIN-16-025", "pe");
	leg1->Draw("SAME");
	jumSun(0,1,40,1);

	drawText("6.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c1, iPeriod, iPos);

	c1->SaveAs(Form("./figs/compare_mid_pT_Jpsi_PR_Sys%d.pdf",isSys));

	TCanvas *c2 = new TCanvas("c2","",900,800);
	c2->cd();
	g_fwdPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_fwdPR->GetXaxis()->CenterTitle();
	g_fwdPR->GetYaxis()->SetTitle("R_{AA}");
	g_fwdPR->GetYaxis()->CenterTitle();
	g_fwdPR->SetTitle();
	g_fwdPR->GetXaxis()->SetLimits(0.,40.);
	g_fwdPR->SetMinimum(0.);
	g_fwdPR->SetMaximum(1.44);

	g_fwdPR->SetMarkerColor(kBlue+2);
	g_fwdPR->SetLineColor(kBlue+2);
	g_fwdPR->SetMarkerStyle(20);
	g_fwdPR->SetMarkerSize(1.4);
	g_fwdPRSys->SetLineColor(kBlue-4);
	g_fwdPRSys->SetFillColorAlpha(kBlue-9,0.40);

	g_fwdPR_old->SetMarkerStyle(25);
	g_fwdPR_old->SetMarkerSize(1.4);
	g_fwdPR_old->SetMarkerColor(15);
	g_fwdPR_old->SetLineColor(15);

	g_fwdPR_oldSys->SetLineColor(15);
	g_fwdPR_oldSys->SetFillColorAlpha(15,0.4);

	//fwdPR_upper->SetLineColor(15);
	//fwdPR_upper->SetLineWidth(2);

	g_fwdPR->Draw("AP");
	//g_fwdPRSys->Draw("5");
	g_fwdPR_old->Draw("P");
	//g_fwdPR_oldSys->Draw("5");
	//fwdPR_upper->Draw("");
	if(isSys) {
		g_fwdPRSys->Draw("5");
		g_fwdPR_oldSys->Draw("5");
	}

	TLegend *leg2 = new TLegend(0.52,0.61,0.76,0.69);
	leg2->SetTextSize(20);
	leg2->SetTextFont(43);
	leg2->SetBorderSize(0);
	leg2->AddEntry(g_fwdPR,"#bf{Prompt J/#psi}, New Data", "pe");
	leg2->AddEntry(g_fwdPR_old,"#bf{Propmt J/#psi} HIN-16-025", "pe");
	leg2->Draw("SAME");
	jumSun(0,1,40,1);

	drawText("3.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c2, iPeriod, iPos);

	c2->SaveAs(Form("./figs/compare_fwd_pT_Jpsi_PR_Sys%d.pdf",isSys));

	TCanvas *c3 = new TCanvas("c3","",900,800);
	c3->cd();

	g_midNP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midNP->GetXaxis()->CenterTitle();
	g_midNP->GetYaxis()->SetTitle("R_{AA}");
	g_midNP->GetYaxis()->CenterTitle();
	g_midNP->SetTitle();
	g_midNP->GetXaxis()->SetLimits(0.,40.);
	g_midNP->SetMinimum(0.);
	g_midNP->SetMaximum(1.44);

	g_midNP->SetMarkerColor(kRed+2);
	g_midNP->SetLineColor(kRed+2);
	g_midNP->SetMarkerStyle(20);
	g_midNP->SetMarkerSize(1.4);
	g_midNPSys->SetLineColor(kRed-4);
	g_midNPSys->SetFillColorAlpha(kRed-9,0.40);
	g_midNP_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midNP_old->GetXaxis()->CenterTitle();
	g_midNP_old->GetYaxis()->SetTitle("R_{AA}");
	g_midNP_old->GetYaxis()->CenterTitle();
	g_midNP_old->SetTitle();
	g_midNP_old->GetXaxis()->SetLimits(0.,40.);
	g_midNP_old->SetMinimum(0.);
	g_midNP_old->SetMaximum(1.44);
	g_midNP_old->SetMarkerStyle(25);
	g_midNP_old->SetMarkerSize(1.4);
	g_midNP_old->SetMarkerColor(15);
	g_midNP_old->SetLineColor(15);

	g_midNP_oldSys->SetLineColor(15);
	g_midNP_oldSys->SetFillColorAlpha(15,0.4);


	g_midNP->Draw("AP");
	g_midNP_old->Draw("P");
	if(isSys) {
	g_midNPSys->Draw("5");
	g_midNP_oldSys->Draw("5");
	}

	TLegend *leg3 = new TLegend(0.50,0.61,0.76,0.69);
	leg3->SetTextSize(20);
	leg3->SetTextFont(43);
	leg3->SetBorderSize(0);
	leg3->AddEntry(g_midNP,    "#bf{Nonprompt J/#psi}, New Data", "pe");
	leg3->AddEntry(g_midNP_old,"#bf{Nonprompt J/#psi}, |y| < 0.6, HIN-16-025", "pe");
	leg3->Draw("SAME");
	jumSun(0,1,50,1);

	drawText("6.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c3, iPeriod, iPos);

	c3->SaveAs(Form("./figs/compare_mid_pT_Jpsi_NP_Sys%d.pdf",isSys));


	TCanvas *c4 = new TCanvas("c4","",900,800);
	c4->cd();
	g_fwdNP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_fwdNP->GetXaxis()->CenterTitle();
	g_fwdNP->GetYaxis()->SetTitle("R_{AA}");
	g_fwdNP->GetYaxis()->CenterTitle();
	g_fwdNP->SetTitle();
	g_fwdNP->GetXaxis()->SetLimits(0.,40.);
	g_fwdNP->SetMinimum(0.);
	g_fwdNP->SetMaximum(1.44);

	g_fwdNP->SetMarkerColor(kRed+2);
	g_fwdNP->SetLineColor(kRed+2);
	g_fwdNP->SetMarkerStyle(20);
	g_fwdNP->SetMarkerSize(1.4);
	g_fwdNPSys->SetLineColor(kRed-4);
	g_fwdNPSys->SetFillColorAlpha(kRed-9,0.40);

	g_fwdNP_old->SetMarkerStyle(25);
	g_fwdNP_old->SetMarkerSize(1.4);
	g_fwdNP_old->SetMarkerColor(15);
	g_fwdNP_old->SetLineColor(15);

	g_fwdNP_oldSys->SetLineColor(15);
	g_fwdNP_oldSys->SetFillColorAlpha(15,0.4);

	//fwdNP_upper->SetLineColor(15);
	//fwdNP_upper->SetLineWidth(2);

	g_fwdNP->Draw("AP");
	//g_fwdNPSys->Draw("5");
	g_fwdNP_old->Draw("P");
	//g_fwdNP_oldSys->Draw("5");
	//fwdNP_upper->Draw("");
	if(isSys) {
		g_fwdNPSys->Draw("5");
		g_fwdNP_oldSys->Draw("5");
	}

	TLegend *leg4 = new TLegend(0.50,0.61,0.76,0.69);
	leg4->SetTextSize(20);
	leg4->SetTextFont(43);
	leg4->SetBorderSize(0);
	leg4->AddEntry(g_fwdNP,"#bf{Nonprompt J/#psi}, New Data", "pe");
	leg4->AddEntry(g_fwdNP_old,"#bf{Nonropmt J/#psi}, 1.8 < |y| < 2.4, HIN-16-025", "pe");
	leg4->Draw("SAME");
	jumSun(0,1,40,1);

	drawText("3.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c4, iPeriod, iPos);

	c4->SaveAs(Form("./figs/compare_fwd_pT_Jpsi_NP_Sys%d.pdf",isSys));
}
