#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"
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

void compare_Npart_Jpsi(bool isSys=false)
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid = new TFile("./roots/RAA_JPsi_midRap_Npart.root"); 
	TFile *f_fwd = new TFile("./roots/RAA_JPsi_forRap_Npart_4Bins.root"); 
	TFile *fSys = new TFile("../syst_summary_Jpsi/syst_roots/total_syst.root");
	//TFile *fOld = new TFile("roots/RAA_PR_psi2S_HIN_16_025_Npart.root");
	TFile *fJpsi_mid_old = new TFile("roots/RAA_PR_JPsi_HIN_16_025_mid_Npart.root");
	TFile *fJpsi_fwd_old = new TFile("roots/RAA_PR_JPsi_HIN_16_025_fwd_Npart.root");
	TFile *fJpsi_NP_old = new TFile("roots/RAA_NP_Jpsi_HIN_16_025_Npart.root");
	TFile *fJpsi_NP_old_fwd = new TFile("roots/RAA_NP_Jpsi_HIN_16_025_fwd_Npart.root");

	TH1D *h_midPR = (TH1D*) f_mid->Get("hRAA_PR");
	TH1D *h_fwdPR = (TH1D*) f_fwd->Get("hRAA_PR");

	TH1D *hSys_mid_PR = (TH1D*) fSys->Get("mid_cent_PR");
	TH1D *hSys_fwd_PR = (TH1D*) fSys->Get("fwd_cent_PR");

	TH1D *h_midNP = (TH1D*) f_mid->Get("hRAA_NP");
	TH1D *h_fwdNP = (TH1D*) f_fwd->Get("hRAA_NP");

	TH1D *hSys_mid_NP = (TH1D*) fSys->Get("mid_cent_NP");
	TH1D *hSys_fwd_NP = (TH1D*) fSys->Get("fwd_cent_NP");

	//TH1F *h_PR_psi2S = (TH1F*) fOld->Get("Table 16/Hist1D_y1");
	//TH1F *h_PR_psi2SErr = (TH1F*) fOld->Get("Table 16/Hist1D_y1_e1");
	//TH1F *h_PR_psi2SSys = (TH1F*) fOld->Get("Table 16/Hist1D_y1_e2plus");


	TH1F *h_PR_Jpsi_mid_old    = (TH1F*) fJpsi_mid_old->Get("Table 15/Hist1D_y1");
	TH1F *h_PR_Jpsi_mid_oldErr = (TH1F*) fJpsi_mid_old->Get("Table 15/Hist1D_y1_e1");
	TH1F *h_PR_Jpsi_mid_oldSys = (TH1F*) fJpsi_mid_old->Get("Table 15/Hist1D_y1_e2plus");

	TH1F *h_PR_Jpsi_fwd_old    = (TH1F*) fJpsi_fwd_old->Get("Table 19/Hist1D_y1");
	TH1F *h_PR_Jpsi_fwd_oldErr = (TH1F*) fJpsi_fwd_old->Get("Table 19/Hist1D_y1_e1");
	TH1F *h_PR_Jpsi_fwd_oldSys = (TH1F*) fJpsi_fwd_old->Get("Table 19/Hist1D_y1_e2plus");

	TH1F *h_NP_Jpsi_mid_old    = (TH1F*) fJpsi_NP_old->Get("Table 30/Hist1D_y1");
	TH1F *h_NP_Jpsi_mid_oldErr = (TH1F*) fJpsi_NP_old->Get("Table 30/Hist1D_y1_e1");
	TH1F *h_NP_Jpsi_mid_oldSys = (TH1F*) fJpsi_NP_old->Get("Table 30/Hist1D_y1_e2plus");

	TH1F *h_NP_Jpsi_fwd_old    = (TH1F*) fJpsi_NP_old->Get("Table 30/Hist1D_y2");
	TH1F *h_NP_Jpsi_fwd_oldErr = (TH1F*) fJpsi_NP_old->Get("Table 30/Hist1D_y2_e1");
	TH1F *h_NP_Jpsi_fwd_oldSys = (TH1F*) fJpsi_NP_old->Get("Table 30/Hist1D_y2_e2plus");

	TH1F *h_NP_Jpsi_fwd_Lowpt_old    = (TH1F*) fJpsi_NP_old_fwd->Get("Table 32/Hist1D_y1");
	TH1F *h_NP_Jpsi_fwd_Lowpt_oldErr = (TH1F*) fJpsi_NP_old_fwd->Get("Table 32/Hist1D_y1_e1");
	TH1F *h_NP_Jpsi_fwd_Lowpt_oldSys = (TH1F*) fJpsi_NP_old_fwd->Get("Table 32/Hist1D_y1_e2plus");


	const int nCentBins=6;
	const int nCentBins_fwd=4;
	const int nBins_fwd_old=3;
	double NpartBin[nCentBins] = {27.12,87.19,131.0,188.2,262.3,356.9};
	double NpartBin_mid_old[nCentBins] = {21.9,86.9,131.4,189.2,264.2,358.8};
	double NpartBin_fwd[nCentBins_fwd+1] = {27.12,109.1,225.2,356.9};
	double NpartBin_fwd_old[nBins_fwd_old] = {32.7,160.3,311.5};
	//double NpartBin_alice[4] = {18,70.74,159.4,309.7};
	//double alice_var[4] = {0.441,0.346,0.299,0.342};
	//double alice_err[4] = {0.115,0.084,0.07,0.063};
	//double alice_sys[4] = {0.104,0.057,0.066,0.061};

	double midPR_new[nCentBins]; 
	double midPR_new_Err[nCentBins]; 
	double fwdPR_new[nCentBins_fwd]; 
	double fwdPR_new_Err[nCentBins_fwd]; 

	double SysPR_mid[nCentBins]; 
	double SysPR_fwd[nCentBins_fwd];

	double midNP_new[nCentBins]; 
	double midNP_new_Err[nCentBins]; 
	double fwdNP_new[nCentBins_fwd]; 
	double fwdNP_new_Err[nCentBins_fwd]; 

	double SysNP_mid[nCentBins]; 
	double SysNP_fwd[nCentBins_fwd];

	double midPR_old[nCentBins];
	double midPR_old_Err[nCentBins];
	double midPR_old_Sys[nCentBins];

	double midPR_Jpsi_old[nCentBins];
	double midPR_Jpsi_old_Err[nCentBins];
	double midPR_Jpsi_old_Sys[nCentBins];

	double fwdPR_Jpsi_old[nBins_fwd_old];
	double fwdPR_Jpsi_old_Err[nBins_fwd_old];
	double fwdPR_Jpsi_old_Sys[nBins_fwd_old];

	double midNP_Jpsi_old[nCentBins];
	double midNP_Jpsi_old_Err[nCentBins];
	double midNP_Jpsi_old_Sys[nCentBins];

	double fwdNP_Jpsi_old[nCentBins];
	double fwdNP_Jpsi_old_Err[nCentBins];
	double fwdNP_Jpsi_old_Sys[nCentBins];

	double fwdNP_Jpsi_Lowpt_old[nCentBins];
	double fwdNP_Jpsi_Lowpt_old_Err[nCentBins];
	double fwdNP_Jpsi_Lowpt_old_Sys[nCentBins];

	double binWidth[nCentBins]={4.3,4.3,4.3,4.3,4.3,4.3};

	for (int i=0; i<nCentBins; i++){
		midPR_new[i] = h_midPR->GetBinContent(i+1);
		midPR_new_Err[i] = h_midPR->GetBinError(i+1);
		SysPR_mid[i] = midPR_new[i]*(hSys_mid_PR->GetBinContent(i+1));
		midNP_new[i] = h_midNP->GetBinContent(i+1);
		midNP_new_Err[i] = h_midNP->GetBinError(i+1);
		SysNP_mid[i] = midNP_new[i]*(hSys_mid_NP->GetBinContent(i+1));
		//cout << "midPR_new : " << midPR_new[i] << "\tSys : " << SysPR_mid[i] << endl;
	}
	for (int i=0; i<nCentBins; i++){
		//midPR_old[i] = h_PR_psi2S->GetBinContent(nCentBins-i-1);
		//midPR_old_Err[i] = h_PR_psi2SErr->GetBinContent(nCentBins-i-1);
		//midPR_old_Sys[i] = h_PR_psi2SSys->GetBinContent(nCentBins-i-1);
		midPR_Jpsi_old[i] 	  = h_PR_Jpsi_mid_old->GetBinContent(nCentBins-i);
		midPR_Jpsi_old_Err[i] = h_PR_Jpsi_mid_oldErr->GetBinContent(nCentBins-i);
		midPR_Jpsi_old_Sys[i] = h_PR_Jpsi_mid_oldSys->GetBinContent(nCentBins-i);
		//cout << "i : " << i << ", nCentBins-i-1 : " << nCentBins-i << " Old Val : " << midPR_old[i] << endl;
		//cout << "i : " << i << ", nCentBins-i : " << nCentBins-i << " Old Val : " << midPR_Jpsi_old[i] << endl;
		midNP_Jpsi_old[i] 	  = h_NP_Jpsi_mid_old->GetBinContent(nCentBins-i);
		midNP_Jpsi_old_Err[i] = h_NP_Jpsi_mid_oldErr->GetBinContent(nCentBins-i);
		midNP_Jpsi_old_Sys[i] = h_NP_Jpsi_mid_oldSys->GetBinContent(nCentBins-i);

		fwdNP_Jpsi_old[i] 	  = h_NP_Jpsi_fwd_old->GetBinContent(nCentBins-i);
		fwdNP_Jpsi_old_Err[i] = h_NP_Jpsi_fwd_oldErr->GetBinContent(nCentBins-i);
		fwdNP_Jpsi_old_Sys[i] = h_NP_Jpsi_fwd_oldSys->GetBinContent(nCentBins-i);

		fwdNP_Jpsi_Lowpt_old[i] 	= h_NP_Jpsi_fwd_Lowpt_old->GetBinContent(nCentBins-i);
		fwdNP_Jpsi_Lowpt_old_Err[i] = h_NP_Jpsi_fwd_Lowpt_oldErr->GetBinContent(nCentBins-i);
		fwdNP_Jpsi_Lowpt_old_Sys[i] = h_NP_Jpsi_fwd_Lowpt_oldSys->GetBinContent(nCentBins-i);
	}
	for (int i=0; i<nBins_fwd_old; i++){
		fwdPR_Jpsi_old[i]     = h_PR_Jpsi_fwd_old->GetBinContent(nBins_fwd_old-i);
		fwdPR_Jpsi_old_Err[i] = h_PR_Jpsi_fwd_oldErr->GetBinContent(nBins_fwd_old-i);
		fwdPR_Jpsi_old_Sys[i] = h_PR_Jpsi_fwd_oldSys->GetBinContent(nBins_fwd_old-i);
		//cout << "i : " << i << ", nCentBins-i-1 : " << nCentBins-i << " Old Val : " << midPR_old[i] << endl;
		//cout << "i : " << i << ", nCentBins-i-1 : " << nCentBins-i << " Old Val : " << fwdPR_Jpsi_old[i] << ", stat. Error : " << fwdPR_Jpsi_old_Err[i] << ", syst. Error : " << fwdPR_Jpsi_old_Sys[i] << endl;
	}
	for (int i=0; i<nCentBins_fwd; i++){
		fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
		fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
		SysPR_fwd[i] = fwdPR_new[i]*(hSys_fwd_PR->GetBinContent(i+1));
		fwdNP_new[i] = h_fwdNP->GetBinContent(i+1);
		fwdNP_new_Err[i] = h_fwdNP->GetBinError(i+1);
		SysNP_fwd[i] = fwdNP_new[i]*(hSys_fwd_NP->GetBinContent(i+1));
		//cout << "fwd new val : " << fwdPR_new[i] << "\tSys : " << SysPR_fwd[i] << endl;
	}

	TGraphErrors *g_midPR;
	TGraphErrors *g_Jpsi_midPR_old;
	TGraphErrors *g_Jpsi_fwdPR_old;
	TGraphErrors *g_fwdPR;
	g_midPR = new TGraphErrors(nCentBins, NpartBin, midPR_new, 0, midPR_new_Err);
	g_Jpsi_midPR_old = new TGraphErrors(nCentBins, NpartBin_mid_old, midPR_Jpsi_old, 0, midPR_Jpsi_old_Err);
	g_Jpsi_fwdPR_old = new TGraphErrors(nBins_fwd_old, NpartBin_fwd_old, fwdPR_Jpsi_old, 0, fwdPR_Jpsi_old_Err);
	g_fwdPR = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, fwdPR_new, 0, fwdPR_new_Err);

	TGraphErrors *g_midPRSys = new TGraphErrors(nCentBins,NpartBin,midPR_new,binWidth,SysPR_mid);
	TGraphErrors *g_Jpsi_midPR_oldSys = new TGraphErrors(nCentBins,NpartBin_mid_old,midPR_Jpsi_old,binWidth,midPR_Jpsi_old_Sys);
	TGraphErrors *g_Jpsi_fwdPR_oldSys = new TGraphErrors(nBins_fwd_old,NpartBin_fwd_old,fwdPR_Jpsi_old,binWidth,fwdPR_Jpsi_old_Sys);
	TGraphErrors *g_fwdPRSys = new TGraphErrors(nCentBins_fwd,NpartBin_fwd,fwdPR_new,binWidth,SysPR_fwd);

	//TGraphErrors *g_aliceSys = new TGraphErrors(4,NpartBin_alice,alice_var,binWidth,alice_sys);
	TGraphErrors *g_midNP;
	TGraphErrors *g_Jpsi_midNP_old;
	TGraphErrors *g_Jpsi_fwdNP_old;
	TGraphErrors *g_Jpsi_fwdNP_Lowpt_old;
	TGraphErrors *g_fwdNP;
	g_midNP = new TGraphErrors(nCentBins, NpartBin, midNP_new, 0, midNP_new_Err);
	g_Jpsi_midNP_old = new TGraphErrors(nCentBins, NpartBin_mid_old, midNP_Jpsi_old, 0, midNP_Jpsi_old_Err);
	g_Jpsi_fwdNP_old = new TGraphErrors(nCentBins, NpartBin_mid_old, fwdNP_Jpsi_old, 0, fwdNP_Jpsi_old_Err);
	g_Jpsi_fwdNP_Lowpt_old = new TGraphErrors(nCentBins, NpartBin_mid_old, fwdNP_Jpsi_Lowpt_old, 0, fwdNP_Jpsi_Lowpt_old_Err);
	g_fwdNP = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, fwdNP_new, 0, fwdNP_new_Err);

	TGraphErrors *g_midNPSys = new TGraphErrors(nCentBins,NpartBin,midNP_new,binWidth,SysNP_mid);
	TGraphErrors *g_Jpsi_midNP_oldSys = new TGraphErrors(nCentBins,NpartBin_mid_old,midNP_Jpsi_old,binWidth,midNP_Jpsi_old_Sys);
	TGraphErrors *g_Jpsi_fwdNP_oldSys = new TGraphErrors(nCentBins,NpartBin_mid_old,fwdNP_Jpsi_old,binWidth,fwdNP_Jpsi_old_Sys);
	TGraphErrors *g_Jpsi_fwdNP_Lowpt_oldSys = new TGraphErrors(nCentBins,NpartBin_mid_old,fwdNP_Jpsi_Lowpt_old,binWidth,fwdNP_Jpsi_Lowpt_old_Sys);
	TGraphErrors *g_fwdNPSys = new TGraphErrors(nCentBins_fwd,NpartBin_fwd,fwdNP_new,binWidth,SysNP_fwd);

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
	g_midPR->GetYaxis()->SetTitle("R_{AA}");
	g_midPR->GetYaxis()->CenterTitle();
	g_midPR->SetTitle();
	g_midPR->GetXaxis()->SetLimits(0.,400.);
	g_midPR->SetMinimum(0.);
	g_midPR->SetMaximum(1.44);

	g_midPR->SetMarkerColor(kBlue+2);
	g_midPR->SetLineColor(kBlue+2);
	g_midPR->SetMarkerStyle(20);
	g_midPR->SetMarkerSize(1.4);
	g_midPRSys->SetLineColor(kBlue-4);
	g_midPRSys->SetFillColorAlpha(kBlue-9,0.40);

	g_Jpsi_midPR_old->GetXaxis()->SetTitle("<N_{Part}>");
	g_Jpsi_midPR_old->GetXaxis()->CenterTitle();
	g_Jpsi_midPR_old->GetYaxis()->SetTitle("R_{AA}");
	g_Jpsi_midPR_old->GetYaxis()->CenterTitle();
	g_Jpsi_midPR_old->SetTitle();
	g_Jpsi_midPR_old->GetXaxis()->SetLimits(0.,400.);
	g_Jpsi_midPR_old->SetMinimum(0.);
	g_Jpsi_midPR_old->SetMaximum(1.44);
	g_Jpsi_midPR_old->SetMarkerStyle(25);
	g_Jpsi_midPR_old->SetMarkerSize(1.4);
	g_Jpsi_midPR_old->SetMarkerColor(15);
	g_Jpsi_midPR_old->SetLineColor(15);

	g_Jpsi_midPR_oldSys->SetLineColor(33);
	g_Jpsi_midPR_oldSys->SetFillColorAlpha(33,0.4);
	g_midPR->Draw("AP");
	g_Jpsi_midPR_old->Draw("P");
	if(isSys==1) {
		g_midPRSys->Draw("5");
		g_Jpsi_midPR_oldSys->Draw("5");
	}

	TLegend *leg1 = new TLegend(0.3,0.6,0.66,0.68);
	leg1->SetTextSize(20);
	leg1->SetTextFont(43);
	leg1->SetBorderSize(0);
	leg1->AddEntry(g_midPR,"#bf{Prompt J/#psi}, New Data");
	//leg1->AddEntry(g_midPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, 6.5 < p_{T} < 30 GeV/c");
	//leg1->AddEntry(g_midPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, #splitline{6.5 < p_{T} < 30 GeV/c}{Cent. 0-100%}");
	leg1->AddEntry(g_Jpsi_midPR_old,"#bf{Prompt J/#psi}, HIN-16-025, 6.5 < p_{T} < 30 GeV/c");
	leg1->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("6.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c1, iPeriod, iPos);

	c1->SaveAs(Form("./figs/compare_PR_mid_Npart_Jpsi.pdf"));

	TCanvas *c2 = new TCanvas("c2","",900,800);
	c2->cd();
	g_fwdPR->GetXaxis()->SetTitle("<N_{Part}>");
	g_fwdPR->GetXaxis()->CenterTitle();
	g_fwdPR->GetYaxis()->SetTitle("R_{AA}");
	g_fwdPR->GetYaxis()->CenterTitle();
	g_fwdPR->SetTitle();
	g_fwdPR->GetXaxis()->SetLimits(0.,400.);
	g_fwdPR->SetMinimum(0.);
	g_fwdPR->SetMaximum(1.44);

	g_fwdPR->SetMarkerColor(kBlue+2);
	g_fwdPR->SetLineColor(kBlue+2);
	g_fwdPR->SetMarkerStyle(20);
	g_fwdPR->SetMarkerSize(1.4);
	g_fwdPRSys->SetLineColor(kBlue-4);
	g_fwdPRSys->SetFillColorAlpha(kBlue-9,0.40);

	g_Jpsi_fwdPR_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_Jpsi_fwdPR_old->GetXaxis()->CenterTitle();
	g_Jpsi_fwdPR_old->GetYaxis()->SetTitle("R_{AA}");
	g_Jpsi_fwdPR_old->GetYaxis()->CenterTitle();
	g_Jpsi_fwdPR_old->SetTitle();
	g_Jpsi_fwdPR_old->GetXaxis()->SetLimits(0.,400.);
	g_Jpsi_fwdPR_old->SetMinimum(0.);
	g_Jpsi_fwdPR_old->SetMaximum(1.44);
	g_Jpsi_fwdPR_old->SetMarkerStyle(25);
	g_Jpsi_fwdPR_old->SetMarkerSize(1.4);
	g_Jpsi_fwdPR_old->SetMarkerColor(15);
	g_Jpsi_fwdPR_old->SetLineColor(15);

	g_Jpsi_fwdPR_oldSys->SetLineColor(33);
	g_Jpsi_fwdPR_oldSys->SetFillColorAlpha(33,0.4);

	//g_alice->GetXaxis()->SetTitle("<N_{Part}>");
	//g_alice->GetXaxis()->CenterTitle();
	//g_alice->GetYaxis()->SetTitle("R_{AA}");
	//g_alice->GetYaxis()->CenterTitle();
	//g_alice->SetTitle();
	//g_alice->GetXaxis()->SetLimits(0.,400.);
	//g_alice->SetMinimum(0.);
	//g_alice->SetMaximum(1.44);

	//g_alice->SetMarkerColor(45);
	//g_alice->SetLineColor(45);
	//g_alice->SetMarkerStyle(33);
	//g_alice->SetMarkerSize(1.8);
	//g_aliceSys->SetLineColor(42);
	//g_aliceSys->SetFillColorAlpha(42,0.40);

	g_fwdPR->Draw("AP");
	//g_fwdPR_old->Draw("P");
	//g_alice->Draw("P");
	g_Jpsi_fwdPR_old->Draw("P");

	if(isSys==1){
		//g_fwdPR_oldSys->Draw("5");
		g_fwdPRSys->Draw("5");
		//g_aliceSys->Draw("5");
		g_Jpsi_fwdPR_oldSys->Draw("5");
	}
	//fwdPR_upper_1->Draw("");
	//fwdPR_upper_2->Draw("");

	TLegend *leg2 = new TLegend(0.3,0.6,0.66,0.68);
	leg2->SetTextSize(20);
	leg2->SetTextFont(43);
	leg2->SetBorderSize(0);
	leg2->AddEntry(g_fwdPR,"#bf{Prompt J/#psi}, New Data");
	//leg2->AddEntry(g_fwdPR_old,"#bf{Propmt #psi(2S)} HIN-16-025, 3 < p_{T} < 30 GeV/c");
	//leg2->AddEntry(g_fwdPR_old,"#bf{Propmt #psi(2S)} HIN-16-025, #splitline{3 < p_{T} < 30 GeV/c}{Cent. 0-100%}");
	leg2->AddEntry(g_Jpsi_fwdPR_old,"#bf{Prompt J/#psi}, HIN-16-025, 3 < p_{T} < 30 GeV/c");
	//leg2->AddEntry(g_alice,"#bf{Inclusive #psi(2S)} ALICE, 2.5 < y < 4,p_{T} < 12 GeV/c");
	leg2->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("3.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c2, iPeriod, iPos);

	c2->SaveAs(Form("./figs/compare_PR_fwd_Npart_Jpsi.pdf"));

	TCanvas *c3 = new TCanvas("c3","",900,800);
	c3->cd();

	g_midNP->GetXaxis()->SetTitle("<N_{Part}>");
	g_midNP->GetXaxis()->CenterTitle();
	g_midNP->GetYaxis()->SetTitle("R_{AA}");
	g_midNP->GetYaxis()->CenterTitle();
	g_midNP->SetTitle();
	g_midNP->GetXaxis()->SetLimits(0.,400.);
	g_midNP->SetMinimum(0.);
	g_midNP->SetMaximum(1.44);

	g_midNP->SetMarkerColor(kRed+2);
	g_midNP->SetLineColor(kRed+2);
	g_midNP->SetMarkerStyle(20);
	g_midNP->SetMarkerSize(1.4);
	g_midNPSys->SetLineColor(kRed-4);
	g_midNPSys->SetFillColorAlpha(kRed-9,0.40);

	g_Jpsi_midNP_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_Jpsi_midNP_old->GetXaxis()->CenterTitle();
	g_Jpsi_midNP_old->GetYaxis()->SetTitle("R_{AA}");
	g_Jpsi_midNP_old->GetYaxis()->CenterTitle();
	g_Jpsi_midNP_old->SetTitle();
	g_Jpsi_midNP_old->GetXaxis()->SetLimits(0.,400.);
	g_Jpsi_midNP_old->SetMinimum(0.);
	g_Jpsi_midNP_old->SetMaximum(1.44);
	g_Jpsi_midNP_old->SetMarkerStyle(25);
	g_Jpsi_midNP_old->SetMarkerSize(1.4);
	g_Jpsi_midNP_old->SetMarkerColor(15);
	g_Jpsi_midNP_old->SetLineColor(15);

	g_Jpsi_midNP_oldSys->SetLineColor(33);
	g_Jpsi_midNP_oldSys->SetFillColorAlpha(33,0.4);
	g_midNP->Draw("AP");
	g_Jpsi_midNP_old->Draw("P");
	if(isSys==1) {
		g_midNPSys->Draw("5");
		g_Jpsi_midNP_oldSys->Draw("5");
	}

	TLegend *leg3 = new TLegend(0.3,0.6,0.66,0.68);
	leg3->SetTextSize(20);
	leg3->SetTextFont(43);
	leg3->SetBorderSize(0);
	leg3->AddEntry(g_midNP,"#bf{Nonrompt J/#psi}, New Data, Cent. 0-90%");
	//leg1->AddEntry(g_midNP_old,"#bf{Prompt #psi(2S)} HIN-16-025, 6.5 < p_{T} < 30 GeV/c");
	//leg1->AddEntry(g_midNP_old,"#bf{Prompt #psi(2S)} HIN-16-025, #splitline{6.5 < p_{T} < 30 GeV/c}{Cent. 0-100%}");
	leg3->AddEntry(g_Jpsi_midNP_old,"#bf{Nonprompt J/#psi}, HIN-16-025, 6.5 < p_{T} < 50 GeV/c, |y| < 0.6");
	leg3->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("6.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c1, iPeriod, iPos);

	c3->SaveAs(Form("./figs/compare_NP_mid_Npart_Jpsi.pdf"));

	TCanvas *c4 = new TCanvas("c4","",900,800);
	c4->cd();
	g_fwdNP->GetXaxis()->SetTitle("<N_{Part}>");
	g_fwdNP->GetXaxis()->CenterTitle();
	g_fwdNP->GetYaxis()->SetTitle("R_{AA}");
	g_fwdNP->GetYaxis()->CenterTitle();
	g_fwdNP->SetTitle();
	g_fwdNP->GetXaxis()->SetLimits(0.,400.);
	g_fwdNP->SetMinimum(0.);
	g_fwdNP->SetMaximum(1.44);

	g_fwdNP->SetMarkerColor(kRed+2);
	g_fwdNP->SetLineColor(kRed+2);
	g_fwdNP->SetMarkerStyle(20);
	g_fwdNP->SetMarkerSize(1.4);
	g_fwdNPSys->SetLineColor(kRed-4);
	g_fwdNPSys->SetFillColorAlpha(kRed-9,0.40);

	g_Jpsi_fwdNP_old->GetXaxis()->SetTitle("<N_{Part}>");
	g_Jpsi_fwdNP_old->GetXaxis()->CenterTitle();
	g_Jpsi_fwdNP_old->GetYaxis()->SetTitle("R_{AA}");
	g_Jpsi_fwdNP_old->GetYaxis()->CenterTitle();
	g_Jpsi_fwdNP_old->SetTitle();
	g_Jpsi_fwdNP_old->GetXaxis()->SetLimits(0.,400.);
	g_Jpsi_fwdNP_old->SetMinimum(0.);
	g_Jpsi_fwdNP_old->SetMaximum(1.44);
	g_Jpsi_fwdNP_old->SetMarkerStyle(25);
	g_Jpsi_fwdNP_old->SetMarkerSize(1.4);
	g_Jpsi_fwdNP_old->SetMarkerColor(15);
	g_Jpsi_fwdNP_old->SetLineColor(15);

	g_Jpsi_fwdNP_oldSys->SetLineColor(33);
	g_Jpsi_fwdNP_oldSys->SetFillColorAlpha(33,0.4);

	g_Jpsi_fwdNP_Lowpt_old->GetXaxis()->SetTitle("<N_{Part}>");
	g_Jpsi_fwdNP_Lowpt_old->GetXaxis()->CenterTitle();
	g_Jpsi_fwdNP_Lowpt_old->GetYaxis()->SetTitle("R_{AA}");
	g_Jpsi_fwdNP_Lowpt_old->GetYaxis()->CenterTitle();
	g_Jpsi_fwdNP_Lowpt_old->SetTitle();
	g_Jpsi_fwdNP_Lowpt_old->GetXaxis()->SetLimits(0.,400.);
	g_Jpsi_fwdNP_Lowpt_old->SetMinimum(0.);
	g_Jpsi_fwdNP_Lowpt_old->SetMaximum(1.44);
	g_Jpsi_fwdNP_Lowpt_old->SetMarkerStyle(27);
	g_Jpsi_fwdNP_Lowpt_old->SetMarkerSize(1.7);
	g_Jpsi_fwdNP_Lowpt_old->SetMarkerColor(15);
	g_Jpsi_fwdNP_Lowpt_old->SetLineColor(15);

	g_Jpsi_fwdNP_Lowpt_oldSys->SetLineColor(33);
	g_Jpsi_fwdNP_Lowpt_oldSys->SetFillColorAlpha(33,0.4);

	g_fwdNP->Draw("AP");
	g_Jpsi_fwdNP_old->Draw("P");
	g_Jpsi_fwdNP_Lowpt_old->Draw("P");

	if(isSys==1){
		g_fwdNPSys->Draw("5");
		g_Jpsi_fwdNP_oldSys->Draw("5");
		g_Jpsi_fwdNP_Lowpt_oldSys->Draw("5");
	}
	TLegend *leg4 = new TLegend(0.28,0.73,0.62,0.8);
	leg4->SetTextSize(20);
	leg4->SetTextFont(43);
	leg4->SetBorderSize(0);
	leg4->AddEntry(g_fwdNP,"#bf{Nonrompt J/#psi}, New Data, Cent. 0-90%");
	leg4->AddEntry(g_Jpsi_fwdNP_Lowpt_old,"#bf{Nonrompt J/#psi}, HIN-16-025, 3 < p_{T} < 6.5 GeV/c, 1.8 < |y| < 2.4");
	leg4->AddEntry(g_Jpsi_fwdNP_old,"#bf{Nonrompt J/#psi}, HIN-16-025, 6.5 < p_{T} < 50 GeV/c, 1.8 < |y| < 2.4");
	leg4->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("3.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c4, iPeriod, iPos);

	c4->SaveAs(Form("./figs/compare_NP_fwd_Npart_Jpsi.pdf"));
}
