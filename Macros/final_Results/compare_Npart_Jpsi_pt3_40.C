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

void compare_Npart_Jpsi_pt3_40(bool isSys=false)
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid = new TFile("./roots/RAA_JPsi_midRap_Npart.root"); 
	TFile *f_fwd = new TFile("./roots/RAA_JPsi_forRap_Npart_pt3_40_4Bins.root"); 
	//TFile *fSys = new TFile("../syst_summary/syst_roots/total_syst.root");
	//TFile *fOld = new TFile("roots/RAA_PR_psi2S_HIN_16_025_Npart.root");
	TFile *fJpsi_mid_old = new TFile("roots/RAA_PR_JPsi_HIN_16_025_mid_Npart.root");
	TFile *fJpsi_fwd_old = new TFile("roots/RAA_PR_JPsi_HIN_16_025_fwd_Npart.root");

	TH1D *h_midPR = (TH1D*) f_mid->Get("hRAA_PR");
	TH1D *h_fwdPR = (TH1D*) f_fwd->Get("hRAA_PR");

	//TH1D *hSys_PR = (TH1D*) fSys->Get("mid_cent_PR");
	//TH1D *hSys_fwd_PR = (TH1D*) fSys->Get("fwd_cent_PR");

	//TH1F *h_PR_psi2S = (TH1F*) fOld->Get("Table 16/Hist1D_y1");
	//TH1F *h_PR_psi2SErr = (TH1F*) fOld->Get("Table 16/Hist1D_y1_e1");
	//TH1F *h_PR_psi2SSys = (TH1F*) fOld->Get("Table 16/Hist1D_y1_e2plus");


	TH1F *h_PR_Jpsi_mid_old    = (TH1F*) fJpsi_mid_old->Get("Table 15/Hist1D_y1");
	TH1F *h_PR_Jpsi_mid_oldErr = (TH1F*) fJpsi_mid_old->Get("Table 15/Hist1D_y1_e1");
	//TH1F *h_PR_Jpsi_mid_oldSys = (TH1F*) fJpsi_mid_old->Get("Table 15/Hist1D_y1_e2plus");

	TH1F *h_PR_Jpsi_fwd_old    = (TH1F*) fJpsi_fwd_old->Get("Table 19/Hist1D_y1");
	TH1F *h_PR_Jpsi_fwd_oldErr = (TH1F*) fJpsi_fwd_old->Get("Table 19/Hist1D_y1_e1");
	//TH1F *h_PR_Jpsi_fwd_oldSys = (TH1F*) fJpsi_fwd_old->Get("Table 19/Hist1D_y1_e2plus");


	const int nCentBins=6;
	const int nCentBins_fwd=3;
	double centBin[nCentBins+1];
	double NpartBin[nCentBins] = {27.12,87.19,131.0,188.2,262.3,356.9};
	double NpartBin_mid_old[nCentBins] = {21.9,86.9,131.4,189.2,264.2,358.8};
	double NpartBin_fwd[nCentBins_fwd+1] = {27.12,109.1,225.2,356.9};
	double NpartBin_fwd_old[nCentBins_fwd] = {32.7,160.3,311.5};
	//double NpartBin_alice[4] = {18,70.74,159.4,309.7};
	//double alice_var[4] = {0.441,0.346,0.299,0.342};
	//double alice_err[4] = {0.115,0.084,0.07,0.063};
	//double alice_sys[4] = {0.104,0.057,0.066,0.061};

	double midPR_new[nCentBins]; 
	double midPR_new_Err[nCentBins]; 
	double fwdPR_new[nCentBins]; 
	double fwdPR_new_Err[nCentBins]; 

	//double SysPR[nCentBins]; 
	//double SysPR_fwd[nCentBins];

	//double midPR_old[nCentBins];
	//double midPR_old_Err[nCentBins];
	//double midPR_old_Sys[nCentBins];

	double midPR_Jpsi_old[nCentBins];
	double midPR_Jpsi_old_Err[nCentBins];
	//double midPR_Jpsi_old_Sys[nCentBins];

	double fwdPR_Jpsi_old[nCentBins];
	double fwdPR_Jpsi_old_Err[nCentBins];
	//double fwdPR_Jpsi_old_Sys[nCentBins];

	//double binWidth[nCentBins]={4.3,4.3,4.3,4.3,4.3,4.3};

	for (int i=0; i<nCentBins; i++){
		midPR_new[i] = h_midPR->GetBinContent(i+1);
		midPR_new_Err[i] = h_midPR->GetBinError(i+1);
		//SysPR[i] = midPR_new[i]*(hSys_PR->GetBinContent(i+1));
		cout << "New Val : " << midPR_new[i]  << endl;
	}
	for (int i=0; i<nCentBins; i++){
		//midPR_old[i] = h_PR_psi2S->GetBinContent(nCentBins-i-1);
		//midPR_old_Err[i] = h_PR_psi2SErr->GetBinContent(nCentBins-i-1);
		//midPR_old_Sys[i] = h_PR_psi2SSys->GetBinContent(nCentBins-i-1);
		midPR_Jpsi_old[i] 	  = h_PR_Jpsi_mid_old->GetBinContent(nCentBins-i);
		midPR_Jpsi_old_Err[i] = h_PR_Jpsi_mid_oldErr->GetBinContent(nCentBins-i);
		//midPR_Jpsi_old_Sys[i] = h_PR_Jpsi_mid_oldSys->GetBinContent(nCentBins-i);
		//cout << "i : " << i << ", nCentBins-i-1 : " << nCentBins-i << " Old Val : " << midPR_old[i] << endl;
		cout << "i : " << i << ", nCentBins-i : " << nCentBins-i << " Old Val : " << midPR_Jpsi_old[i] << endl;
	}
	for (int i=0; i<nCentBins_fwd; i++){
		fwdPR_Jpsi_old[i]     = h_PR_Jpsi_fwd_old->GetBinContent(nCentBins_fwd-i);
		fwdPR_Jpsi_old_Err[i] = h_PR_Jpsi_fwd_oldErr->GetBinContent(nCentBins_fwd-i);
		//fwdPR_Jpsi_old_Sys[i] = h_PR_Jpsi_fwd_oldSys->GetBinContent(nCentBins_fwd-i);
		//cout << "i : " << i << ", nCentBins-i-1 : " << nCentBins-i << " Old Val : " << midPR_old[i] << endl;
		//cout << "i : " << i <<  " Old Val : " << fwdPR_Jpsi_old[i] << ", stat. Error : " << fwdPR_Jpsi_old_Err[i] << ", syst. Error : " << fwdPR_Jpsi_old_Sys[i] << endl;
	}
	for (int i=0; i<nCentBins_fwd+1; i++){
		fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
		fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
		//SysPR_fwd[i] = fwdPR_new[i]*(hSys_fwd_PR->GetBinContent(i+1));
		//cout << "fwd new val : " << fwdPR_new[i] << endl;
	}

	TGraphErrors *g_midPR = new TGraphErrors(nCentBins, NpartBin, midPR_new, 0, midPR_new_Err); 
	//TGraphErrors *g_midPR_old;
	//TGraphErrors *g_Jpsi_midPR;    
	//TGraphErrors *g_Jpsi_midPR_old;
	//TGraphErrors *g_Jpsi_midPR_oldd;
	//TGraphErrors *g_Jpsi_fwdPR;    
	//TGraphErrors *g_Jpsi_fwdPR_old;
	//TGraphErrors *g_fwdPR = new TGraphErrors(nCentBins_fwd+1, NpartBin, fwdPR_new, 0, fwdPR_new_Err);
	TGraphErrors *g_fwdPR = new TGraphErrors(nCentBins_fwd+1, NpartBin_fwd, fwdPR_new, 0, fwdPR_new_Err);
	//TGraphErrors *g_Jpsi_fwdPR_old;
	//TGraphErrors *g_fwdPR_old;
	//TGraphErrors *g_alice;
	TGraphErrors *g_Jpsi_midPR_old = new TGraphErrors(nCentBins, NpartBin_mid_old, midPR_Jpsi_old, 0, midPR_Jpsi_old_Err);
	TGraphErrors *g_Jpsi_fwdPR_old = new TGraphErrors(nCentBins_fwd, NpartBin_fwd_old, fwdPR_Jpsi_old, 0, fwdPR_Jpsi_old_Err);

	//if (isSys == 1)
	//{
	//	g_midPR = new TGraphErrors(nCentBins, NpartBin, midPR_new, 0, midPR_new_Err);
	//	//g_midPR_old = new TGraphErrors(nCentBins - 1, NpartBin_mid_old, midPR_old, 0, midPR_old_Err);
	//	//g_Jpsi_midPR_old = new TGraphErrors(nCentBins, NpartBin_mid_old, midPR_Jpsi_old, 0, midPR_Jpsi_old_Err);
	//	//g_Jpsi_fwdPR_old = new TGraphErrors(nCentBins_fwd, NpartBin_fwd_old, fwdPR_Jpsi_old, 0, fwdPR_Jpsi_old_Err);
	//	g_fwdPR = new TGraphErrors(nCentBins, NpartBin, fwdPR_new, 0, fwdPR_new_Err);
	//	//g_fwdPR_old = new TGraphErrors();
	//	//g_alice = new TGraphErrors(4,NpartBin_alice,alice_var,0,alice_err);
	//}
	//else
	//{
	//	g_midPR = new TGraphErrors(nCentBins, NpartBin, midPR_new, 0, midPR_new_Err);
	//	//g_midPR_old = new TGraphErrors(nCentBins - 1, NpartBin_mid_old, midPR_old, 0, midPR_old_Err);
	//	//g_Jpsi_midPR_old = new TGraphErrors(nCentBins, NpartBin_mid_old, midPR_Jpsi_old, 0, midPR_Jpsi_old_Err);
	//	//g_Jpsi_fwdPR_old = new TGraphErrors(nCentBins_fwd, NpartBin_fwd_old, fwdPR_Jpsi_old, 0, fwdPR_Jpsi_old_Err);
	//	g_fwdPR = new TGraphErrors(nCentBins, NpartBin, fwdPR_new, 0, fwdPR_new_Err);
	//	//g_fwdPR_old = new TGraphErrors();
	//	//g_alice = new TGraphErrors(4,NpartBin_alice,alice_var,0,alice_err);
	//}
	//TGraphErrors *g_midPRSys = new TGraphErrors(nCentBins,NpartBin,midPR_new,binWidth,SysPR);
	//TGraphErrors *g_midPR_oldSys = new TGraphErrors(nCentBins-1,NpartBin_mid_old,midPR_old,binWidth,midPR_old_Sys);
	//TGraphErrors *g_Jpsi_midPR_oldSys = new TGraphErrors(nCentBins,NpartBin_mid_old,midPR_Jpsi_old,binWidth,midPR_Jpsi_old_Sys);
	//TGraphErrors *g_Jpsi_fwdPR_oldSys = new TGraphErrors(3,NpartBin_fwd_old,fwdPR_Jpsi_old,binWidth,fwdPR_Jpsi_old_Sys);
	//TArrow *midPR_upper = new TArrow(358.8,0,358.8,0.138,0.027,"<-|");
	//TGraphErrors *g_fwdPRSys = new TGraphErrors(nCentBins,NpartBin,fwdPR_new,binWidth,SysPR_fwd);
	//TGraphErrors *g_fwdPR_oldSys = new TGraphErrors();

	//TGraphErrors *g_aliceSys = new TGraphErrors(4,NpartBin_alice,alice_var,binWidth,alice_sys);


	//g_fwdPR_old->SetPoint(0,311.5,0.193);
	//g_fwdPR_old->SetPointError(0,0,0.145);
	//g_fwdPR_oldSys->SetPoint(0,311.5,0.193);
	//g_fwdPR_oldSys->SetPointError(0,4.3,0.063);

	//TArrow *fwdPR_upper_1 = new TArrow(160.3, 0, 160.3, 0.290, 0.027, "<-|");
	//TArrow *fwdPR_upper_2 = new TArrow(32.7,0,32.7,0.518,0.027,"<-|");

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

	/*
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
	//g_midPRSys->SetLineColor(kBlue-4);
	//g_midPRSys->SetFillColorAlpha(kBlue-9,0.40);

	//g_midPR_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	//g_midPR_old->GetXaxis()->CenterTitle();
	//g_midPR_old->GetYaxis()->SetTitle("R_{AA}");
	//g_midPR_old->GetYaxis()->CenterTitle();
	//g_midPR_old->SetTitle();
	//g_midPR_old->GetXaxis()->SetLimits(0.,400.);
	//g_midPR_old->SetMinimum(0.);
	//g_midPR_old->SetMaximum(1.44);
	//g_midPR_old->SetMarkerStyle(25);
	//g_midPR_old->SetMarkerSize(1.4);
	//g_midPR_old->SetMarkerColor(15);
	//g_midPR_old->SetLineColor(15);

	//g_midPR_oldSys->SetLineColor(15);
	//g_midPR_oldSys->SetFillColorAlpha(15,0.4);

	g_Jpsi_midPR_old->GetXaxis()->SetTitle("p_{T} (GeV/c)");
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

	//g_Jpsi_midPR_oldSys->SetLineColor(33);
	//g_Jpsi_midPR_oldSys->SetFillColorAlpha(33,0.4);
	   g_midPR_Jpsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	   g_midPR_Jpsi->GetXaxis()->CenterTitle();
	   g_midPR_Jpsi->GetYaxis()->SetTitle("R_{AA}");
	   g_midPR_Jpsi->GetYaxis()->CenterTitle();
	   g_midPR_Jpsi->SetTitle();
	   g_midPR_Jpsi->GetXaxis()->SetLimits(0.,50.);
	   g_midPR_Jpsi->SetMinimum(0.);
	   g_midPR_Jpsi->SetMaximum(1.44);
	   g_midPR_Jpsi->SetMarkerStyle(33);
	   g_midPR_Jpsi->SetMarkerSize(1.4);
	   g_midPR_Jpsi->SetMarkerColor(38);
	   g_midPR_Jpsi->SetLineColor(38);

	   g_midPR_JpsiSys->SetLineColor(38);
	   g_midPR_JpsiSys->SetFillColorAlpha(38,0.4);
	g_midPR->Draw("AP");
	//g_midPR_old->Draw("P");
	g_Jpsi_midPR_old->Draw("P");
	if(isSys==1) {
		//g_midPRSys->Draw("5");
		//g_midPR_oldSys->Draw("5");
		//g_Jpsi_midPR_oldSys->Draw("5");
	}
	//midPR_upper->Draw("");
	//g_midPR_Jpsi->Draw("P");
	//g_midPR_JpsiSys->Draw("5");

	TLegend *leg1 = new TLegend(0.3,0.6,0.66,0.68);
	leg1->SetTextSize(20);
	leg1->SetTextFont(43);
	leg1->SetBorderSize(0);
	leg1->AddEntry(g_midPR,"#bf{Prompt J/#psi}, New Data, Cent. 0-90%");
	//leg1->AddEntry(g_midPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, 6.5 < p_{T} < 30 GeV/c");
	//leg1->AddEntry(g_midPR_old,"#bf{Prompt #psi(2S)} HIN-16-025, #splitline{6.5 < p_{T} < 30 GeV/c}{Cent. 0-100%}");
	leg1->AddEntry(g_Jpsi_midPR_old,"#bf{Prompt J/#psi}, HIN-16-025, 6.5 < p_{T} < 30 GeV/c, Cent. 0-100%");
	leg1->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("6.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c1, iPeriod, iPos);

	c1->SaveAs(Form("./figs/compare_PR_mid_Npart_Jpsi.pdf"));
*/
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
	//g_fwdPRSys->SetLineColor(kBlue-4);
	//g_fwdPRSys->SetFillColorAlpha(kBlue-9,0.40);

	//g_fwdPR_old->SetMarkerStyle(25);
	//g_fwdPR_old->SetMarkerSize(1.4);
	//g_fwdPR_old->SetMarkerColor(15);
	//g_fwdPR_old->SetLineColor(15);

	//g_fwdPR_oldSys->SetLineColor(15);
	//g_fwdPR_oldSys->SetFillColorAlpha(15,0.4);

	//fwdPR_upper_1->SetLineColor(15);
	//fwdPR_upper_1->SetLineWidth(2);
	//fwdPR_upper_2->SetLineColor(15);
	//fwdPR_upper_2->SetLineWidth(2);


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

	//g_Jpsi_fwdPR_oldSys->SetLineColor(33);
	//g_Jpsi_fwdPR_oldSys->SetFillColorAlpha(33,0.4);

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
		//g_fwdPRSys->Draw("5");
		//g_aliceSys->Draw("5");
		//g_Jpsi_fwdPR_oldSys->Draw("5");
	}
	//fwdPR_upper_1->Draw("");
	//fwdPR_upper_2->Draw("");

	TLegend *leg2 = new TLegend(0.3,0.6,0.66,0.68);
	leg2->SetTextSize(20);
	leg2->SetTextFont(43);
	leg2->SetBorderSize(0);
	leg2->AddEntry(g_fwdPR,"#bf{Prompt J/#psi}, New Data, Cent. 0-90%");
	//leg2->AddEntry(g_fwdPR_old,"#bf{Propmt #psi(2S)} HIN-16-025, 3 < p_{T} < 30 GeV/c");
	//leg2->AddEntry(g_fwdPR_old,"#bf{Propmt #psi(2S)} HIN-16-025, #splitline{3 < p_{T} < 30 GeV/c}{Cent. 0-100%}");
	leg2->AddEntry(g_Jpsi_fwdPR_old,"#bf{Prompt J/#psi}, HIN-16-025, 3 < p_{T} < 30 GeV/c, Cent. 0-100%");
	//leg2->AddEntry(g_alice,"#bf{Inclusive #psi(2S)} ALICE, 2.5 < y < 4,p_{T} < 12 GeV/c");
	leg2->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("3 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c2, iPeriod, iPos);

	c2->SaveAs(Form("./figs/compare_PR_fwd_Npart_pt3_40_Jpsi.pdf"));
}
