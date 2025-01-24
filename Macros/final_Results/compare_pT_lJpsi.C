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

void compare_pT_lJpsi()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	TFile *f_mid = new TFile("roots/RAA_JPsi_midRap_pT.root"); 
	TFile *f_fwd = new TFile("roots/RAA_JPsi_forRap_pT.root"); 
	//TFile *fSys = new TFile("../syst_summary/syst_roots/total_syst.root");
	TFile *f_Old_mid = new TFile("roots/RAA_PR_Jpsi_HIN_16_025_mid_pT.root");
	TFile *f_Old_fwd = new TFile("roots/RAA_PR_Jpsi_HIN_16_025_fwd_pT.root");

	TH1D *h_midPR = (TH1D*) f_mid->Get("hRAA_PR");
	TH1D *h_fwdPR = (TH1D*) f_fwd->Get("hRAA_PR");

	//TH1D *hSys_PR = (TH1D*) fSys->Get("mid_pt_PR");
	//TH1D *hSys_fwd_PR = (TH1D*) fSys->Get("fwd_pt_PR");


	const int nPtBins=6;
	const int nPtBins_fwd=4;
	double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,40};
	double ptBin_fwd[nPtBins_fwd+1] = {3.5,5,6.5,12,40};
	double ptBin_old[nPtBins] = {6.5,9,12,15,20,30};
	double ptBin_fwd_old[nPtBins_fwd] = {6.5,12,30};
	double x[nPtBins]; double binWidth[nPtBins]; 
	double x_fwd[nPtBins_fwd]; double binWidth_fwd[nPtBins_fwd];
	double x_old[nPtBins-1]; double binWidth_old[nPtBins-1];
	double x_fwd_old[2]={(6.5+12)/2,(12+30)/2}; double binWidth_fwd_old[2]={(12-6.5)/2,(30-12)/2};

	double midPR_new[nPtBins]; 
	double midPR_new_Err[nPtBins]; 
	//double midPR_new_Sys[nPtBins]; 
	double fwdPR_new[nPtBins_fwd]; 
	double fwdPR_new_Err[nPtBins_fwd]; 
	//double fwdPR_new_Sys[nPtBins_fwd]; 

	//double SysPR[nPtBins]; 
	//double SysPR_fwd[nPtBins_fwd];


	//double midPR_old[nPtBins-1] = {0.109,0.128,0.172,0.142,0.123};
	//double midPR_old_Err[nPtBins-1] = {0.057,0.034,0.048,0.061,0.109};
	//double midPR_old_Sys[nPtBins-1] = {0.018,0.013,0.018,0.020,0.034};
	//double fwdPR_old[2] = {0.135,0.247};
	//double fwdPR_old_Err[2] = {0.079,0.098}; 
	//double fwdPR_old_Sys[2] = {0.038,0.042};

	for (int i=0; i<nPtBins; i++){
		midPR_new[i] = h_midPR->GetBinContent(i+1);
		midPR_new_Err[i] = h_midPR->GetBinError(i+1);
		x[i] = (ptBin[i+1]+ptBin[i])/2;
		x_old[i] = (ptBin_old[i+1]+ptBin_old[i])/2;
		binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		binWidth_old[i] = (ptBin_old[i+1]-ptBin_old[i])/2;
		//SysPR[i] = midPR_new[i]*(hSys_PR->GetBinContent(i+1));
	}
	for (int i=0; i<nPtBins_fwd; i++){
		fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
		fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
		x_fwd[i] = (ptBin_fwd[i+1]+ptBin_fwd[i])/2;
		binWidth_fwd[i] = (ptBin_fwd[i+1]-ptBin_fwd[i])/2;
		//SysPR_fwd[i] = fwdPR_new[i]*(hSys_fwd_PR->GetBinContent(i+1));
	}

	double midPR_lJpsi[nPtBins] = {
		0.299,
		0.302,
		0.322,
		0.341,
		0.388,
		0.420
	};
	
	double fwdPR_lJpsi[nPtBins_fwd] = {
		0.334,
		0.279,
		0.287,
		0.320
	};

	TGraphErrors *g_midPR = new TGraphErrors(nPtBins,x,midPR_new,binWidth,midPR_new_Err);
	TGraphErrors *g_midPR_lJpsi = new TGraphErrors(nPtBins,x,midPR_lJpsi,binWidth,0);
	//TGraphErrors *g_midPRSys = new TGraphErrors(nPtBins,x,midPR_new,binWidth,SysPR);
	//TGraphErrors *g_midPR_old = new TGraphErrors(nPtBins-1,x_old,midPR_old,0,midPR_old_Err);
	//TGraphErrors *g_midPR_oldSys = new TGraphErrors(nPtBins-1,x_old,midPR_old,binWidth_old,midPR_old_Sys);
	TGraphErrors *g_midPR_old = (TGraphErrors*)f_Old_mid->Get("Table 18/Graph1D_y1");

	TGraphErrors *g_fwdPR = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR_new,binWidth_fwd,fwdPR_new_Err);
	TGraphErrors *g_fwdPR_lJpsi = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR_lJpsi,binWidth_fwd,0);
	//TGraphErrors *g_fwdPRSys = new TGraphErrors(nPtBins_fwd,x_fwd,fwdPR_new,binWidth_fwd,SysPR_fwd);
	//TGraphErrors *g_fwdPR_old = new TGraphErrors(2,x_fwd_old,fwdPR_old,0,fwdPR_old_Err);
	//TGraphErrors *g_fwdPR_oldSys = new TGraphErrors(2,x_fwd_old,fwdPR_old,binWidth_fwd_old,fwdPR_old_Sys);
	TGraphErrors *g_fwdPR_old = (TGraphErrors*)f_Old_fwd->Get("Table 22/Graph1D_y1");

	double x_1st = (6.5+3)/2;
	//TArrow *fwdPR_upper = new TArrow(x_1st,0,x_1st,0.395,0.027,"<-|");

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
	//g_midPRSys->SetLineColor(kBlue-4);
	//g_midPRSys->SetFillColorAlpha(kBlue-9,0.40);
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

	g_midPR_lJpsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_midPR_lJpsi->GetXaxis()->CenterTitle();
	g_midPR_lJpsi->GetYaxis()->SetTitle("R_{AA}");
	g_midPR_lJpsi->GetYaxis()->CenterTitle();
	g_midPR_lJpsi->SetTitle();
	g_midPR_lJpsi->GetXaxis()->SetLimits(0.,40.);
	g_midPR_lJpsi->SetMinimum(0.);
	g_midPR_lJpsi->SetMaximum(1.44);
	g_midPR_lJpsi->SetMarkerColor(kAzure-2);
	g_midPR_lJpsi->SetLineColor(kAzure-2);
	g_midPR_lJpsi->SetMarkerStyle(33);
	g_midPR_lJpsi->SetMarkerSize(1.6);
	//g_midPR_oldSys->SetLineColor(15);
	//g_midPR_oldSys->SetFillColorAlpha(15,0.4);


	g_midPR->Draw("AP");
	//g_midPRSys->Draw("5");
	g_midPR_old->Draw("P");
	g_midPR_lJpsi->Draw("P");
	//g_midPR_oldSys->Draw("5");

	TLegend *leg1 = new TLegend(0.50,0.58,0.76,0.69);
	leg1->SetTextSize(20);
	leg1->SetTextFont(43);
	leg1->SetBorderSize(0);
	leg1->AddEntry(g_midPR,"#bf{Prompt} J/#psi, New Data, Cent. 0-90%");
	leg1->AddEntry(g_midPR_lJpsi,"#bf{Prompt} J/#psi, #font[12]{l}_{J/#psi} cut, Cent. 0-90%");
	leg1->AddEntry(g_midPR_old,"#bf{Prompt} J/#psi, HIN-16-025, Cent. 0-100%");
	leg1->Draw("SAME");
	jumSun(0,1,40,1);

	drawText("6.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c1, iPeriod, iPos);

	c1->SaveAs("./figs/compare_mid_pT_Jpsi_PR.pdf");

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
	//g_fwdPRSys->SetLineColor(kBlue-4);
	//g_fwdPRSys->SetFillColorAlpha(kBlue-9,0.40);

	g_fwdPR_old->SetMarkerStyle(25);
	g_fwdPR_old->SetMarkerSize(1.4);
	g_fwdPR_old->SetMarkerColor(15);
	g_fwdPR_old->SetLineColor(15);

	g_fwdPR_lJpsi->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	g_fwdPR_lJpsi->GetXaxis()->CenterTitle();
	g_fwdPR_lJpsi->GetYaxis()->SetTitle("R_{AA}");
	g_fwdPR_lJpsi->GetYaxis()->CenterTitle();
	g_fwdPR_lJpsi->SetTitle();
	g_fwdPR_lJpsi->GetXaxis()->SetLimits(0.,40.);
	g_fwdPR_lJpsi->SetMinimum(0.);
	g_fwdPR_lJpsi->SetMaximum(1.44);

	g_fwdPR_lJpsi->SetMarkerColor(kAzure-2);
	g_fwdPR_lJpsi->SetLineColor(kAzure-2);
	g_fwdPR_lJpsi->SetMarkerStyle(33);
	g_fwdPR_lJpsi->SetMarkerSize(1.6);
	//g_fwdPR_oldSys->SetLineColor(15);
	//g_fwdPR_oldSys->SetFillColorAlpha(15,0.4);

	//fwdPR_upper->SetLineColor(15);
	//fwdPR_upper->SetLineWidth(2);

	g_fwdPR->Draw("AP");
	//g_fwdPRSys->Draw("5");
	g_fwdPR_old->Draw("P");
	g_fwdPR_lJpsi->Draw("P");
	//g_fwdPR_oldSys->Draw("5");
	//fwdPR_upper->Draw("");

	TLegend *leg2 = new TLegend(0.50,0.58,0.76,0.69);
	leg2->SetTextSize(20);
	leg2->SetTextFont(43);
	leg2->SetBorderSize(0);
	leg2->AddEntry(g_fwdPR,"#bf{Prompt} J/#psi, New Data, Cent. 0-90%");
	leg2->AddEntry(g_fwdPR_lJpsi,"#bf{Prompt} J/#psi, #font[12]{l}_{J/#psi} cut, Cent. 0-90%");
	leg2->AddEntry(g_fwdPR_old,"#bf{Propmt} J/#psi HIN-16-025, Cent. 0-100%");
	leg2->Draw("SAME");
	jumSun(0,1,40,1);

	drawText("3.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6< |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(c2, iPeriod, iPos);

	c2->SaveAs("./figs/compare_fwd_pT_Jpsi_PR.pdf");
}
