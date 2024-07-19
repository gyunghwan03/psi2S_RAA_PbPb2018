#include "TFile.h"
#include "TH1.h"
#include "TFile.h"
#include <iostream>
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"

using namespace std;

void draw_syst()
{
	gStyle->SetOptStat(0);
    setTDRStyle();

	int iPeriod = 101;
    int iPos = 33;

	const int nFiles = 9;

	//TString SysName[nFiles] = {"Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction", "Acc", "Eff"};
	//TString SysName[nFiles] = {"Total", "Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction", "Eff"};
	TString SysName[nFiles] = {"Total", "Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction", "TNP", "Eff", "Acc"}; //, "b Fraction", "Eff"};
	//TString SysName[nFiles] = {"Total", "Sig. PDF","Sig. Par","Bkg. PDF", "HF", "b Fraction"}; //, "b Fraction", "Eff"};

	TFile *in_pt[nFiles];
	TFile *in_cent[nFiles];
	TFile *in_Tot = TFile::Open("syst_roots/total_syst.root");

	in_pt[1] = TFile::Open("syst_roots/syst_pt_sigPDF.root");
	in_pt[2] = TFile::Open("syst_roots/syst_pt_sigPAR.root");
	in_pt[3] = TFile::Open("syst_roots/syst_pt_bkgPDF.root");
	in_pt[4] = TFile::Open("syst_roots/syst_pt_HF.root");
	in_pt[5] = TFile::Open("syst_roots/syst_pt_bFrac.root");
	in_pt[6] = TFile::Open("syst_roots/syst_pt_TNP.root");
	in_pt[7] = TFile::Open("syst_roots/syst_pt_eff.root");
	in_pt[8] = TFile::Open("syst_roots/syst_pt_acc.root");

	in_cent[1] = TFile::Open("syst_roots/syst_cent_sigPDF.root");
	in_cent[2] = TFile::Open("syst_roots/syst_cent_sigPAR.root");
	in_cent[3] = TFile::Open("syst_roots/syst_cent_bkgPDF.root");
	in_cent[4] = TFile::Open("syst_roots/syst_cent_HF.root");
	in_cent[5] = TFile::Open("syst_roots/syst_cent_bFrac.root");
	in_cent[6] = TFile::Open("syst_roots/syst_cent_TNP.root");
	in_cent[7] = TFile::Open("syst_roots/syst_cent_eff.root");	
	in_cent[8] = TFile::Open("syst_roots/syst_cent_acc.root");

	TH1D *h_mid_pt_PR[nFiles];
	TH1D *h_mid_pt_NP[nFiles];
	TH1D *h_fwd_pt_PR[nFiles];
	TH1D *h_fwd_pt_NP[nFiles];
	TH1D *h_mid_cent_PR[nFiles];
	TH1D *h_mid_cent_NP[nFiles];
	TH1D *h_fwd_cent_PR[nFiles];
	TH1D *h_fwd_cent_NP[nFiles];

	h_mid_pt_PR[0] = (TH1D*) in_Tot->Get("mid_pt_PR");
	h_mid_pt_NP[0] = (TH1D*) in_Tot->Get("mid_pt_NP");
	h_fwd_pt_PR[0] = (TH1D*) in_Tot->Get("fwd_pt_PR");
	h_fwd_pt_NP[0] = (TH1D*) in_Tot->Get("fwd_pt_NP");
	h_mid_cent_PR[0] = (TH1D*) in_Tot->Get("mid_cent_PR");
	h_mid_cent_NP[0] = (TH1D*) in_Tot->Get("mid_cent_NP");
	h_fwd_cent_PR[0] = (TH1D*) in_Tot->Get("fwd_cent_PR");
	h_fwd_cent_NP[0] = (TH1D*) in_Tot->Get("fwd_cent_NP");
	
	double mid_pt_PR_max = h_mid_pt_PR[0]->GetBinContent(h_mid_pt_PR[0]->GetMaximumBin()) + 0.1;
	double mid_pt_NP_max = h_mid_pt_NP[0]->GetBinContent(h_mid_pt_NP[0]->GetMaximumBin()) + 0.1;
	double fwd_pt_PR_max = h_fwd_pt_PR[0]->GetBinContent(h_fwd_pt_PR[0]->GetMaximumBin()) + 0.1;
	double fwd_pt_NP_max = h_fwd_pt_NP[0]->GetBinContent(h_fwd_pt_NP[0]->GetMaximumBin()) + 0.1;
	double mid_cent_PR_max = h_mid_cent_PR[0]->GetBinContent(h_mid_cent_PR[0]->GetMaximumBin()) + 0.1;
	double mid_cent_NP_max = h_mid_cent_NP[0]->GetBinContent(h_mid_cent_NP[0]->GetMaximumBin()) + 0.1;
	double fwd_cent_PR_max = h_fwd_cent_PR[0]->GetBinContent(h_fwd_cent_PR[0]->GetMaximumBin()) + 0.1;
	double fwd_cent_NP_max = h_fwd_cent_NP[0]->GetBinContent(h_fwd_cent_NP[0]->GetMaximumBin()) + 0.1;

	for (int i=1; i<nFiles-1; i++) { // If Acceptance is included, nFiles-1
		h_mid_pt_PR[i] = (TH1D*) in_pt[i]->Get("mid_PR");
		h_mid_pt_NP[i] = (TH1D*) in_pt[i]->Get("mid_NP");
		h_fwd_pt_PR[i] = (TH1D*) in_pt[i]->Get("fwd_PR");
		h_fwd_pt_NP[i] = (TH1D*) in_pt[i]->Get("fwd_NP");
		h_mid_cent_PR[i] = (TH1D*) in_cent[i]->Get("mid_PR");
		h_mid_cent_NP[i] = (TH1D*) in_cent[i]->Get("mid_NP");
		h_fwd_cent_PR[i] = (TH1D*) in_cent[i]->Get("fwd_PR");
		h_fwd_cent_NP[i] = (TH1D*) in_cent[i]->Get("fwd_NP");

		h_mid_pt_PR[8] = (TH1D*) in_pt[8]->Get("mid"); //Acc
		h_mid_pt_NP[8] = (TH1D*) in_pt[8]->Get("mid"); //Acc
		h_fwd_pt_PR[8] = (TH1D*) in_pt[8]->Get("fwd"); //Acc
		h_fwd_pt_NP[8] = (TH1D*) in_pt[8]->Get("fwd"); //Acc
		h_mid_cent_PR[8] = (TH1D*) in_cent[8]->Get("mid");
		h_mid_cent_NP[8] = (TH1D*) in_cent[8]->Get("mid");
		h_fwd_cent_PR[8] = (TH1D*) in_cent[8]->Get("fwd");
		h_fwd_cent_NP[8] = (TH1D*) in_cent[8]->Get("fwd");
		//cout << h_mid_pt_PR[i]->GetBinContent(1) << endl;
		//cout << h_mid_pt_NP[i]->GetBinContent(1) << endl;
	}

	const int NcentBins = 6;
	double centBin[NcentBins+1] = {0,10,20,30,40,50,90};
	h_mid_cent_PR[6] = new TH1D("mid_PR","mid_PR",NcentBins,centBin);
	h_mid_cent_NP[6] = new TH1D("mid_NP","mid_NP",NcentBins,centBin);
	h_fwd_cent_PR[6] = new TH1D("fwd_PR","fwd_PR",NcentBins,centBin);
	h_fwd_cent_NP[6] = new TH1D("fwd_NP","fwd_NP",NcentBins,centBin);

	for(int i=1; i<NcentBins+1; i++)
	{
		h_mid_cent_PR[6]->SetBinContent(i, ((TH1D*) in_cent[6]->Get("mid_PR"))->GetBinContent(i));
		h_mid_cent_NP[6]->SetBinContent(i, ((TH1D*) in_cent[6]->Get("mid_NP"))->GetBinContent(i));
		h_fwd_cent_PR[6]->SetBinContent(i, ((TH1D*) in_cent[6]->Get("fwd_PR"))->GetBinContent(i));
		h_fwd_cent_NP[6]->SetBinContent(i, ((TH1D*) in_cent[6]->Get("fwd_NP"))->GetBinContent(i));
	}

	TCanvas *c_mid_pt_PR = new TCanvas("c_mid_pt_PR","",700,700);
	TCanvas *c_mid_pt_NP = new TCanvas("c_mid_pt_NP","",700,700);
	TCanvas *c_fwd_pt_PR = new TCanvas("c_fwd_pt_PR","",700,700);
	TCanvas *c_fwd_pt_NP = new TCanvas("c_fwd_pt_NP","",700,700);
	TCanvas *c_mid_cent_PR = new TCanvas("c_mid_cent_PR","",700,700);
	TCanvas *c_mid_cent_NP = new TCanvas("c_mid_cent_NP","",700,700);
	TCanvas *c_fwd_cent_PR = new TCanvas("c_fwd_cent_PR","",700,700);
	TCanvas *c_fwd_cent_NP = new TCanvas("c_fwd_cent_NP","",700,700);

	TLegend *leg_midptPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptPR);
	TLegend *leg_midptNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midptNP);
	TLegend *leg_fwdptPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptPR);
	TLegend *leg_fwdptNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdptNP);

	TLegend *leg_midcentPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentPR);
	TLegend *leg_midcentNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_midcentNP);
	TLegend *leg_fwdcentPR = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentPR);
	TLegend *leg_fwdcentNP = new TLegend(0.58,0.70,0.7,0.90);
    SetLegendStyle(leg_fwdcentNP);

	float pos_x = 0.20;
	float pos_y = 0.87;
	float pos_y_diff = 0.051;
	int text_color = 1;
	float text_size = 23;

	c_mid_pt_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_pt_PR[i]->GetYaxis()->SetRangeUser(0,mid_pt_PR_max);
		h_mid_pt_PR[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_PR[i]->SetLineWidth(2);
		h_mid_pt_PR[i]->SetLineColor(i+1);
		h_mid_pt_PR[i]->Draw("SAME");

		leg_midptPR->AddEntry(h_mid_pt_PR[i],SysName[i]);
		leg_midptPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_PR,iPeriod,iPos);	
	c_mid_pt_PR->SaveAs("./figs/mid_pt_PR.pdf");

	c_mid_pt_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_pt_NP[i]->GetYaxis()->SetRangeUser(0,mid_pt_NP_max);
		h_mid_pt_NP[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_NP[i]->SetLineWidth(2);
		h_mid_pt_NP[i]->SetLineColor(i+1);
		h_mid_pt_NP[i]->Draw("SAME");

		leg_midptNP->AddEntry(h_mid_pt_NP[i],SysName[i]);
		leg_midptNP->Draw("SAME");
	}	
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_pt_NP,iPeriod,iPos);	
	c_mid_pt_NP->SaveAs("./figs/mid_pt_NP.pdf");

	c_fwd_pt_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_pt_PR[i]->GetYaxis()->SetRangeUser(0,fwd_pt_PR_max);
		h_fwd_pt_PR[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_PR[i]->SetLineWidth(2);
		h_fwd_pt_PR[i]->SetLineColor(i+1);
		h_fwd_pt_PR[i]->Draw("SAME");

		leg_fwdptPR->AddEntry(h_fwd_pt_PR[i],SysName[i]);
		leg_fwdptPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_PR,iPeriod,iPos);	
	c_fwd_pt_PR->SaveAs("./figs/fwd_pt_PR.pdf");

	c_fwd_pt_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_pt_NP[i]->GetYaxis()->SetRangeUser(0,fwd_pt_NP_max);
		h_fwd_pt_NP[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_NP[i]->SetLineWidth(2);
		h_fwd_pt_NP[i]->SetLineColor(i+1);
		h_fwd_pt_NP[i]->Draw("SAME");

		leg_fwdptNP->AddEntry(h_fwd_pt_NP[i],SysName[i]);
		leg_fwdptNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_pt_NP,iPeriod,iPos);	
	c_fwd_pt_NP->SaveAs("./figs/fwd_pt_NP.pdf");

	c_mid_cent_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_cent_PR[i]->GetYaxis()->SetRangeUser(0,mid_cent_PR_max);
		h_mid_cent_PR[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_PR[i]->SetLineWidth(2);
		h_mid_cent_PR[i]->SetLineColor(i+1);
		h_mid_cent_PR[i]->Draw("SAME");

		leg_midcentPR->AddEntry(h_mid_cent_PR[i],SysName[i]);
		leg_midcentPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_cent_PR,iPeriod,iPos);	
	c_mid_cent_PR->SaveAs("./figs/mid_cent_PR.pdf");

	c_mid_cent_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_cent_NP[i]->GetYaxis()->SetRangeUser(0,mid_cent_NP_max);
		h_mid_cent_NP[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_NP[i]->SetLineWidth(2);
		h_mid_cent_NP[i]->SetLineColor(i+1);
		h_mid_cent_NP[i]->Draw("SAME");

		leg_midcentNP->AddEntry(h_mid_cent_NP[i],SysName[i]);
		leg_midcentNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("6.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_mid_cent_NP,iPeriod,iPos);	
	c_mid_cent_NP->SaveAs("./figs/mid_cent_NP.pdf");

	c_fwd_cent_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_cent_PR[i]->GetYaxis()->SetRangeUser(0,fwd_cent_PR_max);
		h_fwd_cent_PR[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_PR[i]->SetLineWidth(2);
		h_fwd_cent_PR[i]->SetLineColor(i+1);
		h_fwd_cent_PR[i]->Draw("SAME");

		leg_fwdcentPR->AddEntry(h_fwd_cent_PR[i],SysName[i]);
		leg_fwdcentPR->Draw("SAME");
	}
	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_cent_PR,iPeriod,iPos);	
	c_fwd_cent_PR->SaveAs("./figs/fwd_cent_PR.pdf");

	c_fwd_cent_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_cent_NP[i]->GetYaxis()->SetRangeUser(0,fwd_cent_NP_max);
		h_fwd_cent_NP[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_NP[i]->SetLineWidth(2);
		h_fwd_cent_NP[i]->SetLineColor(i+1);
		h_fwd_cent_NP[i]->Draw("SAME");

		leg_fwdcentNP->AddEntry(h_fwd_cent_NP[i],SysName[i]);
		leg_fwdcentNP->Draw("SAME");
	}
	drawText("NonPrompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	drawText("3.5 < p_{T} < 50 GeV/c", pos_x, pos_y-pos_y_diff*2, text_color, text_size);
    CMS_lumi_v2mass(c_fwd_cent_NP,iPeriod,iPos);	
	c_fwd_cent_NP->SaveAs("./figs/fwd_cent_NP.pdf");




}

