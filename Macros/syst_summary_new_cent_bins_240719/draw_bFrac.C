#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"

using namespace std;

void draw_bFrac()
{
	gStyle->SetOptStat(0);
    setTDRStyle();

	const int nFiles = 4;

	TString SysName[nFiles] = {"Err", "Res", "Bkg","True"};

	TFile *in_pt[nFiles];
	TFile *in_cent[nFiles];
	in_pt[0] = TFile::Open("syst_roots/syst_pt_ctauERR.root");
	in_pt[1] = TFile::Open("syst_roots/syst_pt_ctauRES.root");
	in_pt[2] = TFile::Open("syst_roots/syst_pt_ctauBkg.root");
	in_pt[3] = TFile::Open("syst_roots/syst_pt_ctauTRUE.root");

	in_cent[0] = TFile::Open("syst_roots/syst_cent_ctauERR.root");
	in_cent[1] = TFile::Open("syst_roots/syst_cent_ctauRES.root");
	in_cent[2] = TFile::Open("syst_roots/syst_cent_ctauBkg.root");
	in_cent[3] = TFile::Open("syst_roots/syst_cent_ctauTRUE.root");

	TH1D *h_mid_pt_PR[nFiles];
	TH1D *h_mid_pt_NP[nFiles];
	TH1D *h_fwd_pt_PR[nFiles];
	TH1D *h_fwd_pt_NP[nFiles];
	TH1D *h_mid_cent_PR[nFiles];
	TH1D *h_mid_cent_NP[nFiles];
	TH1D *h_fwd_cent_PR[nFiles];
	TH1D *h_fwd_cent_NP[nFiles];

	for (int i=0; i<nFiles; i++) {
		h_mid_pt_PR[i] = (TH1D*) in_pt[i]->Get("mid_PR");
		h_mid_pt_NP[i] = (TH1D*) in_pt[i]->Get("mid_NP");
		h_fwd_pt_PR[i] = (TH1D*) in_pt[i]->Get("fwd_PR");
		h_fwd_pt_NP[i] = (TH1D*) in_pt[i]->Get("fwd_NP");
		h_mid_cent_PR[i] = (TH1D*) in_cent[i]->Get("mid_PR");
		h_mid_cent_NP[i] = (TH1D*) in_cent[i]->Get("mid_NP");
		h_fwd_cent_PR[i] = (TH1D*) in_cent[i]->Get("fwd_PR");
		h_fwd_cent_NP[i] = (TH1D*) in_cent[i]->Get("fwd_NP");
		//cout << h_mid_pt_PR[i]->GetBinContent(1) << endl;
		//cout << h_mid_pt_NP[i]->GetBinContent(1) << endl;
	}

	TCanvas *c_mid_pt_PR = new TCanvas("c_mid_pt_PR","",700,700);
	TCanvas *c_mid_pt_NP = new TCanvas("c_mid_pt_NP","",700,700);
	TCanvas *c_fwd_pt_PR = new TCanvas("c_fwd_pt_PR","",700,700);
	TCanvas *c_fwd_pt_NP = new TCanvas("c_fwd_pt_NP","",700,700);
	TCanvas *c_mid_cent_PR = new TCanvas("c_mid_cent_PR","",700,700);
	TCanvas *c_mid_cent_NP = new TCanvas("c_mid_cent_NP","",700,700);
	TCanvas *c_fwd_cent_PR = new TCanvas("c_fwd_cent_PR","",700,700);
	TCanvas *c_fwd_cent_NP = new TCanvas("c_fwd_cent_NP","",700,700);

	TLegend *leg_midptPR = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_midptPR);
	TLegend *leg_midptNP = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_midptNP);
	TLegend *leg_fwdptPR = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_fwdptPR);
	TLegend *leg_fwdptNP = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_fwdptNP);

	TLegend *leg_midcentPR = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_midcentPR);
	TLegend *leg_midcentNP = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_midcentNP);
	TLegend *leg_fwdcentPR = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_fwdcentPR);
	TLegend *leg_fwdcentNP = new TLegend(0.78,0.67,0.9,0.87);
    SetLegendStyle(leg_fwdcentNP);

	c_mid_pt_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_pt_PR[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_mid_pt_PR[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_PR[i]->SetLineColor(i+2);
		h_mid_pt_PR[i]->SetLineWidth(2);
		h_mid_pt_PR[i]->Draw("SAME");

		leg_midptPR->AddEntry(h_mid_pt_PR[i],SysName[i]);
		leg_midptPR->Draw("SAME");
	}
	c_mid_pt_PR->SaveAs("./figs/bfrac/mid_pt_PR.pdf");

	c_mid_pt_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_pt_NP[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_mid_pt_NP[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_mid_pt_NP[i]->SetLineColor(i+2);
		h_mid_pt_NP[i]->SetLineWidth(2);
		h_mid_pt_NP[i]->Draw("SAME");

		leg_midptNP->AddEntry(h_mid_pt_NP[i],SysName[i]);
		leg_midptNP->Draw("SAME");
	}	
	c_mid_pt_NP->SaveAs("./figs/bfrac/mid_pt_NP.pdf");

	c_fwd_pt_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_pt_PR[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_fwd_pt_PR[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_PR[i]->SetLineColor(i+2);
		h_fwd_pt_PR[i]->SetLineWidth(2);
		h_fwd_pt_PR[i]->Draw("SAME");

		leg_fwdptPR->AddEntry(h_fwd_pt_PR[i],SysName[i]);
		leg_fwdptPR->Draw("SAME");
	}
	c_fwd_pt_PR->SaveAs("./figs/bfrac/fwd_pt_PR.pdf");

	c_fwd_pt_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_pt_NP[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_fwd_pt_NP[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
		h_fwd_pt_NP[i]->SetLineColor(i+2);
		h_fwd_pt_NP[i]->SetLineWidth(2);
		h_fwd_pt_NP[i]->Draw("SAME");

		leg_fwdptNP->AddEntry(h_fwd_pt_NP[i],SysName[i]);
		leg_fwdptNP->Draw("SAME");
	}
	c_fwd_pt_NP->SaveAs("./figs/bfrac/fwd_pt_NP.pdf");

	c_mid_cent_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_cent_PR[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_mid_cent_PR[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_PR[i]->SetLineColor(i+2);
		h_mid_cent_PR[i]->SetLineWidth(2);
		h_mid_cent_PR[i]->Draw("SAME");

		leg_midcentPR->AddEntry(h_mid_cent_PR[i],SysName[i]);
		leg_midcentPR->Draw("SAME");
	}
	c_mid_cent_PR->SaveAs("./figs/bfrac/mid_cent_PR.pdf");

	c_mid_cent_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_mid_cent_NP[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_mid_cent_NP[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_mid_cent_NP[i]->SetLineColor(i+2);
		h_mid_cent_NP[i]->SetLineWidth(2);
		h_mid_cent_NP[i]->Draw("SAME");

		leg_midcentNP->AddEntry(h_mid_cent_NP[i],SysName[i]);
		leg_midcentNP->Draw("SAME");
	}
	c_mid_cent_NP->SaveAs("./figs/bfrac/mid_cent_NP.pdf");

	c_fwd_cent_PR->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_cent_PR[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_fwd_cent_PR[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_PR[i]->SetLineColor(i+2);
		h_fwd_cent_PR[i]->SetLineWidth(2);
		h_fwd_cent_PR[i]->Draw("SAME");

		leg_fwdcentPR->AddEntry(h_fwd_cent_PR[i],SysName[i]);
		leg_fwdcentPR->Draw("SAME");
	}
	c_fwd_cent_PR->SaveAs("./figs/bfrac/fwd_cent_PR.pdf");

	c_fwd_cent_NP->cd();
	for (int i=0; i<nFiles; i++) {
		h_fwd_cent_NP[i]->GetYaxis()->SetRangeUser(0,0.2);
		h_fwd_cent_NP[i]->GetXaxis()->SetTitle("Centrality (%)");
		h_fwd_cent_NP[i]->SetLineColor(i+2);
		h_fwd_cent_NP[i]->SetLineWidth(2);
		h_fwd_cent_NP[i]->Draw("SAME");

		leg_fwdcentNP->AddEntry(h_fwd_cent_NP[i],SysName[i]);
		leg_fwdcentNP->Draw("SAME");
	}
	c_fwd_cent_NP->SaveAs("./figs/bfrac/fwd_cent_NP.pdf");




}
