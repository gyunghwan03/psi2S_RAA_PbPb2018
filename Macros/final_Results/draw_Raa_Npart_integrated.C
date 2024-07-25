#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
#include "../../Style.h"
#include <RooGaussian.h>
#include <RooFormulaVar.h>
#include <RooCBShape.h>
#include <RooWorkspace.h>
#include <RooChebychev.h>
#include <RooPolynomial.h>
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

void draw_Raa_Npart_integrated(bool isSys = true, double ptHigh = 40)
{
	gStyle->SetOptStat(0);
    setTDRStyle();
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

    TFile *f_mid = new TFile(Form("roots/RAA_psi2S_midRap_Npart_pt6p5_%.f.root",ptHigh));
    TFile *f_fwd = new TFile(Form("roots/RAA_psi2S_forRap_Npart_pt3p5_%.f_4Bins.root",ptHigh));
    TFile *fSys = new TFile("../syst_summary/syst_roots/total_syst.root");

	TH1D *h_midPR = (TH1D*) f_mid->Get("hRAA_PR");
    TH1D *h_fwdPR = (TH1D*) f_fwd->Get("hRAA_PR");

    TH1D *hSys_mid_PR = (TH1D*) fSys->Get("mid_cent_PR");
	TH1D *hSys_fwd_PR = (TH1D*) fSys->Get("fwd_cent_PR");

	TH1D *h_midNP = (TH1D*) f_mid->Get("hRAA_NP");
    TH1D *h_fwdNP = (TH1D*) f_fwd->Get("hRAA_NP");

    TH1D *hSys_mid_NP = (TH1D*) fSys->Get("mid_cent_NP");
	TH1D *hSys_fwd_NP = (TH1D*) fSys->Get("fwd_cent_NP");

	const int nCentBins_mid=8;
    const int nCentBins_fwd=4;
	double NpartBin_mid[nCentBins_mid] = {27.12,87.19,131.0,188.2,241.0,283.6,331.5,382.3};
    double NpartBin_fwd[nCentBins_fwd] = {27.12,109.1,225.2,356.9};
	double binWidth_mid[nCentBins_mid]={4.3,4.4,4.3,4.3,4.3,4.3,4.3,4.3};
	double binWidth_fwd[nCentBins_fwd]={4.3,4.3,4.3,4.3};

	double midPR_new[nCentBins_mid];
    double midPR_new_Err[nCentBins_mid];
    double fwdPR_new[nCentBins_fwd];
    double fwdPR_new_Err[nCentBins_fwd];
	double midNP_new[nCentBins_mid];
    double midNP_new_Err[nCentBins_mid];
    double fwdNP_new[nCentBins_fwd];
    double fwdNP_new_Err[nCentBins_fwd];

	double midPR_Sys[nCentBins_mid];
	double midNP_Sys[nCentBins_mid];
	double fwdPR_Sys[nCentBins_fwd];
	double fwdNP_Sys[nCentBins_fwd];

	for (int i=0; i<nCentBins_mid; i++){
        midPR_new[i] = h_midPR->GetBinContent(i+1);
        midPR_new_Err[i] = h_midPR->GetBinError(i+1);
        midPR_Sys[i] = midPR_new[i]*(hSys_mid_PR->GetBinContent(i+1));
        midNP_new[i] = h_midNP->GetBinContent(i+1);
        midNP_new_Err[i] = h_midNP->GetBinError(i+1);
        midNP_Sys[i] = midNP_new[i]*(hSys_mid_NP->GetBinContent(i+1));
        //    cout << "New Val : " << midPR_new[i]  << endl;
    }
	for (int i=0; i<nCentBins_fwd; i++){
        fwdPR_new[i] = h_fwdPR->GetBinContent(i+1);
        fwdPR_new_Err[i] = h_fwdPR->GetBinError(i+1);
        fwdPR_Sys[i] = fwdPR_new[i]*(hSys_fwd_PR->GetBinContent(i+1));
        fwdNP_new[i] = h_fwdNP->GetBinContent(i+1);
        fwdNP_new_Err[i] = h_fwdNP->GetBinError(i+1);
        fwdNP_Sys[i] = fwdNP_new[i]*(hSys_fwd_NP->GetBinContent(i+1));
        //    cout << "New Val : " << fwdPR_new[i]  << endl;
    }

	double midPR_int = 0.257; double midPR_intErr = 0.0544;
	double midNP_int = 0.247; double midNP_intErr = 0.0522;
	double fwdPR_int = 0.276; double fwdPR_intErr = 0.0294;
	double fwdNP_int = 0.238; double fwdNP_intErr = 0.0254;

	TH1D *hmidPR_int = new TH1D("hmidPR_int", "", 3,0,1);
	TH1D *hmidNP_int = new TH1D("hmidNP_int", "", 3,0,1);

	hmidPR_int->SetBinContent(2, midPR_int);
	hmidPR_int->SetBinError(2, midPR_intErr);
	hmidPR_int->SetMarkerStyle(20);
	hmidPR_int->SetMarkerSize(1.5);
	hmidPR_int->SetLineColor(kBlue+2);
	hmidPR_int->SetMarkerColor(kBlue+2);

	hmidNP_int->SetBinContent(2, midNP_int);
	hmidNP_int->SetBinError(2, midNP_intErr);
	hmidNP_int->SetMarkerStyle(21);
	hmidNP_int->SetMarkerSize(1.5);
	hmidNP_int->SetLineColor(kRed+3);
	hmidNP_int->SetMarkerColor(kRed+3);

	hmidPR_int->GetYaxis()->SetRangeUser(0,1.44);
	hmidNP_int->GetYaxis()->SetRangeUser(0,1.44);
	hmidPR_int->GetYaxis()->SetLabelOffset(0);
	hmidNP_int->GetYaxis()->SetLabelOffset(0);
	hmidPR_int->GetXaxis()->SetLabelOffset(0);
	hmidNP_int->GetXaxis()->SetLabelOffset(0);
	hmidPR_int->GetXaxis()->SetLabelSize(0);
	hmidNP_int->GetXaxis()->SetLabelSize(0);
	hmidPR_int->GetXaxis()->SetTickLength(0);
	hmidNP_int->GetXaxis()->SetTickLength(0);
	hmidPR_int->GetYaxis()->SetLabelSize(0.12);
	hmidNP_int->GetYaxis()->SetLabelSize(0.12);

	TH1D *hfwdPR_int = new TH1D("hfwdPR_int", "", 3,0,1);
	TH1D *hfwdNP_int = new TH1D("hfwdNP_int", "", 3,0,1);

	hfwdPR_int->SetBinContent(2, fwdPR_int);
	hfwdPR_int->SetBinError(2, fwdPR_intErr);
	hfwdPR_int->SetMarkerStyle(20);
	hfwdPR_int->SetMarkerSize(1.5);
	hfwdPR_int->SetLineColor(kBlue+2);
	hfwdPR_int->SetMarkerColor(kBlue+2);

	hfwdNP_int->SetBinContent(2, fwdNP_int);
	hfwdNP_int->SetBinError(2, fwdNP_intErr);
	hfwdNP_int->SetMarkerStyle(21);
	hfwdNP_int->SetMarkerSize(1.5);
	hfwdNP_int->SetLineColor(kRed+3);
	hfwdNP_int->SetMarkerColor(kRed+3);

	hfwdPR_int->GetYaxis()->SetRangeUser(0,1.44);
	hfwdNP_int->GetYaxis()->SetRangeUser(0,1.44);
	hfwdPR_int->GetYaxis()->SetLabelOffset(0);
	hfwdNP_int->GetYaxis()->SetLabelOffset(0);
	hfwdPR_int->GetXaxis()->SetLabelOffset(0);
	hfwdNP_int->GetXaxis()->SetLabelOffset(0);
	hfwdPR_int->GetXaxis()->SetLabelSize(0);
	hfwdNP_int->GetXaxis()->SetLabelSize(0);
	hfwdPR_int->GetXaxis()->SetTickLength(0);
	hfwdNP_int->GetXaxis()->SetTickLength(0);
	hfwdPR_int->GetYaxis()->SetLabelSize(0.12);
	hfwdNP_int->GetYaxis()->SetLabelSize(0.12);

	TGraphErrors *g_midPR;
	TGraphErrors *g_midNP;
	TGraphErrors *g_fwdPR;
	TGraphErrors *g_fwdNP;
	g_midPR = new TGraphErrors(nCentBins_mid, NpartBin_mid, midPR_new, 0, midPR_new_Err);
	g_midNP = new TGraphErrors(nCentBins_mid, NpartBin_mid, midNP_new, 0, midNP_new_Err);
	g_fwdPR = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, fwdPR_new, 0, fwdPR_new_Err);
	g_fwdNP = new TGraphErrors(nCentBins_fwd, NpartBin_fwd, fwdNP_new, 0, fwdNP_new_Err);
	//g_midPR_int = new TGraphErrors(1,1,

	//TGraphErrors *g_midPRSys = new TGraphErrors(nCentBins_mid,NpartBin_mid,midPR_new,binWidth,midPR_Sys);
	//TGraphErrors *g_midNPSys = new TGraphErrors(nCentBins_mid,NpartBin_mid,midNP_new,binWidth,midNP_Sys);
    //TGraphErrors *g_fwdPRSys = new TGraphErrors(nCentBins_fwd,NpartBin_fwd,fwdPR_new,binWidth,fwdPR_Sys);
    //TGraphErrors *g_fwdNPSys = new TGraphErrors(nCentBins_fwd,NpartBin_fwd,fwdNP_new,binWidth,fwdNP_Sys);

	g_midPR->GetXaxis()->SetTitle("<N_{Part}>");
    g_midPR->GetXaxis()->CenterTitle();
    g_midPR->GetYaxis()->SetTitle("R_{AA}");
    g_midPR->SetTitle("");
    g_midPR->GetYaxis()->CenterTitle();
    g_midPR->GetXaxis()->SetLimits(0.,400);
    g_midPR->SetMinimum(0);
    g_midPR->SetMaximum(1.44);
    g_midNP->SetMinimum(0);
    g_midNP->SetMaximum(1.44);

    g_midPR->SetMarkerColor(kBlue+2);
    g_midPR->SetLineColor(kBlue+2);
    g_midPR->SetMarkerStyle(20);
    g_midPR->SetMarkerSize(1.5);
    g_midNP->SetMarkerColor(kRed+3);
    g_midNP->SetLineColor(kRed+3);
    g_midNP->SetMarkerStyle(21);
    g_midNP->SetMarkerSize(1.5);

    //g_midPRSys->SetMinimum(0);
    //g_midPRSys->SetMaximum(1.44);
    //g_midPRSys->SetLineColor(kBlue-4);
    //g_midPRSys->SetFillColorAlpha(kBlue-9,0.40);
    //g_midNPSys->SetLineColor(kRed-4);
    //g_midNPSys->SetFillColorAlpha(kRed-9,0.40);

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
	
    g_fwdPR->SetMarkerColor(kBlue+2);
    g_fwdPR->SetLineColor(kBlue+2);
    g_fwdPR->SetMarkerStyle(20);
    g_fwdPR->SetMarkerSize(1.5);
    g_fwdNP->SetMarkerColor(kRed+3);
    g_fwdNP->SetLineColor(kRed+3);
    g_fwdNP->SetMarkerStyle(21);
    g_fwdNP->SetMarkerSize(1.5);

	TLegend *leg = new TLegend(0.68,0.72,0.8,0.82);
    //SetLegendStyle(leg);
	leg->SetFillColor(0);
	//leg->SetFillStyle(4000);
	leg->SetBorderSize(0);
	leg->SetMargin(0.2);
	leg->SetTextSize(0.029);
	leg->SetTextFont(42);
    leg->AddEntry(g_midPR,"Prompt #psi(2S)", "p");
    leg->AddEntry(g_midNP,"Non Prompt #psi(2S)", "p");

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
	
	g_midPR->Draw("AP");
	g_midNP->Draw("P");
    leg->Draw("SAME");

    jumSun(0,1,400,1);	
	drawText(Form("6.5 < p_{T} < %.f GeV/c",ptHigh), pos_x, pos_y, text_color, text_size);
    drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(pad1_1,iPeriod,iPos);


	pad1_2->cd();
	hmidPR_int->Draw("PE");
	hmidNP_int->Draw("PE SAME");
    drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff, text_color, text_size);
	jumSun(0,1,1,1);

	c1->SaveAs(Form("figs/Integrated_RAA_psi2S_y0_1p6_Npart_Sys%d_pt6p5_%.f.pdf",isSys,ptHigh));


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
	
	g_fwdPR->Draw("AP");
	g_fwdNP->Draw("P");
    leg->Draw("SAME");

    jumSun(0,1,400,1);	
	drawText(Form("3.5 < p_{T} < %.f GeV/c",ptHigh), pos_x, pos_y, text_color, text_size);
    drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
	CMS_lumi_v2mass(pad2_1,iPeriod,iPos);


	pad2_2->cd();
	hfwdPR_int->Draw("PE");
	hfwdNP_int->Draw("PE SAME");
    drawText("Cent. 0-90%", pos_x, pos_y-pos_y_diff, text_color, text_size);
	jumSun(0,1,1,1);

	c2->SaveAs(Form("figs/Integrated_RAA_psi2S_y1p6_2p4_Npart_Sys%d_pt3p5_%.f.pdf",isSys,ptHigh));
}
