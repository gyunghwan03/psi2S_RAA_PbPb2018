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

valErr getYield_PbPb(int i=0);
valErr getYield_pp(int i=0);
valErr getFrac_PbPb(int i=0);
valErr getFrac_pp(int i=0);

void draw_Raa_psi2S_y0_1p6_pT()
{

	//gROOT->SumW2();
	gStyle->SetOptStat(0);
	setTDRStyle();
	int iPeriod = 101;
  	int iPos = 33;

    const int nPtBins=6;
    TFile *fPbPb[nPtBins+1];
    TFile *fpp[nPtBins+1];

	TFile *fEff_PbPbPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20230729.root");
    TFile *fEff_PbPbNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_psi2S_PtW1_tnp1_20230729.root");
    TFile *fEff_ppPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230728.root");
    TFile *fEff_ppNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_nprompt_pp_psi2s_PtW1_tnp1_20230801.root");
    TFile *fAcc_pp = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20230728.root");
    TFile *fAcc_PbPb = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_PbPb_SysUp0_20230728.root");

	TH1D *hEff_PbPbPR = (TH1D*) fEff_PbPbPR -> Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6");
    TH1D *hEff_PbPbNP = (TH1D*) fEff_PbPbNP -> Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy0_1p6");
    TH1D *hEff_ppPR = (TH1D*) fEff_ppPR -> Get("mc_eff_vs_pt_TnP1_PtW1_absy0_1p6");
    TH1D *hEff_ppNP = (TH1D*) fEff_ppNP -> Get("mc_eff_vs_pt_TnP1_PtW1_absy0_1p6");
    TH1D *hAcc_PbPb = (TH1D*) fAcc_PbPb -> Get("hAccPt_2021_midy");
    TH1D *hAcc_pp = (TH1D*) fAcc_pp -> Get("hAccPt_2021_midy");

	Double_t Nmb = 11968044281.;
	Double_t Nmb_err = Nmb*0.01261;
	//Double_t Taa = 5.649; // 0-100%
	Double_t Taa = 6.274; // 0-90%
	Double_t lumi_pp = 3.;
	Double_t lumi_pp_scale = 1e-9;
	Double_t lumi_pp_err = 1*lumi_pp_scale;
	//Double_t Taa_err = 0.123; //0-100%
	Double_t Taa_err = 0.137; //0-90%

    double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,50};
    double fracPP[nPtBins]; double fracPbPb[nPtBins];
    double fracErrPP[nPtBins]; double fracErrPbPb[nPtBins];

    TString kineLabel[nPtBins+1];

    TH1D *hyieldPP_PR = new TH1D("hyieldPP_PR",";p_{T} (GeV/c);",nPtBins, ptBin);
    TH1D *hyieldPbPb_PR = new TH1D("hyieldPbPb_PR",";p_{T} (GeV/c);",nPtBins, ptBin);
    TH1D *hyieldPP_NP = new TH1D("hyieldPP_NP",";p_{T} (GeV/c);",nPtBins, ptBin);
    TH1D *hyieldPbPb_NP = new TH1D("hyieldPbPb_NP",";p_{T} (GeV/c);",nPtBins, ptBin);

	TH1D *hXpp_PR = new TH1D("hXpp_PR",";p_{T} (GeV/c);",nPtBins, ptBin);
	TH1D *hXpp_NP = new TH1D("hXpp_NP",";p_{T} (GeV/c);",nPtBins, ptBin);
	TH1D *hXPbPb_PR = new TH1D("hXPbPb_PR",";p_{T} (GeV/c);",nPtBins, ptBin);
	TH1D *hXPbPb_NP = new TH1D("hXPbPb_NP",";p_{T} (GeV/c);",nPtBins, ptBin);

	TH1D *hRAA_PR = new TH1D("hRAA_PR",";p_{T} (GeV/c);", nPtBins,ptBin);
	TH1D *hRAA_NP = new TH1D("hRAA_NP",";p_{T} (GeV/c);", nPtBins,ptBin);

	double RaaPR[nPtBins]; double RaaPR_err[nPtBins]; double binWidth[nPtBins]; double x[nPtBins];
	double RaaNP[nPtBins]; double RaaNP_err[nPtBins];

	double Xpp_PR[nPtBins]; double Xpp_NP[nPtBins];
	double Xpp_PR_err[nPtBins]; double Xpp_NP_err[nPtBins];
	double XPbPb_PR[nPtBins]; double XPbPb_NP[nPtBins];
	double XPbPb_PR_err[nPtBins]; double XPbPb_NP_err[nPtBins];


	double Taa_Nmb; double step_one[nPtBins]; double step_two[nPtBins]; double RAA[nPtBins];

	for(int i=0;i<nPtBins;i++)
    {
		double weight_ppPR = 1; double weight_ppNP = 1; double weight_PbPbPR = 1; double weight_PbPbNP = 1;
        double eff_ppPR = 1; double eff_ppNP = 1; double eff_PbPbPR = 1; double eff_PbPbNP = 1;
        double acc_pp = 1; double acc_PbPb = 1;


        eff_ppPR=hEff_ppPR->GetBinContent(i+2);
        eff_ppNP=hEff_ppNP->GetBinContent(i+2);
        acc_pp=hAcc_pp->GetBinContent(i+1);
        eff_PbPbPR=hEff_PbPbPR->GetBinContent(i+2);
        eff_PbPbNP=hEff_PbPbNP->GetBinContent(i+2);
        acc_PbPb=hAcc_PbPb->GetBinContent(i+1);

        weight_ppPR=eff_ppPR*acc_pp;
        weight_ppNP=eff_ppNP*acc_pp;
		weight_PbPbPR=eff_PbPbPR*acc_PbPb;
        weight_PbPbNP=eff_PbPbNP*acc_PbPb;

        valErr yieldPP; valErr yieldPbPb; valErr fracPP; valErr fracPbPb;
        yieldPP = getYield_pp(i);
        yieldPbPb = getYield_PbPb(i);
        fracPP = getFrac_pp(i);
        fracPbPb = getFrac_PbPb(i);

        double err1_PP = yieldPP.err;
        double err2_PP = fracPP.err;
        double err_PPNP = yieldPP.val*fracPP.val*sqrt((err1_PP/yieldPP.val)*(err1_PP/yieldPP.val) + (err2_PP/fracPP.val)*(err2_PP/fracPP.val));
        double err_PPPR = yieldPP.val*(1-fracPP.val)*sqrt((err1_PP/yieldPP.val)*(err1_PP/yieldPP.val) + ((err2_PP/(1-fracPP.val))*(err2_PP/(1-fracPP.val))));

        double err1_PbPb = yieldPbPb.err;
        double err2_PbPb = fracPbPb.err;
        double err_PbPbNP = yieldPbPb.val*fracPbPb.val*sqrt((err1_PbPb/yieldPbPb.val)*(err1_PbPb/yieldPbPb.val) + (err2_PbPb/fracPbPb.val)*(err2_PbPb/fracPbPb.val));
        double err_PbPbPR = yieldPbPb.val*(1-fracPbPb.val)*sqrt((err1_PbPb/yieldPbPb.val)*(err1_PbPb/yieldPbPb.val) + ((err2_PbPb/(1-fracPbPb.val))*(err2_PbPb/(1-fracPbPb.val))));
		
		double yieldPP_PR = yieldPP.val*(1-fracPP.val)/weight_ppPR;
		double yieldPP_NP = yieldPP.val*(fracPP.val)/weight_ppNP;
		double yieldPbPb_PR = yieldPbPb.val*(1-fracPbPb.val)/weight_PbPbPR;
		double yieldPbPb_NP = yieldPbPb.val*(fracPbPb.val)/weight_PbPbNP;
		
		double err_PbPbPR_wgt = err_PbPbPR/weight_PbPbPR;
		double err_PbPbNP_wgt = err_PbPbNP/weight_PbPbNP;
		double err_PPPR_wgt = err_PPPR/weight_ppPR;
		double err_PPNP_wgt = err_PPNP/weight_ppNP;

		cout << " " << endl;
		cout << "pt " << ptBin[i] << " - " << ptBin[i+1] << endl;
		cout << "PbPb Yield	: " << yieldPbPb.val << ",	b frac	: " << fracPbPb.val << ",	PR Eff	: " << eff_PbPbPR << ",	NP Eff	: " << eff_PbPbNP<< ",	Acc	: " << acc_PbPb << endl;
		cout << "pp Yield	: " << yieldPP.val << ",	b frac	: " << fracPP.val << ",	PR Eff	: " << eff_ppPR << ",	NP Eff	: " << eff_ppNP << ",	Acc	: " << acc_pp << endl;
        cout << "pp Prompt yield : " << yieldPP_PR << ", PbPb Prompt yield : " << yieldPbPb_PR << ", pp NonPrompt yield : " << yieldPP_NP << ", PbPb NonPrompt yield : " << yieldPbPb_NP << endl;
		cout << "pp Prompt error : " << err_PPPR_wgt << " , PbPb Prompt error : " << err_PbPbPR_wgt << " , pp NonPrompt error : " << err_PPNP_wgt << " , PbPb NonPrompt error : " << err_PbPbNP_wgt << endl;

		hyieldPP_PR -> SetBinContent(i+1,yieldPP.val*(1-fracPP.val));
		hyieldPP_PR -> SetBinError(i+1,err_PPPR_wgt);
		hyieldPbPb_PR -> SetBinContent(i+1,yieldPbPb.val*(1-fracPbPb.val));
		hyieldPbPb_PR -> SetBinError(i+1,err_PbPbPR_wgt);

		hyieldPP_NP -> SetBinContent(i+1,yieldPP.val*(fracPP.val));
		hyieldPP_NP -> SetBinError(i+1,err_PPNP_wgt);
		hyieldPbPb_NP -> SetBinContent(i+1,yieldPbPb.val*(fracPbPb.val));
		hyieldPbPb_NP -> SetBinError(i+1,err_PbPbNP_wgt);

		Xpp_PR[i] = lumi_pp_scale*yieldPP_PR/(lumi_pp*1e+2*(ptBin[i+1]-ptBin[i])*(double)2*(1.6));
		Xpp_PR_err[i] = Xpp_PR[i]*sqrt(TMath::Power(err_PPPR_wgt/yieldPP_PR,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_NP[i] = lumi_pp_scale*yieldPP_NP/(lumi_pp*1e+2*(ptBin[i+1]-ptBin[i])*(double)2*(1.6));
		Xpp_NP_err[i] = lumi_pp_scale*Xpp_NP[i]*sqrt(TMath::Power(err_PPNP_wgt/yieldPP_NP,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));

		XPbPb_PR[i] = yieldPbPb_PR/(Nmb*Taa*(ptBin[i+1]-ptBin[i])*(double)2*(1.6));
		XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(Taa_err/Taa,2) + TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2) + TMath::Power(Nmb_err/Nmb,2));
		XPbPb_NP[i] = yieldPbPb_NP/(Nmb*Taa*(ptBin[i+1]-ptBin[i])*(double)2*(1.6));
		XPbPb_NP_err[i] = XPbPb_NP[i]*sqrt(TMath::Power(Taa_err/Taa,2) + TMath::Power(err_PbPbNP_wgt/yieldPbPb_NP,2) + TMath::Power(Nmb_err/Nmb,2));

		hXpp_PR->SetBinContent(i+1, Xpp_PR[i]);
		hXpp_NP->SetBinContent(i+1, Xpp_NP[i]);
		hXPbPb_PR->SetBinContent(i+1, XPbPb_PR[i]);
		hXPbPb_NP->SetBinContent(i+1, XPbPb_NP[i]);

		hXpp_PR->SetBinError(i+1, Xpp_PR_err[i]);
		hXpp_NP->SetBinError(i+1, Xpp_NP_err[i]);
		hXPbPb_PR->SetBinError(i+1, XPbPb_PR_err[i]);
		hXPbPb_NP->SetBinError(i+1, XPbPb_NP_err[i]);


		//cout << "yield Err NonPrompt : " << err_PbPbNP << "	X-Section Err NP : " << XPbPb_NP_err[i] << endl;

        Taa_Nmb =  Taa*Nmb*2*(2.4-1.6)*(50.-3.);
        step_one[i] = 1./(Taa_Nmb*Xpp_PR[i]);
        step_two[i] = step_one[i];
        RAA[i] = step_two[i]*yieldPbPb_PR;
        //cout << "Taa_Nmb = " << Taa_Nmb << " ,step one = " << step_one[i] << " , step two = " << step_two[i] << ", RAA = " << RAA[i] << endl;
		 cout << "Xpp_PR : " << Xpp_PR[i] << " , Xpp_PR_err : " << Xpp_PR_err[i] << ", XPbPb_PR : " << XPbPb_PR[i] << " , XPbPb_PR_err : " << XPbPb_PR_err[i] << " , Xpp_NP : " << Xpp_NP[i] << " , XPbPb_NP : " << XPbPb_NP[i] << endl;
	}

	hXpp_PR->Sumw2();
	hXpp_NP->Sumw2();
	hXPbPb_PR->Sumw2();
	hXPbPb_NP->Sumw2();


	TH1D *hRaa_pp_PR = (TH1D*)hXpp_PR->Clone("hRaa_pp_PR");
	TH1D *hRaa_PbPb_PR = (TH1D*)hXPbPb_PR->Clone("hRaa_PbPb_PR");
	TH1D *hRaa_pp_NP = (TH1D*)hXpp_NP->Clone("hRaa_pp_NP");
	TH1D *hRaa_PbPb_NP = (TH1D*)hXPbPb_NP->Clone("hRaa_PbPb_NP");

	TCanvas *c2 = new TCanvas("c2","",700,700);
	c2->cd();
	hRaa_PbPb_PR->Divide(hRaa_pp_PR);
	hRaa_PbPb_PR->Sumw2();
	hRaa_PbPb_NP->Divide(hRaa_pp_NP);
	hRaa_PbPb_NP->Sumw2();
	hRaa_PbPb_PR->GetYaxis()->SetRangeUser(0,1.5);
	hRaa_PbPb_PR->SetLineColor(kRed+2);
	hRaa_PbPb_PR->SetMarkerColor(kRed+2);
	hRaa_PbPb_PR->SetMarkerStyle(25);
	hRaa_PbPb_PR->Draw("PE");

	for(int i=0;i<nPtBins;i++)
	{
		RaaPR[i] = hRaa_PbPb_PR->GetBinContent(i+1);
		RaaPR_err[i] = hRaa_PbPb_PR->GetBinError(i+1);
		RaaNP[i] = hRaa_PbPb_NP->GetBinContent(i+1);
		RaaNP_err[i] = hRaa_PbPb_NP->GetBinError(i+1);
		x[i] = (ptBin[i+1]+ptBin[i])/2;
		binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		cout << "pt : " << ptBin[i] << " - " << ptBin[i+1] << ", RAA Value : " << RaaPR[i] << ",	RAA Err : " << RaaPR_err[i] << " x : " << x[i] << ", bin Width : " << binWidth[i] << endl;
		cout << "		, RAA NP : " << RaaNP[i] <<  ",	RaaNP_err : " << RaaNP_err[i] << endl;
		//gRaaNP->SetPoint(i,(ptBin[i+1]-ptBin[i])/2,); 
		//gRaaNP->SetPointError(i,0,Xpp_NP_err); 
		hRAA_PR->SetBinContent(i+1,RaaPR[i]);
		hRAA_PR->SetBinError(i+1,RaaPR_err[i]);
		hRAA_NP->SetBinContent(i+1,RaaNP[i]);
		hRAA_NP->SetBinError(i+1,RaaNP_err[i]);
	}

	//////////////////////////////// Start Plotting ////////////////////////////////
	float pos_x = 0.25;
	float pos_y = 0.85;
	float pos_y_diff = 0.071;
	int text_color = 1;
	float text_size = 25;
	TGraphErrors *gXpp_PR = new TGraphErrors(nPtBins,x,Xpp_PR,binWidth,Xpp_PR_err);
	TGraphErrors *gXpp_NP = new TGraphErrors(nPtBins,x,Xpp_NP,binWidth,Xpp_NP_err);

	TGraphErrors *gXPbPb_PR = new TGraphErrors(nPtBins,x,XPbPb_PR,binWidth,XPbPb_PR_err);
	TGraphErrors *gXPbPb_NP = new TGraphErrors(nPtBins,x,XPbPb_NP,binWidth,XPbPb_NP_err);

	TCanvas *cXPR = new TCanvas("cXPR", "", 700,700);
	cXPR->cd();
	cXPR->SetLogy();
	gXpp_PR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gXpp_PR->GetXaxis()->CenterTitle();
	gXpp_PR->GetYaxis()->SetTitle("#it{B} #times d#sigma/dp_{T} or #it{B} #times (1/T_{AA}N_{MB})dN/dp_{T} (mb/GeV/c)");
	gXpp_PR->GetYaxis()->SetTitleSize(0.04);
	gXpp_PR->GetYaxis()->SetTitleOffset(1.90);
	gXpp_PR->SetTitle("");
	gXpp_PR->GetYaxis()->CenterTitle();
	gXpp_PR->GetXaxis()->SetLimits(0.,50);
	gXpp_PR->SetMinimum(1e-11);
	gXpp_PR->SetMaximum(1e-3);
	gXpp_PR->SetMarkerColor(kBlue+2);
	gXpp_PR->SetLineColor(kBlue+2);
	gXpp_PR->SetMarkerStyle(20);
	gXpp_PR->SetMarkerSize(1.4);
	gXPbPb_PR->SetMarkerColor(kRed+3);
	gXPbPb_PR->SetLineColor(kRed+3);
	gXPbPb_PR->SetMarkerStyle(21);
	gXPbPb_PR->SetMarkerSize(1.4);

	gXpp_PR->Draw("AP");
	gXPbPb_PR->Draw("P");

	TLegend *legPR = new TLegend(0.68,0.72,0.8,0.82);
    SetLegendStyle(legPR);
	legPR->AddEntry(gXpp_PR,"pp", "p");
	legPR->AddEntry(gXPbPb_PR,"PbPb, Cent. 0-90%", "p");
	legPR->Draw("SAME");
	jumSun(0,1,50,1);

	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXPR,iPeriod,iPos);	

	cXPR->SaveAs("figs/CrossSection_PR_y0_1p6_pT.pdf");

	TCanvas *cXNP = new TCanvas("cXNP", "", 700,700);
	cXNP->cd();
	cXNP->SetLogy();
	gXpp_NP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gXpp_NP->GetXaxis()->CenterTitle();
	gXpp_NP->GetYaxis()->SetTitle("#it{B} #times d#sigma/dp_{T} or #it{B} #times (1/T_{AA}N_{MB})dN/dp_{T} (mb/GeV/c)");
	gXpp_NP->GetYaxis()->SetTitleSize(0.04);
	gXpp_NP->GetYaxis()->SetTitleOffset(1.90);
	gXpp_NP->SetTitle("");
	gXpp_NP->GetYaxis()->CenterTitle();
	gXpp_NP->GetXaxis()->SetLimits(0.,50);
	gXpp_NP->SetMinimum(1e-11);
	gXpp_NP->SetMaximum(1e-3);
	gXpp_NP->SetMarkerColor(kBlue+2);
	gXpp_NP->SetLineColor(kBlue+2);
	gXpp_NP->SetMarkerStyle(20);
	gXpp_NP->SetMarkerSize(1.4);
	gXPbPb_NP->SetMarkerColor(kRed+3);
	gXPbPb_NP->SetLineColor(kRed+3);
	gXPbPb_NP->SetMarkerStyle(21);
	gXPbPb_NP->SetMarkerSize(1.4);

	gXpp_NP->Draw("AP");
	gXPbPb_NP->Draw("P");
	
	TLegend *legNP = new TLegend(0.68,0.72,0.8,0.82);
    SetLegendStyle(legNP);
	legNP->AddEntry(gXpp_NP,"pp", "p");
	legNP->AddEntry(gXPbPb_NP,"PbPb, Cent. 0-90%", "p");
	legNP->Draw("SAME");
	jumSun(0,1,50,1);

	drawText("Non Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXNP,iPeriod,iPos);	

	cXNP->SaveAs("figs/CrossSection_NP_y0_1p6_pT.pdf");

	TGraphErrors *gRaaPR = new TGraphErrors(nPtBins,x,RaaPR,binWidth,RaaPR_err);
	TGraphErrors *gRaaNP = new TGraphErrors(nPtBins,x,RaaNP,binWidth,RaaNP_err);
	TCanvas *cRAA = new TCanvas("cRAA", "", 700, 700);
	cRAA->cd();
	gRaaPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gRaaPR->GetXaxis()->CenterTitle();
	gRaaPR->GetYaxis()->SetTitle("R_{AA}");
	gRaaPR->SetTitle("");
	gRaaPR->GetYaxis()->CenterTitle();
	gRaaPR->GetXaxis()->SetLimits(0.,50);
	gRaaPR->SetMinimum(0);
	gRaaPR->SetMaximum(1.44);

	gRaaPR->SetMarkerColor(kBlue+2);
	gRaaPR->SetLineColor(kBlue+2);
	gRaaPR->SetMarkerStyle(20);
	gRaaPR->SetMarkerSize(1.4);
	gRaaNP->SetMarkerColor(kRed+3);
	gRaaNP->SetLineColor(kRed+3);
	gRaaNP->SetMarkerStyle(21);
	gRaaNP->SetMarkerSize(1.4);

	gRaaPR->Draw("AP");
	gRaaNP->Draw("P");
	TLegend *leg = new TLegend(0.68,0.72,0.8,0.82);
    SetLegendStyle(leg);
	leg->AddEntry(gRaaPR,"Prompt #psi(2S)", "p");
	leg->AddEntry(gRaaNP,"Non Prompt #psi(2S)", "p");
	leg->Draw("SAME");
	jumSun(0,1,50,1);

	drawText("Cent. 0-90%", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 1.6", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cRAA,iPeriod,iPos);	

	cRAA->SaveAs("figs/RAA_psi2S_y0_1p6_pT.pdf");

	TFile *outFile = new TFile("./roots/RAA_psi2S_midRap_pT.root","recreate");
	hRAA_PR->Write();
	hRAA_NP->Write();
	outFile->Close();


}

valErr getYield_pp(int i){
    double ptBins[7] = {6.5,9,12,15,20,25,50};
    TString kineLabel[7];
    kineLabel[i] = getKineLabelpp(ptBins[i],ptBins[i+1],0,1.6,0.0);
    TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getYield_PbPb(int i){
    double ptBins[7] = {6.5,9,12,15,20,25,50};
    TString kineLabel[7];
    kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,1.6,0.0,0,180);
    //TFile* inf = new TFile(Form("./psi2S/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TFile* inf = new TFile(Form("../psi2S_230512/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("MassResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getFrac_PbPb(int i) {
    double ptBin[7] = {6.5,9,12,15,20,25,50};
    TString kineLabel[7];
    kineLabel[i] = getKineLabel(ptBin[i],ptBin[i+1],0,1.6,0.0,0,180);
    TFile* inf = new TFile(Form("../psi2S_230512/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getFrac_pp(int i) {
    double ptBin[7] = {6.5,9,12,15,20,25,50};
    TString kineLabel[7];
    kineLabel[i] = getKineLabelpp(ptBin[i],ptBin[i+1],0,1.6,0.0);
    TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
