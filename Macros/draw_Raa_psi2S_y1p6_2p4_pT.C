#include "../rootFitHeaders.h"
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
#include "../Style.h"
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
double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);

void draw_Raa_psi2S_y1p6_2p4_pT()
{

	//gROOT->SumW2();
	gStyle->SetOptStat(0);
	setTDRStyle();
	int iPeriod = 101;
  	int iPos = 33;

    const int nPtBins=3;
    TFile *fPbPb[nPtBins+1];
    TFile *fpp[nPtBins+1];

	TFile *fEff_PbPb = new TFile("../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20230417.root");
    TFile *fEff_pp = new TFile("../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230416.root");
    TFile *fAcc_PbPb = new TFile("../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_PbPb_SysUp0_20230416.root");
    TFile *fAcc_pp = new TFile("../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20230416.root");

	TH1D *hEff_PbPb = (TH1D*) fEff_PbPb -> Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4");
	TH1D *hEff_pp = (TH1D*) fEff_pp -> Get("mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4");
	TH1D *hAcc_PbPb = (TH1D*) fAcc_PbPb -> Get("hAccPt_2021_Fory");
	TH1D *hAcc_pp = (TH1D*) fAcc_pp -> Get("hAccPt_2021_Fory");

	Double_t Nmb = 11968044281.;
	Double_t Nmb_err = Nmb*0.01261;
	//Double_t Taa = 5.649; // 0-100%
	Double_t Taa = 6.274; // 0-90%
	Double_t lumi_pp = 3.002;
	Double_t lumi_pp_scale = 1e-9;
	Double_t lumi_pp_err = 1*lumi_pp_scale;
	Double_t Taa_err = 0.123;

    double ptBin[nPtBins+1] = {3,6.5,12,50};
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

	double RaaPR[nPtBins]; double RaaPR_err[nPtBins]; double binWidth[nPtBins]; double x[nPtBins];
	double RaaNP[nPtBins]; double RaaNP_err[nPtBins];

	double Xpp_PR[nPtBins]; double Xpp_NP[nPtBins];
	double Xpp_PR_err[nPtBins]; double Xpp_NP_err[nPtBins];
	double XPbPb_PR[nPtBins]; double XPbPb_NP[nPtBins];
	double XPbPb_PR_err[nPtBins]; double XPbPb_NP_err[nPtBins];


	for(int i=0;i<nPtBins;i++)
    {
		double weight_pp = 1; double weight_PbPb = 1;
		double eff_pp = 1; double eff_PbPb = 1;
		double acc_pp = 1; double acc_PbPb = 1;


		eff_pp=hEff_pp->GetBinContent(i+2);
		acc_pp=hAcc_pp->GetBinContent(i+1);
		eff_PbPb=hEff_PbPb->GetBinContent(i+2);
		acc_PbPb=hAcc_PbPb->GetBinContent(i+1);

		weight_pp=eff_pp*acc_pp;
		weight_PbPb=eff_PbPb*acc_PbPb;

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

		double yieldPP_PR = yieldPP.val*(1-fracPP.val)/weight_pp;
		double yieldPP_NP = yieldPP.val*(fracPP.val)/weight_pp;
		double yieldPbPb_PR = yieldPbPb.val*(1-fracPbPb.val)/weight_PbPb;
		double yieldPbPb_NP = yieldPbPb.val*(fracPbPb.val)/weight_PbPb;

		double err_PbPbPR_wgt = err_PbPbPR/weight_PbPb;
        double err_PbPbNP_wgt = err_PbPbNP/weight_PbPb;
        double err_PPPR_wgt = err_PPPR/weight_pp;
        double err_PPNP_wgt = err_PPNP/weight_pp;
		 
		cout << "PbPb Yield : " << yieldPbPb.val << ", yield Err : " << err1_PbPb << ", b frac : " << fracPbPb.val << ", frac Err : " << err2_PbPb << ", Eff : " << eff_PbPb << ", Acc : " << acc_PbPb << endl;
        cout << "pp Yield : " << yieldPP.val << ", yield Err : " << err1_PP << ", b frac : " << fracPP.val << ", frac Err : " << err2_PP << ", Eff : " << eff_pp << ", Acc : " << acc_pp << endl;
        cout << "pp Prompt yield : " << yieldPP_PR << ", PbPb Prompt yield : " << yieldPbPb_PR << ", pp NonPrompt yield : " << yieldPP_NP << ", PbPb NonPrompt yield : " << yieldPbPb_NP << endl;
        cout << "pp Prompt error : " << err_PPPR_wgt << " , PbPb Prompt error : " << err_PbPbPR_wgt << " , pp NonPrompt error : " << err_PPNP_wgt << " , PbPb NonPrompt error : " << err_PbPbNP_wgt << endl;

		hyieldPP_PR -> SetBinContent(i+1,yieldPP.val*(1-fracPP.val));
		hyieldPP_PR -> SetBinError(i+1,err_PPPR);
		hyieldPbPb_PR -> SetBinContent(i+1,yieldPbPb.val*(1-fracPbPb.val));
		hyieldPbPb_PR -> SetBinError(i+1,err_PbPbPR);

		hyieldPP_NP -> SetBinContent(i+1,yieldPP.val*(fracPP.val));
		hyieldPP_NP -> SetBinError(i+1,err_PPNP);
		hyieldPbPb_NP -> SetBinContent(i+1,yieldPbPb.val*(fracPbPb.val));
		hyieldPbPb_NP -> SetBinError(i+1,err_PbPbNP);

		Xpp_PR[i] = lumi_pp_scale*yieldPP_PR/(lumi_pp*1e+2*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
		Xpp_PR_err[i] = Xpp_PR[i]*sqrt(TMath::Power(err_PPPR_wgt/yieldPP_PR,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_NP[i] = lumi_pp_scale*yieldPP_NP/(lumi_pp*1e+2*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
		Xpp_NP_err[i] = Xpp_NP[i]*sqrt(TMath::Power(err_PPNP_wgt/yieldPP_NP,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));

		XPbPb_PR[i] = yieldPbPb_PR/(Nmb*Taa*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
		XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(Taa_err/Taa,2) + TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2) + TMath::Power(Nmb_err/Nmb,2));
		XPbPb_NP[i] = yieldPbPb_NP/(Nmb*Taa*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
		XPbPb_NP_err[i] = XPbPb_NP[i]*sqrt(TMath::Power(Taa_err/Taa,2) + TMath::Power(err_PbPbNP_wgt/yieldPbPb_NP,2) + TMath::Power(Nmb_err/Nmb,2));

		hXpp_PR->SetBinContent(i+1, Xpp_PR[i]);
		hXpp_NP->SetBinContent(i+1, Xpp_NP[i]);
		hXPbPb_PR->SetBinContent(i+1, XPbPb_PR[i]);
		hXPbPb_NP->SetBinContent(i+1, XPbPb_NP[i]);

		hXpp_PR->SetBinError(i+1, Xpp_PR_err[i]);
		hXpp_NP->SetBinError(i+1, Xpp_NP_err[i]);
		hXPbPb_PR->SetBinError(i+1, XPbPb_PR_err[i]);
		hXPbPb_NP->SetBinError(i+1, XPbPb_NP_err[i]);

		cout << "i = " << i << " , Xpp_PR_err : " << Xpp_PR_err[i] << " , XPbPb_PR_err : " << XPbPb_PR_err[i] << endl;
	}

	hXpp_PR->Sumw2();
	hXpp_NP->Sumw2();
	hXPbPb_PR->Sumw2();
	hXPbPb_NP->Sumw2();


	TH1D *hRaa_pp_PR = (TH1D*)hXpp_PR->Clone("hRaa_pp_PR");
	TH1D *hRaa_PbPb_PR = (TH1D*)hXPbPb_PR->Clone("hRaa_PbPb_PR");
	TH1D *hRaa_pp_NP = (TH1D*)hXpp_NP->Clone("hRaa_pp_NP");
	TH1D *hRaa_PbPb_NP = (TH1D*)hXPbPb_NP->Clone("hRaa_PbPb_NP");

	hRaa_PbPb_PR->Divide(hRaa_pp_PR);
	hRaa_PbPb_NP->Divide(hRaa_pp_NP);

	for(int i=0;i<4;i++)
	{
		RaaPR[i] = hRaa_PbPb_PR->GetBinContent(i+1);
		RaaPR_err[i] = hRaa_PbPb_PR->GetBinError(i+1);
		RaaNP[i] = hRaa_PbPb_NP->GetBinContent(i+1);
		RaaNP_err[i] = hRaa_PbPb_NP->GetBinError(i+1);
		x[i] = (ptBin[i+1]+ptBin[i])/2;
		binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		cout << "i : " << i << ", RAA Value : " << RaaPR[i] << " , RAA Err : " << RaaPR_err[i] << " x : " << x[i] << ", bin Width : " << binWidth[i] << endl;
		cout << "   , RAA NP : " << RaaNP[i] <<  ", RaaNP_err : " << RaaNP_err[i] << endl;
		//gRaaNP->SetPoint(i,(ptBin[i+1]-ptBin[i])/2,); 
		//gRaaNP->SetPointError(i,0,Xpp_NP_err); 
	}


	//////////////////////////////// Start Plotting ////////////////////////////////
	float pos_x = 0.20;
	float pos_y = 0.85;
	float pos_y_diff = 0.051;
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
	gXpp_PR->GetYaxis()->SetTitle("#it{B} #times d#sigma/dp_{T} or #it{B} #times (1/T_{AA}N_{MB})dN/dp_{T} (nb/GeV/c)");
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

	TLegend *legPR = new TLegend(0.68,0.75,0.8,0.85);
    SetLegendStyle(legPR);
	legPR->AddEntry(gXpp_PR,"pp", "p");
	legPR->AddEntry(gXPbPb_PR,"PbPb, Cent. 0-90%", "p");
	legPR->Draw("SAME");
	jumSun(0,1,50,1);

	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXPR,iPeriod,iPos);	

	cXPR->SaveAs("CrossSection_PR_y1p6_2p4_pT.pdf");

	TCanvas *cXNP = new TCanvas("cXNP", "", 700,700);
	cXNP->cd();
	cXNP->SetLogy();
	gXpp_NP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gXpp_NP->GetXaxis()->CenterTitle();
	gXpp_NP->GetYaxis()->SetTitle("#it{B} #times d#sigma/dp_{T} or #it{B} #times (1/T_{AA}N_{MB})dN/dp_{T} (nb/GeV/c)");
	gXpp_NP->GetYaxis()->SetTitleSize(0.04);
	gXpp_NP->GetYaxis()->SetTitleOffset(1.90);
	gXpp_NP->SetTitle("");
	gXpp_NP->GetYaxis()->CenterTitle();
	gXpp_NP->GetXaxis()->SetLimits(0.,50);
	gXpp_NP->SetMinimum(1e-11);
	gXpp_NP->SetMaximum(1e-4);
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
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXNP,iPeriod,iPos);	

	cXNP->SaveAs("CrossSection_NP_y1p6_2p4_pT.pdf");

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
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cRAA,iPeriod,iPos);	

	cRAA->SaveAs("RAA_psi2S_y1p6_2p4_pT.pdf");

	TFile *f1 = new TFile("test.root","recreate");
	f1->cd();
	hyieldPP_PR->Write();
	hyieldPbPb_PR->Write();
	hyieldPP_NP->Write();
	hyieldPbPb_NP->Write();

	hXpp_PR->Write();
	hXPbPb_PR->Write();
	f1->Close();

}

valErr getYield_pp(int i){
    double ptBins[4] = {3,6.5,12,50};
    TString kineLabel[4];
    kineLabel[i] = getKineLabelpp(ptBins[i],ptBins[i+1],1.6,2.4,0.0);
    TFile* inf = new TFile(Form("./pp_psi2S/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getYield_PbPb(int i){
    double ptBins[4] = {3,6.5,12,50};
    TString kineLabel[4];
    kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
    TFile* inf = new TFile(Form("./psi2S/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getFrac_PbPb(int i) {
    double ptBin[4] = {3,6.5,12,50};
    TString kineLabel[4];
    kineLabel[i] = getKineLabel(ptBin[i],ptBin[i+1],1.6,2.4,0.0,0,180);
    TFile* inf = new TFile(Form("./psi2S/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getFrac_pp(int i) {
    double ptBin[4] = {3,6.5,12,50};
    TString kineLabel[4];
    kineLabel[i] = getKineLabelpp(ptBin[i],ptBin[i+1],1.6,2.4,0.0);
    TFile* inf = new TFile(Form("./pp_psi2S/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
double getAccWeight(TH1D* h, double pt){
    //cout<<"pt: "<<pt<<endl;
    double binN = h->FindBin(pt);
    double weight_ = 1./(h->GetBinContent(binN));
    return weight_;
}
  double getEffWeight(TH1D *h, double pt){
    double binN = h->FindBin(pt);
    double weight_ = 1./(h->GetBinContent(binN));
    return weight_;
}
