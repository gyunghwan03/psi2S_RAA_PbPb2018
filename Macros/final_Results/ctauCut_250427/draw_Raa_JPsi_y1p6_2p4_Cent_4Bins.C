#include "../../../rootFitHeaders.h"
#include "../../../commonUtility.h"
#include "../../../JpsiUtility.h"
#include "../../../cutsAndBin.h"
#include "../../../CMS_lumi_v2mass.C"
#include "../../../tdrstyle.C"
#include "../../../Style.h"
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

valErr getYield_PbPb(int i=0, int isPR=0);
valErr getYield_pp(int i=0, int isPR=0);
valErr getCtauEff_pp(int i=0, int isPR=0);
valErr getCtauEff_PbPb(int i=0, int isPR=0);

void draw_Raa_JPsi_y1p6_2p4_Cent_4Bins(bool isSys=true)
{

    //gROOT->SumW2();
    gStyle->SetOptStat(0);
    setTDRStyle();
    int iPeriod = 101;
    int iPos = 33;

    const int nCentBins=4;
    TFile *fPbPb[nCentBins+1];
    TFile *fpp[nCentBins+1];

	TFile *fSys = new TFile("../../syst_summary/syst_roots/total_syst.root");
    //TFile *fSys = new TFile("../syst_summary/syst_roots/total_syst_NobFrac.root");
	TFile *fEff_PbPbPR = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_Jpsi_PtW1_tnp1_20240805.root");
    TFile *fEff_PbPbNP = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_Jpsi_PtW1_tnp1_20240805.root"); //Must Change
    TFile *fEff_ppPR = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtW1_tnp1_20240805.root");
    TFile *fEff_ppNP = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtW1_tnp1_20240805.root"); //Must Change
    TFile *fAcc_ppPR = new TFile("../../../Eff_Acc/roots/acceptance_Prompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20240805.root");
    TFile *fAcc_PbPbPR = new TFile("../../../Eff_Acc/roots/acceptance_Prompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20240805.root");
    TFile *fAcc_ppNP = new TFile("../../../Eff_Acc/roots/acceptance_NPrompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20240805_NP.root"); //Must Change
    TFile *fAcc_PbPbNP = new TFile("../../../Eff_Acc/roots/acceptance_NPrompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20240805_NP.root"); // Must Change

    TH1D *hEff_PbPbPR = (TH1D*) fEff_PbPbPR -> Get("mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4_2");
    TH1D *hEff_PbPbNP = (TH1D*) fEff_PbPbNP -> Get("mc_eff_vs_cent_TnP1_PtW1_pt_3_to_40_absy1p6_2p4_2");
    TH1D *hEff_ppPR = (TH1D*) fEff_ppPR -> Get("mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4");
    TH1D *hEff_ppNP = (TH1D*) fEff_ppNP -> Get("mc_eff_Integrated_TnP1_PtW1_absy1p6_2p4");
    TH1D *hAcc_PbPbPR = (TH1D*) fAcc_PbPbPR -> Get("hAccPt_2021_Fory_Int");
    TH1D *hAcc_ppPR = (TH1D*) fAcc_ppPR -> Get("hAccPt_2021_Fory_Int");
    TH1D *hAcc_PbPbNP = (TH1D*) fAcc_PbPbNP -> Get("hAccPt_2021_Fory_Int");
    TH1D *hAcc_ppNP = (TH1D*) fAcc_ppNP -> Get("hAccPt_2021_Fory_Int");

	TH1D *hSys_PR = (TH1D*) fSys -> Get("fwd_cent_PR");
    TH1D *hSys_NP = (TH1D*) fSys -> Get("fwd_cent_NP");

    Double_t Nmb = 11968044281.;
    //Double_t Taa = 5.649; // 0-100%
    Double_t lumi_pp = 3.002;
    Double_t lumi_pp_scale = 1e-9;
    Double_t lumi_pp_err = 0.019;
    Double_t Nmb_err = 0.01261;
    //Double_t Taa_err = 0.123; // 0-100%

    //double centBin[nCentBins+1] = {0,20,30,40,50,90};
    //double centBin[nCentBins+1] = {0,10,20,30,40,50,90};
	double centBin[nCentBins+1] = {0,10,30,50,90};
    //double NpartBin[nCentBins+1] = {27.12,87.19,131.0,188.2,309.6};
    //double NpartBin[nCentBins+1] = {27.12,87.19,131.0,188.2,262.3,356.9};
	double NpartBin[nCentBins] = {27.12,109.1,225.2,356.9};
    double fracPP[nCentBins]; double fracPbPb[nCentBins];
    double fracErrPP[nCentBins]; double fracErrPbPb[nCentBins];

    TString kineLabel[nCentBins+1];

    TH1D *hyieldPP_PR = new TH1D("hyieldPP_PR",";Centrality (%);",nCentBins, centBin);
    TH1D *hyieldPbPb_PR = new TH1D("hyieldPbPb_PR",";Centrality (%);",nCentBins, centBin);
    TH1D *hyieldPP_NP = new TH1D("hyieldPP_NP",";Centrality (%);",nCentBins, centBin);
    TH1D *hyieldPbPb_NP = new TH1D("hyieldPbPb_NP",";Centrality (%);",nCentBins, centBin);

    TH1D *hXpp_PR = new TH1D("hXpp_PR",";Centrality (%);",nCentBins, centBin);
    TH1D *hXpp_NP = new TH1D("hXpp_NP",";Centrality (%);",nCentBins, centBin);
    TH1D *hXPbPb_PR = new TH1D("hXPbPb_PR",";Centrality (%);",nCentBins, centBin);
    TH1D *hXPbPb_NP = new TH1D("hXPbPb_NP",";Centrality (%);",nCentBins, centBin);

    TH1D *hRAA_PR = new TH1D("hRAA_PR",";p_{T} (GeV/c);", nCentBins,NpartBin);
    TH1D *hRAA_NP = new TH1D("hRAA_NP",";p_{T} (GeV/c);", nCentBins,NpartBin);

    double RaaPR[nCentBins]; double RaaPR_err[nCentBins]; double binWidth[nCentBins]={4.3,4.3,4.3,4.3}; double x[nCentBins];
    double RaaNP[nCentBins]; double RaaNP_err[nCentBins];
    double SysPR[nCentBins]; double SysNP[nCentBins];

    double cfrac[nCentBins];
	double Taa[nCentBins] = {23.05, 14.39, 3.950, 0.5803};
    double Taa_err[nCentBins] = {0.42, 0.26, 0.132, 0.0288};
    //double Taa[nCentBins] = {23.05, 14.39, 8.798, 5.124, 2.777, 0.5803};
    //double Taa_err[nCentBins] = {0.42, 0.30, 0.219, 0.159, 0.107, 0.0288};


    double Xpp_PR[nCentBins]; double Xpp_NP[nCentBins];
    double Xpp_PR_err[nCentBins]; double Xpp_NP_err[nCentBins];
    double XPbPb_PR[nCentBins]; double XPbPb_NP[nCentBins];
    double XPbPb_PR_err[nCentBins]; double XPbPb_NP_err[nCentBins];

    for(int i=0;i<nCentBins;i++)
    {
        double weight_ppPR = 1; double weight_ppNP = 1; double weight_PbPbPR = 1; double weight_PbPbNP = 1;
        double eff_ppPR = 1; double eff_ppNP = 1; double eff_PbPbPR = 1; double eff_PbPbNP = 1;
        double acc_ppPR = 1; double acc_PbPbPR = 1; double acc_ppNP = 1; double acc_PbPbNP = 1;

        eff_ppPR=hEff_ppPR->GetBinContent(1);
        eff_ppNP=hEff_ppNP->GetBinContent(1);
        acc_ppPR=hAcc_ppPR->GetBinContent(1);
        acc_ppNP=hAcc_ppNP->GetBinContent(1);
        eff_PbPbPR=hEff_PbPbPR->GetBinContent(i+1);
        eff_PbPbNP=hEff_PbPbNP->GetBinContent(i+1);
        acc_PbPbPR=hAcc_PbPbPR->GetBinContent(1);
        acc_PbPbNP=hAcc_PbPbNP->GetBinContent(1);

        weight_ppPR=eff_ppPR*acc_ppPR;
        weight_ppNP=eff_ppNP*acc_ppNP;
        weight_PbPbPR=eff_PbPbPR*acc_PbPbPR;
        weight_PbPbNP=eff_PbPbNP*acc_PbPbNP;

		cout.precision(9);
		cout << " " << endl;
		cout << "Cent " << centBin[i] << " - " << centBin[i+1] << "	Cent Frac : " << cfrac[i] << endl;
		cout << "weight_ppPR :\t" << weight_ppPR << "\tweight_ppNP :\t" << weight_ppNP << "\tweight_PbPbPR :\t" << weight_PbPbPR << "\tweight_PbPbNP :\t" << weight_PbPbNP << endl;
		cout << " " << endl;

		valErr yieldPP; valErr yieldPbPb; valErr ctauEffPP; valErr ctauEffPbPb;

		double yieldPP_PR; double yieldPP_NP; double yieldPbPb_PR; double yieldPbPb_NP;
		double ctauEffPP_PR; double ctauEffPP_NP; double ctauEffPbPb_PR; double ctauEffPbPb_NP;
		double err_PPPR; double err_PPNP; double err_PbPbPR; double err_PbPbNP;

		for (int isPR = 0; isPR < 2; isPR++)
		{
			yieldPP = getYield_pp(i, isPR);
			yieldPbPb = getYield_PbPb(i, isPR);
			ctauEffPP = getCtauEff_pp(i, isPR);
			ctauEffPbPb = getCtauEff_PbPb(i, isPR);
			
			if(isPR == 0)
			{
				ctauEffPP_PR = ctauEffPP.val;
				ctauEffPbPb_PR = ctauEffPbPb.val;
				weight_ppPR   = weight_ppPR * ctauEffPP_PR; 
				weight_PbPbPR = weight_PbPbPR * ctauEffPbPb_PR; 
				yieldPP_PR = yieldPP.val/(weight_ppPR);
				yieldPbPb_PR = yieldPbPb.val/(weight_PbPbPR);
				err_PPPR = yieldPP.err;
				err_PbPbPR = yieldPbPb.err;
			}
			else if(isPR == 1)
			{
				ctauEffPP_NP = ctauEffPP.val;
				ctauEffPbPb_NP = ctauEffPbPb.val;
				weight_ppNP   = weight_ppNP * ctauEffPP_NP; 
				weight_PbPbNP = weight_PbPbNP * ctauEffPbPb_NP; 
				yieldPP_NP = yieldPP.val/(weight_ppNP);
				yieldPbPb_NP = yieldPbPb.val/(weight_PbPbNP);
				//if(i==0) yieldPbPb_NP = 0/(weight_PbPbNP);
				err_PPNP = yieldPP.err;
				err_PbPbNP = yieldPbPb.err;
			}
		}
 
		double err_PbPbPR_wgt = err_PbPbPR/weight_PbPbPR;
        double err_PbPbNP_wgt = err_PbPbNP/weight_PbPbNP;
        double err_PPPR_wgt = err_PPPR/weight_ppPR;
        double err_PPNP_wgt = err_PPNP/weight_ppNP;

        cout << " " << endl;
        cout << "Cent Bin : " << centBin[i] << "-" << centBin[i+1] << endl;
        cout << "Yield pp PR :\t"   << yieldPP_PR*weight_ppPR   << " +/- " << err_PPPR   << "\tEff :\t" << eff_ppPR   << "\tAcc :\t" << acc_ppPR   << "\tctauEff :\t" << ctauEffPP_PR   << "\tweight :\t" << weight_ppPR << endl;
		cout << "Yield PbPb PR :\t" << yieldPbPb_PR*weight_PbPbPR << " +/- " << err_PbPbPR << "\tEff :\t" << eff_PbPbPR << "\tAcc :\t" << acc_PbPbPR << "\tctauEff :\t" << ctauEffPbPb_PR << "\tweight :\t" << weight_PbPbPR << endl;
		cout << "Yield pp NP :\t"   << yieldPP_NP*weight_ppNP   << " +/- " << err_PPNP   << "\tEff :\t" << eff_ppNP   << "\tAcc :\t" << acc_ppNP   << "\tctauEff :\t" << ctauEffPP_NP   << "\tweight :\t" << weight_ppNP << endl;
		cout << "Yield PbPb NP :\t" << yieldPbPb_NP*weight_PbPbNP << " +/- " << err_PbPbNP << "\tEff :\t" << eff_PbPbNP << "\tAcc :\t" << acc_PbPbNP << "\tctauEff :\t" << ctauEffPbPb_NP << "\tweight :\t" << weight_PbPbNP << endl;
        //cout << "pp Prompt yield : " << yieldPP_PR << ", PbPb Prompt yield : " << yieldPbPb_PR << ", pp NonPrompt yield : " << yieldPP_NP << ", PbPb NonPrompt yield : " << yieldPbPb_NP << endl;
        //cout << "pp Prompt error : " << err_PPPR_wgt << " , PbPb Prompt error : " << err_PbPbPR_wgt << " , pp NonPrompt error : " << err_PPNP_wgt << " , PbPb NonPrompt error : " << err_PbPbNP_wgt << endl;

        /*
        hyieldPP_PR -> SetBinContent(i+1,yieldPP.val*(1-fracPP.val));
        hyieldPP_PR -> SetBinError(i+1,err_PPPR_wgt);
        hyieldPbPb_PR -> SetBinContent(i+1,yieldPbPb.val*(1-fracPbPb.val));
        hyieldPbPb_PR -> SetBinError(i+1,err_PbPbPR_wgt);

        hyieldPP_NP -> SetBinContent(i+1,yieldPP.val*(fracPP.val));
        hyieldPP_NP -> SetBinError(i+1,err_PPNP_wgt);
        hyieldPbPb_NP -> SetBinContent(i+1,yieldPbPb.val*(fracPbPb.val));
        hyieldPbPb_NP -> SetBinError(i+1,err_PbPbNP_wgt);
        */

        cfrac[i] = (centBin[i+1]-centBin[i])/90.;

        Xpp_PR[i] = lumi_pp_scale*yieldPP_PR/(lumi_pp*1e+2*(double)(40.-3.5)*2*(double)(2.4-1.6));
		//Xpp_PR_err[i] = lumi_pp_scale*Xpp_PR[i]*sqrt(TMath::Power(err_PPPR_wgt/yieldPP_PR,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_PR_err[i] = Xpp_PR[i]*sqrt(TMath::Power(err_PPPR_wgt/yieldPP_PR,2));// + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_NP[i] = lumi_pp_scale*yieldPP_NP/(lumi_pp*1e+2*(double)(40.-3.5)*2*(double)(2.4-1.6));
		//Xpp_NP_err[i] = lumi_pp_scale*Xpp_NP[i]*sqrt(TMath::Power(err_PPNP_wgt/yieldPP_NP,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_NP_err[i] = Xpp_NP[i]*sqrt(TMath::Power(err_PPNP_wgt/yieldPP_NP,2));// + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));

		XPbPb_PR[i] = yieldPbPb_PR/(Nmb*Taa[i]*(double)(40.-3.5)*(double)2*(2.4-1.6)*cfrac[i]);
		//XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(Taa_err[i]/Taa[i],2) + TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2));// + TMath::Power(Nmb_err/Nmb,2));
		XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2));// + TMath::Power(Nmb_err/Nmb,2));
		XPbPb_NP[i] = yieldPbPb_NP/(Nmb*Taa[i]*(double)(40.-3.5)*(double)2*(2.4-1.6)*cfrac[i]);
		//XPbPb_NP_err[i] = XPbPb_NP[i]*sqrt(TMath::Power(Taa_err[i]/Taa[i],2) + TMath::Power(err_PbPbNP_wgt/yieldPbPb_NP,2));// + TMath::Power(Nmb_err/Nmb,2));
		XPbPb_NP_err[i] = XPbPb_NP[i]*sqrt(TMath::Power(err_PbPbNP_wgt/yieldPbPb_NP,2));// + TMath::Power(Nmb_err/Nmb,2));

        hXpp_PR->SetBinContent(i+1, Xpp_PR[i]);
        hXpp_NP->SetBinContent(i+1, Xpp_NP[i]);
        hXPbPb_PR->SetBinContent(i+1, XPbPb_PR[i]);
        hXPbPb_NP->SetBinContent(i+1, XPbPb_NP[i]);

        hXpp_PR->SetBinError(i+1, Xpp_PR_err[i]);
        hXpp_NP->SetBinError(i+1, Xpp_NP_err[i]);
        hXPbPb_PR->SetBinError(i+1, XPbPb_PR_err[i]);
        hXPbPb_NP->SetBinError(i+1, XPbPb_NP_err[i]);

        //cout << "Cent Frac : " << cfrac[i] << ", Taa : " << Taa[i]  << endl;
        //cout << "Xpp_PR : " << Xpp_PR[i] << " , Xpp_PR_err : " << Xpp_PR_err[i] << ", XPbPb_PR : " << XPbPb_PR[i] << " , XPbPb_PR_err : " << XPbPb_PR_err[i] << endl;

        /*Taa_Nmb[i] =  Taa[i]*Nmb*2*(2.4-1.6)*(50.-3.);
        step_one[i] = (Taa_Nmb[i]*Xpp_PR[i]);
        step_two[i] = step_one[i]*cfrac[i];
        RAA[i] = step_two[i]*yieldPbPb_PR;*/
        //cout << "Taa_Nmb[i] = " << Taa_Nmb[i] << " ,step one = " << step_one[i] << " , step two = " << step_two[i] << ", RAA = " << RAA[i] << endl;
    }

	cout << " " << endl;
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

/*
	for(int i=nCentBins-1;i>=0;i--)
	{
		cout << "nCentBins-1-i = " << nCentBins-1-i << " i = : " << i << endl;
		RaaPR[nCentBins-1-i] = (hRaa_PbPb_PR->GetBinContent(i+1))*cfrac[nCentBins-1-i];
		RaaPR_err[nCentBins-1-i] = (hRaa_PbPb_PR->GetBinError(i+1))*cfrac[nCentBins-1-i];
		RaaNP[nCentBins-1-i] = (hRaa_PbPb_NP->GetBinContent(i+1))*cfrac[nCentBins-1-i];
		RaaNP_err[nCentBins-1-i] = (hRaa_PbPb_NP->GetBinError(i+1))*cfrac[nCentBins-1-i];
	}
*/
	for(int i=nCentBins-1;i>=0;i--)
	{
		//cout << "nCentBins-1-i = " << nCentBins-1-i << " i = : " << i << endl;
		RaaPR[nCentBins-1-i] = (hRaa_PbPb_PR->GetBinContent(i+1));
		RaaPR_err[nCentBins-1-i] = (hRaa_PbPb_PR->GetBinError(i+1));
		RaaNP[nCentBins-1-i] = (hRaa_PbPb_NP->GetBinContent(i+1));
		RaaNP_err[nCentBins-1-i] = (hRaa_PbPb_NP->GetBinError(i+1));
		SysPR[nCentBins-1-i] = RaaPR[nCentBins-1-i]*(hSys_PR->GetBinContent(i+1));
        SysNP[nCentBins-1-i] = RaaNP[nCentBins-1-i]*(hSys_NP->GetBinContent(i+1));

	    cout << " " << endl;
		cout << "Cent. " << centBin[i] << " - " << centBin[i+1] << endl;
        cout << "Raa PR : " << RaaPR[nCentBins-1-i] << " , PR error : " << RaaPR_err[nCentBins-1-i] << endl;
        cout << "Raa NP : " << RaaNP[nCentBins-1-i] << " , NP error : " << RaaNP_err[nCentBins-1-i] << endl;
	}
	cout << " " << endl;

	for (int i=0; i<nCentBins; i++){
        cout << "RAA PR : " << RaaPR[i] << ", RAA NP: " << RaaNP[i] << endl;
        hRAA_PR->SetBinContent(i+1,RaaPR[i]);
        hRAA_NP->SetBinContent(i+1,RaaNP[i]);
        hRAA_PR->SetBinError(i+1,RaaPR_err[i]);
        hRAA_NP->SetBinError(i+1,RaaNP_err[i]);
    }



	//////////////////////////////// Tex Table ////////////////////////////////
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;
	cout.precision(3);
	cout << "\\multicolumn{2}{c}{$1.6 < |y| < 2.4$, $3.5 < \\pt < 40$ (GeV/c)} \\\\ \\hline" << endl;
	for(int i=0;i<nCentBins;i++)
	{
		cout << centBin[i] << " -- " << centBin[i+1] << " & " << RaaPR[i] << " $\\pm$ " << RaaPR_err[i] << " (stat.) $\\pm$ " << SysPR[i] << " (syst.) \\\\" << endl;
	}
	cout << " " << endl;
	cout << "\\multicolumn{2}{c}{$1.6 < |y| < 2.4$, $3.5 < \\pt < 40$ (GeV/c} \\\\ \\hline" << endl;
	for(int i=0;i<nCentBins;i++)
	{
		cout << centBin[i] << " -- " << centBin[i+1] << " & " << RaaNP[i] << " $\\pm$ " << RaaNP_err[i] << " (stat.) $\\pm$ " << SysNP[i] << " (syst.) \\\\" << endl;
	}

	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;



	//////////////////////////////// Start Plotting ////////////////////////////////
	float pos_x = 0.20;
	float pos_y = 0.85;
	float pos_y_diff = 0.051;
	int text_color = 1;
	float text_size = 25;

	TGraphErrors *gXpp_PR = new TGraphErrors(nCentBins,x,Xpp_PR,binWidth,Xpp_PR_err);
	TGraphErrors *gXpp_NP = new TGraphErrors(nCentBins,x,Xpp_NP,binWidth,Xpp_NP_err);

	TGraphErrors *gXPbPb_PR = new TGraphErrors(nCentBins,x,XPbPb_PR,binWidth,XPbPb_PR_err);
	TGraphErrors *gXPbPb_NP = new TGraphErrors(nCentBins,x,XPbPb_NP,binWidth,XPbPb_NP_err);

	TCanvas *cXPR = new TCanvas("cXPR", "", 700,700);
	cXPR->cd();
	cXPR->SetLogy();
	gXpp_PR->GetXaxis()->SetTitle("Centrality (%)");
	//gXpp_PR->GetXaxis()->SetTitle("<N_{Part}>");
	gXpp_PR->GetXaxis()->CenterTitle();
	gXpp_PR->GetYaxis()->SetTitle("#it{B} #times d#sigma/dp_{T} or #it{B} #times (1/T_{AA}N_{MB})dN/dp_{T} (nb/GeV/c)");
	gXpp_PR->GetYaxis()->SetTitleSize(0.04);
	gXpp_PR->GetYaxis()->SetTitleOffset(1.90);
	gXpp_PR->SetTitle("");
	gXpp_PR->GetYaxis()->CenterTitle();
	gXpp_PR->GetXaxis()->SetLimits(0.,100);
	gXpp_PR->SetMinimum(1e-10);
	gXpp_PR->SetMaximum(1e-4);
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
	legPR->AddEntry(gXpp_PR,"pp", "pe");
	legPR->AddEntry(gXPbPb_PR,"PbPb", "pe");
	legPR->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXPR,iPeriod,iPos);	

	cXPR->SaveAs("figs/CrossSection_JPsi_PR_y1p6_2p4_Cent.pdf");

	TCanvas *cXNP = new TCanvas("cXNP", "", 700,700);
	cXNP->cd();
	cXNP->SetLogy();
	gXpp_NP->GetXaxis()->SetTitle("Centrality (%)");
	//gXpp_NP->GetXaxis()->SetTitle("<N_{Part}>");
	gXpp_NP->GetXaxis()->CenterTitle();
	gXpp_NP->GetYaxis()->SetTitle("#it{B} #times d#sigma/dp_{T} or #it{B} #times (1/T_{AA}N_{MB})dN/dp_{T} (nb/GeV/c)");
	gXpp_NP->GetYaxis()->SetTitleSize(0.04);
	gXpp_NP->GetYaxis()->SetTitleOffset(1.90);
	gXpp_NP->SetTitle("");
	gXpp_NP->GetYaxis()->CenterTitle();
	gXpp_NP->GetXaxis()->SetLimits(0.,100);
	gXpp_NP->SetMinimum(1e-10);
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
	legNP->AddEntry(gXpp_NP,"pp", "pe");
	legNP->AddEntry(gXPbPb_NP,"PbPb", "pe");
	legNP->Draw("SAME");
	jumSun(0,1,400,1);

	drawText("Non Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXNP,iPeriod,iPos);	

	cXNP->SaveAs("figs/CrossSection_JPsi_NP_y1p6_2p4_Cent.pdf");

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  RAA  ///////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	TGraphErrors *gRaaPR = new TGraphErrors(nCentBins,NpartBin,RaaPR,0,RaaPR_err);
	TGraphErrors *gRaaNP = new TGraphErrors(nCentBins,NpartBin,RaaNP,0,RaaNP_err);
	TGraphErrors *gSysPR = new TGraphErrors(nCentBins,NpartBin,RaaPR,binWidth,SysPR);
    TGraphErrors *gSysNP = new TGraphErrors(nCentBins,NpartBin,RaaNP,binWidth,SysNP);

	TCanvas *cRAA = new TCanvas("cRAA", "", 800, 700);
	TPad *pad1_1, *pad1_2;
	pad1_1 = new TPad("pad1_1", "", 0.00,0,0.83,1);
    pad1_1 -> SetRightMargin(0);
    pad1_2 = new TPad("pad1_2", "", 0.83,0,1,1);
    pad1_2 -> SetLeftMargin(0);
    pad1_1->SetTicks();
    pad1_2->SetTicks();
	cRAA->cd();
	pad1_2->Draw();
	pad1_1->Draw();

	pad1_1->cd();
	gRaaPR->GetXaxis()->SetTitle("<N_{Part}>");
	gRaaPR->GetXaxis()->CenterTitle();
	gRaaPR->GetYaxis()->SetTitle("R_{AA}");
	gRaaPR->SetTitle("");
	gRaaPR->GetYaxis()->CenterTitle();
	gRaaPR->GetXaxis()->SetLimits(0.,400);
	gRaaPR->SetMinimum(0);
	gRaaPR->SetMaximum(1.44);

	gRaaPR->SetMarkerColor(kBlue+2);
	gRaaPR->SetLineColor(kBlue+2);
	gRaaPR->SetMarkerStyle(20);
	gRaaPR->SetMarkerSize(1.5);
	gRaaNP->SetMarkerColor(kRed+3);
	gRaaNP->SetLineColor(kRed+3);
	gRaaNP->SetMarkerStyle(21);
	gRaaNP->SetMarkerSize(1.5);

	gSysPR->SetMinimum(0);
	gSysPR->SetMaximum(1.44);
	gSysPR->SetLineColor(kBlue-4);
	gSysPR->SetFillColorAlpha(kBlue-9,0.40);
	gSysNP->SetLineColor(kRed-4);
	gSysNP->SetFillColorAlpha(kRed-9,0.40);

	gRaaPR->Draw("AP");
	gRaaNP->Draw("P");
	if(isSys==1){
		gSysPR->Draw("5");
		gSysNP->Draw("5");
	}

	TLegend *leg = new TLegend(0.68,0.72,0.8,0.82);
	SetLegendStyle(leg);
	leg->AddEntry(gRaaPR,"Prompt J/#psi", "pe");
	leg->AddEntry(gRaaNP,"Non Prompt J/#psi", "pe");
	leg->Draw("SAME");
	jumSun(0,1,400,1);

	double glb_Err=0;
	glb_Err = TMath::Sqrt( TMath::Power(Nmb_err,2) + TMath:: Power(lumi_pp_err,2) );

	cout << "Global Uncertainty : " << glb_Err << endl;

	TBox* b_ynerr = new TBox();
	b_ynerr->SetX1(0);
	b_ynerr->SetX2(15);
	b_ynerr->SetY1(1 - glb_Err);
	b_ynerr->SetY2(1 + glb_Err);
	b_ynerr->SetFillColorAlpha(12, 0.8);
	b_ynerr->SetLineWidth(1);
	b_ynerr->SetLineColor(kBlack);
	b_ynerr->Draw("L");

	drawText("3.5 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(pad1_1,iPeriod,iPos);	

	pad1_2->cd();
	double midPR2S_int  = 0.427; double midPR2S_intErr = 0.0021;
	double midNP2S_int  = 0.424; double midNP2S_intErr = 0.0029;
	TH1D *hmidPR_int = new TH1D("hmidPR_int",";;",3,0,1);
	TH1D *hmidNP_int = new TH1D("hmidNP_int",";;",3,0,1);

	hmidPR_int -> SetBinContent(2,midPR2S_int);
	hmidPR_int -> SetBinError(2,midPR2S_intErr);
	hmidNP_int -> SetBinContent(2,midNP2S_int);
	hmidNP_int -> SetBinError(2,midNP2S_intErr);

	hmidPR_int->SetMinimum(0);
	hmidPR_int->SetMaximum(1.44);

	hmidPR_int->SetMarkerColor(kBlue+2);
	hmidPR_int->SetLineColor(kBlue+2);
	hmidPR_int->SetMarkerStyle(20);
	hmidPR_int->SetMarkerSize(1.5);
	hmidNP_int->SetMarkerColor(kRed+3);
	hmidNP_int->SetLineColor(kRed+3);
	hmidNP_int->SetMarkerStyle(21);
	hmidNP_int->SetMarkerSize(1.5);

	hmidPR_int->GetYaxis()->SetLabelOffset(0);
    hmidPR_int->GetXaxis()->SetLabelOffset(0);
    hmidPR_int->GetXaxis()->SetLabelSize(0);
    hmidPR_int->GetXaxis()->SetTickLength(0);
    hmidPR_int->GetYaxis()->SetLabelSize(0.12);
	hmidNP_int->GetYaxis()->SetLabelOffset(0);
    hmidNP_int->GetXaxis()->SetLabelOffset(0);
    hmidNP_int->GetXaxis()->SetLabelSize(0);
    hmidNP_int->GetXaxis()->SetTickLength(0);
    hmidNP_int->GetYaxis()->SetLabelSize(0.12);
	
	hmidPR_int->Draw("PE");
	hmidNP_int->Draw("PE SAME");
	jumSun(0,1,1,1);

	cRAA->SaveAs(Form("figs/RAA_JPsi_y1p6_2p4_Cent_Sys%d_4Bins.pdf",isSys));

	TFile *f1 = new TFile("roots/RAA_JPsi_forRap_Npart_4Bins.root","recreate");
	f1->cd();
	hRAA_PR->Write();
	hRAA_NP->Write();
	f1->Close();

}

valErr getYield_pp(int i, int isPR){
    TString kineLabel;
    kineLabel = getKineLabelpp(3.5,40,1.6,2.4,0.0);
    TString PR = (isPR==0) ? "PRMC" : "NPMC";
	TFile* inf = new TFile(Form("../../Jpsi_L_cut/roots_1S_pp/%s/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", PR.Data(), kineLabel.Data()));
    valErr ret; ret.val = 0; ret.err = 0;
    if (!inf || inf->IsZombie()) {
        cout << "파일을 열 수 없습니다: " << inf->GetName() << endl;
        return ret;
    }
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");
    if (!fitResults) {
        cout << "fitResults 히스토그램이 없습니다! 파일: " << inf->GetName() << endl;
        inf->Close();
        return ret;
    }
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    inf->Close();
    return ret;
}
valErr getYield_PbPb(int i, int isPR){
    double centBin[5] = {0,20,60,100,180};//{0,20,40,60,80,100,180};
    TString kineLabel[5];
	TString PR = (isPR==0) ? "PRMC" : "NPMC";	
    kineLabel[i] = getKineLabel(3.5,40,1.6,2.4,0.0,centBin[i],centBin[i+1]);
    TFile *inf = new TFile(Form("../../Jpsi_L_cut/roots_1S_Pb/%s/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", PR.Data(), kineLabel[i].Data()));
	valErr ret; ret.val = 0; ret.err = 0;
	if (!inf || inf->IsZombie()) {
		cout << "파일을 열 수 없습니다: " << inf->GetName() << endl;
		return ret;
	}
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");
	if (!fitResults) {
		cout << "fitResults 히스토그램이 없습니다! 파일: " << inf->GetName() << endl;
		inf->Close();
		return ret;
	}
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
	inf->Close();
    return ret;
}
valErr getCtauEff_pp(int i, int isPR){
	TString kineLabel;
	kineLabel = getKineLabelpp(3.5,40,1.6,2.4,0.0);
	TString PR;
	if(isPR==0) PR = "PRMC";
	else if(isPR==1) PR = "NPMC";
	TFile* inf = new TFile(Form("../../../ctau_eff/roots_1S_pp/decayL/%s/decay_hist_%s_m2.6-3.5.root", PR.Data(), kineLabel.Data()));
	valErr ret; ret.val = 0; ret.err = 0;
	if (!inf || inf->IsZombie()) {
		cout << "파일을 열 수 없습니다: " << inf->GetName() << endl;
		return ret;
	}
	TH1D* h_ctau = (TH1D*)inf->Get("h_eff");
	if (!h_ctau) {
		cout << "hEff 히스토그램이 없습니다! 파일: " << inf->GetName() << endl;
		inf->Close();
		return ret;
	}
	ret.val = h_ctau->GetBinContent(1);
	ret.err = h_ctau->GetBinError(1);
	inf->Close();
	return ret;
}
valErr getCtauEff_PbPb(int i, int isPR){
	double centBin[5] = {0,20,60,100,180};//{0,20,40,60,80,100,180};
	TString kineLabel[5];
	kineLabel[i] = getKineLabel(3.5,40,1.6,2.4,0.0,centBin[i],centBin[i+1]);
	TString PR;
	if(isPR==0) PR = "PRMC";
	else if(isPR==1) PR = "NPMC";
	TFile* inf = new TFile(Form("../../../ctau_eff/roots_1S_Pb/decayL/%s/decay_hist_%s_m2.6-3.5.root", PR.Data(), kineLabel[i].Data()));
	valErr ret; ret.val = 0; ret.err = 0;
	if (!inf || inf->IsZombie()) {
		cout << "파일을 열 수 없습니다: " << inf->GetName() << endl;
		return ret;
	}
	TH1D* hEff = (TH1D*)inf->Get("h_eff");
	if (!hEff) {
		cout << "hEff 히스토그램이 없습니다! 파일: " << inf->GetName() << endl;
		inf->Close();
		return ret;
	}
	ret.val = hEff->GetBinContent(1);
	ret.err = hEff->GetBinError(1);
	inf->Close();
	return ret;
}