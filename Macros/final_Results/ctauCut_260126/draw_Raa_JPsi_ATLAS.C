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
//valErr getCtauEff_pp(int i=0, int isPR=0);
//valErr getCtauEff_PbPb(int i=0, int isPR=0);
//valErr getYield_PbPb_Integrated(int isPR=0);
//valErr getCtauEff_PbPb_Integrated(int isPR=0);


void draw_Raa_JPsi_ATLAS(bool isSys=false)
{

	//gROOT->SumW2();
	gStyle->SetOptStat(0);
	setTDRStyle();
	int iPeriod = 101;
  	int iPos = 33;

    const int nCentBins=6;
    TFile *fPbPb[nCentBins+1];
    TFile *fpp[nCentBins+1];

	TFile *fSys = new TFile("../../syst_summary_Jpsi/syst_roots/total_syst.root");
	TFile *fEff_PbPbPR = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_ATLAS.root");
	TFile *fEff_PbPbNP = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_JPsi_PtWnomi_tnp1_ctauCut_ATLAS.root");
	//TFile *fEff_PbPbPR = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_Jpsi_PtW1_tnp1_20241126.root");
	//TFile *fEff_PbPbNP = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_Jpsi_PtW1_tnp1_20241126.root"); // Must Change
	//TFile *fEff_PbPbPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_Jpsi_PtW1_tnp1_20240805.root");
    //TFile *fEff_PbPbNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_Jpsi_PtW1_tnp1_20240805.root"); //Must Change
    //TFile *fEff_ppPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtW1_tnp1_20240805.root");
    TFile *fEff_ppPR = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_ATLAS.root");
    TFile *fEff_ppNP = new TFile("../../../Eff_Acc/roots/mc_eff_vs_pt_rap_nprompt_pp_Jpsi_PtWnomi_tnp1_ctauCut_ATLAS.root");
    TFile *fAcc_ppPR = new TFile("../../../Eff_Acc/roots/acceptance_Prompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20241010.root");
    TFile *fAcc_PbPbPR = new TFile("../../../Eff_Acc/roots/acceptance_Prompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20241010.root");
    TFile *fAcc_ppNP   = new TFile("../../../Eff_Acc/roots/acceptance_NonPrompt_JPsi_GenOnly_wgt1_pp_SysUp0_20241014.root"); //Must Change
    TFile *fAcc_PbPbNP = new TFile("../../../Eff_Acc/roots/acceptance_NonPrompt_JPsi_GenOnly_wgt1_pp_SysUp0_20241014.root"); // Must Change
    //TFile *fAcc_ppNP = new TFile("../../Eff_Acc/roots/acceptance_NPrompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20240805_NP.root"); //Must Change
    //TFile *fAcc_PbPbNP = new TFile("../../Eff_Acc/roots/acceptance_NPrompt_Jpsi_GenOnly_wgt1_pp_SysUp0_20240805_NP.root"); // Must Change


    //TH1D *hEff_PbPbPR = (TH1D*) fEff_PbPbPR -> Get("mc_eff_vs_cent_TnP1_PtW0_pt_6p5_to_50_absy0_1p6");
    //TH1D *hEff_PbPbNP = (TH1D*) fEff_PbPbNP -> Get("mc_eff_vs_cent_TnP1_PtW0_pt_6p5_to_50_absy0_1p6");
    //TH1D *hEff_ppPR = (TH1D*) fEff_ppPR -> Get("mc_eff_Integrated_TnP1_PtW0_absy0_1p6");
    //TH1D *hEff_ppNP = (TH1D*) fEff_ppNP -> Get("mc_eff_Integrated_TnP1_PtW0_absy0_1p6");
    TH1D *hEff_PbPbPR = (TH1D*) fEff_PbPbPR -> Get("mc_eff_vs_cent_TnP1_PtW1_pt_9_to_40_absy0_2");
    TH1D *hEff_PbPbNP = (TH1D*) fEff_PbPbNP -> Get("mc_eff_vs_cent_TnP1_PtW1_pt_9_to_40_absy0_2");
    TH1D *hEff_ppPR = (TH1D*) fEff_ppPR -> Get("mc_eff_Integrated_TnP1_PtW1_absy0_2");
    TH1D *hEff_ppNP = (TH1D*) fEff_ppNP -> Get("mc_eff_Integrated_TnP1_PtW1_absy0_2");
    TH1D *hAcc_PbPbPR = (TH1D*) fAcc_PbPbPR -> Get("hAccPt_2021_midy_Int");
    TH1D *hAcc_ppPR = (TH1D*) fAcc_ppPR -> Get("hAccPt_2021_midy_Int");
    TH1D *hAcc_PbPbNP = (TH1D*) fAcc_PbPbNP -> Get("hAccPt_2021_midy_Int");
    TH1D *hAcc_ppNP = (TH1D*) fAcc_ppNP -> Get("hAccPt_2021_midy_Int");

	TH1D *hSys_PR = (TH1D*) fSys -> Get("mid_cent_PR");
    TH1D *hSys_NP = (TH1D*) fSys -> Get("mid_cent_NP");

	Double_t Nmb = 11968044281.*0.8/0.9; //0-80%
	//Double_t Taa = 5.649; // 0-100%
	Double_t lumi_pp = 3.002;
	Double_t lumi_pp_scale = 1e-9;
	Double_t lumi_pp_err = 0.019;
	Double_t Nmb_err = 0.01261;
	//Double_t Taa_err = 0.123; // 0-100%

    double centBin[nCentBins+1] = {0,5,10,20,40,60,80};
	double NpartBin[nCentBins+1] = {7.026,23.52,70.80,159.6,262.3,331.5,382.3};
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

	TH1D *hRAA_PR = new TH1D("hRAA_PR",";Centrality (%);", nCentBins,NpartBin);
	TH1D *hRAA_NP = new TH1D("hRAA_NP",";Centrality (%);", nCentBins,NpartBin);

	double RaaPR[nCentBins]; double RaaPR_err[nCentBins]; double binWidth[nCentBins]={4.3,4.3,4.3,4.3,4.3}; double x[nCentBins];
	double RaaNP[nCentBins]; double RaaNP_err[nCentBins];
	double SysPR[nCentBins]; double SysNP[nCentBins];

	double cfrac[nCentBins];
	double Taa[nCentBins] = {25.7, 20.4, 14.39, 6.961, 2.076, 0.4298};
	double Taa_err[nCentBins] = {0.47, 0.40, 0.3, 0.188, 0.086, 0.0226};


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
		cfrac[i] = (centBin[i+1]-centBin[i])/80.;
		cout.precision(9);
		cout << " " << endl;
		cout << "Cent " << centBin[i] << " - " << centBin[i+1] << "	Cent Frac : " << cfrac[i] << endl;
		cout << "weight_ppPR :\t" << weight_ppPR << "\tweight_ppNP :\t" << weight_ppNP << "\tweight_PbPbPR :\t" << weight_PbPbPR << "\tweight_PbPbNP :\t" << weight_PbPbNP << endl;
		cout << " " << endl;

        valErr yieldPP; valErr yieldPbPb;

		double yieldPP_PR; double yieldPP_NP; double yieldPbPb_PR; double yieldPbPb_NP;
		double err_PPPR; double err_PPNP; double err_PbPbPR; double err_PbPbNP;

		for (int isPR = 0; isPR < 2; isPR++)
		{
			yieldPP = getYield_pp(i, isPR);
			yieldPbPb = getYield_PbPb(i, isPR);
			
			if(isPR == 0)
			{
				yieldPP_PR = yieldPP.val/(weight_ppPR);
				yieldPbPb_PR = yieldPbPb.val/(weight_PbPbPR);
				err_PPPR = yieldPP.err;
				err_PbPbPR = yieldPbPb.err;
			}
			else if(isPR == 1)
			{
				yieldPP_NP = yieldPP.val/(weight_ppNP);
				yieldPbPb_NP = yieldPbPb.val/(weight_PbPbNP);
				err_PPNP = yieldPP.err;
				err_PbPbNP = yieldPbPb.err;
			}
		}

		double err_PbPbPR_wgt = err_PbPbPR/weight_PbPbPR;
        double err_PbPbNP_wgt = err_PbPbNP/weight_PbPbNP;
        double err_PPPR_wgt = err_PPPR/weight_ppPR;
        double err_PPNP_wgt = err_PPNP/weight_ppNP;


		cout << "Yield pp PR :\t"   << yieldPP_PR   << " +/- " << err_PPPR   << "\tEff :\t" << eff_ppPR   << "\tAcc :\t" << acc_ppPR   << "\tweight :\t" << weight_ppPR << endl;
		cout << "Yield PbPb PR :\t" << yieldPbPb_PR << " +/- " << err_PbPbPR << "\tEff :\t" << eff_PbPbPR << "\tAcc :\t" << acc_PbPbPR << "\tweight :\t" << weight_ppNP << endl;
		cout << "Yield pp NP :\t"   << yieldPP_NP   << " +/- " << err_PPNP   << "\tEff :\t" << eff_ppNP   << "\tAcc :\t" << acc_ppNP   << "\tweight :\t" << weight_PbPbPR << endl;
		cout << "Yield PbPb NP :\t" << yieldPbPb_NP << " +/- " << err_PbPbNP << "\tEff :\t" << eff_PbPbNP << "\tAcc :\t" << acc_PbPbNP << "\tweight :\t" << weight_PbPbNP << endl;

		
		hyieldPP_PR -> SetBinContent(i+1,yieldPP_PR);
		//hyieldPP_PR -> SetBinError(i+1,err_PPPR_wgt);
		hyieldPP_PR -> SetBinError(i+1,err_PPPR);
		hyieldPbPb_PR -> SetBinContent(i+1,yieldPbPb_PR);
		//hyieldPbPb_PR -> SetBinError(i+1,err_PbPbPR_wgt);
		hyieldPbPb_PR -> SetBinError(i+1,err_PbPbPR);

		hyieldPP_NP -> SetBinContent(i+1,yieldPP_NP);
		//hyieldPP_NP -> SetBinError(i+1,err_PPNP_wgt);
		hyieldPP_NP -> SetBinError(i+1,err_PPNP);
		hyieldPbPb_NP -> SetBinContent(i+1,yieldPbPb_NP);
		//hyieldPbPb_NP -> SetBinError(i+1,err_PbPbNP_wgt);
		hyieldPbPb_NP -> SetBinError(i+1,err_PbPbNP);
		

		Xpp_PR[i] = lumi_pp_scale*yieldPP_PR/(lumi_pp*1e+2*(double)(40.-9.)*2*(double)(2));
		//Xpp_PR_err[i] = lumi_pp_scale*Xpp_PR[i]*sqrt(TMath::Power(err_PPPR_wgt/yieldPP_PR,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_PR_err[i] = Xpp_PR[i]*sqrt(TMath::Power(err_PPPR_wgt/yieldPP_PR,2));// + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_NP[i] = lumi_pp_scale*yieldPP_NP/(lumi_pp*1e+2*(double)(40.-9.)*2*(double)(2));
		//Xpp_NP_err[i] = lumi_pp_scale*Xpp_NP[i]*sqrt(TMath::Power(err_PPNP_wgt/yieldPP_NP,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
		Xpp_NP_err[i] = Xpp_NP[i]*sqrt(TMath::Power(err_PPNP_wgt/yieldPP_NP,2));// + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));

		XPbPb_PR[i] = yieldPbPb_PR/(Nmb*Taa[i]*(double)(40.-9.)*(double)2*(double)(2)*cfrac[i]);
		//XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(Taa_err[i]/Taa[i],2) + TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2));// + TMath::Power(Nmb_err/Nmb,2));
		XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2));// + TMath::Power(Nmb_err/Nmb,2));
		XPbPb_NP[i] = yieldPbPb_NP/(Nmb*Taa[i]*(double)(40.-9.)*(double)2*(double)(2)*cfrac[i]);
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
/*
	for(int i=0;i<nCentBins;i++)
	{
		RaaPR[i] = hRaa_PbPb_PR->GetBinContent(i+1);
		RaaPR_err[i] = hRaa_PbPb_PR->GetBinError(i+1);
		RaaNP[i] = hRaa_PbPb_NP->GetBinContent(i+1);
		RaaNP_err[i] = hRaa_PbPb_NP->GetBinError(i+1);
		x[i] = (centBin[i+1]+centBin[i])/2;
		binWidth[i] = (centBin[i+1]-centBin[i])/2;
		cout << "i : " << i << ", RAA Value : " << RaaPR[i] << " x : " << x[i] << ", bin Width : " << binWidth[i] << endl;
		cout << "     , RaaNP_err : " << RaaNP_err[i] << endl;
		//gRaaNP->SetPoint(i,(ptBin[i+1]-ptBin[i])/2,); 
		//gRaaNP->SetPointError(i,0,Xpp_NP_err); 
	}
*/
	for(int i=nCentBins-1;i>=0;i--)
    {
        RaaPR[nCentBins-1-i] = hRaa_PbPb_PR->GetBinContent(i+1);
        RaaPR_err[nCentBins-1-i] = hRaa_PbPb_PR->GetBinError(i+1);
        RaaNP[nCentBins-1-i] = hRaa_PbPb_NP->GetBinContent(i+1);
        RaaNP_err[nCentBins-1-i] = hRaa_PbPb_NP->GetBinError(i+1);
		SysPR[nCentBins-1-i] = RaaPR[nCentBins-1-i]*(hSys_PR->GetBinContent(i+1));
        SysNP[nCentBins-1-i] = RaaNP[nCentBins-1-i]*(hSys_NP->GetBinContent(i+1));

		//cout << "i : " << i << ", RAA Value : " << RaaPR[i] << ",   RAA Err : " << RaaPR_err[i] << " x : " << x[i] << ", bin Width : " << binWidth[i] << endl;
        //cout << "   , RAA NP : " << RaaNP[i] <<  ", RaaNP_err : " << RaaNP_err[i] << endl;
		cout << " " << endl;
		cout << "Cent. : " << centBin[i] << " - " << centBin[i+1] << endl;
		cout << "Raa PR : " << RaaPR[nCentBins-1-i] << " , PR error : " << RaaPR_err[nCentBins-1-i] << endl;
		cout << "Raa NP : " << RaaNP[nCentBins-1-i] << " , NP error : " << RaaNP_err[nCentBins-1-i] << endl;

    }

	 cout << " " << endl;
	for (int i=0; i<nCentBins; i++){
		//cout << "RAA PR : " << RaaPR[i] << ", RAA NP: " << RaaNP[i] << endl;
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
	cout << "\\multicolumn{2}{c}{$|y| < 1.6$, $6.5 < \\pt < 40$ (GeV/c)} \\\\ \\hline" << endl;
	for(int i=0;i<nCentBins;i++)
	{
		cout << centBin[i] << " -- " << centBin[i+1] << " & " << hRaa_PbPb_PR->GetBinContent(i+1) << " $\\pm$ " << hRaa_PbPb_PR->GetBinError(i+1) << " (stat.) $\\pm$ " << (hRaa_PbPb_PR->GetBinContent(i+1))*(hSys_PR->GetBinContent(i+1)) << " (syst.) \\\\" << endl;
	}
	cout << " " << endl;
	cout << "\\multicolumn{2}{c}{$|y| < 1.6$, $6.5 < \\pt < 40$ (GeV/c} \\\\ \\hline" << endl;
	for(int i=0;i<nCentBins;i++)
	{
		//cout << centBin[i] << " -- " << centBin[i+1] << " & " << RaaNP[i] << " $\\pm$ " << RaaNP_err[i] << " (stat.) $\\pm$ " << SysNP[i] << " (syst.) \\\\" << endl;
		cout << centBin[i] << " -- " << centBin[i+1] << " & " << hRaa_PbPb_NP->GetBinContent(i+1) << " $\\pm$ " << hRaa_PbPb_NP->GetBinError(i+1) << " (stat.) $\\pm$ " << (hRaa_PbPb_NP->GetBinContent(i+1))*(hSys_NP->GetBinContent(i+1)) << " (syst.) \\\\" << endl;
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
	jumSun(0,1,100,1);

	drawText("Prompt J/#psi", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 2", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXPR,iPeriod,iPos);	

	cXPR->SaveAs("./figs/CrossSection_PR_y0_2_Cent.pdf");

	TCanvas *cXNP = new TCanvas("cXNP", "", 700,700);
	cXNP->cd();
	cXNP->SetLogy();
	gXpp_NP->GetXaxis()->SetTitle("Centrality (%)");
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

	drawText("Non Prompt J/#psi", pos_x, pos_y, text_color, text_size*1.3);
	drawText("|y| < 2", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXNP,iPeriod,iPos);	

	cXNP->SaveAs("figs/CrossSection_NP_y0_2_Cent.pdf");

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
/*
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
	*/
	cRAA->cd();
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

	drawText("9 < p_{T} < 40 GeV/c", pos_x, pos_y, text_color, text_size);
	drawText("|y| < 2", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cRAA,iPeriod,iPos);	

	//pad1_2->cd();
	/*
	// Calculate RAA for full pT range (6.5-40 GeV/c, 0-180%)
	double weight_ppPR_int = 1; double weight_ppNP_int = 1;
	double weight_PbPbPR_int = 1; double weight_PbPbNP_int = 1;
	
	double eff_ppPR = hEff_ppPR->GetBinContent(1);
	double eff_ppNP = hEff_ppNP->GetBinContent(1);
	double acc_ppPR = hAcc_ppPR->GetBinContent(1);
	double acc_ppNP = hAcc_ppNP->GetBinContent(1);
	
	double eff_PbPbPR = hEff_PbPbPR_Integrated->GetBinContent(1);
	double eff_PbPbNP = hEff_PbPbNP_Integrated->GetBinContent(1);
	double acc_PbPbPR = hAcc_PbPbPR->GetBinContent(1);
	double acc_PbPbNP = hAcc_PbPbNP->GetBinContent(1);
	
	valErr yieldPP_PR_int = getYield_pp(0, 0);
	valErr yieldPP_NP_int = getYield_pp(0, 1);
	valErr yieldPbPb_PR_int = getYield_PbPb_Integrated(0);
	valErr yieldPbPb_NP_int = getYield_PbPb_Integrated(1);
	
	//valErr ctauEffPP_PR_int = getCtauEff_pp(0, 0);
	//valErr ctauEffPP_NP_int = getCtauEff_pp(0, 1);
	//valErr ctauEffPbPb_PR_int = getCtauEff_PbPb_Integrated(0);
	//valErr ctauEffPbPb_NP_int = getCtauEff_PbPb_Integrated(1);
	
	weight_ppPR_int = eff_ppPR * acc_ppPR;// * ctauEffPP_PR_int.val;
	weight_ppNP_int = eff_ppNP * acc_ppNP;// * ctauEffPP_NP_int.val;
	weight_PbPbPR_int = eff_PbPbPR * acc_PbPbPR;// * ctauEffPbPb_PR_int.val;
	weight_PbPbNP_int = eff_PbPbNP * acc_PbPbNP;// * ctauEffPbPb_NP_int.val;
	
	double yieldPP_PR_corr = yieldPP_PR_int.val / weight_ppPR_int;
	double yieldPP_NP_corr = yieldPP_NP_int.val / weight_ppNP_int;
	double yieldPbPb_PR_corr = yieldPbPb_PR_int.val / weight_PbPbPR_int;
	double yieldPbPb_NP_corr = yieldPbPb_NP_int.val / weight_PbPbNP_int;
	
	double err_PPPR_int = yieldPP_PR_int.err / weight_ppPR_int;
	double err_PPNP_int = yieldPP_NP_int.err / weight_ppNP_int;
	double err_PbPbPR_int = yieldPbPb_PR_int.err / weight_PbPbPR_int;
	double err_PbPbNP_int = yieldPbPb_NP_int.err / weight_PbPbNP_int;
	
	double Taa_int = 6.274; // 0-90%
	double Xpp_PR_int = lumi_pp_scale * yieldPP_PR_corr / (lumi_pp*1e+2*(double)(40.-6.5)*2*(double)(1.6-0));
	double Xpp_NP_int = lumi_pp_scale * yieldPP_NP_corr / (lumi_pp*1e+2*(double)(40.-6.5)*2*(double)(1.6-0));
	double XPbPb_PR_int = yieldPbPb_PR_corr / (Nmb*Taa_int*(double)(40.-6.5)*(double)2*(1.6-0));
	double XPbPb_NP_int = yieldPbPb_NP_corr / (Nmb*Taa_int*(double)(40.-6.5)*(double)2*(1.6-0));
	
	double Xpp_PR_int_err = Xpp_PR_int * sqrt(TMath::Power(err_PPPR_int/yieldPP_PR_corr,2));
	double Xpp_NP_int_err = Xpp_NP_int * sqrt(TMath::Power(err_PPNP_int/yieldPP_NP_corr,2));
	double XPbPb_PR_int_err = XPbPb_PR_int * sqrt(TMath::Power(err_PbPbPR_int/yieldPbPb_PR_corr,2));
	double XPbPb_NP_int_err = XPbPb_NP_int * sqrt(TMath::Power(err_PbPbNP_int/yieldPbPb_NP_corr,2));
	
	double RaaPR_int = XPbPb_PR_int / Xpp_PR_int;
	double RaaNP_int = XPbPb_NP_int / Xpp_NP_int;
	double RaaPR_int_err = RaaPR_int * sqrt(TMath::Power(XPbPb_PR_int_err/XPbPb_PR_int,2) + TMath::Power(Xpp_PR_int_err/Xpp_PR_int,2));
	double RaaNP_int_err = RaaNP_int * sqrt(TMath::Power(XPbPb_NP_int_err/XPbPb_NP_int,2) + TMath::Power(Xpp_NP_int_err/Xpp_NP_int,2));
	
	cout << " " << endl;
	cout << "Integrated RAA (6.5 < pT < 40 GeV/c, 0-180%):" << endl;
	cout << "Prompt: " << RaaPR_int << " +/- " << RaaPR_int_err << endl;
	cout << "NonPrompt: " << RaaNP_int << " +/- " << RaaNP_int_err << endl;
	cout << " " << endl;
	
	TH1D *hmidPR_int = new TH1D("hmidPR_int",";;",3,0,1);
	TH1D *hmidNP_int = new TH1D("hmidNP_int",";;",3,0,1);

	hmidPR_int -> SetBinContent(2,RaaPR_int);
	hmidPR_int -> SetBinError(2,RaaPR_int_err);
	hmidNP_int -> SetBinContent(2,RaaNP_int);
	hmidNP_int -> SetBinError(2,RaaNP_int_err);

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
	hmidNP_int->Draw("PE SAME");*/
	jumSun(0,1,1,1);
	drawText("Cent. 0-90%", 0.12, 0.76, text_color, 20);

	cRAA->SaveAs(Form("./figs/RAA_JPsi_y0_2_ATLAS_Npart_Sys%d.pdf",isSys));

	TFile *outFile = new TFile("./roots/RAA_JPsi_ATLAS_Comp_Npart.root","recreate");
	hRAA_PR->Write();
	hRAA_NP->Write();
	//hmidPR_int->Write();
	//hmidNP_int->Write();

	outFile->Close();

}

valErr getYield_pp(int i, int isPR){
    TString kineLabel;
    kineLabel = getKineLabelpp(9,40,0,2,0.0);
	TString PR;
	if(isPR==0) PR = "PRMC";
	else if(isPR==1) PR = "NPMC";
    TFile* inf = new TFile(Form("../../Jpsi_L_cut/roots_1S_pp/%s/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", PR.Data(), kineLabel.Data()));
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");

    valErr ret; 
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getYield_PbPb(int i, int isPR){
    double centBin[7] = {0,10,20,40,80,120,160};
    TString kineLabel[7];
	TString PR;
	if(isPR==0) PR = "PRMC";
	else if(isPR==1) PR = "NPMC";
    kineLabel[i] = getKineLabel(9,40,0,2,0.0,centBin[i],centBin[i+1]);
	TFile* inf = new TFile(Form("../../Jpsi_L_cut/roots_1S_Pb/%s/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", PR.Data(), kineLabel[i].Data()));
    //TFile* inf = new TFile(Form("./Jpsi/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
/*
valErr getYield_PbPb_Integrated(int isPR){
	TString kineLabel;
	TString PR;
	if(isPR==0) PR = "PRMC";
	else if(isPR==1) PR = "NPMC";
	kineLabel = getKineLabel(9,40,0,2,0.0,0,180);
	TFile* inf = new TFile(Form("../../Jpsi_L_cut/roots_1S_Pb/%s/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", PR.Data(), kineLabel.Data()));
	TH1D* fitResults = (TH1D*)inf->Get("fitResults");
	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
valErr getCtauEff_PbPb_Integrated(int isPR){
	TString kineLabel;
	TString PR;
	if(isPR==0) PR = "PRMC";
	else if(isPR==1) PR = "NPMC";
	kineLabel = getKineLabel(9,40,0,2,0.0,0,180);
	TFile* f_ctau = new TFile(Form("../../../ctau_eff/roots_1S_Pb/decayL/%s/decay_hist_%s_m2.6-3.5.root",PR.Data(), kineLabel.Data()));
	TH1D* h_ctau = (TH1D*)f_ctau->Get("h_eff");
	valErr ret;
	ret.val = h_ctau->GetBinContent(1);
	ret.err = h_ctau->GetBinError(1);
	return ret;
}
*/