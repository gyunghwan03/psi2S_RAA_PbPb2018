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
using namespace TMath;

valErr getYield_PbPb(int i=0, double ptHigh=40);
valErr getYield_pp(int i=0, double ptHigh=40);
valErr getFrac_PbPb(int i=0, double ptHigh=40);
valErr getFrac_pp(int i=0, double ptHigh=40);
double getAccWeight(TH1D* h = 0, double pt = 0);
double getEffWeight(TH1D* h = 0, double pt = 0);

void draw_Raa_psi2S_y1p6_2p4_pT(bool isSys = true, double ptHigh = 40)
{
	//gROOT->SumW2();
    gStyle->SetOptStat(0);
    setTDRStyle();
    int iPeriod = 101;
    int iPos = 33;

    const int nPtBins=4;
    TFile *fPbPb[nPtBins+1];
    TFile *fpp[nPtBins+1];

    //TFile *fEff_PbPbPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20230602.root");
    //TFile *fEff_PbPbNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nonprompt_pbpb_psi2s_PtW1_tnp1_20230602.root");
    //TFile *fEff_ppPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230718.root");
    //TFile *fEff_ppNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_nprompt_pp_psi2S_PtW1_tnp1_20230801.root");
    //TFile *fAcc_pp = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20230728.root");
    //TFile *fAcc_PbPb = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_PbPb_SysUp0_20230728.root");
    //TFile *fSys = new TFile("../syst_summary/syst_roots/total_syst.root");
    //TFile *fSys = new TFile("../syst_summary/syst_roots/total_syst_NobFrac.root");
    //TFile *fEff_PbPbPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20231226.root");
    //TFile *fEff_PbPbNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_nprompt_pbpb_psi2s_PtW1_tnp1_20231226.root");
    //TFile *fEff_ppPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20231226.root");
    //TFile *fEff_ppNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_nprompt_pp_psi2S_PtW1_tnp1_20231226.root");
    //TFile *fAcc_ppPR = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20231226.root");
    //TFile *fAcc_PbPbPR = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_PbPb_SysUp0_20231226.root");
    //TFile *fAcc_ppNP = new TFile("../../Eff_Acc/roots/acceptance_NonPrompt_psi2s_GenOnly_wgt1_pp_SysUp0_20240514.root");
    //TFile *fAcc_PbPbNP = new TFile("../../Eff_Acc/roots/acceptance_NonPrompt_psi2s_GenOnly_wgt1_PbPb_SysUp0_20240515.root");
	TFile *fEff_PbPbPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20240807_ppWfunc.root");
    TFile *fEff_PbPbNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_cent_0_to_180_rap_prompt_pbpb_psi2s_PtW1_tnp1_20240807_ppWfunc.root");
    TFile *fEff_ppPR = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_prompt_pp_psi2s_PtW1_tnp1_20230728.root");
    TFile *fEff_ppNP = new TFile("../../Eff_Acc/roots/mc_eff_vs_pt_rap_nprompt_pp_psi2s_PtW1_tnp1_20230801.root");
    TFile *fAcc_ppPR = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20230728.root");
    TFile *fAcc_PbPbPR = new TFile("../../Eff_Acc/roots/acceptance_Prompt_psi2s_GenOnly_wgt1_pp_SysUp0_20230728.root");
    TFile *fAcc_ppNP = new TFile("../../Eff_Acc/roots/acceptance_NonPrompt_psi2s_GenOnly_wgt1_pp_SysUp0_20240514.root");
    TFile *fAcc_PbPbNP = new TFile("../../Eff_Acc/roots/acceptance_NonPrompt_psi2s_GenOnly_wgt1_pp_SysUp0_20240514.root");
    TFile *fSys = new TFile("../syst_summary_psi2S/syst_roots/total_syst.root");

    TH1D *hEff_PbPbPR = (TH1D*) fEff_PbPbPR -> Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4");
    TH1D *hEff_PbPbNP = (TH1D*) fEff_PbPbNP -> Get("mc_eff_vs_pt_TnP1_PtW1_cent_0_to_180_absy1p6_2p4");
    TH1D *hEff_ppPR = (TH1D*) fEff_ppPR -> Get("mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4");
    TH1D *hEff_ppNP = (TH1D*) fEff_ppNP -> Get("mc_eff_vs_pt_TnP1_PtW1_absy1p6_2p4");
    TH1D *hAcc_PbPbPR = (TH1D*) fAcc_PbPbPR -> Get("hAccPt_2021_Fory");
    TH1D *hAcc_PbPbNP = (TH1D*) fAcc_PbPbNP -> Get("hAccPt_2021_Fory");
    TH1D *hAcc_ppPR = (TH1D*) fAcc_ppPR -> Get("hAccPt_2021_Fory");
    TH1D *hAcc_ppNP = (TH1D*) fAcc_ppNP -> Get("hAccPt_2021_Fory");

	TH1D *hSys_PR = (TH1D*) fSys -> Get("mid_pt_PR");
    TH1D *hSys_PR_pb = (TH1D *)fSys->Get("mid_pt_PR_pb");
    TH1D *hSys_PR_pp = (TH1D *)fSys->Get("mid_pt_PR_pp");
    TH1D *hSys_NP = (TH1D *)fSys->Get("mid_pt_NP");
    TH1D *hSys_NP_pb = (TH1D *)fSys->Get("mid_pt_NP_pb");
    TH1D *hSys_NP_pp = (TH1D *)fSys->Get("mid_pt_NP_pp");

    Double_t Nmb = 11968044281.;
    Double_t Nmb_err = 0.01261;
    //Double_t Taa = 5.649; // 0-100%
    Double_t Taa = 6.274; // 0-90%
    Double_t lumi_pp = 3.002;
    Double_t lumi_pp_scale = 1e-9;
    Double_t lumi_pp_err = 1*lumi_pp_scale;
    Double_t Taa_err = 0.137;


    double ptBin[nPtBins+1] = {3.5,6.5,9,12,ptHigh};
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
	double SysPR[nPtBins]; double SysNP[nPtBins];
    double SysPR_pb[nPtBins]; double SysNP_pb[nPtBins];
    double SysPR_pp[nPtBins]; double SysNP_pp[nPtBins];

    double Xpp_PR[nPtBins]; double Xpp_NP[nPtBins];
    double Xpp_PR_err[nPtBins]; double Xpp_NP_err[nPtBins];
    double XPbPb_PR[nPtBins]; double XPbPb_NP[nPtBins];
    double XPbPb_PR_err[nPtBins]; double XPbPb_NP_err[nPtBins];


    for(int i=0;i<nPtBins;i++)
    {
        double weight_ppPR = 1; double weight_ppNP = 1; double weight_PbPbPR = 1; double weight_PbPbNP = 1;
        double eff_ppPR = 1; double eff_ppNP = 1; double eff_PbPbPR = 1; double eff_PbPbNP = 1;
        double acc_ppPR = 1; double acc_PbPbPR = 1; double acc_ppNP = 1; double acc_PbPbNP = 1;


        eff_ppPR=hEff_ppPR->GetBinContent(i+2);
        eff_ppNP=hEff_ppNP->GetBinContent(i+2);
        acc_ppPR=hAcc_ppPR->GetBinContent(i+1);
        acc_ppNP=hAcc_ppNP->GetBinContent(i+1);
        eff_PbPbPR=hEff_PbPbPR->GetBinContent(i+2);
        eff_PbPbNP=hEff_PbPbNP->GetBinContent(i+2);
        acc_PbPbPR=hAcc_PbPbPR->GetBinContent(i+1);
        acc_PbPbNP=hAcc_PbPbNP->GetBinContent(i+1);

        weight_ppPR=eff_ppPR*acc_ppPR;
        weight_ppNP=eff_ppNP*acc_ppNP;
        weight_PbPbPR=eff_PbPbPR*acc_PbPbPR;
        weight_PbPbNP=eff_PbPbNP*acc_PbPbNP;

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
        cout << "pT " << ptBin[i] << "-" << ptBin[i+1] << endl;
        cout << "PbPb Yield : " << yieldPbPb.val << ", yield Err : " << err1_PbPb << ", b frac : " << fracPbPb.val << ", frac Err : " << err2_PbPb << ",    Eff PR : " << eff_PbPbPR << ",  Eff NP : " << eff_PbPbNP << ", Acc PR : " << acc_PbPbPR << ",   Acc NP :    " << acc_PbPbNP << endl;
        cout << "pp Yield   : " << yieldPP.val << ", yield Err : " << err1_PP << ", b frac : " << fracPP.val << ", frac Err : " << err2_PP << ",    Eff PR : " << eff_ppPR << ",    Eff NP : " << eff_ppNP <<  ", Acc PR : " << acc_ppPR  << ", Acc NP :    " << acc_ppNP << endl;
        //cout << "pp Prompt yield : " << yieldPP_PR << ", PbPb Prompt yield : " << yieldPbPb_PR << ", pp NonPrompt yield : " << yieldPP_NP << ", PbPb NonPrompt yield : " << yieldPbPb_NP << endl;
        //cout << "pp Prompt error : " << err_PPPR_wgt << " , PbPb Prompt error : " << err_PbPbPR_wgt << " , pp NonPrompt error : " << err_PPNP_wgt << " , PbPb NonPrompt error : " << err_PbPbNP_wgt << endl;

        hyieldPP_PR -> SetBinContent(i+1,yieldPP_PR);
        hyieldPP_PR -> SetBinError(i+1,err_PPPR_wgt);
        hyieldPbPb_PR -> SetBinContent(i+1,yieldPbPb_PR);
        hyieldPbPb_PR -> SetBinError(i+1,err_PbPbPR_wgt);

        hyieldPP_NP -> SetBinContent(i+1,yieldPP_NP);
        hyieldPP_NP -> SetBinError(i+1,err_PPNP_wgt);
        hyieldPbPb_NP -> SetBinContent(i+1,yieldPbPb_NP);
        hyieldPbPb_NP -> SetBinError(i+1,err_PbPbNP_wgt);

        cout << "PbPb PR Yield  : " << yieldPbPb.val*(1-fracPbPb.val) << ", PbPb NP Yield   : " << yieldPbPb.val*fracPbPb.val << endl;
        cout << "pp PR Yield    : " << yieldPP.val*(1-fracPP.val) << ", pp NP Yield : " << yieldPP.val*fracPP.val << endl;

        Xpp_PR[i] = lumi_pp_scale*yieldPP_PR/(lumi_pp*1e+2*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
        Xpp_PR_err[i] = Xpp_PR[i]*sqrt(TMath::Power(err_PPPR_wgt/yieldPP_PR,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));
        Xpp_NP[i] = lumi_pp_scale*yieldPP_NP/(lumi_pp*1e+2*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
        Xpp_NP_err[i] = lumi_pp_scale*Xpp_NP[i]*sqrt(TMath::Power(err_PPNP_wgt/yieldPP_NP,2) + TMath::Power(lumi_pp_err/(lumi_pp*1e+2),2));

        XPbPb_PR[i] = yieldPbPb_PR/(Nmb*Taa*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
        //XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(Taa_err/Taa,2) + TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2) + TMath::Power(Nmb_err,2));
        XPbPb_PR_err[i] = XPbPb_PR[i]*sqrt(TMath::Power(Taa_err/Taa,2) + TMath::Power(err_PbPbPR_wgt/yieldPbPb_PR,2));// + TMath::Power(Nmb_err,2));
        XPbPb_NP[i] = yieldPbPb_NP/(Nmb*Taa*(ptBin[i+1]-ptBin[i])*(double)2*(2.4-1.6));
        XPbPb_NP_err[i] = XPbPb_NP[i]*sqrt(TMath::Power(Taa_err/Taa,2) + TMath::Power(err_PbPbNP_wgt/yieldPbPb_NP,2));// + TMath::Power(Nmb_err,2));

        hXpp_PR->SetBinContent(i+1, Xpp_PR[i]);
        hXpp_NP->SetBinContent(i+1, Xpp_NP[i]);
        hXPbPb_PR->SetBinContent(i+1, XPbPb_PR[i]);
        hXPbPb_NP->SetBinContent(i+1, XPbPb_NP[i]);

        hXpp_PR->SetBinError(i+1, Xpp_PR_err[i]);
        hXpp_NP->SetBinError(i+1, Xpp_NP_err[i]);
        hXPbPb_PR->SetBinError(i+1, XPbPb_PR_err[i]);
        hXPbPb_NP->SetBinError(i+1, XPbPb_NP_err[i]);

        //cout << "i = " << i << " , Xpp_PR_err : " << Xpp_PR_err[i] << " , XPbPb_PR_err : " << XPbPb_PR_err[i] << endl;
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

    for(int i=0;i<nPtBins;i++)
    {
        RaaPR[i] = hRaa_PbPb_PR->GetBinContent(i+1);
        //RaaPR_err[i] = hRaa_PbPb_PR->GetBinError(i+1);
        RaaPR_err[i] = RaaPR[i]*sqrt(TMath::Power( (hyieldPP_PR->GetBinError(i+1)) / (hyieldPP_PR->GetBinContent(i+1)) ,2) + TMath::Power( (hyieldPbPb_PR->GetBinError(i+1)) / (hyieldPbPb_PR->GetBinContent(i+1)),2));
        RaaNP[i] = hRaa_PbPb_NP->GetBinContent(i+1);
        RaaNP_err[i] = RaaNP[i]*sqrt(TMath::Power( (hyieldPP_NP->GetBinError(i+1)) / (hyieldPP_NP->GetBinContent(i+1)) ,2) + TMath::Power( (hyieldPbPb_NP->GetBinError(i+1)) / (hyieldPbPb_NP->GetBinContent(i+1)),2));
        //RaaNP_err[i] = hRaa_PbPb_NP->GetBinError(i+1);
        x[i] = (ptBin[i+1]+ptBin[i])/2;
        binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
        cout << "i : " << i << ", RAA PR : " << RaaPR[i] << ",   RaaPR Err : " << RaaPR_err[i] << " x : " << x[i] << ", bin Width : " << binWidth[i] << endl;
        cout << "       RAA NP : " << RaaNP[i] <<  ",   RaaNP_err : " << RaaNP_err[i] << endl;
        //gRaaNP->SetPoint(i,(ptBin[i+1]-ptBin[i])/2,);
        //gRaaNP->SetPointError(i,0,Xpp_NP_err);
		SysPR[i] = RaaPR[i]*(hSys_PR->GetBinContent(i+1));
        SysNP[i] = RaaNP[i]*(hSys_NP->GetBinContent(i+1));
        SysPR_pb[i] = XPbPb_PR[i] * (hSys_PR_pb->GetBinContent(i + 1));
        SysNP_pb[i] = XPbPb_NP[i] * (hSys_NP_pb->GetBinContent(i + 1));
        SysPR_pp[i] = Xpp_PR[i] * (hSys_PR_pp->GetBinContent(i + 1));
        SysNP_pp[i] = Xpp_NP[i] * (hSys_NP_pp->GetBinContent(i + 1));
        hRAA_PR->SetBinContent(i+1,RaaPR[i]);
        hRAA_PR->SetBinError(i+1,RaaPR_err[i]);
        hRAA_NP->SetBinContent(i+1,RaaNP[i]);
        hRAA_NP->SetBinError(i+1,RaaNP_err[i]);
    }


	//////////////////////////////// Tex Table ////////////////////////////////
	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;
	cout << "\\multicolumn{2}{c}{$1.6 < |y| < 2.4, 0-90\\%$} \\\\ " << endl;
	cout.precision(3);
	for(int i=0;i<nPtBins;i++)
	{
		cout << ptBin[i] << " -- " << ptBin[i+1] << " & " << RaaPR[i] << " $\\pm$ " << RaaPR_err[i] << " (stat.) $\\pm$ " << SysPR[i] << " (syst.) \\\\" << endl;
	}

	cout << " " << endl;
	cout << "\\multicolumn{2}{c}{$1.6 < |y| < 2.4, 0-90\\%$} \\\\ " << endl;
	cout.precision(3);
	for(int i=0;i<nPtBins;i++)
	{
		cout << ptBin[i] << " -- " << ptBin[i+1] << " & " << RaaNP[i] << " $\\pm$ " << RaaNP_err[i] << " (stat.) $\\pm$ " << SysNP[i] << " (syst.) \\\\" << endl;
	}

	cout << " " << endl;
	cout << " " << endl;
	cout << " " << endl;

	cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Cross Section  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"  << endl;
	cout << "\\multicolumn{2}{c}{$|y| < 1.6, 0-90\\%$} \\\\ " << endl;
	for(int i=0;i<nPtBins;i++)
	{
		cout << ptBin[i] << " -- " << ptBin[i+1] << " & " << Xpp_PR[i] << " $\\pm$ " << Xpp_PR_err[i] << " (stat.) $\\pm$ " << SysPR_pp[i] << " (syst.) \\\\" << endl;
	}
	cout << " " << endl;
	for(int i=0;i<nPtBins;i++)
	{
		cout << ptBin[i] << " -- " << ptBin[i+1] << " & " << Xpp_NP[i] << " $\\pm$ " << Xpp_NP_err[i] << " (stat.) $\\pm$ " << SysNP_pp[i] << " (syst.) \\\\" << endl;
	}
	cout << " " << endl;
	cout << "\\multicolumn{2}{c}{$|y| < 1.6, 0-90\\%$} \\\\ " << endl;
	for(int i=0;i<nPtBins;i++)
	{
		cout << ptBin[i] << " -- " << ptBin[i+1] << " & " << XPbPb_PR[i] << " $\\pm$ " << XPbPb_PR_err[i] << " (stat.) $\\pm$ " << SysPR_pb[i] << " (syst.) \\\\" << endl;
	}
	cout << " " << endl;
	for(int i=0;i<nPtBins;i++)
	{
		cout << ptBin[i] << " -- " << ptBin[i+1] << " & " << XPbPb_NP[i] << " $\\pm$ " << XPbPb_NP_err[i] << " (stat.) $\\pm$ " << SysNP_pb[i] << " (syst.) \\\\" << endl;
	}
	//////////////////////////////// Start Plotting ////////////////////////////////
	float pos_x = 0.20;
	float pos_y = 0.85;
	float pos_y_diff = 0.051;
	int text_color = 1;
	float text_size = 25;

    TGraphErrors *gXpp_PR;
    TGraphErrors *gXpp_NP;
    TGraphErrors *gXPbPb_PR;
    TGraphErrors *gXPbPb_NP;
    if (isSys == 1)
    {
        gXpp_PR = new TGraphErrors(nPtBins, x, Xpp_PR, 0, Xpp_PR_err);
        gXpp_NP = new TGraphErrors(nPtBins, x, Xpp_NP, 0, Xpp_NP_err);

        gXPbPb_PR = new TGraphErrors(nPtBins, x, XPbPb_PR, 0, XPbPb_PR_err);
        gXPbPb_NP = new TGraphErrors(nPtBins, x, XPbPb_NP, 0, XPbPb_NP_err);
    }
    else
    {
        gXpp_PR = new TGraphErrors(nPtBins, x, Xpp_PR, binWidth, Xpp_PR_err);
        gXpp_NP = new TGraphErrors(nPtBins, x, Xpp_NP, binWidth, Xpp_NP_err);

        gXPbPb_PR = new TGraphErrors(nPtBins, x, XPbPb_PR, binWidth, XPbPb_PR_err);
        gXPbPb_NP = new TGraphErrors(nPtBins, x, XPbPb_NP, binWidth, XPbPb_NP_err);
    }

    TGraphErrors *gSysPR_pb = new TGraphErrors(nPtBins, x, XPbPb_PR, binWidth, SysPR_pb);
    TGraphErrors *gSysNP_pb = new TGraphErrors(nPtBins, x, XPbPb_NP, binWidth, SysNP_pb);
    TGraphErrors *gSysPR_pp = new TGraphErrors(nPtBins, x, Xpp_PR, binWidth, SysPR_pp);
    TGraphErrors *gSysNP_pp = new TGraphErrors(nPtBins, x, Xpp_NP, binWidth, SysNP_pp);

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
	gXpp_PR->GetXaxis()->SetLimits(0.,ptHigh);
	gXpp_PR->SetMinimum(1e-11);
	gXpp_PR->SetMaximum(1e-5);
	gXpp_PR->SetMarkerColor(kBlue+2);
	gXpp_PR->SetLineColor(kBlue+2);
	gXpp_PR->SetMarkerStyle(20);
	gXpp_PR->SetMarkerSize(1.2);
	gXPbPb_PR->SetMarkerColor(kRed+3);
    gXPbPb_PR->SetLineColor(kRed+3);
    gXPbPb_PR->SetMarkerStyle(21);
    gXPbPb_PR->SetMarkerSize(1.2);

    gSysPR_pb->GetXaxis()->SetLimits(0., ptHigh);
    gSysPR_pb->SetMinimum(0);
    gSysPR_pb->SetMaximum(1.44);
    gSysPR_pb->SetLineColor(kRed - 4);
    gSysPR_pb->SetFillColorAlpha(kRed - 9, 0.40);

    gSysPR_pp->GetXaxis()->SetLimits(0., ptHigh);
    gSysPR_pp->SetMinimum(0);
    gSysPR_pp->SetMaximum(1.44);
    gSysPR_pp->SetLineColor(kBlue - 4);
    gSysPR_pp->SetFillColorAlpha(kBlue - 9, 0.40);

    gXpp_PR->Draw("AP");
	gXPbPb_PR->Draw("P");
    if (isSys == 1)
    {
        gSysPR_pb->Draw("5");
        gSysPR_pp->Draw("5");
    }


	TLegend *legPR = new TLegend(0.66,0.75,0.8,0.85);
    SetLegendStyle(legPR);
	legPR->AddEntry(gXpp_PR,"pp", "pe");
	legPR->AddEntry(gXPbPb_PR,"PbPb, Cent. 0-90%", "pe");
	legPR->Draw("SAME");
	jumSun(0,1,ptHigh,1);

	drawText("Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXPR,iPeriod,iPos);	

	cXPR->SaveAs(Form("./figs/CrossSection_PR_y1p6_2p4_pT_%.f.pdf",ptHigh));

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
	gXpp_NP->GetXaxis()->SetLimits(0.,ptHigh);
	gXpp_NP->SetMinimum(1e-11);
	gXpp_NP->SetMaximum(1e-5);
	gXpp_NP->SetMarkerColor(kBlue+2);
	gXpp_NP->SetLineColor(kBlue+2);
	gXpp_NP->SetMarkerStyle(20);
	gXpp_NP->SetMarkerSize(1.2);
	gXPbPb_NP->SetMarkerColor(kRed+3);
	gXPbPb_NP->SetLineColor(kRed+3);
	gXPbPb_NP->SetMarkerStyle(21);
	gXPbPb_NP->SetMarkerSize(1.2);

    gSysNP_pb->GetXaxis()->SetLimits(0., ptHigh);
    gSysNP_pb->SetMinimum(0);
    gSysNP_pb->SetMaximum(1.44);
    gSysNP_pb->SetLineColor(kRed - 4);
    gSysNP_pb->SetFillColorAlpha(kRed - 9, 0.40);

    gSysNP_pp->GetXaxis()->SetLimits(0., ptHigh);
    gSysNP_pp->SetMinimum(0);
    gSysNP_pp->SetMaximum(1.44);
    gSysNP_pp->SetLineColor(kBlue - 4);
    gSysNP_pp->SetFillColorAlpha(kBlue - 9, 0.40);

    gXpp_NP->Draw("AP");
	gXPbPb_NP->Draw("P");
    if (isSys == 1)
    {
        gSysNP_pb->Draw("5");
        gSysNP_pp->Draw("5");
    }

    TLegend *legNP = new TLegend(0.66,0.72,0.8,0.82);
    SetLegendStyle(legNP);
	legNP->AddEntry(gXpp_NP,"pp", "pe");
	legNP->AddEntry(gXPbPb_NP,"PbPb, Cent. 0-90%", "pe");
	legNP->Draw("SAME");
	jumSun(0,1,ptHigh,1);

	drawText("Non Prompt #psi(2S)", pos_x, pos_y, text_color, text_size*1.3);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cXNP,iPeriod,iPos);	

	cXNP->SaveAs(Form("./figs/CrossSection_NP_y1p6_2p4_pT_%.f.pdf",ptHigh));

    TGraphErrors *gRaaPR;
    TGraphErrors *gRaaNP;

    if(isSys==1) {
	gRaaPR = new TGraphErrors(nPtBins,x,RaaPR,0,RaaPR_err);
	gRaaNP = new TGraphErrors(nPtBins,x,RaaNP,0,RaaNP_err);
    }
    else {
	gRaaPR = new TGraphErrors(nPtBins,x,RaaPR,binWidth,RaaPR_err);
	gRaaNP = new TGraphErrors(nPtBins,x,RaaNP,binWidth,RaaNP_err);
    }
    TGraphErrors *gSysPR = new TGraphErrors(nPtBins,x,RaaPR,binWidth,SysPR);
    TGraphErrors *gSysNP = new TGraphErrors(nPtBins,x,RaaNP,binWidth,SysNP);
	TCanvas *cRAA = new TCanvas("cRAA", "", 700, 700);
	cRAA->cd();
	gRaaPR->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gRaaPR->GetXaxis()->CenterTitle();
	gRaaPR->GetYaxis()->SetTitle("R_{AA}");
	gRaaPR->SetTitle("");
	gRaaPR->GetYaxis()->CenterTitle();
	gRaaPR->GetXaxis()->SetLimits(0.,ptHigh);
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

	gSysPR->SetLineColor(kBlue-4);
    gSysPR->SetFillColorAlpha(kBlue-9,0.40);
    gSysNP->SetLineColor(kRed-4);
    gSysNP->SetFillColorAlpha(kRed-9,0.40);

	gRaaPR->Draw("AP");
	gRaaNP->Draw("P");
    if (isSys==1){
        gSysPR->Draw("5");
        gSysNP->Draw("5");
    }

	TLegend *leg = new TLegend(0.66,0.72,0.8,0.82);
    SetLegendStyle(leg);
	leg->AddEntry(gRaaPR,"Prompt #psi(2S)", "pe");
	leg->AddEntry(gRaaNP,"Non Prompt #psi(2S)", "pe");
	leg->Draw("SAME");
	jumSun(0,1,ptHigh,1);

	double glb_Err=0;
	glb_Err = TMath::Sqrt( TMath::Power(Nmb_err,2) + TMath::Power(0.022,2) + TMath:: Power(0.02,2) );
	
	cout << "Global Uncertainty : " << glb_Err << endl;

	TBox* b_ynerr = new TBox();
	b_ynerr->SetX1(0);
	b_ynerr->SetX2(2);
	b_ynerr->SetY1(1 - glb_Err);
	b_ynerr->SetY2(1 + glb_Err);
	b_ynerr->SetFillColorAlpha(12, 0.8);
	b_ynerr->SetLineWidth(1);
	b_ynerr->SetLineColor(kBlack);
	b_ynerr->Draw("L");

	drawText("Cent. 0-90%", pos_x, pos_y, text_color, text_size);
	drawText("1.6 < |y| < 2.4", pos_x, pos_y-pos_y_diff, text_color, text_size);
    CMS_lumi_v2mass(cRAA,iPeriod,iPos);	

	cRAA->SaveAs(Form("./figs/RAA_psi2S_y1p6_2p4_pT_%.f_Sys%d.pdf",ptHigh,isSys));

	TFile *f1 = new TFile(Form("./roots/RAA_psi2S_forRap_pT_%.f.root",ptHigh),"recreate");
	f1->cd();
	hRAA_PR->Write();
	hRAA_NP->Write();
	f1->Close();

}

valErr getYield_pp(int i, double ptHigh){
    double ptBins[5] = {3.5,6.5,9,12,ptHigh};
    TString kineLabel[5];
    kineLabel[i] = getKineLabelpp(ptBins[i],ptBins[i+1],1.6,2.4,0.0);
    TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    //TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/backup_231122/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("fitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getYield_PbPb(int i, double ptHigh){
    double ptBins[5] = {3.5,6.5,9,12,ptHigh};
    TString kineLabel[5];
    kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
    //TFile* inf = new TFile(Form("./psi2S/roots/2DFit_No_Weight/Mass/Mass_FixedFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TFile* inf = new TFile(Form("../psi2S_230512/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("MassResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getFrac_PbPb(int i, double ptHigh) {
    double ptBin[5] = {3.5,6.5,9,12,ptHigh};
    TString kineLabel[5];
    kineLabel[i] = getKineLabel(ptBin[i],ptBin[i+1],1.6,2.4,0.0,0,180);
    TFile* inf = new TFile(Form("../psi2S_230512/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

    valErr ret;
    ret.val = fitResults->GetBinContent(1);
    ret.err = fitResults->GetBinError(1);
    return ret;
}
valErr getFrac_pp(int i, double ptHigh) {
    double ptBin[5] = {3.5,6.5,9,12,ptHigh};
    TString kineLabel[5];
    kineLabel[i] = getKineLabelpp(ptBin[i],ptBin[i+1],1.6,2.4,0.0);
    TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
    //TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/backup_231122/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
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
