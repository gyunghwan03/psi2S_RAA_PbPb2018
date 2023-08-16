#include <iostream>
#include "../rootFitHeaders.h"
#include "../commonUtility.h"
#include "../JpsiUtility.h"
#include "../cutsAndBin.h"
#include "../CMS_lumi_v2mass.C"
#include "../tdrstyle.C"
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

valErr getYield_Jpsi(int i=0);
valErr getYield_psi2S(int i=0);
valErr getFrac_Jpsi(int i=0);
valErr getFrac_psi2S(int i=0);
void DoubleRatio_Charmonia()
{
	gStyle->SetOptStat(0);

	TString DATE="221116";
	TFile *fJpsi[7];
	TFile *fpsi2S[7];

	double ptBins[7] = {3,6.5,9,12,15,20,50};
	double fracJpsi[7]; double fracErrJpsi[7]; double fracpsi2S[7]; double fracErrpsi2S[7];

	TString kineLabel[7];

	TH1D* hJpsiPR = new TH1D("hJPsiPR", ";p_{T};", 6, ptBins);
	TH1D* hJpsiNP = new TH1D("hJPsiNP", ";p_{T};", 6, ptBins);
	TH1D* hpsi2SPR = new TH1D("hpsi2SPR", ";p_{T};", 6, ptBins);
	TH1D* hpsi2SNP = new TH1D("hpsi2SNP", ";p_{T};", 6, ptBins);

	for(int i=0;i<6;i++)
	{
		valErr yieldJpsi; valErr yieldPsi2S; valErr fracJpsi; valErr fracPsi2S;
		yieldJpsi = getYield_Jpsi(i);
		yieldPsi2S = getYield_psi2S(i);
		fracJpsi = getFrac_Jpsi(i);
		fracPsi2S = getFrac_psi2S(i);

		double err1_Jpsi = yieldJpsi.err;
		double err2_Jpsi = fracJpsi.err;
		double err_JpsiPR = yieldJpsi.val*fracJpsi.val*sqrt((err1_Jpsi/yieldJpsi.val)*(err1_Jpsi/yieldJpsi.val) + (err2_Jpsi/fracJpsi.val)*(err2_Jpsi/fracJpsi.val));
		double err_JpsiNP = yieldJpsi.val*(1-fracJpsi.val)*sqrt((err1_Jpsi/yieldJpsi.val)*(err1_Jpsi/yieldJpsi.val)) + ((err2_Jpsi/(1-fracJpsi.val))*(err2_Jpsi/(1-fracJpsi.val)));

		double err1_Psi2S = yieldPsi2S.err;
		double err2_Psi2S = fracPsi2S.err;
		double err_Psi2SPR = yieldPsi2S.val*fracPsi2S.val*sqrt((err1_Psi2S/yieldPsi2S.val)*(err1_Psi2S/yieldPsi2S.val) + (err2_Psi2S/fracPsi2S.val)*(err2_Psi2S/fracPsi2S.val));
		double err_Psi2SNP = yieldPsi2S.val*(1-fracPsi2S.val)*sqrt((err1_Psi2S/yieldPsi2S.val)*(err1_Psi2S/yieldPsi2S.val)) + ((err2_Psi2S/(1-fracPsi2S.val))*(err2_Psi2S/(1-fracPsi2S.val)));

		hJpsiPR->SetBinContent(i+1,yieldJpsi.val*(1-fracJpsi.val));
		hJpsiPR->SetBinError(i+1,err_JpsiPR);
		hJpsiNP->SetBinContent(i+1,yieldJpsi.val*(fracJpsi.val));
		hJpsiNP->SetBinError(i+1,err_JpsiNP);
		hpsi2SPR->SetBinContent(i+1,yieldPsi2S.val*(1-fracPsi2S.val));
		hpsi2SPR->SetBinError(i+1,err_Psi2SPR);
		hpsi2SNP->SetBinContent(i+1,yieldPsi2S.val*(fracPsi2S.val));
		hpsi2SNP->SetBinError(i+1,err_Psi2SNP);
		cout << "pT	" << ptBins[i] << "	-	" << ptBins[i+1] << "	Jpsi Yield	:	" << yieldJpsi.val <<  " +/- " << yieldJpsi.err << "		frac	:	" << fracJpsi.val  << " +/- " << fracJpsi.err << endl;
		cout << "				" << "psi 2S Yield	:	" << yieldPsi2S.val << " +/- " << yieldPsi2S.err << "		frac	:	" << fracPsi2S.val << " +/- " << fracPsi2S.err <<  endl;
	}

	TH1D* hPR; TH1D* hNP;
	hpsi2SPR->Divide(hJpsiPR);
	hpsi2SNP->Divide(hJpsiNP);

	int text_size=17;
	double text_x = 0.18;
	double text_y = 0.8;
	TCanvas *c1 = new TCanvas("c1","",550,520);
	c1->cd();
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 0.98, 1.0);
	pad1->SetTicks(1,1);
	pad1->Draw();pad1->cd();
	hpsi2SPR->SetMarkerColor(kRed+2);
	hpsi2SPR->SetLineColor(kRed+2);
	hpsi2SPR->SetMarkerStyle(20);
	hpsi2SNP->SetMarkerColor(kBlue+2);
	hpsi2SNP->SetLineColor(kBlue+2);
	hpsi2SNP->SetMarkerStyle(21);

	hpsi2SPR->GetXaxis()->SetLabelSize(0);
	hpsi2SPR->GetXaxis()->SetTitleSize(0);
	hpsi2SPR->GetXaxis()->CenterTitle();
	hpsi2SPR->GetYaxis()->SetRangeUser(0,0.1);
	hpsi2SPR->Draw("PE");

	hpsi2SNP->GetXaxis()->SetLabelSize(0);
	hpsi2SNP->GetXaxis()->SetTitleSize(0);
	hpsi2SNP->GetXaxis()->CenterTitle();
	hpsi2SNP->GetYaxis()->SetRangeUser(0,0.1);
	hpsi2SNP->Draw("PE SAME");

	TLegend *leg = new TLegend(text_x+0.45,text_y-0.2,text_x+0.65,text_y);
	leg->SetTextSize(text_size);
	leg->SetTextFont(43);
	leg->SetBorderSize(0);
	leg->AddEntry(hpsi2SPR,"Prompt","PE");
	leg->AddEntry(hpsi2SNP,"Non Prompt","PE");
	leg->Draw("SAME");

	TPad *pad_A_2 = new TPad("pad_A_2", "pad_A_2", 0, 0.001, 0.98, 0.32);
	c1->cd();
	pad_A_2->Draw();
	pad_A_2->cd();
	pad_A_2->SetTopMargin(0); // Upper and lower plot are joined
	pad_A_2->SetBottomMargin(0.67);
	pad_A_2->SetBottomMargin(0.4);
	pad_A_2->SetFillStyle(4000);
	pad_A_2->SetFrameFillStyle(4000);
	pad_A_2->SetTicks(1,1);

	TH1D* hPull= (TH1D*)hpsi2SNP->Clone("TMP");
	hPull->Divide(hpsi2SPR);
	hPull->SetMarkerSize(0.8);
	hPull->SetMarkerColor(kBlack);
	hPull->SetLineColor(kBlack);
	hPull->SetTitle("");
	hPull->SetTitleSize(0);
	hPull->GetYaxis()->SetTitleOffset(0.3) ;
	hPull->GetYaxis()->SetTitle("Pull") ;
	hPull->GetYaxis()->SetTitleSize(0.08) ;
	hPull->GetYaxis()->SetLabelSize(0.08) ;
	hPull->GetYaxis()->SetRangeUser(-2.8,2.8);
	hPull->GetYaxis()->CenterTitle();

	hPull->GetXaxis()->SetTitleOffset(1.3) ;
	hPull->GetXaxis()->SetLabelOffset(0.04) ;
	hPull->GetXaxis()->SetLabelSize(0.08) ;
	hPull->GetXaxis()->SetTitleSize(0.08) ;
	hPull->GetXaxis()->CenterTitle();

	hPull->GetYaxis()->SetTickSize(0.04);
	hPull->GetYaxis()->SetNdivisions(404);
	hPull->GetXaxis()->SetTickSize(0.03);
	hPull->Draw() ;

	c1->SaveAs("Charmonia_DR_pT.pdf");
}

valErr getYield_Jpsi(int i){
	double ptBins[7] = {3,6.5,9,12,15,20,50};
	TString kineLabel[7];
	if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
	else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
	TFile* inf = new TFile(Form("./Jpsi/roots/2DFit_221116/Mass/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("fitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
valErr getYield_psi2S(int i){
	double ptBins[7] = {3,6.5,9,12,15,20,50};
	TString kineLabel[7];
	if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
	else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
	TFile* inf = new TFile(Form("./psi2S/roots/2DFit_221116/Mass/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("fitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
valErr getFrac_Jpsi(int i) {
	double ptBins[7] = {3,6.5,9,12,15,20,50};
	TString kineLabel[7]; 
	if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
	else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
	TFile* inf = new TFile(Form("./Jpsi/roots/2DFit_221116/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
valErr getFrac_psi2S(int i) {
	double ptBins[7] = {3,6.5,9,12,15,20,50};
	TString kineLabel[7]; 
	if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
	else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
	TFile* inf = new TFile(Form("./psi2S/roots/2DFit_221116/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
/*
   double getFracErr_Jpsi(int i) {
   double ptBins[7] = {3,6.5,9,12,15,20,50};
   TString kineLabel[7]; 
   if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
   else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
   TFile* inf = new TFile(Form("./Jpsi/roots/2DFit_221116/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
   TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");
   double frac;
   frac = fitResults->GetBinError(1);
   return frac;
   }
   double getFrac_psi2S(int i) {
   double ptBins[7] = {3,6.5,9,12,15,20,50};
   TString kineLabel[7]; 
   if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
   else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
   TFile* inf = new TFile(Form("./psi2S/roots/2DFit_221116/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
   TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");
   double frac;
   frac = fitResults->GetBinContent(1);
   return frac;
   }
   double getFracErr_psi2S(int i) {
   double ptBins[7] = {3,6.5,9,12,15,20,50};
   TString kineLabel[7]; 
   if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
   else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
   TFile* inf = new TFile(Form("./psi2S/roots/2DFit_221116/Final/2DFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
   TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");
   double frac;
   frac = fitResults->GetBinError(1);
   return frac;
   }
   double getYield_Jpsi(int i=0) {
   double ptBins[7] = {3,6.5,9,12,15,20,50};
   TString kineLabel[7]; 
   if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
   else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
   TFile* inf = new TFile(Form("./Jpsi/roots/2DFit_221116/Mass/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
   TH1D* fitResults = (TH1D*)inf->Get("fitResults");
   double yield_Jpsi;
   yield_Jpsi = fitResults->GetBinContent(1);
   return yield_Jpsi;
   }
   double getYield_psi2S(int i=0) {
   double ptBins[7] = {3,6.5,9,12,15,20,50};
   TString kineLabel[7]; 
   if(i==0)  kineLabel[i]= getKineLabel(ptBins[i],ptBins[i+1],1.6,2.4,0.0,0,180);
   else kineLabel[i] = getKineLabel(ptBins[i],ptBins[i+1],0,2.4,0.0,0,180);
   TFile* inf = new TFile(Form("./psi2S/roots/2DFit_221116/Mass/Mass_FixedFitResult_%s_PRw_Effw1_Accw1_PtW1_TnP1.root", kineLabel[i].Data()));
   TH1D* fitResults = (TH1D*)inf->Get("fitResults");
   double yield_psi2S;
   yield_psi2S = fitResults->GetBinContent(1);
   return yield_psi2S;
   }*/
