#include "../../rootFitHeaders.h"
#include "../../commonUtility.h"
#include "../../JpsiUtility.h"
#include "../../cutsAndBin.h"
#include "../../CMS_lumi_v2mass.C"
#include "../../tdrstyle.C"
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
#include <iostream>

using namespace std;
using namespace RooFit;

valErr getFrac_psi2S_PbPb(int i=0);
valErr getFrac_psi2S_pp(int i=0);
void draw_Bfraction_y0_1p6()
{
	gStyle->SetOptStat(0);
	setTDRStyle();

	int iPeriod = 101;
	int iPos = 33;

	const int nPtBins=6;
	TFile *fPbPb[nPtBins+1];
	TFile *fpp[nPtBins+1];

	double ptBin[nPtBins+1] = {6.5,9,12,15,20,25,50};
	double fracPP[nPtBins]; double fracPbPb[nPtBins];
	double fracErrPP[nPtBins]; double fracErrPbPb[nPtBins];
	double x[nPtBins]; double binWidth[nPtBins];

	double ptBin_BPH[10]={6.5,8.,9.,10.,11.,12.,13.5,15.,18.,30.};
	double var_BPH[9] = {0.326,0.350,0.335,0.378,0.394,0.431,0.456,0.487,0.572};
	double err_BPH[9] = {0.041,0.025,0.021,0.023,0.024,0.022,0.025,0.024,0.025};
	double x_BPH[9]; double binWidth_BPH[9];

	double BfracPP[nPtBins]; double BfracPbPb[nPtBins];

	TString kineLabel[nPtBins+1];

	TH1D *hfracPP = new TH1D("hfracPP",";p_{T} (GeV/c);non-prompt fraction",nPtBins, ptBin);
	TH1D *hfracPbPb = new TH1D("hfracPbPb",";p_{T} (GeV/c);non-prompt fraction",nPtBins, ptBin);


	for(int i=0;i<nPtBins;i++)
	{
		x[i] = (ptBin[i+1]+ptBin[i])/2;
        binWidth[i] = (ptBin[i+1]-ptBin[i])/2;
		valErr fracPP; valErr fracPbPb;
		fracPP = getFrac_psi2S_pp(i);
		fracPbPb = getFrac_psi2S_PbPb(i);

		BfracPP[i] = fracPP.val; BfracPbPb[i] = fracPbPb.val;

		fracErrPP[i] = fracPP.err;
		fracErrPbPb[i] = fracPbPb.err;
		
		hfracPP->SetBinContent(i+1,fracPP.val);
		hfracPP->SetBinError(i+1,fracPP.err);
		hfracPbPb->SetBinContent(i+1,fracPbPb.val);
		hfracPbPb->SetBinError(i+1,fracPbPb.err);

		//gfracPP->SetPoint(i,(ptBin[i]+ptBin[i+1])/2, fracPP.val);
		//gfracPP->SetPointError(i,i,fracPP.err);
		//gfracPbPb->SetPoint(i,(ptBin[i]+ptBin[i+1])/2, fracPbPb.val);
		//gfracPbPb->SetPointError(i,i,fracPbPb.err);

		cout << "pT	" << ptBin[i] << "	-	" << ptBin[i+1] << "	fraction PP :	" << fracPP.val << " +/- " << fracPP.err << endl;
		cout << "				fraction PbPb :	" << fracPbPb.val << " +/- " << fracPbPb.err << endl;
	}

	for(int i=0; i<9; i++)
	{
		x_BPH[i] = (ptBin_BPH[i+1]+ptBin_BPH[i])/2;
		binWidth_BPH[i] = (ptBin_BPH[i+1]-ptBin_BPH[i])/2;
	}

	TGraphErrors* gfracPP = new TGraphErrors(nPtBins,x,BfracPP,binWidth,fracErrPP); 
	TGraphErrors* gfracPbPb = new TGraphErrors(nPtBins,x,BfracPbPb, binWidth,fracErrPbPb);
	TGraphErrors* gfracBPH = new TGraphErrors(9,x_BPH,var_BPH,binWidth_BPH,err_BPH);

	int text_size=17;
	double text_x = 0.18;
	double text_y = 0.8;
	double y_diff = 0.08;
	TCanvas *c1 = new TCanvas("c1","",550,520);
	c1->cd();
	gfracPP->SetMarkerColor(kBlue+2);
	gfracPP->SetLineColor(kBlue+2);
	gfracPP->SetMarkerStyle(20);
	gfracPbPb->SetMarkerColor(kRed+2);
	gfracPbPb->SetLineColor(kRed+2);
	gfracPbPb->SetMarkerStyle(21);
	gfracBPH->SetMarkerColor(15);
	gfracBPH->SetLineColor(15);
	gfracBPH->SetMarkerStyle(27);
	gfracBPH->SetMarkerSize(1.4);
	
	gfracPP->SetTitle("");
	gfracPP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gfracPP->GetYaxis()->SetTitle("non-propmt fraction");
	gfracPP->GetXaxis()->CenterTitle();
	gfracPP->GetYaxis()->CenterTitle();
	gfracPbPb->GetXaxis()->CenterTitle();
	gfracPbPb->GetYaxis()->CenterTitle();
	gfracPP->GetYaxis()->SetRangeUser(0,0.8);
	gfracPbPb->GetYaxis()->SetRangeUser(0,0.8);
	gfracPP->GetXaxis()->SetLimits(0,50);
	gfracBPH->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	gfracBPH->GetYaxis()->SetTitle("non-propmt fraction");
	gfracBPH->GetXaxis()->CenterTitle();
	gfracBPH->GetYaxis()->CenterTitle();
	gfracBPH->GetXaxis()->SetLimits(0,50);
	gfracBPH->GetYaxis()->SetRangeUser(0,0.75);

	gfracPP->Draw("AP");
	gfracPbPb->Draw("P");
	gfracBPH->Draw("P");

	TLegend *leg = new TLegend(0.2,0.9,0.3,0.8);
	leg->SetTextSize(text_size);
	leg->SetTextFont(43);
	leg->SetBorderSize(0);
	leg->AddEntry(gfracPP,"pp 5.02 TeV, |y|<1.6","PE");
	leg->AddEntry(gfracPbPb,"PbPb 5.02 TeV, |y|<1.6","PE");
	leg->AddEntry(gfracBPH,"pp 7 TeV (BPH-10-014), |y| < 1.2","P");
	leg->Draw("SAME");

	CMS_lumi_v2mass(c1,iPeriod,iPos);

	c1->SaveAs("figs/Bfraction_y0_1p6.pdf");


}
valErr getFrac_psi2S_PbPb(int i) {
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
valErr getFrac_psi2S_pp(int i) {
	double ptBin[7] = {6.5,9,12,15,20,25,50};
	TString kineLabel[7];
	kineLabel[i] = getKineLabelpp(ptBin[i],ptBin[i+1],0,1.6,0.0);
	TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/2DFit_No_Weight/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
	//TFile* inf = new TFile(Form("../pp_psi2S_230512/roots/backup_231122/Final/2DFitResult_%s_PRw_Effw0_Accw0_PtW0_TnP0.root", kineLabel[i].Data()));
	TH1D* fitResults = (TH1D*)inf->Get("2DfitResults");

	valErr ret;
	ret.val = fitResults->GetBinContent(1);
	ret.err = fitResults->GetBinError(1);
	return ret;
}
